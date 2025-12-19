import os
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Model, Structure, Selection, NeighborSearch
from pymol import cmd


# Assumes TMscore is in the path
TMSCORE_PATH = 'TMscore'
DEFAULT_CUTOFF = 6
DEFAULT_INTERACTION_THRESHOLD = 3
DEFAULT_INTERFACE_BB_CUTOFF = 8


# -------- Helper functions --------
def get_residues_id(structure, chain_id):
    """
    Returns the residues of chain_id chain in the structure.
    Ignores het residues.
    """
    model = list(structure.get_models())[0]
    chain = model[chain_id]
    residues = []
    for r in chain.child_list:
        hetflag, resseq, icode = r.get_id()
        if hetflag.strip() != '':
            continue
        residues.append(r)
    return residues


def extract_chains(structure, chain_ids, model_id=0):
    """
    Given a Biopython PDB structure, return a new structure that has just
    the chains specified in `chain_ids` in the order specified. 

    Assumes that there is only one model (it uses the first one) unless
    model is specified as other.
    
    Return: a Biopython structure
    """
    # Extract chains
    chains_list = []
    chains_dict = {}
    model = structure[model_id]
    for chain in chain_ids:
        chains_dict[chain] = model[chain]
        chains_list.append(model[chain])

    # Initialize new model
    new_model = Model.Model(model_id)
    new_model.child_list = chains_list
    new_model.child_dict = chains_dict

    # Initialize new structure
    new_structure = Structure.Structure(structure.get_id())
    new_structure.child_list = [new_model]
    new_structure.child_dict = {model_id: new_model}

    return new_structure


def find_interface_residues(struct, peptide_chain=0, protein_chain=1, cutoff=DEFAULT_CUTOFF):
    """
    Returns a list of residue indices that participate in a binding interface
    between the chains in struct. Peptide assumed to be 0th chain,
    protein assumed to be 1st chain.

    Interaction between two residues is defined as having
    at least one atom `cutoff` distance away from each other.
    """
    # Find interface residues for predicted structure
    if peptide_chain == 0 and protein_chain == 1:
        peptide_residues = get_residues_id(struct, 0)
        protein_residues = get_residues_id(struct, 1)
    else:
        peptide_residues = get_residues_id(struct, peptide_chain)
        protein_residues = get_residues_id(struct, protein_chain)

    prot_interface_residues = set()
    pep_interface_residues = set()

    for prot_res in protein_residues:
        for pep_res in peptide_residues:
            for prot_atom in prot_res.get_atoms():
                for pep_atom in pep_res.get_atoms():
                    if prot_atom - pep_atom < cutoff:
                        prot_interface_residues.add(prot_res)
                        pep_interface_residues.add(pep_res)

    prot_res_ids = sorted([s.get_full_id()[3][1] for s in prot_interface_residues])
    pep_res_ids = sorted([s.get_full_id()[3][1] for s in pep_interface_residues])

    return prot_res_ids, pep_res_ids


def run_tmalign(path1, path2, verbose=False):
    """
    Runs tmalign on path1 and path2, both of which should be PDB files.
    """
    command = f'{TMSCORE_PATH} {path1} {path2} -seq 2> /dev/null'
    output = os.popen(command).read().strip()
    lines = output.split('\n')
    if verbose:
        print(lines)
    return lines[-4:-1]


def extract_pdb(pdb_id, assembly, tmpdir='tmp', pdbdir='/mnt/shared3/PDB/cif', pdbformat='cif'):
    """
    Extracts PDB file from .gz to tmpdir.
    
    Returns the path to the extracted file, or False if the file doesn't exist.
    """
    hash_dir = pdb_id[1:3]
    dest = os.path.join(tmpdir, f'{pdb_id}.{pdbformat}')
    path = os.path.join(pdbdir, f'{hash_dir}/{pdb_id}-assembly{assembly}.{pdbformat}.gz')
    if os.path.isfile(path):
        command = f'gunzip -c {path} > {dest}'
        os.system(command)
        return dest
    else:
        return False

    
def delete_path(path):
    """
    Removes this file if it exists.
    """
    if os.path.isfile(path):
        os.remove(path)

        
def delete_paths(arr):
    """
    Deletes a list of paths.
    """
    for path in arr:
        delete_path(path)
            
            
def get_tm_mapping(tm_alignment):
    """
    Returns mapping between sequence indices given a
    structure-based residue alignment.
    """
    mapping = {}
    seq1 = tm_alignment[0]
    close = tm_alignment[1]
    seq2 = tm_alignment[2]
    seq1_index = 0
    seq2_index = 0

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-' and close[i] == ':':
            mapping[seq1_index] = seq2_index

        if seq1[i] != '-':
            seq1_index += 1
        if seq2[i] != '-':
            seq2_index += 1
    return mapping


def get_seqs_from_tm(tm_alignment):
    """
    Returns original sequences given tm_alignment.
    """
    return tm_alignment[0].replace('-', ''), tm_alignment[2].replace('-', '')


def save_structure(struct, path):
    """
    Saves structure to path using Biopython's PDBIO.
    """
    io = PDBIO()
    io.set_structure(struct)
    io.save(path)
    
    
def all_equal(iterator):
    """
    Check if all elements in an iterator are equal.
    """
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)


def clean(ser):
    """
    Cleans series of boolean values.
    """
    if ser == 'NaN' or ser == 0 or ser == 'False' or (type(ser) != str and np.isnan(ser)):
        return False
    else:
        return True

    
def calc_residue_dist(residue_one, residue_two):
    """
    Returns the C-alpha distance between two residues.
    """
    if 'CA' not in residue_one or 'CA' not in residue_two:
        return np.inf
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two):
    """
    Returns a matrix of C-alpha distances between two chains.
    """
    mat = np.zeros((len(chain_one), len(chain_two)))
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            mat[row, col] = calc_residue_dist(residue_one, residue_two)
    return mat


# -------- Main functions --------
def get_binding_site_hit(query_id, 
                         query_chain, 
                         peptide_chain,
                         target_id, 
                         target_chain, 
                         target_bioassembly,
                         query_path=None,
                         query_chain_only_path=None,
                         interaction_threshold=DEFAULT_INTERACTION_THRESHOLD,
                         tmpdir='./tmp'):
    """
    Analyzes binding site overlap between query and target protein structures.
    
    REQUIRED:
    query_id: PDB ID of query
    query_chain: receptor/protein chain of query
    peptide_chain: peptide chain of query
    target_id: PDB ID of target
    target_chain: receptor/protein chain of target
    target_bioassembly: bioassembly ID of target (this is often 1)
    
    OPTIONAL:
    interaction_threshold: number of contacts to consider two chains to be interacting
    query_path: path to query PDB file
    query_chain_only_path: path to query PDB file containing only the receptor/protein chain
    tmpdir: directory for temporary files
    
    RETURN:
    A dictionary with keys:
        binding_site_overlap: number of overlapping residues between the two binding sites
                             If 0, this binding site is totally novel
                             If a small number (<4), there may be some overlap but it might be flanking residues
                             If a large number, it is likely the "same" binding site
        binding_partner_chains: IDs of binding site partners
        other_chain_ids: IDs of other chains in the target structure
        total_chain_count: total # of chains in the target
        error: whether an error was encountered
        msg: error message, if applicable
    """
    target_id = target_id.lower()

    output_dict = {
        'binding_site_overlap': 0,
        'binding_partner_chains': set(),
        'other_chain_ids': set(), 
        'total_chain_count': 0,
        'error': True,
        'msg': ''
    }
    
    try:
        # Load query structure
        delete_query_tmp = False
        delete_chain_tmp = False
        
        if query_path is None or query_chain_only_path is None:
            if not query_path:
                query_path = extract_pdb(query_id, assembly=1, tmpdir=tmpdir)
                delete_query_tmp = True

            query_chain_only_path = query_path.replace(query_id, f'{query_id}_{query_chain}')
            delete_chain_tmp = True

        model = 0
        if '.cif' in query_path:
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        query_struct = parser.get_structure('_', query_path)
        query_chain_only = extract_chains(query_struct, [query_chain], model_id=model)
        save_structure(query_chain_only, query_chain_only_path)

        query_res = get_residues_id(query_struct, query_chain)

        # maps from residue IDs to positions in sequence
        mapping_from_query_resids = dict(zip([r.id[1] for r in query_res], range(len(query_res))))
        receptor_interface_res = [
            mapping_from_query_resids[r] 
            for r in find_interface_residues(
                query_struct,
                peptide_chain=peptide_chain,
                protein_chain=query_chain
            )[0]
        ]

        if 'pdb' in target_chain:
            if delete_query_tmp:
                delete_path(query_path)
            if delete_chain_tmp:
                delete_path(query_chain_only_path)
            output_dict['msg'] = 'cant_handle_pdb_in_name'
            return output_dict

        elif 'MODEL' in target_chain:
            model = int(target_chain.split('MODEL_')[1].split('_')[0]) - 1
            target_chain = target_chain.split('_')[-1]
            
        # Get target
        try:
            target_path = extract_pdb(target_id, assembly=target_bioassembly, tmpdir=tmpdir)
            if not target_path:
                if delete_query_tmp:
                    delete_path(query_path)
                if delete_chain_tmp:
                    delete_path(query_chain_only_path)
                output_dict['msg'] = 'no_pdb_file_target'
                return output_dict
        except (IOError, OSError) as e:
            if delete_query_tmp:
                delete_path(query_path)
            if delete_chain_tmp:
                delete_path(query_chain_only_path)
            output_dict['msg'] = f'no_pdb_file: {str(e)}'
            return output_dict

        if '.cif' in target_path:
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)            

        target_struct = parser.get_structure('_', target_path)

        # Check number of chains
        target_chains = [c.id for c in target_struct[model].get_chains()]
        total_chains = len(target_chains)
        
        if total_chains == 1:
            delete_path(target_path)
            if delete_query_tmp:
                delete_path(query_path)
            if delete_chain_tmp:
                delete_path(query_chain_only_path)
            output_dict['error'] = False
            return output_dict

        if total_chains == 0:
            delete_path(target_path)
            if delete_query_tmp:
                delete_path(query_path)
            if delete_chain_tmp:
                delete_path(query_chain_only_path)
            output_dict['msg'] = 'no_chains'
            return output_dict

        # Save chain-only structure for use with TMalign
        target_chain_only = extract_chains(target_struct, [target_chain], model_id=model)
        target_chain_only_path = target_path.replace('.pdb', f'_{target_chain}.pdb')
        save_structure(target_chain_only, target_chain_only_path)

        tmalignment = run_tmalign(query_chain_only_path, target_chain_only_path, verbose=False)
        seq1 = tmalignment[0]
        close = tmalignment[1]
        seq2 = tmalignment[2]
        
        if not all_equal([len(seq1), len(close), len(seq2)]):
            delete_paths([target_path, target_chain_only_path])
            if delete_query_tmp:
                delete_path(query_path)
            if delete_chain_tmp:
                delete_path(query_chain_only_path)
            output_dict['msg'] = 'bad_alignment'
            return output_dict

        tmmapping = get_tm_mapping(tmalignment)
        seqs = get_seqs_from_tm(tmalignment)

        target_res = get_residues_id(target_struct, target_chain)
        mapping_to_target_resid = dict(zip(
            range(len(target_res)),
            [r.id[1] for r in target_res]
        ))

        potential_target_interface_res = []
        for r in receptor_interface_res:
            if r in tmmapping:
                potential_target_interface_res.append(mapping_to_target_resid[tmmapping[r]])

        atoms = Selection.unfold_entities(target_struct, 'A')
        ns = NeighborSearch(atoms)
        target_bs_interactions = 0
        interacting_chains = set()

        for r in potential_target_interface_res:
            try:
                target_res_coord = list(target_struct.get_models())[model][target_chain][r]['CA'].get_coord()
            except (KeyError, IndexError):
                continue
                
            close_res = ns.search(target_res_coord, DEFAULT_CUTOFF, 'R')
            close_chains = set([res.get_parent().id for res in close_res])
            close_chains.discard(target_chain)

            if close_chains:
                interacting_chains = interacting_chains | close_chains
                target_bs_interactions += 1

        other_chains = set()
        if len(target_chains) > 2 and len(interacting_chains) > 0:
            # Check for other protein-protein interactions
            potential_other_chains = set(target_chains) - interacting_chains - {target_chain}
            for other_chain in potential_other_chains:
                chain_one = target_struct[model][target_chain]
                chain_two = target_struct[model][other_chain]
                distances = calc_dist_matrix(chain_one, chain_two)
                if sum(np.array(distances).flatten() < DEFAULT_CUTOFF) >= interaction_threshold:
                    other_chains.add(other_chain)

        # Delete temp files
        delete_paths([target_path, target_chain_only_path])
        if delete_query_tmp:
            delete_path(query_path)
        if delete_chain_tmp:
            delete_path(query_chain_only_path)

        output_dict = {
            'binding_site_overlap': target_bs_interactions,
            'binding_partner_chains': interacting_chains,
            'other_chain_ids': other_chains, 
            'total_chain_count': total_chains,
            'error': False,
            'msg': ''
        }
        return output_dict
        
    except Exception as e:
        output_dict['msg'] = f'error: {str(e)}'
        return output_dict


def calculate_ligand_rms(query_id, 
                         query_receptor_chain,
                         query_ligand_chain,
                         target_id, 
                         target_receptor_chain,
                         target_ligand_chain, 
                         target_bioassembly,
                         query_path=None,
                         interface_bb_cutoff=DEFAULT_INTERFACE_BB_CUTOFF,
                         tmpdir='./tmp',
                         verbose=False):
    """
    Calculate RMSD between aligned ligands from query and target structures.
    
    REQUIRED:
    query_id: PDB ID of query
    query_receptor_chain: receptor/protein chain of query
    query_ligand_chain: peptide chain of query
    target_id: PDB ID of target
    target_receptor_chain: receptor/protein chain of target
    target_ligand_chain: peptide chain of target
    target_bioassembly: bioassembly ID of target (this is often 1)
    
    OPTIONAL:
    query_path: path to query PDB file
    interface_bb_cutoff: backbone cutoff distance for interface residues
    tmpdir: path to tmp dir
    verbose: whether to do extra debugging
    
    RETURN:
    A dictionary with keys:
        rmsd: RMSD between aligned ligands
        align_len_atoms: length of alignment in atoms
        rmsd_error: whether an error was encountered
        rmsd_msg: error msg, if applicable
    """
    output_dict = {
        'rmsd': np.nan,
        'align_len_atoms': 0,
        'rmsd_error': True,
        'rmsd_msg': ''
    }
    
    delete_query_tmp = False
    if query_path is None:
        query_path = extract_pdb(query_id, assembly=1, tmpdir=tmpdir)
        if not query_path:
            output_dict['rmsd_msg'] = 'no_pdb_file_query'
            return output_dict
        delete_query_tmp = True

    target_path = extract_pdb(target_id, assembly=target_bioassembly, tmpdir=tmpdir)
    
    if not target_path:
        if delete_query_tmp:
            delete_path(query_path)
        output_dict['rmsd_msg'] = 'no_pdb_file_target'
        return output_dict

    try:
        cmd.delete('all')
        
        if not verbose:
            cmd.feedback("disable", "all", "warnings") 
            cmd.feedback("disable", "all", "actions")
            cmd.feedback("disable", "all", "details")

        cmd.load(target_path, 'raw_target')
        cmd.create('target', f'raw_target and chain {target_receptor_chain}+{target_ligand_chain}')
        cmd.delete('raw_target')

        cmd.load(query_path, 'query')
        cmd.align(f"query and chain {query_receptor_chain}", f"target and chain {target_receptor_chain}")

        cmd.select('query_lig', f'br. (query and chain {query_ligand_chain} and bb.) within {interface_bb_cutoff} of (query and chain {query_receptor_chain} and bb.)')
        cmd.select('target_lig', f'br. (target and chain {target_ligand_chain} and bb.) within {interface_bb_cutoff} of (target and chain {target_receptor_chain} and bb.)')

        atom_length_target_lig = cmd.count_atoms('target_lig')

        if atom_length_target_lig == 0:
            output_dict['rmsd_msg'] = 'no_ligand_atoms_to_align'
            return output_dict

        cmd.super("query_lig", "target_lig", cycles=0, transform=0, object="aln")
        rmsd = cmd.rms_cur("query_lig & aln", "target_lig & aln", matchmaker=-1)
        raw_aln = cmd.get_raw_alignment('aln')

        output_dict = {
            'rmsd': rmsd,
            'align_len_atoms': len(raw_aln),
            'rmsd_error': False,
            'rmsd_msg': ''
        }
        return output_dict
        
    except Exception as e:
        output_dict['rmsd_msg'] = f'error_aligning: {str(e)}'
        return output_dict
    finally:
        # Cleanup temp files
        delete_path(target_path)
        if delete_query_tmp:
            delete_path(query_path)
