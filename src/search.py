import argparse
import math
import os
import pickle
import sys

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.PDB import *
from Bio.Seq import Seq
from tqdm import tqdm
import utils


# Create the parser
parser = argparse.ArgumentParser(
    description='Aligns and searches for peptide binding sites based on foldseek/mmseqs target matches.'
)

# Positional argument (required)
parser.add_argument(
    'input_dir',
    help='Directory containing foldseek and mmseqs outputs. Results will be written to {input_dir}.'
)

# Parse the arguments
args = parser.parse_args()


def main():
    # Use the arguments
    input_dir = args.input_dir

    seq_hits = pd.read_csv(os.path.join(input_dir, 'mmseqs/outputs/results'), sep='\s+')
    seq_hits['target_chain'] = seq_hits['target'].str.split('_').str[-1]
    seq_hits['target_pdb_id'] = seq_hits['target'].str.split('_').str[0]
    seq_hits['query_chain'] = seq_hits['query'].str.split('_').str[-1].str[0]
    seq_hits['query_id'] = seq_hits['query'].str.split('_').str[0] + '_' + seq_hits['query_chain']
    seq_hits['query_receptor_chain'] = seq_hits['query'].str.split('_').str[1].str[0]
    seq_hits['query_peptide_chain'] = seq_hits['query'].str.split('_').str[1].str[1]
    seq_hits['target_id'] = seq_hits['target']

    struct_hits = pd.read_csv(os.path.join(input_dir, 'foldseek/outputs/results'), sep='\s+')
    struct_hits = struct_hits[struct_hits['evalue'] >= 0.6] # Filter for TM-score
    struct_hits['assembly'] = struct_hits['target'].str.split('-').str[1].str.split('.').str[0].str.split('assembly').str[-1].astype(int)
    struct_hits['target_pdb_id'] = struct_hits['target'].str.split('-').str[0]
    struct_hits['target_chain'] = struct_hits['target'].str.split('_').str[-1]
    struct_hits['query_chain'] = struct_hits['query'].str.split('_').str[-1]
    struct_hits['query_receptor_chain'] = struct_hits['query'].str.split('_').str[1].str[0]
    struct_hits['query_peptide_chain'] = struct_hits['query'].str.split('_').str[1].str[1]
    struct_hits['query_id'] = struct_hits['query'].str.split('_').str[0] + '_' + struct_hits['query_chain']
    struct_hits['target_id'] = struct_hits['target_pdb_id'] + '_' + struct_hits['target_chain']

    hits = pd.concat([seq_hits, struct_hits])
    hits = hits[['query_id', 
                  'target_id', 
                  'query_receptor_chain', 
                  'query_peptide_chain',
                  'query_chain',
                  'target_pdb_id',
                  'target_chain', 
                  'assembly']]

    for col in ['query_id', 
                  'target_id', 
                  'query_receptor_chain', 
                  'query_peptide_chain',
                  'query_chain',
                  'target_pdb_id',
                  'target_chain']:
        hits[col] = hits[col].astype(str)

    hits['query_full_id'] = hits['query_id'].str.split('_').str[0] + '_' + \
                            hits['query_receptor_chain'] + hits['query_peptide_chain']
    hits['assembly'] = hits['assembly'].fillna(1).astype(int)
    hits = hits.drop_duplicates()

    bs_output_root = os.path.join(input_dir, 'binding_sites')
    if not os.path.isdir(bs_output_root):
        os.mkdir(bs_output_root)

    output_root = os.path.join(bs_output_root, 'outputs')
    if not os.path.isdir(output_root):
        os.mkdir(output_root)

    main_output = []

    for query_full_id in hits['query_full_id'].unique():
        pdb_id = query_full_id[:4]
        receptor_chain = query_full_id[-2]
        peptide_chain = query_full_id[-1]
        query_id = pdb_id + '_' + receptor_chain
        
        hit_rows = hits[hits['query_id'] == query_id]
        hit_rows = hit_rows.to_dict(orient='records')
        
        output_csv_path = os.path.join(output_root, f'{pdb_id}_{receptor_chain}.csv')
        
        for row in tqdm(hit_rows, desc=query_full_id):
            target_id = row['target_pdb_id']
            target_chain = str(row['target_chain'])
                
            assembly = row['assembly']
            if not assembly:
                assembly = 1

            if '-' in target_chain:
                target_chain = target_chain.split('-')[0]
            
            # Note that we are NOT skipping instances where the target is the same as the query
            try:
                query_path = os.path.join(input_dir, 'foldseek/inputs'+f'/{query_full_id}.pdb')
                output = utils.get_binding_site_hit(pdb_id, receptor_chain, peptide_chain, 
                                                      target_id, target_chain,
                                                      target_bioassembly=int(assembly),
                                                      query_path=query_path,
                                                      tmpdir=os.path.join(input_dir, 'tmp')
                                                     )

                row.update(output)

                if row['binding_site_overlap'] > 0 and not row['error']:
                    binding_partner_chains = row['binding_partner_chains']
                    for binding_partner in binding_partner_chains:
                        rmsd = utils.calculate_ligand_rms(query_id=pdb_id, 
                                                          query_receptor_chain=receptor_chain,
                                                          query_ligand_chain=peptide_chain,
                                                          target_id=target_id, 
                                                          target_receptor_chain=target_chain,
                                                          target_ligand_chain=binding_partner, 
                                                          target_bioassembly=int(assembly),
                                                          query_path=query_path,
                                                          tmpdir=os.path.join(input_dir, 'tmp'),
                                                          verbose=False
                                                         )
                        new_row = row.copy()
                        new_row['target_peptide_chain'] = binding_partner
                        new_row.pop('binding_partner_chains', None)
                        new_row['other_chain_ids'] = '|'.join(new_row['other_chain_ids'])

                        new_row.update(rmsd)
                        main_output.append(new_row)
            except:
                row['major_error'] = True
                
    summary_df = pd.DataFrame(main_output).drop(columns=['query_chain', 'target_chain'])
    summary_df.to_csv(os.path.join(bs_output_root, 'summary.csv'), index=None)


if __name__ == '__main__':
    main()
