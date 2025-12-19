

# peptide_binding_site_search

This repo contains code for running a peptide binding site search, as described in Figure 2 of:

Guan L, Keating AE. Training bias and sequence alignments shape protein–peptide docking by AlphaFold and related methods. _Protein Science_. 2025; 34(11):e70331. [https://doi.org/10.1002/pro.70331](https://doi.org/10.1002/pro.70331)


## Prerequisites
Software:

 - [foldseek](https://github.com/steineggerlab/foldseek)
 - [mmseqs](https://github.com/soedinglab/MMseqs2)
 - biopython
 - numpy
 - pandas
 - tqdm
 - [pymol](https://pymol.org/conda/) python API
 - [TMalign](https://aideepmed.com/TM-score/)

Databases:

 - A copy of the [PDB](https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/), as mmCIF bioassembly files divided into hash directories. An rsync script can be found [here](https://files.wwpdb.org/pub/pdb/software/rsyncPDB.sh).
 - A foldseek database of the PDB. You can download a pre-generated one: 
 `foldseek databases PDB pdb tmp`.


 - A copy of PDB [SEQRES](https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz) records.


## Running the script

**Overview**: The `run.sh` script will run mmseqs and foldseek to search for sequence- and structure-based matches for the target chains. It then runs the script to count binding site overlap and calculate peptide RMSD, outputting the results to `DIR/binding_sites` .

**Setup**: Please inspect the top of `run.sh` and see that the paths are consistent with your machine. Also, make sure `TMscore` is in the PATH environment variable (or change the TMscore path in `src/utils.py`). 


The script will expect the input directory to look like this:

```
├── DIR 
│    ├── foldseek
│        └── inputs
│            ├── 2ych_AB.pdb
│            ├── ...
│            └── 8qm0_AE.pdb  
│    └── mmseqs 
│        └── inputs
│            └── queries.fa  
```
`queries.fa` should have headers that match the format used for foldseek inputs, i.e.,`{PDB_ID}_{RECEPTOR_CHAIN}{PEPTIDE_CHAIN}`.

**Run**: `bash ./run.sh`


