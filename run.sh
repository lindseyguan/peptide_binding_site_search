#!/bin/bash

# ------- Foldseek parameters -------
PDB_PATH="/mnt/shared3/PDB/cif" # CHANGE ME: Path to divided PDB mmCIF bioassembly files (e.g., 2ych-assembly1.cif.gz)
FOLDSEEK_DATA_DIR="./foldseek/pdb/pdb" # CHANGE ME: Path to foldseek db created from PDB
FOLDSEEK_EXEC="foldseek" # CHANGE ME: Path to foldseek executable
BASE_OUTPUT_PATH="./example" # CHANGE ME: Base path for outputs (e.g., "../example")

PBS_REPO=$(pwd) # CHANGE ME, OPTIONALLY: Path to peptide_binding_site_search repo.

# ------- mmseqs parameters -------
PDB_SEQRES_PATH="/mnt/shared3/PDB/pdb_seqres.fa" # CHANGE ME: Path to pdb_seqres file in fasta format (see www.rcsb.org/downloads/fasta)
MMSEQS_EXEC="mmseqs" # CHANGE ME: Path to mmseqs executable

FOLDSEEK_INPUT_DIR="${BASE_OUTPUT_PATH}/foldseek/inputs"
FOLDSEEK_OUTPUT_DIR="${BASE_OUTPUT_PATH}/foldseek/outputs"
FOLDSEEK_OUTPUT_PATH="${FOLDSEEK_OUTPUT_DIR}/results"

FASTA_INPUT_PATH="${BASE_OUTPUT_PATH}/mmseqs/inputs/queries.fa"
MMSEQS_OUTPUT_DIR="${BASE_OUTPUT_PATH}/mmseqs/outputs"
MMSEQS_OUTPUT_PATH="${MMSEQS_OUTPUT_DIR}/results"

TMP_DIR="${BASE_OUTPUT_PATH}/tmp"

OUTPUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

# ------- Create output directories -------
echo "Creating output directories..."
mkdir -p "${FOLDSEEK_OUTPUT_DIR}"
mkdir -p "${MMSEQS_OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

# ------- Run! -------
echo "Running foldseek..."

${FOLDSEEK_EXEC} easy-search ${FOLDSEEK_INPUT_DIR} ${FOLDSEEK_DATA_DIR} ${FOLDSEEK_OUTPUT_PATH} ${TMP_DIR} \
		--exhaustive-search --alignment-type 1 -c 0.8 \
		--cov-mode 2 --tmscore-threshold 0.6 \
		--format-output ${OUTPUT_FORMAT}

echo "Done running foldseek!"

# this will add the output columns to the output file, for convenience.
echo -e "${OUTPUT_FORMAT//,/\\t}" | cat - ${FOLDSEEK_OUTPUT_PATH} > ${FOLDSEEK_OUTPUT_PATH}.tmp
mv ${FOLDSEEK_OUTPUT_PATH}.tmp ${FOLDSEEK_OUTPUT_PATH}

echo "Running mmseqs..."

${MMSEQS_EXEC} easy-search ${FASTA_INPUT_PATH} ${PDB_SEQRES_PATH} ${MMSEQS_OUTPUT_PATH} ${TMP_DIR} \
		-c 0.8 --cov-mode 2 \
		--format-output ${OUTPUT_FORMAT}

echo -e "${OUTPUT_FORMAT//,/\\t}" | cat - ${MMSEQS_OUTPUT_PATH} > ${MMSEQS_OUTPUT_PATH}.tmp
mv ${MMSEQS_OUTPUT_PATH}.tmp ${MMSEQS_OUTPUT_PATH}

echo "Done running mmseqs!"

python ${PBS_REPO}/src/search.py ${BASE_OUTPUT_PATH}
