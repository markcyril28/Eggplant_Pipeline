#!/bin/bash
# Download SlDMP from NCBI by locus ID (Zhong et al., 2022b)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Solanum lycopersicum"
LOCUS_LIST=("Solyc05g007920")
GENE_NAME_LIST=("SlDMP3")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
