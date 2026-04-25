#!/bin/bash
# Download StDMP from NCBI by locus ID (J. Zhang et al., 2022)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Solanum tuberosum"
LOCUS_LIST=("Soltu.DM.05G005100")
GENE_NAME_LIST=("StDMP")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
