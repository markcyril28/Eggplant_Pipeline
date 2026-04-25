#!/bin/bash
# Download BjuDMP1-4 from NCBI by gene symbol (Chu et al., 2025)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Brassica juncea"
SYMBOL_LIST=("BjuDMP1" "BjuDMP2" "BjuDMP3" "BjuDMP4")
GENE_NAME_LIST=("BjuDMP1" "BjuDMP2" "BjuDMP3" "BjuDMP4")

for ((i=0; i<${#SYMBOL_LIST[@]}; i++)); do
    ncbi_fetch_by_symbol "${SYMBOL_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
