#!/bin/bash
# Download ClDMP3 from NCBI by locus ID (Chen et al., 2023)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Citrullus lanatus"
LOCUS_LIST=("Cla97C06G121370")
GENE_NAME_LIST=("ClDMP3")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
