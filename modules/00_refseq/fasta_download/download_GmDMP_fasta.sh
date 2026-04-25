#!/bin/bash
# Download GmDMP1 + GmDMP2 from NCBI by locus ID (Zhong et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Glycine max"
LOCUS_LIST=("Glyma.18G097400" "Glyma.18G098300")
GENE_NAME_LIST=("GmDMP1" "GmDMP2")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
