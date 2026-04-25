#!/bin/bash
# Download OsDMP1 + OsDMP3 from NCBI by locus ID (Liang et al., 2025)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Oryza sativa"
LOCUS_LIST=("LOC_Os08g01530" "LOC_Os01g29240")
GENE_NAME_LIST=("OsDMP1" "OsDMP3")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
