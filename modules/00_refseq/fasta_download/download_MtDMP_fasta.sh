#!/bin/bash
# Download MtDMP8 + MtDMP9 from NCBI by locus ID (N. Wang et al., 2022)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Medicago truncatula"
LOCUS_LIST=("Medtr7g010890" "Medtr5g044580")
GENE_NAME_LIST=("MtDMP8" "MtDMP9")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
