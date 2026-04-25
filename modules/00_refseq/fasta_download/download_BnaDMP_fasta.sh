#!/bin/bash
# Download BnaDMP 4 paralogs from NCBI (Y. Li et al., 2022)
# Loci are Darmor-bzh v4.1 IDs; resolved via NCBI search.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Brassica napus"
LOCUS_LIST=(
    "BnaA03g55920D"
    "BnaC03g03890D"
    "BnaA04g09480D"
    "BnaC04g31700D"
)
GENE_NAME_LIST=(
    "BnaDMP_A03"
    "BnaDMP_C03"
    "BnaDMP_A04"
    "BnaDMP_C04"
)

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
