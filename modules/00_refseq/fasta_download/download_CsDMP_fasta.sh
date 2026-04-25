#!/bin/bash
# Download CsDMP from NCBI by locus ID (Yin et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Cucumis sativus"
LOCUS_LIST=("CsaV3_1G028660")
GENE_NAME_LIST=("CsDMP")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
