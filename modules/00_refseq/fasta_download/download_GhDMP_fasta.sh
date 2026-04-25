#!/bin/bash
# Download GhDMPa + GhDMPd from NCBI by locus ID (Long et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Gossypium hirsutum"
LOCUS_LIST=("Gh_A11G3045" "Gh_D11G0735")
GENE_NAME_LIST=("GhDMPa" "GhDMPd")

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
