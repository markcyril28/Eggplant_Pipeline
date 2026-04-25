#!/bin/bash
# Download NtDMP1-3 from NCBI by gene symbol (X. Zhang et al., 2022)
# Per-paralog locus IDs not published in the table; query by symbol.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Nicotiana tabacum"
SYMBOL_LIST=("NtDMP1" "NtDMP2" "NtDMP3")
GENE_NAME_LIST=("NtDMP1" "NtDMP2" "NtDMP3")

for ((i=0; i<${#SYMBOL_LIST[@]}; i++)); do
    ncbi_fetch_by_symbol "${SYMBOL_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
