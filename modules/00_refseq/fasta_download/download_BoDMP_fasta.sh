#!/bin/bash
# Download BoC03.DMP9 from NCBI by gene symbol (Zhao et al., 2022)
# No published per-paralog locus ID in the HIR table; query by symbol.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Brassica oleracea"
SYMBOL_LIST=("BoDMP9" "DMP9")
GENE_NAME_LIST=("BoC03_DMP9" "BoDMP9_alt")

for ((i=0; i<${#SYMBOL_LIST[@]}; i++)); do
    ncbi_fetch_by_symbol "${SYMBOL_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"
done
