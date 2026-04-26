#!/bin/bash
# Download GmDMP1 + GmDMP2 from NCBI by locus ID (Zhong et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Glycine max"
# Glyma.18G IDs are Phytozome v2.1 (Wm82.a2) identifiers; NCBI Gene DB cross-references
# them. Phytozome Gmax_Wm82.a2.v1 is the authoritative source if NCBI misses.
LOCUS_LIST=("Glyma.18G097400" "Glyma.18G098300")
GENE_NAME_LIST=("GmDMP1" "GmDMP2")
PHYTOZOME_ORGANISM="Gmax"

source "$SCRIPT_DIR/_lib_phytozome.sh"

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    out_file="${GENE_NAME_LIST[i]}_${LOCUS_LIST[i]}.fasta"
    if ncbi_fetch_via_gene_db "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        continue
    fi
    echo "    [INFO] NCBI gene DB miss for ${LOCUS_LIST[i]}; trying Phytozome"
    phytozome_fetch_gene_cds "$PHYTOZOME_ORGANISM" "${GENE_NAME_LIST[i]}" \
        "${LOCUS_LIST[i]//./\\.}|${GENE_NAME_LIST[i]}" "$out_file" || true
done
