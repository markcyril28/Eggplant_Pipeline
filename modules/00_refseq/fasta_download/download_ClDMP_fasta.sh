#!/bin/bash
# Download ClDMP3 from NCBI by locus ID (Chen et al., 2023)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Citrullus lanatus"
# Cla97C06G121370 is a CucurBit Genomics DB locus ID; not in NCBI nuccore.
# NCBI Gene DB indexes it via cross-reference; fall back to Phytozome if needed.
LOCUS_LIST=("Cla97C06G121370")
GENE_NAME_LIST=("ClDMP3")
PHYTOZOME_ORGANISM="Clanatus"
PHYTOZOME_PATTERN="Cla97C06G121370|ClDMP3|DMP|DUF679"

source "$SCRIPT_DIR/_lib_phytozome.sh"

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    out_file="${GENE_NAME_LIST[i]}_${LOCUS_LIST[i]}.fasta"
    if ncbi_fetch_via_gene_db "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        continue
    fi
    echo "    [INFO] NCBI gene DB miss for ${LOCUS_LIST[i]}; trying Phytozome"
    phytozome_fetch_gene_cds "$PHYTOZOME_ORGANISM" "${GENE_NAME_LIST[i]}" \
        "$PHYTOZOME_PATTERN" "$out_file" || true
done
