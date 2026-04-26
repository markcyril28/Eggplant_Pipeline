#!/bin/bash
# Download CsDMP from NCBI by locus ID (Yin et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Cucumis sativus"
# CsaV3_1G028660 is a Phytozome/CucurBit Genomics v3 locus ID; not in NCBI nuccore.
# NCBI Gene DB may carry it as a cross-reference; Phytozome is the definitive source.
LOCUS_LIST=("CsaV3_1G028660")
GENE_NAME_LIST=("CsDMP")
PHYTOZOME_ORGANISM="Csativus"
PHYTOZOME_PATTERN="CsaV3_1G028660|CsDMP|DMP|DUF679"

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
