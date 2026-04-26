#!/bin/bash
# Download OsDMP1 + OsDMP3 from NCBI by locus ID (Liang et al., 2025)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"
source "$SCRIPT_DIR/_lib_phytozome.sh"

ORGANISM="Oryza sativa"
# LOC_Os IDs are MSU7 / RAP-DB identifiers. NCBI Gene DB carries them as
# cross-references; Phytozome Osativa (v7.0) is the authoritative source.
LOCUS_LIST=("LOC_Os08g01530" "LOC_Os01g29240")
GENE_NAME_LIST=("OsDMP1" "OsDMP3")
PHYTOZOME_ORGANISM="Osativa"

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    out_file="${GENE_NAME_LIST[i]}_${LOCUS_LIST[i]}.fasta"
    if ncbi_fetch_via_gene_db "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        continue
    fi
    echo "    [INFO] NCBI gene DB miss for ${LOCUS_LIST[i]}; trying Phytozome"
    phytozome_fetch_gene_cds "$PHYTOZOME_ORGANISM" "${GENE_NAME_LIST[i]}" \
        "${LOCUS_LIST[i]}|${GENE_NAME_LIST[i]}" "$out_file" || true
done
