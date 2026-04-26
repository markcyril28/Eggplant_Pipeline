#!/bin/bash
# Download StDMP from NCBI by locus ID (J. Zhang et al., 2022)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"
source "$SCRIPT_DIR/_lib_phytozome.sh"

ORGANISM="Solanum tuberosum"
# Soltu.DM.05G005100 is a PGSC locus ID; not indexed in NCBI nuccore directly.
# NCBI Gene DB cross-references it under the gene entry. Phytozome Stuberosum_448_v6.1
# is the definitive source if the gene DB link returns no nuccore records.
LOCUS_LIST=("Soltu.DM.05G005100")
GENE_NAME_LIST=("StDMP")
PHYTOZOME_ORGANISM="Stuberosum"

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    out_file="${GENE_NAME_LIST[i]}_${LOCUS_LIST[i]}.fasta"
    if ncbi_fetch_via_gene_db "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        continue
    fi
    echo "    [INFO] NCBI gene DB miss for ${LOCUS_LIST[i]}; trying Phytozome"
    phytozome_fetch_gene_cds "$PHYTOZOME_ORGANISM" "${GENE_NAME_LIST[i]}" \
        "${LOCUS_LIST[i]//./\\.}|${GENE_NAME_LIST[i]}|DMP|DUF679" "$out_file" || true
done
