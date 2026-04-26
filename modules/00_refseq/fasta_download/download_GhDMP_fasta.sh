#!/bin/bash
# Download GhDMPa + GhDMPd from NCBI by locus ID (Long et al., 2024)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"
source "$SCRIPT_DIR/_lib_phytozome.sh"

ORGANISM="Gossypium hirsutum"
# Gh_A11G / Gh_D11G are CottonGen identifiers (Long et al., 2024) — not in NCBI nuccore directly.
# Phytozome Ghirsutum (TM-1 v3.1) mirrors the same loci but uses the "Gohir.A/D11G..." naming
# with 6-digit suffixes (CottonGen 4-digit IDs are padded with "00"):
#   CottonGen Gh_A11G3045  =  Phytozome Gohir.A11G304500
#   CottonGen Gh_D11G0735  =  Phytozome Gohir.D11G073500
LOCUS_LIST=("Gh_A11G3045" "Gh_D11G0735")
PHYTOZOME_LOCUS_LIST=("Gohir.A11G304500" "Gohir.D11G073500")
GENE_NAME_LIST=("GhDMPa" "GhDMPd")
PHYTOZOME_ORGANISM="Ghirsutum"

for ((i=0; i<${#LOCUS_LIST[@]}; i++)); do
    out_file="${GENE_NAME_LIST[i]}_${LOCUS_LIST[i]}.fasta"
    if ncbi_fetch_via_gene_db "${LOCUS_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        continue
    fi
    echo "    [INFO] NCBI gene DB miss for ${LOCUS_LIST[i]}; trying Phytozome (${PHYTOZOME_LOCUS_LIST[i]})"
    phytozome_fetch_gene_cds "$PHYTOZOME_ORGANISM" "${GENE_NAME_LIST[i]}" \
        "${PHYTOZOME_LOCUS_LIST[i]//./\\.}" "$out_file" || true
done
