#!/bin/bash
# Download BoC03.DMP9 candidates for Brassica oleracea (Zhao et al., 2022).
#
# Source paper: Zhao et al. (2022) Plant Biotechnology J., 20:2495-2505.
#   "In vivo maternal haploid induction based on genome editing of DMP in
#    Brassica oleracea". HIR 0.41-2.35%.
#
# The Zhao paper does not publish a per-paralog locus ID; it identifies 15
# DMP-like proteins in the B. oleracea genome and selects "BoC03.DMP9" /
# "BoC04.DMP9" by chromosome (C03/C04) and phylogenetic similarity to ZmDMP /
# AtDMP9. Since the paper supplies no locus ID, we resolve via NCBI Gene DB:
# of the 13 DUF679-domain genes annotated in the B. oleracea HDEM RefSeq
# assembly, exactly two map to chromosome C03:
#
#   LOC106333853  XM_013772244.1  chr C3  2007555..2008366    -> primary BoC03.DMP9 candidate
#   LOC106333617  XM_013772041.1  chr C3  56005807..56006620  -> alt C3 candidate
#
# These replace earlier symbol-search hits LOC106304616 (chr C7) and
# LOC106300314 (chr C1), which were not on chromosome C03 and therefore could
# not be the published BoC03.DMP9.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Brassica oleracea"

# Direct accession fetch (nuccore search; bypasses NCBI Gene DB symbol lookup
# which fails because B. oleracea DMP genes are annotated as "uncharacterized
# protein" with only the pfam05078/DUF679 domain attribution).
ACCESSION_LIST=(
    "XM_013772244.1"  # LOC106333853, chr C3 -> BoC03.DMP9 primary
    "XM_013772041.1"  # LOC106333617, chr C3 -> BoC03.DMP9 alternative
)
GENE_NAME_LIST=(
    "BoC03_DMP9"
    "BoC03_DMP9_alt"
)

for ((i=0; i<${#ACCESSION_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${ACCESSION_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM" || true
done
