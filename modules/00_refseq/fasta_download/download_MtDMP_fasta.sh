#!/bin/bash
# Download MtDMP8 + MtDMP9 from NCBI by RefSeq mRNA accession (N. Wang et al., 2022).
#
# Source paper: N. Wang et al. (2022) Plant Biotechnology J., 20:22-24.
#   "In planta haploid induction by genome editing of DMP in the model legume
#    Medicago truncatula". HIR 0.29-0.82%.
#
# The cited Mt4.0v1 locus IDs (Medtr7g010890, Medtr5g044580) are not registered
# in NCBI nuccore; searching by them fell through to whole-chromosome WGS
# records (~50 MB each). NCBI Gene has them under MtrunA17r5.0-ANR aliases:
#
#   Medtr7g010890  ->  LOC11409144 (MtrunA17_Chr7g0217511)  XM_003621193.2  protein DMP9 (chr 7)
#   Medtr5g044580  ->  LOC11431589 (MtrunA17_Chr5g0419161)  XM_003614037.1  protein DMP9 (chr 5)
#
# The Wang paper labels them MtDMP8 / MtDMP9 by phylogenetic analogy with
# AtDMP8 / AtDMP9; NCBI annotates both as "protein DMP9" because they cluster
# in the AtDMP9 clade.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Medicago truncatula"
ACCESSION_LIST=(
    "XM_003621193.2"  # LOC11409144 / MtrunA17_Chr7g0217511 / Medtr7g010890
    "XM_003614037.1"  # LOC11431589 / MtrunA17_Chr5g0419161 / Medtr5g044580
)
GENE_NAME_LIST=(
    "MtDMP8_Medtr7g010890_LOC11409144"
    "MtDMP9_Medtr5g044580_LOC11431589"
)

for ((i=0; i<${#ACCESSION_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${ACCESSION_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM" || true
done
