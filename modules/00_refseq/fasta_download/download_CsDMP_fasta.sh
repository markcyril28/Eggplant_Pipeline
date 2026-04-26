#!/bin/bash
# Download CsDMP from NCBI by RefSeq mRNA accession (Yin et al., 2024).
#
# Source paper: Yin et al. (2024) Plant Physiology, 194:1282-1290.
#   "Mutating the maternal haploid inducer gene CsDMP in cucumber produces
#    haploids in planta". HIR 0.09-0.40%.
#
# Cited locus is CsaV3_1G028660, a Phytozome/CucurBit Genomics v3 ID on
# chromosome 1. Not registered in NCBI nuccore; previous symbol search
# returned LOC101205404 on chromosome 3 (DMP2-like), which is the wrong gene.
# The single chromosome-1 DUF679 candidate in NCBI's annotation is:
#
#   LOC101222659  XM_004146681.2  chr 1  15616822..15617762  -> CsDMP (translates to DMP9)
#
# NCBI Gene page explicitly states "translates to the protein DMP9" for this
# locus; chromosome and position match the CsaV3_1G028660 mapping.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Cucumis sativus"
ACCESSION_LIST=("XM_004146681.2")
GENE_NAME_LIST=("CsDMP_CsaV3_1G028660_LOC101222659")

for ((i=0; i<${#ACCESSION_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${ACCESSION_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM" || true
done
