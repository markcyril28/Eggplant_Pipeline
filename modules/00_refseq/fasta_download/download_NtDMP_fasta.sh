#!/bin/bash
# Download NtDMP candidates for Nicotiana tabacum (X. Zhang et al., 2022).
#
# Source paper: X. Zhang et al. (2022) Planta, 255:98.
#   "Haploid induction in allotetraploid tobacco using DMPs mutation". HIR 1.52-1.75%.
#
# The Zhang paper names three CRISPR-targeted homologs "NtDMP1", "NtDMP2",
# "NtDMP3" but does not publish per-paralog locus IDs; the canonical sequences
# are listed only in the corresponding patent CN117285610A (SEQ ID NO: 5/7/9).
# Because we cannot machine-fetch the patent sequence listing, we download the
# full N. tabacum DMP-family panel from NCBI's RefSeq annotation (allowing
# downstream HMMER/BLAST stages to identify the three NtDMP1-3 candidates):
#
#   LOC107816576  XM_016642301.2  protein DMP3-like
#   LOC107767355  XM_016586331.1  protein DMP2-like
#   LOC107760505  XM_016578562.2  protein DMP4-like
#   LOC107783066  XM_016604032.2  protein DMP9-like  (canonical AtDMP9 ortholog -- most likely true HI gene)
#   LOC107772236  XM_016591727.2  protein DMP7-like
#   LOC107824786  protein DMP6-like
#   LOC107815137  protein DMP6-like
#
# The DMP9-like and DMP7-like homologs (added in this version) are the closest
# orthologs of the canonical AtDMP8/AtDMP9 haploid-inducer pair and are the
# most plausible true identities of the Zhang paper's NtDMP1-3.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Nicotiana tabacum"

ACCESSION_LIST=(
    "XM_016642301.2"  # LOC107816576, DMP3-like  (legacy NtDMP1)
    "XM_016586331.1"  # LOC107767355, DMP2-like  (legacy NtDMP2)
    "XM_016578562.2"  # LOC107760505, DMP4-like  (legacy NtDMP3)
    "XM_016604032.2"  # LOC107783066, DMP9-like  (canonical HI candidate)
    "XM_016591727.2"  # LOC107772236, DMP7-like
)
GENE_NAME_LIST=(
    "NtDMP3like"
    "NtDMP2like"
    "NtDMP4like"
    "NtDMP9like"
    "NtDMP7like"
)

for ((i=0; i<${#ACCESSION_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${ACCESSION_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM" || true
done

# Also try Gene DB symbol search in case future RefSeq annotation registers
# NtDMP1/2/3 directly. Harmless if it returns no hits.
SYMBOL_LIST=("NtDMP1" "NtDMP2" "NtDMP3")
for sym in "${SYMBOL_LIST[@]}"; do
    ncbi_fetch_via_gene_db "$sym" "$sym" "$ORGANISM" || true
done
