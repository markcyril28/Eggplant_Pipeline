#!/bin/bash
# Download BjuDMP1-4 candidates for Brassica juncea (Chu et al., 2025).
#
# Source paper: Chu et al. (2025) Horticulture Research, 12(7):uhaf094.
#   "In vivo maternal haploid induction in Brassica juncea". HIR 0.64-1.51%.
#
# Data-availability gap: the Chu paper releases BjuDMP1-4 sequences only via a
# Chinese cloud-storage link (https://kdocs.cn/l/cpGGcwjHItGV) which requires
# manual sign-in and cannot be machine-fetched. NCBI Gene contains zero
# DUF679/DMP entries for "Brassica juncea"[Organism] (verified 2026-04-26).
# Until the Chu sequences can be obtained manually and dropped into this dir,
# we proxy via the B. juncea AA subgenome donor (Brassica rapa, NCBI annotated)
# using the four DMP9/DMP7/DMP2/DMP3 paralogs most likely to match the Chu
# nomenclature, plus the second B. rapa DMP9 paralog so the HMMER stage has the
# full canonical-HI candidate panel:
#
#   LOC103855829  XM_009132872.3  B. rapa  chr A03  protein DMP9   -> BjuDMP1 (canonical HI)
#   LOC103863904  XM_*            B. rapa  chr A04  protein DMP9   -> BjuDMP1_alt (canonical HI)
#   LOC103861835  XM_009139540.3  B. rapa  chr A03  protein DMP7   -> BjuDMP2
#   LOC103835095  XM_009111198.3  B. rapa  chr A01  protein DMP2   -> BjuDMP3
#   LOC103861987  XM_033288643.1  B. rapa  chr A01  protein DMP3   -> BjuDMP4
#
# B. nigra (BB subgenome donor) has no DUF679/DMP genes in NCBI Gene, so no BB
# proxy is available. To replace these AA-only proxies with the true Chu
# sequences, manually download the kdocs amino-acid file and place per-paralog
# FASTAs here named BjuDMP1.fasta..BjuDMP4.fasta (then re-run with overwrite=false).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

# AA-subgenome (Brassica rapa) canonical DMP candidates as proxies.
ORGANISM_PROXY="Brassica rapa"
ACCESSION_LIST=(
    "XM_009132872.3"  # LOC103855829, B. rapa chrA03, DMP9
    "XM_009139540.3"  # LOC103861835, B. rapa chrA03, DMP7
    "XM_009111198.3"  # LOC103835095, B. rapa chrA01, DMP2
    "XM_033288643.1"  # LOC103861987, B. rapa chrA01, DMP3
)
GENE_NAME_LIST=(
    "BjuDMP1_proxy_BraDMP9_A03"
    "BjuDMP2_proxy_BraDMP7_A03"
    "BjuDMP3_proxy_BraDMP2_A01"
    "BjuDMP4_proxy_BraDMP3_A01"
)

for ((i=0; i<${#ACCESSION_LIST[@]}; i++)); do
    ncbi_fetch_by_locus "${ACCESSION_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM_PROXY" || true
done

# Add the second B. rapa DMP9 paralog (chr A04) for full canonical-HI panel.
ncbi_fetch_via_gene_db "103863904" "BjuDMP1alt_proxy_BraDMP9_A04" "$ORGANISM_PROXY" || true

# Optional: try a direct B. juncea symbol search; harmless when it returns nothing.
ORGANISM="Brassica juncea"
SYMBOL_LIST=("BjuDMP1" "BjuDMP2" "BjuDMP3" "BjuDMP4")
for sym in "${SYMBOL_LIST[@]}"; do
    ncbi_fetch_via_gene_db "$sym" "$sym" "$ORGANISM" || true
done
