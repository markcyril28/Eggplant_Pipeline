#!/bin/bash
# Download StDMP (potato DMP; AtDMP8/9 haploid-induction clade) from NCBI RefSeq.
# -----------------------------------------------------------------------------
# Corrected 2026-06-07. The previous version resolved the PGSC/Phytozome locus
# Soltu.DM.05G005100 and fetched the DM v6.1 primary-transcript model, which is an
# over-extended fusion (529 aa = a 228-aa DUF679 DMP + ~301 aa of unrelated
# C-terminal sequence) and produced a long, misleading branch in the DMP tree.
# The clean curated single-DUF679 CDS for the same gene is the NCBI RefSeq model
# below (228 aa). Cross-references (UniProt M1DKT1): Soltu.DM.05G005100 =
# PGSC0003DMG400040190 = LOC102591030 = XP_006355155.1 (Zhang et al. 2022).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Solanum tuberosum"
ACC="XM_006355093.1"        # RefSeq mRNA (CDS -> protein XP_006355155.1)
GENE_NAME="StDMP"
out_file="${GENE_NAME}_${ACC%.*}.fasta"

if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
    echo "    [SKIP] $out_file exists"
    exit 0
fi
echo ">>> $GENE_NAME  ($ACC)  in  $ORGANISM  [RefSeq CDS]"
url="$(_eutils_url efetch.fcgi "db=nuccore&id=${ACC}&rettype=fasta_cds_na&retmode=text")"
tmp="$(mktemp)"
wget -qO "$tmp" --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$url"
# Rewrite the lcl|... CDS header to a label-first header (first token = GENE_NAME)
{
  echo ">${GENE_NAME} ${ACC} Solanum tuberosum DMP CDS [refseq_protein=XP_006355155.1 locus=Soltu.DM.05G005100 gene=LOC102591030 uniprot=M1DKT1]"
  grep -v '^>' "$tmp"
} > "$out_file"
rm -f "$tmp"
[[ -s "$out_file" ]] && echo "    -> $out_file" || { echo "    [ERROR] efetch failed for $ACC" >&2; exit 1; }
