#!/bin/bash
# Download GhDMPa + GhDMPd (cotton DMP; AtDMP8/9 haploid-induction clade) from NCBI RefSeq.
# -----------------------------------------------------------------------------
# Corrected 2026-06-07. The previous version mapped the CottonGen IDs (Long et al.,
# 2024) Gh_A11G3045 / Gh_D11G0735 to Phytozome by padding the 4-digit suffix with
# "00" (-> Gohir.A11G304500 / Gohir.D11G073500). That assumption was wrong: the
# fetched models were not the DMP gene (829 aa / 431 aa; zero 12-mer overlap with
# any of the 47 NCBI cotton DUF679 proteins), producing long spurious tree branches.
# The correct clean single-DUF679 homeologs are the NCBI RefSeq "protein DMP9" pair
# on chromosomes A11 (A-subgenome = GhDMPa) and D11 (D-subgenome = GhDMPd), 218 aa
# each, confirmed against Long et al. 2024 (~59.3% identity to AtDMP8/9).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"

ORGANISM="Gossypium hirsutum"
# index 0 = GhDMPa (A-subgenome, chr A11); index 1 = GhDMPd (D-subgenome, chr D11)
ACC_LIST=("XM_016854833.2" "XM_016839743.2")            # RefSeq mRNA (CDS)
PROT_LIST=("XP_016710322.1" "XP_016695232.1")           # encoded protein
LOCUS_LIST=("LOC107924398" "LOC107911807")
SUBGENOME_LIST=("A-subgenome chrA11 Gh_A11G3045" "D-subgenome chrD11 Gh_D11G0735")
GENE_NAME_LIST=("GhDMPa" "GhDMPd")

for ((i=0; i<${#ACC_LIST[@]}; i++)); do
    acc="${ACC_LIST[i]}"; gene="${GENE_NAME_LIST[i]}"
    out_file="${gene}_${acc%.*}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"; continue
    fi
    echo ">>> $gene  ($acc)  in  $ORGANISM  [RefSeq CDS]"
    url="$(_eutils_url efetch.fcgi "db=nuccore&id=${acc}&rettype=fasta_cds_na&retmode=text")"
    tmp="$(mktemp)"
    wget -qO "$tmp" --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$url"
    {
      echo ">${gene} ${acc} Gossypium hirsutum protein DMP9 ${SUBGENOME_LIST[i]} CDS [refseq_protein=${PROT_LIST[i]} gene=${LOCUS_LIST[i]}]"
      grep -v '^>' "$tmp"
    } > "$out_file"
    rm -f "$tmp"
    [[ -s "$out_file" ]] && echo "    -> $out_file" || echo "    [ERROR] efetch failed for $acc" >&2
    sleep "$EUTILS_DELAY"
done
