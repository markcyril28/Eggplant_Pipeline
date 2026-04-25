#!/bin/bash
# ============================================================================
# Translate nucleotide FASTA to protein using seqkit (frame +1).
# Consolidates the former translate_Bull.sh, translate_Debernardi.sh, and
# translate_combined.sh which differed only in filename prefix.
#
# Usage:
#   bash translate_nuc_to_prot.sh <PREFIX>
#   bash translate_nuc_to_prot.sh Bull
#   bash translate_nuc_to_prot.sh Debernardi
#   bash translate_nuc_to_prot.sh Combined
# ============================================================================

set -euo pipefail

PREFIX="${1:?Usage: $0 <PREFIX>  (e.g. Bull, Debernardi, Combined)}"

NucSeq=(
    "${PREFIX}_Phylo_GIF_NucSeq_with_SmelGIF_e-value_1e-5.fasta"
    "${PREFIX}_Phylo_GRF_NucSeq_with_SmelGRF_e-value_1e-5.fasta"
    "${PREFIX}_Phylo_GRF_NucSeq_with_SmelGRF_e-value_1e-10.fasta"
)

ProtSeq=(
    "${PREFIX}_Phylo_GIF_ProtSeq_with_SmelGIF_e-value_1e-5.fasta"
    "${PREFIX}_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-5.fasta"
    "${PREFIX}_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-10.fasta"
)

for i in "${!NucSeq[@]}"; do
    >"${ProtSeq[$i]}"
    seqkit translate -f 1 "${NucSeq[$i]}" > "${ProtSeq[$i]}"
done
