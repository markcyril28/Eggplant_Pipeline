#!/bin/bash
set -euo pipefail

NucSeq=(
    Debernardi_Phylo_GIF_NucSeq_with_SmelGIF_e-value_1e-5.fasta
    Debernardi_Phylo_GRF_NucSeq_with_SmelGRF_e-value_1e-5.fasta
    Debernardi_Phylo_GRF_NucSeq_with_SmelGRF_e-value_1e-10.fasta
)

ProtSeq=(
    Debernardi_Phylo_GIF_ProtSeq_with_SmelGIF_e-value_1e-5.fasta
    Debernardi_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-5.fasta
    Debernardi_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-10.fasta
)

for i in "${!NucSeq[@]}"; do
    >"${ProtSeq[$i]}"
    seqkit translate -f 1 "${NucSeq[$i]}" > "${ProtSeq[$i]}"
done
