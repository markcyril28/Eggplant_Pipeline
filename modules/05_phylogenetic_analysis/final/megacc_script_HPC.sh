#!/bin/bash
set -euo pipefail

: << 'OFF'
# sudo dpkg -i mega-cc_11.0.13-1_amd64.deb

# Nucleotide
megacc -a infer_ML_nucleotide_100.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_Nuc_seq.fas -o Nuc_ML_tree_100.nwks

megacc -a infer_ML_nucleotide_1000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_Nuc_seq.fas -o Nuc_ML_tree_1000.nwks

megacc -a infer_ML_nucleotide_10000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_Nuc_seq.fas -o Nuc_ML_tree_10000.nwks

megacc -a infer_ML_nucleotide_100000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_Nuc_seq.fas -o Nuc_ML_tree_100000.nwks


# Amino Acid
megacc -a infer_ML_amino_acid_100.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_100.nwks

megacc -a infer_ML_amino_acid_1000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_1000.nwks

megacc -a infer_ML_amino_acid_10000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_10000.nwks

megacc -a infer_ML_amino_acid_100000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_100000.nwks

OFF

megacc -a infer_ML_amino_acid_JTT_10000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_JTT_10000.nwks
