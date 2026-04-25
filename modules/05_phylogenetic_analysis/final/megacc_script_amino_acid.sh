#!/bin/bash
set -euo pipefail

megacc -a infer_ML_amino_acid_100.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_100.nwks

megacc -a infer_ML_amino_acid_1000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_1000.nwks

megacc -a infer_ML_amino_acid_10000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_10000.nwks

megacc -a infer_ML_amino_acid_100000.mao -d All_SmelDMPs_Closest_DMPs_and_Outgroups_AA_seq.fas -o AA_ML_tree_100000.nwks

: << 'VERSION1'
megacc -a infer_ML_nucleotide.mao -d JRO.fas -o JRO_ML_tree.nwk -n JRO_bootstrap \
	-s K2P+G --GAMMA_CATEGORIES 5 -b 10000 --cpu 3
VERSION1
