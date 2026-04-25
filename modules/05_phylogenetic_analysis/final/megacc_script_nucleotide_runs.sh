#!/bin/bash
set -euo pipefail

megacc -a infer_ML_nucleotide.mao -d JRO.fas -o JRO_ML_tree_100.nwks

: << 'VERSION1'
megacc -a infer_ML_nucleotide.mao -d JRO.fas -o JRO_ML_tree.nwk -n JRO_bootstrap \
	-s K2P+G --GAMMA_CATEGORIES 5 -b 10000 --cpu 3
VERSION1
