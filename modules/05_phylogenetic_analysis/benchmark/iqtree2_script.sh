#!/bin/bash
set -euo pipefail

iqtree2 -s JRO.fas -st DNA -m K2P+G -bb 1000 -alrt 1000 -nt 4 -rcluster 5 --bnni --prefix JRO_iqtree2_bb_100

: << 'EXPLANATION'
Explanation:

    -s JRO.fas → input FASTA file.

    -st DNA → specify that sequences are nucleotides (DNA).

    -m K2P+G → Kimura 2-Parameter model with Gamma-distributed rate variation.

    -bb 1000 → perform 1000 ultrafast bootstrap replicates.

    -alrt 1000 → (optional, but recommended) approximate likelihood-ratio test with 1000 replicates for branch support.

    -nt 64 → use 64 CPU cores.

    -gammacat 5 → set the number of discrete Gamma categories to 5.

EXPLANATION
