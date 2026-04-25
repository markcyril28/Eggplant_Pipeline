#!/bin/bash
set -euo pipefail

iqtree -s JRO.fas -st DNA -m K2P+G -bb 1000 -alrt 1000 \
	-nt 4 \
	-rcluster 5 --bnni --prefix JRO_iqtree1_bb_1000
