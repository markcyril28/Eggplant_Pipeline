#!/bin/bash
set -euo pipefail

timestamp=$(date +"%Y%m%d_%H%M%S")
tar -czf "3_PHYLOGENETIC_TREE_RESULTS_HPC_${timestamp}.tar.gz" 3_PHYLOGENETIC_TREE_RESULTS/
