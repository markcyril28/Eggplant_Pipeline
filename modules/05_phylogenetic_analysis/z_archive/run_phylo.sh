#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

input_group="h_Final"

mkdir -p logs
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TIME_LOG="logs/phylo_pipeline_${TIMESTAMP}_${input_group}_script_time_log.log"

# Run with timing
/usr/bin/time -v bash "$SCRIPT_DIR/generate_Alignment_and_Phylo.sh" --group "$input_group" \
    --alignment TRUE --phylo TRUE 2>> "$TIME_LOG"

#cd /mnt/c/Users/admon/Pipeline/4_Phylogenetic_Anaylsis
