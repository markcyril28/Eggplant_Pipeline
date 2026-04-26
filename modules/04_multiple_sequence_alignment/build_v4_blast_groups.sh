#!/bin/bash
# ============================================================================
# Build pre-curated MSA Selected/v4_BLAST_Groups_bitscore<N> input FASTAs
# ============================================================================
# For each anchor SmelDMP paralog, this script writes:
#   <out_dir>/<paralog>_NUCLEOTIDE_Sequence.fasta   (anchor + blastn hits >= threshold)
#   <out_dir>/<paralog>_AMINO_ACID_Sequence.fasta   (anchor + blastp hits >= threshold)
#
# Configuration (read from config/<GROUP>/04_multiple_sequence_alignment.toml):
#   [blast_groups_msa]
#     bitscore_threshold = 200
#     paralog_ids        = [...]
#     output_dir         = "04_MSA/merged_input/Selected/v4_BLAST_Groups_bitscore200"
#
# Source CSVs are auto-discovered as
#   III_RESULT/<GROUP>/02_BLAST_Alignment/<GENOME>/curated_results/merged_blast{n,p}_*_plant_only.csv
#
# Usage:
#   bash build_v4_blast_groups.sh                         # DMP, GPE001970, defaults from config
#   bash build_v4_blast_groups.sh --gene-group DMP --bitscore 250
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

GENE_GROUP="${GENE_GROUP:-DMP}"
GENOME="${GENOME:-GPE001970_SMEL5}"
BITSCORE="${BITSCORE:-}"
TYPES="${TYPES:-nucleotide,amino_acid}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gene-group)   GENE_GROUP="$2"; shift 2 ;;
        --genome)       GENOME="$2"; shift 2 ;;
        --bitscore)     BITSCORE="$2"; shift 2 ;;
        --types)        TYPES="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,/^# =====/p' "${BASH_SOURCE[0]}" | sed 's/^# *//'
            exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

ARGS=(--pipeline-dir "$PIPELINE_DIR" --gene-group "$GENE_GROUP" --genome "$GENOME" --types "$TYPES")
[[ -n "$BITSCORE" ]] && ARGS+=(--bitscore "$BITSCORE")

python3 "$SCRIPT_DIR/build_v4_blast_groups.py" "${ARGS[@]}"
