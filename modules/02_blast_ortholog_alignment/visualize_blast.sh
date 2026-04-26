#!/bin/bash
# ============================================================================
# Module: BLASTn Visualization
# ============================================================================
# Thin wrapper that invokes visualize_blast_results.py for one curated-results
# directory.  Called by 02_blast_alignment.sh after BLASTn CSVs are generated.
#
# Usage:
#   bash visualize_blast.sh \
#       --results-dir /path/to/curated_results \
#       --top-n 10 \
#       --figures heatmap,heatmap_evalue,lollipop
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RESULTS_DIR=""
GENE_GROUP=""
TOP_N="10"
FIGURES="heatmap,heatmap_evalue,lollipop"
COLORMAP="RdYlGn"
FIGURE_DPI="150"
SAVE_DPI="300"
HEATMAP_VMIN="65.0"
HEATMAP_VMAX="100.0"
HEATMAP_W_SCALE="1.20"
HEATMAP_H_SCALE="0.92"
LOLLIPOP_NCOLS="2"
LOLLIPOP_X_PAD="1.60"
LOLLIPOP_DOT_SIZE="100"
LOLLIPOP_DOT_SIZE_HI="150"
HI_STEM_COLOR="#c7920a"
STEM_COLOR="#d1d5db"
HI_LABELS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir)          RESULTS_DIR="$2";        shift 2 ;;
        --gene-group)           GENE_GROUP="$2";         shift 2 ;;
        --top-n)                TOP_N="$2";              shift 2 ;;
        --figures)              FIGURES="$2";            shift 2 ;;
        --colormap)             COLORMAP="$2";           shift 2 ;;
        --figure-dpi)           FIGURE_DPI="$2";         shift 2 ;;
        --save-dpi)             SAVE_DPI="$2";           shift 2 ;;
        --heatmap-vmin)         HEATMAP_VMIN="$2";       shift 2 ;;
        --heatmap-vmax)         HEATMAP_VMAX="$2";       shift 2 ;;
        --heatmap-w-scale)      HEATMAP_W_SCALE="$2";    shift 2 ;;
        --heatmap-h-scale)      HEATMAP_H_SCALE="$2";    shift 2 ;;
        --lollipop-ncols)       LOLLIPOP_NCOLS="$2";     shift 2 ;;
        --lollipop-x-pad)       LOLLIPOP_X_PAD="$2";     shift 2 ;;
        --lollipop-dot-size)    LOLLIPOP_DOT_SIZE="$2";  shift 2 ;;
        --lollipop-dot-size-hi) LOLLIPOP_DOT_SIZE_HI="$2"; shift 2 ;;
        --hi-stem-color)        HI_STEM_COLOR="$2";      shift 2 ;;
        --stem-color)           STEM_COLOR="$2";         shift 2 ;;
        --hi-labels)            HI_LABELS="$2";          shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$RESULTS_DIR" ]] && { echo "ERROR: --results-dir is required"; exit 1; }
[[ ! -d "$RESULTS_DIR" ]] && { echo "ERROR: directory does not exist: $RESULTS_DIR"; exit 1; }

_GENE_GROUP_ARG=()
[[ -n "$GENE_GROUP" ]] && _GENE_GROUP_ARG=(--gene-group "$GENE_GROUP")

_HI_LABELS_ARG=()
[[ -n "$HI_LABELS" ]] && _HI_LABELS_ARG=(--hi-labels "$HI_LABELS")

python3 "$SCRIPT_DIR/visualize_blast_results.py" \
    --results-dir "$RESULTS_DIR" \
    "${_GENE_GROUP_ARG[@]}" \
    --top-n "$TOP_N" \
    --figures "$FIGURES" \
    --colormap "$COLORMAP" \
    --figure-dpi "$FIGURE_DPI" \
    --save-dpi "$SAVE_DPI" \
    --heatmap-vmin "$HEATMAP_VMIN" \
    --heatmap-vmax "$HEATMAP_VMAX" \
    --heatmap-w-scale "$HEATMAP_W_SCALE" \
    --heatmap-h-scale "$HEATMAP_H_SCALE" \
    --lollipop-ncols "$LOLLIPOP_NCOLS" \
    --lollipop-x-pad "$LOLLIPOP_X_PAD" \
    --lollipop-dot-size "$LOLLIPOP_DOT_SIZE" \
    --lollipop-dot-size-hi "$LOLLIPOP_DOT_SIZE_HI" \
    --hi-stem-color "$HI_STEM_COLOR" \
    --stem-color "$STEM_COLOR" \
    "${_HI_LABELS_ARG[@]}"
