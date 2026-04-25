#!/bin/bash

# Combined Motif and Function Heatmap Pipeline
set -euo pipefail

# Toggles
OVERWRITE=true
OPTION="C"  # A: All, B: Filter, C: Custom list

# Paths
RAW_DIR=""
MATRIX_DIR=""
OUTPUT_DIR=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$SCRIPT_DIR"

# Config
GENES=()
DPI=900
COLOR_PALETTE=""
ROW_FONT=10
COL_FONT=8
CELL_FONT=8
CELL_SIZE=1
LABEL_FONT=14
CELL_BORDER_COLOR="white"
CELL_BORDER_WIDTH=0.5
COLUMN_ROTATION=22.5
LEGEND_HEIGHT=4
GENE_LABEL="Genes of Interest"

# Option C: Custom function list 
SELECTED_FUNCTIONS_V2=(
    #"core promoter element around -30 of transcription start"
    #"common CAE in promoter and enhancer regions"
    "CARE involved in the MeJA-responsiveness"
    #"MYB binding site involved in drought-inducibility"
    #"MYB binding site involved in light responsiveness"
    #"gibberellin-responsive element"
    #"CAE involved in gibberellin-responsiveness"
    "CARE related to meristem specific activation"
    "CARE related to meristem expression"
    "CARE involved in auxin responsiveness"
    #"part of an auxin-responsive element"
    # "60K protein binding site"
    # "MYBHv1 binding site"
    # "auxin-responsive element"
    # "binding site of AT-rich DNA binding protein (ATBP-1)"
    # "CAE involved in defense and stress responsiveness"
    # "CAE involved in light responsiveness"
    # "CAE involved in low-temperature responsiveness"
    # "CAE involved in salicylic acid responsiveness"
    # "CAE involved in the abscisic acid responsiveness"
    "CARE"
    "CARE essential for the anaerobic induction"
    "CARE involved in circadian control"
    "CARE involved in light responsiveness"
    "CARE involved in zein metabolism regulation"
    # "cis-regulatory element involved in endosperm expression"
    # "element for maximal elicitor-mediated activation (2copies)"
    # "light responsive element"
    # "part of a conserved DNA module involved in phytohormone and stress responsive expression (CMA3)"
    # "part of a conserved DNA module involved in light responsiveness"
    #"part of a light responsive element"
    # "part of a light responsive element"
    # "part of a light responsive module"
    # "part of a module for light response"
    # "part of gapA in (gapA-CMA1) involved with light responsiveness"
    # "protein binding site"
)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-dir) RAW_DIR="$2"; shift 2 ;;
        --matrix-dir) MATRIX_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --overwrite)
            if [[ $# -ge 2 && "$2" != --* ]]; then
                OVERWRITE="$2"; shift 2
            else
                OVERWRITE=true; shift
            fi
            ;;
        --overwrite=*)
            OVERWRITE="${1#*=}"
            shift
            ;;
        --genes)
            IFS=',' read -ra GENES <<< "$2"; shift 2
            ;;
        --dpi)              DPI="$2"; shift 2 ;;
        --color-palette)    COLOR_PALETTE="$2"; shift 2 ;;
        --row-font)         ROW_FONT="$2"; shift 2 ;;
        --col-font)         COL_FONT="$2"; shift 2 ;;
        --cell-font)        CELL_FONT="$2"; shift 2 ;;
        --cell-size)        CELL_SIZE="$2"; shift 2 ;;
        --label-font)       LABEL_FONT="$2"; shift 2 ;;
        --cell-border-color) CELL_BORDER_COLOR="$2"; shift 2 ;;
        --cell-border-width) CELL_BORDER_WIDTH="$2"; shift 2 ;;
        --column-rotation)  COLUMN_ROTATION="$2"; shift 2 ;;
        --legend-height)    LEGEND_HEIGHT="$2"; shift 2 ;;
        --gene-label)       GENE_LABEL="$2"; shift 2 ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

BASE_DIR="$(pwd)"
[[ -z "$RAW_DIR" ]] && RAW_DIR="$BASE_DIR/../PlantCARE_Results"
[[ -z "$MATRIX_DIR" ]] && MATRIX_DIR="$BASE_DIR/05_Matrices"
[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$BASE_DIR/07_Heatmaps_Combined_Motif_Function"

if [[ ${#GENES[@]} -eq 0 ]]; then
    # Auto-detect gene names from .tab files in RAW_DIR
    mapfile -t GENES < <(find "$RAW_DIR" -maxdepth 1 -name "*.tab" -printf "%f\n" 2>/dev/null | sed 's/\.tab$//' | sort)
    if [[ ${#GENES[@]} -eq 0 ]]; then
        echo "Error: No genes specified and no .tab files found in $RAW_DIR" >&2
        exit 1
    fi
    echo "Auto-detected ${#GENES[@]} genes from $(basename "$RAW_DIR")"
fi

mkdir -p "$OUTPUT_DIR"

echo "=== Combined Motif-Function Heatmap Pipeline ==="
echo "Option: $OPTION"
echo "Number of selected functions: ${#SELECTED_FUNCTIONS_V2[@]}"

# Generate matrix
if [ "$OPTION" == "C" ] && [ ${#SELECTED_FUNCTIONS_V2[@]} -gt 0 ]; then
    echo "Using custom function list..."
    python3 "$MODULE_DIR/create_combined_motif_function_matrix.py" \
        --raw_dir "$RAW_DIR" \
        --function_matrix "$MATRIX_DIR/plantcare_heatmap_matrix_v2_function.tsv" \
        --output_dir "$OUTPUT_DIR" \
        --genes "${GENES[@]}" \
        --option "C" \
        --functions "${SELECTED_FUNCTIONS_V2[@]}"
else
    echo "Using option $OPTION (all functions)..."
    python3 "$MODULE_DIR/create_combined_motif_function_matrix.py" \
        --raw_dir "$RAW_DIR" \
        --function_matrix "$MATRIX_DIR/plantcare_heatmap_matrix_v2_function.tsv" \
        --output_dir "$OUTPUT_DIR" \
        --genes "${GENES[@]}" \
        --option "$OPTION"
fi

# Generate heatmap
Rscript "$MODULE_DIR/create_combined_heatmap.R" \
    "$OUTPUT_DIR/combined_matrix.tsv" \
    "$OUTPUT_DIR" \
    "$DPI" \
    "$COLOR_PALETTE" \
    "$ROW_FONT" "$COL_FONT" "$CELL_FONT" "$CELL_SIZE" \
    "$LABEL_FONT" "$CELL_BORDER_COLOR" "$CELL_BORDER_WIDTH" \
    "$COLUMN_ROTATION" "$LEGEND_HEIGHT" "$GENE_LABEL"

echo "=== Pipeline Complete ==="
echo "Output: $OUTPUT_DIR"
