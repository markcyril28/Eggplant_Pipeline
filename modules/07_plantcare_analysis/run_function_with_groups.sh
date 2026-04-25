#!/bin/bash
# Heatmap with Functions Organized by Major Groups

set -euo pipefail

# ========================================================================
# TOGGLES
# ========================================================================
OVERWRITE=true
OPTION="C"  # A: All, B: Filter, C: Custom list

# ========================================================================
# OPERATIONS - Toggle by commenting/uncommenting
# ========================================================================
OPERATIONS=(
    "create_matrix"
    #"create_heatmap"
    "create_heatmap_horizontal"
)

# ========================================================================
# PATHS
# ========================================================================
RAW_DIR=""
OUTPUT_DIR=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$SCRIPT_DIR"

# ========================================================================
# CONFIGURATION
# ========================================================================
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

# ========================================================================
# OPTION C: CUSTOM SELECTIONS
# ========================================================================

# List of Major Groups to include
SELECTED_MAJOR_GROUPS=(
    # "Core Promoter Elements"
    #"General CARE"
    "Phytohormone-Responsive Elements"
    "Growth and Development-Related Elements"
    "Stress-Responsive Elements"
    "Light-Responsive Elements"
)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-dir) RAW_DIR="$2"; shift 2 ;;
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
[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$BASE_DIR/08_Heatmaps_Function_Groups"

if [[ ${#GENES[@]} -eq 0 ]]; then
    # Auto-detect gene names from .tab files in RAW_DIR
    mapfile -t GENES < <(find "$RAW_DIR" -maxdepth 1 -name "*.tab" -printf "%f\n" 2>/dev/null | sed 's/\.tab$//' | sort)
    if [[ ${#GENES[@]} -eq 0 ]]; then
        echo "Error: No genes specified and no .tab files found in $RAW_DIR" >&2
        exit 1
    fi
    echo "Auto-detected ${#GENES[@]} genes from $(basename "$RAW_DIR")"
fi

# ========================================================================
# MAIN LOGIC
# ========================================================================

mkdir -p "$OUTPUT_DIR"

echo "=== Function with Major Groups Heatmap Pipeline ==="
echo "Option: $OPTION"

# Create matrix
if [[ " ${OPERATIONS[@]} " =~ " create_matrix " ]]; then
    echo "--- Creating Matrix ---"
    
    MATRIX_OUTPUT="$OUTPUT_DIR/function_groups_matrix.tsv"
    
    if [ "$OVERWRITE" = true ] || [ ! -f "$MATRIX_OUTPUT" ]; then
        if [ "$OPTION" == "C" ] && [ ${#SELECTED_MAJOR_GROUPS[@]} -gt 0 ]; then
            echo "Using selected major groups: ${#SELECTED_MAJOR_GROUPS[@]}"
            python3 "$MODULE_DIR/create_function_with_groups_matrix.py" \
                --raw_dir "$RAW_DIR" \
                --output_dir "$OUTPUT_DIR" \
                --genes "${GENES[@]}" \
                --option "C" \
                --groups "${SELECTED_MAJOR_GROUPS[@]}"
        else
            echo "Using all functions and groups"
            python3 "$MODULE_DIR/create_function_with_groups_matrix.py" \
                --raw_dir "$RAW_DIR" \
                --output_dir "$OUTPUT_DIR" \
                --genes "${GENES[@]}" \
                --option "$OPTION"
        fi
        echo "Matrix created"
    else
        echo "Matrix exists, skipping"
    fi
fi

# Create heatmap
if [[ " ${OPERATIONS[@]} " =~ " create_heatmap " ]]; then
    echo "--- Creating Heatmap ---"
    
    MATRIX_FILE="$OUTPUT_DIR/function_groups_matrix.tsv"
    HEATMAP_OUTPUT="$OUTPUT_DIR/function_groups_heatmap.png"
    
    if [ ! -f "$MATRIX_FILE" ]; then
        echo "Error: Matrix not found at $MATRIX_FILE"
        exit 1
    fi
    
    if [ "$OVERWRITE" = true ] || [ ! -f "$HEATMAP_OUTPUT" ]; then
        Rscript "$MODULE_DIR/create_function_groups_heatmap.R" \
            "$MATRIX_FILE" \
            "$OUTPUT_DIR" \
            "$DPI" \
            "$COLOR_PALETTE" \
            "$ROW_FONT" "$COL_FONT" "$CELL_FONT" "$CELL_SIZE" \
            "$LABEL_FONT" "$CELL_BORDER_COLOR" "$CELL_BORDER_WIDTH" \
            "$COLUMN_ROTATION" "$LEGEND_HEIGHT" "$GENE_LABEL"
        if [ -f "$HEATMAP_OUTPUT" ]; then
            echo "Heatmap created"
        else
            echo "Heatmap skipped (empty or missing metadata/matrix)"
        fi
    else
        echo "Heatmap exists, skipping"
    fi
fi

# Create heatmap with horizontal labels
if [[ " ${OPERATIONS[@]} " =~ " create_heatmap_horizontal " ]]; then
    echo "--- Creating Heatmap with Horizontal Labels ---"
    
    MATRIX_FILE="$OUTPUT_DIR/function_groups_matrix.tsv"
    HEATMAP_OUTPUT="$OUTPUT_DIR/function_groups_heatmap_horizontal.png"
    
    if [ ! -f "$MATRIX_FILE" ]; then
        echo "Error: Matrix not found at $MATRIX_FILE"
        exit 1
    fi
    
    if [ "$OVERWRITE" = true ] || [ ! -f "$HEATMAP_OUTPUT" ]; then
        Rscript "$MODULE_DIR/create_function_groups_heatmap_horizontal.R" \
            "$MATRIX_FILE" \
            "$OUTPUT_DIR" \
            "$DPI" \
            "$COLOR_PALETTE" \
            "$ROW_FONT" "$COL_FONT" "$CELL_FONT" "$CELL_SIZE" \
            "$LABEL_FONT" "$CELL_BORDER_COLOR" "$CELL_BORDER_WIDTH" \
            "$COLUMN_ROTATION" "$LEGEND_HEIGHT" "$GENE_LABEL"
        if [ -f "$HEATMAP_OUTPUT" ]; then
            echo "Heatmap with horizontal labels created"
        else
            echo "Heatmap with horizontal labels skipped (empty or missing metadata/matrix)"
        fi
    else
        echo "Heatmap exists, skipping"
    fi
fi

echo "=== Pipeline Complete ==="
echo "Output: $OUTPUT_DIR"
