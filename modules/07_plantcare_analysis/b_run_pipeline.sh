#!/bin/bash

# Main pipeline script for PlantCARE analysis
# Supports Version 1 (Motif Name) and Version 2 (Function) heatmaps

set -euo pipefail

# Default toggles
RUN_ALL=false
RUN_POST_PROCESSING=false
RUN_MATRIX_GENERATION=false
RUN_HEATMAP_V1=false
RUN_HEATMAP_V2=true
OVERWRITE=true

# ========================================================================
# FILTERING OPTIONS - Toggle by commenting/uncommenting
# ========================================================================

# Option A: All Genes and All Columns (Default - no filtering)
#FILTER_GENES=false
#FILTER_COLUMNS=false

# Option B: Filter specific genes and/or columns
# Uncomment lines below to enable filtering

#FILTER_GENES=true
FILTER_COLUMNS=true

# --- Gene Filtering ---
# List genes to include in heatmap (comment out unwanted genes)
SELECTED_GENES=(
    "SmelGRF05_970"
    "SmelGRF08_140"
    "SmelGIF11_070"
    "SmelGIF11_650"
    "SmelGIF11_790"
)

# --- Column Filtering for Version 1 (Motif Names) ---
# Uncomment specific motifs to include (only these will be shown)
SELECTED_MOTIFS_V1=(
    # "CAAT-box"
    # "TATA-box"
    # Add more motifs as needed
)

# --- Column Filtering for Version 2 (Functions) ---
# Uncomment specific functions to include (only these will be shown)
SELECTED_FUNCTIONS_V2=(
    # "common CAE, in promoter and enhancer regions"
    # "part of a conserved DNA module light responsiveness"
    # "60K protein binding site"
    "MYB binding site drought-inducibility"
    "MYB binding site light responsiveness"
    # "MYBHv1 binding site"
    # "auxin-responsive element"
    # "binding site of AT-rich DNA binding protein (ATBP-1)"
    # "CAE, defense and stress responsiveness"
    # "CAE, gibberellin-responsiveness"
    # "CAE, light responsiveness"
    # "CAE, low-temperature responsiveness"
    # "CAE, salicylic acid responsiveness"
    # "CAE, the abscisic acid responsiveness"
    # "CARE"
    # "CARE essential for the anaerobic induction"
    "CARE auxin responsiveness"
    # "CARE circadian control"
    # "CARE light responsiveness"
    "CARE the MeJA-responsiveness"
    # "CARE zein metabolism regulation"
    "CARE related to meristem expression"
    "CARE related to meristem specific activation"
    # "cis-regulatory element endosperm expression"
    "common CAE in promoter and enhancer regions"
    "core promoter element -30 of transcription start"
    # "element for maximal elicitor-mediated activation (2copies)"
    "gibberellin-responsive element"
    # "light responsive element"
    # "part of a conserved DNA module array (CMA3)"
    # "part of a conserved DNA module light responsiveness"
    "part of a light response element"
    # "part of a light responsive element"
    # "part of a light responsive module"
    # "part of a module for light response"
    # "part of an auxin-responsive element"
    # "part of gapA in (gapA-CMA1) involved with light responsiveness"
    # "protein binding site"
)

# --- Configuration ---
# Heatmap settings
GENE_LABEL="Genes of Interest"
MOTIF_LABEL_V1="CAREs (CARE)"
FUNCTION_LABEL_V2="Functions"
CELL_SIZE=1
MIN_FREQ=0
DPI=900
ROW_FONT=12
COL_FONT=10
LABEL_FONT=14
CELL_FONT_SIZE=10

BASE_DIR=$(pwd)

# Directories
RAW_DIR="$BASE_DIR/03_PlantCARE_Results"
POST_PROCESSED_DIR="$BASE_DIR/04_Post_Processed"
MATRIX_DIR="$BASE_DIR/05_Matrices"
HEATMAP_DIR="$BASE_DIR/06_Heatmaps"
TBTOOLS_HEATMAP_DIR="$BASE_DIR/06_Heatmaps_TBTools"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$SCRIPT_DIR"

# Scripts
POST_PROCESS_SCRIPT="$MODULE_DIR/post_process_plantcare.py"
MATRIX_SCRIPT="$MODULE_DIR/plantCARE_to_matrix.py"
HEATMAP_SCRIPT="$MODULE_DIR/visualize_plantCARE_matrix.R"

# ========================================================================

# --- Functions ---

# Function to run post-processing
run_post_processing() {
    echo "--- Running Post-processing ---"
    if [ "$OVERWRITE" = true ] || [ ! -d "$POST_PROCESSED_DIR" ] || [ -z "$(ls -A "$POST_PROCESSED_DIR" 2>/dev/null)" ]; then
        python3 "$POST_PROCESS_SCRIPT" -i "$RAW_DIR" -o "$POST_PROCESSED_DIR"
        echo "Post-processing complete"
    else
        echo "Post-processing skipped (output exists)"
    fi
}

# Function to run matrix generation
run_matrix_generation() {
    echo "--- Running Matrix Generation ---"
    if [ "$OVERWRITE" = true ] || [ ! -d "$MATRIX_DIR" ] || [ -z "$(ls -A "$MATRIX_DIR" 2>/dev/null)" ]; then
        python3 "$MATRIX_SCRIPT" -i "$POST_PROCESSED_DIR" -o "$MATRIX_DIR"
        echo "Matrix generation complete"
    else
        echo "Matrix generation skipped (output exists)"
    fi
}

# Function to run heatmap visualization Version 1 (Motif Names)
run_heatmap_v1() {
    echo "--- Running Heatmap Visualization V1 (Motif Names) ---"
    MATRIX_V1="$MATRIX_DIR/plantcare_heatmap_matrix_v1_motif.tsv"
    OUTPUT_V1="$HEATMAP_DIR/plantcare_heatmap_v1_motif.png"
    
    if [ ! -f "$MATRIX_V1" ]; then
        echo "Error: Matrix V1 not found at $MATRIX_V1"
        return 1
    fi
    
    # Apply filtering if enabled
    MATRIX_TO_USE="$MATRIX_V1"
    if [ "$FILTER_GENES" = true ] || [ "$FILTER_COLUMNS" = true ]; then
        echo "Applying filters..."
        FILTERED_MATRIX="$MATRIX_DIR/filtered_v1_motif.tsv"

        local -a FILTER_ARGS=( python3 "$MODULE_DIR/filter_matrix.py" -i "$MATRIX_V1" -o "$FILTERED_MATRIX" )

        if [ "$FILTER_GENES" = true ] && [ ${#SELECTED_GENES[@]} -gt 0 ]; then
            GENE_LIST=$(IFS=,; echo "${SELECTED_GENES[*]}")
            FILTER_ARGS+=( --genes "$GENE_LIST" )
            echo "Filtering genes: $GENE_LIST"
        fi

        if [ "$FILTER_COLUMNS" = true ] && [ ${#SELECTED_MOTIFS_V1[@]} -gt 0 ]; then
            COL_LIST=$(IFS=,; echo "${SELECTED_MOTIFS_V1[*]}")
            FILTER_ARGS+=( --columns "$COL_LIST" )
            echo "Filtering motifs: $COL_LIST"
        fi

        "${FILTER_ARGS[@]}"
        MATRIX_TO_USE="$FILTERED_MATRIX"
    fi
    
    if [ "$OVERWRITE" = true ] || [ ! -f "$OUTPUT_V1" ]; then
        Rscript "$HEATMAP_SCRIPT" \
            -i "$MATRIX_TO_USE" \
            -o "$OUTPUT_V1" \
            --gene_label "$GENE_LABEL" \
            --column_label "$MOTIF_LABEL_V1" \
            --cell_size "$CELL_SIZE" \
            --min_freq "$MIN_FREQ" \
            --dpi "$DPI" \
            --row_font "$ROW_FONT" \
            --col_font "$COL_FONT" \
            --label_font "$LABEL_FONT" \
            --cell_font_size "$CELL_FONT_SIZE"
        echo "Heatmap V1 complete"
    else
        echo "Heatmap V1 skipped (output exists)"
    fi
}

# Function to run heatmap visualization Version 2 (Functions)
run_heatmap_v2() {
    echo "--- Running Heatmap Visualization V2 (Functions) ---"
    MATRIX_V2="$MATRIX_DIR/plantcare_heatmap_matrix_v2_function.tsv"
    OUTPUT_V2="$HEATMAP_DIR/plantcare_heatmap_v2_function.png"
    
    if [ ! -f "$MATRIX_V2" ]; then
        echo "Warning: Matrix V2 not found at $MATRIX_V2. Skipping."
        return 0
    fi
    
    # Apply filtering if enabled
    MATRIX_TO_USE="$MATRIX_V2"
    if [ "$FILTER_GENES" = true ] || [ "$FILTER_COLUMNS" = true ]; then
        echo "Applying filters..."
        FILTERED_MATRIX="$MATRIX_DIR/filtered_v2_function.tsv"

        local -a FILTER_ARGS=( python3 "$MODULE_DIR/filter_matrix.py" -i "$MATRIX_V2" -o "$FILTERED_MATRIX" )

        if [ "$FILTER_GENES" = true ] && [ ${#SELECTED_GENES[@]} -gt 0 ]; then
            GENE_LIST=$(IFS=,; echo "${SELECTED_GENES[*]}")
            FILTER_ARGS+=( --genes "$GENE_LIST" )
            echo "Filtering genes: $GENE_LIST"
        fi

        if [ "$FILTER_COLUMNS" = true ] && [ ${#SELECTED_FUNCTIONS_V2[@]} -gt 0 ]; then
            COL_LIST=$(IFS=,; echo "${SELECTED_FUNCTIONS_V2[*]}")
            FILTER_ARGS+=( --columns "$COL_LIST" )
            echo "Filtering functions: $COL_LIST"
        fi

        "${FILTER_ARGS[@]}"
        MATRIX_TO_USE="$FILTERED_MATRIX"
    fi
    
    if [ "$OVERWRITE" = true ] || [ ! -f "$OUTPUT_V2" ]; then
        Rscript "$HEATMAP_SCRIPT" \
            -i "$MATRIX_TO_USE" \
            -o "$OUTPUT_V2" \
            --gene_label "$GENE_LABEL" \
            --column_label "$FUNCTION_LABEL_V2" \
            --cell_size "$CELL_SIZE" \
            --min_freq "$MIN_FREQ" \
            --dpi "$DPI" \
            --row_font "$ROW_FONT" \
            --col_font "$COL_FONT" \
            --label_font "$LABEL_FONT" \
            --cell_font_size "$CELL_FONT_SIZE"
        echo "Heatmap V2 complete"
    else
        echo "Heatmap V2 skipped (output exists)"
    fi
}

# --- Main Logic ---

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        all) RUN_ALL=true; shift ;;
        post-process) RUN_POST_PROCESSING=true; shift ;;
        matrix) RUN_MATRIX_GENERATION=true; shift ;;
        heatmap-v1) RUN_HEATMAP_V1=true; shift ;;
        heatmap-v2) RUN_HEATMAP_V2=true; shift ;;
        heatmap) RUN_HEATMAP_V1=true; RUN_HEATMAP_V2=true; shift ;;
        --overwrite) OVERWRITE=true; shift ;;
        -h|--help) 
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  all              Run all steps"
            echo "  post-process     Run post-processing only"
            echo "  matrix           Run matrix generation only"
            echo "  heatmap-v1       Run heatmap V1 (Motif Names) only"
            echo "  heatmap-v2       Run heatmap V2 (Functions) only"
            echo "  heatmap          Run both heatmap versions"
            echo "  --overwrite      Overwrite existing outputs"
            echo "  -h, --help       Show this help message"
            echo ""
            echo "Filtering Options:"
            echo "  Edit the script to enable filtering by setting:"
            echo "    FILTER_GENES=true       # Filter specific genes"
            echo "    FILTER_COLUMNS=true     # Filter specific motifs/functions"
            echo "  Then specify genes/motifs/functions in:"
            echo "    SELECTED_GENES array"
            echo "    SELECTED_MOTIFS_V1 array (for Version 1)"
            echo "    SELECTED_FUNCTIONS_V2 array (for Version 2)"
            exit 0
            ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

# If no specific step is selected, run all
if [ "$RUN_ALL" = false ] && [ "$RUN_POST_PROCESSING" = false ] && \
   [ "$RUN_MATRIX_GENERATION" = false ] && [ "$RUN_HEATMAP_V1" = false ] && \
   [ "$RUN_HEATMAP_V2" = false ]; then
    RUN_ALL=true
fi

# Create directories
mkdir -p "$POST_PROCESSED_DIR"
mkdir -p "$MATRIX_DIR"
mkdir -p "$HEATMAP_DIR"
mkdir -p "$TBTOOLS_HEATMAP_DIR"

# Execute pipeline steps
if [ "$RUN_ALL" = true ] || [ "$RUN_POST_PROCESSING" = true ]; then
    run_post_processing
fi

if [ "$RUN_ALL" = true ] || [ "$RUN_MATRIX_GENERATION" = true ]; then
    run_matrix_generation
fi

if [ "$RUN_ALL" = true ] || [ "$RUN_HEATMAP_V1" = true ]; then
    run_heatmap_v1
fi

if [ "$RUN_ALL" = true ] || [ "$RUN_HEATMAP_V2" = true ]; then
    run_heatmap_v2
fi

echo "=== Pipeline Finished ==="
