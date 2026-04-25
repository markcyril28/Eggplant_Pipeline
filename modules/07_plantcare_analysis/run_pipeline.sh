#!/bin/bash
# Module: PlantCARE Promoter Analysis Pipeline
# Usage: bash run_pipeline.sh --raw-dir <dir> --outdir <dir> --module-dir <dir> \
#        [--steps post-process,matrix,heatmap] [--overwrite] \
#        [--genes "Gene1,Gene2"] [--functions "func1,func2"]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../logging/logging_utils.sh"

# Defaults
RAW_DIR=""
OUTPUT_DIR="."
MODULE_DIR="$SCRIPT_DIR"   # default: scripts live alongside this file
STEPS="post-process,matrix,heatmap"
OVERWRITE=false
SELECTED_GENES=""
SELECTED_FUNCTIONS=""
DPI=900
ROW_FONT=12
COL_FONT=10
LABEL_FONT=14
CELL_FONT=10
CELL_SIZE=1
MIN_FREQ=0
COLUMN_ROTATION=22.5
CELL_BORDER_COLOR="white"
CELL_BORDER_WIDTH=0.5
LEGEND_HEIGHT=4
GENE_LABEL="Genes of Interest"
COLOR_PALETTE=""
CPU=$(nproc)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --raw-dir)          RAW_DIR="$2"; shift 2 ;;
        --outdir)           OUTPUT_DIR="$2"; shift 2 ;;
        --module-dir)       MODULE_DIR="$2"; shift 2 ;;
        --steps)            STEPS="$2"; shift 2 ;;
        --overwrite)        OVERWRITE=true; shift ;;
        --genes)            SELECTED_GENES="$2"; shift 2 ;;
        --functions)        SELECTED_FUNCTIONS="$2"; shift 2 ;;
        --dpi)              DPI="$2"; shift 2 ;;
        --threads)          CPU="$2"; shift 2 ;;
        --row-font)         ROW_FONT="$2"; shift 2 ;;
        --col-font)         COL_FONT="$2"; shift 2 ;;
        --label-font)       LABEL_FONT="$2"; shift 2 ;;
        --cell-font)        CELL_FONT="$2"; shift 2 ;;
        --cell-size)        CELL_SIZE="$2"; shift 2 ;;
        --min-freq)         MIN_FREQ="$2"; shift 2 ;;
        --column-rotation)  COLUMN_ROTATION="$2"; shift 2 ;;
        --cell-border-color) CELL_BORDER_COLOR="$2"; shift 2 ;;
        --cell-border-width) CELL_BORDER_WIDTH="$2"; shift 2 ;;
        --legend-height)    LEGEND_HEIGHT="$2"; shift 2 ;;
        --gene-label)       GENE_LABEL="$2"; shift 2 ;;
        --color-palette)    COLOR_PALETTE="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$RAW_DIR" ]]   && { log_error "Missing --raw-dir (PlantCARE raw results)"; exit 1; }
[[ ! -d "$RAW_DIR" ]] && { log_error "Raw PlantCARE directory does not exist: $RAW_DIR"; exit 1; }

POST_PROCESSED_DIR="$OUTPUT_DIR/04_Post_Processed"
MATRIX_DIR="$OUTPUT_DIR/05_Matrices"
HEATMAP_DIR="$OUTPUT_DIR/06_Heatmaps"

mkdir -p "$POST_PROCESSED_DIR" "$MATRIX_DIR" "$HEATMAP_DIR"

should_run() { [[ ",${STEPS}," == *",${1},"* ]]; }

# Post-processing
if should_run "post-process"; then
    log_step "PlantCARE Post-processing"
    if [[ "$OVERWRITE" == true ]] || [[ -z "$(ls -A "$POST_PROCESSED_DIR" 2>/dev/null)" ]]; then
        python3 "$MODULE_DIR/post_process_plantcare.py" -i "$RAW_DIR" -o "$POST_PROCESSED_DIR"
        log_info "Post-processing complete"
    else
        log_info "Post-processing skipped (output exists)"
    fi
fi

# Matrix generation
if should_run "matrix"; then
    log_step "Matrix Generation"
    if [[ "$OVERWRITE" == true ]] || [[ -z "$(ls -A "$MATRIX_DIR" 2>/dev/null)" ]]; then
        python3 "$MODULE_DIR/plantCARE_to_matrix.py" -i "$POST_PROCESSED_DIR" -o "$MATRIX_DIR"
        log_info "Matrix generation complete"
    else
        log_info "Matrix generation skipped (output exists)"
    fi
fi

# --- Heatmap generation (3 independent steps run in parallel) ---
HEATMAP_PIDS=()

# Heatmap
if should_run "heatmap"; then
    (
        log_step "Heatmap Generation"

        # Helper: generate a single heatmap from a matrix file
        generate_heatmap() {
            local matrix_file="$1" output_file="$2" col_label="$3" filter_tag="$4"
            if [[ ! -f "$matrix_file" ]]; then
                log_warn "Matrix not found: $matrix_file"
                return 0
            fi

            local matrix_to_use="$matrix_file"

            # Apply filtering if requested
            if [[ -n "$SELECTED_GENES" || -n "$SELECTED_FUNCTIONS" ]]; then
                local filtered="$MATRIX_DIR/filtered_${filter_tag}.tsv"
                local -a filter_args=( python3 "$MODULE_DIR/filter_matrix.py" -i "$matrix_file" -o "$filtered" )
                [[ -n "$SELECTED_GENES" ]]     && filter_args+=( --genes "$SELECTED_GENES" )
                [[ -n "$SELECTED_FUNCTIONS" ]] && filter_args+=( --columns "$SELECTED_FUNCTIONS" )
                "${filter_args[@]}"
                matrix_to_use="$filtered"
            fi

            if [[ "$OVERWRITE" == true ]] || [[ ! -f "$output_file" ]]; then
                local -a r_cmd=(
                    Rscript "$MODULE_DIR/visualize_plantCARE_matrix.R"
                    -i "$matrix_to_use" -o "$output_file"
                    --gene_label "$GENE_LABEL"
                    --column_label "$col_label"
                    --cell_size "$CELL_SIZE" --min_freq "$MIN_FREQ"
                    --dpi "$DPI" --row_font "$ROW_FONT" --col_font "$COL_FONT"
                    --label_font "$LABEL_FONT" --cell_font_size "$CELL_FONT"
                    --column_rotation "$COLUMN_ROTATION"
                    --cell_border_color "$CELL_BORDER_COLOR"
                    --cell_border_width "$CELL_BORDER_WIDTH"
                    --legend_height "$LEGEND_HEIGHT"
                )
                [[ -n "$COLOR_PALETTE" ]] && r_cmd+=(--color_palette "$COLOR_PALETTE")
                "${r_cmd[@]}"
                log_info "Heatmap -> $output_file"
            else
                log_info "Heatmap skipped (exists): $output_file"
            fi
        }

        # Version 1: Motif-based heatmap
        generate_heatmap \
            "$MATRIX_DIR/plantcare_heatmap_matrix_v1_motif.tsv" \
            "$HEATMAP_DIR/plantcare_heatmap_v1_motif.png" \
            "Motifs" "v1_motif"

        # Version 2: Function-based heatmap
        generate_heatmap \
            "$MATRIX_DIR/plantcare_heatmap_matrix_v2_function.tsv" \
            "$HEATMAP_DIR/plantcare_heatmap_v2_function.png" \
            "Functions" "v2_function"

        # TBTools: Motif-based heatmap (TBTools-compatible format)
        generate_heatmap \
            "$MATRIX_DIR/plantcare_tbtools_matrix.tsv" \
            "$HEATMAP_DIR/plantcare_heatmap_tbtools_motif.png" \
            "Motifs" "tbtools_motif"
    ) &
    HEATMAP_PIDS+=($!)
fi

# Combined motif-function heatmap
if should_run "combined-heatmap"; then
    (
        log_step "Combined Motif-Function Heatmap"
        COMBINED_DIR="$OUTPUT_DIR/07_Heatmaps_Combined_Motif_Function"
        mkdir -p "$COMBINED_DIR"
        COMBINED_SCRIPT="$SCRIPT_DIR/c_run_function_with_motif_heatmap.sh"
        if [[ -f "$COMBINED_SCRIPT" ]]; then
            combined_cmd=(
                bash "$COMBINED_SCRIPT"
                --raw-dir "$RAW_DIR"
                --matrix-dir "$MATRIX_DIR"
                --output-dir "$COMBINED_DIR"
                --overwrite "$OVERWRITE"
                --dpi "$DPI"
                --color-palette "$COLOR_PALETTE"
                --row-font "$ROW_FONT"
                --col-font "$COL_FONT"
                --cell-font "$CELL_FONT"
                --cell-size "$CELL_SIZE"
                --label-font "$LABEL_FONT"
                --cell-border-color "$CELL_BORDER_COLOR"
                --cell-border-width "$CELL_BORDER_WIDTH"
                --column-rotation "$COLUMN_ROTATION"
                --legend-height "$LEGEND_HEIGHT"
                --gene-label "$GENE_LABEL"
            )
            [[ -n "$SELECTED_GENES" ]] && combined_cmd+=(--genes "$SELECTED_GENES")
            "${combined_cmd[@]}"
            log_info "Combined motif-function heatmap complete"
        else
            log_warn "Combined heatmap script not found: $COMBINED_SCRIPT"
        fi
    ) &
    HEATMAP_PIDS+=($!)
fi

# Function with major groups heatmap
if should_run "groups-heatmap"; then
    (
        log_step "Function Groups Heatmap"
        GROUPS_DIR="$OUTPUT_DIR/08_Heatmaps_Function_Groups"
        mkdir -p "$GROUPS_DIR"
        GROUPS_SCRIPT="$SCRIPT_DIR/run_function_with_groups.sh"
        if [[ -f "$GROUPS_SCRIPT" ]]; then
            groups_cmd=(
                bash "$GROUPS_SCRIPT"
                --raw-dir "$RAW_DIR"
                --output-dir "$GROUPS_DIR"
                --overwrite "$OVERWRITE"
                --dpi "$DPI"
                --color-palette "$COLOR_PALETTE"
                --row-font "$ROW_FONT"
                --col-font "$COL_FONT"
                --cell-font "$CELL_FONT"
                --cell-size "$CELL_SIZE"
                --label-font "$LABEL_FONT"
                --cell-border-color "$CELL_BORDER_COLOR"
                --cell-border-width "$CELL_BORDER_WIDTH"
                --column-rotation "$COLUMN_ROTATION"
                --legend-height "$LEGEND_HEIGHT"
                --gene-label "$GENE_LABEL"
            )
            [[ -n "$SELECTED_GENES" ]] && groups_cmd+=(--genes "$SELECTED_GENES")
            "${groups_cmd[@]}"
            log_info "Function groups heatmap complete"
        else
            log_warn "Function groups script not found: $GROUPS_SCRIPT"
        fi
    ) &
    HEATMAP_PIDS+=($!)
fi

# Wait for all heatmap jobs
FAILED=0
for pid in "${HEATMAP_PIDS[@]+${HEATMAP_PIDS[@]}}"; do
    wait "$pid" || ((FAILED++))
done
(( FAILED > 0 )) && { log_error "$FAILED heatmap job(s) failed"; exit 1; }

log_step "PlantCARE pipeline complete"
teardown_logging
