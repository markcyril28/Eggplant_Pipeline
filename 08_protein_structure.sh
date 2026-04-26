#!/bin/bash
# ============================================================================
# Program 8: Protein Structure Analysis
# ============================================================================
# Edit [protein_structure].gene_groups in 08_protein_structureCONFIG.toml, then run:
#   bash h_protein_structure.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# All settings including gene_groups are loaded from 08_protein_structureCONFIG.toml.
# Edit [protein_structure].gene_groups there to select which groups to process.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# Load GENE_GROUPS from shared config (before the per-group loop)
SHARED_CONFIG="$PIPELINE_DIR/08_protein_structureCONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" protein_structure gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: protein_structure.gene_groups is empty in 08_protein_structureCONFIG.toml" >&2
    exit 1
fi

TEMP_FILES=()
cleanup_all() {
    rm -f "${TEMP_FILES[@]+"${TEMP_FILES[@]}"}" 2>/dev/null || true
    safe_teardown_logging
}
trap cleanup_all EXIT

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# Resolve config: split directory or monolithic file
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_protstruct_cfg_XXXXXX.toml")
    TEMP_FILES+=("$CONFIG_FILE")
    python3 "$MODULES/utils/merge_toml.py" \
        "$PIPELINE_DIR/08_protein_structureCONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/h_protein_structure_analysis.toml" > "$CONFIG_FILE"
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

MACHINE=$(get_toml protein_structure machine 2>/dev/null || echo "Local")
CPU=$(get_toml protein_structure compute "$MACHINE" threads 2>/dev/null || nproc)
MAX_PARALLEL=$(get_toml protein_structure compute "$MACHINE" max_parallel 2>/dev/null || echo "2")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
OVERWRITE=$(get_toml protein_structure overwrite 2>/dev/null || echo "true")

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"

PROT_DIR="$BASE_DIR/08_Protein_Structure"
mkdir -p "$PROT_DIR"

setup_logging
log_step "Protein Structure Analysis: $GENE_GROUP"

NUC_FASTA_REL=$(get_toml protein_structure input_nucleotide_fasta 2>/dev/null || true)
if [[ -n "$NUC_FASTA_REL" ]]; then
    NUC_FASTA="$PIPELINE_DIR/$NUC_FASTA_REL"
else
    NUC_FASTA=""
fi
NUC_BASE=$(basename "${NUC_FASTA:-unnamed}" | sed 's/\.\(fa\|fasta\|fna\)$//')
PROT_OUTPUT="$PROT_DIR/${NUC_BASE}_polypeptide.fa"

if [[ -n "$NUC_FASTA" && -f "$NUC_FASTA" ]]; then
    bash "$MODULES/08_protein_structure_analysis/translate.sh" \
        --input "$NUC_FASTA" \
        --output "$PROT_OUTPUT" \
        --threads "$CPU"
elif [[ -z "$NUC_FASTA" ]]; then
    log_warn "No input_nucleotide_fasta configured for $GENE_GROUP — skipping translation"
else
    log_warn "Nucleotide FASTA not found: $NUC_FASTA"
    log_info "Place your nucleotide FASTA at the path above, then rerun."
fi

# ── Gene-specific structure processing (Extract → CIF→PDB → HEADER → Render)
COLOR_CONFIG="$PIPELINE_DIR/$(get_toml protein_structure color_config 2>/dev/null || echo "config/colors_config/protein_structure_colors.toml")"
STRUCT_PROCESSOR="$MODULES/08_protein_structure_analysis/gene_specific_structure_processer/process_structures.sh"

# Find run directories: direct children of $PROT_DIR that contain
# zip files OR AlphaFold3_Results/, Protein_Structures/, or Protein_Structure_* folders.
RUN_DIRS=()
for d in "$PROT_DIR"/*/; do
    d="${d%/}"
    [[ -d "$d" ]] || continue
    if compgen -G "$d/*.zip" >/dev/null 2>&1 || \
       find "$d" -maxdepth 1 -type d \( -name "AlphaFold3_Results" -o -name "Protein_Structures" -o -name "Protein_Structure_*" \) -print -quit 2>/dev/null | grep -q .; then
        RUN_DIRS+=("$d")
    fi
done

if (( ${#RUN_DIRS[@]} > 0 )); then
    log_step "Processing ${#RUN_DIRS[@]} run directory(ies)"
    # Divide threads among concurrent run dirs
    THREADS_PER_RUN=$(( CPU / (${#RUN_DIRS[@]} < MAX_PARALLEL ? ${#RUN_DIRS[@]} : MAX_PARALLEL) ))
    (( THREADS_PER_RUN < 1 )) && THREADS_PER_RUN=1
    wait_for_slot() { while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done; }

    # Parse operation lists from TOML → derive DO_* env vars from (DO_TAG) markers.
    # op_on "DO_TAG" "$list" → "true" if the tag appears in list, else "false".
    op_on() { [[ "$2" == *"($1)"* ]] && echo "true" || echo "false"; }
    _ops15=$(get_toml protein_structure operation_1_5 2>/dev/null || true)
    _ops69=$(get_toml protein_structure operation_6_9 2>/dev/null || true)
    # Ops 1-5: default all-true when list is absent (backward compat); honour tag when list present.
    if [[ -z "$_ops15" ]]; then
        export DO_EXTRACT=true DO_COPY_MODEL0=true DO_CIF_TO_PDB=true DO_UPDATE_HEADER=true DO_RENDER=true
    else
        export DO_EXTRACT=$(op_on "DO_EXTRACT" "$_ops15")
        export DO_COPY_MODEL0=$(op_on "DO_COPY_MODEL0" "$_ops15")
        export DO_CIF_TO_PDB=$(op_on "DO_CIF_TO_PDB" "$_ops15")
        export DO_UPDATE_HEADER=$(op_on "DO_UPDATE_HEADER" "$_ops15")
        export DO_RENDER=$(op_on "DO_RENDER" "$_ops15")
    fi
    # Ops 6-9: default false when list is absent.
    export DO_EXTRACT_METRICS=$(op_on "DO_EXTRACT_METRICS" "$_ops69")
    export DO_STRUCT_ALIGN=$(op_on "DO_STRUCT_ALIGN" "$_ops69")
    export DO_COMPARE_RENDER=$(op_on "DO_COMPARE_RENDER" "$_ops69")
    export DO_COMPARE_REPORT=$(op_on "DO_COMPARE_REPORT" "$_ops69")
    # Backgrounds: read list from TOML (newline-separated → space-separated for CLI).
    BACKGROUNDS=$(get_toml protein_structure backgrounds 2>/dev/null | tr '\n' ' ' | sed 's/ $//' || echo "black white")
    # Color versions: read list from TOML; empty string means "all versions from color config".
    COLOR_VERSIONS=$(get_toml protein_structure color_versions 2>/dev/null | tr '\n' ' ' | sed 's/ $//' || echo "")
    # Op 8 sub-types: space-separated; empty string falls back to "overlay" default.
    COMPARE_RENDER_TYPES=$(get_toml protein_structure operation_8_types 2>/dev/null | tr '\n' ' ' | sed 's/ $//' || echo "")
    # Multi-model genes: genes that also get extra models processed.
    MULTI_MODEL_GENES=$(get_toml protein_structure multi_model_genes 2>/dev/null | tr '\n' ' ' | sed 's/ $//' || echo "")
    MULTI_MODEL_NUMBERS=$(get_toml protein_structure multi_model_numbers 2>/dev/null | tr '\n' ' ' | sed 's/ $//' || echo "1")

    STRUCT_PIDS=()
    for RUN_DIR in "${RUN_DIRS[@]}"; do
        wait_for_slot
        log_info "  → $(basename "$RUN_DIR") (${THREADS_PER_RUN} threads)"
        bash "$STRUCT_PROCESSOR" \
            --input-dir "$RUN_DIR" \
            --color-config "$COLOR_CONFIG" \
            --gene-group "$GENE_GROUP" \
            --overwrite "$OVERWRITE" \
            --threads "$THREADS_PER_RUN" \
            --backgrounds "$BACKGROUNDS" \
            --color-versions "$COLOR_VERSIONS" \
            --compare-render-types "$COMPARE_RENDER_TYPES" \
            --multi-model-genes "$MULTI_MODEL_GENES" \
            --multi-model-numbers "$MULTI_MODEL_NUMBERS" &
        STRUCT_PIDS+=("$!")
    done
    for pid in "${STRUCT_PIDS[@]}"; do wait "$pid"; done
else
    log_warn "No structure data found under $PROT_DIR"
    log_info "Place AlphaFold3 zip files or AlphaFold3_Results/ folders under $PROT_DIR/<run>/"
fi

# ── Confidence + Interaction composite ───────────────────────────────────────
COMPOSE_ENABLED=$(get_toml protein_structure compose enabled 2>/dev/null || echo "false")
if [[ "$COMPOSE_ENABLED" == "true" && "$GENE_GROUP" == "DMP" ]]; then
    COMPOSE_BG=$(get_toml protein_structure compose background 2>/dev/null || echo "black")
    COMPOSE_VERSION=$(get_toml protein_structure compose interaction_version 2>/dev/null || echo "interaction")
    COMPOSE_SUBDIR=$(get_toml protein_structure compose output_subdir 2>/dev/null || echo "Combined_Confidence_Interaction")
    log_step "Composing confidence + interaction panels (version: $COMPOSE_VERSION, bg: $COMPOSE_BG)"
    conda run -n egg python3 \
        "$MODULES/08_protein_structure_analysis/compose_confidence_interaction.py" \
        --background "$COMPOSE_BG" \
        --interaction-version "$COMPOSE_VERSION" \
        --output-subdir "$COMPOSE_SUBDIR"
fi

log_step "Protein structure analysis complete: $GENE_GROUP"
teardown_logging

done
