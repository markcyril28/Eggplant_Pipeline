#!/bin/bash
# ============================================================================
# Program 9: CRISPR Off-Target Analysis Pipeline
# ============================================================================
# Operations (configured via [crispr].operations in TOML):
#   score_filter    — filter raw sgRNA scoring results by on-target score
#   blast_offtarget — BLAST-based genome-wide off-target search
#   cas_offinder    — Cas-OFFinder GPU-accelerated off-target search
#   report          — summary report + publication-quality visualizations
#
# Substage layout under 09_CRISPR_Off-Target_Analysis/{genome}/:
#   01_Raw_Scoring_Results_from_CRISPR-P_V2_0/
#   02_Filtered_Score_0.5/     (High + Moderate)
#   03_Filtered_Score_0.7/     (High only)
#   04_Off_Target_BLAST/       (genome-scoped)
#   05_Off_Target_Cas-OFFinder/ (genome-scoped)
#   06_Summary_Report/         (text report, guide_summary.csv, plots)
#
# Edit 9_crispr_v1CONFIG.toml to change gene groups, compute settings, and
# operations, then run:
#   bash i_crispr_v1.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# GENE_GROUPS, CPU, MAX_PARALLEL, OVERWRITE, and OPERATIONS are all loaded
# from 9_crispr_v1CONFIG.toml [pipeline] — edit gene_groups there.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# Load GENE_GROUPS from shared config (read before the per-group loop)
SHARED_CONFIG="$PIPELINE_DIR/9_crispr_v1CONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 9_crispr_v1CONFIG.toml" >&2
    exit 1
fi

wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done
}

# ─────────────────────────────────────────────────────────────
# Checks if an operation is in the OPERATIONS array
# ─────────────────────────────────────────────────────────────
op_enabled() {
    local op="$1"
    for o in "${OPERATIONS[@]}"; do [[ "$o" == "$op" ]] && return 0; done
    return 1
}

TMP_CONFIG_FILES=()

cleanup_tmp_configs() {
    local cfg
    for cfg in "${TMP_CONFIG_FILES[@]:-}"; do
        [[ -n "$cfg" && -f "$cfg" ]] && rm -f "$cfg"
    done
}

trap 'teardown_logging; cleanup_tmp_configs' EXIT

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# ── Config resolution ──────────────────────────────────────────
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_crispr_cfg_XXXXXX.toml")
    cat "$PIPELINE_DIR/9_crispr_v1CONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/09_crispr_analysis.toml" > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
BLAST_OPTIMAL_THREADS=$(get_toml pipeline compute "$MACHINE" blast_optimal_threads 2>/dev/null || echo "4")
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "$(( CPU / BLAST_OPTIMAL_THREADS ))")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
THREADS_PER_JOB=$BLAST_OPTIMAL_THREADS
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"
GENOME="$PIPELINE_DIR/$(get_toml reference eggplant_v4_1_genome)"

setup_logging

CRISPR_ENABLED=$(get_toml crispr enabled 2>/dev/null || echo "false")
if [[ "$CRISPR_ENABLED" != "True" && "$CRISPR_ENABLED" != "true" ]]; then
    log_info "CRISPR pipeline not enabled for $GENE_GROUP. Skipping."
    teardown_logging
    continue
fi

# Read operations list; fall back to running all if absent
mapfile -t OPERATIONS < <(get_toml crispr operations 2>/dev/null || printf '%s\n' "score_filter" "blast_offtarget" "cas_offinder" "report")

CRISPR_DIR="$BASE_DIR/09_CRISPR_Off-Target_Analysis"
mkdir -p "$CRISPR_DIR"

# Genome subdirectories to process
mapfile -t GENOME_NAMES < <(get_toml crispr genomes names 2>/dev/null || echo "GPE001970_SMEL5")

log_step "CRISPR Off-Target Analysis Pipeline: $GENE_GROUP"
log_info "Operations: ${OPERATIONS[*]}"

# ============================================================================
# Operation 1: Score-based filtering of raw sgRNA scoring results
# ============================================================================
if op_enabled "score_filter"; then
    log_step "Operation 1/4: Score filtering"

    mapfile -t SCORE_THRESHOLDS < <(get_toml crispr score_thresholds 2>/dev/null || printf '%s\n' 0.5 0.7)

    for genome_name in "${GENOME_NAMES[@]}"; do
        RAW_DIR="$CRISPR_DIR/${genome_name}/01_Raw_Scoring_Results_from_CRISPR-P_V2_0"
        if [[ ! -d "$RAW_DIR" ]]; then
            log_warn "Raw results not found: $RAW_DIR — skipping filtering for $genome_name"
            continue
        fi

        # Materialize tab-delimited siblings of any CRISPR-P v2.0 *.csv exports
        # so downstream stages always see a consistent TSV format. Idempotent:
        # only regenerates when the .tsv is missing or older than the .csv.
        python3 "$MODULES/09_crispr_analysis/csv_to_tsv.py" "$RAW_DIR" \
            | while IFS= read -r line; do log_info "$line"; done

        # Dynamic substage numbers: 0.5 → 02, 0.7 → 03, etc.
        n=2
        for threshold in "${SCORE_THRESHOLDS[@]}"; do
            SUBSTAGE=$(printf '%02d' "$n")
            FILTER_DIR="$CRISPR_DIR/${genome_name}/${SUBSTAGE}_Filtered_Score_${threshold}"
            log_info "  [$genome_name] score >= $threshold -> ${SUBSTAGE}_Filtered_Score_${threshold}"
            python3 "$MODULES/09_crispr_analysis/filter_scores.py" \
                --input-dir "$RAW_DIR" \
                --output-dir "$FILTER_DIR" \
                --threshold "$threshold"
            (( n++ ))
        done
    done
fi

# ============================================================================
# Operation 2: BLAST off-target analysis
# ============================================================================
if op_enabled "blast_offtarget"; then
    log_step "Operation 2/4: BLAST off-target analysis"

    # Read BLAST accuracy parameters from config (with most-accurate defaults)
    BLAST_WORD_SIZE=$(get_toml crispr blast word_size        2>/dev/null || echo "7")
    BLAST_EVALUE=$(get_toml crispr blast evalue              2>/dev/null || echo "10000")
    BLAST_DUST=$(get_toml crispr blast dust                  2>/dev/null || echo "false")
    BLAST_UNGAPPED=$(get_toml crispr blast ungapped          2>/dev/null || echo "true")
    BLAST_REWARD=$(get_toml crispr blast reward              2>/dev/null || echo "1")
    BLAST_PENALTY=$(get_toml crispr blast penalty            2>/dev/null || echo "-1")
    BLAST_MAX_TARGETS=$(get_toml crispr blast max_target_seqs 2>/dev/null || echo "10000")
    BLAST_QCOV=$(get_toml crispr blast qcov_hsp_perc        2>/dev/null || echo "80")
    BLAST_MM=$(get_toml crispr blast max_mismatches          2>/dev/null || echo "4")

    # Convert dust boolean to BLAST flag value
    [[ "$BLAST_DUST" == "false" || "$BLAST_DUST" == "False" ]] && BLAST_DUST="no" || BLAST_DUST="yes"

    mapfile -t GRNA_FASTAS < <(get_toml crispr offtarget "$GENE_GROUP" grna_fastas 2>/dev/null || true)

    for genome_name in "${GENOME_NAMES[@]}"; do
        BLAST_DIR="$CRISPR_DIR/${genome_name}/04_Off_Target_BLAST"
        BLAST_PIDS=()
        for grna in "${GRNA_FASTAS[@]}"; do
            GRNA_FULL="$PIPELINE_DIR/$grna"
            [[ -f "$GRNA_FULL" ]] || { log_warn "gRNA FASTA not found: $GRNA_FULL"; continue; }
            wait_for_slot
            bash "$MODULES/09_crispr_analysis/off_target_blast.sh" \
                --grna-fasta     "$GRNA_FULL" \
                --genome         "$GENOME" \
                --outdir         "$BLAST_DIR" \
                --threads        "$THREADS_PER_JOB" \
                --max-mismatches "$BLAST_MM" \
                --word-size      "$BLAST_WORD_SIZE" \
                --evalue         "$BLAST_EVALUE" \
                --dust           "$BLAST_DUST" \
                --ungapped       "$BLAST_UNGAPPED" \
                --reward         "$BLAST_REWARD" \
                --penalty        "$BLAST_PENALTY" \
                --max-target-seqs "$BLAST_MAX_TARGETS" \
                --qcov-hsp-perc  "$BLAST_QCOV" &
            BLAST_PIDS+=("$!")
        done
        for pid in "${BLAST_PIDS[@]}"; do wait "$pid"; done
    done
fi

# ============================================================================
# Operation 3: Cas-OFFinder off-target analysis
# ============================================================================
if op_enabled "cas_offinder"; then
    CASOFF_ENABLED=$(get_toml crispr cas_offinder enabled 2>/dev/null || echo "false")
    if [[ "$CASOFF_ENABLED" == "true" || "$CASOFF_ENABLED" == "True" ]]; then
        log_step "Operation 3/4: Cas-OFFinder off-target analysis"

        PAM=$(get_toml crispr cas_offinder pam             2>/dev/null || echo "NNNNNNNNNNNNNNNNNNNNNGG")
        DEVICE=$(get_toml crispr cas_offinder device          2>/dev/null || echo "G0")
        CASOFF_MM=$(get_toml crispr cas_offinder max_mismatches 2>/dev/null || echo "4")

        mapfile -t GRNA_FASTAS < <(get_toml crispr offtarget "$GENE_GROUP" grna_fastas 2>/dev/null || true)

        for genome_name in "${GENOME_NAMES[@]}"; do
            CASOFF_DIR="$CRISPR_DIR/${genome_name}/05_Off_Target_Cas-OFFinder"
            CASOFF_PIDS=()
            for grna in "${GRNA_FASTAS[@]}"; do
                GRNA_FULL="$PIPELINE_DIR/$grna"
                [[ -f "$GRNA_FULL" ]] || continue
                wait_for_slot
                bash "$MODULES/09_crispr_analysis/cas_offinder.sh" \
                    --grna-fasta  "$GRNA_FULL" \
                    --genome      "$GENOME" \
                    --outdir      "$CASOFF_DIR" \
                    --pam         "$PAM" \
                    --max-mismatches "$CASOFF_MM" \
                    --device      "$DEVICE" &
                CASOFF_PIDS+=("$!")
            done
            for pid in "${CASOFF_PIDS[@]}"; do wait "$pid"; done
        done
    else
        log_info "Cas-OFFinder not enabled — skipping."
    fi
fi

# ============================================================================
# Operation 4: Summary report + visualizations
# ============================================================================
if op_enabled "report"; then
    log_step "Operation 4/4: Summary report & visualizations"

    REPORT_DPI=$(get_toml crispr report dpi    2>/dev/null || echo "300")
    REPORT_FMT=$(get_toml crispr report format 2>/dev/null || echo "png")
    REPORT_Y_MAX_CAP=$(get_toml crispr report y_max_cap 2>/dev/null || echo "0")
    REPORT_OFFT_STRICT=$(get_toml crispr report offtarget_strict    2>/dev/null || echo "0")
    REPORT_OFFT_MOD=$(get_toml crispr report offtarget_moderate     2>/dev/null || echo "10")

    mapfile -t SCORE_THRESHOLDS < <(get_toml crispr score_thresholds 2>/dev/null || printf '%s\n' 0.5 0.7)

    # Scatter/table thresholds derive from [crispr].score_thresholds so the
    # filtering operation and the visualization stay in sync. Convention:
    # smallest = moderate (orange line), largest = strict (green line).
    SCORE_THRESHOLDS_SORTED=($(printf '%s\n' "${SCORE_THRESHOLDS[@]}" | sort -g))
    REPORT_SCORE_MOD="${SCORE_THRESHOLDS_SORTED[0]}"
    REPORT_SCORE_STRICT="${SCORE_THRESHOLDS_SORTED[-1]}"

    for genome_name in "${GENOME_NAMES[@]}"; do
        GENOME_DIR="$CRISPR_DIR/${genome_name}"
        REPORT_DIR="$GENOME_DIR/06_Summary_Report"
        BLAST_DIR="$GENOME_DIR/04_Off_Target_BLAST"
        CASOFF_DIR="$GENOME_DIR/05_Off_Target_Cas-OFFinder"

        # Build optional args for off-target dirs
        REPORT_EXTRAS=()
        [[ -d "$BLAST_DIR" ]]  && REPORT_EXTRAS+=(--blast-dir "$BLAST_DIR")
        [[ -d "$CASOFF_DIR" ]] && REPORT_EXTRAS+=(--casoff-dir "$CASOFF_DIR")

        log_info "  [$genome_name] Generating summary report..."
        python3 "$MODULES/09_crispr_analysis/generate_report.py" \
            --crispr-dir  "$GENOME_DIR" \
            --genome      "$genome_name" \
            --output-dir  "$REPORT_DIR" \
            --score-thresholds "${SCORE_THRESHOLDS[@]}" \
            "${REPORT_EXTRAS[@]}"

        # Generate plots from the summary CSV
        SUMMARY_CSV="$REPORT_DIR/guide_summary.csv"
        if [[ -f "$SUMMARY_CSV" ]]; then
            log_info "  [$genome_name] Generating visualizations (dpi=$REPORT_DPI, fmt=$REPORT_FMT)..."
            python3 "$MODULES/09_crispr_analysis/plot_crispr_results.py" \
                --summary-csv       "$SUMMARY_CSV" \
                --output-dir        "$REPORT_DIR" \
                --dpi               "$REPORT_DPI" \
                --format            "$REPORT_FMT" \
                --y-max-cap         "$REPORT_Y_MAX_CAP" \
                --score-strict      "$REPORT_SCORE_STRICT" \
                --score-moderate    "$REPORT_SCORE_MOD" \
                --offtarget-strict  "$REPORT_OFFT_STRICT" \
                --offtarget-moderate "$REPORT_OFFT_MOD"
        fi
    done
fi

log_step "CRISPR Off-Target Analysis Pipeline complete: $GENE_GROUP"
teardown_logging

done
