#!/bin/bash
# ============================================================================
# Program 9 v3: CRISPR KO Prediction Pipeline — PLANT-ONLY
# ============================================================================
# Plant-first redesign of the v2 pipeline. Uses CRISPR-P v2.0 (Liu 2017) as
# the SOLE on-target guide source and enforces plant-biology rules (U6/U3
# Pol III compatibility, GC window, plant long-3'UTR NMD) as mandatory.
#
# Module scripts under modules/09_crispr_analysis/v2/ are SHARED with v2 —
# v3 expresses plant-first behaviour only through this orchestrator, the
# TOML config 09_crispr_v3CONFIG.toml, and the group override
# config/{GROUP}/09_crispr_analysis_v3.toml.
#
# Mode:
#   mode = "crisprp_only"   (only supported mode in v3)
#     Symlinks pre-existing CRISPR-P v2.0 raw result TSVs into the arm's
#     01_Raw_Input/ folder, renames the score column to "crisprP_score",
#     then runs stages [1b] plant_filter, [3] off-target curation,
#     [4] indel outcome, [5] mutant transcripts, [6] protein consequence,
#     [7] plant NMD (50-nt + ≥350-nt 3'UTR), [8] composite ranking, and
#     [9] comparison scatter.
#
# Output tree:
#   09_CRISPR_v3/{genome}/
#     01_Raw_Input/                    ← symlinked CRISPR-P v2.0 results
#     01b_Plant_Filter/                ← TTTT / GC / 5' G|A filter (mandatory)
#     02_Plant_Rescorers/              ← DeepCRISPR / CRISPR-Local (pending)
#     03_OffTargets_CFD/               ← CFD score + paralog filter
#     04_Indels_inDelphi_Lindel/       ← indel OUTCOMES (mammalian, with caveat)
#     05_Transcripts_Biopython/        ← Biopython mutant CDS/mRNA
#     06_Protein_Biopython_Pfam/       ← Biopython + Pfam domain TSV
#     07_NMD_plantRules/               ← 50-nt EJC + plant ≥350-nt 3'UTR
#     08_Ranking_Composite/            ← plant-weighted KO score
#     09_Guide_Scatter/                ← on-target vs off-target scatter
#
# Edit 09_crispr_v3CONFIG.toml then run:
#   bash i_crispr_v3_pipeline.sh
# ============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"
# v3 SHARES v2's Python modules — bug fixes propagate automatically.
V2_MOD="$MODULES/09_crispr_analysis/v2"
PROJECT_ROOT="$PIPELINE_DIR"

# WSL2/NTFS workaround (see v2 orchestrator for rationale).
safe_mkdir() {
    local d
    for d in "$@"; do
        mkdir -p "$d" 2>/dev/null || [[ -d "$d" ]] || {
            echo "ERROR: Cannot create directory: $d" >&2
            return 1
        }
    done
}

safe_mkdir "$PROJECT_ROOT/logs"
source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# ─── Load GENE_GROUPS from v3 shared config ──────────────────────────────────
SHARED_CONFIG="$PIPELINE_DIR/09_crispr_v3CONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 09_crispr_v3CONFIG.toml" >&2
    exit 1
fi

# ─── Helpers ─────────────────────────────────────────────────────────────────
wait_for_slot() {
    local limit="${1:-$MAX_PARALLEL}"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

op_enabled() {
    local op="$1"
    for o in "${OPERATIONS[@]}"; do [[ "$o" == "$op" ]] && return 0; done
    return 1
}

# ensure_indelphi_env / preflight — identical to v2 (same conda env name).
ensure_indelphi_env() {
    local env_name="$1"
    [[ -z "$env_name" ]] && return 0
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env_name"; then
        log_info "[inDelphi] Legacy env '$env_name' already exists."
        return 0
    fi
    log_info "[inDelphi] Creating legacy conda env '$env_name' (scikit-learn=0.20.0, Python 3.7)..."
    # Pin channels to match setup_conda_crispr_v3.sh Step 5 so auto-creation and
    # scripted creation produce bit-identical envs (prevents sklearn build drift
    # that would break inDelphi's pickled model weights).
    if ! conda create -y -n "$env_name" \
            -c conda-forge -c bioconda --strict-channel-priority \
            python=3.7 scikit-learn=0.20.0 numpy=1.19 pandas=1.1; then
        log_warn "[inDelphi] 'conda create' failed for '$env_name'. inDelphi will be skipped."
        return 0
    fi
    log_info "[inDelphi] Legacy env '$env_name' ready (inDelphi loaded from bundled tools dir)."
}

preflight_check_indelphi_sklearn() {
    local env_name="$1"
    [[ -z "$env_name" ]] && return 0
    local actual_ver
    actual_ver=$(conda run -n "$env_name" python3 -c \
        "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo "UNAVAILABLE")
    if [[ "$actual_ver" == "0.18.1" || "$actual_ver" == "0.20.0" ]]; then
        log_info "[pre-flight] inDelphi legacy env '$env_name': scikit-learn $actual_ver — OK"
    else
        log_error "[pre-flight] inDelphi legacy env '$env_name' has scikit-learn $actual_ver — INCOMPATIBLE (requires 0.18.1 or 0.20.0 exactly). Stage 04 will fall back to Lindel-only."
    fi
}

TMP_CONFIG_FILES=()
cleanup_tmp_configs() {
    local cfg
    for cfg in "${TMP_CONFIG_FILES[@]:-}"; do
        [[ -n "$cfg" && -f "$cfg" ]] && rm -f "$cfg"
    done
}
trap 'cleanup_tmp_configs; safe_teardown_logging' EXIT

# ─────────────────────────────────────────────────────────────────────────────
# run_plant_stages ARM_DIR GENOME_NAME [_MP]
#   Runs v3 stages [1b]–[9] on an arm root directory. Stage [1] is NOT
#   included because v3 sources guides from pre-existing CRISPR-P v2.0
#   raw results linked into ARM_DIR/01_Raw_Input/ by the caller.
#
#   Stage module scripts are the SAME as v2 (V2_MOD path); only the config
#   keys and output folder differ.
# ─────────────────────────────────────────────────────────────────────────────
run_plant_stages() {
    local ARM_DIR="$1"
    local GENOME_NAME="$2"
    local _MP="${3:-$MAX_PARALLEL}"

    # ── [1b] Plant sgRNA pre-filter (MANDATORY in v3) ───────────────────────
    # Always runs regardless of the `operations` list — v3 contract.
    local STAGE_DIR="$ARM_DIR/01b_Plant_Filter"
    safe_mkdir "$STAGE_DIR"
    local PIDS=()
    local -a _S1_TSVS=()
    # v3 stage-1 source directory is 01_Raw_Input/ (symlinked CRISPR-P TSVs).
    for tsv in "$ARM_DIR/01_Raw_Input/"*.tsv "$ARM_DIR/01_Raw_Input/"*.csv; do
        [[ -f "$tsv" ]] && _S1_TSVS+=("$tsv")
    done
    if [[ ${#_S1_TSVS[@]} -eq 0 ]]; then
        log_warn "[01b_plant] No CRISPR-P raw TSV/CSV under $ARM_DIR/01_Raw_Input/ — pipeline aborted for this arm."
        return 1
    fi
    for tsv in "${_S1_TSVS[@]}"; do
        wait_for_slot "$_MP"
        $CONDA_RUN python3 "$V2_MOD/01b_plant_sgrna_filter.py" \
            --input                  "$tsv" \
            --outdir                 "$STAGE_DIR" \
            --termination-run-length "$PLANT_TERM_RUN" \
            --termination-action     "$PLANT_TERM_ACTION" \
            --promoter-type          "$PLANT_PROMOTER_TYPE" \
            --promoter-action        "$PLANT_PROMOTER_ACTION" \
            --gc-min                 "$PLANT_GC_MIN" \
            --gc-max                 "$PLANT_GC_MAX" \
            --enabled                "true" \
            "$OVERWRITE_FLAG" &
        PIDS+=("$!")
    done
    for pid in "${PIDS[@]}"; do wait "$pid"; done

    # ── [2] Plant-trained rescorers (optional; pending tool bundling) ───────
    # v3 does NOT dispatch mammalian rescorers. When plant_scorer.enabled AND
    # plant_predictors is non-empty, this stage runs DeepCRISPR / CRISPR-Local
    # into 02_Plant_Rescorers/. Otherwise it is skipped.
    if op_enabled "rescore_ontarget" && \
       [[ "${PLANT_SCORER_ENABLED:-false}" == "true" ]] && \
       (( ${#PLANT_PREDICTORS[@]} > 0 )); then
        local PR_DIR="$ARM_DIR/02_Plant_Rescorers"
        safe_mkdir "$PR_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/01b_Plant_Filter/"*.plant_filtered.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            local -a _PR_ARGS=(
                --input          "$tsv"
                --outdir         "$PR_DIR"
                --predictors     "${PLANT_PREDICTORS[@]}"
                --flag-threshold "$FLAG_THRESH"
                --workers        "$INNER_WORKERS"
                "$OVERWRITE_FLAG"
            )
            # Thread DeepCRISPR's TF-1.x env + CRISPR-Local's py2 env to
            # the module so each scorer's dispatch path can subprocess into
            # its correct runtime.
            [[ -n "${DEEPCRISPR_ENV:-}" ]] && \
                _PR_ARGS+=(--deepcrispr-env "$DEEPCRISPR_ENV")
            [[ -n "${CRISPR_LOCAL_ENV:-}" ]] && \
                _PR_ARGS+=(--crispr-local-env "$CRISPR_LOCAL_ENV")
            $CONDA_RUN python3 "$V2_MOD/02_rescore_ontarget.py" "${_PR_ARGS[@]}" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    else
        log_info "[02_rescore] Plant-trained rescore skipped (no plant_predictors bundled)."
    fi

    # ── [3] Curate off-targets ───────────────────────────────────────────────
    if op_enabled "curate_offtargets"; then
        local STAGE_DIR="$ARM_DIR/03_OffTargets_CFD"
        # v3 off-target source: prefer 02_Plant_Rescorers output when present,
        # otherwise fall back to 01b_Plant_Filter output.
        local OFFTARGET_TSV=""
        for _ot_candidate in \
                "$ARM_DIR/01_Raw_Input/offtargets/offtargets.tsv" \
                "$ARM_DIR/01_Raw_Input/offtargets.tsv" \
                "$ARM_DIR/01_Raw_Input/"*offtarget*.tsv \
                "$ARM_DIR/01_Raw_Input/"*offtarget*.csv; do
            [[ -f "$_ot_candidate" ]] && { OFFTARGET_TSV="$_ot_candidate"; break; }
        done
        [[ -z "$OFFTARGET_TSV" ]] && OFFTARGET_TSV="$ARM_DIR/01_Raw_Input/offtargets.tsv"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        # Prefer rescored output; fall back to plant-filtered
        local -a _S2_OUT=()
        for tsv in "$ARM_DIR/02_Plant_Rescorers/"*.rescored.tsv; do
            [[ -f "$tsv" ]] && _S2_OUT+=("$tsv")
        done
        if [[ ${#_S2_OUT[@]} -eq 0 ]]; then
            for tsv in "$ARM_DIR/01b_Plant_Filter/"*.plant_filtered.tsv; do
                [[ -f "$tsv" ]] && _S2_OUT+=("$tsv")
            done
        fi
        for tsv in "${_S2_OUT[@]:-}"; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            $CONDA_RUN python3 "$V2_MOD/03_curate_offtargets.py" \
                --input                  "$tsv" \
                --offtarget-tsv          "$OFFTARGET_TSV" \
                --outdir                 "$STAGE_DIR" \
                --paralog-patterns       "${PARALOG_PATTERNS[@]}" \
                "${PARALOG_ID_ARGS[@]}" \
                --paralog-hit-threshold  "$PARALOG_HIT_THRESH" \
                --cfd-sum-threshold      "$CFD_SUM_THRESH" \
                "$OVERWRITE_FLAG" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [4] Indel outcome prediction (mammalian, with plant caveat) ──────────
    if op_enabled "predict_indels"; then
        local STAGE_DIR="$ARM_DIR/04_Indels_inDelphi_Lindel"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/03_OffTargets_CFD/"*.curated.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            local _stem; _stem="$(basename "$tsv")"; _stem="${_stem%%.*}"
            local INDEL_ARGS=(
                --input                  "$tsv"
                --outdir                 "$STAGE_DIR"
                --predictors             "${INDEL_PREDICTORS[@]}"
                --indelphi-cell-type     "$INDELPHI_CELL"
                --frameshift-threshold   "$FS_THRESH"
                --top-outcomes           "$TOP_OUTCOMES"
                --indelphi-conda-env     "$CONDA_ENV_INDELPHI"
                --workers                "$INNER_WORKERS"
                --genome-fasta           "$GENOME"
                --gene-group             "$GENE_GROUP"
                --gene-id                "$_stem"
                "$OVERWRITE_FLAG"
            )
            [[ -n "${GTF_FOR_STAGE:-}" && -f "$GTF_FOR_STAGE" ]] && \
                INDEL_ARGS+=(--gtf "$GTF_FOR_STAGE")
            $CONDA_RUN python3 "$V2_MOD/04_predict_indels.py" "${INDEL_ARGS[@]}" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [5] Rebuild transcripts ──────────────────────────────────────────────
    if op_enabled "rebuild_transcripts"; then
        local STAGE_DIR="$ARM_DIR/05_Transcripts_Biopython"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/04_Indels_inDelphi_Lindel/"*.indels.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            local _stem; _stem="$(basename "$tsv")"; _stem="${_stem%%.*}"
            local REBUILD_ARGS=(
                --input        "$tsv"
                --outdir       "$STAGE_DIR"
                --genome-fasta "$GENOME"
                --gene-group   "$GENE_GROUP"
                --gene-id      "$_stem"
                --top-indels   "$TOP_INDELS"
                "$OVERWRITE_FLAG"
            )
            [[ -n "${GTF_FOR_STAGE:-}" && -f "$GTF_FOR_STAGE" ]] && \
                REBUILD_ARGS+=(--gtf "$GTF_FOR_STAGE")
            $CONDA_RUN python3 "$V2_MOD/05_rebuild_transcripts.py" "${REBUILD_ARGS[@]}" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [6] Protein consequence ──────────────────────────────────────────────
    if op_enabled "protein_consequence"; then
        local STAGE_DIR="$ARM_DIR/06_Protein_Biopython_Pfam"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/05_Transcripts_Biopython/"*.transcripts.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            local PROTEIN_ARGS=(
                --input              "$tsv"
                --outdir             "$STAGE_DIR"
                --flag-domain-hits   "$FLAG_DOMAINS"
                "$OVERWRITE_FLAG"
            )
            [[ -n "${DOMAIN_TSV:-}" && -f "$DOMAIN_TSV" ]] && \
                PROTEIN_ARGS+=(--domain-tsv "$DOMAIN_TSV")
            $CONDA_RUN python3 "$V2_MOD/06_protein_consequence.py" "${PROTEIN_ARGS[@]}" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [6b] ESMFold local structure prediction (single-GPU, serial) ─────────
    # One worker call processes every flagged sequence across all *.protein.tsv
    # so the ~8 GB ESMFold weights load only once. Skipped unless backend=local.
    if op_enabled "fold_protein"; then
        if [[ "$ESMFOLD_BACKEND" == "local" ]]; then
            local STAGE_DIR="$ARM_DIR/06_Protein_Biopython_Pfam"
            local PDB_DIR="$STAGE_DIR/esmfold_structures"
            safe_mkdir "$PDB_DIR"
            $CONDA_RUN python3 "$V2_MOD/06b_run_esmfold.py" \
                --protein-tsv-dir    "$STAGE_DIR" \
                --transcripts-dir    "$ARM_DIR/05_Transcripts_Biopython" \
                --outdir             "$PDB_DIR" \
                --esmfold-env        "$ESMFOLD_ENV" \
                --chunk-size         "$ESMFOLD_CHUNK_SIZE" \
                --max-protein-length "$ESMFOLD_MAX_LEN" \
                "$OVERWRITE_FLAG" || \
                log_warn "[06b_esmfold] Folding stage failed for $ARM_DIR (continuing)"
        else
            log_info "[06b_esmfold] esmfold_backend=$ESMFOLD_BACKEND (not 'local') — skipping local fold."
        fi
    fi

    # ── [7] Plant-aware NMD (50-nt + ≥350-nt 3'UTR) ──────────────────────────
    if op_enabled "predict_nmd"; then
        local STAGE_DIR="$ARM_DIR/07_NMD_plantRules"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/06_Protein_Biopython_Pfam/"*.protein.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            local _stem; _stem="$(basename "$tsv")"; _stem="${_stem%%.*}"
            local NMD_ARGS=(
                --input                    "$tsv"
                --outdir                   "$STAGE_DIR"
                --gene-group               "$GENE_GROUP"
                --gene-id                  "$_stem"
                --ptc-distance-threshold   "$PTC_DIST"
                --long-3utr-threshold      "$LONG_3UTR_THRESH"
                "$OVERWRITE_FLAG"
            )
            [[ -n "${GTF:-}" && -f "$GTF" ]] && NMD_ARGS+=(--gtf "$GTF")
            $CONDA_RUN python3 "$V2_MOD/07_predict_nmd.py" "${NMD_ARGS[@]}" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [8] Rank guides (plant-biased weights from [crispr_v3.rank_guides]) ──
    if op_enabled "rank_guides"; then
        local STAGE_DIR="$ARM_DIR/08_Ranking_Composite"
        local PER_GENE_DIR="$STAGE_DIR/per_gene"
        safe_mkdir "$PER_GENE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/07_NMD_plantRules/"*.nmd.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            $CONDA_RUN python3 "$V2_MOD/08_rank_guides.py" \
                --input         "$tsv" \
                --outdir        "$PER_GENE_DIR" \
                --weights       "$WEIGHTS_JSON" \
                --top-n         "$TOP_N" \
                --output-format "$OUT_FMT" \
                --dpi           "$REPORT_DPI" \
                --format        "$REPORT_FMT" \
                --cfd-sum-cap   "$CFD_SUM_CAP" \
                "$OVERWRITE_FLAG" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done

        # Aggregate
        $CONDA_RUN python3 "$V2_MOD/08_rank_guides.py" \
            --aggregate-dir "$PER_GENE_DIR" \
            --outdir        "$STAGE_DIR" \
            --weights       "$WEIGHTS_JSON" \
            --top-n         "$TOP_N" \
            --output-format "$OUT_FMT" \
            --dpi           "$REPORT_DPI" \
            --format        "$REPORT_FMT" \
            --cfd-sum-cap   "$CFD_SUM_CAP" \
            "$OVERWRITE_FLAG" || \
            log_warn "[08_rank] Aggregate step failed for $ARM_DIR"
    fi

    # ── [9] Guide comparison scatter (plant-only scorer variants) ────────────
    if op_enabled "comparison_scatter"; then
        local STAGE_DIR="$ARM_DIR/09_Guide_Scatter"
        local RANKED_DIR="$ARM_DIR/08_Ranking_Composite/per_gene"
        safe_mkdir "$STAGE_DIR"
        if compgen -G "$RANKED_DIR/*.ranked.tsv" > /dev/null 2>&1 || \
           compgen -G "$RANKED_DIR/*_ranked_guides.tsv" > /dev/null 2>&1; then
            _SCATTER_EXTRA=()
            [[ -n "$SCATTER_VARIANTS_ARG" ]] && \
                _SCATTER_EXTRA=(--score-variants "$SCATTER_VARIANTS_ARG")
            $CONDA_RUN python3 "$V2_MOD/09_comparison_scatter.py" \
                --input-dir       "$RANKED_DIR" \
                --outdir          "$STAGE_DIR" \
                --top-n           "$SCATTER_TOP_N" \
                --y-max-cap       "$SCATTER_Y_MAX_CAP" \
                --tier-high       "$SCATTER_TIER_HIGH" \
                --tier-moderate   "$SCATTER_TIER_MOD" \
                --x-score-column  "$SCATTER_X_COLUMN" \
                --x-axis-high     "$SCATTER_X_HIGH" \
                --x-axis-moderate "$SCATTER_X_MOD" \
                --x-axis-max      "$SCATTER_X_MAX" \
                --dpi             "$REPORT_DPI" \
                --format          "$REPORT_FMT" \
                "${_SCATTER_EXTRA[@]}" \
                "$OVERWRITE_FLAG" || \
                log_warn "[09_scatter] Comparison scatter failed for $ARM_DIR"
        else
            log_warn "[09_scatter] No ranked TSVs under $RANKED_DIR — skipping."
        fi
    fi
}

# =============================================================================
# Logging setup
# =============================================================================
setup_logging

# Track which inDelphi envs have already been provisioned this run, so the
# per-group loop below does not re-check the same env name multiple times
# when several gene groups share it.
_INDELPHI_ENVS_SEEN=""

# =============================================================================
# Per gene-group loop
# =============================================================================
for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# ── Config resolution ────────────────────────────────────────────────────────
MERGE_TOML="$MODULES/utils/merge_toml.py"
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
if [[ -d "$CONFIG_DIR" && -f "$CONFIG_DIR/09_crispr_analysis_v3.toml" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_crispr_v3_cfg_XXXXXX.toml")
    python3 "$MERGE_TOML" \
        "$PIPELINE_DIR/09_crispr_v3CONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/09_crispr_analysis_v3.toml" \
        > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$SHARED_CONFIG"
    log_warn "No group override config/${GENE_GROUP}/09_crispr_analysis_v3.toml — using shared config only."
fi

# ── Compute settings ─────────────────────────────────────────────────────────
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "4")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
INNER_WORKERS=$(( CPU / MAX_PARALLEL ))
(( INNER_WORKERS < 1 )) && INNER_WORKERS=1
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")
[[ "$OVERWRITE" == "true" || "$OVERWRITE" == "True" ]] && OVERWRITE_FLAG="--overwrite" || OVERWRITE_FLAG="--no-overwrite"

_base_dir_rel=$(get_toml general base_dir 2>/dev/null || echo "")
[[ -z "$_base_dir_rel" ]] && _base_dir_rel="III_RESULT/$GENE_GROUP"
BASE_DIR="$PIPELINE_DIR/$_base_dir_rel"

# ── Genome / GTF / domain paths ──────────────────────────────────────────────
_default_genome_key=$(get_toml crispr_v3 rebuild_transcripts genome_key 2>/dev/null || echo "")
_default_gtf_key=$(get_toml    crispr_v3 rebuild_transcripts gtf_key    2>/dev/null || echo "")
_default_genome_rel=$( [[ -n "$_default_genome_key" ]] && get_toml reference "$_default_genome_key" 2>/dev/null || echo "")
DEFAULT_GENOME="${_default_genome_rel:+$PIPELINE_DIR/$_default_genome_rel}"
DEFAULT_GTF_REL=$( [[ -n "$_default_gtf_key" ]] && get_toml reference "$_default_gtf_key" 2>/dev/null || echo "")
DEFAULT_GTF="${DEFAULT_GTF_REL:+$PIPELINE_DIR/$DEFAULT_GTF_REL}"

DOMAIN_TSV_KEY=$(get_toml crispr_v3 protein_consequence domain_tsv_key 2>/dev/null || echo "")
DOMAIN_TSV=""
if [[ -n "$DOMAIN_TSV_KEY" ]]; then
    DOMAIN_TSV_REL=$(get_toml reference "$DOMAIN_TSV_KEY" 2>/dev/null || echo "")
    [[ -n "$DOMAIN_TSV_REL" ]] && DOMAIN_TSV="$PIPELINE_DIR/$DOMAIN_TSV_REL"
fi

CONDA_ENV=$(get_toml crispr_v3 conda_env 2>/dev/null || echo "crispr_v2")
CONDA_RUN="conda run -n $CONDA_ENV --no-capture-output"
CONDA_ENV_INDELPHI=$(get_toml crispr_v3 conda_env_indelphi 2>/dev/null || echo "")

ENABLED=$(get_toml crispr_v3 enabled 2>/dev/null || echo "false")
if [[ "$ENABLED" != "True" && "$ENABLED" != "true" ]]; then
    log_info "CRISPR v3 pipeline not enabled for $GENE_GROUP. Skipping."
    continue
fi

# Provision / preflight-check the inDelphi legacy env using the MERGED
# config value so a group-level override of conda_env_indelphi is honoured
# (the previous one-time setup above the loop read the shared TOML only).
# Runs AFTER the ENABLED gate so a disabled group does not pay the cost of
# creating / probing the legacy sklearn 0.20 env. Dedup by env name so
# several enabled groups sharing the same env do not re-check it.
if [[ -n "$CONDA_ENV_INDELPHI" && " $_INDELPHI_ENVS_SEEN " != *" $CONDA_ENV_INDELPHI "* ]]; then
    ensure_indelphi_env "$CONDA_ENV_INDELPHI" || true
    preflight_check_indelphi_sklearn "$CONDA_ENV_INDELPHI"
    _INDELPHI_ENVS_SEEN="$_INDELPHI_ENVS_SEEN $CONDA_ENV_INDELPHI"
fi

# ── Mode ─────────────────────────────────────────────────────────────────────
MODE=$(get_toml crispr_v3 mode 2>/dev/null || echo "crisprp_only")
if [[ "$MODE" != "crisprp_only" ]]; then
    log_warn "v3 only supports mode='crisprp_only' (got '$MODE'); forcing."
    MODE="crisprp_only"
fi

# ── Operations, with MANDATORY plant_filter re-injection ─────────────────────
mapfile -t OPERATIONS < <(
    get_toml crispr_v3 operations 2>/dev/null || \
    printf '%s\n' plant_filter curate_offtargets predict_indels \
                   rebuild_transcripts protein_consequence fold_protein \
                   predict_nmd rank_guides comparison_scatter
)
# v3 contract: plant_filter is non-optional. Re-inject even if the user
# commented it out in the TOML.
_has_pf=0
for _op in "${OPERATIONS[@]}"; do [[ "$_op" == "plant_filter" ]] && _has_pf=1; done
(( _has_pf == 0 )) && OPERATIONS=("plant_filter" "${OPERATIONS[@]}")

CRISPR_DIR="$BASE_DIR/09_CRISPR_v3"
safe_mkdir "$CRISPR_DIR"

mapfile -t GENOME_NAMES < <(
    get_toml crispr_v3 genomes names 2>/dev/null || echo "GPE001970_SMEL5"
)

log_step "CRISPR v3 KO Prediction Pipeline (PLANT-ONLY): $GENE_GROUP"
log_info "Mode:          $MODE"
log_info "Operations:    ${OPERATIONS[*]}"
log_info "Conda env:     $CONDA_ENV"
log_info "CPU / MAX_PARALLEL / INNER_WORKERS: $CPU / $MAX_PARALLEL / $INNER_WORKERS"

# ── Plant sgRNA pre-filter params ────────────────────────────────────────────
PLANT_TERM_RUN=$(get_toml       crispr_v3 plant_filter termination_run_length 2>/dev/null || echo "4")
PLANT_TERM_ACTION=$(get_toml    crispr_v3 plant_filter termination_action     2>/dev/null || echo "reject")
PLANT_PROMOTER_TYPE=$(get_toml  crispr_v3 plant_filter promoter_type          2>/dev/null || echo "U6")
PLANT_PROMOTER_ACTION=$(get_toml crispr_v3 plant_filter promoter_action       2>/dev/null || echo "flag")
PLANT_GC_MIN=$(get_toml         crispr_v3 plant_filter gc_min                 2>/dev/null || echo "0.30")
PLANT_GC_MAX=$(get_toml         crispr_v3 plant_filter gc_max                 2>/dev/null || echo "0.70")

# ── Plant-trained rescorers ──────────────────────────────────────────────────
PLANT_SCORER_ENABLED=$(get_toml crispr_v3 plant_scorer enabled 2>/dev/null || echo "false")
mapfile -t PLANT_PREDICTORS < <(get_toml crispr_v3 plant_scorer plant_predictors 2>/dev/null || true)
# Isolated TF-1.x env for DeepCRISPR (empty ⇒ direct-import attempt only).
DEEPCRISPR_ENV=$(get_toml crispr_v3 plant_scorer deepcrispr_env 2>/dev/null || echo "")
# Python-2 env for CRISPR-Local's rs2_score_calculator.py (Python 2).
CRISPR_LOCAL_ENV=$(get_toml crispr_v3 plant_scorer crispr_local_env 2>/dev/null || echo "")
FLAG_THRESH=$(get_toml crispr_v3 rescore_ontarget flag_threshold 2>/dev/null || echo "0.3")

# ── Paralog / off-target curation ────────────────────────────────────────────
PARALOG_HIT_THRESH=$(get_toml   crispr_v3 curate_offtargets paralog_hit_threshold 2>/dev/null || echo "1")
CFD_SUM_THRESH=$(get_toml       crispr_v3 curate_offtargets cfd_sum_threshold     2>/dev/null || echo "0.2")
# Read patterns as a proper TOML array (one-per-line via mapfile) so a
# multi-word regex like "^Smel(?:DMP|HAP2)\\d+$" is preserved intact.
# The previous `read -ra <<< string` form word-split on whitespace and
# would have broken any pattern containing a space.
mapfile -t PARALOG_PATTERNS < <(get_toml crispr_v3 curate_offtargets paralog_patterns 2>/dev/null || true)
(( ${#PARALOG_PATTERNS[@]} == 0 )) && PARALOG_PATTERNS=("DMP")
mapfile -t PARALOG_IDS < <(get_toml crispr_v3 curate_offtargets paralog_gene_ids 2>/dev/null || true)
PARALOG_ID_ARGS=()
[[ ${#PARALOG_IDS[@]} -gt 0 ]] && PARALOG_ID_ARGS=(--paralog-gene-ids "${PARALOG_IDS[@]}")

# ── Indel outcome prediction ─────────────────────────────────────────────────
INDELPHI_CELL=$(get_toml  crispr_v3 predict_indels indelphi_cell_type  2>/dev/null || echo "mESC")
FS_THRESH=$(get_toml      crispr_v3 predict_indels frameshift_threshold 2>/dev/null || echo "0.5")
TOP_OUTCOMES=$(get_toml   crispr_v3 predict_indels top_outcomes         2>/dev/null || echo "5")
# Read predictors as a proper TOML array (one-per-line via mapfile). The
# previous `read -ra <<< "$RAW"` form only consumed the first line of
# parse_toml.py's newline-separated output, silently dropping every
# predictor after the first (e.g. `["inDelphi", "Lindel"]` collapsed to
# just "inDelphi" and Lindel was never dispatched to stage [4]).
mapfile -t INDEL_PREDICTORS < <(get_toml crispr_v3 predict_indels predictors 2>/dev/null || printf '%s\n' inDelphi Lindel)
(( ${#INDEL_PREDICTORS[@]} == 0 )) && INDEL_PREDICTORS=(inDelphi Lindel)

# ── Transcript rebuild ───────────────────────────────────────────────────────
GTF_KEY=$(get_toml crispr_v3 rebuild_transcripts gtf_key 2>/dev/null || echo "")
GTF_FOR_STAGE_FROM_KEY=""
if [[ -n "$GTF_KEY" ]]; then
    GTF_REL2=$(get_toml reference "$GTF_KEY" 2>/dev/null || echo "")
    [[ -n "$GTF_REL2" ]] && GTF_FOR_STAGE_FROM_KEY="$PIPELINE_DIR/$GTF_REL2"
fi
TOP_INDELS=$(get_toml crispr_v3 rebuild_transcripts top_indels 2>/dev/null || echo "3")

FLAG_DOMAINS=$(get_toml crispr_v3 protein_consequence flag_domain_hits 2>/dev/null || echo "true")

# ── ESMFold local-fold backend (stage [6b]) ──────────────────────────────────
# Config keys live under [crispr_v3.protein_consequence] alongside the
# domain-flag knobs, since folding consumes the structure_flag column.
ESMFOLD_BACKEND=$(get_toml    crispr_v3 protein_consequence esmfold_backend     2>/dev/null || echo "local")
ESMFOLD_ENV=$(get_toml        crispr_v3 protein_consequence esmfold_env         2>/dev/null || echo "esmfold")
ESMFOLD_CHUNK_SIZE=$(get_toml crispr_v3 protein_consequence esmfold_chunk_size  2>/dev/null || echo "64")
ESMFOLD_MAX_LEN=$(get_toml    crispr_v3 protein_consequence esmfold_max_length  2>/dev/null || echo "400")

# ── NMD ──────────────────────────────────────────────────────────────────────
PTC_DIST=$(get_toml         crispr_v3 predict_nmd ptc_distance_threshold 2>/dev/null || echo "50")
# v3 default: plant long-3'UTR rule ON (350 nt, Kerényi 2008).
LONG_3UTR_THRESH=$(get_toml crispr_v3 predict_nmd long_3utr_threshold    2>/dev/null || echo "350")

# ── Composite KO ranking ─────────────────────────────────────────────────────
_w_ontarget=$(get_toml  crispr_v3 rank_guides weights w_ontarget   2>/dev/null || echo "0.20")
_w_frameshift=$(get_toml crispr_v3 rank_guides weights w_frameshift 2>/dev/null || echo "0.30")
_w_nmd=$(get_toml       crispr_v3 rank_guides weights w_nmd        2>/dev/null || echo "0.25")
_w_offtarget=$(get_toml crispr_v3 rank_guides weights w_offtarget  2>/dev/null || echo "0.15")
_w_domain=$(get_toml    crispr_v3 rank_guides weights w_domain     2>/dev/null || echo "0.10")
WEIGHTS_JSON=$(printf '{"w_ontarget":%s,"w_frameshift":%s,"w_nmd":%s,"w_offtarget":%s,"w_domain":%s}' \
    "$_w_ontarget" "$_w_frameshift" "$_w_nmd" "$_w_offtarget" "$_w_domain")

TOP_N=$(get_toml       crispr_v3 rank_guides top_n         2>/dev/null || echo "10")
OUT_FMT=$(get_toml     crispr_v3 rank_guides output_format 2>/dev/null || echo "tsv")
CFD_SUM_CAP=$(get_toml crispr_v3 rank_guides cfd_sum_cap   2>/dev/null || echo "5.0")
REPORT_DPI=$(get_toml  crispr_v3 report dpi                2>/dev/null || echo "600")
REPORT_FMT=$(get_toml  crispr_v3 report format             2>/dev/null || echo "png")

# ── Scatter ──────────────────────────────────────────────────────────────────
SCATTER_TOP_N=$(get_toml     crispr_v3 comparison_scatter top_n         2>/dev/null || echo "3")
SCATTER_Y_MAX_CAP=$(get_toml crispr_v3 comparison_scatter y_max_cap     2>/dev/null || echo "100")
SCATTER_TIER_HIGH=$(get_toml crispr_v3 comparison_scatter tier_high     2>/dev/null || echo "0.7")
SCATTER_TIER_MOD=$(get_toml  crispr_v3 comparison_scatter tier_moderate 2>/dev/null || echo "0.5")
SCATTER_X_COLUMN=$(get_toml  crispr_v3 comparison_scatter x_score_column  2>/dev/null || echo "crisprP_score")
SCATTER_X_HIGH=$(get_toml    crispr_v3 comparison_scatter x_axis_high     2>/dev/null || echo "40")
SCATTER_X_MOD=$(get_toml     crispr_v3 comparison_scatter x_axis_moderate 2>/dev/null || echo "30")
SCATTER_X_MAX=$(get_toml     crispr_v3 comparison_scatter x_axis_max      2>/dev/null || echo "100")
mapfile -t _SCATTER_VARIANTS < <(get_toml crispr_v3 comparison_scatter score_variants 2>/dev/null || true)
SCATTER_VARIANTS_ARG=""
if (( ${#_SCATTER_VARIANTS[@]} > 0 )); then
    SCATTER_VARIANTS_ARG=$(IFS='|'; echo "${_SCATTER_VARIANTS[*]}")
fi

# ── CRISPR-P v2.0 raw dir template ───────────────────────────────────────────
CRISPRP_RAW_TEMPLATE=$(get_toml crispr_v3 crispr_p raw_dir_template 2>/dev/null || \
    echo "III_RESULT/{GROUP}/09_CRISPR_Off-Target_Analysis/{GENOME}/01_Raw_Scoring_Results_from_CRISPR-P_V2_0")

# =============================================================================
# Per-genome loop
# =============================================================================
for genome_name in "${GENOME_NAMES[@]}"; do

GENOME_BASE="$CRISPR_DIR/${genome_name}"

# Per-genome FASTA / GTF resolution
_genome_key=$(echo "${genome_name}" | tr '.' '_')
GENOME_REL=$(get_toml reference "${_genome_key}_genome" 2>/dev/null || echo "")
if [[ -n "$GENOME_REL" ]]; then
    GENOME="$PIPELINE_DIR/$GENOME_REL"
else
    GENOME="$DEFAULT_GENOME"
fi
GTF_REL_LOOP=$(get_toml reference "${_genome_key}_gtf" 2>/dev/null || echo "")
if [[ -n "$GTF_REL_LOOP" ]]; then
    GTF="$PIPELINE_DIR/$GTF_REL_LOOP"
else
    GTF="$DEFAULT_GTF"
fi
GTF_FOR_STAGE="${GTF_FOR_STAGE_FROM_KEY:-$GTF}"

# ── Single arm (crisprp_only mode) ───────────────────────────────────────────
ARM="$GENOME_BASE/crisprp_only"
safe_mkdir "$ARM"
log_step "[$genome_name / crisprp_only]"

# ── [1] Link CRISPR-P v2.0 raw results into 01_Raw_Input/ ────────────────────
CRISPRP_RAW="${CRISPRP_RAW_TEMPLATE/\{GROUP\}/$GENE_GROUP}"
CRISPRP_RAW="${CRISPRP_RAW/\{GENOME\}/$genome_name}"
CRISPRP_RAW="$PIPELINE_DIR/$CRISPRP_RAW"

STAGE1_DIR="$ARM/01_Raw_Input"
safe_mkdir "$STAGE1_DIR"

if [[ -d "$CRISPRP_RAW" ]]; then
    log_info "[1/8] Linking CRISPR-P v2.0 raw results from: $CRISPRP_RAW"
    for f in "$CRISPRP_RAW/"*.tsv "$CRISPRP_RAW/"*.csv; do
        [[ -f "$f" ]] || continue
        target="$STAGE1_DIR/$(basename "$f")"
        ln -sfn "$f" "$target"
    done
else
    log_error "[1/8] CRISPR-P v2.0 raw dir not found: $CRISPRP_RAW"
    log_error "      v3 requires manually-generated CRISPR-P v2.0 results. Place them at the"
    log_error "      path above (or override crispr_v3.crispr_p.raw_dir_template) and re-run."
    continue
fi

# ── [1b]–[9] downstream stages ───────────────────────────────────────────────
log_info "[1b-9] Running plant-only downstream stages"
run_plant_stages "$ARM" "$genome_name"
log_info "[$genome_name] crisprp_only complete"

# ── Manual-action reminders ──────────────────────────────────────────────────
log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "MANUAL ACTIONS before downstream interpretation:"
log_info ""
log_info "  [A] Verify CRISPR-P v2.0 raw inputs were current (rerun web tool if stale):"
log_info "      $CRISPRP_RAW"
log_info ""
log_info "  [B] AlphaFold3 (optional — only if higher-confidence structures needed)"
log_info "      ESMFold local PDBs are now produced automatically by stage [6b]:"
log_info "         $ARM/06_Protein_Biopython_Pfam/esmfold_structures/"
log_info "      For AlphaFold3 (web-only), submit FASTAs from"
log_info "         $ARM/06_Protein_Biopython_Pfam/"
log_info "      to https://alphafoldserver.com and place PDBs in:"
log_info "         $ARM/06_Protein_Biopython_Pfam/alphafold3_structures/"
log_info ""
log_info "  [C] Review ranked-guide table before experimental validation:"
log_info "      $ARM/08_Ranking_Composite/"
log_info "      Cross-check top guides against published SmelDMP / AtDMP off-target data."
log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

done  # end genome_name loop

log_step "CRISPR v3 Pipeline complete: $GENE_GROUP  (mode=$MODE)"
log_info "Root output: $CRISPR_DIR"

done  # end GENE_GROUP loop
