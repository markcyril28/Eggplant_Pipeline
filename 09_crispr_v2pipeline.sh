#!/bin/bash
# ============================================================================
# Program 9 v2: CRISPR KO Prediction Pipeline
# ============================================================================
# Eight-stage pipeline for guide RNA knockout efficacy prediction.
# Two mutually exclusive modes produce separate output folder trees:
#
#   mode = crispor_only
#     Runs all 8 stages on CRISPOR-designed guides only.
#     Output tree:
#       09_CRISPR_v2/{genome}/crispor_only/
#         01_Design_CRISPOR/               ← CRISPOR
#           logs/        ← CRISPOR stdout/stderr per gene
#           offtargets/  ← CRISPOR --offtargets table
#         01b_Plant_Filter/                ← TTTT / GC / 5' G|A sgRNA pre-filter
#         02_Rescore_crisprOn_DeepSpCas9/
#           non_plant_scorers/  ← mammalian/human models: CRISPRon, DeepSpCas9
#           plant_scorers/      ← plant-trained models: DeepCRISPR, CRISPR-Local
#                                  (populated only when [crispr_v2.plant_scorer]
#                                   .enabled = true AND plant_predictors non-empty)
#         03_OffTargets_CFD/               ← CFD score + paralog filter
#         04_Indels_inDelphi_Lindel/       ← inDelphi, Lindel
#         05_Transcripts_Biopython/        ← Biopython / pyfaidx
#         06_Protein_Biopython_Pfam/       ← Biopython + Pfam domain TSV
#         07_NMD_50ntRule/                 ← 50-nt EJC + plant ≥350-nt 3'UTR
#         08_Ranking_Composite/            ← weighted composite KO score
#           per_gene/    ← per-gene ranked tables and plots
#           ranked_guides.tsv, top10_guides.tsv, top10_ko_score.png (aggregate)
#
#   mode = comparison
#     Runs all 8 stages in two parallel arms, then compares them.
#     Output tree:
#       09_CRISPR_v2/{genome}/crispor/          ← CRISPOR arm (stages 1-8)
#       09_CRISPR_v2/{genome}/crispr_p_v2/      ← CRISPR-P v2.0 arm (stages 2-8;
#                                                   stage 1 = raw v1 input)
#       09_CRISPR_v2/{genome}/comparison/       ← compare_tools.py report
#
# Edit 09_crispr_v2CONFIG.toml then run:
#   bash i_crispr_v2_pipeline.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# All settings loaded from 09_crispr_v2CONFIG.toml — edit there.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"
V2_MOD="$MODULES/09_crispr_analysis/v2"

PROJECT_ROOT="$PIPELINE_DIR"

# WSL2/NTFS workaround: mkdir -p uses lstat internally and fails with EEXIST
# on Windows reparse-point directories even though the path is a valid dir.
# safe_mkdir falls back to [[ -d ]] (which follows links) before erroring.
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

# ─── Load GENE_GROUPS from shared config ─────────────────────────────────────
SHARED_CONFIG="$PIPELINE_DIR/09_crispr_v2CONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 09_crispr_v2CONFIG.toml" >&2
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

# ensure_indelphi_env ENV_NAME
#   Creates the legacy conda env for inDelphi (scikit-learn=0.20.0) if absent.
#   No-op when ENV_NAME is empty.
#   IMPORTANT: pin scikit-learn=0.20.0 exactly — 0.20.x patch releases (e.g.
#   0.20.4) fail inDelphi's hard version assertion in init_model().
ensure_indelphi_env() {
    local env_name="$1"
    [[ -z "$env_name" ]] && return 0

    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env_name"; then
        log_info "[inDelphi] Legacy env '$env_name' already exists."
        return 0
    fi

    log_info "[inDelphi] Creating legacy conda env '$env_name' (scikit-learn=0.20.0, Python 3.7)..."
    if ! conda create -y -n "$env_name" python=3.7 scikit-learn=0.20.0 numpy pandas; then
        log_warn "[inDelphi] 'conda create' failed for '$env_name'. inDelphi will be skipped."
        return 0
    fi

    # inDelphi has no PyPI/setup.py release. The bundled copy at
    # modules/09_crispr_analysis/v2/tools/inDelphi/ is added to PYTHONPATH at
    # runtime by 04_predict_indels.py, so no pip install step is needed here.
    log_info "[inDelphi] Legacy env '$env_name' ready (inDelphi loaded from bundled tools dir)."
}

# preflight_check_indelphi_sklearn ENV_NAME
#   Validates that the legacy conda env has a scikit-learn version inDelphi
#   supports (0.18.1 or 0.20.0 exactly). Logs an actionable ERROR if not,
#   so the run fails fast rather than silently degrading at Stage 04.
#   Non-fatal: the pipeline continues and Stage 04 falls back to Lindel-only.
preflight_check_indelphi_sklearn() {
    local env_name="$1"
    [[ -z "$env_name" ]] && return 0
    local actual_ver
    actual_ver=$(conda run -n "$env_name" python3 -c \
        "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo "UNAVAILABLE")
    if [[ "$actual_ver" == "0.18.1" || "$actual_ver" == "0.20.0" ]]; then
        log_info "[pre-flight] inDelphi legacy env '$env_name': scikit-learn $actual_ver — OK"
    else
        log_error "[pre-flight] inDelphi legacy env '$env_name' has scikit-learn $actual_ver — INCOMPATIBLE (requires 0.18.1 or 0.20.0 exactly). Stage 04 will fall back to Lindel-only. Fix with: conda install -n $env_name scikit-learn=0.20.0 --yes"
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
# run_8_stages ARM_DIR GENOME_NAME
#   Runs stages [2]–[8] on an arm root directory.
#
# Args:
#   ARM_DIR     — root of the arm (e.g. .../crispor_only or .../crispr_p_v2)
#   GENOME_NAME — genome subdirectory name (for logging)
#
# Required outer-scope variables (set before calling):
#   CONDA_RUN, MAX_PARALLEL, OVERWRITE
#   PREDICTORS, FLAG_THRESH                         (stage 2)
#   PARALOG_PATTERNS, PARALOG_ID_ARGS,
#     PARALOG_HIT_THRESH, CFD_SUM_THRESH            (stage 3)
#   INDEL_PREDICTORS, INDELPHI_CELL,
#     FS_THRESH, TOP_OUTCOMES                       (stage 4)
#   GENOME, GENE_GROUP, GTF_FOR_STAGE, TOP_INDELS   (stage 5)
#   DOMAIN_TSV, FLAG_DOMAINS                        (stage 6)
#   GTF, PTC_DIST                                   (stage 7)
#   WEIGHTS_JSON, TOP_N, OUT_FMT,
#     REPORT_DPI, REPORT_FMT                        (stage 8)
# ─────────────────────────────────────────────────────────────────────────────
run_8_stages() {
    local ARM_DIR="$1"
    local GENOME_NAME="$2"
    # Optional third arg overrides MAX_PARALLEL for this arm.
    # Comparison mode passes MAX_PARALLEL/2 so two concurrent arms stay within budget.
    local _MP="${3:-$MAX_PARALLEL}"

    # ── [1b] Plant sgRNA pre-filter (TTTT, GC, 5' G/A) ────────────────────────
    # Applied BEFORE rescoring so the on-target ranking reflects only guides
    # that can be expressed in a plant U6/U3 vector. Downstream stages prefer
    # its *.plant_filtered.tsv output over raw *.filtered.tsv when present.
    if op_enabled "plant_filter"; then
        local STAGE_DIR="$ARM_DIR/01b_Plant_Filter"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        local -a _S1_TSVS=()
        for tsv in "$ARM_DIR/01_"*/*.filtered.tsv; do
            [[ -f "$tsv" ]] && _S1_TSVS+=("$tsv")
        done
        if [[ ${#_S1_TSVS[@]} -eq 0 ]]; then
            log_warn "[01b_plant] No *.filtered.tsv under $ARM_DIR/01_* — plant filter skipped."
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
                --enabled                "$PLANT_FILTER_ENABLED" \
                "$OVERWRITE_FLAG" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [2] Rescore ───────────────────────────────────────────────────────────
    # Output layout (separates training-species provenance):
    #   02_Rescore_crisprOn_DeepSpCas9/
    #     non_plant_scorers/   ← mammalian/human-trained: CRISPRon, DeepSpCas9
    #     plant_scorers/       ← plant-trained: DeepCRISPR, CRISPR-Local
    #                            (populated only when [crispr_v2.plant_scorer]
    #                             .enabled = true AND its predictors list is
    #                             non-empty; otherwise created but left empty)
    # Stage [3] globs non_plant_scorers/ first, then plant_scorers/, then the
    # legacy top-level directory for backwards compatibility.
    if op_enabled "rescore_ontarget"; then
        local STAGE_DIR="$ARM_DIR/02_Rescore_crisprOn_DeepSpCas9"
        local NP_DIR="$STAGE_DIR/non_plant_scorers"
        local PL_DIR="$STAGE_DIR/plant_scorers"
        safe_mkdir "$STAGE_DIR" "$NP_DIR" "$PL_DIR"
        local PIDS=()
        # Stage-2 input preference order:
        #   1. *.plant_filtered.tsv  (stage 1b output — plant-biology filter applied)
        #   2. *.filtered.tsv        (stage 1 output — raw CRISPOR min-score filter)
        # Legacy *_filtered.tsv form also accepted for partially-run prior outputs.
        # IMPORTANT: do NOT fall back to *.tsv — 01_Design_CRISPOR/ also
        # contains *.design.tsv (raw, unfiltered CRISPOR output). Feeding that
        # into stage 2 would silently bypass the min-score filter.
        local -a STAGE1_TSVS=()
        for tsv in "$ARM_DIR/01b_Plant_Filter/"*.plant_filtered.tsv; do
            [[ -f "$tsv" ]] && STAGE1_TSVS+=("$tsv")
        done
        if [[ ${#STAGE1_TSVS[@]} -eq 0 ]]; then
            for tsv in "$ARM_DIR/01_"*/*.filtered.tsv "$ARM_DIR/01_"*/*_filtered.tsv; do
                [[ -f "$tsv" ]] && STAGE1_TSVS+=("$tsv")
            done
        fi
        if [[ ${#STAGE1_TSVS[@]} -eq 0 ]]; then
            log_warn "[02_rescore] No plant_filtered/filtered TSV found under $ARM_DIR — stage 2 skipped for this arm."
        fi
        for tsv in "${STAGE1_TSVS[@]}"; do
            wait_for_slot "$_MP"
            # Resolve the ORIGINAL target FASTA that CRISPOR received in stage
            # [1] for this gene stem. DeepSpCas9 needs the 4 nt upstream + 3
            # nt downstream flank that CRISPOR's 23-nt targetSeq does not
            # emit — see 02_rescore_ontarget.py::build_deepspcas9_context.
            # We glob across TARGET_FASTAS rather than trying to regenerate
            # the path so configs using `target_fastas` (absolute list) and
            # `target_fasta_template` both work.
            local _stem; _stem="$(basename "$tsv")"; _stem="${_stem%%.*}"
            local _target_fa=""
            for _cand in "${TARGET_FASTAS[@]:-}"; do
                local _full="$PIPELINE_DIR/$_cand"
                local _cand_stem; _cand_stem="$(basename "$_cand")"
                _cand_stem="${_cand_stem%.fasta}"; _cand_stem="${_cand_stem%.fa}"
                _cand_stem="${_cand_stem%.fna}"; _cand_stem="${_cand_stem%.fas}"
                if [[ "$_cand_stem" == "$_stem" && -f "$_full" ]]; then
                    _target_fa="$_full"
                    break
                fi
            done
            # Mammalian/human-trained rescorers → non_plant_scorers/
            local -a _S2_ARGS=(
                --input          "$tsv"
                --outdir         "$NP_DIR"
                --predictors     "${PREDICTORS[@]}"
                --flag-threshold "$FLAG_THRESH"
                --workers        "$INNER_WORKERS"
                "$OVERWRITE_FLAG"
            )
            [[ -n "$_target_fa" ]] && _S2_ARGS+=(--target-fasta "$_target_fa")
            $CONDA_RUN python3 "$V2_MOD/02_rescore_ontarget.py" "${_S2_ARGS[@]}" &
            PIDS+=("$!")

            # Plant-trained rescorers → plant_scorers/
            # Only runs when the user opts in via [crispr_v2.plant_scorer] and
            # lists ≥1 plant predictor. The module's plant-scorer dispatch is
            # a placeholder until DeepCRISPR / CRISPR-Local wrappers are
            # bundled — see [crispr_v2.plant_scorer] notes in the TOML.
            if [[ "${PLANT_SCORER_ENABLED:-false}" == "true" && ${#PLANT_PREDICTORS[@]} -gt 0 ]]; then
                wait_for_slot "$_MP"
                local -a _S2_PLANT_ARGS=(
                    --input          "$tsv"
                    --outdir         "$PL_DIR"
                    --predictors     "${PLANT_PREDICTORS[@]}"
                    --flag-threshold "$FLAG_THRESH"
                    --workers        "$INNER_WORKERS"
                    "$OVERWRITE_FLAG"
                )
                [[ -n "$_target_fa" ]] && _S2_PLANT_ARGS+=(--target-fasta "$_target_fa")
                $CONDA_RUN python3 "$V2_MOD/02_rescore_ontarget.py" "${_S2_PLANT_ARGS[@]}" &
                PIDS+=("$!")
            fi
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [3] Curate off-targets ────────────────────────────────────────────────
    if op_enabled "curate_offtargets"; then
        local STAGE_DIR="$ARM_DIR/03_OffTargets_CFD"
        # Find the off-target TSV in whichever stage-1 directory exists
        # (01_Design_CRISPOR for CRISPOR arm, 01_Raw_Input for CRISPR-P arm).
        # Preferred location: 01_*/offtargets/offtargets.tsv (new layout);
        # legacy fallback: 01_*/offtargets.tsv and variants.
        local OFFTARGET_TSV=""
        for _ot_candidate in \
                "$ARM_DIR/01_"*/offtargets/offtargets.tsv \
                "$ARM_DIR/01_"*/offtargets.tsv \
                "$ARM_DIR/01_"*/off_targets.tsv \
                "$ARM_DIR/01_"*/*offtarget*.tsv; do
            [[ -f "$_ot_candidate" ]] && { OFFTARGET_TSV="$_ot_candidate"; break; }
        done
        [[ -z "$OFFTARGET_TSV" ]] && OFFTARGET_TSV="$ARM_DIR/01_Design_CRISPOR/offtargets/offtargets.tsv"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        # Stage-[2] output lives under non_plant_scorers/ (mammalian models)
        # and optionally plant_scorers/ (when [crispr_v2.plant_scorer].enabled
        # is true). Read from both, with the legacy top-level path kept as a
        # third fallback so pre-reorg outputs still flow through stage [3].
        local -a _S2_OUT_TSVS=()
        for tsv in "$ARM_DIR/02_Rescore_crisprOn_DeepSpCas9/non_plant_scorers/"*.rescored.tsv \
                   "$ARM_DIR/02_Rescore_crisprOn_DeepSpCas9/plant_scorers/"*.rescored.tsv \
                   "$ARM_DIR/02_Rescore_crisprOn_DeepSpCas9/"*.rescored.tsv; do
            [[ -f "$tsv" ]] && _S2_OUT_TSVS+=("$tsv")
        done
        # De-duplicate by gene stem: if both non-plant and plant outputs exist
        # for the same gene, prefer non-plant (listed first above) — the
        # composite KO score in stage [8] already merges all scorer columns
        # via its mean-of-all-scorers aggregation.
        local -A _seen_stems=()
        local -a _S2_DEDUP=()
        for _p in "${_S2_OUT_TSVS[@]:-}"; do
            local _st; _st="$(basename "$_p")"; _st="${_st%%.*}"
            if [[ -z "${_seen_stems[$_st]:-}" ]]; then
                _seen_stems[$_st]=1
                _S2_DEDUP+=("$_p")
            fi
        done
        for tsv in "${_S2_DEDUP[@]:-}"; do
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

    # ── [4] Indel prediction ──────────────────────────────────────────────────
    if op_enabled "predict_indels"; then
        local STAGE_DIR="$ARM_DIR/04_Indels_inDelphi_Lindel"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/03_OffTargets_CFD/"*.curated.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            $CONDA_RUN python3 "$V2_MOD/04_predict_indels.py" \
                --input                  "$tsv" \
                --outdir                 "$STAGE_DIR" \
                --predictors             "${INDEL_PREDICTORS[@]}" \
                --indelphi-cell-type     "$INDELPHI_CELL" \
                --frameshift-threshold   "$FS_THRESH" \
                --top-outcomes           "$TOP_OUTCOMES" \
                --indelphi-conda-env     "$CONDA_ENV_INDELPHI" \
                --workers                "$INNER_WORKERS" \
                "$OVERWRITE_FLAG" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
    fi

    # ── [5] Rebuild transcripts ───────────────────────────────────────────────
    if op_enabled "rebuild_transcripts"; then
        local STAGE_DIR="$ARM_DIR/05_Transcripts_Biopython"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/04_Indels_inDelphi_Lindel/"*.indels.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            # Derive the target gene id from the filename stem (e.g.
            # SMEL5_01g008730.indels.tsv -> SMEL5_01g008730) so load_cds_coords
            # can match GTF gene_id exactly. Without this, the module falls back
            # to a family-name regex that never matches SMEL5/GPE001970 ids.
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

    # ── [6] Protein consequence ───────────────────────────────────────────────
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

    # ── [7] NMD prediction ────────────────────────────────────────────────────
    if op_enabled "predict_nmd"; then
        local STAGE_DIR="$ARM_DIR/07_NMD_50ntRule"
        safe_mkdir "$STAGE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/06_Protein_Biopython_Pfam/"*.protein.tsv; do
            [[ -f "$tsv" ]] || continue
            wait_for_slot "$_MP"
            # Derive target gene id from filename stem (see stage [5] comment).
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

    # ── [8] Rank guides ───────────────────────────────────────────────────────
    # Per-gene tables + plots go under per_gene/ to keep the stage root clean;
    # a follow-up aggregate call emits ranked_guides / top_guides / plot at
    # the stage root by concatenating the per-gene ranked TSVs.
    if op_enabled "rank_guides"; then
        local STAGE_DIR="$ARM_DIR/08_Ranking_Composite"
        local PER_GENE_DIR="$STAGE_DIR/per_gene"
        safe_mkdir "$PER_GENE_DIR"
        local PIDS=()
        for tsv in "$ARM_DIR/07_NMD_50ntRule/"*.nmd.tsv; do
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

        # Aggregate across all per-gene ranked tables -> stage root
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

    # ── [9] Guide comparison scatter (v1-style on-target vs off-target plot) ──
    # Aggregates all per-gene ranked tables from [8] into a single scatter
    # (on-target score x, off-target count y, gene-coloured, tier-shaped).
    # Runs once per arm — not per gene — so no parallel loop here.
    if op_enabled "comparison_scatter"; then
        local STAGE_DIR="$ARM_DIR/09_Guide_Scatter"
        local RANKED_DIR="$ARM_DIR/08_Ranking_Composite/per_gene"
        safe_mkdir "$STAGE_DIR"
        # Stage-8 emits {gene}.ranked.tsv; accept the legacy _ranked_guides.tsv
        # form too so older runs still plot cleanly.
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
            log_warn "[09_scatter] No ranked TSVs (*.ranked.tsv / *_ranked_guides.tsv) in $RANKED_DIR — skipping."
        fi
    fi
}

# =============================================================================
# Logging — set up once before the loop; teardown in EXIT trap
# =============================================================================
setup_logging

# =============================================================================
# One-time setup: ensure inDelphi legacy env exists before any work starts
# =============================================================================
_INDELPHI_ENV_GLOBAL=$(python3 "$TOML_PARSER" "$SHARED_CONFIG" crispr_v2 conda_env_indelphi 2>/dev/null || echo "")
ensure_indelphi_env "$_INDELPHI_ENV_GLOBAL" || true
preflight_check_indelphi_sklearn "$_INDELPHI_ENV_GLOBAL"

# =============================================================================
# Per gene-group loop
# =============================================================================
for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# ── Config resolution ─────────────────────────────────────────────────────────
# Use merge_toml.py (deep merge) rather than cat to avoid duplicate TOML table
# errors when both the shared config and the group override define sub-tables
# under the same [crispr_v2.*] parent (e.g. [crispr_v2.design_score]).
MERGE_TOML="$MODULES/utils/merge_toml.py"
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_crispr_v2_cfg_XXXXXX.toml")
    python3 "$MERGE_TOML" \
        "$PIPELINE_DIR/09_crispr_v2CONFIG.toml" \
        "$CONFIG_DIR/00_common_${GENE_GROUP}.toml" \
        "$CONFIG_DIR/09_crispr_analysis_v2_${GENE_GROUP}.toml" \
        > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

# ── Compute settings ──────────────────────────────────────────────────────────
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "4")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
# Per-process worker count for intra-module parallelism (row-level threading
# inside each Python module). Formula: total_cores / concurrent_gene_processes.
# Ensures the product MAX_PARALLEL * INNER_WORKERS stays at or below CPU.
INNER_WORKERS=$(( CPU / MAX_PARALLEL ))
(( INNER_WORKERS < 1 )) && INNER_WORKERS=1
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")
[[ "$OVERWRITE" == "true" || "$OVERWRITE" == "True" ]] && OVERWRITE_FLAG="--overwrite" || OVERWRITE_FLAG="--no-overwrite"

_base_dir_rel=$(get_toml general base_dir 2>/dev/null || echo "")
if [[ -z "$_base_dir_rel" ]]; then
    _base_dir_rel="III_RESULT/$GENE_GROUP"
    log_info "general.base_dir not in config — defaulting to $_base_dir_rel"
fi
BASE_DIR="$PIPELINE_DIR/$_base_dir_rel"

# Default genome FASTA / GTF — resolved via the keys declared in
# [crispr_v2.rebuild_transcripts] so a single TOML edit switches genomes.
_default_genome_key=$(get_toml crispr_v2 rebuild_transcripts genome_key 2>/dev/null || echo "")
_default_gtf_key=$(get_toml crispr_v2 rebuild_transcripts gtf_key 2>/dev/null || echo "")
_default_genome_rel=$( [[ -n "$_default_genome_key" ]] && get_toml reference "$_default_genome_key" 2>/dev/null || echo "")
DEFAULT_GENOME="${_default_genome_rel:+$PIPELINE_DIR/$_default_genome_rel}"
DEFAULT_GTF_REL=$( [[ -n "$_default_gtf_key" ]] && get_toml reference "$_default_gtf_key" 2>/dev/null || echo "")
DEFAULT_GTF="${DEFAULT_GTF_REL:+$PIPELINE_DIR/$DEFAULT_GTF_REL}"

DOMAIN_TSV_KEY=$(get_toml crispr_v2 protein_consequence domain_tsv_key 2>/dev/null || echo "")
DOMAIN_TSV=""
if [[ -n "$DOMAIN_TSV_KEY" ]]; then
    DOMAIN_TSV_REL=$(get_toml reference "$DOMAIN_TSV_KEY" 2>/dev/null || echo "")
    [[ -n "$DOMAIN_TSV_REL" ]] && DOMAIN_TSV="$PIPELINE_DIR/$DOMAIN_TSV_REL"
fi

CONDA_ENV=$(get_toml crispr_v2 conda_env 2>/dev/null || echo "crispr_v2")
CONDA_RUN="conda run -n $CONDA_ENV --no-capture-output"
CONDA_ENV_INDELPHI=$(get_toml crispr_v2 conda_env_indelphi 2>/dev/null || echo "")

ENABLED=$(get_toml crispr_v2 enabled 2>/dev/null || echo "false")
if [[ "$ENABLED" != "True" && "$ENABLED" != "true" ]]; then
    log_info "CRISPR v2 pipeline not enabled for $GENE_GROUP. Skipping."
    continue
fi

# ── Mode ──────────────────────────────────────────────────────────────────────
MODE=$(get_toml crispr_v2 mode 2>/dev/null || echo "crispor_only")
if [[ "$MODE" != "crispor_only" && "$MODE" != "comparison" ]]; then
    log_warn "Unknown mode '$MODE'; falling back to crispor_only."
    MODE="crispor_only"
fi

mapfile -t OPERATIONS < <(
    get_toml crispr_v2 operations 2>/dev/null || \
    printf '%s\n' design_score plant_filter rescore_ontarget curate_offtargets \
                   predict_indels rebuild_transcripts protein_consequence \
                   predict_nmd rank_guides
)

CRISPR_DIR="$BASE_DIR/09_CRISPR_v2"
safe_mkdir "$CRISPR_DIR"

mapfile -t GENOME_NAMES < <(
    get_toml crispr_v2 genomes names 2>/dev/null || echo "GPE001970_SMEL5"
)

log_step "CRISPR v2 KO Prediction Pipeline: $GENE_GROUP"
log_info "Mode:          $MODE"
log_info "Operations:    ${OPERATIONS[*]}"
log_info "Conda env:     $CONDA_ENV"
log_info "CPU / MAX_PARALLEL / INNER_WORKERS: $CPU / $MAX_PARALLEL / $INNER_WORKERS"

# Report whether the Azimuth/Fusi scorer is patched for sklearn >= 1.0.
# The patch (try/except fallback to 0.0) is applied by setup_conda_crispr_v2.sh.
_EFFSCORE="$V2_MOD/tools/crispor/crisporEffScores.py"
if [[ -f "$_EFFSCORE" ]]; then
    if grep -q "Azimuth/Fusi score skipped" "$_EFFSCORE"; then
        # Permanent expected state (sklearn >=1.0 incompatible with bundled
        # Azimuth pickles); log at INFO so it doesn't pollute error/warn logs.
        log_info "Azimuth/Fusi (fusi) scoring is DISABLED: crisporEffScores.py is patched for scikit-learn >= 1.0 compatibility. Fusi scores will be 0; all other CRISPOR scores are unaffected."
    else
        log_info "Azimuth/Fusi scoring: enabled (sklearn version compatible)."
    fi
fi
unset _EFFSCORE

# ── Shared stage config (loaded once, used by run_8_stages) ──────────────────
# ── Plant sgRNA pre-filter (stage [1b]) ───────────────────────────────────────
# All knobs live under [crispr_v2.plant_filter]. Defaults here mirror the
# TOML defaults so the pipeline still behaves correctly for older configs.
PLANT_FILTER_ENABLED=$(get_toml crispr_v2 plant_filter enabled               2>/dev/null || echo "true")
PLANT_TERM_RUN=$(get_toml       crispr_v2 plant_filter termination_run_length 2>/dev/null || echo "4")
PLANT_TERM_ACTION=$(get_toml    crispr_v2 plant_filter termination_action     2>/dev/null || echo "reject")
PLANT_PROMOTER_TYPE=$(get_toml  crispr_v2 plant_filter promoter_type          2>/dev/null || echo "U6")
PLANT_PROMOTER_ACTION=$(get_toml crispr_v2 plant_filter promoter_action       2>/dev/null || echo "flag")
PLANT_GC_MIN=$(get_toml         crispr_v2 plant_filter gc_min                 2>/dev/null || echo "0.30")
PLANT_GC_MAX=$(get_toml         crispr_v2 plant_filter gc_max                 2>/dev/null || echo "0.70")

PREDICTORS_RAW=$(get_toml crispr_v2 rescore_ontarget predictors    2>/dev/null || echo "crisprOn DeepSpCas9")
FLAG_THRESH=$(get_toml crispr_v2 rescore_ontarget flag_threshold    2>/dev/null || echo "0.3")
read -ra PREDICTORS <<< "$PREDICTORS_RAW"

# ── Plant-trained on-target rescorers ([crispr_v2.plant_scorer]) ─────────────
# Kept as a separate array from PREDICTORS so the two groups land in
# distinct subdirectories under stage [2]. Disabled by default until
# DeepCRISPR / CRISPR-Local wrappers are bundled in modules/09_*/v2/tools/.
PLANT_SCORER_ENABLED=$(get_toml crispr_v2 plant_scorer enabled 2>/dev/null || echo "false")
mapfile -t PLANT_PREDICTORS < <(get_toml crispr_v2 plant_scorer plant_predictors 2>/dev/null || true)

PARALOG_PATTERNS_RAW=$(get_toml crispr_v2 curate_offtargets paralog_patterns 2>/dev/null || echo "DMP HAP2 GCS1")
PARALOG_HIT_THRESH=$(get_toml crispr_v2 curate_offtargets paralog_hit_threshold 2>/dev/null || echo "1")
CFD_SUM_THRESH=$(get_toml crispr_v2 curate_offtargets cfd_sum_threshold        2>/dev/null || echo "0.2")
# Max mismatches for CRISPOR's --mm in stage 1. Pulled from the curation section
# because it governs the off-target enumeration that both stages 1 and 3 share.
MAX_MISMATCHES=$(get_toml crispr_v2 curate_offtargets max_mismatches             2>/dev/null || echo "4")
read -ra PARALOG_PATTERNS <<< "$PARALOG_PATTERNS_RAW"
mapfile -t PARALOG_IDS < <(get_toml crispr_v2 curate_offtargets paralog_gene_ids 2>/dev/null || true)
PARALOG_ID_ARGS=()
[[ ${#PARALOG_IDS[@]} -gt 0 ]] && PARALOG_ID_ARGS=(--paralog-gene-ids "${PARALOG_IDS[@]}")

INDELPHI_CELL=$(get_toml crispr_v2 predict_indels indelphi_cell_type  2>/dev/null || echo "HEK293")
FS_THRESH=$(get_toml crispr_v2 predict_indels frameshift_threshold     2>/dev/null || echo "0.5")
TOP_OUTCOMES=$(get_toml crispr_v2 predict_indels top_outcomes          2>/dev/null || echo "5")
INDEL_PRED_RAW=$(get_toml crispr_v2 predict_indels predictors          2>/dev/null || echo "inDelphi Lindel")
read -ra INDEL_PREDICTORS <<< "$INDEL_PRED_RAW"

GTF_KEY=$(get_toml crispr_v2 rebuild_transcripts gtf_key 2>/dev/null || echo "")
GTF_FOR_STAGE_FROM_KEY=""
if [[ -n "$GTF_KEY" ]]; then
    GTF_REL2=$(get_toml reference "$GTF_KEY" 2>/dev/null || echo "")
    [[ -n "$GTF_REL2" ]] && GTF_FOR_STAGE_FROM_KEY="$PIPELINE_DIR/$GTF_REL2"
fi
# GTF_FOR_STAGE is finalised inside the genome loop (falls back to per-genome $GTF)
TOP_INDELS=$(get_toml crispr_v2 rebuild_transcripts top_indels 2>/dev/null || echo "3")

FLAG_DOMAINS=$(get_toml crispr_v2 protein_consequence flag_domain_hits 2>/dev/null || echo "true")

PTC_DIST=$(get_toml crispr_v2 predict_nmd ptc_distance_threshold 2>/dev/null || echo "50")
# Plant NMD long-3'UTR rule (Kerényi 2008) — 0 disables the rule entirely and
# reverts to mammalian-only NMD prediction. Default 350 nt matches the Kerényi
# threshold reported for A. thaliana.
LONG_3UTR_THRESH=$(get_toml crispr_v2 predict_nmd long_3utr_threshold 2>/dev/null || echo "350")

# Weights JSON for rank_guides — read each key individually via parse_toml.py
# (merge_toml.py expands inline-tables to [section.subsection] headers, so a
# regex on the merged config file can never match the inline form).
_w_ontarget=$(get_toml  crispr_v2 rank_guides weights w_ontarget   2>/dev/null || echo "0.30")
_w_frameshift=$(get_toml crispr_v2 rank_guides weights w_frameshift 2>/dev/null || echo "0.25")
_w_nmd=$(get_toml       crispr_v2 rank_guides weights w_nmd        2>/dev/null || echo "0.20")
_w_offtarget=$(get_toml crispr_v2 rank_guides weights w_offtarget  2>/dev/null || echo "0.15")
_w_domain=$(get_toml    crispr_v2 rank_guides weights w_domain     2>/dev/null || echo "0.10")
WEIGHTS_JSON=$(printf '{"w_ontarget":%s,"w_frameshift":%s,"w_nmd":%s,"w_offtarget":%s,"w_domain":%s}' \
    "$_w_ontarget" "$_w_frameshift" "$_w_nmd" "$_w_offtarget" "$_w_domain")

TOP_N=$(get_toml crispr_v2 rank_guides top_n       2>/dev/null || echo "10")
OUT_FMT=$(get_toml crispr_v2 rank_guides output_format 2>/dev/null || echo "tsv")
# Absolute CFD-sum cap for reproducible off-target penalty normalisation.
# Replaces the old per-batch max-normalisation so a guide's composite score
# no longer depends on which siblings were in the run.
CFD_SUM_CAP=$(get_toml crispr_v2 rank_guides cfd_sum_cap 2>/dev/null || echo "5.0")
REPORT_DPI=$(get_toml crispr_v2 report dpi          2>/dev/null || echo "600")
REPORT_FMT=$(get_toml crispr_v2 report format       2>/dev/null || echo "png")

# ─── Stage-9 guide comparison scatter ────────────────────────────────────────
# Every knob below maps 1:1 to a key in [crispr_v2.comparison_scatter].
# Top-N guides per gene to label + tabulate.
SCATTER_TOP_N=$(get_toml crispr_v2 comparison_scatter top_n 2>/dev/null || echo "3")
# Off-target y-axis upper bound (0 = auto-scale; any positive number clamps
# the axis and renders over-cap guides at the top edge with "(>cap)" label).
SCATTER_Y_MAX_CAP=$(get_toml crispr_v2 comparison_scatter y_max_cap 2>/dev/null || echo "0")
# Tier cut-offs on composite_ko_score — control marker shape (▲/●/■).
SCATTER_TIER_HIGH=$(get_toml crispr_v2 comparison_scatter tier_high      2>/dev/null || echo "0.7")
SCATTER_TIER_MOD=$(get_toml  crispr_v2 comparison_scatter tier_moderate  2>/dev/null || echo "0.5")
# Single-variant fallback (used when score_variants is empty).
SCATTER_X_COLUMN=$(get_toml   crispr_v2 comparison_scatter x_score_column  2>/dev/null || echo "Moreno-Mateos-Score")
SCATTER_X_HIGH=$(get_toml     crispr_v2 comparison_scatter x_axis_high     2>/dev/null || echo "40")
SCATTER_X_MOD=$(get_toml      crispr_v2 comparison_scatter x_axis_moderate 2>/dev/null || echo "30")
SCATTER_X_MAX=$(get_toml      crispr_v2 comparison_scatter x_axis_max      2>/dev/null || echo "100")
# Multi-variant scatter spec: one plot set per scoring column. Entries in the
# TOML array are already pipe-separated column:high:mod:max:suffix strings —
# pipe-join them again here so the whole spec goes through one CLI flag.
mapfile -t _SCATTER_VARIANTS < <(get_toml crispr_v2 comparison_scatter score_variants 2>/dev/null || true)
SCATTER_VARIANTS_ARG=""
if (( ${#_SCATTER_VARIANTS[@]} > 0 )); then
    SCATTER_VARIANTS_ARG=$(IFS='|'; echo "${_SCATTER_VARIANTS[*]}")
fi

# ─── CRISPOR stage-1 params ──────────────────────────────────────────────────
PAM=$(get_toml crispr_v2 design_score pam              2>/dev/null || echo "NGG")
CRISPOR_GENOME=$(get_toml crispr_v2 design_score crispor_genome 2>/dev/null || echo "melongena")
MIN_SCORE=$(get_toml crispr_v2 design_score min_score  2>/dev/null || echo "20")
BATCH_SIZE=$(get_toml crispr_v2 design_score batch_size 2>/dev/null || echo "100")

# Ordered list of score columns for the stage-1 filter. Loaded as a bash
# array, then pipe-joined when passed to 01_design_score.sh (column names
# contain spaces and quotes, so ' ' / ',' separators aren't safe).
mapfile -t SCORE_COLUMNS < <(get_toml crispr_v2 design_score score_columns 2>/dev/null || true)
if [[ ${#SCORE_COLUMNS[@]} -eq 0 ]]; then
    SCORE_COLUMNS=(
        "Moreno-Mateos-Score"
        "Doench '16-Score"
        "Doench-RuleSet3-Score"
        "DoenchScore"
        "on_target_score"
    )
fi
SCORE_COLUMNS_ARG=$(IFS='|'; echo "${SCORE_COLUMNS[*]}")

# Target nucleotide FASTAs for CRISPOR de-novo guide design.
# Resolution (in order of precedence):
#   1. target_fastas   (explicit list) — used as-is for all genomes
#   2. target_fasta_template — {GENOME} substituted per genome inside the loop
#   3. legacy grna_fastas — deprecated fallback (logs a warning)
mapfile -t TARGET_FASTAS_EXPLICIT < <(get_toml crispr_v2 design_score target_fastas 2>/dev/null || true)
TARGET_FASTA_TEMPLATE=$(get_toml crispr_v2 design_score target_fasta_template 2>/dev/null || echo "")
if [[ ${#TARGET_FASTAS_EXPLICIT[@]} -eq 0 && -z "$TARGET_FASTA_TEMPLATE" ]]; then
    mapfile -t TARGET_FASTAS_EXPLICIT < <(get_toml crispr_v2 design_score grna_fastas 2>/dev/null || true)
    if [[ ${#TARGET_FASTAS_EXPLICIT[@]} -gt 0 ]]; then
        log_warn "[design_score] Using legacy 'grna_fastas' key. Rename to 'target_fastas'"
        log_warn "[design_score] and point at gene/CDS nucleotide FASTAs, not pre-designed guides."
    fi
fi

# ─── CRISPR-P v2.0 raw dir template (comparison mode) ───────────────────────
CRISPRP_RAW_TEMPLATE=$(get_toml crispr_v2 crispr_p raw_dir_template 2>/dev/null || \
    echo "III_RESULT/{GROUP}/09_CRISPR_Off-Target_Analysis/{GENOME}/01_Raw_Scoring_Results_from_CRISPR-P_V2_0")

# ─── Comparison config ───────────────────────────────────────────────────────
COMP_METRICS_RAW=$(get_toml crispr_v2 comparison metrics 2>/dev/null || \
    echo "doench2016OnTarget crisprOn_score best_fs_frac cfd_sum composite_ko_score")
read -ra COMP_METRICS <<< "$COMP_METRICS_RAW"
JACCARD_TOP_N=$(get_toml crispr_v2 comparison jaccard_top_n 2>/dev/null || echo "10")
COMP_DPI=$(get_toml crispr_v2 comparison dpi    2>/dev/null || echo "$REPORT_DPI")
COMP_FMT=$(get_toml crispr_v2 comparison format 2>/dev/null || echo "$REPORT_FMT")

# =============================================================================
# MODE ROUTING
# =============================================================================

for genome_name in "${GENOME_NAMES[@]}"; do

GENOME_BASE="$CRISPR_DIR/${genome_name}"

# ── Per-genome FASTA / GTF resolution ────────────────────────────────────────
# Try TOML key reference.<genome_name>_genome; fall back to DEFAULT_GENOME.
_genome_key=$(echo "${genome_name}" | tr '.' '_')  # dots → underscores for TOML keys
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
# Finalise GTF_FOR_STAGE: prefer the explicit gtf_key, then fall back to per-genome GTF
GTF_FOR_STAGE="${GTF_FOR_STAGE_FROM_KEY:-$GTF}"

# ── Per-genome TARGET_FASTAS resolution ──────────────────────────────────────
# If an explicit list was provided, use it verbatim for every genome.
# Otherwise, substitute {GENOME} in the template.
TARGET_FASTAS=()
if [[ ${#TARGET_FASTAS_EXPLICIT[@]} -gt 0 ]]; then
    TARGET_FASTAS=("${TARGET_FASTAS_EXPLICIT[@]}")
elif [[ -n "$TARGET_FASTA_TEMPLATE" ]]; then
    TARGET_FASTAS=("${TARGET_FASTA_TEMPLATE//\{GENOME\}/$genome_name}")
fi

# ─────────────────────────────────────────────────────────────────────────────
# MODE: crispor_only
# ─────────────────────────────────────────────────────────────────────────────
if [[ "$MODE" == "crispor_only" ]]; then

    ARM="$GENOME_BASE/crispor_only"
    safe_mkdir "$ARM"
    log_step "[$genome_name / crispor_only]"

    # ── [1] CRISPOR design + scoring ─────────────────────────────────────────
    if op_enabled "design_score"; then
        log_info "[1/8] CRISPOR design + scoring"
        STAGE_DIR="$ARM/01_Design_CRISPOR"
        safe_mkdir "$STAGE_DIR"

        # Prepare CRISPOR genome dir (BWA index + 2bit) — idempotent.
        # Must complete before parallel CRISPOR calls.
        log_info "[1/8] Preparing CRISPOR genome: $CRISPOR_GENOME -> $GENOME"
        $CONDA_RUN bash "$V2_MOD/00_prepare_crispor_genome.sh" \
            --genome-fasta "$GENOME" \
            --name         "$CRISPOR_GENOME"

        # Per-job BWA thread count — saturate cores without oversubscribing when
        # MAX_PARALLEL CRISPOR calls run concurrently (each runs its own bwa aln).
        CRISPOR_THREADS=$(( CPU / MAX_PARALLEL ))
        (( CRISPOR_THREADS < 1 )) && CRISPOR_THREADS=1

        PIDS=()
        for target in "${TARGET_FASTAS[@]}"; do
            TARGET_FULL="$PIPELINE_DIR/$target"
            [[ -f "$TARGET_FULL" ]] || { log_warn "Target FASTA not found: $TARGET_FULL"; continue; }
            wait_for_slot "$MAX_PARALLEL"
            $CONDA_RUN bash "$V2_MOD/01_design_score.sh" \
                --grna-fasta      "$TARGET_FULL" \
                --genome-fasta    "$GENOME" \
                --outdir          "$STAGE_DIR" \
                --pam             "$PAM" \
                --crispor-genome  "$CRISPOR_GENOME" \
                --min-score       "$MIN_SCORE" \
                --score-columns   "$SCORE_COLUMNS_ARG" \
                --batch-size      "$BATCH_SIZE" \
                --threads         "$CRISPOR_THREADS" \
                --max-mismatches  "$MAX_MISMATCHES" \
                --overwrite       "$OVERWRITE" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
        log_info "[1/8] Done"
    fi

    # ── [2]–[8] ──────────────────────────────────────────────────────────────
    log_info "[2-8/8] Running downstream stages on crispor_only arm"
    run_8_stages "$ARM" "$genome_name"
    log_info "[$genome_name] crispor_only complete"

    # ── Manual action reminders ───────────────────────────────────────────────
    # Expected informational notices, not warnings — use log_info so they stay
    # out of logs/error_warn_logs/ on every run.
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_info "MANUAL ACTIONS REQUIRED before proceeding to downstream analysis:"
    log_info ""
    log_info "  [A] AlphaFold3 (web server — cannot be automated)"
    log_info "      1. Collect mutant protein FASTA sequences from:"
    log_info "         $ARM/06_Protein_Biopython_Pfam/"
    log_info "      2. Submit each to https://alphafoldserver.com"
    log_info "      3. Download predicted PDB structures and place them under:"
    log_info "         $ARM/06_Protein_Biopython_Pfam/alphafold3_structures/"
    log_info "      Note: free tier allows ~20 jobs/day; batch submissions may"
    log_info "      require an institutional account."
    log_info ""
    log_info "  [B] ESMFold (only if esmfold_backend = 'api' and API is rate-limited)"
    log_info "      The pipeline calls the public API automatically, but if you"
    log_info "      hit rate limits, resubmit failed sequences manually at:"
    log_info "      https://esmatlas.com/resources?action=fold"
    log_info "      and place the resulting PDBs in the same structures folder."
    log_info ""
    log_info "  [C] CRISPR-P v2.0 (comparison mode only — currently mode=$MODE)"
    log_info "      If switching to mode=comparison, run CRISPR-P v2.0 manually"
    log_info "      for $GENE_GROUP / $genome_name and place results under:"
    log_info "      III_RESULT/$GENE_GROUP/09_CRISPR_Off-Target_Analysis/$genome_name/"
    log_info "      01_Raw_Scoring_Results_from_CRISPR-P_V2_0/"
    log_info ""
    log_info "  [D] Review ranked guide table before experimental validation:"
    log_info "      $ARM/08_Ranking_Composite/"
    log_info "      Confirm top guides against published DMP off-target data."
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# ─────────────────────────────────────────────────────────────────────────────
# MODE: comparison
# ─────────────────────────────────────────────────────────────────────────────
elif [[ "$MODE" == "comparison" ]]; then

    ARM_C="$GENOME_BASE/crispor"          # CRISPOR arm
    ARM_P="$GENOME_BASE/crispr_p_v2"      # CRISPR-P v2.0 arm
    ARM_COMP="$GENOME_BASE/comparison"    # comparison report
    safe_mkdir "$ARM_C" "$ARM_P" "$ARM_COMP"

    log_step "[$genome_name / comparison]"

    # ── CRISPOR ARM [1]: design + scoring ────────────────────────────────────
    if op_enabled "design_score"; then
        log_info "[1/8] CRISPOR arm — design + scoring"
        STAGE_DIR="$ARM_C/01_Design_CRISPOR"
        safe_mkdir "$STAGE_DIR"

        # Prepare CRISPOR genome dir (BWA index + 2bit) — idempotent.
        log_info "[1/8] Preparing CRISPOR genome: $CRISPOR_GENOME -> $GENOME"
        $CONDA_RUN bash "$V2_MOD/00_prepare_crispor_genome.sh" \
            --genome-fasta "$GENOME" \
            --name         "$CRISPOR_GENOME"

        # Per-job BWA thread count — see crispor_only branch for rationale.
        CRISPOR_THREADS=$(( CPU / MAX_PARALLEL ))
        (( CRISPOR_THREADS < 1 )) && CRISPOR_THREADS=1

        PIDS=()
        for target in "${TARGET_FASTAS[@]}"; do
            TARGET_FULL="$PIPELINE_DIR/$target"
            [[ -f "$TARGET_FULL" ]] || { log_warn "Target FASTA not found: $TARGET_FULL"; continue; }
            wait_for_slot "$MAX_PARALLEL"
            $CONDA_RUN bash "$V2_MOD/01_design_score.sh" \
                --grna-fasta      "$TARGET_FULL" \
                --genome-fasta    "$GENOME" \
                --outdir          "$STAGE_DIR" \
                --pam             "$PAM" \
                --crispor-genome  "$CRISPOR_GENOME" \
                --min-score       "$MIN_SCORE" \
                --score-columns   "$SCORE_COLUMNS_ARG" \
                --batch-size      "$BATCH_SIZE" \
                --threads         "$CRISPOR_THREADS" \
                --max-mismatches  "$MAX_MISMATCHES" \
                --overwrite       "$OVERWRITE" &
            PIDS+=("$!")
        done
        for pid in "${PIDS[@]}"; do wait "$pid"; done
        log_info "[1/8] CRISPOR arm — design done"
    fi

    # ── CRISPR-P v2.0 ARM [1]: copy/link raw v1 results ─────────────────────
    # Replace {GROUP} and {GENOME} placeholders in the template path
    CRISPRP_RAW="${CRISPRP_RAW_TEMPLATE/\{GROUP\}/$GENE_GROUP}"
    CRISPRP_RAW="${CRISPRP_RAW/\{GENOME\}/$genome_name}"
    CRISPRP_RAW="$PIPELINE_DIR/$CRISPRP_RAW"

    CRISPRP_STAGE1="$ARM_P/01_Raw_Input"
    safe_mkdir "$CRISPRP_STAGE1"

    if [[ -d "$CRISPRP_RAW" ]]; then
        log_info "[1/8] CRISPR-P v2.0 arm — linking raw results from: $CRISPRP_RAW"
        # Materialize tab-delimited siblings of any *.csv exports so the
        # symlink loop below picks up TSVs first and downstream stages stay on
        # a single delimiter. Idempotent: only regenerates when missing/stale.
        python3 "$MODULES/09_crispr_analysis/csv_to_tsv.py" "$CRISPRP_RAW" \
            | while IFS= read -r line; do log_info "$line"; done
        # Symlink each TSV/CSV into the arm's stage-1 folder so downstream
        # stages can glob for files the same way as the CRISPOR arm.
        for f in "$CRISPRP_RAW/"*.tsv "$CRISPRP_RAW/"*.csv; do
            [[ -f "$f" ]] || continue
            target="$CRISPRP_STAGE1/$(basename "$f")"
            [[ -e "$target" ]] || ln -s "$f" "$target"
        done
    else
        log_warn "[1/8] CRISPR-P v2.0 raw dir not found: $CRISPRP_RAW"
        log_warn "      Stages [2]-[8] for the CRISPR-P arm will be skipped."
    fi

    # ── CRISPOR ARM [2]–[8] and CRISPR-P ARM [2]–[8] — run concurrently ──────
    # Each arm gets half the slot budget so combined parallelism stays within MAX_PARALLEL.
    _ARM_PARALLEL=$(( MAX_PARALLEL / 2 ))
    (( _ARM_PARALLEL < 1 )) && _ARM_PARALLEL=1

    log_info "[2-8/8] CRISPOR arm — downstream stages (concurrent, _MP=$_ARM_PARALLEL)"
    ( trap - EXIT; run_8_stages "$ARM_C" "$genome_name" "$_ARM_PARALLEL" ) &
    _ARM_C_PID=$!

    if [[ -d "$CRISPRP_STAGE1" ]] && compgen -G "$CRISPRP_STAGE1/"'*.tsv' > /dev/null 2>&1; then
        log_info "[2-8/8] CRISPR-P v2.0 arm — downstream stages (concurrent, _MP=$_ARM_PARALLEL)"
        ( trap - EXIT; run_8_stages "$ARM_P" "$genome_name" "$_ARM_PARALLEL" ) &
        _ARM_P_PID=$!
    else
        log_warn "CRISPR-P v2.0 arm skipped (no input TSVs found in $CRISPRP_STAGE1)"
        _ARM_P_PID=""
    fi

    wait "$_ARM_C_PID"
    [[ -n "$_ARM_P_PID" ]] && wait "$_ARM_P_PID"

    # ── [C] Comparison report ─────────────────────────────────────────────────
    if op_enabled "compare_tools"; then
        log_step "[C] Comparison report"

        # Find the ranked-guides TSV from each arm
        CRISPOR_RANKED=$(find "$ARM_C/08_Ranking_Composite" -name "ranked_guides.*" 2>/dev/null | head -1 || true)
        CRISPRP_RANKED=$(find "$ARM_P/08_Ranking_Composite" -name "ranked_guides.*" 2>/dev/null | head -1 || true)

        if [[ -z "$CRISPOR_RANKED" ]]; then
            log_warn "CRISPOR ranked table not found — skipping comparison."
        elif [[ -z "$CRISPRP_RANKED" ]]; then
            log_warn "CRISPR-P v2.0 ranked table not found — skipping comparison."
        else
            $CONDA_RUN python3 "$V2_MOD/compare_tools.py" \
                --crispor-ranked  "$CRISPOR_RANKED" \
                --crisprp-ranked  "$CRISPRP_RANKED" \
                --outdir          "$ARM_COMP" \
                --metrics         "${COMP_METRICS[@]}" \
                --jaccard-top-n   "$JACCARD_TOP_N" \
                --dpi             "$COMP_DPI" \
                --format          "$COMP_FMT" \
                "$OVERWRITE_FLAG"
            log_info "[C] Comparison report -> $ARM_COMP"
        fi
    fi

    log_info "[$genome_name] comparison complete"
    log_info "  CRISPOR arm:      $ARM_C"
    log_info "  CRISPR-P v2.0 arm: $ARM_P"
    log_info "  Comparison:       $ARM_COMP"

    # ── Manual action reminders ───────────────────────────────────────────────
    # Expected informational notices, not warnings — use log_info so they stay
    # out of logs/error_warn_logs/ on every run.
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_info "MANUAL ACTIONS REQUIRED before proceeding to downstream analysis:"
    log_info ""
    log_info "  [A] AlphaFold3 (web server — cannot be automated)"
    log_info "      1. Collect mutant protein FASTA sequences from:"
    log_info "         $ARM_C/06_Protein_Biopython_Pfam/"
    log_info "      2. Submit each to https://alphafoldserver.com"
    log_info "      3. Download predicted PDB structures and place them under:"
    log_info "         $ARM_C/06_Protein_Biopython_Pfam/alphafold3_structures/"
    log_info "      Note: free tier allows ~20 jobs/day; batch submissions may"
    log_info "      require an institutional account."
    log_info ""
    log_info "  [B] ESMFold (only if esmfold_backend = 'api' and API is rate-limited)"
    log_info "      If you hit rate limits, resubmit failed sequences manually at:"
    log_info "      https://esmatlas.com/resources?action=fold"
    log_info "      and place the resulting PDBs in the same structures folder."
    log_info ""
    log_info "  [C] Review comparison report and ranked guide tables:"
    log_info "      CRISPOR:    $ARM_C/08_Ranking_Composite/"
    log_info "      CRISPR-P:   $ARM_P/08_Ranking_Composite/"
    log_info "      Comparison: $ARM_COMP"
    log_info "      Confirm top guides against published DMP off-target data."
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

fi   # end mode branch

done  # end genome_name loop

log_step "CRISPR v2 Pipeline complete: $GENE_GROUP  (mode=$MODE)"
log_info "Root output: $CRISPR_DIR"

done  # end GENE_GROUP loop
# teardown_logging is called by the EXIT trap (line 78)
