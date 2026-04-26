#!/bin/bash
# ============================================================================
# Stage 05: Phylogenetic Analysis
# ============================================================================
# Comment in/out the gene groups below, then run:
#   bash 05_phylo.sh
# ============================================================================

set -euo pipefail

# Ensure conda environment is activated
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "egg" ]]; then
    echo "Activating conda environment: egg"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate egg
fi

# ===================== IMPORTANT VARIABLES =====================
# Comment in/out the gene groups to process:
GENE_GROUPS=(
    "DMP"
    #"GRF_GIF"
    #"PLA"
)

# CPU, MAX_PARALLEL, OVERWRITE, OPERATIONS, and PHYLO_SOFTWARE are loaded from
# 05_phyloCONFIG.toml — edit that file to change them.

# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

should_run() { [[ " ${OPERATIONS[@]} " =~ " $1 " ]]; }

load_operations_from_config() {
    local ops_str
    ops_str=$(get_toml pipeline operations 2>/dev/null || true)
    if [[ -n "$ops_str" ]]; then
        # parse_toml.py outputs one item per line; mapfile -t splits on newlines
        mapfile -t OPERATIONS <<< "$ops_str"
    else
        OPERATIONS=("build_tree" "visualize_tree" "compare_trees")
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

normalize_phylo_software() {
    local raw="$1"
    local normalized
    normalized=$(echo "$raw" | tr '[:lower:]' '[:upper:]' | tr '-' '_')

    case "$normalized" in
        MEGACC|MEGA_CC) echo "MEGA_CC" ;;
        IQTREE2|IQ_TREE2) echo "IQTREE2" ;;
        RAXML|RAXML_NG|RAXMLNG) echo "RAXML" ;;
        *) return 1 ;;
    esac
}

# Wait until a parallel slot opens before launching the next job
wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 0.5
    done
}

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# Resolve config: deep-merge shared defaults + group overrides
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
MERGE_TOML="$MODULES/utils/merge_toml.py"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_phylo_cfg_XXXXXX.toml")
    python3 "$MERGE_TOML" \
        "$PIPELINE_DIR/05_phyloCONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/05_phylogenetic_analysis.toml" > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

MACHINE=$(get_toml pipeline machine)
CPU=$(get_toml pipeline compute "$MACHINE" threads)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel)
OVERWRITE=$(get_toml pipeline overwrite)
load_operations_from_config

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"
MEGACC_CONFIG_MAO="$PIPELINE_DIR/$(get_toml phylogenetics megacc config_file 2>/dev/null || get_toml phylogenetics config_file 2>/dev/null || echo "")"
_megacc_nuc=$(get_toml phylogenetics megacc config_file_nucleotide 2>/dev/null || echo "")
_megacc_aa=$(get_toml phylogenetics megacc config_file_protein 2>/dev/null || echo "")
MEGACC_CONFIG_NUC="$( [[ -n "$_megacc_nuc" ]] && echo "$PIPELINE_DIR/$_megacc_nuc" || echo "" )"
MEGACC_CONFIG_AA="$( [[ -n "$_megacc_aa" ]] && echo "$PIPELINE_DIR/$_megacc_aa" || echo "" )"
unset _megacc_nuc _megacc_aa

mapfile -t _raw_software < <(get_toml phylogenetics software 2>/dev/null || true)
PHYLO_SOFTWARE=()
for sw in "${_raw_software[@]}"; do
    [[ -z "$sw" ]] && continue
    canonical=$(normalize_phylo_software "$sw") || {
        log_error "Unsupported phylogeny software: $sw (use: megacc, iqtree2, or raxml)"
        exit 1
    }
    PHYLO_SOFTWARE+=("$canonical")
done
unset _raw_software
if (( ${#PHYLO_SOFTWARE[@]} == 0 )); then
    log_error "No phylogenetics software configured. Set [phylogenetics].software in the TOML config."
    exit 1
fi

IQTREE2_BOOTSTRAP=$(get_toml phylogenetics iqtree2 bootstrap 2>/dev/null || get_toml phylogenetics bootstrap 2>/dev/null || echo "5000")
IQTREE2_ALRT=$(get_toml phylogenetics iqtree2 alrt 2>/dev/null || get_toml phylogenetics alrt 2>/dev/null || echo "5000")
IQTREE2_MODEL=$(get_toml phylogenetics iqtree2 model 2>/dev/null || echo "MFP+MERGE")
IQTREE2_MODEL_NUC=$(get_toml phylogenetics iqtree2 nucleotide_model 2>/dev/null || echo "")
IQTREE2_MODEL_AA=$(get_toml phylogenetics iqtree2 protein_model 2>/dev/null || echo "")
IQTREE2_ALLNNI=$(get_toml phylogenetics iqtree2 allnni 2>/dev/null || echo "true")
IQTREE2_POLYTOMY=$(get_toml phylogenetics iqtree2 polytomy 2>/dev/null || echo "false")
IQTREE2_SAFE=$(get_toml phylogenetics iqtree2 safe 2>/dev/null || echo "true")
IQTREE2_BNNI=$(get_toml phylogenetics iqtree2 bnni 2>/dev/null || echo "true")
IQTREE2_REDO=$(get_toml phylogenetics iqtree2 redo 2>/dev/null || echo "true")
IQTREE2_FAST=$(get_toml phylogenetics iqtree2 fast 2>/dev/null || echo "false")
IQTREE2_PERS=$(get_toml phylogenetics iqtree2 pers 2>/dev/null || echo "0.05")
RAXML_BINARY=$(get_toml phylogenetics raxml binary 2>/dev/null || echo "raxml-ng")
RAXML_MODEL_NUC=$(get_toml phylogenetics raxml nucleotide model 2>/dev/null || echo "GTR+FC+R4")
RAXML_MODEL_AA=$(get_toml phylogenetics raxml protein model 2>/dev/null || echo "LG+FC+R4")
RAXML_SEED=$(get_toml phylogenetics raxml seed 2>/dev/null || echo "12345")
RAXML_MODE=$(get_toml phylogenetics raxml mode 2>/dev/null || echo "all")
RAXML_BS_TREES=$(get_toml phylogenetics raxml bs_trees 2>/dev/null || echo "5000")
RAXML_SEARCH_REPLICATES=$(get_toml phylogenetics raxml search_replicates 2>/dev/null || echo "50")
RAXML_REDO=$(get_toml phylogenetics raxml redo 2>/dev/null || echo "true")

MSA_DIR_NAME=$(get_toml output_dirs msa 2>/dev/null || echo "04_MSA")
ALIGN_DIR="$BASE_DIR/$MSA_DIR_NAME"
PHYLO_INPUT_PATTERN=$(get_toml phylogenetics input_pattern 2>/dev/null || echo "*.fas")

# Read input_subdirs array; fall back to legacy input_subdir scalar if absent.
# parse_toml.py emits one item per line, so use mapfile (read -ra would only
# capture the first line and silently drop the rest).
mapfile -t PHYLO_INPUT_SUBDIRS < <(get_toml phylogenetics input_subdirs 2>/dev/null || true)
if (( ${#PHYLO_INPUT_SUBDIRS[@]} == 0 )); then
    mapfile -t PHYLO_INPUT_SUBDIRS < <(get_toml phylogenetics input_subdir 2>/dev/null || true)
fi

# Validate every configured subdir and build the list of real paths
ALIGN_INPUT_DIRS=()
for _subdir in "${PHYLO_INPUT_SUBDIRS[@]:-}"; do
    [[ -z "$_subdir" ]] && continue
    _dir="$ALIGN_DIR/$_subdir"
    if [[ -d "$_dir" ]]; then
        ALIGN_INPUT_DIRS+=("$_dir")
    else
        log_warn "Phylo input subdir not found, skipping: $_dir"
    fi
done
unset _subdir _dir

# Read input_files array (direct .fas file references, relative to ALIGN_DIR).
# Use mapfile because parse_toml.py emits one entry per line.
mapfile -t _file_list < <(get_toml phylogenetics input_files 2>/dev/null || true)
ALIGN_INPUT_FILES=()
for _f in "${_file_list[@]:-}"; do
    [[ -z "$_f" ]] && continue
    _fp="$ALIGN_DIR/$_f"
    if [[ -f "$_fp" ]]; then
        ALIGN_INPUT_FILES+=("$_fp")
    else
        log_warn "Phylo input file not found, skipping: $_fp"
    fi
done
unset _file_list _f _fp

PHYLO_DIR_NAME=$(get_toml output_dirs phylogenetics 2>/dev/null || echo "05_Phylogenetics")
PHYLO_DIR="$BASE_DIR/$PHYLO_DIR_NAME"
mkdir -p "$PHYLO_DIR"

if (( MAX_PARALLEL < 1 )); then
    log_error "pipeline.compute.${MACHINE}.max_parallel must be >= 1 (current: $MAX_PARALLEL)"
    exit 1
fi

setup_logging
log_step "Phylogenetic Analysis: $GENE_GROUP"
log_info "Operations: ${OPERATIONS[*]}"
log_info "Program-level parallelism: max ${MAX_PARALLEL} software runs at once"

# ======================== Build Trees ========================
if should_run "build_tree"; then
    if (( ${#ALIGN_INPUT_DIRS[@]} == 0 && ${#ALIGN_INPUT_FILES[@]} == 0 )); then
        log_error "No valid alignment inputs found: set phylogenetics.input_subdirs and/or phylogenetics.input_files"
        exit 1
    fi

    # Dynamic threading per software (z_docs/Rules.md):
    #   optimal_threads >= CPU  →  1 job at a time, full CPU per job
    #   optimal_threads  < CPU  →  CPU/optimal_threads parallel jobs
    # Dynamic software-level parallelization:
    #   run up to MAX_PARALLEL phylogeny programs in parallel.
    SOFTWARE_PIDS=()
# Divide total CPU budget across concurrent software programs to avoid oversubscription
num_sw=${#PHYLO_SOFTWARE[@]}
active_sw=$(( num_sw < MAX_PARALLEL ? num_sw : MAX_PARALLEL ))
CPU_PER_SOFTWARE=$(( CPU / active_sw ))
(( CPU_PER_SOFTWARE < 1 )) && CPU_PER_SOFTWARE=1
log_info "CPU budget per software: $CPU_PER_SOFTWARE (total=$CPU, concurrent=$active_sw)"

for software in "${PHYLO_SOFTWARE[@]}"; do
    wait_for_slot "$MAX_PARALLEL"

    (
        case "$software" in
            MEGA_CC)
                optimal=$(get_toml phylogenetics megacc optimal_threads 2>/dev/null || get_toml phylogenetics mega_cc optimal_threads 2>/dev/null || echo "$CPU_PER_SOFTWARE")
                ;;
            IQTREE2)
                optimal=$(get_toml phylogenetics iqtree2 optimal_threads 2>/dev/null || echo "$CPU_PER_SOFTWARE")
                ;;
            RAXML)
                optimal=$(get_toml phylogenetics raxml optimal_threads 2>/dev/null || echo "$CPU_PER_SOFTWARE")
                ;;
            *)
                software_lower=$(echo "$software" | tr 'A-Z' 'a-z')
                optimal=$(get_toml phylogenetics "$software_lower" optimal_threads 2>/dev/null || echo "$CPU_PER_SOFTWARE")
                ;;
        esac

        # Cap optimal_threads to this software's CPU share to prevent oversubscription
        (( optimal > CPU_PER_SOFTWARE )) && optimal=$CPU_PER_SOFTWARE

        if (( optimal >= CPU_PER_SOFTWARE )); then
            max_parallel=1
            threads_per_job=$CPU_PER_SOFTWARE
        else
            max_parallel=$(( CPU_PER_SOFTWARE / optimal ))
            threads_per_job=$optimal
        fi

        log_info "$software: ${threads_per_job} threads/job, max ${max_parallel} parallel"
        log_info "Input scan: ${#ALIGN_INPUT_DIRS[@]} dir(s) + ${#ALIGN_INPUT_FILES[@]} direct file(s) under $ALIGN_DIR ($PHYLO_INPUT_PATTERN)"

        matched_count=0
        submitted_count=0
        skipped_empty_count=0
        inner_pids=()
        while IFS= read -r aligned; do
            matched_count=$(( matched_count + 1 ))
            [[ -f "$aligned" ]] || continue
            if [[ ! -s "$aligned" ]]; then
                log_warn "$software: skipping empty alignment file: $aligned"
                skipped_empty_count=$(( skipped_empty_count + 1 ))
                continue
            fi
            wait_for_slot "$max_parallel"

            # Mirror the full MSA folder layout under PHYLO_DIR. Pass the relative
            # parent of the aligned file (relative to ALIGN_DIR, basename stripped)
            # as --subpath so each tree lands beside its source alignment:
            #   <ALIGN_DIR>/<genome>/<output_subdir>/<set>/<METHOD>_aligned/foo.fas
            #     →  <PHYLO_DIR>/<genome>/<output_subdir>/<set>/<METHOD>_aligned/<software>/foo_<software>.<ext>
            _rel_path="${aligned#$ALIGN_DIR/}"
            _subpath="$(dirname "$_rel_path")"
            [[ "$_subpath" == "." ]] && _subpath=""

            extra_args=()
            config_args=()
            [[ -n "$_subpath" ]] && extra_args+=(--subpath "$_subpath")
            if [[ "$software" == "MEGA_CC" ]]; then
                config_args+=(--config "$MEGACC_CONFIG_MAO")
                [[ -n "$MEGACC_CONFIG_NUC" ]] && extra_args+=(--megacc-config-nuc "$MEGACC_CONFIG_NUC")
                [[ -n "$MEGACC_CONFIG_AA"  ]] && extra_args+=(--megacc-config-aa  "$MEGACC_CONFIG_AA")
            fi

            if [[ "$software" == "IQTREE2" ]]; then
                extra_args+=(--bootstrap "$IQTREE2_BOOTSTRAP")
                extra_args+=(--alrt "$IQTREE2_ALRT")
                extra_args+=(--model "$IQTREE2_MODEL")
                [[ -n "$IQTREE2_MODEL_NUC" ]] && extra_args+=(--iqtree2-model-nuc "$IQTREE2_MODEL_NUC")
                [[ -n "$IQTREE2_MODEL_AA"  ]] && extra_args+=(--iqtree2-model-aa  "$IQTREE2_MODEL_AA")
                extra_args+=(--allnni "$IQTREE2_ALLNNI")
                extra_args+=(--polytomy "$IQTREE2_POLYTOMY")
                extra_args+=(--safe "$IQTREE2_SAFE")
                extra_args+=(--bnni "$IQTREE2_BNNI")
                extra_args+=(--fast "$IQTREE2_FAST")
                extra_args+=(--pers "$IQTREE2_PERS")
                extra_args+=(--redo "$IQTREE2_REDO")
            fi

            if [[ "$software" == "RAXML" ]]; then
                extra_args+=(--raxml-binary "$RAXML_BINARY")
                extra_args+=(--raxml-model  "$RAXML_MODEL_NUC")  # generic fallback = nucleotide model
                extra_args+=(--raxml-model-nuc "$RAXML_MODEL_NUC")
                extra_args+=(--raxml-model-aa  "$RAXML_MODEL_AA")
                extra_args+=(--raxml-seed "$RAXML_SEED")
                extra_args+=(--raxml-mode "$RAXML_MODE")
                extra_args+=(--raxml-bs-trees "$RAXML_BS_TREES")
                extra_args+=(--raxml-search-replicates "$RAXML_SEARCH_REPLICATES")
                extra_args+=(--raxml-redo "$RAXML_REDO")
            fi

            bash "$MODULES/05_phylogenetic_analysis/build_tree.sh" \
                --input "$aligned" \
                --software "$software" \
                --outdir "$PHYLO_DIR" \
                "${config_args[@]}" \
                --threads "$threads_per_job" \
                "${extra_args[@]}" &
            inner_pids+=("$!")
            submitted_count=$(( submitted_count + 1 ))
        done < <({ [[ ${#ALIGN_INPUT_DIRS[@]} -gt 0 ]] && find "${ALIGN_INPUT_DIRS[@]}" -type f -name "$PHYLO_INPUT_PATTERN" 2>/dev/null; printf '%s\n' "${ALIGN_INPUT_FILES[@]:-}"; } | sort -u)

        if (( matched_count == 0 )); then
            log_error "$software: no alignment files found in ${ALIGN_INPUT_DIRS[*]:-$ALIGN_DIR} matching $PHYLO_INPUT_PATTERN"
            exit 1
        fi
        if (( submitted_count == 0 )); then
            log_error "$software: no non-empty alignment files to process (matched=$matched_count, skipped_empty=$skipped_empty_count)"
            exit 1
        fi
        if (( skipped_empty_count > 0 )); then
            log_warn "$software: skipped $skipped_empty_count empty file(s); processing $submitted_count file(s)"
        else
            log_info "$software: processing $submitted_count file(s)"
        fi

        failed=0
        for p in "${inner_pids[@]}"; do
            wait "$p" || failed=$(( failed + 1 ))
        done
        if (( failed > 0 )); then
            log_error "$software: $failed tree build(s) failed"
            exit 1
        fi
    ) &
    SOFTWARE_PIDS+=("$!")
done

    for pid in "${SOFTWARE_PIDS[@]}"; do
        wait "$pid"
    done
else
    log_info "Skipping tree building (not in OPERATIONS)"
fi

# ======================== Tree Visualization ========================
if should_run "visualize_tree"; then
    VIZ_ENABLED=$(get_toml visualization enabled 2>/dev/null || echo "true")
    if [[ "$VIZ_ENABLED" == "true" ]]; then
        log_step "Tree Visualization: $GENE_GROUP"
        bash "$MODULES/05_phylogenetic_analysis/visualize_tree.sh" \
            --treedir "$PHYLO_DIR" \
            --outdir "$PHYLO_DIR" \
            --config "$CONFIG_FILE" \
            --threads "$CPU" \
            --overwrite "$OVERWRITE"
        log_step "Tree visualization complete: $GENE_GROUP"
    else
        log_info "Tree visualization disabled (visualization.enabled = false)"
    fi
else
    log_info "Skipping tree visualization (not in OPERATIONS)"
fi

# ======================== Tree Comparison (Tanglegram) ========================
if should_run "compare_trees"; then
    log_step "Tree Comparison (IQ-TREE2 vs RAxML): $GENE_GROUP"
    bash "$MODULES/05_phylogenetic_analysis/compare_trees.sh" \
        --treedir  "$PHYLO_DIR" \
        --outdir   "$PHYLO_DIR" \
        --config   "$CONFIG_FILE" \
        --threads  "$CPU" \
        --overwrite "$OVERWRITE"
    log_step "Tree comparison complete: $GENE_GROUP"
else
    log_info "Skipping tree comparison (not in OPERATIONS)"
fi

# ======================== Combined Bootstrap Tree ========================
if should_run "combined_bootstrap_tree"; then
    log_step "Combined Bootstrap Tree (IQ-TREE2 UFBoot + RAxML-NG BS): $GENE_GROUP"
    bash "$MODULES/05_phylogenetic_analysis/combined_bootstrap_tree.sh" \
        --treedir  "$PHYLO_DIR" \
        --outdir   "$PHYLO_DIR" \
        --config   "$CONFIG_FILE" \
        --threads  "$CPU" \
        --overwrite "$OVERWRITE"
    log_step "Combined bootstrap tree complete: $GENE_GROUP"
else
    log_info "Skipping combined bootstrap tree (not in OPERATIONS)"
fi

log_step "Phylogenetic analysis complete: $GENE_GROUP"
teardown_logging

done
