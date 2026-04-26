#!/bin/bash
# ==============================================================================
# 03_gea.sh — Gene Expression Analysis orchestrator
#
# Wraps RNA-SeqMAP's three runners (run_alignment_script.sh,
# run_post_processing.sh, run_concordance_analysis.sh) under a single config
# (03_geaCONFIG.toml) and mirrors their figure/CSV/TSV outputs into
# III_RESULT/{GENE_GROUP}/03_Gene_Expression_Analysis/RNA-SeqMAP/.
#
# RNA-SeqMAP is a sibling pipeline (separate workspace, separate `gea` conda
# env). This orchestrator is the only entry point — do not edit RNA-SeqMAP's
# run_*_script.sh files from this codebase. Each runner's hardcoded config
# array (CONFIG_FILES / PIPELINE_CONFIGS / CONCORDANCE_CONFIGS) is patched
# in a temp copy at runtime; the originals stay untouched.
# ==============================================================================
set -euo pipefail

# ---- Paths --------------------------------------------------------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$PIPELINE_DIR/03_geaCONFIG.toml}"
TOML_PARSER="$PIPELINE_DIR/modules/utils/parse_toml.py"
LOGGING_LIB="$PIPELINE_DIR/modules/logging/logging_utils.sh"

[[ -f "$CONFIG_FILE" ]]   || { echo "[ERROR] Config missing: $CONFIG_FILE" >&2; exit 1; }
[[ -f "$TOML_PARSER" ]]   || { echo "[ERROR] parse_toml.py missing: $TOML_PARSER" >&2; exit 1; }
[[ -f "$LOGGING_LIB" ]]   || { echo "[ERROR] logging_utils.sh missing: $LOGGING_LIB" >&2; exit 1; }

# Project-standard 6-component logging: full logs, time/space/combined CSVs,
# error+warning logs, and a software catalog. Set PROJECT_ROOT before sourcing
# so log paths resolve to absolute project-relative paths even if the working
# directory changes (e.g., when child runners cd into RNA-SeqMAP/).
export PROJECT_ROOT="$PIPELINE_DIR"
export STAGE_NAME="03_gea"
# shellcheck disable=SC1090
source "$LOGGING_LIB"

get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# ---- Load configuration -------------------------------------------------------
mapfile -t GENE_GROUPS < <(get_toml pipeline gene_groups)
mapfile -t OPERATIONS  < <(get_toml pipeline operations)
MACHINE="$(get_toml pipeline machine)"
CPU="$(get_toml pipeline compute "$MACHINE" threads)"
MAX_PARALLEL="$(get_toml pipeline compute "$MACHINE" max_parallel)"
CONDA_ENV="$(get_toml pipeline conda_env)"
DRY_RUN="${DRY_RUN:-$(get_toml pipeline dry_run)}"
REGENERATE="${REGENERATE:-$(get_toml pipeline regenerate 2>/dev/null || echo true)}"

RNASEQ_SUBDIR="$(get_toml rnaseq_map rnaseq_dir)"

RESULT_ROOT_REL="$(get_toml output result_root)"
OUTPUT_SUBDIR="$(get_toml output subdir)"
COPY_MODE="$(get_toml output copy_mode)"
MIRROR_SOURCE="$(get_toml output mirror_source)"
mapfile -t MIRROR_EXTENSIONS       < <(get_toml output mirror_extensions)
mapfile -t MIRROR_INCLUDE_PATTERNS < <(get_toml output mirror_include_patterns 2>/dev/null || true)
mapfile -t MIRROR_EXCLUDE_PATTERNS < <(get_toml output mirror_exclude_patterns 2>/dev/null || true)

# Export for child processes (RNA-SeqMAP runner reads THREADS/JOBS via TOML, but
# downstream tools that respect OMP/OpenMP envs benefit from a sane default).
export THREADS="$CPU" JOBS="$MAX_PARALLEL"

RNASEQ_DIR="$PIPELINE_DIR/$RNASEQ_SUBDIR"
RNASEQ_RESULTS_ROOT="$RNASEQ_DIR/$MIRROR_SOURCE"
RESULT_ROOT="$PIPELINE_DIR/$RESULT_ROOT_REL"

[[ -d "$RNASEQ_DIR" ]]    || { echo "[ERROR] RNA-SeqMAP not found: $RNASEQ_DIR" >&2; exit 1; }

# Map TOML operation name -> [section] in 03_geaCONFIG.toml. Hyphens in
# operation names map to underscores in TOML section names (TOML-friendly).
op_to_section() {
    case "$1" in
        Alignment)             echo "alignment" ;;
        Post-Processing)       echo "post_processing" ;;
        Concordance_Analysis)  echo "concordance_analysis" ;;
        *) echo "" ;;
    esac
}

# Track every patched runner across all operations so the EXIT trap can clean
# them all (previously only the last operation's temp file was removed).
_TMP_RUNNERS=()
cleanup() {
    local f
    for f in "${_TMP_RUNNERS[@]}"; do
        [[ -f "$f" ]] && rm -f "$f"
    done
    # Use safe_teardown_logging (bounded, hang-proof) rather than the bare
    # teardown_logging — see logging_utils.sh:safe_teardown_logging() for the
    # rationale. Bare teardown_logging can leave a tee in a process
    # substitution waiting on a pipe that never sees EOF, which on WSL2
    # surfaces as either a hang or a "unexpected EOF while looking for
    # matching `\"'" error attributed to the parent script.
    if type -t safe_teardown_logging &>/dev/null; then
        safe_teardown_logging
    elif type -t teardown_logging &>/dev/null; then
        teardown_logging
    fi
}
trap cleanup EXIT

# ---- Conda activation ---------------------------------------------------------
activate_conda_env() {
    local env_name="$1"
    if command -v conda >/dev/null 2>&1; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        if conda env list | awk '{print $1}' | grep -qx "$env_name"; then
            conda activate "$env_name"
            log_info "Activated conda env: $env_name"
        else
            log_warn "Conda env '$env_name' not found. Run: bash $RNASEQ_DIR/setup_conda_gea.sh"
            log_warn "Continuing with current shell environment."
        fi
    else
        log_warn "conda not on PATH; relying on system tools for STAR/samtools."
    fi
}

# ---- Build a patched RNA-SeqMAP runner ---------------------------------------
# Each RNA-SeqMAP runner hardcodes its config list inside an array assignment
# (CONFIG_FILES, PIPELINE_CONFIGS, or CONCORDANCE_CONFIGS). Rather than edit
# the runners, we materialize a temp copy with that array replaced by the
# operation-specific configs from 03_geaCONFIG.toml, then run the temp copy
# with PROJECT_ROOT pointed at RNA-SeqMAP/.
patch_runner() {
    local src="$1" out="$2" array_name="$3"; shift 3
    local cfgs=("$@")
    awk -v arr="$array_name" -v cfgs="${cfgs[*]}" '
        BEGIN {
            n = split(cfgs, items, " ")
            in_block = 0; replaced = 0
        }
        { sub(/\r$/, "") }                              # normalize CRLF -> LF
        $0 ~ "^"arr"=\\(" {
            print arr "=("
            for (i = 1; i <= n; i++) print "    \"" items[i] "\""
            print ")"
            in_block = 1
            replaced = 1
            next
        }
        in_block {
            if ($0 ~ /^\)/) { in_block = 0 }
            next
        }
        { print }
        END {
            if (!replaced) {
                print "[03_gea] ERROR: " arr "=( marker not found in runner" > "/dev/stderr"
                exit 2
            }
        }
    ' "$src" > "$out"
    chmod +x "$out"
}

# ---- Dispatch a single operation ---------------------------------------------
# Reads runner / config_array / configs from the [<section>] block in
# 03_geaCONFIG.toml, patches the runner, and executes it under
# PROJECT_ROOT=RNA-SeqMAP/.
run_operation() {
    local op="$1"
    local section
    section="$(op_to_section "$op")"
    if [[ -z "$section" ]]; then
        log_warn "Unknown operation '$op' — skipping (expected: Alignment, Post-Processing, Concordance_Analysis)"
        return 0
    fi

    local runner_name array_name
    runner_name="$(get_toml "$section" runner)"
    array_name="$(get_toml "$section" config_array)"
    local -a op_configs=()
    mapfile -t op_configs < <(get_toml "$section" configs)

    local runner_path="$RNASEQ_DIR/$runner_name"
    [[ -x "$runner_path" ]] || { log_error "Runner not executable: $runner_path"; return 1; }
    if [[ ${#op_configs[@]} -eq 0 ]]; then
        log_warn "No configs defined under [$section].configs — skipping $op"
        return 0
    fi

    # Verify every config exists inside RNA-SeqMAP/
    local cfg
    for cfg in "${op_configs[@]}"; do
        [[ -f "$RNASEQ_DIR/$cfg" ]] || { log_error "Config missing: $RNASEQ_DIR/$cfg"; return 1; }
    done

    # Prerequisite checks: some operations consume artefacts produced by an
    # earlier operation. Fail fast with a clear message rather than letting
    # the runner crash deep inside its own logic.
    case "$op" in
        Post-Processing)
            local _align_dir="$RNASEQ_DIR/II_RESULTS/2_ALIGNMENT_RESULTs"
            if [[ ! -d "$_align_dir" ]]; then
                log_error "Post-Processing requires alignment outputs at $_align_dir"
                log_error "  -> Enable \"Alignment\" in [pipeline].operations first, or"
                log_error "     copy 2_ALIGNMENT_RESULTs/ in from a prior run."
                return 1
            fi
            ;;
        Concordance_Analysis)
            local _pp_dir="$RNASEQ_DIR/II_RESULTS/3_POST_PROC"
            if [[ ! -d "$_pp_dir" ]]; then
                log_error "Concordance_Analysis requires post-processing outputs at $_pp_dir"
                log_error "  -> Enable \"Post-Processing\" in [pipeline].operations first."
                return 1
            fi
            ;;
    esac

    log_step "Operation: $op"
    log_info "  Runner:  $runner_name"
    log_info "  Array:   $array_name"
    log_info "  Configs: ${op_configs[*]}"

    local tmp_runner
    tmp_runner="$(mktemp -t "${runner_name%.sh}.XXXXXX.sh")"
    _TMP_RUNNERS+=("$tmp_runner")
    patch_runner "$runner_path" "$tmp_runner" "$array_name" "${op_configs[@]}"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY_RUN=true — patched runner: $tmp_runner"
        log_info "Patched array block:"
        grep -A "$((${#op_configs[@]} + 1))" "^${array_name}=" "$tmp_runner" || true
        return 0
    fi

    log_info "Launching $runner_name under PROJECT_ROOT=$RNASEQ_DIR"
    # Wrap in run_with_space_time_log so each operation gets its own time/space
    # row in logs/space_time_logs/*.csv; --output points at II_RESULTS so the
    # post-run delta size is recorded. `env PROJECT_ROOT=...` only sets the
    # var for the child process — the orchestrator's own PROJECT_ROOT is not
    # mutated, so no save/restore is needed.
    run_with_space_time_log --output "$RNASEQ_DIR/II_RESULTS" \
        env PROJECT_ROOT="$RNASEQ_DIR" bash "$tmp_runner"
}

# ---- Run all enabled operations in order -------------------------------------
run_rnaseqmap() {
    if [[ ${#OPERATIONS[@]} -eq 0 ]]; then
        log_warn "No operations enabled in [pipeline].operations — nothing to run."
        return 0
    fi
    local op
    for op in "${OPERATIONS[@]}"; do
        run_operation "$op"
    done
}

# ---- Mirror only figures + tabular outputs into III_RESULT/{GROUP}/ ----------
# Walks RNA-SeqMAP/<MIRROR_SOURCE> and copies *only* files matching the
# extensions in [output].mirror_extensions (figures, CSVs, TSVs by default).
# BAMs, BAM indices, _STARtmp/, genome indexes, FASTQ, and log files are
# deliberately skipped — they live only in the RNA-SeqMAP workspace.
#
# Relative paths under MIRROR_SOURCE are preserved so files don't collide
# (e.g., 2_ALIGNMENT_RESULTs/<fasta_tag>/star/.../ReadsPerGene.tsv stays
# distinct from 3_POST_PROC/<group>/star/.../ReadsPerGene.tsv).
mirror_outputs() {
    if [[ ! -d "$RNASEQ_RESULTS_ROOT" ]]; then
        log_warn "No RNA-SeqMAP output at $RNASEQ_RESULTS_ROOT — nothing to mirror."
        return 0
    fi
    if [[ ${#MIRROR_EXTENSIONS[@]} -eq 0 ]]; then
        log_warn "No mirror_extensions configured — nothing to mirror."
        return 0
    fi

    # Build a single find expression: ( -iname *.png -o -iname *.svg -o ... )
    local -a find_expr=( -type f \( )
    local first=1 ext
    for ext in "${MIRROR_EXTENSIONS[@]}"; do
        # Normalize: accept ".png" or "png"; emit "*.png"
        local pat="*.${ext#.}"
        if (( first )); then
            find_expr+=( -iname "$pat" ); first=0
        else
            find_expr+=( -o -iname "$pat" )
        fi
    done
    find_expr+=( \) )

    # Returns 0 (mirror) / 1 (skip) for the relative path in $1.
    # Implements the two-stage filter from [output].mirror_{include,exclude}_patterns:
    #   - Skip if any exclude pattern matches.
    #   - Skip if include list is non-empty AND no include pattern matches.
    # Empty include list = match everything (then excludes still apply).
    _path_passes_filters() {
        local rel="$1" pat
        for pat in "${MIRROR_EXCLUDE_PATTERNS[@]}"; do
            [[ -n "$pat" && "$rel" =~ $pat ]] && return 1
        done
        if [[ ${#MIRROR_INCLUDE_PATTERNS[@]} -eq 0 ]]; then
            return 0
        fi
        for pat in "${MIRROR_INCLUDE_PATTERNS[@]}"; do
            [[ -n "$pat" && "$rel" =~ $pat ]] && return 0
        done
        return 1
    }

    local group dest src rel target count skipped
    for group in "${GENE_GROUPS[@]}"; do
        dest="$RESULT_ROOT/$group/$OUTPUT_SUBDIR"
        mkdir -p "$dest"
        log_info "Mirroring outputs (${MIRROR_EXTENSIONS[*]}) -> $dest  (mode=$COPY_MODE)"
        if [[ ${#MIRROR_INCLUDE_PATTERNS[@]} -gt 0 ]]; then
            log_info "  Include filters: ${MIRROR_INCLUDE_PATTERNS[*]}"
        fi
        if [[ ${#MIRROR_EXCLUDE_PATTERNS[@]} -gt 0 ]]; then
            log_info "  Exclude filters: ${MIRROR_EXCLUDE_PATTERNS[*]}"
        fi

        count=0; skipped=0
        while IFS= read -r -d '' src; do
            rel="${src#$RNASEQ_RESULTS_ROOT/}"
            if ! _path_passes_filters "$rel"; then
                ((skipped++)) || true
                continue
            fi
            target="$dest/$rel"
            mkdir -p "$(dirname "$target")"
            # Replace any existing target so reruns don't accumulate stale files
            # (and don't accidentally write through an old hardlink).
            rm -f "$target"
            # -p preserves mode + mtime so the mirrored file timestamps reflect
            # when RNA-SeqMAP produced each artefact, not when the mirror ran.
            if [[ "$COPY_MODE" == "link" ]]; then
                cp -l "$src" "$target" 2>/dev/null || cp -p "$src" "$target"
            else
                cp -p "$src" "$target"
            fi
            ((count++)) || true
        done < <(find "$RNASEQ_RESULTS_ROOT" "${find_expr[@]}" -print0)

        log_info "Mirrored $count file(s) under $dest  (skipped $skipped by filter)"
    done
}

# ---- Main --------------------------------------------------------------------
main() {
    # Initialize 6-component logging. Pass "true" to clear previous logs;
    # default keeps prior runs for comparison. Sets up dual-output (terminal
    # + log files) via process substitution.
    setup_logging "${CLEAR_LOGS:-false}"

    log_step "Eggplant_Pipeline GEA orchestrator -> RNA-SeqMAP (STAR / SmelDMP)"
    log_info "Config:     $CONFIG_FILE"
    log_info "GENE_GROUPS=${GENE_GROUPS[*]}  MACHINE=$MACHINE  CPU=$CPU  MAX_PARALLEL=$MAX_PARALLEL"
    log_info "OPERATIONS: ${OPERATIONS[*]:-<none>}"
    log_info "REGENERATE=$REGENERATE  DRY_RUN=$DRY_RUN  COPY_MODE=$COPY_MODE  MIRROR_EXTENSIONS=${MIRROR_EXTENSIONS[*]}"

    # Capture tool versions once per run for reproducibility audit
    # (logs/software_catalogs/*.csv). Skipped under DRY_RUN to keep dry runs cheap.
    if [[ "$DRY_RUN" != "true" ]] && type -t catalog_all_software &>/dev/null; then
        catalog_all_software 2>/dev/null || log_warn "catalog_all_software failed (non-fatal)"
        # catalog_all_software installs a RETURN trap that references its local
        # _ver_tmpdir; under `set -u` that trap re-fires with the var out of
        # scope when main() later returns, crashing with "unbound variable".
        # Clear the leaked trap to break that chain.
        trap - RETURN
    fi

    if [[ "$REGENERATE" == "true" ]]; then
        activate_conda_env "$CONDA_ENV"
        run_rnaseqmap
    else
        log_info "REGENERATE=false -> skipping RNA-SeqMAP runners; mirroring existing $RNASEQ_RESULTS_ROOT only."
    fi

    log_step "Mirroring outputs"
    mirror_outputs

    log_step "GEA pipeline complete"
    log_info "Logs:  $LOG_DIR/"
    log_info "       $TIME_DIR/  $SPACE_DIR/  $SPACE_TIME_DIR/"
    log_info "       $ERROR_WARN_DIR/  $SOFTWARE_CATALOG_DIR/"
}

main "$@"
