#!/bin/bash
# ============================================================================
# Program 3: Multiple Sequence Alignment
# ============================================================================
# Reads [alignment] config to determine methods and input sets, then runs
# align.sh for every (input_set × method) combination.
#
# Input sources (referenced in-place from earlier pipeline steps):
#   - HMMER-identified sequences : 01_Identification → 02_BLAST_Alignment
#   - Cross-species query FASTAs: ortholog_blast.query_fastas
#
# Each input_set is a labeled combination of (scope × genome × sequence_type).
# Outputs go to 04_MSA/<top_folder>/All_Result/<set_name>/<METHOD>_aligned/ so
# basenames never collide and genome-specific sets are separated at the first folder.
# Selected/curated results live in 04_MSA/<top_folder>/Selected_Result/v1_Full/.
#
# Threading strategy — goal: 100% CPU saturation across all alignment jobs.
#   Rule: parallel_jobs = CPU / optimal_threads  (floored, min 1)
#         threads_per_job = optimal_threads (or CPU when optimal >= CPU)
#
#   Multi-threaded tools (CLUSTALO, MAFFT, MUSCLE v5):
#     optimal_threads = CPU → 1 job at full CPU → CPU threads/job, 1 parallel
#   Single-threaded tools (MUSCLE v3, PROBCONS, CLUSTALW):
#     optimal_threads = 1  → CPU parallel jobs  → 1 thread/job, CPU parallel
#
#   Both cases give CPU × 1 = $(nproc) total threads in use.
#   Controlled by optimal_threads in each [alignment.<method>] config block.
#
# Usage:
#   bash d_msa_alignment.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# Comment in/out the gene groups to process:
GENE_GROUPS=(
    "DMP"
    #"GRF_GIF"
    #"PLA"
)

# All other settings (overwrite, threads, parallelism) are loaded from
# 04_msa_alignmentCONFIG.toml [pipeline.compute.$MACHINE] section.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

# ---------------------------------------------------------------------------
# wait_for_slot <max_parallel>
#   Block until fewer than max_parallel background jobs are running.
# ---------------------------------------------------------------------------
wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 0.5
    done
}

# ---------------------------------------------------------------------------
# resolve_top_folder <set_name>
#   Route set outputs to a genome-specific first-level folder when possible.
# ---------------------------------------------------------------------------
resolve_top_folder() {
    local set_name="$1"
    case "$set_name" in
        *smel_v4_1*) echo "Solanum_melongena_v4.1" ;;
        *gpe001970*|*selected_v1*|*selected_v2*) echo "GPE001970_SMEL5" ;;
        *) echo "shared" ;;
    esac
}

# ---------------------------------------------------------------------------
# resolve_output_subdir <set_name> <config_file>
#   Determine which sub-folder under <top_folder>/ to write results into.
#   Priority: explicit output_subdir in config → name-based pattern → All_Result
#   Naming rules (applied to set_name):
#     *selected*v2* | *v2_reduced*          → Selected_Result/v2_Reduced
#     *selected*                             → Selected_Result/v1_Full
#     (default)                              → All_Result
# ---------------------------------------------------------------------------
resolve_output_subdir() {
    local set_name="$1"
    # 1. Honour explicit config override
    local explicit
    explicit=$(get_toml alignment input_sets "$set_name" output_subdir 2>/dev/null) || explicit=""
    if [[ -n "$explicit" ]]; then
        echo "$explicit"
        return
    fi
    # 2. Infer from set name
    case "$set_name" in
        *selected*v2*|*v2_reduced*) echo "Selected_Result/v2_Reduced" ;;
        *selected*)                  echo "Selected_Result/v1_Full"    ;;
        *)                           echo "All_Result"                 ;;
    esac
}

# ---------------------------------------------------------------------------
# Accumulate temp files for cleanup (handles multiple gene groups correctly)
# ---------------------------------------------------------------------------
TEMP_FILES=()
cleanup_all() {
    # Note: do NOT pre-kill jobs -rp here -- that list includes the logging tee
    # process substitutions, and killing them before teardown_logging restores
    # FDs would lose buffered log writes. safe_teardown_logging closes pipe
    # write-ends first (giving tees clean EOF), then SIGTERMs leftover children
    # (orphaned aligner subprocesses) at the end.
    rm -f "${TEMP_FILES[@]+"${TEMP_FILES[@]}"}" 2>/dev/null || true
    safe_teardown_logging
}
trap cleanup_all EXIT

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# --- Resolve config: split directory or monolithic file --------------------
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_msa_cfg_XXXXXX.toml")
    TEMP_FILES+=("$CONFIG_FILE")
    cat "$PIPELINE_DIR/04_msa_alignmentCONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/04_multiple_sequence_alignment.toml" \
        "$CONFIG_DIR/02_blast_ortholog_alignment.toml" > "$CONFIG_FILE"
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "$CPU")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")
IFS=' ' read -ra ALIGNMENT_METHODS <<< "$(get_toml alignment methods)"

MSA_DIR_NAME=$(get_toml output_dirs msa 2>/dev/null || echo "04_MSA")
ALIGN_DIR="$BASE_DIR/$MSA_DIR_NAME"
MERGE_DIR="$ALIGN_DIR/merged_input/All"
mkdir -p "$ALIGN_DIR" "$MERGE_DIR"

setup_logging
log_step "Multiple Sequence Alignment: $GENE_GROUP (OVERWRITE=$OVERWRITE)"

# =========================================================================
# Discover input sets and build (set_name, fasta_path) pairs
# =========================================================================
# The active_input_sets list in [alignment] controls which sets run.
# Comment in/out entries there — no per-set "enabled" flag needed.
# Arrays: SET_NAMES[i] and SET_FASTAS[i] are parallel arrays
SET_NAMES=()
SET_FASTAS=()

IFS=' ' read -ra ACTIVE_SETS <<< "$(get_toml alignment active_input_sets 2>/dev/null || true)"
if [[ ${#ACTIVE_SETS[@]} -eq 0 ]]; then
    log_error "alignment.active_input_sets is empty or missing — nothing to align"
    teardown_logging
    continue
fi
log_info "Active input sets: ${#ACTIVE_SETS[@]}"

# --- Handle OVERWRITE: clean only the active sets' previous outputs --------
if [[ "$OVERWRITE" == "true" ]]; then
    log_info "OVERWRITE=true — removing outputs for active input sets only"
    for _set in "${ACTIVE_SETS[@]}"; do
        _set=$(echo "$_set" | tr 'A-Z' 'a-z')
        _top_folder="$(resolve_top_folder "$_set")"
        _output_subdir="$(resolve_output_subdir "$_set")"
        _outdir="$ALIGN_DIR/$_top_folder/$_output_subdir/$_set"
        if [[ -d "$_outdir" ]]; then
            log_info "  Removing: $_top_folder/$_output_subdir/$_set"
            rm -rf "$_outdir"
        fi
        _merged="$MERGE_DIR/${_set}.fa"
        [[ -f "$_merged" ]] && { log_info "  Removing merged: ${_set}.fa"; rm -f "$_merged"; }
    done
fi

for set_name in "${ACTIVE_SETS[@]}"; do
    set_name=$(echo "$set_name" | tr 'A-Z' 'a-z')

    is_merge=$(get_toml alignment input_sets "$set_name" merge 2>/dev/null) || is_merge="false"
    mapfile -t fastas_raw < <(get_toml alignment input_sets "$set_name" fastas 2>/dev/null) || fastas_raw=()

    if [[ "$is_merge" == "true" ]]; then
        # --- Merge mode: concatenate listed FASTAs + query_fastas into one file
        merged="$MERGE_DIR/${set_name}.fa"
        if [[ ! -s "$merged" ]]; then
            log_info "Merging: $set_name"
            > "$merged"
            # Safely iterate over fastas_raw (may be empty with set -u)
            for rel in "${fastas_raw[@]+"${fastas_raw[@]}"}"; do
                abs="$BASE_DIR/$rel"
                if [[ -s "$abs" ]]; then
                    cat "$abs" >> "$merged"
                    # Ensure trailing newline so next file's > header starts on a new line
                    [[ -s "$merged" ]] && sed -i -e '$a\' "$merged"
                else
                    log_warn "  Missing: $rel"
                fi
            done
            query_fastas_key=$(get_toml alignment input_sets "$set_name" query_fastas_key 2>/dev/null) || query_fastas_key=""
            if [[ -z "$query_fastas_key" ]]; then
                if [[ "$set_name" == *"amino_acid"* ]]; then
                    query_fastas_key=$(get_toml msa_query_selection amino_acid_query_fastas_key 2>/dev/null) || query_fastas_key="query_protein_fastas"
                else
                    query_fastas_key=$(get_toml msa_query_selection nucleotide_query_fastas_key 2>/dev/null) || query_fastas_key="query_fastas"
                fi
            fi

            mapfile -t QUERY_FASTAS < <(get_toml ortholog_blast "$query_fastas_key" 2>/dev/null) || QUERY_FASTAS=()
            if [[ ${#QUERY_FASTAS[@]} -eq 0 && "$query_fastas_key" == "query_protein_fastas" ]]; then
                log_warn "  ortholog_blast.query_protein_fastas is empty — falling back to query_fastas"
                mapfile -t QUERY_FASTAS < <(get_toml ortholog_blast query_fastas 2>/dev/null) || QUERY_FASTAS=()
            fi
            # Safely iterate over QUERY_FASTAS (may be empty with set -u)
            for qf in "${QUERY_FASTAS[@]+"${QUERY_FASTAS[@]}"}"; do
                # Note: query FASTAs are relative to PIPELINE_DIR (II_INPUTS/), not BASE_DIR
                qf_abs="$PIPELINE_DIR/$qf"
                if [[ -s "$qf_abs" ]]; then
                    cat "$qf_abs" >> "$merged"
                    # Ensure trailing newline so next file's > header starts on a new line
                    [[ -s "$merged" ]] && sed -i -e '$a\' "$merged"
                else
                    log_warn "  Query FASTA missing: $qf"
                fi
            done
            before_dedup=$(grep -c '^>' "$merged" 2>/dev/null || echo 0)

            # --- Deduplicate: remove sequences whose header ID (first word) repeats
            python3 -c "
import sys
fpath = sys.argv[1]
seqs, seen, dupes = [], set(), 0
with open(fpath) as fh:
    hdr, lines = None, []
    for line in fh:
        if line.startswith('>'):
            if hdr: seqs.append((hdr, lines))
            hdr, lines = line, []
        else:
            lines.append(line)
    if hdr: seqs.append((hdr, lines))
with open(fpath, 'w') as fh:
    for h, ls in seqs:
        sid = h.strip().split()[0]
        if sid in seen:
            dupes += 1; continue
        seen.add(sid)
        fh.write(h)
        for l in ls: fh.write(l)
if dupes: print(f'  Removed {dupes} duplicate sequence(s)')
" "$merged"
            after_dedup=$(grep -c '^>' "$merged" 2>/dev/null || echo 0)
            log_info "  -> $after_dedup sequences (dedup removed $((before_dedup - after_dedup)))"
        fi
        [[ -s "$merged" ]] && { SET_NAMES+=("$set_name"); SET_FASTAS+=("$merged"); }
    else
        # --- Reference mode: point to each file in-place
        # Safely iterate over fastas_raw (may be empty with set -u)
        for rel in "${fastas_raw[@]+"${fastas_raw[@]}"}"; do
            abs="$BASE_DIR/$rel"
            if [[ -s "$abs" ]]; then
                SET_NAMES+=("$set_name"); SET_FASTAS+=("$abs")
            else
                log_warn "Missing: $rel"
            fi
        done
    fi
done

# Count unique set names safely (handle empty array with set -u)
if [[ ${#SET_NAMES[@]} -gt 0 ]]; then
    unique_sets=$(echo "${SET_NAMES[@]}" | tr ' ' '\n' | sort -u | wc -l)
else
    unique_sets=0
fi
log_info "Input sets resolved: ${#SET_FASTAS[@]} FASTA files across ${unique_sets} sets"

# --- Early exit if no valid input sets found ------------------------------
if [[ ${#SET_FASTAS[@]} -eq 0 ]]; then
    log_error "No valid input sets found — skipping alignment for $GENE_GROUP"
    teardown_logging
    continue
fi

# =========================================================================
# Run alignments — dynamic threading per method
# =========================================================================
# Threading strategy (z_docs/Rules.md):
#   optimal_threads >= CPU  →  1 job at a time, full CPU per job
#   optimal_threads  < CPU  →  CPU/optimal_threads parallel jobs

for method in "${ALIGNMENT_METHODS[@]}"; do
    method_lower=$(echo "$method" | tr 'A-Z' 'a-z')
    optimal=$(get_toml alignment "$method_lower" optimal_threads 2>/dev/null) || optimal=1

    if (( optimal >= CPU )); then
        max_parallel=1
        threads_per_job=$CPU
    else
        max_parallel=$(( CPU / optimal ))
        threads_per_job=$optimal
    fi
    # Cap at MAX_PARALLEL if set
    if (( max_parallel > MAX_PARALLEL )); then
        max_parallel=$MAX_PARALLEL
    fi
    total_cpu=$(( threads_per_job * max_parallel ))

    log_step "Running $method"
    log_info "  Threading: ${threads_per_job} threads/job × ${max_parallel} parallel jobs = ${total_cpu}/${CPU} CPUs used"

    ALIGN_PIDS=()
    for i in "${!SET_FASTAS[@]}"; do
        set_name="${SET_NAMES[$i]}"
        fasta="${SET_FASTAS[$i]}"
        top_folder="$(resolve_top_folder "$set_name")"
        output_subdir="$(resolve_output_subdir "$set_name")"
        outdir="$ALIGN_DIR/$top_folder/$output_subdir/$set_name"
        mkdir -p "$outdir"

        wait_for_slot "$max_parallel"
        log_info "  Launching: $set_name ($(basename "$fasta"))"
        bash "$MODULES/04_multiple_sequence_alignment/align.sh" \
            --input  "$fasta" \
            --method "$method" \
            --outdir "$outdir" \
            --threads "$threads_per_job" \
            --config "$CONFIG_FILE" &
        ALIGN_PIDS+=($!)
    done

    log_info "Waiting for ${#ALIGN_PIDS[@]} $method job(s) to complete..."
    FAILED=0
    for pid in "${ALIGN_PIDS[@]+"${ALIGN_PIDS[@]}"}"; do
        wait "$pid" || (( ++FAILED ))
    done
    if (( FAILED > 0 )); then
        log_error "$FAILED $method alignment job(s) failed"
        exit 1
    fi
    log_info "All $method jobs completed successfully"
done

log_step "Alignment complete: $GENE_GROUP"
teardown_logging

done
