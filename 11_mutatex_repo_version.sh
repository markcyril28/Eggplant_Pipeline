#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# CONFIGURATION — edit variables below for easy per-run overrides.
# All other settings (CPU/threads, parallelism, force-rerun) come from TOML.
# Edit 11_mutatex_repo_versionCONFIG.toml to change those settings.
# ==============================================================================

GENE_GROUP="DMP-HAP2"   # Gene group to process (matches config/PPI/ and III_RESULT/ paths)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/modules/PPI/config_parser.sh"

CONFIG_DIR="${SCRIPT_DIR}/config/PPI"
load_config "${SCRIPT_DIR}/11_mutatex_repo_versionCONFIG.toml"

# ------------------------------------------------------------------------------
# Select compute profile (Local vs HPC) - same pattern as 05_phyloCONFIG.toml.
# Edit `machine` in the TOML to switch profiles; do not edit this script.
# ------------------------------------------------------------------------------
MACHINE="$(toml_get "pipeline.machine" "Local")"

# ------------------------------------------------------------------------------
# Bulk-load scalar config (single Python3 invocation instead of N forks).
# Compute keys (threads, max_parallel_pdbs) come from the selected machine profile.
# ------------------------------------------------------------------------------
eval "$(toml_get_bulk \
    _cfg_env=environment.conda_env:PPI \
    "_cfg_results=paths.results_dir:III_RESULT/${GENE_GROUP}/11_PPI_MutateX" \
    "_cfg_threads=pipeline.compute.${MACHINE}.threads:0" \
    _cfg_nruns=settings.nruns:3 \
    _cfg_force=settings.force_rerun:true \
    _cfg_cleanup=cleanup.enabled:true \
    _cfg_interval=cleanup.worker_interval_seconds:30 \
    _cfg_debug=settings.debug.debug_mode:true \
    _cfg_verbose=settings.debug.verbose_logging:true \
    "_cfg_max_parallel_pdbs=pipeline.compute.${MACHINE}.max_parallel_pdbs:0" \
    _cfg_timeout_hours=settings.per_pdb_timeout_hours:48 \
    _cfg_pp_workers=settings.post_processing_workers:0 \
)"

ENV_NAME="${ENV_NAME:-$_cfg_env}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/$_cfg_results}"
THREADS_CFG="$_cfg_threads"
NRUNS="$_cfg_nruns"
FORCE_RERUN="${FORCE_RERUN:-$_cfg_force}"
CLEANUP_ENABLED="${CLEANUP_ENABLED:-$_cfg_cleanup}"
CLEANUP_WORKER_INTERVAL="${CLEANUP_WORKER_INTERVAL:-$_cfg_interval}"
DEBUG_MODE="${DEBUG_MODE:-$_cfg_debug}"
VERBOSE_LOGGING="${VERBOSE_LOGGING:-$_cfg_verbose}"
MAX_PARALLEL_PDBS_CFG="$_cfg_max_parallel_pdbs"
PER_PDB_TIMEOUT_HOURS="$_cfg_timeout_hours"
PP_WORKERS_CFG="$_cfg_pp_workers"

# ------------------------------------------------------------------------------
# 3. Input PDB Files
#    Listed under [structures] active = [...] in a_mutatex.toml
# ------------------------------------------------------------------------------
PDB_FILES_CONFIG=()
while IFS= read -r rel_path; do
    [[ -n "$rel_path" ]] && PDB_FILES_CONFIG+=("$rel_path")
done < <(toml_get_array "structures.active")

# ==============================================================================
# SCRIPT SETUP - DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU'RE DOING
# ==============================================================================

# Source logging utilities
source "$SCRIPT_DIR/modules/logging/logging_utils.sh"
source "$SCRIPT_DIR/modules/logging/gpu_utils.sh"
source "$SCRIPT_DIR/modules/PPI/mutatex_utils.sh"

# Set up logging with absolute paths (needed because script changes directories)
printf -v RUN_ID 'mutatex_%(%Y%m%d_%H%M%S)T' -1
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
ERROR_WARN_DIR="$SCRIPT_DIR/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$SCRIPT_DIR/logs/software_catalogs"
GPU_LOG_DIR="$SCRIPT_DIR/logs/gpu_log"
CONSOLE_LOG_DIR="$SCRIPT_DIR/logs/console_logs"
setup_logging

# Robust cleanup: kill background jobs, teardown logging, remove temp files
cleanup_exit() {
    local bg_pids; bg_pids=$(jobs -rp 2>/dev/null)
    [[ -n "$bg_pids" ]] && kill $bg_pids 2>/dev/null
    wait 2>/dev/null
    teardown_logging
    [[ -d "${STATUS_DIR:-}" ]] && rm -rf "$STATUS_DIR"
}
trap 'cleanup_exit' EXIT

# Create console log directory and file for full output capture
mkdir -p "$CONSOLE_LOG_DIR"
CONSOLE_LOG="$CONSOLE_LOG_DIR/console_${RUN_ID}.log"
TROUBLESHOOT_LOG="$CONSOLE_LOG_DIR/troubleshoot_${RUN_ID}.log"
touch "$CONSOLE_LOG" "$TROUBLESHOOT_LOG"

# ==============================================================================
# TROUBLESHOOTING HELPER FUNCTIONS
# ==============================================================================

log_debug() {
  [[ "$DEBUG_MODE" != "true" ]] && return 0
  local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
  local msg="[$ts] [DEBUG] $*"
  echo "$msg"
  echo "$msg" >> "$TROUBLESHOOT_LOG"
}

log_troubleshoot() {
  local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
  echo "[$ts] [TROUBLESHOOT] $*" >> "$TROUBLESHOOT_LOG"
}

# Generic section dumper — writes a titled block to the troubleshoot log
# Usage: _dump_section "TITLE" "key1: val1" "key2: val2" ...
_dump_section() {
  local title="$1"; shift
  local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
  {
    echo "[$ts] [TROUBLESHOOT] === $title ==="
    for line in "$@"; do echo "[$ts] [TROUBLESHOOT] $line"; done
    echo "[$ts] [TROUBLESHOOT] === END $title ==="
  } >> "$TROUBLESHOOT_LOG"
}

dump_environment() {
  _dump_section "ENVIRONMENT DUMP" \
    "Date: $(printf '%(%Y-%m-%d %H:%M:%S %Z)T' -1)" \
    "Hostname: $(hostname)" \
    "User: $(whoami)" \
    "Working Directory: $(pwd)" \
    "Script Directory: $SCRIPT_DIR" \
    "Shell: $SHELL" \
    "Bash Version: $BASH_VERSION" \
    "PATH: $PATH" \
    "PYTHONPATH: ${PYTHONPATH:-not set}" \
    "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:-not set}"
}

dump_system_info() {
  _dump_section "SYSTEM INFO" \
    "OS: $(uname -a)" \
    "CPU cores: $(nproc 2>/dev/null || echo 'unknown')" \
    "Memory: $(free -h 2>/dev/null | grep Mem || echo 'unknown')" \
    "Disk space (script dir): $(df -h "$SCRIPT_DIR" 2>/dev/null | tail -1 || echo 'unknown')"
}

dump_conda_info() {
  local lines=(
    "Conda env: $ENV_NAME"
    "Conda base: ${CONDA_BASE:-not set}"
    "Active env: ${CONDA_DEFAULT_ENV:-none}"
  )
  if command -v conda &>/dev/null; then
    lines+=(
      "Conda version: $(conda --version 2>&1)"
      "Python: $(python --version 2>&1)"
      "Python path: $(which python 2>&1)"
    )
  fi
  _dump_section "CONDA INFO" "${lines[@]}"
}

dump_foldx_info() {
  _dump_section "FOLDX CONFIGURATION" \
    "FoldX binary: $FOLDX_BINARY" \
    "FoldX exists: $(test -f "$FOLDX_BINARY" && echo 'yes' || echo 'NO - MISSING!')" \
    "FoldX executable: $(test -x "$FOLDX_BINARY" && echo 'yes' || echo 'NO - NOT EXECUTABLE!')" \
    "Rotabase: $ROTABASE" \
    "Rotabase exists: $(test -f "$ROTABASE" && echo 'yes' || echo 'NO - MISSING!')"
}

dump_pdb_info() {
  local pdb_file="$1"
  local lines=("Exists: $(test -f "$pdb_file" && echo 'yes' || echo 'NO - MISSING!')")
  if [[ -f "$pdb_file" ]]; then
    lines+=(
      "Size: $(du -h "$pdb_file" 2>/dev/null | cut -f1)"
      "Lines: $(wc -l < "$pdb_file" 2>/dev/null || echo 'unknown')"
      "Chains: $(grep ^ATOM "$pdb_file" 2>/dev/null | cut -c22 | sort -u | tr '\n' ' ' || echo 'unknown')"
    )
  fi
  _dump_section "PDB FILE INFO: $pdb_file" "${lines[@]}"
}

log_mutatex_invocation() {
  local pdb="$1" outdir="$2" threads="${3:-$THREADS_PER_PDB}"
  local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
  _dump_section "MUTATEX INVOCATION" \
    "Timestamp: $ts" \
    "Input PDB: $pdb" \
    "Output dir: $outdir" \
    "Parallel FoldX jobs: $threads" \
    "N-runs: $NRUNS" \
    "FoldX: $FOLDX_BINARY" \
    "Command: mutatex $pdb -p $threads -n $NRUNS -x $FOLDX_BINARY -b $ROTABASE -f suite5 -R $REPAIR_TEMPLATE -M $MUTATE_TEMPLATE -I $INTERFACE_TEMPLATE -B -v -l"
}

log_run_result() {
  local pdb="$1" exit_status="$2" outdir="$3"
  local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
  local lines=(
    "PDB: $pdb"
    "Exit status: $exit_status"
    "Timestamp: $ts"
  )
  if [[ -d "$outdir" ]]; then
    lines+=(
      "Output dir size: $(du -sh "$outdir" 2>/dev/null | cut -f1 || echo 'unknown')"
      "Output files count: $(find "$outdir" -type f 2>/dev/null | wc -l || echo 'unknown')"
    )
  fi
  _dump_section "RUN RESULT" "${lines[@]}"
}

# Initialize troubleshooting session
_dump_section "MUTATEX TROUBLESHOOTING SESSION START" "RUN_ID: $RUN_ID"

log_info "Console log: $CONSOLE_LOG"
log_info "Troubleshoot log: $TROUBLESHOOT_LOG"

log_step "Starting MutateX batch processing"

# Activate the conda environment
CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Dump initial troubleshooting info
if [[ "$DEBUG_MODE" == "true" ]]; then
  dump_environment
  dump_system_info
  dump_conda_info
  log_debug "Conda environment activated: $ENV_NAME"
fi

# Source mutatex environment hints if available
ENV_HINT_FILE="${ENV_HINT_FILE:-$SCRIPT_DIR/mutatex.env}"
if [[ -f "$ENV_HINT_FILE" ]]; then
  log_info "Loading environment hints from $ENV_HINT_FILE"
  source "$ENV_HINT_FILE"
fi

# Build full paths for PDB files from config
# Structure paths are relative to input_base (same convention as GROMACS pipeline)
INPUT_BASE="$SCRIPT_DIR/$(toml_get "paths.input_base" "III_RESULT/${GENE_GROUP}/08_Protein_Structure/GPE001970_SMEL5/AlphaFold3_Results")"
PDB_FILES=()
for pdb in "${PDB_FILES_CONFIG[@]}"; do
  PDB_FILES+=("$INPUT_BASE/$pdb")
done

# FoldX configuration (from TOML, paths relative to SCRIPT_DIR)
FOLDX_BINARY="$SCRIPT_DIR/$(toml_get "paths.foldx_binary" "modules/PPI/foldx_20270131_5.1")"
ROTABASE="$SCRIPT_DIR/$(toml_get "paths.rotabase" "modules/PPI/rotabase.txt")"

# Templates (from TOML, paths relative to SCRIPT_DIR)
REPAIR_TEMPLATE="$SCRIPT_DIR/$(toml_get "paths.templates.repair" "modules/PPI/repair_runfile_template.txt")"
MUTATE_TEMPLATE="$SCRIPT_DIR/$(toml_get "paths.templates.mutate" "modules/PPI/mutate_runfile_template.txt")"
INTERFACE_TEMPLATE="$SCRIPT_DIR/$(toml_get "paths.templates.interface" "modules/PPI/interface_runfile_template.txt")"

log_info "Using FoldX: $FOLDX_BINARY"
log_info "Using templates from: $SCRIPT_DIR/modules/"

# ==============================================================================
# AUTO-DETECT CPU AND COMPUTE PARALLELISM
#
# FoldX is inherently single-threaded (hardcoded). mutatex -p controls how many
# parallel single-threaded FoldX jobs run simultaneously. To saturate the machine,
# we auto-detect available CPUs, then determine:
#   1) How many PDBs to process in parallel (PARALLEL_PDBS)
#   2) How many FoldX jobs each mutatex instance runs (THREADS_PER_PDB)
#
# Formula: TOTAL_CPUS = PARALLEL_PDBS × THREADS_PER_PDB
# This adapts automatically: 12-core laptop, 64-core server, etc.
# ==============================================================================
TOTAL_CPUS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
(( TOTAL_CPUS < 1 )) && TOTAL_CPUS=1
NUM_PDBS=${#PDB_FILES[@]}

log_info "Compute profile: ${MACHINE} (threads=${THREADS_CFG}, max_parallel_pdbs=${MAX_PARALLEL_PDBS_CFG})"

if (( NUM_PDBS == 0 )); then
  log_error "No active PDB files found. Check [structures].active in $CONFIG_DIR/mutatex/a_mutatex.toml"
  _dump_section "CONFIGURATION ERROR" \
    "Reason: [structures].active resolved to an empty list" \
    "Input base: $INPUT_BASE" \
    "Config file: $CONFIG_DIR/mutatex/a_mutatex.toml"
  exit 1
fi

# FoldX optimal threads = 1 (single-threaded by design — do not change)
FOLDX_OPTIMAL_THREADS=1

if (( THREADS_CFG > 0 )); then
  # User explicitly set threads in TOML → use as parallel FoldX jobs, 1 PDB at a time
  THREADS_PER_PDB=$THREADS_CFG
  PARALLEL_PDBS=1
  log_info "Config override: ${THREADS_PER_PDB} FoldX jobs, 1 PDB at a time"
elif (( MAX_PARALLEL_PDBS_CFG > 0 )); then
  # User explicitly set max parallel PDBs
  PARALLEL_PDBS=$MAX_PARALLEL_PDBS_CFG
  (( PARALLEL_PDBS > NUM_PDBS )) && PARALLEL_PDBS=$NUM_PDBS
  (( PARALLEL_PDBS < 1 )) && PARALLEL_PDBS=1
  THREADS_PER_PDB=$(( TOTAL_CPUS / PARALLEL_PDBS ))
  (( THREADS_PER_PDB < 1 )) && THREADS_PER_PDB=1
  log_info "Config: ${PARALLEL_PDBS} parallel PDBs × ${THREADS_PER_PDB} FoldX jobs (${TOTAL_CPUS} CPUs)"
else
  # Full auto-detect: saturate all CPUs aggressively.
  # FoldX is single-threaded → maximize parallelization at every level.
  #
  # Wall time of a (P, T) layout ≈ ceil(N/P) / T (one wave-time per wave).
  # Subject to P*T ≤ TOTAL_CPUS, T ≥ MIN_FOLDX_JOBS_PER_PDB, 1 ≤ P ≤ N.
  # Equivalently, minimize work = ceil(N/P) * P (waves * concurrency);
  # smaller work means fewer "wasted" CPU-slots in the tail wave.
  #
  # Greedy max-P (the previous heuristic) loses badly when N does not
  # divide TOTAL_CPUS/T_min cleanly: e.g. N=7 PDBs on 12 CPUs gave
  # P=6×T=2 (wall=2 unit-times) instead of optimal P=1×T=12 (wall=7/12).
  # Scaling from laptops (12-core) to servers (64+ core) with no hardcoded cap.
  MIN_FOLDX_JOBS_PER_PDB=2
  P_MAX_BY_T=$(( TOTAL_CPUS / MIN_FOLDX_JOBS_PER_PDB ))
  (( P_MAX_BY_T < 1 )) && P_MAX_BY_T=1
  P_MAX=$NUM_PDBS
  (( P_MAX > P_MAX_BY_T )) && P_MAX=$P_MAX_BY_T

  # Seed with P=1, T=TOTAL_CPUS (always valid; tiebreak prefers larger T).
  PARALLEL_PDBS=1
  THREADS_PER_PDB=$TOTAL_CPUS
  BEST_WORK=$(( NUM_PDBS ))   # ceil(N/1) * 1

  for (( _p = 2; _p <= P_MAX; _p++ )); do
    _t=$(( TOTAL_CPUS / _p ))
    (( _t < MIN_FOLDX_JOBS_PER_PDB )) && continue
    _waves=$(( (NUM_PDBS + _p - 1) / _p ))
    _work=$(( _waves * _p ))
    # Strictly less work wins; ties keep the larger-T (smaller-P) seed.
    if (( _work < BEST_WORK )); then
      PARALLEL_PDBS=$_p; THREADS_PER_PDB=$_t; BEST_WORK=$_work
    fi
  done
  WAVES=$(( (NUM_PDBS + PARALLEL_PDBS - 1) / PARALLEL_PDBS ))
  log_info "Auto: ${TOTAL_CPUS} CPUs, ${NUM_PDBS} PDBs → ${PARALLEL_PDBS} parallel × ${THREADS_PER_PDB} FoldX jobs (${WAVES} waves)"
fi

# Per-PDB timeout to prevent infinite hangs (0 = no timeout)
PER_PDB_TIMEOUT_SECS=$(( PER_PDB_TIMEOUT_HOURS * 3600 ))
if (( PER_PDB_TIMEOUT_SECS > 0 )); then
  log_info "Per-PDB timeout: ${PER_PDB_TIMEOUT_HOURS}h"
else
  log_info "Per-PDB timeout: disabled"
fi

log_info "Runs per mutation: $NRUNS"

# Dump FoldX configuration for troubleshooting
if [[ "$DEBUG_MODE" == "true" ]]; then
  dump_foldx_info
  _dump_section "PARALLELISM CONFIG" \
    "Compute profile: $MACHINE" \
    "Total CPUs: $TOTAL_CPUS" \
    "FoldX optimal threads: $FOLDX_OPTIMAL_THREADS (single-threaded)" \
    "Parallel PDBs: $PARALLEL_PDBS" \
    "FoldX jobs per PDB: $THREADS_PER_PDB" \
    "Effective CPU usage: $(( PARALLEL_PDBS * THREADS_PER_PDB )) / $TOTAL_CPUS" \
    "Per-PDB timeout: ${PER_PDB_TIMEOUT_HOURS}h (${PER_PDB_TIMEOUT_SECS}s)" \
    "Num PDBs to process: $NUM_PDBS"
  log_debug "Configuration validated. Starting processing."
fi

# ==============================================================================
# PARALLEL PDB PROCESSING
# ==============================================================================

# Status tracking via temp directory (works across subshells, unlike arrays)
STATUS_DIR=$(mktemp -d)

# Slot limiter: wait until a background slot is free (parameterized).
wait_for_slot() {
  local limit="$1"
  while (( $(jobs -rp | wc -l) >= limit )); do
    sleep 0.5
  done
}

# Process a single PDB file — designed to run as a background job via &
# Usage: process_single_pdb <pdb_path> <threads_per_pdb>
process_single_pdb() {
  local pdb="$1"
  local threads="$2"
  local base_name="${pdb##*/}"; base_name="${base_name%.pdb}"

  # FoldX filename safety: replace dots with underscores to prevent FoldX truncation
  local safe_base_name="${base_name//./_}"
  local pdb_to_use="$pdb"

  if [[ "$safe_base_name" != "$base_name" ]]; then
    log_warn "[$base_name] PDB filename contains dots → safe copy: ${safe_base_name}.pdb"
    local safe_pdb_dir="$RESULTS_DIR/${base_name}/safe_input"
    mkdir -p "$safe_pdb_dir"
    local safe_pdb="$safe_pdb_dir/${safe_base_name}.pdb"
    cp "$pdb" "$safe_pdb"
    pdb_to_use="$safe_pdb"
  fi

  local outdir="$RESULTS_DIR/${base_name}/nruns_${NRUNS}"
  local completion_marker="$outdir/.completed"

  # Handle existing runs: skip or clean up
  if [[ -d "$outdir" ]]; then
    if [[ "$FORCE_RERUN" == "true" ]]; then
      log_warn "[$base_name] FORCE_RERUN → removing previous run"
      rm -rf "$outdir"
    else
      if [[ -f "$completion_marker" ]]; then
        log_warn "[$base_name] Skipping (already completed)"
      else
        log_warn "[$base_name] Incomplete run. Set FORCE_RERUN=true to restart."
      fi
      echo "SKIPPED" > "$STATUS_DIR/${base_name}.status"
      return 0
    fi
  fi

  log_step "[$base_name] Processing ($threads parallel FoldX jobs)"

  if [[ "$DEBUG_MODE" == "true" ]]; then
    dump_pdb_info "$pdb_to_use"
  fi

  mkdir -p "$outdir"
  local pdb_log="$CONSOLE_LOG_DIR/console_${RUN_ID}_${base_name}.log"

  log_mutatex_invocation "$pdb_to_use" "$outdir" "$threads"

  # Start background cleanup worker (prevents 196GB disk usage per PDB)
  local cleanup_pid=""
  if [[ "$CLEANUP_ENABLED" == "true" ]]; then
    log_info "[$base_name] Cleanup worker started (interval=${CLEANUP_WORKER_INTERVAL}s)"
    cleanup_pid=$(start_cleanup_worker "$outdir/mutations" "$CLEANUP_WORKER_INTERVAL")
  fi

  # Run mutatex from the output directory (mutatex operates in cwd)
  # Each background job gets its own subshell for cd, so it won't affect others
  local exit_status=0
  if (( PER_PDB_TIMEOUT_SECS > 0 )); then
    ( cd "$outdir" && exec timeout --signal=TERM --kill-after=120 "$PER_PDB_TIMEOUT_SECS" \
        mutatex "$pdb_to_use" \
          -p "$threads" \
          -n "$NRUNS" \
          -x "$FOLDX_BINARY" \
          -b "$ROTABASE" \
          -f suite5 \
          -R "$REPAIR_TEMPLATE" \
          -M "$MUTATE_TEMPLATE" \
          -I "$INTERFACE_TEMPLATE" \
          -B -v -l ) > "$pdb_log" 2>&1
    exit_status=$?
  else
    ( cd "$outdir" && mutatex "$pdb_to_use" \
        -p "$threads" \
        -n "$NRUNS" \
        -x "$FOLDX_BINARY" \
        -b "$ROTABASE" \
        -f suite5 \
        -R "$REPAIR_TEMPLATE" \
        -M "$MUTATE_TEMPLATE" \
        -I "$INTERFACE_TEMPLATE" \
        -B -v -l ) > "$pdb_log" 2>&1
    exit_status=$?
  fi

  # Stop cleanup worker (must happen even on failure/timeout)
  if [[ -n "$cleanup_pid" ]]; then
    stop_cleanup_worker "$cleanup_pid"
  fi

  log_run_result "$pdb_to_use" "$exit_status" "$outdir"

  if [[ $exit_status -eq 0 ]]; then
    log_info "[$base_name] MutateX completed successfully"

    # Generate PyMOL visualization using ddg2pdb
    generate_pymol_visualization "$outdir" "$base_name" >> "$pdb_log" 2>&1

    touch "$completion_marker"
    echo "SUCCESS" > "$STATUS_DIR/${base_name}.status"
  elif [[ $exit_status -eq 124 ]]; then
    log_error "[$base_name] TIMED OUT after ${PER_PDB_TIMEOUT_HOURS}h — killed"
    rm -f "$completion_marker"
    echo "FAILED" > "$STATUS_DIR/${base_name}.status"
    _dump_section "TIMEOUT for $base_name" \
      "Timeout: ${PER_PDB_TIMEOUT_HOURS}h" \
      "Last 30 lines:" "$(tail -30 "$pdb_log" 2>/dev/null || echo 'no log')"
  else
    log_error "[$base_name] MutateX failed (exit $exit_status)"
    rm -f "$completion_marker"
    echo "FAILED" > "$STATUS_DIR/${base_name}.status"
    _dump_section "FAILURE CONTEXT for $base_name" \
      "Exit code: $exit_status" \
      "Last 30 lines:" "$(tail -30 "$pdb_log" 2>/dev/null || echo 'no log')"
  fi

  # Post-run cleanup to save disk space (196GB+ per run without this)
  if [[ "$CLEANUP_ENABLED" == "true" && -d "$outdir/mutations" ]]; then
    log_info "[$base_name] Cleaning up intermediate files..."
    cleanup_mutations_dir "$outdir/mutations"
  fi

  log_step "[$base_name] Done"
}

# ==============================================================================
# DISPATCH: Process all PDBs in parallel (capped at PARALLEL_PDBS)
# ==============================================================================
printf -v START_TIME '%(%s)T' -1

log_step "Dispatching ${NUM_PDBS} PDBs (max ${PARALLEL_PDBS} in parallel, ${THREADS_PER_PDB} FoldX jobs each)"

for pdb in "${PDB_FILES[@]}"; do
  wait_for_slot "$PARALLEL_PDBS"
  process_single_pdb "$pdb" "$THREADS_PER_PDB" &
done

# Wait for all PDB jobs to finish
wait

log_step "All PDB jobs completed"

# ==============================================================================
# COLLECT RESULTS FROM STATUS FILES
# ==============================================================================
declare -a PROCESSED_PDBS=()
declare -a FAILED_PDBS=()
declare -a SKIPPED_PDBS=()

for status_file in "$STATUS_DIR"/*.status; do
  [[ -f "$status_file" ]] || continue
  local_base="${status_file##*/}"; local_base="${local_base%.status}"
  case "$(< "$status_file")" in
    SUCCESS) PROCESSED_PDBS+=("$local_base") ;;
    FAILED)  FAILED_PDBS+=("$local_base") ;;
    SKIPPED) SKIPPED_PDBS+=("$local_base") ;;
  esac
done

log_step "Batch processing complete"

# Calculate total runtime (printf builtin avoids date subprocess)
printf -v END_TIME '%(%s)T' -1
TOTAL_RUNTIME=$((END_TIME - START_TIME))
RUNTIME_HOURS=$((TOTAL_RUNTIME / 3600))
RUNTIME_MINS=$(((TOTAL_RUNTIME % 3600) / 60))
RUNTIME_SECS=$((TOTAL_RUNTIME % 60))

# Generate troubleshooting summary
summary_lines=(
  "Total runtime: ${RUNTIME_HOURS}h ${RUNTIME_MINS}m ${RUNTIME_SECS}s"
  "Total PDBs configured: ${#PDB_FILES[@]}"
  "Parallel PDBs: $PARALLEL_PDBS"
  "FoldX jobs per PDB: $THREADS_PER_PDB"
  "Successfully processed: ${#PROCESSED_PDBS[@]}"
  "Failed: ${#FAILED_PDBS[@]}"
  "Skipped (already complete): ${#SKIPPED_PDBS[@]}"
)
[[ ${#PROCESSED_PDBS[@]} -gt 0 ]] && summary_lines+=("Processed PDBs: ${PROCESSED_PDBS[*]}")
[[ ${#FAILED_PDBS[@]} -gt 0 ]]   && summary_lines+=("Failed PDBs: ${FAILED_PDBS[*]}")
[[ ${#SKIPPED_PDBS[@]} -gt 0 ]]  && summary_lines+=("Skipped PDBs: ${SKIPPED_PDBS[*]}")
summary_lines+=(
  "Final disk space: $(df -h "$SCRIPT_DIR" 2>/dev/null | tail -1 || echo 'unknown')"
  "Results directory size: $(du -sh "$RESULTS_DIR" 2>/dev/null | cut -f1 || echo 'unknown')"
)
_dump_section "MUTATEX PROCESSING SUMMARY" "${summary_lines[@]}"

# Print summary to console
log_step "Processing Summary"
log_info "Total runtime: ${RUNTIME_HOURS}h ${RUNTIME_MINS}m ${RUNTIME_SECS}s"
log_info "Processed: ${#PROCESSED_PDBS[@]} | Failed: ${#FAILED_PDBS[@]} | Skipped: ${#SKIPPED_PDBS[@]}"
log_info "Parallelism: ${PARALLEL_PDBS} PDBs × ${THREADS_PER_PDB} FoldX jobs (${TOTAL_CPUS} CPUs detected)"
if [[ ${#FAILED_PDBS[@]} -gt 0 ]]; then
  log_warn "Failed PDBs: ${FAILED_PDBS[*]}"
  log_warn "Check troubleshooting log for details: $TROUBLESHOOT_LOG"
fi

# ==============================================================================
# POST-PROCESSING: Run extraction + visualization in parallel
# (Independent outputs: critical_residues/ vs visualizations/ — safe to overlap)
# Both scripts run simultaneously; each gets internal parallelism workers.
# CPU budget: TOTAL_CPUS / 2 per script (two scripts sharing the machine).
# ==============================================================================

# Auto-determine post-processing workers (or use TOML override)
if (( PP_WORKERS_CFG > 0 )); then
  PP_WORKERS=$PP_WORKERS_CFG
else
  PP_WORKERS=$(( TOTAL_CPUS / 2 ))
  (( PP_WORKERS < 1 )) && PP_WORKERS=1
fi
log_info "Post-processing: 2 scripts in parallel × ${PP_WORKERS} workers each"

SCRIPTS_DIR="$SCRIPT_DIR/modules/PPI/scripts"

python "$SCRIPTS_DIR/c_mutatex_extract_all_critical_residues.py" \
  --run-results "$RESULTS_DIR" \
  --inputs "$INPUT_BASE" \
  --workers "$PP_WORKERS" >> "$TROUBLESHOOT_LOG" 2>&1 &
PID_EXTRACT=$!

MUTATEX_WORKERS=$PP_WORKERS python "$SCRIPTS_DIR/b_mutatex_generate_mutatex_visualizations.py" \
  "$RESULTS_DIR" >> "$TROUBLESHOOT_LOG" 2>&1 &
PID_VIS=$!

wait "$PID_EXTRACT"; PYTHON_EXIT1=$?
log_debug "c_mutatex_extract_all_critical_residues.py exit code: $PYTHON_EXIT1"
[[ $PYTHON_EXIT1 -ne 0 ]] && log_warn "Critical residue extraction failed (exit $PYTHON_EXIT1)"

wait "$PID_VIS"; PYTHON_EXIT2=$?
log_debug "b_mutatex_generate_mutatex_visualizations.py exit code: $PYTHON_EXIT2"
[[ $PYTHON_EXIT2 -ne 0 ]] && log_warn "Visualization generation failed (exit $PYTHON_EXIT2)"

log_step "All logs saved"
log_info "Full log: $LOG_FILE"
log_info "Console log (full output): $CONSOLE_LOG"
log_info "Per-PDB logs: $CONSOLE_LOG_DIR/console_${RUN_ID}_*.log"
log_info "Troubleshooting log: $TROUBLESHOOT_LOG"
