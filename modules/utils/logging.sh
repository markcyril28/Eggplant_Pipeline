#!/bin/bash
set -euo pipefail
# ============================================================================
# Common Logging Utilities for Eggplant Pipeline
# Source this file in any module or orchestrator script.
#
# Usage:
#   source "$(dirname "${BASH_SOURCE[0]}")/modules/utils/logging.sh"
#   setup_logging "my_pipeline"
#   log_info "Starting work..."
# ============================================================================

# Colors
readonly _RED='\033[0;31m'
readonly _GREEN='\033[0;32m'
readonly _YELLOW='\033[1;33m'
readonly _BLUE='\033[0;34m'
readonly _PURPLE='\033[0;35m'
readonly _CYAN='\033[0;36m'
readonly _WHITE='\033[1;37m'
readonly _NC='\033[0m'  # No Color

_LOG_FILE=""

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# Colored logging functions
log() {
    local lvl="$1"
    shift
    case "$lvl" in
        "ERROR")   printf "${_RED}[%s] [%s]${_NC} %s\n" "$(timestamp)" "$lvl" "$*"; ;;
        "WARN")    printf "${_YELLOW}[%s] [%s]${_NC} %s\n" "$(timestamp)" "$lvl" "$*"; ;;
        "INFO")    printf "${_GREEN}[%s] [%s]${_NC} %s\n" "$(timestamp)" "$lvl" "$*"; ;;
        "STEP")    printf "${_CYAN}[%s] [%s]${_NC} %s\n" "$(timestamp)" "$lvl" "$*"; ;;
        *)         printf "${_WHITE}[%s] [%s]${_NC} %s\n" "$(timestamp)" "$lvl" "$*"; ;;
    esac
}

log_info()  { log INFO  "$@"; }
log_warn()  { log WARN  "$@"; }
log_error() { log ERROR "$@"; }
log_step()  { log STEP "============================== $* =============================="; }

setup_logging() {
    local name="${1:-pipeline}"
    # Default logs to the pipeline root when available so logs are not written
    # to an unexpected caller working directory.
    local default_log_dir="${PIPELINE_DIR:-$PWD}/logs"
    local log_dir="${2:-$default_log_dir}"
    local run_id
    run_id="$(date +%Y%m%d_%H%M%S)"

    mkdir -p "$log_dir"
    _LOG_FILE="$log_dir/${name}_${run_id}.log"

    # Save original stdout/stderr so we can restore them on exit.
    # Restoring breaks the exec > >(...) pipe cleanly and lets tee
    # receive EOF, preventing the script from hanging after 'wait'.
    exec 3>&1 4>&2
    exec > >(tee -a "$_LOG_FILE") 2>&1
    _LOGGING_SETUP=true
    log_info "Log file: $_LOG_FILE"
}

teardown_logging() {
    if [[ "${_LOGGING_SETUP:-}" == "true" ]]; then
        exec 1>&3 2>&4 3>&- 4>&-
        _LOGGING_SETUP=false
        # Give tee a moment to flush remaining output
        sleep 0.2
    fi
}

run_timed() {
    local start end elapsed
    start=$(date +%s)
    "$@"
    local rc=$?
    end=$(date +%s)
    elapsed=$((end - start))
    log_info "Elapsed: $(date -u -d @${elapsed} +%H:%M:%S 2>/dev/null || printf '%02d:%02d:%02d' $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60)))"
    return $rc
}
