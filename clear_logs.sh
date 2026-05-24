#!/bin/bash
# clear_logs.sh — Remove contents of logs/ (z_logs_prompt/ always preserved).
#
# Modes:
#   (default)        Remove everything under logs/ except z_logs_prompt/.
#   --errors-only    Remove only logs from runs that contain ERROR entries
#                    (keeps clean runs). A "run" is identified by its
#                    pipeline_YYYYMMDD_HHMMSS prefix in error_warn_logs/.
#   -n, --dry-run    Print what would be deleted; do not delete or prompt.
#   -y, --yes        Skip the confirmation prompt.
#   -h, --help       Show this help.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGS_DIR="$SCRIPT_DIR/logs"
KEEP_DIR="z_logs_prompt"

# --- User-configurable defaults (edit here to change behavior without flags) ---
MODE="errors"       # "all" = clear everything except z_logs_prompt/
                 # "errors" = remove only error-run files (same as --errors-only)
DRY_RUN=false    # true = print what would be deleted, no changes (same as -n)
ASSUME_YES=false # true = skip confirmation prompt (same as -y)

usage() {
    sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --errors-only) MODE="errors" ;;
        -n|--dry-run)  DRY_RUN=true ;;
        -y|--yes)      ASSUME_YES=true ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
    shift
done

if [[ ! -d "$LOGS_DIR" ]]; then
    echo "ERROR: logs directory not found: $LOGS_DIR" >&2
    exit 1
fi

confirm_or_exit() {
    $ASSUME_YES && return 0
    read -r -p "Proceed? [y/N] " confirm
    if [[ "${confirm,,}" != "y" ]]; then
        echo "Aborted."
        exit 0
    fi
}

if [[ "$MODE" == "all" ]]; then
    echo "The following will be deleted from: $LOGS_DIR"
    echo "--------------------------------------------"
    find "$LOGS_DIR" -maxdepth 1 ! -name "$(basename "$LOGS_DIR")" ! -name "$KEEP_DIR" -print
    echo "--------------------------------------------"
    echo "Directory preserved: $KEEP_DIR"
    echo ""
    $DRY_RUN && { echo "[dry-run] no changes made."; exit 0; }
    confirm_or_exit
    find "$LOGS_DIR" -maxdepth 1 ! -name "$(basename "$LOGS_DIR")" ! -name "$KEEP_DIR" -exec rm -rf {} +
    echo "Done. logs/ cleared (z_logs_prompt preserved)."
    exit 0
fi

# --- --errors-only mode --------------------------------------------------
# Identify runs with ERRORs by scanning every error_warn_logs* directory.
# A run is keyed by its 'pipeline_YYYYMMDD_HHMMSS' prefix.

mapfile -t ERROR_FILES < <(
    find "$LOGS_DIR" -maxdepth 2 -type f \
        -path "*/error_warn_logs*/pipeline_*_errors_warnings.log" 2>/dev/null
)

ERROR_PREFIXES=()
for f in "${ERROR_FILES[@]}"; do
    # A run had errors if its error_warn log contains at least one ERROR line.
    if grep -q -E '\bERROR\b' "$f" 2>/dev/null; then
        base="$(basename "$f")"
        prefix="${base%_errors_warnings.log}"   # -> pipeline_YYYYMMDD_HHMMSS
        ERROR_PREFIXES+=("$prefix")
    fi
done

if [[ ${#ERROR_PREFIXES[@]} -eq 0 ]]; then
    echo "No runs with ERROR entries found. Nothing to delete."
    exit 0
fi

# De-duplicate prefixes.
mapfile -t ERROR_PREFIXES < <(printf '%s\n' "${ERROR_PREFIXES[@]}" | sort -u)

echo "Runs with errors (${#ERROR_PREFIXES[@]}):"
printf '  %s\n' "${ERROR_PREFIXES[@]}"
echo ""

# Collect every file across logs/ (excluding z_logs_prompt/) whose name
# starts with one of the error prefixes.
TO_DELETE=()
for prefix in "${ERROR_PREFIXES[@]}"; do
    while IFS= read -r -d '' match; do
        TO_DELETE+=("$match")
    done < <(
        find "$LOGS_DIR" -mindepth 2 -type f \
            -name "${prefix}*" \
            -not -path "$LOGS_DIR/$KEEP_DIR/*" \
            -print0 2>/dev/null
    )
done

if [[ ${#TO_DELETE[@]} -eq 0 ]]; then
    echo "No matching files to delete."
    exit 0
fi

echo "Files to delete (${#TO_DELETE[@]}):"
echo "--------------------------------------------"
printf '%s\n' "${TO_DELETE[@]}"
echo "--------------------------------------------"
echo "Directory preserved: $KEEP_DIR"
echo "Clean runs (no ERROR entries) preserved."
echo ""

$DRY_RUN && { echo "[dry-run] no changes made."; exit 0; }
confirm_or_exit

for f in "${TO_DELETE[@]}"; do
    rm -f "$f"
done

echo "Done. Removed ${#TO_DELETE[@]} files from ${#ERROR_PREFIXES[@]} error-containing runs."
