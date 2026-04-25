#!/bin/bash
# Module: Merge BLAST CSV Results
# Usage: bash merge_blast_csv.sh --input-dir <dir> --output <file.csv> [--include-pattern <pat>] [--exclude-pattern <pat>]
#
# --include-pattern  Only merge CSVs whose filename matches this grep pattern
# --exclude-pattern  Exclude CSVs whose filename matches this grep pattern
# (both patterns are applied to the basename of each CSV file)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

INPUT_DIR=""
OUTPUT_FILE=""
INCLUDE_PATTERN=""
EXCLUDE_PATTERN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)        INPUT_DIR="$2";        shift 2 ;;
        --output)           OUTPUT_FILE="$2";       shift 2 ;;
        --include-pattern)  INCLUDE_PATTERN="$2";   shift 2 ;;
        --exclude-pattern)  EXCLUDE_PATTERN="$2";   shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" ]]   && { log_error "Missing --input-dir"; exit 1; }
[[ -z "$OUTPUT_FILE" ]] && { log_error "Missing --output";    exit 1; }

mkdir -p "$(dirname "$OUTPUT_FILE")"

mapfile -t all_csv_files < <(find "$INPUT_DIR" -type f -name "*.csv" 2>/dev/null | sort)

# Apply include/exclude filters on filename basenames
csv_files=()
for f in "${all_csv_files[@]}"; do
    bn=$(basename "$f")
    [[ -n "$INCLUDE_PATTERN" ]] && ! echo "$bn" | grep -q "$INCLUDE_PATTERN" && continue
    [[ -n "$EXCLUDE_PATTERN" ]] &&   echo "$bn" | grep -q "$EXCLUDE_PATTERN" && continue
    csv_files+=("$f")
done

if [[ ${#csv_files[@]} -eq 0 ]]; then
    log_info "No CSV files matched in $INPUT_DIR (include='$INCLUDE_PATTERN' exclude='$EXCLUDE_PATTERN')"
    exit 0
fi

# Merge: keep header from first file, skip headers from rest
# Then sort by Subject ID (col 1), E-value ascending (col 3), Percent Identity descending (col 4)
{
    awk 'NR==1' "${csv_files[@]}"
    awk 'FNR>1' "${csv_files[@]}" | sort -t',' -k1,1 -k3,3g -k4,4gr | uniq
} > "$OUTPUT_FILE"

log_info "Merged ${#csv_files[@]} CSV(s) -> $(basename "$OUTPUT_FILE") (sorted by Subject ID, E-value, %Identity)"
