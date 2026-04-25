#!/bin/bash
# ==============================================================================
# enforce_fasta_format.sh — Scan and enforce FASTA formatting rules
# ==============================================================================
# Rules enforced:
#   1. Sequence lines wrapped at 80 characters max
#   2. Exactly one blank line between the last sequence line and the next header
#   3. No leading blank lines before the first header
#   4. No trailing blank lines after the last sequence
#
# Usage:
#   bash enforce_fasta_format.sh [--fix] [--dir <path>] [--dry-run] [--quiet]
#                                [--max-size <MB>] [--include-refs]
#                                [--parallel <N>]
#
#   --fix           Rewrite non-compliant files in-place (default: scan only)
#   --dir  <path>   Root directory to scan (default: current script directory)
#   --dry-run       Show what --fix would do without writing
#   --quiet         Only print files with violations (suppress OK files)
#   --max-size <MB> Skip files larger than this many MB (default: 100)
#   --include-refs  Also scan I_RefSeqs/ and II_INPUTS/ (skipped by default)
#   --parallel <N>  Max parallel file checks (default: 12)
# ==============================================================================
set -euo pipefail

# ── Configurable variables ────────────────────────────────────────────────────
LINE_WIDTH=80
FASTA_EXTENSIONS=("fa" "fasta" "faa" "fna" "fas" "pep")
SKIP_DIRS=("z_archive" ".git" "node_modules" "__pycache__")
IMMUTABLE_DIRS=("I_RefSeqs" "II_INPUTS")
MAX_SIZE_MB=100
MAX_PARALLEL=12
STAGE_THRESHOLD_MB=5   # Files larger than this are staged to /tmp before processing
# ──────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR"
PROJECT_ROOT="$PIPELINE_DIR"
MODULES="$PIPELINE_DIR/modules"

source "$MODULES/logging/logging_utils.sh"

FIX=false
DRY_RUN=false
QUIET=false
INCLUDE_REFS=false
SCAN_DIR="$SCRIPT_DIR"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fix)           FIX=true;           shift ;;
        --dry-run)       DRY_RUN=true;       shift ;;
        --quiet)         QUIET=true;         shift ;;
        --include-refs)  INCLUDE_REFS=true;  shift ;;
        --dir)           SCAN_DIR="$2";      shift 2 ;;
        --max-size)      MAX_SIZE_MB="$2";   shift 2 ;;
        --parallel)      MAX_PARALLEL="$2";  shift 2 ;;
        -h|--help)
            sed -n '2,/^# ====/{ /^# ====/d; s/^# \?//p }' "$0"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Build find command ───────────────────────────────────────────────────────
find_fasta_files() {
    local prune_args=()
    for skip in "${SKIP_DIRS[@]}"; do
        prune_args+=(-name "$skip" -o)
    done
    if ! $INCLUDE_REFS; then
        for immutable in "${IMMUTABLE_DIRS[@]}"; do
            prune_args+=(-name "$immutable" -o)
        done
    fi
    # Remove trailing -o
    unset 'prune_args[-1]'

    local ext_args=()
    local first=true
    for ext in "${FASTA_EXTENSIONS[@]}"; do
        if $first; then
            ext_args+=(-name "*.${ext}")
            first=false
        else
            ext_args+=(-o -name "*.${ext}")
        fi
    done

    local max_size_bytes
    max_size_bytes=$(( MAX_SIZE_MB * 1024 * 1024 ))

    find "$SCAN_DIR" \
        \( "${prune_args[@]}" \) -prune \
        -o \( "${ext_args[@]}" \) -type f -size -"${max_size_bytes}"c -print
}

# ── AWK: check a FASTA file for violations (fast, no rewrite) ────────────────
check_fasta_awk() {
    awk -v width="$LINE_WIDTH" '
    BEGIN {
        violations = 0
        blank_count = 0
        found_first_header = 0
    }
    /^[[:space:]]*$/ {
        blank_count++
        next
    }
    /^>/ {
        if (!found_first_header) {
            found_first_header = 1
            if (blank_count > 0) {
                print "Line " NR ": leading blank lines before first header"
                violations++
            }
        } else {
            if (blank_count == 0) {
                print "Line " NR ": missing blank line before header"
                violations++
            } else if (blank_count > 1) {
                print "Line " NR ": " blank_count " blank lines before header (expected 1)"
                violations++
            }
        }
        blank_count = 0
        next
    }
    {
        # Sequence line
        if (blank_count > 0 && found_first_header) {
            print "Line " NR ": unexpected blank line(s) within sequence block"
            violations++
        }
        blank_count = 0
        if (length($0) > width) {
            print "Line " NR ": sequence line " length($0) " chars (max " width ")"
            violations++
        }
    }
    END {
        if (blank_count > 0 && found_first_header) {
            print "EOF: trailing blank line(s) after last sequence"
            violations++
        }
        exit (violations > 0) ? 1 : 0
    }
    ' "$1"
}

# ── AWK: reformat a FASTA file in-place ─────────────────────────────────────
# Accepts an optional second argument for the source file to read from
# (used when input has been staged to /tmp). Output always replaces $1.
format_fasta() {
    local input_file="$1"
    local source_file="${2:-$input_file}"
    local tmp_file
    tmp_file=$(mktemp "${TMPDIR:-/tmp}/fasta_fmt_out_XXXXXX")

    awk -v width="$LINE_WIDTH" '
    BEGIN { seq = ""; header = ""; first = 1 }
    /^>/ {
        if (header != "") {
            if (!first) printf "\n"
            print header
            if (seq != "") {
                for (i = 1; i <= length(seq); i += width)
                    print substr(seq, i, width)
            }
            first = 0
        }
        header = $0
        seq = ""
        next
    }
    /^[[:space:]]*$/ { next }
    {
        gsub(/[[:space:]]/, "")
        seq = seq $0
    }
    END {
        if (header != "") {
            if (!first) printf "\n"
            print header
            if (seq != "")
                for (i = 1; i <= length(seq); i += width)
                    print substr(seq, i, width)
        }
    }
    ' "$source_file" > "$tmp_file" \
        && cp "$tmp_file" "$input_file" \
        && rm -f "$tmp_file" \
        || { rm -f "$tmp_file"; return 1; }
}

# ── Parallelization helpers ──────────────────────────────────────────────────
wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.1; done
}

RESULT_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fasta_fmt_XXXXXX")
cleanup() {
    # Flush logging FIRST so tee can drain its pipe buffer before we kill jobs.
    teardown_logging 2>/dev/null
    local pids
    pids=$(jobs -rp 2>/dev/null) && [[ -n "$pids" ]] && kill $pids 2>/dev/null || true
    wait 2>/dev/null || true
    rm -rf "$RESULT_DIR"
}
trap cleanup EXIT

# ── Stage large files to /tmp for faster I/O on slow mounts (WSL /mnt/c) ─────
stage_threshold_bytes=$(( STAGE_THRESHOLD_MB * 1024 * 1024 ))

stage_file_if_large() {
    local file="$1"
    local file_size
    file_size=$(stat -c%s "$file" 2>/dev/null || echo 0)
    if (( file_size > stage_threshold_bytes )); then
        local staged
        staged=$(mktemp "${TMPDIR:-/tmp}/fasta_stage_XXXXXX")
        cp "$file" "$staged"
        echo "$staged"
    fi
}

process_file() {
    local file="$1"
    local idx="$2"
    local rel_path="${file#"${SCAN_DIR}/"}"
    local out_file="$RESULT_DIR/${idx}.out"
    local stat_file="$RESULT_DIR/${idx}.stat"

    # Stage large files to /tmp to avoid slow I/O on network/WSL mounts
    local staged_file=""
    staged_file=$(stage_file_if_large "$file")
    local work_file="${staged_file:-$file}"

    if $FIX && ! $DRY_RUN; then
        # ── Single-pass: reformat unconditionally, compare checksums ──────
        local orig_hash
        orig_hash=$(md5sum "$work_file" | awk '{print $1}')

        format_fasta "$file" "$work_file"   # reads work_file, writes to file

        local new_hash
        new_hash=$(md5sum "$file" | awk '{print $1}')

        if [[ "$orig_hash" != "$new_hash" ]]; then
            # Something changed — run checker on the ORIGINAL to report what was wrong
            set +e
            local violations_output
            violations_output=$(check_fasta_awk "$work_file" 2>&1)
            set -e
            local violation_count=0
            [[ -n "$violations_output" ]] && violation_count=$(echo "$violations_output" | wc -l)
            {
                log_warn "FAIL [${violation_count}] ${rel_path}"
                while IFS= read -r v; do
                    log_warn "  → $v"
                done <<< "$violations_output"
                log_info "  ✓ Fixed: ${rel_path}"
            } > "$out_file" 2>&1
            echo "FAIL ${violation_count} true" > "$stat_file"
        else
            if ! $QUIET; then
                log_info "OK   ${rel_path}" > "$out_file" 2>&1
            else
                : > "$out_file"
            fi
            echo "OK 0 false" > "$stat_file"
        fi
    else
        # ── Scan-only or dry-run: just check ─────────────────────────────
        set +e
        local violations_output
        violations_output=$(check_fasta_awk "$work_file" 2>&1)
        local awk_exit=$?
        set -e

        local violation_count=0
        if [[ -n "$violations_output" ]]; then
            violation_count=$(echo "$violations_output" | wc -l)
        fi

        if (( awk_exit != 0 )) && (( violation_count > 0 )); then
            {
                log_warn "FAIL [${violation_count}] ${rel_path}"
                while IFS= read -r v; do
                    log_warn "  → $v"
                done <<< "$violations_output"
                if $DRY_RUN; then
                    log_info "  (would fix): ${rel_path}"
                fi
            } > "$out_file" 2>&1
            echo "FAIL ${violation_count} false" > "$stat_file"
        else
            if ! $QUIET; then
                log_info "OK   ${rel_path}" > "$out_file" 2>&1
            else
                : > "$out_file"
            fi
            echo "OK 0 false" > "$stat_file"
        fi
    fi

    # Clean up staged copy
    [[ -n "$staged_file" ]] && rm -f "$staged_file"
}

# ── Main ─────────────────────────────────────────────────────────────────────

# ── Auto-tune parallelism for slow I/O mounts ───────────────────────────────
# WSL /mnt/ paths use 9p filesystem with poor parallel I/O; cap at 4.
if [[ "$SCAN_DIR" == /mnt/* ]] && (( MAX_PARALLEL > 4 )); then
    MAX_PARALLEL=4
fi

setup_logging
log_step "FASTA Format Enforcement"
log_info "Scan directory : ${SCAN_DIR}"
log_info "Mode           : $(if $FIX; then echo 'FIX'; elif $DRY_RUN; then echo 'DRY-RUN'; else echo 'SCAN'; fi)"
log_info "Line width     : ${LINE_WIDTH}"
log_info "Max file size  : ${MAX_SIZE_MB} MB"
log_info "Include refs   : $($INCLUDE_REFS && echo 'yes' || echo 'no (I_RefSeqs/, II_INPUTS/ skipped)')"
log_info "Stage threshold: ${STAGE_THRESHOLD_MB} MB (files above this staged to /tmp)"
log_info "Parallel jobs  : ${MAX_PARALLEL}$(if [[ "$SCAN_DIR" == /mnt/* ]]; then echo ' (WSL auto-capped)'; fi)"

mapfile -t fasta_files < <(find_fasta_files 2>/dev/null | sort)

if (( ${#fasta_files[@]} == 0 )); then
    log_warn "No FASTA files found."
    teardown_logging
    exit 0
fi

log_info "Found ${#fasta_files[@]} FASTA file(s) to check (max ${MAX_PARALLEL} parallel)."
echo ""

# Launch parallel file checks (collect PIDs so we wait only for these,
# not the background tee process spawned by setup_logging — bare `wait`
# would deadlock because tee blocks on its stdin pipe).
_job_pids=()
for idx in "${!fasta_files[@]}"; do
    wait_for_slot
    process_file "${fasta_files[$idx]}" "$idx" &
    _job_pids+=($!)
done
wait "${_job_pids[@]}" 2>/dev/null || true

# Aggregate results (preserving original file order for deterministic output)
total_files=${#fasta_files[@]}
violated_files=0
fixed_files=0
total_violations=0

for (( i=0; i<total_files; i++ )); do
    [[ -s "$RESULT_DIR/${i}.out" ]] && cat "$RESULT_DIR/${i}.out"
    if [[ -f "$RESULT_DIR/${i}.stat" ]]; then
        read -r status vcount was_fixed < "$RESULT_DIR/${i}.stat"
        if [[ "$status" == "FAIL" ]]; then
            (( violated_files++ )) || true
            (( total_violations += vcount )) || true
            [[ "$was_fixed" == "true" ]] && (( fixed_files++ )) || true
        fi
    fi
done

# ── Summary ──────────────────────────────────────────────────────────────────
log_step "Summary"
log_info "Files scanned     : ${total_files}"
log_info "Files with issues : ${violated_files}"
log_info "Total violations  : ${total_violations}"
if $FIX && ! $DRY_RUN; then
    log_info "Files fixed       : ${fixed_files}"
fi

if (( violated_files > 0 )) && ! $FIX; then
    log_warn "Run with --fix to auto-correct all violations."
    log_warn "Use --include-refs to also check I_RefSeqs/ and II_INPUTS/."
    teardown_logging
    exit 1
fi

teardown_logging
exit 0
