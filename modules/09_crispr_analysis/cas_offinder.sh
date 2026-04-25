#!/bin/bash
# Module: CRISPR Off-Target Analysis (Cas-OFFinder)
# Usage: bash cas_offinder.sh --grna-fasta <fasta> --genome <ref.fa> --outdir <dir> \
#            [--pam "NNNNNNNNNNNNNNNNNNNNNGG"] [--max-mismatches 3] [--device G0]
#
# Converts gRNA FASTA to Cas-OFFinder input format, runs cas-offinder, outputs results.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../logging/logging_utils.sh" 2>/dev/null || true

# Defaults
GRNA_FASTA=""
GENOME=""
OUTPUT_DIR="."
PAM="NNNNNNNNNNNNNNNNNNNNNGG"   # SpCas9 NGG PAM (23-nt pattern)
MAX_MISMATCHES=3
DEVICE="G0"                      # G0 = first GPU, C = CPU-only

while [[ $# -gt 0 ]]; do
    case "$1" in
        --grna-fasta)      GRNA_FASTA="$2"; shift 2 ;;
        --genome)          GENOME="$2"; shift 2 ;;
        --outdir)          OUTPUT_DIR="$2"; shift 2 ;;
        --pam)             PAM="$2"; shift 2 ;;
        --max-mismatches)  MAX_MISMATCHES="$2"; shift 2 ;;
        --device)          DEVICE="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$GRNA_FASTA" ]] && { echo "ERROR: Missing --grna-fasta"; exit 1; }
[[ -z "$GENOME" ]]     && { echo "ERROR: Missing --genome"; exit 1; }

if ! command -v cas-offinder &>/dev/null; then
    echo "WARNING: cas-offinder not found in PATH. Skipping Cas-OFFinder analysis."
    exit 0
fi

BASENAME=$(basename "$GRNA_FASTA" .fasta)
BASENAME=$(basename "$BASENAME" .fa)

CASOFF_DIR="$OUTPUT_DIR/cas_offinder"
mkdir -p "$CASOFF_DIR/output"

# Build Cas-OFFinder input file from gRNA FASTA
# Format:
#   Line 1: path to genome FASTA
#   Line 2: PAM pattern
#   Lines 3+: gRNA_sequence<space>max_mismatches
INPUT_FILE="$CASOFF_DIR/${BASENAME}_input.txt"
{
    echo "$GENOME"
    echo "$PAM"
    # Extract sequences from FASTA (skip headers), append PAM-matching suffix + mismatch count
    while IFS= read -r line; do
        [[ "$line" =~ ^'>' ]] && continue
        [[ -z "$line" ]] && continue
        seq=$(echo "$line" | tr -d '[:space:]')
        echo "${seq} ${MAX_MISMATCHES}"
    done < "$GRNA_FASTA"
} > "$INPUT_FILE"

OUTPUT_FILE="$CASOFF_DIR/output/${BASENAME}_output.txt"

# Run Cas-OFFinder with GPU; fall back to CPU on OpenCL build failure (common in WSL2)
run_cas_offinder() {
    local device="$1"
    echo "Running Cas-OFFinder: $BASENAME (device=$device, mismatches<=$MAX_MISMATCHES)"
    cas-offinder "$INPUT_FILE" "$device" "$OUTPUT_FILE" 2>&1
}

CASOFF_LOG=$(run_cas_offinder "$DEVICE") || true
echo "$CASOFF_LOG"

if echo "$CASOFF_LOG" | grep -q "clBuildProgram Failed" && [[ "$DEVICE" != "C" ]]; then
    echo "WARNING: OpenCL GPU build failed for $BASENAME — retrying with CPU (C)"
    rm -f "$OUTPUT_FILE"
    CASOFF_LOG=$(run_cas_offinder "C") || true
    echo "$CASOFF_LOG"
fi

if [[ -f "$OUTPUT_FILE" ]]; then
    TOTAL=$(wc -l < "$OUTPUT_FILE")
    echo "Cas-OFFinder hits: $TOTAL -> $OUTPUT_FILE"
else
    echo "WARNING: No output produced for $BASENAME"
fi
