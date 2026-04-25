#!/bin/bash
# Module: Nucleotide to Protein Translation
# Usage: bash translate.sh --input <nucleotide.fa> --output <protein.fa>

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

INPUT_FASTA=""
OUTPUT_FASTA=""
MAX_PARALLEL=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)   INPUT_FASTA="$2"; shift 2 ;;
        --output)  OUTPUT_FASTA="$2"; shift 2 ;;
        --threads) MAX_PARALLEL="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_FASTA" ]]  && { log_error "Missing --input"; exit 1; }
[[ -z "$OUTPUT_FASTA" ]] && { log_error "Missing --output"; exit 1; }

mkdir -p "$(dirname "$OUTPUT_FASTA")"

# Count sequences to decide whether parallel split is worthwhile
SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA" || echo 0)

if (( SEQ_COUNT <= MAX_PARALLEL || MAX_PARALLEL <= 1 )); then
    # Few sequences or single-thread: translate directly
    log_info "Translating: $(basename "$INPUT_FASTA") ($SEQ_COUNT seqs, 1 job)"
    transeq -sequence "$INPUT_FASTA" -outseq "$OUTPUT_FASTA"
else
    # Split FASTA into chunks, translate in parallel, merge
    CHUNK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/transeq_chunks_XXXXXX")
    trap "rm -rf '$CHUNK_DIR'" EXIT

    log_info "Translating: $(basename "$INPUT_FASTA") ($SEQ_COUNT seqs, $MAX_PARALLEL parallel jobs)"

    # Split into MAX_PARALLEL chunks using awk (split on '>' headers)
    awk -v n="$MAX_PARALLEL" -v dir="$CHUNK_DIR" '
        /^>/ { chunk = (chunk % n) + 1 }
        { print >> (dir "/chunk_" chunk ".fa") }
    ' "$INPUT_FASTA"

    # Translate each chunk in parallel
    PIDS=()
    for chunk in "$CHUNK_DIR"/chunk_*.fa; do
        [[ -f "$chunk" ]] || continue
        out="${chunk%.fa}_prot.fa"
        transeq -sequence "$chunk" -outseq "$out" &
        PIDS+=($!)
    done

    FAILED=0
    for pid in "${PIDS[@]}"; do
        wait "$pid" || ((FAILED++))
    done
    (( FAILED > 0 )) && { log_error "$FAILED translation chunk(s) failed"; exit 1; }

    # Merge translated chunks
    cat "$CHUNK_DIR"/chunk_*_prot.fa > "$OUTPUT_FASTA"
fi

log_info "Translation complete."
