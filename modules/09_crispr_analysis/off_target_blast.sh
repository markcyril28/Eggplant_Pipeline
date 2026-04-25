#!/bin/bash
# Module: CRISPR Off-Target Analysis (BLAST-based)
# Usage: bash off_target_blast.sh --grna-fasta <fasta> --genome <ref.fa> --outdir <dir> \
#            [--threads N] [--max-mismatches N]
#
# BLAST parameters are tuned for SpCas9 gRNA off-target detection:
#   - ungapped alignment (CRISPR off-targets are mismatches, not indels)
#   - no dust masking (short guides shouldn't be masked)
#   - high e-value (short queries need permissive threshold)
#   - query coverage filter (reject short partial matches)
#   - high max_target_seqs (capture all off-target loci)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../logging/logging_utils.sh" 2>/dev/null || true

# ── Defaults (most accurate for gRNA off-target detection) ─────
NUM_THREADS=4
MAX_MISMATCHES=4
GENOME=""
GRNA_FASTA=""
OUTPUT_DIR="."

# BLAST accuracy parameters
WORD_SIZE=7
EVALUE=10000
DUST="no"
UNGAPPED=true
REWARD=1
PENALTY=-1
MAX_TARGET_SEQS=10000
QCOV_HSP_PERC=80

# ── Arg parsing ────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --grna-fasta)       GRNA_FASTA="$2";      shift 2 ;;
        --genome)           GENOME="$2";           shift 2 ;;
        --outdir)           OUTPUT_DIR="$2";       shift 2 ;;
        --threads)          NUM_THREADS="$2";      shift 2 ;;
        --max-mismatches)   MAX_MISMATCHES="$2";   shift 2 ;;
        --word-size)        WORD_SIZE="$2";        shift 2 ;;
        --evalue)           EVALUE="$2";           shift 2 ;;
        --dust)             DUST="$2";             shift 2 ;;
        --ungapped)         UNGAPPED="$2";         shift 2 ;;
        --reward)           REWARD="$2";           shift 2 ;;
        --penalty)          PENALTY="$2";          shift 2 ;;
        --max-target-seqs)  MAX_TARGET_SEQS="$2";  shift 2 ;;
        --qcov-hsp-perc)    QCOV_HSP_PERC="$2";   shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$GRNA_FASTA" ]] && { log_error "Missing --grna-fasta"; exit 1; }
[[ -z "$GENOME" ]]     && { log_error "Missing --genome"; exit 1; }

# Build BLAST DB on native Linux filesystem to avoid mmap errors on WSL2 DrvFS (/mnt/*)
GENOME_HASH=$(md5sum "$GENOME" 2>/dev/null | cut -d' ' -f1 || echo "fallback")
DB_TMPDIR="${TMPDIR:-/tmp}/blast_db_${GENOME_HASH}"
mkdir -p "$DB_TMPDIR"
DB_NAME="$DB_TMPDIR/genome_for_gRNAs"

# Create or reuse BLAST database (flock prevents parallel races)
LOCK_FILE="$DB_TMPDIR/.build.lock"
(
    flock -x 200
    if [[ ! -f "${DB_NAME}.nsq" ]]; then
        log_info "Creating BLAST database in $DB_TMPDIR ..."
        makeblastdb -in "$GENOME" -dbtype nucl -out "$DB_NAME"
    else
        log_info "Reusing existing BLAST database: $DB_NAME"
    fi
) 200>"$LOCK_FILE"

BASENAME=$(basename "$GRNA_FASTA" .fasta)
BASENAME=$(basename "$BASENAME" .fa)
RESULT_DIR="$OUTPUT_DIR/results/$BASENAME"
mkdir -p "$RESULT_DIR"

log_step "Off-target analysis: $BASENAME"
log_info "  Parameters: word_size=$WORD_SIZE evalue=$EVALUE dust=$DUST ungapped=$UNGAPPED"
log_info "  reward=$REWARD penalty=$PENALTY max_target_seqs=$MAX_TARGET_SEQS qcov=$QCOV_HSP_PERC%"

# ── Build BLAST command ────────────────────────────────────────
BLAST_COMMON=(
    -query "$GRNA_FASTA"
    -db "$DB_NAME"
    -task blastn-short
    -word_size "$WORD_SIZE"
    -evalue "$EVALUE"
    -dust "$DUST"
    -reward "$REWARD"
    -penalty "$PENALTY"
    -max_target_seqs "$MAX_TARGET_SEQS"
    -qcov_hsp_perc "$QCOV_HSP_PERC"
    -num_threads "$NUM_THREADS"
)

# Ungapped alignment (critical: CRISPR off-targets are mismatches, not indels)
if [[ "$UNGAPPED" == "true" ]]; then
    BLAST_COMMON+=( -ungapped )
fi

# Tabular output — all hits
blastn "${BLAST_COMMON[@]}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore" \
    > "$RESULT_DIR/${BASENAME}_all.txt"

# Full pairwise alignments for manual inspection
blastn "${BLAST_COMMON[@]}" \
    -out "$RESULT_DIR/${BASENAME}_alignments.txt"

# ── Post-filter by mismatch count ─────────────────────────────
awk -v mm="$MAX_MISMATCHES" '$5 <= mm' "$RESULT_DIR/${BASENAME}_all.txt" | \
    sort -k4,4nr -k5,5n > "$RESULT_DIR/${BASENAME}_filtered.txt"

# Split filtered results into per-guide files
awk -v dir="$RESULT_DIR" '{print > dir"/"$1"_filtered.txt"}' "$RESULT_DIR/${BASENAME}_filtered.txt"

TOTAL=$(wc -l < "$RESULT_DIR/${BASENAME}_all.txt")
FILTERED=$(wc -l < "$RESULT_DIR/${BASENAME}_filtered.txt")
log_info "Total hits: $TOTAL | Filtered (mismatches <= $MAX_MISMATCHES): $FILTERED"
log_info "Results -> $RESULT_DIR"
