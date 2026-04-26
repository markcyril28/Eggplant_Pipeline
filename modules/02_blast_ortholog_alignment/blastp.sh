#!/bin/bash
# Module: BLASTp Search
# Usage: bash blastp.sh --query <fasta> --db <db_path> --evalue <val> --word-size <ws> --outdir <dir> [--threads N]
#        bash blastp.sh --config <config.toml> --query <fasta> --db <db_path> --evalue <val> --word-size <ws> --outdir <dir> [--threads N]
#        Additional flags: --max-target-seqs N --max-hsps N --qcov-hsp-perc N --matrix <name> --seg <yes|no>

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

# Default values for standalone mode
NUM_THREADS=4
E_VALUE="1e-15"
WORD_SIZE=2
MAX_TARGET_SEQS=500
MAX_HSPS=1
QCOV_HSP_PERC=0
MATRIX="BLOSUM62"
SEG="yes"
OUTPUT_DIR="."
QUERY_FASTA=""
DB_PATH=""
GENOME_NAME=""
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)           CONFIG_FILE="$2"; shift 2 ;;
        --query)            QUERY_FASTA="$2"; shift 2 ;;
        --db)               DB_PATH="$2"; shift 2 ;;
        --evalue)           E_VALUE="$2"; shift 2 ;;
        --word-size)        WORD_SIZE="$2"; shift 2 ;;
        --max-target-seqs)  MAX_TARGET_SEQS="$2"; shift 2 ;;
        --max-hsps)         MAX_HSPS="$2"; shift 2 ;;
        --qcov-hsp-perc)    QCOV_HSP_PERC="$2"; shift 2 ;;
        --matrix)           MATRIX="$2"; shift 2 ;;
        --seg)              SEG="$2"; shift 2 ;;
        --outdir)           OUTPUT_DIR="$2"; shift 2 ;;
        --threads)          NUM_THREADS="$2"; shift 2 ;;
        --genome-name)      GENOME_NAME="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Load configuration from TOML if provided
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    # Import the TOML parser function
    TOML_PARSER="$SCRIPT_DIR/../utils/parse_toml.py"
    
    # Override defaults with config values if not set via command line
    if [[ -z "$NUM_THREADS" || "$NUM_THREADS" == "4" ]]; then
        NUM_THREADS=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "threads" 2>/dev/null || echo "$NUM_THREADS")
    fi
    
    if [[ -z "$E_VALUE" || "$E_VALUE" == "1e-15" ]]; then
        E_VALUE=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "e_value" 2>/dev/null || echo "$E_VALUE")
    fi
    
    if [[ -z "$WORD_SIZE" || "$WORD_SIZE" == "2" ]]; then
        WORD_SIZE=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "word_size" 2>/dev/null || echo "$WORD_SIZE")
    fi

    if [[ "$MAX_TARGET_SEQS" == "500" ]]; then
        MAX_TARGET_SEQS=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "max_target_seqs" 2>/dev/null || echo "$MAX_TARGET_SEQS")
    fi

    if [[ "$MAX_HSPS" == "1" ]]; then
        MAX_HSPS=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "max_hsps" 2>/dev/null || echo "$MAX_HSPS")
    fi

    if [[ "$QCOV_HSP_PERC" == "0" ]]; then
        QCOV_HSP_PERC=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "qcov_hsp_perc" 2>/dev/null || echo "$QCOV_HSP_PERC")
    fi

    if [[ "$MATRIX" == "BLOSUM62" ]]; then
        MATRIX=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "matrix" 2>/dev/null || echo "$MATRIX")
    fi

    if [[ "$SEG" == "yes" ]]; then
        SEG=$(python3 "$TOML_PARSER" "$CONFIG_FILE" "ortholog_blast" "blastp_params" "seg" 2>/dev/null || echo "$SEG")
    fi
fi

[[ -z "$QUERY_FASTA" ]] && { log_error "Missing --query"; exit 1; }
[[ -z "$DB_PATH" ]]     && { log_error "Missing --db"; exit 1; }

DB_BASENAME=$(basename "$DB_PATH")
QUERY_BASENAME=$(basename "$QUERY_FASTA" .fasta)
QUERY_BASENAME=$(basename "$QUERY_BASENAME" .fa)
TIMESTAMP=$(date +%F_%H-%M)

# Determine genome-specific output directory
if [[ -n "$GENOME_NAME" ]]; then
    RESULT_DIR="$OUTPUT_DIR/$GENOME_NAME/blastp_results/ev_${E_VALUE}_ws_${WORD_SIZE}_${DB_BASENAME}_${TIMESTAMP}"
else
    # Fallback to original structure if no genome name provided
    RESULT_DIR="$OUTPUT_DIR/blastp_results/ev_${E_VALUE}_ws_${WORD_SIZE}_${DB_BASENAME}_${TIMESTAMP}"
fi

mkdir -p "$RESULT_DIR"

ARCHIVE="$RESULT_DIR/blastp_${DB_BASENAME}_VS_${QUERY_BASENAME}_${E_VALUE}_${WORD_SIZE}.asn"
CSV_OUT="$RESULT_DIR/blastp_${DB_BASENAME}_VS_${QUERY_BASENAME}_${E_VALUE}_${WORD_SIZE}.csv"
TXT_OUT="$RESULT_DIR/blastp_${DB_BASENAME}_VS_${QUERY_BASENAME}_${E_VALUE}_${WORD_SIZE}.txt"

log_info "BLASTp: $QUERY_BASENAME vs $DB_BASENAME (e=$E_VALUE, ws=$WORD_SIZE, matrix=$MATRIX)"

# Build optional flags
OPTIONAL_FLAGS=()
if [[ "$MAX_TARGET_SEQS" -gt 0 ]]; then
    OPTIONAL_FLAGS+=(-max_target_seqs "$MAX_TARGET_SEQS")
fi
if [[ "$MAX_HSPS" -gt 0 ]]; then
    OPTIONAL_FLAGS+=(-max_hsps "$MAX_HSPS")
fi
if [[ "$QCOV_HSP_PERC" -gt 0 ]]; then
    OPTIONAL_FLAGS+=(-qcov_hsp_perc "$QCOV_HSP_PERC")
fi

# Single BLAST run: save as archive, then derive CSV and TXT
blastp -query "$QUERY_FASTA" -db "$DB_PATH" \
    -evalue "$E_VALUE" -word_size "$WORD_SIZE" \
    -matrix "$MATRIX" \
    -seg "$SEG" \
    -num_threads "$NUM_THREADS" \
    "${OPTIONAL_FLAGS[@]}" \
    -outfmt 11 -out "$ARCHIVE"

echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Score" > "$CSV_OUT"
blast_formatter -archive "$ARCHIVE" \
    -outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score" >> "$CSV_OUT"

blast_formatter -archive "$ARCHIVE" -outfmt 0 > "$TXT_OUT"

rm -f "$ARCHIVE"

log_info "BLASTp complete -> $RESULT_DIR"
