#!/bin/bash
# Module: BLAST Database Creation
# Usage: bash makeblastdb.sh --input <fasta> --dbtype <nucl|prot> --outdir <dir>

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

DB_TYPE="nucl"
OUTPUT_DIR="."
INPUT_FASTA=""
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)   INPUT_FASTA="$2"; shift 2 ;;
        --dbtype)  DB_TYPE="$2"; shift 2 ;;
        --outdir)  OUTPUT_DIR="$2"; shift 2 ;;
        --force)   FORCE=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_FASTA" ]] && { log_error "Missing --input"; exit 1; }

FASTA_BASE=$(basename "$INPUT_FASTA")
FASTA_BASE="${FASTA_BASE%.*}"
DB_DIR="$OUTPUT_DIR/${FASTA_BASE}_db"
DB_PATH="$DB_DIR/$FASTA_BASE"

mkdir -p "$DB_DIR"

# Check existing database
EXT=$([[ "$DB_TYPE" == "nucl" ]] && echo "nsq" || echo "psq")
if [[ -f "${DB_PATH}.${EXT}" && "$FORCE" != true ]]; then
    log_info "Database exists: $DB_PATH (use --force to recreate)"
    echo "$DB_PATH"
    exit 0
fi

log_info "Creating $DB_TYPE database: $FASTA_BASE"
makeblastdb -in "$INPUT_FASTA" -dbtype "$DB_TYPE" -out "$DB_PATH"
log_info "Database created: $DB_PATH"
echo "$DB_PATH"
