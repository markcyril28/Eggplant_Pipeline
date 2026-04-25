#!/bin/bash
# ============================================================================
# Module: 00_prepare_crispor_genome.sh
# Prepare a CRISPOR-compatible genome directory from a plain FASTA.
# ============================================================================
# CRISPOR needs each target genome under tools/crispor/genomes/{NAME}/ with:
#   {NAME}.fa                       (renamed copy of the input FASTA)
#   {NAME}.fa.{amb,ann,bwt,pac,sa}  (BWA index)
#   {NAME}.2bit                     (UCSC 2bit for efficiency scoring)
#   {NAME}.sizes                    (chromosome sizes)
#   genomeInfo.tab                  (minimal metadata row)
#
# This script is idempotent — re-running only rebuilds missing files.
#
# Usage:
#   bash 00_prepare_crispor_genome.sh \
#       --genome-fasta <path> \
#       --name         melongena \
#       [--crispor-dir <path>]   # default: sibling tools/crispor
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CRISPOR_DIR="$SCRIPT_DIR/tools/crispor"

# Source the shared logger so stdout lines get [YYYY-MM-DD HH:MM:SS] [INFO]
# prefixes matching the orchestrator log format. Non-fatal if absent.
source "$SCRIPT_DIR/../../utils/logging.sh" 2>/dev/null || true
if ! declare -F log_info >/dev/null 2>&1; then
    log_info() { printf '[%(%Y-%m-%d %H:%M:%S)T] [INFO] %s\n' -1 "$*"; }
fi

GENOME_FASTA=""
NAME=""
CRISPOR_DIR="$DEFAULT_CRISPOR_DIR"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome-fasta) GENOME_FASTA="$2"; shift 2 ;;
        --name)         NAME="$2";          shift 2 ;;
        --crispor-dir)  CRISPOR_DIR="$2";   shift 2 ;;
        *) echo "Unknown flag: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$GENOME_FASTA" ]] && { echo "ERROR: --genome-fasta required" >&2; exit 1; }
[[ -z "$NAME"         ]] && { echo "ERROR: --name required"         >&2; exit 1; }
[[ -f "$GENOME_FASTA" ]] || { echo "ERROR: genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }

BIN_DIR="$CRISPOR_DIR/bin/Linux-x86_64"
BWA="$BIN_DIR/bwa"
FATOTWOBIT="$BIN_DIR/faToTwoBit"
TWOBITINFO="$BIN_DIR/twoBitInfo"

for tool in "$BWA" "$FATOTWOBIT" "$TWOBITINFO"; do
    [[ -x "$tool" ]] || { echo "ERROR: required tool not executable: $tool" >&2; exit 1; }
done

GENOME_DIR="$CRISPOR_DIR/genomes/$NAME"
TARGET_FA="$GENOME_DIR/${NAME}.fa"
TARGET_2BIT="$GENOME_DIR/${NAME}.2bit"
TARGET_SIZES="$GENOME_DIR/${NAME}.sizes"
INFO_TAB="$GENOME_DIR/genomeInfo.tab"

mkdir -p "$GENOME_DIR"

# ─── 1. Copy/link FASTA into the genome dir ──────────────────────────────────
if [[ ! -f "$TARGET_FA" ]]; then
    log_info "[prepare_crispor_genome] Linking FASTA: $GENOME_FASTA -> $TARGET_FA"
    ln -sf "$(readlink -f "$GENOME_FASTA")" "$TARGET_FA"
fi

# ─── 2. BWA index ────────────────────────────────────────────────────────────
need_bwa=0
for ext in amb ann bwt pac sa; do
    [[ -f "$TARGET_FA.$ext" ]] || need_bwa=1
done
if (( need_bwa )); then
    log_info "[prepare_crispor_genome] Building BWA index (this takes a few minutes)..."
    "$BWA" index "$TARGET_FA"
else
    log_info "[prepare_crispor_genome] BWA index present — skipping."
fi

# ─── 3. 2bit + sizes ─────────────────────────────────────────────────────────
if [[ ! -f "$TARGET_2BIT" ]]; then
    log_info "[prepare_crispor_genome] Building 2bit: $TARGET_2BIT"
    "$FATOTWOBIT" "$TARGET_FA" "$TARGET_2BIT"
fi
if [[ ! -f "$TARGET_SIZES" ]]; then
    log_info "[prepare_crispor_genome] Writing chrom sizes: $TARGET_SIZES"
    "$TWOBITINFO" "$TARGET_2BIT" "$TARGET_SIZES"
fi

# ─── 4. Minimal genomeInfo.tab ───────────────────────────────────────────────
if [[ ! -f "$INFO_TAB" ]]; then
    log_info "[prepare_crispor_genome] Writing genomeInfo.tab"
    {
        printf '#name\tdescription\tnibPath\torganism\tdefaultPos\tactive\torderKey\tgenome\tscientificName\thtmlPath\thgNearOk\thgPbOk\tsourceName\ttaxId\tserver\n'
        printf '%s\tCustom genome %s\t/gbdb/%s\tCustom\t\t1\t1\t%s\tCustom\t\t0\t0\tLocal\t0\tlocal\n' \
            "$NAME" "$NAME" "$NAME" "$NAME"
    } > "$INFO_TAB"
fi

log_info "[prepare_crispor_genome] Genome ready: $GENOME_DIR"
