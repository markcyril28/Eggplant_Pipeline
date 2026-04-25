#!/bin/bash
# ============================================================================
# Module: 01_design_score.sh
# Stage [1] — gRNA design + on/off-target scoring via CRISPOR
# ============================================================================
# Usage (orchestrated):
#   bash 01_design_score.sh \
#       --grna-fasta  <path>   \
#       --genome-fasta <path>  \
#       --outdir      <path>   \
#       --pam         NGG      \
#       --crispor-genome melongena \
#       --min-score   20       \
#       --batch-size  100      \
#       --threads     4        \
#       --overwrite   true
#
# Standalone (hardcoded defaults for development):
#   bash 01_design_score.sh
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../utils/logging.sh" 2>/dev/null || true

# ─── Defaults (standalone / dev mode) ────────────────────────────────────────
GRNA_FASTA=""
GENOME_FASTA=""
OUTDIR="."
PAM="NGG"
CRISPOR_GENOME="melongena"
MIN_SCORE=20
BATCH_SIZE=100
THREADS=4
# Max mismatches for CRISPOR's internal bwa off-target search. Overridden by
# the TOML key crispr_v2.curate_offtargets.max_mismatches via --max-mismatches.
# CRISPOR supports values 0–6; 4 is the upstream default.
MAX_MISMATCHES=4
OVERWRITE=true
# Pipe-separated ordered list of CRISPOR on-target score columns the filter
# will try. The first column with a numeric value per guide is compared
# against MIN_SCORE (guides with NA / NotEnoughFlankSeq / missing scores in
# ALL listed columns pass through automatically). Override via --score-columns
# or the TOML key crispr_v2.design_score.score_columns.
# Default keeps Moreno-Mateos first because the Azimuth/Fusi-based
# "Doench '16-Score" returns 0.0 for every guide under sklearn >= 1.0.
SCORE_COLUMNS_STR="Moreno-Mateos-Score|Doench '16-Score|Doench-RuleSet3-Score|DoenchScore|on_target_score"

# ─── Parse flags ─────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --grna-fasta)      GRNA_FASTA="$2";      shift 2 ;;
        --genome-fasta)    GENOME_FASTA="$2";     shift 2 ;;
        --outdir)          OUTDIR="$2";           shift 2 ;;
        --pam)             PAM="$2";              shift 2 ;;
        --crispor-genome)  CRISPOR_GENOME="$2";   shift 2 ;;
        --min-score)       MIN_SCORE="$2";        shift 2 ;;
        --score-columns)   SCORE_COLUMNS_STR="$2"; shift 2 ;;
        --batch-size)      BATCH_SIZE="$2";       shift 2 ;;
        --threads)         THREADS="$2";          shift 2 ;;
        --max-mismatches)  MAX_MISMATCHES="$2";   shift 2 ;;
        --overwrite)       OVERWRITE="$2";        shift 2 ;;
        *) echo "Unknown flag: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$GRNA_FASTA" ]]   && { echo "ERROR: --grna-fasta required" >&2; exit 1; }
[[ -f "$GRNA_FASTA" ]]   || { echo "ERROR: gRNA FASTA not found: $GRNA_FASTA" >&2; exit 1; }

# WSL2/NTFS: mkdir -p may fail on Windows reparse-point dirs; fall back to [[ -d ]].
for _d in "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/offtargets"; do
    mkdir -p "$_d" 2>/dev/null || [[ -d "$_d" ]] || { echo "ERROR: Cannot create $_d" >&2; exit 1; }
done

# ─── Output paths ────────────────────────────────────────────────────────────
# Strip both .fasta and .fa extensions so the gene stem is clean
# (e.g. SMEL5_01g008730.fa  →  SMEL5_01g008730).
GENE_NAME=$(basename "$GRNA_FASTA")
GENE_NAME="${GENE_NAME%.fasta}"
GENE_NAME="${GENE_NAME%.fa}"
GENE_NAME=$(echo "$GENE_NAME" | tr ' ' '_')
OUTFILE="$OUTDIR/${GENE_NAME}.design.tsv"
FILTERED="$OUTDIR/${GENE_NAME}.filtered.tsv"

# Guard on the filtered output (the final deliverable of this stage).
# If only *.design.tsv exists because the filter step previously crashed,
# we must still re-run; guarding on $OUTFILE would silently skip that repair.
if [[ "$OVERWRITE" != "true" && -f "$FILTERED" ]]; then
    log_info "[01_design_score] Output exists, skipping: $FILTERED" 2>/dev/null || \
        echo "[01_design_score] Skipping (overwrite=false): $FILTERED"
    exit 0
fi

log_info "[01_design_score] Running CRISPOR on $(basename "$GRNA_FASTA") ..." 2>/dev/null || \
    echo "[01_design_score] Running CRISPOR on $(basename "$GRNA_FASTA") ..."

# ─── Locate CRISPOR ──────────────────────────────────────────────────────────
# Search order:
#   1. Bundled copy alongside this module (tools/crispor/crispor.py)
#   2. conda env opt/ directory (installed by setup_conda_crispr_v2.sh)
#   3. sys.prefix opt/ (handles 'conda run' where CONDA_PREFIX is unset)
#   4. PATH (if installed as an executable)
BUNDLED_CRISPOR="$SCRIPT_DIR/tools/crispor/crispor.py"
CONDA_CRISPOR="${CONDA_PREFIX:-}/opt/crispor/crispor.py"
SYS_PREFIX_CRISPOR="$(python3 -c 'import sys; print(sys.prefix)' 2>/dev/null)/opt/crispor/crispor.py"
if [[ -f "$BUNDLED_CRISPOR" ]]; then
    CRISPOR_CMD="$BUNDLED_CRISPOR"
elif [[ -f "$CONDA_CRISPOR" ]]; then
    CRISPOR_CMD="$CONDA_CRISPOR"
elif [[ -f "$SYS_PREFIX_CRISPOR" ]]; then
    CRISPOR_CMD="$SYS_PREFIX_CRISPOR"
elif command -v crispor.py &>/dev/null; then
    CRISPOR_CMD="$(command -v crispor.py)"
else
    echo "ERROR: crispor.py not found. Run setup_conda_crispr_v2.sh first." >&2
    echo "  Searched: $BUNDLED_CRISPOR" >&2
    echo "  Searched: $CONDA_CRISPOR" >&2
    echo "  Searched: $SYS_PREFIX_CRISPOR" >&2
    exit 1
fi

# ─── Run CRISPOR in batch mode ───────────────────────────────────────────────
# CRISPOR accepts a FASTA of target sequences (one per guide + 30 bp flanks).
# Output: tab-separated file with Doench2016 and other scores.
# NOTE: CRISPOR runs in CGI mode by default and prints HTML error pages to
# stdout with exit 0 on dependency/runtime failures. We detect this by
# checking (a) that the output TSV was created and (b) that the log does not
# start with a 'Content-type:' HTML response.
CRISPOR_LOG="$OUTDIR/logs/${GENE_NAME}.crispor.log"
OFFTARGET_TSV="$OUTDIR/offtargets/offtargets.tsv"
# CLI reference (bundled crispor.py):
#   positional:  org inFile guideOutFile
#   flags:       -p/--pam, -o/--offtargets, --mm (max mismatches),
#                --skipAlign (input not aligned to genome)
# Note: older scripts used --batchSize / --offtargetMaxMismatches — those
# flags exist only in the web-server fork of CRISPOR and are unsupported here.
set +e
# BWA thread count for this CRISPOR invocation. CRISPOR 's internal `bwa aln`
# call was patched to honour CRISPOR_BWA_THREADS (see crispor.py). The pipeline
# orchestrator divides total cores by MAX_PARALLEL before calling this module,
# so $THREADS here is the per-job allowance — no further division needed.
export CRISPOR_BWA_THREADS="$THREADS"
python3 "$CRISPOR_CMD" \
    "$CRISPOR_GENOME" \
    "$GRNA_FASTA" \
    "$OUTFILE" \
    --pam "$PAM" \
    --offtargets "$OFFTARGET_TSV" \
    --mm "$MAX_MISMATCHES" \
    > "$CRISPOR_LOG" 2>&1
CRISPOR_RC=$?
set -e

if [[ $CRISPOR_RC -ne 0 ]] || [[ ! -s "$OUTFILE" ]] || head -1 "$CRISPOR_LOG" | grep -qi "^Content-type:"; then
    echo "ERROR: CRISPOR failed for $(basename "$GRNA_FASTA")" >&2
    echo "  Exit code: $CRISPOR_RC" >&2
    echo "  Log:       $CRISPOR_LOG" >&2
    echo "  --- Last 20 lines of log ---" >&2
    tail -20 "$CRISPOR_LOG" >&2
    echo "  -----------------------------" >&2
    exit 1
fi

log_info "[01_design_score] Raw CRISPOR output: $OUTFILE" 2>/dev/null || true

# ─── Filter by minimum on-target score ───────────────────────────────────────
# Heredoc emits only the integer guide count on stdout so the orchestrator log
# stays timestamped; the human-readable line is emitted via log_info below.
# Export into the heredoc's env so column names that contain quotes / spaces
# survive without heredoc-expansion escaping issues.
export _FILTER_INFILE="$OUTFILE"
export _FILTER_OUTFILE="$FILTERED"
export _FILTER_MIN_SCORE="$MIN_SCORE"
export _FILTER_SCORE_COLUMNS="$SCORE_COLUMNS_STR"
KEPT_COUNT=$(python3 - <<'PYEOF'
import csv
import os

infile   = os.environ["_FILTER_INFILE"]
outfile  = os.environ["_FILTER_OUTFILE"]
min_sc   = float(os.environ["_FILTER_MIN_SCORE"])
# Pipe-separated ordered column list from the TOML key
# crispr_v2.design_score.score_columns — first numeric column wins.
score_cols = [c for c in os.environ["_FILTER_SCORE_COLUMNS"].split("|") if c]

kept = 0
with open(infile, newline="") as fh, open(outfile, "w", newline="") as out:
    reader = csv.DictReader(fh, delimiter="\t")
    writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t")
    writer.writeheader()
    for row in reader:
        sc_val = None
        for col in score_cols:
            if col in row and row[col] not in ("", "NA", "N/A", "NotEnoughFlankSeq"):
                try:
                    sc_val = float(row[col])
                    break
                except ValueError:
                    pass
        # Guides with no numeric value in any configured column pass through;
        # numeric guides must meet the threshold.
        if sc_val is None or sc_val >= min_sc:
            writer.writerow(row)
            kept += 1

print(kept)
PYEOF
)
unset _FILTER_INFILE _FILTER_OUTFILE _FILTER_MIN_SCORE _FILTER_SCORE_COLUMNS

log_info "[01_design_score] Filtered guides: ${KEPT_COUNT} pass min_score=${MIN_SCORE} -> $FILTERED" 2>/dev/null || \
    echo "[01_design_score] Filtered guides: ${KEPT_COUNT} pass min_score=${MIN_SCORE} -> $FILTERED"
log_info "[01_design_score] Done: $(basename "$GRNA_FASTA")" 2>/dev/null || \
    echo "[01_design_score] Done: $(basename "$GRNA_FASTA")"
