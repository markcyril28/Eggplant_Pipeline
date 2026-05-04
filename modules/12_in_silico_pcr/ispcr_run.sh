#!/bin/bash
# ============================================================================
# Module 12: UCSC isPcr wrapper for one primer set against one genome
# ----------------------------------------------------------------------------
# Reference: Jim Kent, UCSC Genome Browser - kent source utilities
# Distribution: https://hgdownload.soe.ucsc.edu/admin/exe/
#
# isPcr usage:
#   isPcr [options] database query.tsv output.fa
# query.tsv columns (no header):  name<TAB>fwd_primer<TAB>rev_primer
#
# Options used:
#   -minSize=N    minimum amplicon size
#   -maxSize=N    maximum amplicon size
#   -minPerfect=N minimum perfect match at 3' end
#   -minGood=N    minimum 'good' match (allowing some MM)
#   -tileSize=N   index tile size (must match .ooc if used)
#   -ooc=FILE     pre-built repeat tile mask (built by build_indices.sh)
#   -flipReverse  reverse-complement reverse primer in output
#   -out=fa       FASTA output (default)
# ============================================================================
set -euo pipefail

PRIMERS_TSV="" SET_NAME="" GENOME="" GENOME_NAME="" INDEX_DIR="" OUTDIR=""
ISPCR_BIN="" MIN_SIZE="75" MAX_SIZE="4000" MIN_PERFECT="15" MIN_GOOD="15"
TILE_SIZE="11" FLIP_REVERSE="true" OVERWRITE="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --primers-tsv)   PRIMERS_TSV="$2"; shift 2 ;;
        --set-name)      SET_NAME="$2"; shift 2 ;;
        --genome)        GENOME="$2"; shift 2 ;;
        --genome-name)   GENOME_NAME="$2"; shift 2 ;;
        --index-dir)     INDEX_DIR="$2"; shift 2 ;;
        --outdir)        OUTDIR="$2"; shift 2 ;;
        --ispcr-bin)     ISPCR_BIN="$2"; shift 2 ;;
        --min-size)      MIN_SIZE="$2"; shift 2 ;;
        --max-size)      MAX_SIZE="$2"; shift 2 ;;
        --min-perfect)   MIN_PERFECT="$2"; shift 2 ;;
        --min-good)      MIN_GOOD="$2"; shift 2 ;;
        --tile-size)     TILE_SIZE="$2"; shift 2 ;;
        --flip-reverse)  FLIP_REVERSE="$2"; shift 2 ;;
        --overwrite)     OVERWRITE="$2"; shift 2 ;;
        *) echo "[ispcr_run] unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -f "$PRIMERS_TSV" ]] || { echo "[ispcr_run] primers TSV missing: $PRIMERS_TSV" >&2; exit 1; }
[[ -f "$GENOME" ]]      || { echo "[ispcr_run] genome FASTA missing: $GENOME" >&2; exit 1; }
[[ -x "$ISPCR_BIN" ]]   || { echo "[ispcr_run] isPcr binary missing/not exec: $ISPCR_BIN" >&2; exit 1; }
mkdir -p "$OUTDIR"

out_fa="$OUTDIR/${SET_NAME}__${GENOME_NAME}.fa"
if [[ -s "$out_fa" && "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
    echo "[ispcr_run] $out_fa exists - skip"
    exit 0
fi

# ── Convert pipeline TSV (with header) → isPcr 3-col TSV (no header) ──────
ispcr_in=$(mktemp "${TMPDIR:-/tmp}/${SET_NAME}_${GENOME_NAME}_XXXXXX.tsv")
trap 'rm -f "$ispcr_in"' EXIT

python3 - "$PRIMERS_TSV" "$ispcr_in" <<'PY'
import sys, csv
from collections import defaultdict
src, dst = sys.argv[1], sys.argv[2]
pairs = defaultdict(dict)
with open(src, newline='') as fin:
    reader = csv.reader(fin, delimiter='\t')
    header = next((r for r in reader if r and not r[0].startswith('#')), None)
    if header is None:
        sys.exit("ERROR: no header row found in TSV")
    cols = {c.lower(): i for i, c in enumerate(header)}
    for need in ('primer_id', 'direction', 'sequence'):
        if need not in cols:
            sys.exit(f"primers TSV missing required column: {need}")
    for row in reader:
        if not row or row[0].startswith('#'):
            continue
        pid = row[cols['primer_id']].strip()
        direction = row[cols['direction']].strip().upper()
        seq = row[cols['sequence']].strip().upper()
        if pid and direction in ('F', 'R') and seq:
            pairs[pid][direction] = seq
with open(dst, 'w') as fout:
    for pid, dirs in pairs.items():
        if 'F' in dirs and 'R' in dirs:
            fout.write(f"{pid}\t{dirs['F']}\t{dirs['R']}\n")
        else:
            missing = 'R' if 'F' in dirs else 'F'
            print(f"  [WARN] {pid}: missing {missing} primer - skipped", file=sys.stderr)
if not any('F' in d and 'R' in d for d in pairs.values()):
    sys.exit("ERROR: no complete primer pairs found in TSV")
PY

OPTS=(
    -minSize="$MIN_SIZE"
    -maxSize="$MAX_SIZE"
    -minPerfect="$MIN_PERFECT"
    -minGood="$MIN_GOOD"
    -tileSize="$TILE_SIZE"
    -out=fa
)
[[ "$FLIP_REVERSE" == "true" || "$FLIP_REVERSE" == "True" ]] && OPTS+=(-flipReverse)

ooc_file="$INDEX_DIR/${GENOME_NAME}.ooc"
[[ -f "$ooc_file" ]] && OPTS+=(-ooc="$ooc_file")

echo "[ispcr_run] isPcr ${OPTS[*]}  $GENOME  $ispcr_in  $out_fa"
"$ISPCR_BIN" "${OPTS[@]}" "$GENOME" "$ispcr_in" "$out_fa"

# Emit a flat TSV alongside the FASTA for easy merging.
out_tsv="${out_fa%.fa}.tsv"
python3 - "$out_fa" "$out_tsv" "$GENOME_NAME" <<'PY'
import sys, re, os
src, dst, genome = sys.argv[1], sys.argv[2], sys.argv[3]
hdr_re = re.compile(r'^>(\S+)\s+(\S+):(\d+)([+-])(\d+)\s+(\d+)bp\s+(\S+)\s+(\S+)')
with open(dst, 'w') as fout:
    fout.write("genome\tprimer_id\tchrom\tstart\tstrand\tend\tsize\tforward\treverse\n")
    if not os.path.exists(src):
        print("  [INFO] isPcr produced no amplicons (output FASTA absent)", file=sys.stderr)
    else:
        with open(src) as fh:
            for line in fh:
                if not line.startswith('>'):
                    continue
                m = hdr_re.match(line.strip())
                if not m:
                    continue
                pid, chrom, start, strand, end, size, fwd, rev = m.groups()
                fout.write(f"{genome}\t{pid}\t{chrom}\t{start}\t{strand}\t{end}\t{size}\t{fwd}\t{rev}\n")
PY
echo "[ispcr_run] wrote $out_tsv"
