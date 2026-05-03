#!/bin/bash
# ============================================================================
# Module 12: MFEprimer-3.0 wrapper for one primer set against one genome
# ----------------------------------------------------------------------------
# Reference: Wang et al. 2019, Nucleic Acids Research, doi:10.1093/nar/gkz351
# Usage (called by 12_In_Silico_PCR.sh):
#   bash mfeprimer_run.sh --primers-tsv <tsv> --set-name <name> \
#       --index-dir <dir> --genome-name <name> --outdir <dir> [tunables]
#
# Input TSV columns: primer_id, forward, reverse, expected_size (optional)
# Converted internally to FASTA pairs (<id>_F / <id>_R) per MFEprimer convention.
# ============================================================================
set -euo pipefail

PRIMERS_TSV="" SET_NAME="" INDEX_DIR="" GENOME_NAME="" OUTDIR=""
MIN_AMP="75" MAX_AMP="1000" MIS_END="3" INDEX_K="9"
DIVALENT="1.5" MONOVALENT="50" DNTP="0.25" OLIGO="50"
THREADS="4" OVERWRITE="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --primers-tsv)  PRIMERS_TSV="$2"; shift 2 ;;
        --set-name)     SET_NAME="$2"; shift 2 ;;
        --index-dir)    INDEX_DIR="$2"; shift 2 ;;
        --genome-name)  GENOME_NAME="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --min-amplicon) MIN_AMP="$2"; shift 2 ;;
        --max-amplicon) MAX_AMP="$2"; shift 2 ;;
        --mis-end)      MIS_END="$2"; shift 2 ;;
        --divalent)     DIVALENT="$2"; shift 2 ;;
        --monovalent)   MONOVALENT="$2"; shift 2 ;;
        --dntp)         DNTP="$2"; shift 2 ;;
        --oligo)        OLIGO="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --index-k)      INDEX_K="$2"; shift 2 ;;
        --overwrite)    OVERWRITE="$2"; shift 2 ;;
        *) echo "[mfeprimer_run] unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -f "$PRIMERS_TSV" ]] || { echo "[mfeprimer_run] primers TSV missing: $PRIMERS_TSV" >&2; exit 1; }
mkdir -p "$OUTDIR"

# Resolve mfeprimer binary: prefer modules/12_in_silico_pcr/bin/ (download_mfeprimer.sh
# installs here); fall back to PATH for users with a manual install.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MFE_BIN="$SCRIPT_DIR/bin/mfeprimer"
if [[ ! -x "$MFE_BIN" ]]; then
    if command -v mfeprimer &>/dev/null; then
        MFE_BIN=$(command -v mfeprimer)
    else
        echo "[mfeprimer_run] WARN: mfeprimer not found." >&2
        echo "  Install with: bash $SCRIPT_DIR/download_mfeprimer.sh" >&2
        echo "  Skipping $SET_NAME / $GENOME_NAME" >&2
        exit 0
    fi
fi

# Locate the .ufm-indexed FASTA inside INDEX_DIR (built by build_indices.sh)
ufm_target=$(find "$INDEX_DIR" -maxdepth 1 -name '*.ufm' | head -1 || true)
if [[ -z "$ufm_target" ]]; then
    echo "[mfeprimer_run] no .ufm index in $INDEX_DIR — run index_genomes first" >&2
    exit 1
fi
indexed_fa="${ufm_target%.ufm}"

out_prefix="$OUTDIR/${SET_NAME}__${GENOME_NAME}"
out_json="${out_prefix}.json"
out_tsv="${out_prefix}.tsv"

if [[ -s "$out_json" && "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
    echo "[mfeprimer_run] $out_json exists — skip"
    exit 0
fi

# ── Convert TSV → FASTA primer pairs (MFEprimer expects <id>_F / <id>_R) ──
primers_fa=$(mktemp "${TMPDIR:-/tmp}/${SET_NAME}_${GENOME_NAME}_XXXXXX.fa")
trap 'rm -f "$primers_fa"' EXIT

python3 - "$PRIMERS_TSV" "$primers_fa" <<'PY'
import sys, csv
src, dst = sys.argv[1], sys.argv[2]
with open(src, newline='') as fin, open(dst, 'w') as fout:
    reader = csv.reader(fin, delimiter='\t')
    header = next(reader)
    cols = {c.lower(): i for i, c in enumerate(header)}
    for need in ('primer_id', 'forward', 'reverse'):
        if need not in cols:
            sys.exit(f"primers TSV missing required column: {need}")
    for row in reader:
        if not row or row[0].startswith('#'):
            continue
        pid = row[cols['primer_id']].strip()
        fwd = row[cols['forward']].strip().upper()
        rev = row[cols['reverse']].strip().upper()
        if not (pid and fwd and rev):
            continue
        fout.write(f">{pid}_F\n{fwd}\n>{pid}_R\n{rev}\n")
PY

echo "[mfeprimer_run] $MFE_BIN spec  primers=$(grep -c '^>' "$primers_fa") set=$SET_NAME genome=$GENOME_NAME"
"$MFE_BIN" spec \
    -i "$primers_fa" \
    -d "$indexed_fa" \
    -k "$INDEX_K" \
    -S "$MAX_AMP" -s "$MIN_AMP" \
    --misEnd "$MIS_END" \
    --divalent "$DIVALENT" --monovalent "$MONOVALENT" \
    --dntp "$DNTP" --oligo "$OLIGO" \
    -c "$THREADS" \
    -j -o "$out_prefix"

# mfeprimer writes <out_prefix> (text), <out_prefix>.json, and <out_prefix>.html
# Promote a flat TSV summary for downstream merge.
if [[ -f "$out_json" ]]; then
    python3 - "$out_json" "$out_tsv" "$GENOME_NAME" <<'PY'
import json, sys, csv
src, dst, genome = sys.argv[1], sys.argv[2], sys.argv[3]
with open(src) as fh:
    data = json.load(fh)
amps = data.get('AmpList') or data.get('amplicons') or []
fields = ['genome','primer_id','chrom','start','end','size','tm_f','tm_r','dg_f','dg_r','product_gc']
with open(dst, 'w', newline='') as fout:
    w = csv.writer(fout, delimiter='\t')
    w.writerow(fields)
    for a in amps:
        f = a.get('F', {}); r = a.get('R', {})
        # PrimerID is stored like "<id>_F" — strip the suffix to recover the pair name
        pid_raw = f.get('ID', '') or r.get('ID', '')
        pid = pid_raw.rsplit('_', 1)[0] if pid_raw.endswith(('_F','_R')) else pid_raw
        w.writerow([
            genome, pid,
            a.get('Hit', {}).get('ID', ''),
            a.get('Start', ''), a.get('End', ''), a.get('ProdSize', ''),
            f.get('Tm', ''), r.get('Tm', ''),
            f.get('Dg', ''), r.get('Dg', ''),
            a.get('GC', ''),
        ])
PY
    echo "[mfeprimer_run] wrote $out_tsv"
fi
