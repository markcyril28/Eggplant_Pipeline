#!/bin/bash
# ============================================================================
# Module 12: MFEprimer-3.0 wrapper for one primer set against one genome
# ----------------------------------------------------------------------------
# Reference: Wang et al. 2019, Nucleic Acids Research, doi:10.1093/nar/gkz351
# Usage (called by 12_In_Silico_PCR.sh):
#   bash mfeprimer_run.sh --primers-tsv <tsv> --set-name <name> \
#       --index-dir <dir> --genome-name <name> --outdir <dir> [tunables]
#
# Input TSV columns: primer_id, direction, sequence, expected_size (optional)
#   direction = F or R; one row per primer; pairs combined internally into FASTA.
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

# Locate the indexed FASTA in INDEX_DIR. MFEprimer v4.x writes <fasta>.primerqc
# (the k-mer database); v3.x writes <fasta>.uni; pre-v3 writes <fasta>.ufm.
# Check all three. The FASTA itself is the value passed to -d.
idx_target=$(find "$INDEX_DIR" -maxdepth 1 \
                  \( -name '*.primerqc' -o -name '*.uni' -o -name '*.ufm' \) \
                  ! -name '*.fai' | head -1 || true)
if [[ -z "$idx_target" ]]; then
    echo "[mfeprimer_run] no MFEprimer index (.primerqc/.uni/.ufm) in $INDEX_DIR — run index_genomes first" >&2
    exit 1
fi
# Strip whichever extension the index file uses to recover the FASTA path
indexed_fa="${idx_target%.primerqc}"
indexed_fa="${indexed_fa%.uni}"
indexed_fa="${indexed_fa%.ufm}"

out_prefix="$OUTDIR/${SET_NAME}__${GENOME_NAME}"
out_json="${out_prefix}.json"
out_tsv="${out_prefix}.tsv"

if [[ -s "$out_json" && "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
    echo "[mfeprimer_run] $out_json exists — skip"
    exit 0
fi
# mfeprimer spec refuses to run if any of its output files already exist.
# Different mfeprimer versions emit different suffixes (text, .json, .html,
# .spec.tsv, etc.), so wipe everything matching the prefix to stay general.
# Preserve the .tsv summary we wrote ourselves only on skip; on overwrite,
# the run will regenerate it.
if [[ "$OVERWRITE" == "true" || "$OVERWRITE" == "True" ]]; then
    # shellcheck disable=SC2115
    rm -f "$out_prefix" "$out_prefix".*
fi

# ── Convert TSV → FASTA primer pairs (MFEprimer expects <id>_F / <id>_R) ──
primers_fa=$(mktemp "${TMPDIR:-/tmp}/${SET_NAME}_${GENOME_NAME}_XXXXXX.fa")
trap 'rm -f "$primers_fa"' EXIT

python3 - "$PRIMERS_TSV" "$primers_fa" <<'PY'
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
            fout.write(f">{pid}_F\n{dirs['F']}\n>{pid}_R\n{dirs['R']}\n")
        else:
            missing = 'R' if 'F' in dirs else 'F'
            print(f"  [WARN] {pid}: missing {missing} primer — skipped", file=sys.stderr)
if not any('F' in d and 'R' in d for d in pairs.values()):
    sys.exit("ERROR: no complete primer pairs found in TSV")
PY

# If the index is a symlink to a /tmp staged file, follow it back so we run
# mfeprimer in the same dir the index actually lives in. Same cross-device
# rename hazard as build_indices: mfeprimer spec also uses tempfile rename.
real_indexed_fa=$(readlink -f "$indexed_fa")
spec_run_dir=$(dirname "$real_indexed_fa")
spec_relative_db=$(basename "$real_indexed_fa")

echo "[mfeprimer_run] cd $spec_run_dir && $MFE_BIN spec  primers=$(grep -c '^>' "$primers_fa") set=$SET_NAME genome=$GENOME_NAME GOMAXPROCS=$THREADS -c $THREADS"
# GOMAXPROCS caps Go OS threads; -c bounds mfeprimer's worker pool; -m caps
# RAM usage to leave headroom for the orchestrator on WSL2.
( cd "$spec_run_dir" \
  && GOMAXPROCS="$THREADS" "$MFE_BIN" spec \
      -i "$primers_fa" \
      -d "$spec_relative_db" \
      -k "$INDEX_K" \
      -S "$MAX_AMP" -s "$MIN_AMP" \
      --misEnd "$MIS_END" \
      --diva "$DIVALENT" --mono "$MONOVALENT" \
      --dntp "$DNTP" --oligo "$OLIGO" \
      -c "$THREADS" \
      -j -o "$out_prefix" )

# mfeprimer writes <out_prefix> (text), <out_prefix>.json, and <out_prefix>.html
# Promote a flat TSV summary for downstream merge.
# Schema notes (MFEprimer 3.x AmpList entry):
#   amp.F.Seq.ID / amp.R.Seq.ID  - primer IDs (named "<pair>_F" / "<pair>_R" by us)
#   amp.Hid / amp.Chr            - chromosome / hit contig
#   amp.Size                     - product size (bp)
#   amp.F.Start / amp.R.End      - amplicon span on the contig
#   amp.GC                       - product GC %
# A chimeric amplicon (F from pair X paired with R from pair Y) is itself a
# primer cross-talk off-target, so we keep both primer IDs and a derived
# pair_name = X if both match, else "X__Y" (sorted) for cross-pair hits.
if [[ -f "$out_json" ]]; then
    out_pair_summary="${out_prefix}.bands_per_pair.tsv"
    python3 - "$out_json" "$out_tsv" "$out_pair_summary" "$GENOME_NAME" <<'PY'
import json, sys, csv
from collections import defaultdict
src, dst, pair_dst, genome = sys.argv[1:5]
with open(src) as fh:
    data = json.load(fh)
amps = (data if isinstance(data, list)
        else data.get('AmpList') or data.get('amplicons') or data.get('Amps') or [])

def stem(pid):
    return pid.rsplit('_', 1)[0] if pid.endswith(('_F', '_R')) else pid

def pair_label(fs, rs):
    if not fs and not rs:
        return ''
    if fs and rs and fs == rs:
        return fs
    if fs and rs:
        return '__'.join(sorted([fs, rs]))
    return fs or rs

fields = ['genome', 'primer_id', 'f_primer_id', 'r_primer_id',
          'is_chimeric', 'chrom', 'start', 'end', 'size',
          'tm_f', 'tm_r', 'dg_f', 'dg_r', 'product_gc', 'ppc']

# Track band counts per pair label (clean + chimeric separately).
band_count = defaultdict(lambda: {'clean': 0, 'chimeric': 0, 'chroms': set()})

with open(dst, 'w', newline='') as fout:
    w = csv.writer(fout, delimiter='\t')
    w.writerow(fields)
    for a in amps:
        f = a.get('F', {})
        r = a.get('R', {})
        f_pid = (f.get('Seq') or {}).get('ID', '') or f.get('ID', '')
        r_pid = (r.get('Seq') or {}).get('ID', '') or r.get('ID', '')
        f_stem = stem(f_pid)
        r_stem = stem(r_pid)
        chimeric = bool(f_stem and r_stem and f_stem != r_stem)
        pair = pair_label(f_stem, r_stem)
        chrom = a.get('Hid') or a.get('Chr') or ''
        f_start = f.get('Start', '')
        r_end = r.get('End', '')
        size = a.get('Size', '') or a.get('ProdSize', '')
        w.writerow([
            genome, pair, f_pid, r_pid,
            'yes' if chimeric else 'no',
            chrom, f_start, r_end, size,
            f.get('Tm', ''), r.get('Tm', ''),
            f.get('Dg', ''), r.get('Dg', ''),
            a.get('GC', ''), a.get('PPC', ''),
        ])
        if pair:
            key = 'chimeric' if chimeric else 'clean'
            band_count[pair][key] += 1
            if chrom:
                band_count[pair]['chroms'].add(str(chrom))

# Per-pair band-count rollup. Lists every input pair from the FASTA so
# zero-band pairs are visible too. Pull pair names from the PrimerList block.
all_pairs = set()
for p in (data.get('PrimerList') or []):
    pid = (p.get('Seq') or {}).get('ID', '') or p.get('ID', '')
    if pid:
        all_pairs.add(stem(pid))
all_pairs.update(band_count.keys())
# Drop chimeric "X__Y" composites from the per-pair view so the row count
# matches the number of input primer pairs; cross-talk is reflected in
# each pair's chimeric_bands tally below.
single_pairs = sorted(p for p in all_pairs if '__' not in p)

with open(pair_dst, 'w', newline='') as pf:
    pw = csv.writer(pf, delimiter='\t')
    pw.writerow(['genome', 'pair_name', 'clean_bands', 'chimeric_bands',
                 'total_bands', 'chrom_count', 'chroms'])
    for pair in single_pairs:
        clean = band_count[pair]['clean']
        # A chimeric amplicon involves two pairs; credit both pairs.
        chim = sum(v['chimeric'] for k, v in band_count.items()
                   if '__' in k and pair in k.split('__'))
        chroms = set(band_count[pair]['chroms'])
        for k, v in band_count.items():
            if '__' in k and pair in k.split('__'):
                chroms.update(v['chroms'])
        pw.writerow([genome, pair, clean, chim, clean + chim,
                     len(chroms), ';'.join(sorted(chroms))])
PY
    echo "[mfeprimer_run] wrote $out_tsv"
    echo "[mfeprimer_run] wrote $out_pair_summary"
fi
