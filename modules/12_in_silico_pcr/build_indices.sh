#!/bin/bash
# ============================================================================
# Module 12: Build per-genome indices for MFEprimer (.ufm) and isPcr (.ooc)
# ----------------------------------------------------------------------------
# References:
#   MFEprimer-3.0   https://github.com/quwubin/MFEprimer-3.0
#   isPcr (UCSC)    https://hgdownload.soe.ucsc.edu/admin/exe/
#
# Both indices are built once per genome and reused across primer sets.
# ============================================================================
set -euo pipefail

GENOME="" GENOME_NAME="" OUTDIR="" KMER="9" OOC_TILE="11" OOC_REPEAT="1024"
ISPCR_BIN="" OVERWRITE="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)       GENOME="$2"; shift 2 ;;
        --genome-name)  GENOME_NAME="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --kmer)         KMER="$2"; shift 2 ;;
        --ooc-tile)     OOC_TILE="$2"; shift 2 ;;
        --ooc-repeat)   OOC_REPEAT="$2"; shift 2 ;;
        --ispcr-bin)    ISPCR_BIN="$2"; shift 2 ;;
        --overwrite)    OVERWRITE="$2"; shift 2 ;;
        *) echo "[build_indices] unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -f "$GENOME" ]] || { echo "[build_indices] genome FASTA missing: $GENOME" >&2; exit 1; }
mkdir -p "$OUTDIR"

genome_basename=$(basename "$GENOME")
ufm_link="$OUTDIR/${genome_basename}"
ufm_index="${ufm_link}.ufm"
ooc_file="$OUTDIR/${GENOME_NAME}.ooc"

# ── MFEprimer index (build under OUTDIR via a symlink so .ufm sits in OUTDIR)
if command -v mfeprimer &>/dev/null; then
    if [[ -f "$ufm_index" && "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
        echo "[build_indices] [$GENOME_NAME] MFEprimer index exists — skip"
    else
        ln -sf "$GENOME" "$ufm_link"
        echo "[build_indices] [$GENOME_NAME] mfeprimer index -i $ufm_link -k $KMER"
        mfeprimer index -i "$ufm_link" -k "$KMER" -f
    fi
else
    echo "[build_indices] WARN: mfeprimer not found in PATH — skipping MFEprimer index" >&2
fi

# ── isPcr .ooc (only if binary present) ────────────────────────────────────
if [[ -x "$ISPCR_BIN" ]]; then
    blat_bin="$(dirname "$ISPCR_BIN")/blat"
    if [[ -f "$ooc_file" && "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
        echo "[build_indices] [$GENOME_NAME] isPcr .ooc exists — skip"
    elif [[ -x "$blat_bin" ]]; then
        echo "[build_indices] [$GENOME_NAME] blat .ooc tile=$OOC_TILE repMatch=$OOC_REPEAT"
        "$blat_bin" "$GENOME" /dev/null /dev/null \
            -tileSize="$OOC_TILE" -repMatch="$OOC_REPEAT" -makeOoc="$ooc_file" || true
    else
        echo "[build_indices] [$GENOME_NAME] blat not next to isPcr — .ooc skipped (isPcr will run without it, slower)"
    fi
else
    echo "[build_indices] [$GENOME_NAME] isPcr binary absent — skipping .ooc"
fi
