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
ISPCR_BIN="" OVERWRITE="true" THREADS="1" ENGINE="both"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)       GENOME="$2"; shift 2 ;;
        --genome-name)  GENOME_NAME="$2"; shift 2 ;;
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --kmer)         KMER="$2"; shift 2 ;;
        --ooc-tile)     OOC_TILE="$2"; shift 2 ;;
        --ooc-repeat)   OOC_REPEAT="$2"; shift 2 ;;
        --ispcr-bin)    ISPCR_BIN="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --engine)       ENGINE="$2"; shift 2 ;;     # mfeprimer | ispcr | both
        --overwrite)    OVERWRITE="$2"; shift 2 ;;
        *) echo "[build_indices] unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -f "$GENOME" ]] || { echo "[build_indices] genome FASTA missing: $GENOME" >&2; exit 1; }
mkdir -p "$OUTDIR"

genome_basename=$(basename "$GENOME")
ufm_link="$OUTDIR/${genome_basename}"
# MFEprimer v3.x writes <fasta>.ufm; v4.x writes <fasta>.uni alongside .fai/.2bit.
# Treat either extension as proof the index exists.
ooc_file="$OUTDIR/${GENOME_NAME}.ooc"
mfe_index_present() {
    [[ -f "${ufm_link}.uni" || -f "${ufm_link}.ufm" ]]
}

# ── MFEprimer index (build under OUTDIR via a symlink so .ufm sits in OUTDIR)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MFE_BIN="$SCRIPT_DIR/bin/mfeprimer"
[[ -x "$MFE_BIN" ]] || MFE_BIN=$(command -v mfeprimer 2>/dev/null || echo "")

if [[ "$ENGINE" == "mfeprimer" || "$ENGINE" == "both" ]]; then
if [[ -n "$MFE_BIN" ]]; then
    if mfe_index_present && [[ "$OVERWRITE" != "true" && "$OVERWRITE" != "True" ]]; then
        echo "[build_indices] [$GENOME_NAME] MFEprimer index exists — skip"
    else
        # Stage indexing OFF /mnt/c/. NTFS-via-WSL2 is 10-100x slower than
        # ext4 for the many small writes mfeprimer makes; on /mnt/c/ even a
        # 41 MB transcript FASTA can take 10+ minutes and look "stuck."
        # We index in /tmp/ (Linux tmpfs/ext4) then symlink results back.
        if [[ "$GENOME" == /mnt/* ]]; then
            stage_dir="/tmp/eggplant_mfe_index/$GENOME_NAME"
            mkdir -p "$stage_dir"
            staged_fa="$stage_dir/$(basename "$GENOME")"
            if [[ ! -f "$staged_fa" || "$(stat -c %s "$staged_fa" 2>/dev/null || echo 0)" != "$(stat -c %s "$GENOME")" ]]; then
                echo "[build_indices] [$GENOME_NAME] staging FASTA → $stage_dir (Linux fs)"
                cp -f "$GENOME" "$staged_fa"
            fi
            target_fa="$staged_fa"
        else
            ln -sf "$GENOME" "$ufm_link"
            target_fa="$ufm_link"
        fi

        # Run mfeprimer index with FULL output capture so silent failures
        # become visible. Without this, set -e propagates a bare non-zero
        # exit and the orchestrator dies with no error printed.
        log_file="${target_fa}.index.log"
        echo "[build_indices] [$GENOME_NAME] $MFE_BIN index -i $target_fa -k $KMER  (GOMAXPROCS=$THREADS)"
        echo "[build_indices] [$GENOME_NAME] live output → $log_file"
        set +e
        GOMAXPROCS="$THREADS" "$MFE_BIN" index -i "$target_fa" -k "$KMER" -f 2>&1 | tee "$log_file"
        rc=${PIPESTATUS[0]}
        set -e
        if (( rc != 0 )); then
            echo "[build_indices] [$GENOME_NAME] ERROR: mfeprimer index exited $rc" >&2
            echo "[build_indices] last 30 lines of $log_file:" >&2
            tail -n 30 "$log_file" >&2 || true
            exit "$rc"
        fi

        # Stage results back into OUTDIR (only when staging was used)
        if [[ "$target_fa" == /tmp/* ]]; then
            for f in "$target_fa" "${target_fa}.uni" "${target_fa}.ufm" "${target_fa}.fai" "${target_fa}.2bit"; do
                [[ -e "$f" ]] && ln -sf "$f" "$OUTDIR/$(basename "$f")"
            done
            echo "[build_indices] [$GENOME_NAME] indexed; symlinks placed in $OUTDIR"
        fi
    fi
else
    echo "[build_indices] WARN: mfeprimer not found — install via: bash $SCRIPT_DIR/download_mfeprimer.sh" >&2
fi
fi

# ── isPcr .ooc (only if binary present) ────────────────────────────────────
if [[ "$ENGINE" == "ispcr" || "$ENGINE" == "both" ]]; then
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
fi
