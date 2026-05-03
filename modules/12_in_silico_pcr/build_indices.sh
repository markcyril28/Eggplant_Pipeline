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
    # MFEprimer v4.x writes <fasta>.primerqc (the big k-mer hash);
    # v3.x writes <fasta>.uni;  pre-v3 writes <fasta>.ufm.
    # Treat any of them as proof the index exists.
    [[ -f "${ufm_link}.primerqc" || -f "${ufm_link}.uni" || -f "${ufm_link}.ufm" ]]
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
        # Stage on Linux ext4 if input lives on the Windows-mounted /mnt/.
        # MFEprimer index writes a tempfile in CWD then atomically renames
        # it onto the destination. If CWD is on /mnt/c (9p/NTFS) and the
        # destination is on /tmp (ext4), Linux rename() fails with
        # "invalid cross-device link" — so we MUST cd into the stage dir.
        # We ALSO normalize the FASTA to one-line-per-sequence on staging:
        # mfeprimer v4.x rejects FASTAs with inconsistent line widths inside
        # a single sequence with the message "different line length in
        # sequence: <id>" — and worse, exits 0 anyway, so a bad FASTA
        # silently produces no .uni index. Single-line sequences are
        # trivially uniform and avoid the whole problem.
        if [[ "$GENOME" == /mnt/* ]]; then
            stage_dir="/tmp/eggplant_mfe_index/$GENOME_NAME"
            mkdir -p "$stage_dir"
            staged_fa="$stage_dir/$(basename "$GENOME")"
            src_size=$(stat -c %s "$GENOME")
            dst_size=$(stat -c %s "$staged_fa" 2>/dev/null || echo 0)
            # Always re-stage if missing OR markedly different size from
            # source — normalization shrinks files modestly so we tolerate
            # up to 10% shrink as "already normalized."
            need_stage=1
            if [[ -s "$staged_fa" ]]; then
                ratio=$(( dst_size * 100 / src_size ))
                (( ratio >= 85 && ratio <= 105 )) && need_stage=0
            fi
            if (( need_stage )); then
                echo "[build_indices] [$GENOME_NAME] normalizing + staging FASTA ($((src_size/1024/1024)) MB) → $stage_dir"
                # awk: each sequence on a single line, removes any whitespace.
                awk 'BEGIN{seq=""}
                     /^>/{if(seq!="") print seq; print; seq=""; next}
                     {gsub(/[[:space:]]/,"",$0); seq=seq $0}
                     END{if(seq!="") print seq}' "$GENOME" > "$staged_fa"
                if [[ ! -s "$staged_fa" ]]; then
                    echo "[build_indices] [$GENOME_NAME] ERROR: FASTA normalization produced empty file" >&2
                    exit 1
                fi
                normalized_size=$(stat -c %s "$staged_fa")
                echo "[build_indices] [$GENOME_NAME] normalized: $((normalized_size/1024/1024)) MB"
            fi
            run_dir="$stage_dir"
            relative_fa="$(basename "$staged_fa")"
            target_fa="$staged_fa"
        else
            ln -sf "$GENOME" "$ufm_link"
            run_dir="$OUTDIR"
            relative_fa="$(basename "$ufm_link")"
            target_fa="$ufm_link"
        fi

        # Memory cap: mfeprimer defaults to 70% of host RAM. On WSL2 with
        # a 16 GB ceiling that is ~11 GB and OOM-kills under tee buffering.
        # 50% (~8 GB) leaves headroom for the orchestrator and tee.
        log_file="${target_fa}.index.log"
        echo "[build_indices] [$GENOME_NAME] cd $run_dir && $MFE_BIN index -i $relative_fa -k $KMER -c $THREADS -m 50 -f"
        echo "[build_indices] [$GENOME_NAME] live output → $log_file"
        set +e
        ( cd "$run_dir" \
          && GOMAXPROCS="$THREADS" "$MFE_BIN" index \
              -i "$relative_fa" -k "$KMER" -c "$THREADS" -m 50 -f \
        ) 2>&1 | tee "$log_file"
        rc=${PIPESTATUS[0]}
        set -e
        # mfeprimer v4.x can return 0 even on indexing failure (e.g. when the
        # FASTA has variable line widths). Verify a real index file exists
        # and is non-empty regardless of exit code.
        index_ok=0
        for ext in primerqc uni ufm; do
            [[ -s "${target_fa}.${ext}" ]] && index_ok=1
        done
        if (( rc != 0 )) || (( index_ok == 0 )); then
            echo "[build_indices] [$GENOME_NAME] ERROR: mfeprimer index failed (rc=$rc, no .primerqc/.uni/.ufm)" >&2
            echo "[build_indices] tee log ($log_file):" >&2
            tail -n 40 "$log_file" >&2 || true
            if [[ -s "${target_fa}.log" ]]; then
                echo "[build_indices] mfeprimer internal log (${target_fa}.log):" >&2
                tail -n 20 "${target_fa}.log" >&2
            fi
            exit 1
        fi

        # Symlink the indexed FASTA + companion files back into OUTDIR so
        # mfeprimer_run.sh discovers them via the standard search regardless
        # of where indexing happened. v4.x produces .primerqc, .primerqc.fai,
        # .fai, .json; older versions produce .uni or .ufm.
        if [[ "$target_fa" == /tmp/* ]]; then
            for f in "$target_fa" "${target_fa}.primerqc" "${target_fa}.primerqc.fai" \
                     "${target_fa}.uni" "${target_fa}.ufm" "${target_fa}.fai" \
                     "${target_fa}.json" "${target_fa}.2bit"; do
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
