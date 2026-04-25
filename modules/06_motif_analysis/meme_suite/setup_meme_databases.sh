#!/bin/bash
# ============================================================================
# Module: MEME Suite Motif Database Setup
# ============================================================================
# Extracts the bundled motif_databases.12.27.tgz to the configured location
# (default: 2_INPUTS/meme_motif_databases/).
#
# Run once after cloning the repository:
#   bash modules/06_motif_analysis/meme_suite/setup_meme_databases.sh
#
# With a custom target directory:
#   bash setup_meme_databases.sh --outdir /path/to/databases
#
# With explicit tarball path:
#   bash setup_meme_databases.sh --tarball /path/to/motif_databases.12.27.tgz
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../logging/logging_utils.sh"

# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DEFAULT_TARBALL="$SCRIPT_DIR/../assets/motif_databases.12.27.tgz"
DEFAULT_OUTDIR="$PIPELINE_DIR/2_INPUTS/meme_motif_databases"
OVERWRITE=false
# ===============================================================

TARBALL="$DEFAULT_TARBALL"
OUTDIR="$DEFAULT_OUTDIR"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tarball)  TARBALL="$2";  shift 2 ;;
        --outdir)   OUTDIR="$2";   shift 2 ;;
        --overwrite) OVERWRITE=true; shift ;;
        *) log_error "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate tarball
# ---------------------------------------------------------------------------
if [[ ! -f "$TARBALL" ]]; then
    log_error "Motif database tarball not found: $TARBALL"
    log_error "Expected: modules/06_motif_analysis/motif_databases.12.27.tgz"
    exit 1
fi

# ---------------------------------------------------------------------------
# Skip if already extracted and not overwriting
# ---------------------------------------------------------------------------
MARKER="$OUTDIR/.extracted"
if [[ -f "$MARKER" && "$OVERWRITE" != true ]]; then
    log_info "Motif databases already extracted to: $OUTDIR"
    log_info "Use --overwrite to re-extract."
    exit 0
fi

# ---------------------------------------------------------------------------
# Extract
# ---------------------------------------------------------------------------
mkdir -p "$OUTDIR"
log_step "Extracting motif databases"
log_info "Source: $TARBALL"
log_info "Target: $OUTDIR"

# Log tarball size before extraction
tarball_size_mb=$(awk "BEGIN { printf \"%.1f\", $(stat -c%s "$TARBALL" 2>/dev/null || stat -f%z "$TARBALL" 2>/dev/null || echo 0) / 1048576 }")
log_info "Tarball size: ${tarball_size_mb} MB"

# The tarball contains a 'motif_databases/' top-level directory;
# strip it so files land directly in OUTDIR.
# --warning=no-unknown-keyword suppresses macOS xattr metadata warnings
# (LIBARCHIVE.xattr.com.apple.quarantine, etc.) from tarballs created on macOS.
tar_start=$SECONDS
if ! tar -xzf "$TARBALL" \
    --strip-components=1 \
    --warning=no-unknown-keyword \
    -C "$OUTDIR" 2>&1; then
    log_error "tar extraction failed for: $TARBALL"
    exit 1
fi
tar_elapsed=$(( SECONDS - tar_start ))

touch "$MARKER"

# Log extraction summary
file_count=$(find "$OUTDIR" -type f | wc -l)
dir_size_mb=$(awk "BEGIN { printf \"%.1f\", $(du -sk "$OUTDIR" 2>/dev/null | cut -f1) / 1024 }")

log_step "Motif databases extracted successfully"
log_info "Elapsed: ${tar_elapsed}s | Files: ${file_count} | Size: ${dir_size_mb} MB"
log_info "Plant-relevant databases:"
log_info "  JASPAR: $OUTDIR/JASPAR/"
log_info "  ARABD:  $OUTDIR/ARABD/"
log_info ""
log_info "Key plant databases ready for FIMO / TOMTOM:"
log_info "  JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme"
log_info "  ARABD/ArabidopsisDAPv1.meme"
