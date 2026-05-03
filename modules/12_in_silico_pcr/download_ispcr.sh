#!/bin/bash
# ============================================================================
# Download UCSC isPcr (+ blat for .ooc generation) into modules/12_in_silico_pcr/bin/
# ----------------------------------------------------------------------------
# Reference: https://hgdownload.soe.ucsc.edu/admin/exe/
# License:   Free for non-commercial use (kent source)
#
# Usage:
#   bash modules/12_in_silico_pcr/download_ispcr.sh
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"
mkdir -p "$BIN_DIR"

# Detect platform; UCSC publishes Linux x86_64 and macOS x86_64 binaries.
case "$(uname -s)" in
    Linux*)   PLATFORM="linux.x86_64" ;;
    Darwin*)  PLATFORM="macOSX.x86_64" ;;
    *)
        echo "[download_ispcr] No native UCSC binary for $(uname -s)." >&2
        echo "  On Windows, run inside WSL (Linux subsystem)." >&2
        exit 1 ;;
esac

BASE="https://hgdownload.soe.ucsc.edu/admin/exe/${PLATFORM}"

for tool in isPcr blat; do
    target="$BIN_DIR/$tool"
    if [[ -x "$target" ]]; then
        echo "[download_ispcr] $tool already present: $target"
        continue
    fi
    echo "[download_ispcr] Fetching $tool from $BASE/$tool"
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$target" "$BASE/$tool"
    elif command -v curl &>/dev/null; then
        curl -fL -o "$target" "$BASE/$tool"
    else
        echo "[download_ispcr] Need wget or curl on PATH." >&2
        exit 1
    fi
    chmod +x "$target"
done

echo
echo "Installed:"
ls -lh "$BIN_DIR"
echo
echo "Verify:"
"$BIN_DIR/isPcr" 2>&1 | head -3 || true
