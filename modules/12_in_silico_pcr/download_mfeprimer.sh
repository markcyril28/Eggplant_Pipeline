#!/bin/bash
# ============================================================================
# Download MFEprimer-3.0 binary into modules/12_in_silico_pcr/bin/
# ----------------------------------------------------------------------------
# Reference: Wang et al. 2019, Nucleic Acids Research, doi:10.1093/nar/gkz351
# Source:    https://github.com/quwubin/MFEprimer-3.0/releases
# License:   GPL-3.0
#
# MFEprimer is distributed as a precompiled Go binary on GitHub Releases.
# It is NOT on bioconda or PyPI, so installation is by direct download.
#
# Usage:
#   bash modules/12_in_silico_pcr/download_mfeprimer.sh            # latest tag
#   bash modules/12_in_silico_pcr/download_mfeprimer.sh v3.3.1     # pin version
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"
mkdir -p "$BIN_DIR"

REPO="quwubin/MFEprimer-3.0"
PINNED_VERSION="${1:-}"

# ── Detect platform ────────────────────────────────────────────────────────
case "$(uname -s)" in
    Linux*)   OS="linux" ;;
    Darwin*)  OS="darwin" ;;
    *)
        echo "[download_mfeprimer] No native binary for $(uname -s)." >&2
        echo "  On Windows, run inside WSL (Linux subsystem)." >&2
        exit 1 ;;
esac

case "$(uname -m)" in
    x86_64|amd64)  ARCH="amd64" ;;
    aarch64|arm64) ARCH="arm64" ;;
    *)
        echo "[download_mfeprimer] Unsupported architecture: $(uname -m)" >&2
        exit 1 ;;
esac

# ── Resolve release tag ────────────────────────────────────────────────────
if [[ -n "$PINNED_VERSION" ]]; then
    TAG="$PINNED_VERSION"
else
    if command -v curl &>/dev/null; then
        TAG=$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" \
              | grep -E '"tag_name":' | head -1 | sed -E 's/.*"([^"]+)".*/\1/')
    elif command -v wget &>/dev/null; then
        TAG=$(wget -qO- "https://api.github.com/repos/${REPO}/releases/latest" \
              | grep -E '"tag_name":' | head -1 | sed -E 's/.*"([^"]+)".*/\1/')
    fi
    if [[ -z "${TAG:-}" ]]; then
        TAG="v4.2.4"
        echo "[download_mfeprimer] Could not resolve latest tag; falling back to $TAG"
    fi
fi
VERSION="${TAG#v}"
# NOTE: The repo is named MFEprimer-3.0 for historical reasons but distributes
# v4.x binaries since 2022. The release tag drives the asset filename.
echo "[download_mfeprimer] target: MFEprimer ${TAG} (${OS}-${ARCH})"

# ── Download + extract ─────────────────────────────────────────────────────
ASSET="mfeprimer-${VERSION}-${OS}-${ARCH}.gz"
URL="https://github.com/${REPO}/releases/download/${TAG}/${ASSET}"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

target="$BIN_DIR/mfeprimer"
if [[ -x "$target" ]]; then
    have_ver=$("$target" version 2>&1 | head -1 || echo "")
    echo "[download_mfeprimer] existing binary: $target ($have_ver)"
fi

echo "[download_mfeprimer] fetching $URL"
if command -v curl &>/dev/null; then
    curl -fL -o "$TMP/$ASSET" "$URL"
elif command -v wget &>/dev/null; then
    wget -q --show-progress -O "$TMP/$ASSET" "$URL"
else
    echo "[download_mfeprimer] need curl or wget on PATH." >&2
    exit 1
fi

gunzip -c "$TMP/$ASSET" > "$target"
chmod +x "$target"

echo
# v4.x exposes version as a subcommand; fall back to the help banner if the
# subcommand was renamed in a future release.
"$target" version 2>&1 | head -3 || "$target" --help 2>&1 | head -3 || true
echo
echo "Installed: $target"
echo "Add to PATH (optional):  export PATH=\"$BIN_DIR:\$PATH\""
