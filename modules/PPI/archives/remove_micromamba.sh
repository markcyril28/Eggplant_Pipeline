#!/bin/bash
set -euo pipefail

echo "=== Removing Micromamba ==="

# 1. Remove micromamba binary
MICROMAMBA_BIN="${HOME}/miniconda3/bin/micromamba"
if [ -f "$MICROMAMBA_BIN" ]; then
    rm "$MICROMAMBA_BIN"
    echo "Removed binary: $MICROMAMBA_BIN"
else
    echo "Micromamba binary not found at $MICROMAMBA_BIN"
fi

# 2. Remove micromamba environments and caches
MICROMAMBA_DIR="$HOME/.micromamba"
if [ -d "$MICROMAMBA_DIR" ]; then
    rm -rf "$MICROMAMBA_DIR"
    echo "Removed Micromamba environment/cache folder: $MICROMAMBA_DIR"
else
    echo "No Micromamba environment/cache folder found at $MICROMAMBA_DIR"
fi

# 3. Remove micromamba shell hooks from .bashrc and .zshrc
for shellfile in "$HOME/.bashrc" "$HOME/.zshrc"; do
    if [ -f "$shellfile" ]; then
        sed -i.bak '/micromamba/d' "$shellfile"
        echo "Removed Micromamba hooks from $shellfile (backup saved as $shellfile.bak)"
    fi
done

echo "=== Micromamba removal complete ==="
