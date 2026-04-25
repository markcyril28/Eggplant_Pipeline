#!/bin/bash
# Install MEGA-CC 12.0.14 and common alignment tools

set -euo pipefail

# Update and fix broken dependencies
sudo apt-get update && sudo apt-get install -f -y

# Use conda/bioconda for alignment tools
if command -v mamba >/dev/null 2>&1; then
    PKG_MGR="mamba"
elif command -v conda >/dev/null 2>&1; then
    PKG_MGR="conda"
else
    echo "Error: conda or mamba not found. Install Miniconda/Mambaforge and rerun." >&2
    exit 1
fi

$PKG_MGR install -y -c conda-forge -c bioconda \
    muscle clustalw mafft t-coffee probcons iqtree
