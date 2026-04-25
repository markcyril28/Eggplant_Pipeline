#!/bin/bash
set -euo pipefail

# ========================================
# GRF-GIF Project Conda Environment Setup
# ========================================
# Creates a conda environment with all dependencies needed for
# GRF-GIF gene family analysis including BLAST and Python tools.

# Configuration
readonly ENV_NAME="blast"
readonly PYTHON_VERSION="3.9"

echo "=============================== GRF-GIF Environment Setup ==============================="
echo "Environment: $ENV_NAME"
echo "Python: $PYTHON_VERSION"
echo ""

# Check if environment already exists
if conda env list | grep -q "^$ENV_NAME "; then
    echo "Environment '$ENV_NAME' already exists."
    read -p "Update existing environment? [y/N]: " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Updating existing environment..."
    else
        echo "Setup cancelled."
        exit 0
    fi
else
    echo "Creating new environment..."
fi

# Create and setup environment
echo "Creating conda environment..."
conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y

echo "Installing packages..."
conda install -n "$ENV_NAME" -y \
    -c bioconda \
    -c conda-forge \
    blast \
    pandas \
    numpy 

echo ""
echo "=============================== Setup Complete ==============================="
echo "Environment '$ENV_NAME' created successfully!"
echo ""
echo "Usage:"
echo "  Activate:   conda activate $ENV_NAME"
echo "  Deactivate: conda deactivate"
echo "  Remove:     conda env remove -n $ENV_NAME"
echo ""
echo "Installed tools: BLAST+, Python $PYTHON_VERSION, pandas, numpy"