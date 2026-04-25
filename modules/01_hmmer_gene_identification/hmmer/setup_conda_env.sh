#!/bin/bash
set -euo pipefail

# ================================================================
# Conda Environment Setup for HMMER Gene Identification Workflow
# ================================================================
# This script creates and configures a conda environment with all
# required dependencies for the HMMER-based gene identification pipeline.
# ================================================================

ENV_NAME="hmmer"

echo "============================================================"
echo "Setting up Conda Environment: $ENV_NAME"
echo "============================================================"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "[ERROR] Conda is not installed or not in PATH."
    echo "        Please install Miniconda or Anaconda first:"
    echo "        https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Remove existing environment if it exists (optional - comment out to skip)
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[INFO] Environment '$ENV_NAME' already exists."
    read -p "Do you want to remove and recreate it? (y/n): " response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        echo "[INFO] Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
    else
        echo "[INFO] Keeping existing environment. Activating..."
        echo ""
        echo "To activate the environment, run:"
        echo "    conda activate $ENV_NAME"
        exit 0
    fi
fi

echo "[INFO] Creating new conda environment: $ENV_NAME"
conda create -n "$ENV_NAME" python=3.10 -y

echo "[INFO] Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "[INFO] Installing bioinformatics tools from bioconda..."
conda install -y -c bioconda -c conda-forge \
    hmmer \
    cd-hit \
    mafft \
    trimal \
    wget

echo "[INFO] Verifying installations..."
echo ""
echo "------------------------------------------------------------"
echo "Installed Tool Versions:"
echo "------------------------------------------------------------"

# Check HMMER
if command -v hmmsearch &> /dev/null; then
    echo "HMMER:    $(hmmsearch -h | head -2 | tail -1)"
else
    echo "HMMER:    [NOT INSTALLED]"
fi

# Check CD-HIT
if command -v cd-hit &> /dev/null; then
    echo "CD-HIT:   $(cd-hit -h 2>&1 | head -1)"
else
    echo "CD-HIT:   [NOT INSTALLED]"
fi

# Check MAFFT
if command -v mafft &> /dev/null; then
    echo "MAFFT:    $(mafft --version 2>&1)"
else
    echo "MAFFT:    [NOT INSTALLED]"
fi

# Check trimAl
if command -v trimal &> /dev/null; then
    echo "trimAl:   $(trimal --version 2>&1 | head -1)"
else
    echo "trimAl:   [NOT INSTALLED]"
fi

# Check wget
if command -v wget &> /dev/null; then
    echo "wget:     $(wget --version | head -1)"
else
    echo "wget:     [NOT INSTALLED]"
fi

echo "------------------------------------------------------------"
echo ""
echo "============================================================"
echo "[SUCCESS] Environment '$ENV_NAME' is ready!"
echo "============================================================"
echo ""
echo "To activate this environment, run:"
echo "    conda activate $ENV_NAME"
echo ""
echo "To run the HMMER workflow:"
echo "    conda activate $ENV_NAME"
echo "    bash HMMER_script_v2.sh"
echo ""
echo "To deactivate when done:"
echo "    conda deactivate"
echo ""
