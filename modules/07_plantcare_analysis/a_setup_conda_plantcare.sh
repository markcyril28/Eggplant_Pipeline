#!/bin/bash
set -euo pipefail

# Setup script for eggplant pipeline conda environment
# Creates environment and installs required packages with skip-if-installed logic

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

echo "=== eggplant Pipeline Setup ==="

# Check if conda is installed
if ! command_exists conda; then
    echo "Error: Conda is not installed. Please install conda first."
    exit 1
fi

# Check if the environment already exists
if conda env list | grep -q "^egg "; then
    echo "Conda environment 'egg' already exists."
    read -p "Do you want to update packages? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Setup cancelled."
        exit 0
    fi
else
    echo "Creating conda environment 'egg' with Python 3.8..."
    conda create -n egg python=3.8 -y
fi

# Activate the environment
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate egg

# Install Python packages
echo "--- Python Packages ---"
for pkg in pandas numpy; do
    if python -c "import $pkg" 2>/dev/null; then
        echo "$pkg already installed"
    else
        echo "Installing $pkg..."
        pip install $pkg
    fi
done

# Install bioinformatics tools
echo "--- Bioinformatics Tools ---"
for tool in bedtools seqtk; do
    if command_exists $tool; then
        echo "$tool already installed"
    else
        echo "Installing $tool..."
        conda install -c bioconda $tool -y
    fi
done

# Install R
echo "--- R and R Packages ---"
if ! command_exists R; then
    echo "Installing R..."
    conda install -c conda-forge r-base r-essentials -y
else
    echo "R already installed ($(R --version | head -n1))"
fi

# Install R packages
echo "Installing R packages..."
Rscript -e "
options(repos = c(CRAN = 'https://cran.r-project.org'))

# Function to check and install package
install_if_missing <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf('Installing %s...\\n', pkg))
        if (bioc) {
            if (!requireNamespace('BiocManager', quietly = TRUE)) {
                install.packages('BiocManager')
            }
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else {
            install.packages(pkg)
        }
    } else {
        cat(sprintf('%s already installed\\n', pkg))
    }
}

# Install packages
install_if_missing('BiocManager')
install_if_missing('ComplexHeatmap', bioc = TRUE)
install_if_missing('circlize')
install_if_missing('argparse')
install_if_missing('ggplot2')
install_if_missing('pheatmap')
install_if_missing('RColorBrewer')
install_if_missing('tidyr')
install_if_missing('dplyr')
install_if_missing('grid')

cat('\\n=== All R packages installed ===\\n')
"

echo ""
echo "=== Setup Complete ==="
echo "Conda environment 'egg' is ready."
echo "To activate: conda activate egg"
