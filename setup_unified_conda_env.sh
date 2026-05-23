#!/bin/bash
# Hello World
# ============================================================================
# Unified Conda Environment Setup for the Eggplant Pipeline
# ============================================================================
# This script creates a single conda environment with ALL dependencies needed
# across the entire Eggplant_Pipeline project, consolidating from:
#
#   - HMMER Gene Identification    (formerly: "hmmer"     env)
#   - BLAST Ortholog Analysis      (formerly: "blast"     env)
#   - PlantCARE Analysis           (formerly: "plantcare" env)
#   - Phylogenetic Analysis        (formerly: conda/local installs)
#   - Protein Structure Analysis   (PyMOL + BioPython)
#   - CRISPR Off-Target Analysis   (Cas-OFFinder)
#
# Usage:
#   bash setup_unified_conda_env.sh            # Create environment
#   bash setup_unified_conda_env.sh --remove   # Remove and recreate
#   bash setup_unified_conda_env.sh --update   # Update existing environment
#
# After setup, activate with:
#   conda activate egg
# ============================================================================

set -euo pipefail

# ======================== Configuration ========================
ENV_NAME="egg"
PYTHON_VERSION="3.10"

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

print_step() {
    echo -e "${GREEN}[STEP]${NC} $1"
}

print_info() {
    echo -e "${YELLOW}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[OK]${NC} $1"
}

# ======================== Parse Arguments ========================
MODE="create"
if [[ "${1:-}" == "--remove" ]]; then
    MODE="remove"
elif [[ "${1:-}" == "--update" ]]; then
    MODE="update"
fi

# ======================== Pre-flight Checks ========================
print_header "Eggplant Pipeline - Unified Conda Environment Setup"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    print_error "Conda is not installed or not in PATH."
    echo "        Please install Miniconda or Anaconda first:"
    echo "        https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check for mamba (faster dependency resolution)
if command -v mamba &> /dev/null; then
    PKG_MGR="mamba"
    print_info "Using mamba for faster package installation."
else
    PKG_MGR="conda"
    print_info "Using conda (install mamba for faster setup: conda install -n base -c conda-forge mamba)"
fi

# ======================== Handle Existing Environment ========================
if [[ "$MODE" == "update" ]]; then
    if ! conda env list | grep -q "^${ENV_NAME} "; then
        print_error "Environment '${ENV_NAME}' does not exist. Cannot update."
        echo "        Run without flags first to create the environment:"
        echo "        bash setup_unified_conda_env.sh"
        exit 1
    fi
    print_info "Updating existing environment '${ENV_NAME}'..."

    # Source conda for activation in script context
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"

    print_success "Environment activated for update."

elif conda env list | grep -q "^${ENV_NAME} "; then
    if [[ "$MODE" == "remove" ]]; then
        print_info "Removing existing environment '${ENV_NAME}'..."
        conda env remove -n "$ENV_NAME" -y
    else
        print_info "Environment '${ENV_NAME}' already exists."
        echo ""
        echo "  [u] Update packages in the existing environment"
        echo "  [r] Remove and recreate from scratch"
        echo "  [n] Cancel (keep as-is)"
        echo ""
        read -p "Choose an option (u/r/n): " response
        case "${response,,}" in
            u)
                MODE="update"
                source "$(conda info --base)/etc/profile.d/conda.sh"
                conda activate "$ENV_NAME"
                print_success "Environment activated for update."
                ;;
            r)
                print_info "Removing existing environment..."
                conda env remove -n "$ENV_NAME" -y
                ;;
            *)
                echo ""
                echo "No changes made. To activate the existing environment:"
                echo "    conda activate $ENV_NAME"
                exit 0
                ;;
        esac
    fi
fi

# ======================== Create Environment ========================
if [[ "$MODE" != "update" ]]; then
    print_header "Step 1/6: Creating Conda Environment"
    print_step "Creating environment '${ENV_NAME}' with Python ${PYTHON_VERSION}..."

    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y

    # Source conda for activation in script context
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"

    print_success "Environment created and activated."
else
    print_header "Step 1/6: Skipped (updating existing environment)"
fi

# ======================== Configure Channels ========================
if [[ "$MODE" != "update" ]]; then
    print_header "Step 2/6: Configuring Conda Channels"

    conda config --env --add channels defaults 2>/dev/null || true
    conda config --env --add channels bioconda
    conda config --env --add channels conda-forge
    conda config --env --set channel_priority strict

    print_success "Channels configured: conda-forge > bioconda > defaults (strict priority)"
else
    print_header "Step 2/6: Skipped (channels already configured)"
fi

# ======================== Install All Conda Packages ========================
# All packages are installed in a single call so dependency resolution runs
# only once instead of once per tool group (the primary speed bottleneck).
if [[ "$MODE" == "update" ]]; then
    print_header "Step 3/6: Updating All Conda Packages"
else
    print_header "Step 3/6: Installing All Conda Packages"
fi

print_step "Resolving and installing all conda packages in a single pass..."
CONDA_PACKAGES=(
    # Gene identification
    hmmer cd-hit blast gffread
    # Sequence alignment
    mafft muscle clustalo clustalw t-coffee probcons trimal
    # Phylogenetics
    iqtree fasttree raxml-ng modeltest-ng
    # Genomic utilities
    samtools bedtools seqtk
    # Motif analysis
    meme
    # CRISPR off-target
    cas-offinder
    # General utilities
    wget pandoc
    # Protein structure
    pymol-open-source emboss
    # R base + only the CRAN packages we need (r-essentials is bloated)
    r-base r-ggplot2 r-pheatmap r-rcolorbrewer r-argparse r-tidyr r-dplyr r-circlize
    # Phylogenetic tree visualization (ggtree + dependencies)
    bioconductor-ggtree bioconductor-treeio
    # Tree comparison: tanglegram (phytools), normalized RF distance (ape), fallback panels (patchwork)
    r-ape r-phytools r-patchwork
)

if [[ "$MODE" == "update" ]]; then
    # 'update --all' upgrades every installed package; works for both conda and mamba
    $PKG_MGR update -y --all \
        -c conda-forge -c bioconda --strict-channel-priority
    # Also ensure any newly added packages are present
    $PKG_MGR install -y \
        -c conda-forge -c bioconda --strict-channel-priority \
        "${CONDA_PACKAGES[@]}"
else
    $PKG_MGR install -y \
        -c conda-forge -c bioconda --strict-channel-priority \
        "${CONDA_PACKAGES[@]}"
fi

print_success "All conda packages installed."

# ======================== Install Python Packages ========================
if [[ "$MODE" == "update" ]]; then
    print_header "Step 4/6: Updating Python Packages"
    PIP_UPDATE_FLAG="--upgrade"
else
    print_header "Step 4/6: Installing Python Packages"
    PIP_UPDATE_FLAG=""
fi

print_step "Installing Python scientific libraries..."
pip install --quiet $PIP_UPDATE_FLAG \
    pandas \
    numpy \
    biopython \
    tomli \
    markdown \
    weasyprint

print_success "All Python packages installed."

# ======================== Install R and R Packages ========================
if [[ "$MODE" == "update" ]]; then
    print_header "Step 5/6: Updating R and R Packages"
else
    print_header "Step 5/6: Installing R and R Packages"
fi

R_UPDATE_MODE="$MODE"
print_step "Installing Bioconductor R packages (CRAN packages already installed via conda)..."
Rscript -e "
options(repos = c(CRAN = 'https://cran.r-project.org'))
update_mode <- '$R_UPDATE_MODE' == 'update'

if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')

# Resolve the Bioconductor release compatible with the installed R version.
# BiocManager may default to the latest release (e.g. 3.22) which can require
# a newer R than is present; map explicitly to avoid that error.
r_ver <- paste0(R.Version()\$major, '.', floor(as.numeric(R.Version()\$minor)))
bioc_ver <- switch(r_ver,
    '4.5' = , '4.6' = '3.22',
    '4.4' = '3.20',
    '4.3' = '3.18',
    '4.2' = '3.16'
    # NULL for unrecognized R versions: let BiocManager auto-select
)
if (!is.null(bioc_ver)) {
    cat(sprintf('  R %s -> Bioconductor %s\n', r_ver, bioc_ver))
} else {
    cat(sprintf('  R %s not in version table; BiocManager will auto-select\n', r_ver))
}

# Align any already-installed Bioconductor packages to the target version.
# This handles the case where a prior run installed packages at a newer
# Bioconductor release (e.g. 3.22); BiocManager requires this downgrade
# step before it will install further packages at the older version.
if (!is.null(bioc_ver)) {
    cat(sprintf('  Aligning installed Bioconductor packages to version %s...\n', bioc_ver))
    BiocManager::install(version = bioc_ver, ask = FALSE, force = TRUE)
}

bioc_pkgs <- c('ComplexHeatmap')
missing <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0 || update_mode) {
    cat(sprintf('  Installing/updating: %s\n', paste(if (update_mode) bioc_pkgs else missing, collapse = ', ')))
    pkgs_to_install <- if (update_mode) bioc_pkgs else missing
    if (!is.null(bioc_ver)) {
        BiocManager::install(pkgs_to_install, version = bioc_ver, update = update_mode, ask = FALSE)
    } else {
        BiocManager::install(pkgs_to_install, update = update_mode, ask = FALSE)
    }
} else {
    cat('  ComplexHeatmap already installed\n')
}
cat('All R packages ready.\n')
"

print_success "R environment configured."

# ======================== Wire In Silico PCR bin/ into env PATH ============
# MFEprimer + isPcr are downloaded as standalone binaries into
# modules/12_in_silico_pcr/bin/. Put that directory on PATH automatically
# whenever the user runs `conda activate egg`.
print_header "Step 5.5/6: Wiring In Silico PCR binaries into conda env PATH"

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PCR_BIN_DIR="$PIPELINE_DIR/modules/12_in_silico_pcr/bin"
ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "$ACTIVATE_DIR/in_silico_pcr_path.sh" <<EOF
#!/bin/bash
# Auto-generated by setup_unified_conda_env.sh — prepends stage 12 binaries
# (MFEprimer, isPcr, blat) to PATH so they run from any working directory.
export _EGG_PCR_PATH_BACKUP="\$PATH"
export PATH="$PCR_BIN_DIR:\$PATH"
EOF

cat > "$DEACTIVATE_DIR/in_silico_pcr_path.sh" <<'EOF'
#!/bin/bash
# Auto-generated by setup_unified_conda_env.sh — restore original PATH on
# `conda deactivate`.
if [[ -n "${_EGG_PCR_PATH_BACKUP:-}" ]]; then
    export PATH="$_EGG_PCR_PATH_BACKUP"
    unset _EGG_PCR_PATH_BACKUP
fi
EOF

chmod +x "$ACTIVATE_DIR/in_silico_pcr_path.sh" "$DEACTIVATE_DIR/in_silico_pcr_path.sh"

# Also export it for the remainder of THIS script run, so verification picks
# up MFEprimer/isPcr if the user already downloaded the binaries.
export PATH="$PCR_BIN_DIR:$PATH"

print_success "Activated:  PATH now includes $PCR_BIN_DIR"
print_info    "Persistent: takes effect on every  conda activate $ENV_NAME"
print_info    "Manual:     export PATH=\"$PCR_BIN_DIR:\$PATH\""

# ======================== Verify Installations ========================
print_header "Step 6/6: Verifying Installations"
print_info "Checking tool versions..."

echo ""
echo "--- Bioinformatics Tools ---"

declare -A TOOLS=(
    ["HMMER"]="hmmsearch -h | head -2 | tail -1"
    ["CD-HIT"]="cd-hit -h 2>&1 | head -1"
    ["BLAST+"]="blastn -version | head -1"
    ["MAFFT"]="mafft --version 2>&1"
    ["MUSCLE"]="muscle --version 2>&1 | head -1"
    ["Clustal Omega"]="clustalo --version 2>&1"
    ["ClustalW"]="clustalw -help 2>&1 | head -1"
    ["T-Coffee"]="t_coffee -version 2>&1 | head -1"
    ["PROBCONS"]="probcons 2>&1 | head -1"
    ["trimAl"]="trimal -version 2>&1 | head -1"
    ["IQ-TREE"]="iqtree --version 2>&1 | head -1"
    ["RAxML-NG"]="raxml-ng --version 2>&1 | head -1"
    ["ModelTest-NG"]="modeltest-ng --version 2>&1 | head -1"
    ["FastTree"]="FastTree 2>&1 | head -1"
    ["SAMtools"]="samtools --version 2>&1 | head -1"
    ["BEDtools"]="bedtools --version 2>&1"
    ["SeqTK"]="seqtk 2>&1 | head -1"
    ["Cas-OFFinder"]="cas-offinder 2>&1 | head -1"
    ["MFEprimer"]="mfeprimer version 2>&1 | head -1"
    ["isPcr (UCSC)"]="isPcr 2>&1 | head -1"
    ["MEME"]="meme --version 2>&1"
    ["gffread"]="gffread --version 2>&1 | head -1"
    ["EMBOSS transeq"]="transeq -version 2>&1 | head -1"
    ["wget"]="wget --version 2>&1 | head -1"
)

# Run all version checks in parallel for speed
TMPDIR_VERIFY=$(mktemp -d)
for tool in "${!TOOLS[@]}"; do
    (
        set +eo pipefail
        version=$(eval "${TOOLS[$tool]}" 2>/dev/null)
        [[ -z "$version" ]] && version="[NOT FOUND]"
        printf "  %-25s %s\n" "$tool:" "$version" > "$TMPDIR_VERIFY/$tool"
    ) &
done

# Python + R checks in parallel too
(
    echo "" > "$TMPDIR_VERIFY/_python"
    echo "--- Python Packages ---" >> "$TMPDIR_VERIFY/_python"
    python -c "
import pandas; print(f'  pandas:    {pandas.__version__}')
import numpy; print(f'  numpy:     {numpy.__version__}')
import Bio; print(f'  biopython: {Bio.__version__}')
" >> "$TMPDIR_VERIFY/_python" 2>/dev/null || echo "  [Some Python packages not found]" >> "$TMPDIR_VERIFY/_python"
    python -c "import pymol; print(f'  pymol:     available')" >> "$TMPDIR_VERIFY/_python" 2>/dev/null || echo "  pymol:     [available via pymol command]" >> "$TMPDIR_VERIFY/_python"
) &

(
    echo "" > "$TMPDIR_VERIFY/_r"
    echo "--- R Packages ---" >> "$TMPDIR_VERIFY/_r"
    Rscript -e "
pkgs <- c('ComplexHeatmap', 'pheatmap', 'ggplot2',
          'circlize', 'RColorBrewer', 'argparse', 'tidyr', 'dplyr',
          'ggtree', 'treeio', 'ape', 'phytools', 'patchwork')
for (pkg in pkgs) {
    v <- tryCatch(as.character(packageVersion(pkg)), error = function(e) 'NOT FOUND')
    cat(sprintf('  %-20s %s\n', paste0(pkg, ':'), v))
}
" >> "$TMPDIR_VERIFY/_r" 2>/dev/null || echo "  [Some R packages not found]" >> "$TMPDIR_VERIFY/_r"
) &

wait

# Print collected results
for tool in "${!TOOLS[@]}"; do
    cat "$TMPDIR_VERIFY/$tool" 2>/dev/null
done
cat "$TMPDIR_VERIFY/_python"
cat "$TMPDIR_VERIFY/_r"
rm -rf "$TMPDIR_VERIFY"

# ======================== Summary ========================
print_header "Setup Complete!"

cat << 'EOF'

  Environment Name:  egg
  Python Version:    3.10

  -------------------------------------------------------
  CONSOLIDATED TOOLS BY PIPELINE MODULE
  -------------------------------------------------------

  Gene Identification:
    hmmer, cd-hit, blast+, gffread

  Sequence Alignment:
    mafft, muscle, clustalo, clustalw, t-coffee, probcons

  Alignment Trimming:
    trimal

  Phylogenetic Analysis:
    iqtree, raxml-ng, fasttree, modeltest-ng
    R: ggtree, treeio, ape, phytools, patchwork
    (MEGA-CC requires separate .deb install - see note below)

  Promoter / Regulatory Analysis (PlantCARE):
    samtools, bedtools, seqtk
    R: ComplexHeatmap, circlize, argparse, tidyr, dplyr

  Protein Structure & Visualization:
    pymol-open-source, emboss (transeq), biopython
    (AlphaFold3 runs via Google DeepMind server - not local)
    (SWISS-MODEL and STRING are web-based tools)

  CRISPR Off-Target:
    cas-offinder

  In Silico PCR (stage 12):
    Both engines distribute precompiled binaries (no conda/pip packages).
    Fetch them with:
        bash modules/12_in_silico_pcr/download_mfeprimer.sh
        bash modules/12_in_silico_pcr/download_ispcr.sh

  General Utilities:
    wget, pandas, numpy, tomli (TOML config parser for Python <3.11)

  -------------------------------------------------------
  TOOLS NOT INCLUDED (require separate installation):
  -------------------------------------------------------

  * MEGA-CC 12.0.14:
      Install via .deb package (Linux only):
      sudo dpkg -i 1_CONFIG_FILES/mega-cc_12.0.14-1_amd64_beta.deb
      Or via Docker:
      docker build -t megacc -f 9_RESULTS/GRF_GIF/3_BLAST_Alignment_and_Phylogenetic_Analysis/Dockerfile .

  * In silico PCR engines (stage 12):
      MFEprimer-3.0 (Go binary, GitHub releases):
          bash modules/12_in_silico_pcr/download_mfeprimer.sh
      UCSC isPcr (kent source, free for non-commercial use):
          bash modules/12_in_silico_pcr/download_ispcr.sh
      Both target Linux/macOS; on Windows run inside WSL.

  * Web-based / External tools (no local install needed):
      - AlphaFold3 Server (alphafoldserver.com)
      - SWISS-MODEL (swissmodel.expasy.org)
      - STRING (string-db.org)
      - PlantCARE (bioinformatics.psb.ugent.be/webtools/plantcare)
      - MEME Suite (meme-suite.org)
      - Gene Structure Display Server (gsds.gao-lab.org)
      - MG2C (mg2c.iask.in)

  -------------------------------------------------------

  To activate:    conda activate egg
  To deactivate:  conda deactivate
  To update:      bash setup_unified_conda_env.sh --update
  To remove:      conda env remove -n egg

EOF
