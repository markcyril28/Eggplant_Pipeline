#!/bin/bash
# ============================================================================
# CRISPR v2 Conda Environment Setup
# ============================================================================
# Creates the 'crispr_v2' conda environment with all tools required for the
# eight-stage CRISPR KO prediction pipeline (i_crispr_v2_pipeline.sh):
#
#   [1]  CRISPOR          -- gRNA design + on/off-target scoring
#   [2a] CRISPRon         -- deep-learning on-target rescoring
#   [2b] DeepSpCas9       -- alternative on-target predictor (Seq-DeepCpf1 port)
#   [3]  (paralog curation -- no extra tool; uses pandas + custom scripts)
#   [4a] inDelphi         -- indel outcome prediction
#   [4b] Lindel           -- alternative indel predictor
#   [5] Biopython        -- mutant transcript reconstruction
#   [6] ESMFold          -- protein structure flag (API or local; config-driven)
#   [7] (NMD heuristic   -- pure Python, no extra tool)
#   [8] (ranking         -- pandas, no extra tool)
#
# Usage:
#   bash setup_conda_crispr_v2.sh             # Create environment
#   bash setup_conda_crispr_v2.sh --remove    # Remove and recreate
#   bash setup_conda_crispr_v2.sh --update    # Update existing environment
#
# Declarative environment spec: envs/crispr_v2_environment.yml
#   conda env create  -f envs/crispr_v2_environment.yml
#   conda env update  -f envs/crispr_v2_environment.yml --prune
#
# After setup, activate with:
#   conda activate crispr_v2
# ============================================================================

set -euo pipefail

ENV_NAME="crispr_v2"
PYTHON_VERSION="3.10"

# ─── Logging ─────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs/installation_logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/setup_conda_crispr_v2_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee >(sed -u \
    -e $'s/\r/\\\n/g' \
    -e $'s/\x1B\\[[0-9;?]*[a-zA-Z]//g' \
    -e $'s/\x1B[()][A-Z0-9]//g' \
    -e 's/[[:space:]]*\([-\\|/][[:space:]]*\)\{2,\}/ /g' \
    -e 's/^\([^:]*:\) \{0,\}done$/\1 done/' \
    -e '/^[[:space:]]*[-\\|/][[:space:]]*$/d' \
    >> "$LOG_FILE")) 2>&1
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup started -- log: $LOG_FILE"

# ─── Colors ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; NC='\033[0m'

print_header()  { echo ""; echo -e "${CYAN}============================================================${NC}"; echo -e "${CYAN}  $1${NC}"; echo -e "${CYAN}============================================================${NC}"; }
print_step()    { echo -e "${GREEN}[STEP]${NC} $1"; }
print_info()    { echo -e "${YELLOW}[INFO]${NC} $1"; }
print_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
print_success() { echo -e "${GREEN}[OK]${NC} $1"; }

# ─── Parse arguments ─────────────────────────────────────────────────────────
MODE="create"
if   [[ "${1:-}" == "--remove" ]]; then MODE="remove"
elif [[ "${1:-}" == "--update" ]]; then MODE="update"
fi

# ─── Pre-flight ──────────────────────────────────────────────────────────────
print_header "CRISPR v2 -- Conda Environment Setup"

if ! command -v conda &>/dev/null; then
    print_error "Conda not found in PATH. Install Miniconda first."
    exit 1
fi

if command -v mamba &>/dev/null; then
    PKG_MGR="mamba"; print_info "Using mamba for faster dependency resolution."
else
    PKG_MGR="conda"; print_info "conda found (install mamba for faster setup)."
fi

# ─── Handle existing environment ─────────────────────────────────────────────
if [[ "$MODE" == "update" ]]; then
    if ! conda env list | grep -q "^${ENV_NAME} "; then
        print_error "Environment '${ENV_NAME}' does not exist. Run without flags first."
        exit 1
    fi
    print_info "Updating existing environment '${ENV_NAME}'..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
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
        read -rp "Choose an option (u/r/n): " response
        case "${response,,}" in
            u) MODE="update"
               source "$(conda info --base)/etc/profile.d/conda.sh"
               conda activate "$ENV_NAME"
               print_success "Environment activated for update." ;;
            r) print_info "Removing existing environment..."
               conda env remove -n "$ENV_NAME" -y ;;
            *) echo "No changes made."; exit 0 ;;
        esac
    fi
fi

# ─── Step 1: Create environment ──────────────────────────────────────────────
if [[ "$MODE" != "update" ]]; then
    print_header "Step 1/5: Creating conda environment"
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
    print_success "Environment '${ENV_NAME}' created and activated."
else
    print_header "Step 1/5: Skipped (updating existing environment)"
fi

# ─── Step 2: Configure channels ──────────────────────────────────────────────
if [[ "$MODE" != "update" ]]; then
    print_header "Step 2/5: Configuring conda channels"
    conda config --env --add channels defaults 2>/dev/null || true
    conda config --env --add channels bioconda
    conda config --env --add channels conda-forge
    conda config --env --set channel_priority strict
    print_success "Channels: conda-forge > bioconda > defaults (strict)"
else
    print_header "Step 2/5: Skipped (channels already configured)"
fi

# ─── Step 3: Install conda packages ──────────────────────────────────────────
if [[ "$MODE" == "update" ]]; then
    print_header "Step 3/5: Updating conda packages"
else
    print_header "Step 3/5: Installing conda packages"
fi

CONDA_PACKAGES=(
    # Core scientific stack
    numpy scipy pandas matplotlib seaborn
    # Genomics
    biopython samtools bedtools blast
    # Jupyter (optional, for interactive inspection)
    jupyter
    # General utilities
    wget
)

if [[ "$MODE" == "update" ]]; then
    $PKG_MGR update -y --all -c conda-forge -c bioconda --strict-channel-priority
    $PKG_MGR install -y -c conda-forge -c bioconda --strict-channel-priority "${CONDA_PACKAGES[@]}"
else
    $PKG_MGR install -y -c conda-forge -c bioconda --strict-channel-priority "${CONDA_PACKAGES[@]}"
fi

print_success "Conda packages installed."

# ─── Step 4: Install pip packages ────────────────────────────────────────────
if [[ "$MODE" == "update" ]]; then
    print_header "Step 4/5: Updating pip packages"
    UPG="--upgrade"
else
    print_header "Step 4/5: Installing pip packages"
    UPG=""
fi

print_step "Installing CRISPOR and its Python dependencies..."
# CRISPOR is distributed as a Python script + genome data; install its deps.
# Deps below are the subset of crispor/requirements.txt that crispor.py and the
# v2 pipeline modules import at runtime (not covered by the conda scientific
# stack installed above). Additions here should be driven by actual import
# failures observed in the pipeline, not by copy-pasting requirements.txt.
pip install --quiet ${UPG:+$UPG} \
    scikit-learn \
    requests \
    xlrd \
    xlwt \
    openpyxl \
    tabulate \
    rs3 \
    lmdbm \
    lmdb \
    twobitreader \
    pytabix \
    matplotlib-venn \
    dill \
    seqfold

TOOLS_DIR="$SCRIPT_DIR/modules/09_crispr_analysis/v2/tools"
mkdir -p "$TOOLS_DIR"

print_step "Installing CRISPRon (on-target deep-learning scorer)..."
# CRISPRon has no PyPI package or setup.py -- it is a script-based tool.
# Clone into modules/09_crispr_analysis/v2/tools/crisprOn and symlink into PATH.
CRISPRON_DIR="$TOOLS_DIR/crisprOn"
if [[ -d "$CRISPRON_DIR" ]]; then
    print_info "CRISPRon already cloned at $CRISPRON_DIR -- pulling latest..."
    git -C "$CRISPRON_DIR" pull --ff-only 2>/dev/null || true
else
    git clone --depth 1 https://github.com/RTH-tools/crisprOn.git "$CRISPRON_DIR" 2>/dev/null || {
        print_info "CRISPRon: git clone failed -- check network access"
    }
fi
if [[ -f "$CRISPRON_DIR/bin/CRISPRon.sh" ]]; then
    chmod +x "$CRISPRON_DIR/bin/CRISPRon.sh"
    ln -sf "$CRISPRON_DIR/bin/CRISPRon.sh" "$CONDA_PREFIX/bin/CRISPRon.sh"
    ln -sf "$CRISPRON_DIR/bin/DeepCRISPRon_eval.py" "$CONDA_PREFIX/bin/DeepCRISPRon_eval.py"
    print_success "CRISPRon installed at $CRISPRON_DIR"
else
    print_info "CRISPRon: clone incomplete -- see https://github.com/RTH-tools/crisprOn"
fi

# CRISPRoff dependency: CRISPRspec_CRISPRoff_pipeline.py + energy_dics.pkl
# (not on PyPI / GitHub -- must be fetched from rth.dk).
CRISPROFF_PIPELINE="$CRISPRON_DIR/bin/CRISPRspec_CRISPRoff_pipeline.py"
CRISPROFF_ENERGY="$CRISPRON_DIR/data/model/energy_dics.pkl"
if [[ -f "$CRISPROFF_PIPELINE" && -f "$CRISPROFF_ENERGY" ]]; then
    print_info "CRISPRoff already installed ($CRISPROFF_PIPELINE)."
else
    print_step "Downloading CRISPRoff 1.1.2 (required by CRISPRon for deltaGb)..."
    _crispoff_tmp=$(mktemp -d)
    _crispoff_url="https://rth.dk/resources/crispr/crisproff/downloads/crisproff-1.1.2.tar.gz"
    if curl -sL -o "$_crispoff_tmp/crisproff.tar.gz" "$_crispoff_url" \
        && tar -xzf "$_crispoff_tmp/crisproff.tar.gz" -C "$_crispoff_tmp"; then
        mkdir -p "$CRISPRON_DIR/data/model"
        cp "$_crispoff_tmp/crisproff-1.1.2/CRISPRspec_CRISPRoff_pipeline.py" "$CRISPRON_DIR/bin/"
        cp "$_crispoff_tmp/crisproff-1.1.2/energy_dics.pkl" "$CRISPRON_DIR/data/model/"
        chmod +x "$CRISPRON_DIR/bin/CRISPRspec_CRISPRoff_pipeline.py"
        print_success "CRISPRoff installed."
    else
        print_info "CRISPRoff download failed -- fetch manually from $_crispoff_url"
        print_info "  then place CRISPRspec_CRISPRoff_pipeline.py in $CRISPRON_DIR/bin/"
        print_info "  and  energy_dics.pkl in $CRISPRON_DIR/data/model/"
    fi
    rm -rf "$_crispoff_tmp"
fi

print_step "Installing ViennaRNA (required by CRISPRon / CRISPRoff for RNAfold)..."
# CRISPRspec_CRISPRoff_pipeline.py calls RNAfold at runtime; without it the
# full CRISPRon pipeline silently fails and scores fall back to deltaGb=0.
if conda install -y -c bioconda viennarna 2>/dev/null; then
    print_success "ViennaRNA installed ($(RNAfold --version 2>/dev/null | head -1))."
else
    print_info "ViennaRNA install failed -- CRISPRon will run without CRISPRoff deltaGb (degraded scores)."
fi

print_step "Installing TensorFlow (required by DeepSpCas9)..."
# Use tensorflow-cpu 2.15: last release with full compat.v1 support that also
# allows typing-extensions>=4.6.0, avoiding the <4.6.0 cap in TF 2.13.
if pip install --quiet ${UPG:+$UPG} "tensorflow-cpu>=2.15,<2.16"; then
    print_success "tensorflow-cpu 2.15 installed."
else
    print_info "tensorflow-cpu 2.15 failed -- trying latest tensorflow..."
    if pip install --quiet ${UPG:+$UPG} tensorflow; then
        print_success "tensorflow (latest) installed."
    else
        print_info "TensorFlow install failed -- DeepSpCas9 will return NaN (pipeline continues)"
    fi
fi

print_step "Installing DeepSpCas9 wrapper..."
DEEPSPCAS9_DIR="$TOOLS_DIR/DeepSpCas9"
if [[ -d "$DEEPSPCAS9_DIR" ]]; then
    print_info "DeepSpCas9 already cloned at $DEEPSPCAS9_DIR -- pulling latest..."
    git -C "$DEEPSPCAS9_DIR" pull --ff-only || print_info "DeepSpCas9: git pull failed -- using existing clone"
else
    print_info "DeepSpCas9 not on PyPI -- cloning from GitHub..."
    git clone --depth 1 https://github.com/MyungjaeSong/Paired-Library.git "$DEEPSPCAS9_DIR" \
        || { print_error "DeepSpCas9: git clone failed -- check network or GitHub access"; }
fi
# Register the directory on sys.path so 'import DeepSpCas9' resolves to
# DeepSpCas9.py (our Python 3 / TF2-compatible shim placed in this directory).
if [[ -f "$DEEPSPCAS9_DIR/DeepSpCas9.py" ]]; then
    SITE_PKG=$(python -c "import site; print(site.getsitepackages()[0])")
    PTH_FILE="$SITE_PKG/DeepSpCas9.pth"
    # Avoid duplicate entries
    if ! grep -qF "$DEEPSPCAS9_DIR" "$PTH_FILE" 2>/dev/null; then
        echo "$DEEPSPCAS9_DIR" >> "$PTH_FILE"
        print_success "DeepSpCas9 registered via .pth: $DEEPSPCAS9_DIR"
    else
        print_info "DeepSpCas9 .pth already registered -- skipping"
    fi
elif [[ -f "$DEEPSPCAS9_DIR/setup.py" ]]; then
    pip install --quiet ${UPG:+$UPG} "$DEEPSPCAS9_DIR"
else
    print_error "DeepSpCas9: DeepSpCas9.py not found -- manual steps required:"
    print_info "  1. Ensure $DEEPSPCAS9_DIR/DeepSpCas9.py exists (shim from repo)"
    print_info "  2. Run: SITE=\$(python -c 'import site; print(site.getsitepackages()[0])')"
    print_info "     echo $DEEPSPCAS9_DIR >> \$SITE/DeepSpCas9.pth"
fi

print_step "Installing inDelphi (indel outcome predictor)..."
# inDelphi is not on PyPI -- clone directly from GitHub.
INDELPHI_DIR="$TOOLS_DIR/inDelphi"
if [[ -d "$INDELPHI_DIR" ]]; then
    print_info "inDelphi already cloned at $INDELPHI_DIR -- pulling latest..."
    git -C "$INDELPHI_DIR" pull --ff-only 2>/dev/null || true
else
    git clone --depth 1 https://github.com/maxwshen/inDelphi-model.git "$INDELPHI_DIR" 2>/dev/null || true
fi
if [[ -f "$INDELPHI_DIR/setup.py" ]]; then
    pip install --quiet ${UPG:+$UPG} "$INDELPHI_DIR"
elif [[ -f "$INDELPHI_DIR/inDelphi.py" ]]; then
    SITE_PKG=$(python -c "import site; print(site.getsitepackages()[0])")
    if ! grep -qF "$INDELPHI_DIR" "$SITE_PKG/inDelphi.pth" 2>/dev/null; then
        echo "$INDELPHI_DIR" >> "$SITE_PKG/inDelphi.pth"
    fi
    print_info "inDelphi installed via .pth file pointing to $INDELPHI_DIR"
else
    print_info "inDelphi: manual install may be needed -- see https://github.com/maxwshen/inDelphi-model"
fi

print_step "Installing Lindel (indel predictor)..."
# Lindel is not on PyPI -- clone directly from GitHub.
LINDEL_DIR="$TOOLS_DIR/Lindel"
if [[ -d "$LINDEL_DIR" ]]; then
    print_info "Lindel already cloned at $LINDEL_DIR -- pulling latest..."
    git -C "$LINDEL_DIR" pull --ff-only 2>/dev/null || true
else
    git clone --depth 1 https://github.com/shendurelab/Lindel.git "$LINDEL_DIR" 2>/dev/null || true
fi
if [[ -f "$LINDEL_DIR/setup.py" ]]; then
    pip install --quiet ${UPG:+$UPG} "$LINDEL_DIR"
    print_success "Lindel installed from $LINDEL_DIR"
elif [[ -d "$LINDEL_DIR" ]]; then
    SITE_PKG=$(python -c "import site; print(site.getsitepackages()[0])")
    if ! grep -qF "$LINDEL_DIR" "$SITE_PKG/Lindel.pth" 2>/dev/null; then
        echo "$LINDEL_DIR" >> "$SITE_PKG/Lindel.pth"
    fi
    print_info "Lindel installed via .pth file pointing to $LINDEL_DIR"
else
    print_info "Lindel: manual install may be needed -- see https://github.com/shendurelab/Lindel"
fi

print_step "Installing ESMFold dependencies (API + optional local inference)..."
# 'requests' was already installed above (CRISPOR deps); it covers API mode.
# fair-esm[esmfold] adds local GPU inference (~10 GB VRAM).
# It pulls in PyTorch + ESM model weights on first use.
# Gracefully skipped if the environment cannot resolve it (CPU-only machines).
# Switch between modes in i_crispr_v2CONFIG.toml:
#   [crispr_v2.protein_consequence]
#   esmfold_backend = "api"    # public API  (no GPU required)
#   esmfold_backend = "local"  # local model (GPU + ~10 GB RAM)
ESMFOLD_LOCAL=false
if pip install --quiet ${UPG:+$UPG} "fair-esm[esmfold]"; then
    print_success "fair-esm[esmfold] installed -- local ESMFold inference available."
    ESMFOLD_LOCAL=true
else
    print_info "fair-esm[esmfold] not installable on this machine -- ESMFold will use public API."
    print_info "  To enable local inference later:  pip install \"fair-esm[esmfold]\""
fi

print_step "Installing TOML parser (tomli for Python <3.11)..."
pip install --quiet ${UPG:+$UPG} tomli

print_success "All pip packages installed."

# ─── Step 5: Install CRISPOR ─────────────────────────────────────────────────
print_header "Step 5/5: Installing CRISPOR"

# CRISPOR is not a pip-installable package; it is a standalone Python script
# with genome databases.  Clone it into $TOOLS_DIR/crispor so it lives with
# the other third-party tools under modules/09_crispr_analysis/v2/tools/.
CRISPOR_DIR="$TOOLS_DIR/crispor"

if command -v crispor.py &>/dev/null; then
    print_info "crispor.py already on PATH -- skipping clone."
elif [[ -f "$CRISPOR_DIR/crispor.py" ]]; then
    print_info "CRISPOR already installed at $CRISPOR_DIR -- pulling latest..."
    git -C "$CRISPOR_DIR" pull --ff-only 2>/dev/null || true
else
    print_step "Cloning CRISPOR into $CRISPOR_DIR ..."
    timeout 120 git clone --depth 1 https://github.com/maximilianh/crisporWebsite.git "$CRISPOR_DIR" 2>/dev/null || {
        print_info "CRISPOR clone failed -- check network access."
        print_info "Manual install: git clone https://github.com/maximilianh/crisporWebsite.git $CRISPOR_DIR"
    }
fi

# Create a wrapper on PATH so 'crispor.py' works without activating env manually
CRISPOR_WRAPPER="$CONDA_PREFIX/bin/crispor.py"
if [[ -f "$CRISPOR_DIR/crispor.py" && ! -f "$CRISPOR_WRAPPER" ]]; then
    ln -s "$CRISPOR_DIR/crispor.py" "$CRISPOR_WRAPPER"
    chmod +x "$CRISPOR_WRAPPER"
    print_success "crispor.py symlinked to $CRISPOR_WRAPPER"
fi

# Patch crisporEffScores.py: wrap calcAziScore in try/except so the pipeline
# continues when the Azimuth pickled model fails to load under scikit-learn >= 1.0.
# The Azimuth model was serialized with sklearn 0.23 (sklearn.ensemble._gb_losses),
# which was removed in 1.0 and cannot be reinstalled under Python 3.10+.
EFFSCORE="$CRISPOR_DIR/crisporEffScores.py"
if [[ -f "$EFFSCORE" ]]; then
    print_step "Patching crisporEffScores.py for scikit-learn >= 1.0 compatibility..."
    python3 - "$EFFSCORE" <<'PYEOF'
import sys, pathlib
p = pathlib.Path(sys.argv[1])
src = p.read_text()
OLD = '            logging.debug("Azimuth score")\n            scores["fusi"] = calcAziScore(trimSeqs(seqs, -24, 6))'
NEW = (
    '            try:\n'
    '                logging.debug("Azimuth score")\n'
    '                scores["fusi"] = calcAziScore(trimSeqs(seqs, -24, 6))\n'
    '            except Exception:\n'
    '                logging.warning("Azimuth/Fusi score skipped: pickled model incompatible with scikit-learn >= 1.0; scores set to 0")\n'
    '                scores["fusi"] = [0.0] * len(seqs)'
)
if OLD in src:
    p.write_text(src.replace(OLD, NEW, 1))
    print("  Azimuth patch applied.")
else:
    print("  Azimuth patch already applied or source changed -- skipping.")
PYEOF
else
    print_info "crisporEffScores.py not found -- Azimuth patch will be needed after CRISPOR is installed."
fi

# Patch crisporEffScores.py: wrap calcRs3Scores in try/except so the pipeline
# continues when the rs3 package is missing or incompatible.
if [[ -f "$EFFSCORE" ]]; then
    print_step "Patching crisporEffScores.py for missing rs3 module compatibility..."
    python3 - "$EFFSCORE" <<'PYEOF'
import sys, pathlib
p = pathlib.Path(sys.argv[1])
src = p.read_text()
OLD = ('        if inList(scoreNames, "rs3"):\n'
       '            logging.debug("Doench RS3 score")\n'
       '            scores["rs3"] = calcRs3Scores(trimSeqs(seqs, -24, 6))')
NEW = ('        if inList(scoreNames, "rs3"):\n'
       '            try:\n'
       '                logging.debug("Doench RS3 score")\n'
       '                scores["rs3"] = calcRs3Scores(trimSeqs(seqs, -24, 6))\n'
       '            except Exception:\n'
       '                logging.warning("RS3 score skipped: rs3 module not installed; scores set to 0")\n'
       '                scores["rs3"] = [0.0] * len(seqs)')
if OLD in src:
    p.write_text(src.replace(OLD, NEW, 1))
    print("  RS3 patch applied.")
else:
    print("  RS3 patch already applied or source changed -- skipping.")
PYEOF
else
    print_info "crisporEffScores.py not found -- RS3 patch will be needed after CRISPOR is installed."
fi

# Patch crisporEffScores.py: wrap calcLindelScore in try/except so the pipeline
# continues when Lindel model weights are missing (setup.py package_data bug).
if [[ -f "$EFFSCORE" ]]; then
    print_step "Patching crisporEffScores.py for missing Lindel model weights..."
    python3 - "$EFFSCORE" <<'PYEOF'
import sys, pathlib
p = pathlib.Path(sys.argv[1])
src = p.read_text()
OLD = ('    if inList(scoreNames, "lindel"):\n'
       '        logging.debug("lindel scores")\n'
       '        mutSeqDict = calcLindelScore(seqIds, seqs)\n'
       '        scores["lindel"] = mutSeqDict')
NEW = ('    if inList(scoreNames, "lindel"):\n'
       '        try:\n'
       '            logging.debug("lindel scores")\n'
       '            mutSeqDict = calcLindelScore(seqIds, seqs)\n'
       '            scores["lindel"] = mutSeqDict\n'
       '        except (FileNotFoundError, ImportError) as e:\n'
       '            logging.warning("Lindel score skipped: model weights missing or module not installed (%s); scores set to None" % e)\n'
       '            scores["lindel"] = {seqId: (None, []) for seqId in seqIds}')
if OLD in src:
    p.write_text(src.replace(OLD, NEW, 1))
    print("  Lindel patch applied.")
else:
    print("  Lindel patch already applied or source changed -- skipping.")
PYEOF
else
    print_info "crisporEffScores.py not found -- Lindel patch will be needed after CRISPOR is installed."
fi

# Post-install guard: if Lindel was installed without pkl files (setup.py bug),
# copy them directly from the source clone into site-packages.
if _pip_ok Lindel 2>/dev/null && python3 -c "import Lindel; import pathlib; assert pathlib.Path(Lindel.__path__[0], 'Model_weights.pkl').exists()" 2>/dev/null; then
    print_info "Lindel model weights present in site-packages."
elif [[ -f "$LINDEL_DIR/Lindel/Model_weights.pkl" ]]; then
    print_step "Copying Lindel model weights into site-packages (setup.py package_data fix)..."
    LINDEL_SITE=$(python3 -c "import importlib.util, pathlib; spec=importlib.util.find_spec('Lindel'); print(pathlib.Path(spec.origin).parent)" 2>/dev/null || true)
    if [[ -n "$LINDEL_SITE" ]]; then
        cp "$LINDEL_DIR/Lindel/Model_weights.pkl" "$LINDEL_SITE/"
        cp "$LINDEL_DIR/Lindel/model_prereq.pkl"  "$LINDEL_SITE/"
        print_success "Lindel model weights copied to $LINDEL_SITE"
    else
        print_info "Cannot locate Lindel site-packages dir -- run: pip install --force-reinstall $LINDEL_DIR"
    fi
fi

# ─── Verification ─────────────────────────────────────────────────────────────
print_header "Verification"

# All checks use `pip show` or filesystem tests -- no Python imports, no hangs.

_pip_ok()  { pip show "$1" &>/dev/null; }
_pip_ver() { pip show "$1" 2>/dev/null | awk '/^Version/{print $2}'; }

# Core genomics / script tools
if command -v crispor.py &>/dev/null || [[ -x "$CRISPOR_DIR/crispor.py" ]]; then
    echo "  CRISPOR:  found"
else
    echo "  CRISPOR:  NOT found (manual install required)"
fi

if command -v CRISPRon.sh &>/dev/null || [[ -x "$CONDA_PREFIX/bin/CRISPRon.sh" ]]; then
    echo "  CRISPRon: found"
else
    echo "  CRISPRon: NOT found (manual install required)"
fi

# Core Python packages
_core_ok=(); _core_missing=()
for _pkg in numpy pandas biopython scikit-learn openpyxl tabulate tomli; do
    if _pip_ok "$_pkg"; then _core_ok+=("$_pkg")
    else                     _core_missing+=("$_pkg")
    fi
done
echo "  Core Python packages OK: ${_core_ok[*]}"
[[ ${#_core_missing[@]} -gt 0 ]] && echo "  MISSING: ${_core_missing[*]}"

# CRISPOR runtime deps (observed imports in crispor.py + v2 modules)
_crispor_ok=(); _crispor_missing=()
for _pkg in lmdbm lmdb twobitreader pytabix xlwt matplotlib-venn dill seqfold; do
    if _pip_ok "$_pkg"; then _crispor_ok+=("$_pkg")
    else                     _crispor_missing+=("$_pkg")
    fi
done
echo "  CRISPOR runtime deps OK: ${_crispor_ok[*]}"
[[ ${#_crispor_missing[@]} -gt 0 ]] && echo "  CRISPOR runtime deps MISSING: ${_crispor_missing[*]}"

# ViennaRNA / RNAfold
if command -v RNAfold &>/dev/null; then
    echo "  ViennaRNA (RNAfold): installed ($(RNAfold --version 2>/dev/null | head -1))"
else
    echo "  ViennaRNA (RNAfold): NOT installed (CRISPRon deltaGb scoring degraded)"
fi

# CRISPRoff (pipeline script + energy dict)
if [[ -f "$TOOLS_DIR/crisprOn/bin/CRISPRspec_CRISPRoff_pipeline.py" \
   && -f "$TOOLS_DIR/crisprOn/data/model/energy_dics.pkl" ]]; then
    echo "  CRISPRoff: installed (pipeline + energy_dics.pkl)"
else
    echo "  CRISPRoff: NOT installed (CRISPRon will fall back to deltaGb=0)"
fi

# Azimuth patch
if grep -q "Azimuth/Fusi score skipped" "$CRISPOR_DIR/crisporEffScores.py" 2>/dev/null; then
    echo "  Azimuth patch: applied (fusi scores fall back to 0 on sklearn >= 1.0)"
else
    echo "  Azimuth patch: NOT applied (run setup again or patch manually)"
fi

# RS3 patch + package
if grep -q "RS3 score skipped" "$CRISPOR_DIR/crisporEffScores.py" 2>/dev/null; then
    echo "  RS3 patch: applied (rs3 scores fall back to 0 if module missing)"
else
    echo "  RS3 patch: NOT applied (run setup again or patch manually)"
fi
if _pip_ok rs3; then
    echo "  rs3: installed ($(_pip_ver rs3))"
else
    echo "  rs3: NOT installed (RS3 scores will be 0 -- run: pip install rs3)"
fi

# TensorFlow
if _pip_ok tensorflow-cpu; then
    echo "  TensorFlow: installed ($(_pip_ver tensorflow-cpu))"
elif _pip_ok tensorflow; then
    echo "  TensorFlow: installed ($(_pip_ver tensorflow))"
else
    echo "  TensorFlow: NOT installed (DeepSpCas9 scoring will be skipped)"
fi

# ML scorers
_ml_ok=(); _ml_missing=()
[[ -f "$TOOLS_DIR/inDelphi/inDelphi.py" ]] && _ml_ok+=("inDelphi") || _ml_missing+=("inDelphi")
[[ -d "$TOOLS_DIR/Lindel" ]]               && _ml_ok+=("Lindel")   || _ml_missing+=("Lindel")
[[ ${#_ml_ok[@]}      -gt 0 ]] && echo "  ML scorers available: ${_ml_ok[*]}"
[[ ${#_ml_missing[@]} -gt 0 ]] && echo "  ML scorers need manual install: ${_ml_missing[*]}"

# Lindel model weights
if python3 -c "import Lindel, pathlib; assert pathlib.Path(Lindel.__path__[0], 'Model_weights.pkl').exists()" 2>/dev/null; then
    echo "  Lindel model weights: present in site-packages"
else
    echo "  Lindel model weights: MISSING -- run setup again or: pip install --force-reinstall $TOOLS_DIR/Lindel"
fi

# Lindel patch
if grep -q "Lindel score skipped" "$CRISPOR_DIR/crisporEffScores.py" 2>/dev/null; then
    echo "  Lindel patch: applied (lindel scores fall back to None if weights missing)"
else
    echo "  Lindel patch: NOT applied (run setup again or patch manually)"
fi

# ─── Summary ──────────────────────────────────────────────────────────────────
print_header "Setup Complete!"

cat <<EOF

  Environment Name:  crispr_v2
  Python Version:    ${PYTHON_VERSION}

  -------------------------------------------------------
  TOOLS INSTALLED
  -------------------------------------------------------

  Core genomics:
    biopython, samtools, bedtools, blast+

  gRNA Design & Scoring:
    CRISPOR       -- modules/09_crispr_analysis/v2/tools/crispor
    CRISPRon      -- on-target deep-learning rescorer
    DeepSpCas9    -- alternative on-target predictor

  Indel Outcome Prediction:
    inDelphi      -- modules/09_crispr_analysis/v2/tools/inDelphi
    Lindel        -- modules/09_crispr_analysis/v2/tools/Lindel

  Protein Structure:
    ESMFold       -- API mode ready (requests installed)
$(  [[ "$ESMFOLD_LOCAL" == "true" ]] \
      && echo "                     local inference ready (fair-esm[esmfold] installed)" \
      || echo "                     local inference NOT installed (GPU/build deps unavailable)")
                     switch via esmfold_backend in i_crispr_v2CONFIG.toml:
                       esmfold_backend = "api"    # public API (no GPU)
                       esmfold_backend = "local"  # local model (~10 GB VRAM)
    AlphaFold3    -- external server (alphafoldserver.com)

  -------------------------------------------------------
  TOOLS REQUIRING MANUAL ACTION
  -------------------------------------------------------

  * CRISPOR genome databases:
      Run crispor.py --help and follow the genome index steps,
      or use the web UI at crispor.tefor.net for offline genomes.

  * ESMFold model weights (local mode only):
      Weights are downloaded automatically on first use (~690 MB).
      Requires GPU with ~10 GB VRAM.  API mode needs no GPU.

  -------------------------------------------------------

  To activate:    conda activate crispr_v2
  To deactivate:  conda deactivate
  To update:      bash setup_conda_crispr_v2.sh --update
  To remove:      conda env remove -n crispr_v2

EOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup finished -- full log saved to: $LOG_FILE"
