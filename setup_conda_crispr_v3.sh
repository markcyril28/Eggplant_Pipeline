#!/bin/bash
# ============================================================================
# CRISPR v3 Conda Environment Setup (PLANT-ONLY)
# ============================================================================
# Creates a plant-focused conda environment for the v3 CRISPR KO prediction
# pipeline (i_crispr_v3_pipeline.sh). Unlike setup_conda_crispr_v2.sh, this
# script DOES NOT install mammalian-trained scorer tools (CRISPOR, CRISPRon,
# DeepSpCas9, RS3, TensorFlow, ViennaRNA, CRISPRoff) because v3 does not run
# them. The footprint is therefore much smaller (~2 GB vs ~10 GB for v2).
#
#  Stages covered:
#   [1]  Guide source      -- CRISPR-P v2.0 raw CSVs (external web tool; no
#                             conda dep -- user supplies the CSVs manually).
#   [1b] Plant sgRNA filter -- pure Python (numpy + stdlib).
#   [2]  Plant rescorers    -- OPTIONAL; DeepCRISPR / CRISPR-Local clones
#                             enabled with --with-plant-scorers. Dispatch
#                             code in 02_rescore_ontarget.py must be extended
#                             separately before these become functional.
#   [3]  Off-target curation -- pandas + stdlib.
#   [4]  Indel outcome       -- inDelphi (legacy sklearn 0.20 env) + Lindel.
#   [5]  Mutant transcripts  -- Biopython.
#   [6]  Protein consequence -- Biopython (domain-hit + structure_flag).
#   [6b] Local ESMFold fold  -- SEPARATE 'esmfold' conda env (PyTorch 2.8 +
#                                cu128 + openfold). Opt in with
#                                --with-esmfold-local (WSL + NVIDIA GPU only).
#   [7]  Plant NMD           -- pure Python.
#   [8]  Composite ranking   -- pandas, matplotlib, optional openpyxl.
#   [9]  Comparison scatter  -- matplotlib.
#
# Usage:
#   bash setup_conda_crispr_v3.sh                      # Create env (default)
#   bash setup_conda_crispr_v3.sh --remove             # Remove and recreate
#   bash setup_conda_crispr_v3.sh --update             # Update existing env
#   bash setup_conda_crispr_v3.sh --with-plant-scorers # Clone DeepCRISPR +
#                                                       CRISPR-Local (pending
#                                                       dispatch code)
#   bash setup_conda_crispr_v3.sh --with-esmfold-local # Provision the separate
#                                                       'esmfold' conda env for
#                                                       stage [6b] local fold
#                                                       (WSL + NVIDIA GPU only)
#   bash setup_conda_crispr_v3.sh --use-v3-env          # Create / target a
#                                                       separate crispr_v3 env
#                                                       instead of the default
#                                                       shared crispr_v2 env
#                                                       (reserved for future
#                                                       v3-only divergence)
#
# DEFAULT ENV: Despite the filename, this script installs into the shared
# `crispr_v2` / `crispr_v2_indelphi` envs by default — v3 uses the same
# v2 Python modules, so env sharing is the simpler and currently-supported
# path. The `--use-v3-env` flag is wired and ready for when v3 needs its
# own dependency footprint.
#
# NOTE: Stage [6b] local ESMFold uses a SEPARATE conda env named `esmfold`
# (PyTorch 2.8 + cu128 + openfold). Opt in with --with-esmfold-local.
# The build applies the WSL-specific CUDA patches and openfold source
# patches documented in memories/repo/esmfold_local_env_setup.md. If the
# env is missing, stage [6b] skips folding and the rest of the pipeline
# continues.
# ============================================================================

set -euo pipefail

# Default env names: share the v2 envs. v3 currently reuses v2 modules
# verbatim, so a separate env is unnecessary. --use-v3-env opt-in below
# switches these to crispr_v3 / crispr_v3_indelphi for future divergence.
ENV_NAME="crispr_v2"
ENV_INDELPHI="crispr_v2_indelphi"
PYTHON_VERSION="3.10"

WITH_PLANT_SCORERS=false
WITH_ESMFOLD_LOCAL=false
USE_V3_ENV=false
MODE="create"

# Pinned versions for the separate 'esmfold' env (stage [6b] local fold).
# These are the exact versions that have been verified to build under
# WSL + RTX 5050 (Blackwell, sm_120, CUDA 12.8). Changing any of them
# may reintroduce the deepspeed/pytorch-lightning/openfold build issues
# documented in memories/repo/esmfold_local_env_setup.md.
ESMFOLD_ENV_NAME="esmfold"
ESMFOLD_PYTHON_VERSION="3.10"
ESMFOLD_TORCH_VERSION="2.8.0"
ESMFOLD_CUDA_INDEX="https://download.pytorch.org/whl/cu128"
ESMFOLD_CUDA_CHANNEL="nvidia/label/cuda-12.8.0"
ESMFOLD_OPENFOLD_COMMIT="4b41059"
ESMFOLD_TORCH_CUDA_ARCH_LIST="12.0"   # Blackwell sm_120; adjust per GPU

# ─── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --remove)               MODE="remove"; shift ;;
        --update)               MODE="update"; shift ;;
        --with-plant-scorers)   WITH_PLANT_SCORERS=true; shift ;;
        --with-esmfold-local)   WITH_ESMFOLD_LOCAL=true; shift ;;
        --use-v3-env)           USE_V3_ENV=true; ENV_NAME="crispr_v3"
                                ENV_INDELPHI="crispr_v3_indelphi"; shift ;;
        -h|--help)
            # Print the leading comment block (everything from the first line
            # up to — but not including — the first non-comment, non-shebang
            # line), then stop. Avoids the brittle `head -N` cap that used to
            # truncate any new usage section added past the 40/55-line mark.
            awk '
                NR==1 && /^#!/ { next }
                /^#/ { sub(/^# ?/, ""); print; next }
                { exit }
            ' "$0"
            exit 0 ;;
        *) echo "Unknown flag: $1" >&2; exit 1 ;;
    esac
done

# ─── Logging ─────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs/installation_logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/setup_conda_crispr_v3_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee >(sed -u \
    -e $'s/\r/\\\n/g' \
    -e $'s/\x1B\\[[0-9;?]*[a-zA-Z]//g' \
    -e $'s/\x1B[()][A-Z0-9]//g' \
    -e 's/[[:space:]]*\([-\\|/][[:space:]]*\)\{2,\}/ /g' \
    -e '/^[[:space:]]*[-\\|/][[:space:]]*$/d' \
    >> "$LOG_FILE")) 2>&1
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup started -- log: $LOG_FILE"

# ─── Colors ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; NC='\033[0m'

print_header()  { echo; echo -e "${CYAN}============================================================${NC}"; echo -e "${CYAN}  $1${NC}"; echo -e "${CYAN}============================================================${NC}"; }
print_step()    { echo -e "${GREEN}[STEP]${NC} $1"; }
print_info()    { echo -e "${YELLOW}[INFO]${NC} $1"; }
print_error()   { echo -e "${RED}[ERROR]${NC} $1"; }
print_success() { echo -e "${GREEN}[OK]${NC} $1"; }

# ─── Pre-flight ──────────────────────────────────────────────────────────────
print_header "CRISPR v3 -- Plant-Only Conda Environment Setup"
print_info "Main env:     $ENV_NAME"
print_info "Legacy env:   $ENV_INDELPHI  (inDelphi sklearn 0.20 compatibility)"
print_info "Python:       $PYTHON_VERSION"
[[ "$WITH_PLANT_SCORERS" == "true" ]] && print_info "Plant scorers: DeepCRISPR + CRISPR-Local clones ENABLED"
[[ "$WITH_ESMFOLD_LOCAL" == "true" ]] && print_info "ESMFold local: '$ESMFOLD_ENV_NAME' env build ENABLED (WSL + NVIDIA GPU required)"
[[ "$USE_V3_ENV"         == "true" ]] && print_info "Using SEPARATE '$ENV_NAME' env (opt-in; default is shared crispr_v2)"

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
        print_error "Environment '${ENV_NAME}' does not exist. Run without --update first."
        exit 1
    fi
    print_info "Updating existing environment '${ENV_NAME}'..."
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
elif conda env list | grep -q "^${ENV_NAME} "; then
    if [[ "$MODE" == "remove" ]]; then
        print_info "Removing existing environment '${ENV_NAME}'..."
        conda env remove -n "$ENV_NAME" -y
    elif [[ "$USE_V3_ENV" != "true" && "$ENV_NAME" == "crispr_v2" ]]; then
        # Default path: v2 env already exists (created by setup_conda_crispr_v2.sh
        # or a prior v3 run). Auto-update without prompting — the user did not
        # ask for a fresh env, so adding any missing plant deps is the least
        # surprising behaviour.
        print_info "Reusing existing crispr_v2 env (adding any missing plant deps)."
        MODE="update"
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$ENV_NAME"
    else
        print_info "Environment '${ENV_NAME}' already exists."
        # Non-interactive runs (CI, agents, `bash ... < /dev/null`, nohup) must
        # not block on stdin. When no TTY is attached, fall through to update
        # mode so the script finishes without human input. The user can still
        # force a fresh env with --remove, or cancel by passing --update which
        # exits early above when the env is missing.
        if [[ ! -t 0 ]]; then
            print_info "No TTY detected — defaulting to UPDATE mode (use --remove to recreate)."
            MODE="update"
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "$ENV_NAME"
            print_success "Environment activated for update."
        else
            echo
            echo "  [u] Update packages in the existing environment"
            echo "  [r] Remove and recreate from scratch"
            echo "  [n] Cancel (keep as-is)"
            echo
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

# ─── Step 3: Install conda packages (plant-only footprint) ───────────────────
# Notable absences vs v2:
#   - NO tensorflow-cpu      (DeepSpCas9 not used in v3)
#   - NO viennarna / RNAfold (CRISPRon/CRISPRoff not used)
#   - NO blast               (CRISPOR dep; v3 doesn't run CRISPOR)
print_header "Step 3/5: Installing conda packages (plant-only)"
CONDA_PACKAGES=(
    # Core scientific stack
    numpy scipy pandas matplotlib seaborn
    # Genomics
    biopython samtools bedtools
    # Utilities
    wget curl
)

# Narrow-scope update policy: we NEVER run `$PKG_MGR update --all` here
# because the default env is the shared crispr_v2 env (see header). A blast
# update would upgrade every v2 package, risking a silent break of stage
# [2] mammalian rescorers / stage [1] CRISPOR in v2 runs. Instead we only
# `install` the v3 plant-only set: `conda install` bumps those specific
# packages to the latest version satisfying the channel constraints when
# they already exist, and installs them fresh when they don't.
if [[ "$MODE" == "update" ]]; then
    print_info "Update mode: refreshing only the plant-only package list (${#CONDA_PACKAGES[@]} pkgs);"
    print_info "  'update --all' is deliberately avoided to protect the shared crispr_v2 env."
fi
$PKG_MGR install -y -c conda-forge -c bioconda --strict-channel-priority "${CONDA_PACKAGES[@]}"

print_success "Conda packages installed."

# ─── Step 4: Install pip packages (plant-only pipeline deps) ─────────────────
print_header "Step 4/5: Installing pip packages"
# Same narrow-scope policy as Step 3: in update mode we upgrade ONLY the
# explicit plant-only pip packages (via `--upgrade`), never a blanket
# `pip install --upgrade *`. Anything else in the shared env is left alone.
UPG=""
[[ "$MODE" == "update" ]] && UPG="--upgrade"

# Notable absences vs v2:
#   - NO rs3, xlrd, xlwt, matplotlib-venn, dill, seqfold, lmdbm, lmdb,
#     twobitreader, pytabix   (all CRISPOR runtime deps)
pip install --quiet ${UPG} \
    scikit-learn \
    openpyxl \
    tabulate \
    tomli \
    requests

print_success "Core pip packages installed."

# ─── Step 4b: Tool clones (stage [4] indel prediction) ───────────────────────
# inDelphi and Lindel are INDEL OUTCOME predictors (distinct from efficacy
# scorers). v3 retains them with an explicit plant-NHEJ caveat (see
# [crispr_v3.predict_indels] in i_crispr_v3CONFIG.toml). No plant-trained
# indel-outcome model currently exists.
TOOLS_DIR="$SCRIPT_DIR/modules/09_crispr_analysis/v2/tools"
mkdir -p "$TOOLS_DIR"

print_step "Installing inDelphi (indel outcome predictor)..."
# inDelphi is not on PyPI -- clone + .pth registration. The module's pickled
# model weights require scikit-learn 0.18.1 or 0.20.0, which is why a
# separate legacy conda env ($ENV_INDELPHI) is created below.
INDELPHI_DIR="$TOOLS_DIR/inDelphi"
if [[ -d "$INDELPHI_DIR" ]]; then
    print_info "inDelphi already cloned at $INDELPHI_DIR -- pulling latest..."
    git -C "$INDELPHI_DIR" pull --ff-only 2>/dev/null || true
else
    git clone --depth 1 https://github.com/maxwshen/inDelphi-model.git "$INDELPHI_DIR" 2>/dev/null \
        || print_info "inDelphi: git clone failed -- check network."
fi
if [[ -f "$INDELPHI_DIR/inDelphi.py" ]]; then
    SITE_PKG=$(python -c "import site; print(site.getsitepackages()[0])")
    if ! grep -qF "$INDELPHI_DIR" "$SITE_PKG/inDelphi.pth" 2>/dev/null; then
        echo "$INDELPHI_DIR" >> "$SITE_PKG/inDelphi.pth"
    fi
    print_info "inDelphi registered via .pth: $INDELPHI_DIR"
else
    print_info "inDelphi.py not present -- manual install may be required."
fi

print_step "Installing Lindel (indel predictor)..."
LINDEL_DIR="$TOOLS_DIR/Lindel"
if [[ -d "$LINDEL_DIR" ]]; then
    print_info "Lindel already cloned at $LINDEL_DIR -- pulling latest..."
    git -C "$LINDEL_DIR" pull --ff-only 2>/dev/null || true
else
    git clone --depth 1 https://github.com/shendurelab/Lindel.git "$LINDEL_DIR" 2>/dev/null \
        || print_info "Lindel: git clone failed -- check network."
fi
if [[ -f "$LINDEL_DIR/setup.py" ]]; then
    pip install --quiet ${UPG} "$LINDEL_DIR" || print_info "Lindel pip install failed."
    # Post-install guard: Lindel's setup.py package_data sometimes drops the
    # pkl weights from the wheel. Copy them directly from the source clone.
    LINDEL_SITE=$(python3 -c "import importlib.util, pathlib; spec=importlib.util.find_spec('Lindel'); print(pathlib.Path(spec.origin).parent)" 2>/dev/null || true)
    if [[ -n "$LINDEL_SITE" && -f "$LINDEL_DIR/Lindel/Model_weights.pkl" ]]; then
        [[ -f "$LINDEL_SITE/Model_weights.pkl" ]] || cp "$LINDEL_DIR/Lindel/Model_weights.pkl" "$LINDEL_SITE/"
        [[ -f "$LINDEL_SITE/model_prereq.pkl" ]]  || cp "$LINDEL_DIR/Lindel/model_prereq.pkl"  "$LINDEL_SITE/"
        print_info "Lindel model weights verified in $LINDEL_SITE"
    fi
elif [[ -d "$LINDEL_DIR" ]]; then
    SITE_PKG=$(python -c "import site; print(site.getsitepackages()[0])")
    if ! grep -qF "$LINDEL_DIR" "$SITE_PKG/Lindel.pth" 2>/dev/null; then
        echo "$LINDEL_DIR" >> "$SITE_PKG/Lindel.pth"
    fi
    print_info "Lindel registered via .pth file: $LINDEL_DIR"
else
    print_info "Lindel: manual install may be required."
fi

# ─── Step 4c: Plant-trained scorers (OPT-IN via --with-plant-scorers) ────────
# DeepCRISPR and CRISPR-Local are plant-trained on-target scorers configured
# under [crispr_v3.plant_scorer] in i_crispr_v3CONFIG.toml. They require
# dispatch code to be added to 02_rescore_ontarget.py before they become
# functional. Cloning here is purely so the tool source is on disk; the
# pipeline will log "Unknown predictor" until dispatch branches exist.
if [[ "$WITH_PLANT_SCORERS" == "true" ]]; then
    print_step "Installing plant scorers (DeepCRISPR + CRISPR-Local)..."

    # ─── DeepCRISPR (mammalian-trained CNN; commonly used as a plant baseline)
    DEEPCRISPR_DIR="$TOOLS_DIR/DeepCRISPR"
    if [[ -f "$DEEPCRISPR_DIR/run_examples.py" ]]; then
        print_info "DeepCRISPR already installed at $DEEPCRISPR_DIR"
    else
        print_info "Downloading DeepCRISPR ZIP snapshot (git clone fails on some networks)..."
        mkdir -p "$TOOLS_DIR/_downloads"
        if timeout 600 curl -sSL \
                -o "$TOOLS_DIR/_downloads/DeepCRISPR.zip" \
                https://github.com/bm2-lab/DeepCRISPR/archive/refs/heads/master.zip \
                && unzip -q "$TOOLS_DIR/_downloads/DeepCRISPR.zip" -d "$TOOLS_DIR/_downloads/"; then
            rm -rf "$DEEPCRISPR_DIR"
            mv "$TOOLS_DIR/_downloads/DeepCRISPR-master" "$DEEPCRISPR_DIR"
            print_success "DeepCRISPR installed at $DEEPCRISPR_DIR"
        else
            print_info "DeepCRISPR download failed -- see DeepCRISPR/README.md for manual steps."
        fi
    fi
    # Copy the v3 inference helper into the DeepCRISPR folder if the v3 tree
    # includes it (setup-script colocation pattern).
    if [[ -f "$DEEPCRISPR_DIR/deepcrispr_infer.py" ]]; then
        print_info "DeepCRISPR inference helper present: $DEEPCRISPR_DIR/deepcrispr_infer.py"
    fi

    # ─── CRISPR-Local (real GitHub repo, discovered April 2026)
    CRISPR_LOCAL_DIR="$TOOLS_DIR/CRISPR-Local"
    # Accept pre-existing ZIP-extracted folder names too, so a user who
    # unzipped manually doesn't need to rename.
    for _cl_alt in "$TOOLS_DIR/CRISPR-Local" "$TOOLS_DIR/CRISPR-Local-master" \
                   "$TOOLS_DIR/CRISPR-Local.new"; do
        if [[ -f "$_cl_alt/RD-build.pl" || -f "$_cl_alt/PL-search.pl" ]]; then
            CRISPR_LOCAL_DIR="$_cl_alt"
            break
        fi
    done
    if [[ -f "$CRISPR_LOCAL_DIR/RD-build.pl" ]]; then
        print_info "CRISPR-Local already installed at $CRISPR_LOCAL_DIR"
    else
        print_info "Cloning CRISPR-Local from sunjiamin0824/CRISPR-Local..."
        if timeout 300 git -c http.postBuffer=524288000 clone --depth 1 \
                https://github.com/sunjiamin0824/CRISPR-Local.git \
                "$TOOLS_DIR/CRISPR-Local.new" 2>/dev/null; then
            # Prefer the canonical lowercase name, but Windows sometimes
            # blocks the rename — v3 dispatch accepts both.
            if mv "$TOOLS_DIR/CRISPR-Local.new" "$TOOLS_DIR/CRISPR-Local" 2>/dev/null; then
                CRISPR_LOCAL_DIR="$TOOLS_DIR/CRISPR-Local"
            else
                CRISPR_LOCAL_DIR="$TOOLS_DIR/CRISPR-Local.new"
                print_info "CRISPR-Local rename blocked by OS; dispatch accepts the .new name."
            fi
            print_success "CRISPR-Local installed at $CRISPR_LOCAL_DIR"
        else
            print_info "CRISPR-Local git clone failed -- see CRISPR-Local/README.md for manual steps."
        fi
    fi

    # ─── Isolated runtime envs (idempotent: skip if already present)
    if [[ -f "$DEEPCRISPR_DIR/run_examples.py" ]]; then
        if conda env list | grep -q "^crispr_v3_deepcrispr "; then
            print_info "crispr_v3_deepcrispr env already exists -- skipping."
        else
            print_step "Creating crispr_v3_deepcrispr env (Python 3.6 + TF 1.3 + sonnet 1.9)..."
            if conda create -y -n crispr_v3_deepcrispr python=3.6 numpy pandas -y \
                    && conda run -n crispr_v3_deepcrispr pip install \
                        tensorflow==1.3.0 dm-sonnet==1.9 2>/dev/null; then
                print_success "crispr_v3_deepcrispr ready. Set in i_crispr_v3CONFIG.toml:"
                print_info "  [crispr_v3.plant_scorer].deepcrispr_env = \"crispr_v3_deepcrispr\""
            else
                print_info "crispr_v3_deepcrispr env creation failed (TF 1.3 wheels unavailable on"
                print_info "  some newer platforms). Scorer will return NaN until this env exists."
            fi
        fi
    fi
    if [[ -f "$CRISPR_LOCAL_DIR/Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py" ]]; then
        if conda env list | grep -q "^crispr_v3_crispr_local "; then
            print_info "crispr_v3_crispr_local env already exists -- skipping."
        else
            print_step "Creating crispr_v3_crispr_local env (Python 2.7 for Rule Set 2)..."
            if conda create -y -n crispr_v3_crispr_local python=2.7 \
                    scikit-learn=0.17 pandas biopython; then
                print_success "crispr_v3_crispr_local ready. Set in i_crispr_v3CONFIG.toml:"
                print_info "  [crispr_v3.plant_scorer].crispr_local_env = \"crispr_v3_crispr_local\""
            else
                print_info "crispr_v3_crispr_local env creation failed. Scorer will return NaN."
            fi
        fi
    fi
else
    print_info "Plant-scorer install skipped. Use --with-plant-scorers to enable."
fi

# ─── Step 4d: ESMFold local inference — separate 'esmfold' env ───────────────
# Stage [6b] subprocess-invokes a dedicated conda env (PyTorch 2.8 + cu128 +
# openfold). The build sequence below codifies the recipe in
# memories/repo/esmfold_local_env_setup.md. All steps are best-effort: any
# failure logs a warning and leaves the env in whatever state it reached —
# stage [6b] will simply skip folding if the env cannot load ESMFold.
install_esmfold_env() {
    local env="$ESMFOLD_ENV_NAME"

    # Platform gate: openfold's nvcc build expects WSL libcuda.so at the WSL
    # path. Native Windows bash has no NVIDIA driver stack accessible to
    # conda, and macOS has no CUDA. We permit bare Linux + WSL; others skip.
    if [[ -z "${WSL_DISTRO_NAME:-}" ]] && [[ "$(uname -s)" != "Linux" ]]; then
        print_info "ESMFold build requires Linux or WSL (uname=$(uname -s)); skipping."
        return 0
    fi

    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$env"; then
        print_info "ESMFold env '$env' already exists — skipping create."
    else
        print_step "Creating '$env' env (Python $ESMFOLD_PYTHON_VERSION)..."
        conda create -y -n "$env" python="$ESMFOLD_PYTHON_VERSION" \
            || { print_error "Failed to create '$env' env."; return 1; }
    fi

    local env_prefix
    env_prefix=$(conda env list 2>/dev/null | awk -v n="$env" '$1==n {print $NF; exit}')
    if [[ -z "$env_prefix" || ! -d "$env_prefix" ]]; then
        print_error "Could not resolve prefix for '$env' env."
        return 1
    fi

    print_step "Installing CUDA 12.8 toolchain into '$env' (cuda-nvcc + libcu{blas,sparse,solver,rand,fft}-dev + libnvjitlink-dev)..."
    conda install -y -n "$env" -c "$ESMFOLD_CUDA_CHANNEL" \
        cuda-nvcc cuda-cudart-dev \
        libcublas-dev libcusparse-dev libcusolver-dev libcurand-dev libcufft-dev \
        libnvjitlink-dev \
        || print_info "CUDA toolchain install reported errors — openfold compile may fail."

    print_step "Installing PyTorch $ESMFOLD_TORCH_VERSION (+cu128) via pip..."
    conda run -n "$env" pip install --quiet \
        "torch==$ESMFOLD_TORCH_VERSION" --index-url "$ESMFOLD_CUDA_INDEX" \
        || { print_error "PyTorch install failed."; return 1; }

    print_step "Installing pinned fair-esm / deepspeed / pytorch-lightning / setuptools..."
    # openfold imports `pkg_resources` (removed in setuptools 81+), uses
    # `pytorch_lightning.utilities.seed.seed_everything` (removed in PL 2.x),
    # and older deepspeed (<0.13) referenced the removed `torch._six` module.
    conda run -n "$env" pip install --quiet \
        "fair-esm[esmfold]==2.0.0" \
        "deepspeed>=0.18,<0.19" \
        "pytorch-lightning==1.9.5" \
        "setuptools<81" \
        dllogger@git+https://github.com/NVIDIA/dllogger.git \
        || print_info "One or more esmfold deps failed to install."

    # ── Openfold build from source with two setup.py patches ─────────────────
    local build_dir
    build_dir=$(mktemp -d -t openfold_build_XXXXXX)
    print_step "Cloning openfold ($ESMFOLD_OPENFOLD_COMMIT) into $build_dir..."
    if ! git clone --quiet https://github.com/aqlaboratory/openfold.git "$build_dir" 2>/dev/null; then
        print_info "openfold clone failed — skipping build."
        rm -rf "$build_dir"
        return 0
    fi
    ( cd "$build_dir" && git checkout --quiet "$ESMFOLD_OPENFOLD_COMMIT" ) \
        || print_info "openfold commit checkout failed (using default branch)."

    print_step "Patching openfold/setup.py (GPU arch list + C++17)..."
    # Patch 1: replace the hardcoded compute_capabilities set with Blackwell
    # (12.0) + current datacentre (8.0, 9.0) and allow TORCH_CUDA_ARCH_LIST
    # override so the build doesn't invoke get_nvidia_cc() (fails under WSL).
    python3 - "$build_dir/setup.py" <<'PY' || print_info "setup.py patch (arch) failed."
import re, sys, pathlib
p = pathlib.Path(sys.argv[1])
s = p.read_text()
# Replace the static set of (major, minor) compute capabilities with the
# modern triple. Tolerates whitespace variation between the original lines.
s_new = re.sub(
    r"compute_capabilities\s*=\s*set\(\s*\[\s*\(.*?\)\s*\]\s*\)",
    "compute_capabilities = {(8, 0), (9, 0), (12, 0)}",
    s, count=1, flags=re.DOTALL,
)
# Inject TORCH_CUDA_ARCH_LIST override before any get_nvidia_cc() call.
if "TORCH_CUDA_ARCH_LIST" not in s_new:
    s_new = s_new.replace(
        "get_nvidia_cc()",
        "(os.environ.get('TORCH_CUDA_ARCH_LIST') or get_nvidia_cc())",
        1,
    )
    header = "\n".join(s_new.splitlines()[:30])
    if not re.search(r"^\s*import os\b", header, re.MULTILINE):
        s_new = "import os\n" + s_new
p.write_text(s_new)
PY
    # Patch 2: bump C++ standard (PyTorch 2.8 headers need C++17).
    sed -i 's/-std=c++14/-std=c++17/g' "$build_dir/setup.py" 2>/dev/null || true

    print_step "Building openfold (--no-build-isolation; this takes ~10-20 min)..."
    (
        cd "$build_dir"
        CUDA_HOME="$env_prefix" \
        PATH="$env_prefix/bin:$PATH" \
        LD_LIBRARY_PATH="/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}" \
        CPATH="$env_prefix/targets/x86_64-linux/include:${CPATH:-}" \
        TORCH_CUDA_ARCH_LIST="$ESMFOLD_TORCH_CUDA_ARCH_LIST" \
        conda run -n "$env" pip install --no-build-isolation --quiet .
    ) || print_info "openfold build failed — stage [6b] will be disabled until fixed."

    # Runtime shim: deepspeed 0.13+ removed `deepspeed.utils.is_initialized`.
    # Openfold's primitives.py still calls it at import time; replace with
    # a getattr fallback that returns False (we don't use deepspeed anyway).
    print_step "Applying deepspeed.is_initialized shim to openfold/model/primitives.py..."
    local primitives
    primitives=$(conda run -n "$env" python -c \
        "import openfold.model.primitives as m; print(m.__file__)" 2>/dev/null || echo "")
    if [[ -n "$primitives" && -f "$primitives" ]]; then
        sed -i \
            's|deepspeed\.utils\.is_initialized()|getattr(deepspeed.utils, "is_initialized", lambda: False)()|g' \
            "$primitives" \
            || print_info "primitives.py patch (sed) failed."
    else
        print_info "openfold.model.primitives not importable — shim skipped."
    fi

    # Activation hook: sets CUDA_HOME + LD_LIBRARY_PATH + CPATH on env activate
    # so deepspeed imports succeed (deepspeed probes CUDA_HOME at module load).
    print_step "Writing conda activation hook ($env_prefix/etc/conda/activate.d/cuda_env.sh)..."
    mkdir -p "$env_prefix/etc/conda/activate.d"
    cat > "$env_prefix/etc/conda/activate.d/cuda_env.sh" <<EOF
export CUDA_HOME="\$CONDA_PREFIX"
case ":\${LD_LIBRARY_PATH:-}:" in
    *":/usr/lib/wsl/lib:"*) ;;
    *) export LD_LIBRARY_PATH="/usr/lib/wsl/lib:\${LD_LIBRARY_PATH:-}" ;;
esac
case ":\${CPATH:-}:" in
    *":\$CONDA_PREFIX/targets/x86_64-linux/include:"*) ;;
    *) export CPATH="\$CONDA_PREFIX/targets/x86_64-linux/include:\${CPATH:-}" ;;
esac
EOF

    # Pre-download model weights with wget -c to avoid the silent-truncation
    # failure described in the memory file (torch.hub.load_state_dict_from_url
    # renames a partial download to final if the connection drops).
    print_step "Pre-downloading ESMFold weights (~8.5 GB total) to ~/.cache/torch/hub/checkpoints/..."
    local ckpt_dir="$HOME/.cache/torch/hub/checkpoints"
    mkdir -p "$ckpt_dir"
    local base="https://dl.fbaipublicfiles.com/fair-esm/models"
    for f in esmfold_3B_v1.pt esm2_t36_3B_UR50D.pt esm2_t36_3B_UR50D-contact-regression.pt; do
        if [[ ! -s "$ckpt_dir/$f" ]]; then
            wget -c --tries=3 -q --show-progress -O "$ckpt_dir/$f" "$base/$f" \
                || print_info "Weight download failed: $f (stage [6b] will retry via torch.hub)."
        else
            print_info "Weight already present: $f"
        fi
    done

    # Verification: try loading ESMFold. Skip the forward pass (would require
    # an allocated GPU and take too long for a setup script).
    print_step "Verifying esm.pretrained.esmfold_v1() loads..."
    if conda run -n "$env" python -c \
            "import esm; m = esm.pretrained.esmfold_v1(); print('ESMFold load OK')" 2>/dev/null; then
        print_success "ESMFold env '$env' is ready."
    else
        print_info "ESMFold env '$env' built but esmfold_v1() failed to import. "
        print_info "  Check memories/repo/esmfold_local_env_setup.md for the manual fix list."
    fi

    rm -rf "$build_dir"
    return 0
}

if [[ "$WITH_ESMFOLD_LOCAL" == "true" ]]; then
    print_header "Step 4d/5: Provisioning 'esmfold' env (stage [6b] local fold)"
    install_esmfold_env || print_info "ESMFold env provisioning did not complete cleanly."
else
    print_info "ESMFold local inference skipped. Stage [6b] will log a WARN"
    print_info "  and continue without folding. Re-run with --with-esmfold-local"
    print_info "  (WSL + NVIDIA GPU only) to provision the separate 'esmfold' env."
fi

# ─── Step 5: Create legacy inDelphi conda env ────────────────────────────────
# inDelphi's pickled model weights were serialized under scikit-learn 0.18.1
# / 0.20.0; they fail to load under modern sklearn. A separate Python 3.7
# env ($ENV_INDELPHI) isolates the legacy dependency.
print_header "Step 5/5: Creating legacy env for inDelphi ($ENV_INDELPHI)"
if conda env list | grep -q "^${ENV_INDELPHI} "; then
    print_info "Legacy env '${ENV_INDELPHI}' already exists -- skipping."
else
    print_step "Creating ${ENV_INDELPHI} (Python 3.7, scikit-learn=0.20.0)..."
    if conda create -y -n "$ENV_INDELPHI" \
            -c conda-forge -c bioconda --strict-channel-priority \
            python=3.7 scikit-learn=0.20.0 numpy=1.19 pandas=1.1; then
        print_success "Legacy env '${ENV_INDELPHI}' created."
    else
        print_info "Legacy env creation failed -- stage [4] will fall back to Lindel only."
    fi
fi

# ─── Verification ─────────────────────────────────────────────────────────────
print_header "Verification"

_pip_ok()  { pip show "$1" &>/dev/null; }
_pip_ver() { pip show "$1" 2>/dev/null | awk '/^Version/{print $2}'; }

# Core scientific stack
_core_ok=(); _core_missing=()
for _pkg in numpy pandas biopython scikit-learn matplotlib openpyxl tabulate tomli requests; do
    if _pip_ok "$_pkg"; then _core_ok+=("$_pkg")
    else                     _core_missing+=("$_pkg")
    fi
done
echo "  Core pip packages OK:    ${_core_ok[*]}"
[[ ${#_core_missing[@]} -gt 0 ]] && echo "  MISSING: ${_core_missing[*]}"

# Genomics via conda
for _bin in samtools bedtools; do
    if command -v "$_bin" &>/dev/null; then
        echo "  $_bin: $( $_bin --version 2>&1 | head -1 )"
    else
        echo "  $_bin: NOT installed"
    fi
done

# Stage [4] tools
if [[ -f "$TOOLS_DIR/inDelphi/inDelphi.py" ]]; then
    echo "  inDelphi: cloned at $TOOLS_DIR/inDelphi"
else
    echo "  inDelphi: NOT cloned -- stage [4] will log 'unavailable' and skip"
fi
if python3 -c "import Lindel" 2>/dev/null; then
    echo "  Lindel:   importable"
    if python3 -c "import Lindel, pathlib; assert pathlib.Path(Lindel.__path__[0], 'Model_weights.pkl').exists()" 2>/dev/null; then
        echo "  Lindel model weights: present"
    else
        echo "  Lindel model weights: MISSING -- re-run setup or reinstall from $TOOLS_DIR/Lindel"
    fi
else
    echo "  Lindel:   NOT importable -- stage [4] will fall back to NA"
fi

# Plant scorers
if [[ "$WITH_PLANT_SCORERS" == "true" ]]; then
    [[ -d "$TOOLS_DIR/DeepCRISPR" ]]    && echo "  DeepCRISPR clone:   present (dispatch code still needed)" \
                                         || echo "  DeepCRISPR clone:   MISSING"
    [[ -d "$TOOLS_DIR/CRISPR-Local" ]]  && echo "  CRISPR-Local clone: present (dispatch code still needed)" \
                                         || echo "  CRISPR-Local clone: MISSING"
fi

# ESMFold (separate env — provisioned via --with-esmfold-local or manually)
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$ESMFOLD_ENV_NAME"; then
    _esm_ok=$(conda run -n "$ESMFOLD_ENV_NAME" python -c \
        "import esm; esm.pretrained.esmfold_v1(); print('OK')" 2>/dev/null || echo "FAIL")
    if [[ "$_esm_ok" == "OK" ]]; then
        echo "  $ESMFOLD_ENV_NAME env: ready (esmfold_v1 loads)"
    else
        echo "  $ESMFOLD_ENV_NAME env: present but esmfold_v1 failed to load"
        echo "    See memories/repo/esmfold_local_env_setup.md for the fix list."
    fi
else
    echo "  $ESMFOLD_ENV_NAME env: NOT created -- stage [6b] will skip local folding"
    echo "    Re-run with --with-esmfold-local to provision it (WSL + NVIDIA GPU only)."
fi

# Legacy inDelphi env
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_INDELPHI"; then
    _indelphi_sklearn=$(conda run -n "$ENV_INDELPHI" python -c "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo "UNAVAILABLE")
    echo "  $ENV_INDELPHI: present (sklearn=$_indelphi_sklearn)"
    if [[ "$_indelphi_sklearn" != "0.20.0" && "$_indelphi_sklearn" != "0.18.1" ]]; then
        echo "  WARNING: inDelphi legacy env has sklearn $_indelphi_sklearn (expected 0.20.0 or 0.18.1)"
    fi
else
    echo "  $ENV_INDELPHI: NOT created -- inDelphi unavailable; Lindel-only fallback active"
fi

# ─── Summary ──────────────────────────────────────────────────────────────────
print_header "Setup Complete!"

cat <<EOF

  Main environment:    $ENV_NAME
  Legacy env:          $ENV_INDELPHI
  Python:              $PYTHON_VERSION

  -------------------------------------------------------
  PLANT-ONLY TOOLS INSTALLED
  -------------------------------------------------------

  Core genomics:
    biopython, samtools, bedtools

  Stage [4] indel outcome predictors:
    inDelphi   -- $TOOLS_DIR/inDelphi
                  (runs under $ENV_INDELPHI for sklearn 0.20 compat)
    Lindel     -- $TOOLS_DIR/Lindel

  Stage [6b] ESMFold local:
$(  [[ "$WITH_ESMFOLD_LOCAL" == "true" ]] \
      && echo "    Provisioned separate '$ESMFOLD_ENV_NAME' env (PyTorch 2.8 + cu128 + openfold)" \
      || echo "    NOT provisioned -- re-run with --with-esmfold-local (WSL + NVIDIA GPU)")
    Recipe + troubleshooting: memories/repo/esmfold_local_env_setup.md

  -------------------------------------------------------
  TOOLS DELIBERATELY NOT INSTALLED (vs v2)
  -------------------------------------------------------
  CRISPOR, CRISPRon, DeepSpCas9, TensorFlow, ViennaRNA, CRISPRoff, rs3 -- all
  mammalian-scorer infrastructure. v3 does not run these.

  -------------------------------------------------------
  MANUAL ACTIONS REQUIRED BEFORE RUNNING v3
  -------------------------------------------------------

  [1] Provide CRISPR-P v2.0 raw scoring CSVs:
        Run https://crispr.hzau.edu.cn/CRISPR2/ on each target gene FASTA
        and place the output CSVs in the directory resolved by
        [crispr_v3.crispr_p].raw_dir_template in i_crispr_v3CONFIG.toml.

  [2] Annotation format:
        The stage-[5]/[7] parsers accept BOTH GTF ('gene_id "X"'; 'transcript_id "Y"')
        and GFF3 ('ID=X;Parent=Y'; Helixer-style '.N' isoform suffixes are
        stripped to derive gene_id) attribute styles. No conversion step is
        required for the bundled GPE001970.gff. If your workflow already
        produces GTF, point [reference].gpe001970_v5_gtf at it in
        config/DMP/00_common.toml; otherwise leave the .gff path as-is.

  [3] (Optional) Add Pfam domain TSV for stage [6]:
        [reference].dmp_domain_tsv = "II_INPUTS/DMP/Pfam/DMP_domains.tsv"
        in config/DMP/00_common.toml. Without it, stage [6] will log an
        empty domain dict and c_domain will always be 0.

  [4] (Optional) Enable plant-trained scorers:
        Re-run with --with-plant-scorers, then add dispatch branches for
        'DeepCRISPR' / 'CRISPR-Local' to
        modules/09_crispr_analysis/v2/02_rescore_ontarget.py.

  -------------------------------------------------------

  To activate:   conda activate $ENV_NAME
  To deactivate: conda deactivate
  To update:     bash setup_conda_crispr_v3.sh --update
  To remove:     conda env remove -n $ENV_NAME && conda env remove -n $ENV_INDELPHI

EOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setup finished -- full log: $LOG_FILE"
