#!/bin/bash

#===============================================================================
# CONDA ENVIRONMENT SETUP: PRODIGY
#===============================================================================
# Creates/updates the conda environment used by Stage 14
# (14_interaction_Domain_Mapping.sh) for the interface_analysis and
# prodigy_dg operations:
#   - freesasa        : per-residue SASA / BSA for interface_analysis
#   - Biopython       : structure I/O (CIF / PDB) and chain handling
#   - PRODIGY         : protein-protein binding-affinity predictor (Xue 2016)
#
# Env name MUST match [tools].prodigy_conda_env and [conda_envs].interface_analysis /
# [conda_envs].prodigy_dg in 14_interaction_Domain_MappingCONFIG.toml (default: 'prodigy').
#
# Usage:
#   bash setup_conda_prodigy.sh                # create or update
#   bash setup_conda_prodigy.sh --update       # update existing env
#   bash setup_conda_prodigy.sh --restart      # remove and recreate
#   bash setup_conda_prodigy.sh --dry-run      # show actions without executing
#
# References:
#   PRODIGY    Xue 2016         10.1093/bioinformatics/btw514
#   freesasa   Mitternacht 2016 10.12688/f1000research.7931.1
#   Biopython  Cock 2009        10.1093/bioinformatics/btp163
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

ENV_NAME="prodigy"
PYTHON_VERSION="3.10"
CHANNELS="-c conda-forge -c bioconda"

# Conda-installable packages (everything that has a stable conda recipe).
# Versions pinned to a tested combination; bump only after re-testing
# Stage 14 interface_analysis / prodigy_dg on a known complex.
CONDA_PACKAGES=(
    "biopython=1.83"
    "freesasa=2.2.1"
    "numpy=1.26.4"
    "pandas=2.2.2"
    "scipy=1.13.1"
)

# pip-only packages (prodigy-prot is not on conda channels with reliable
# version pinning; the upstream PyPI release is the canonical distribution).
PIP_PACKAGES=(
    "prodigy-prot==2.1.4"
)

UPDATE_MODE=false
ENV_RESTART_MODE=false
DRY_RUN=false

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

for arg in "$@"; do
    case "$arg" in
        --update)   UPDATE_MODE=true ;;
        --restart)  ENV_RESTART_MODE=true ;;
        --dry-run)  DRY_RUN=true ;;
        -h|--help)
            sed -n '3,25p' "$0"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $arg" >&2
            echo "Valid: --update, --restart, --dry-run, --help" >&2
            exit 1
            ;;
    esac
done

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

log_info()  { echo "[INFO] $*"; }
log_warn()  { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }

check_command() { command -v "$1" &>/dev/null; }

run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] Would execute: $*"
    else
        "$@"
    fi
}

#===============================================================================
# DETECT PACKAGE MANAGER (prefer mamba > conda)
#===============================================================================

if check_command mamba; then
    PKG_MGR="mamba"
    log_info "Using mamba for faster installation"
else
    PKG_MGR="conda"
    log_warn "mamba not found; falling back to conda."
    log_warn "Install mamba for faster setup:  conda install -c conda-forge mamba"
fi

#===============================================================================
# CONDA INITIALIZATION
#===============================================================================

log_info "Initializing conda..."

if [[ -z "${CONDA_EXE:-}" ]]; then
    if   [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        log_error "Cannot find conda installation"
        exit 1
    fi
fi

if [[ "${WF_MANAGED_ENV:-}" != "true" ]]; then
    eval "$(conda shell.bash hook)"
fi

#===============================================================================
# ENVIRONMENT CREATION / UPDATE
#===============================================================================

log_info "========================================"
log_info "Setting up environment: $ENV_NAME"
log_info "========================================"
log_info "Conda packages: ${#CONDA_PACKAGES[@]}  (${CONDA_PACKAGES[*]})"
log_info "Pip packages:   ${#PIP_PACKAGES[@]}  (${PIP_PACKAGES[*]})"

if ${PKG_MGR} env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    log_info "Environment '${ENV_NAME}' already exists."

    if [[ "$ENV_RESTART_MODE" == true ]]; then
        log_info "Removing existing environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" env remove -n "${ENV_NAME}" -y
        log_info "Recreating environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y
    elif [[ "$UPDATE_MODE" == true ]]; then
        log_info "Update mode: reinstalling conda packages into existing env..."
    else
        log_info "Use --update to refresh packages, --restart to recreate."
        log_info "Continuing to (re-)install packages into existing env."
    fi
else
    log_info "Creating new environment '${ENV_NAME}'..."
    run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y
fi

#===============================================================================
# INSTALL CONDA PACKAGES
#===============================================================================

log_info "Installing conda packages into '${ENV_NAME}'..."
if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY RUN] Would execute: ${PKG_MGR} install -n ${ENV_NAME} ${CHANNELS} -y ${CONDA_PACKAGES[*]}"
else
    "${PKG_MGR}" install -n "${ENV_NAME}" ${CHANNELS} -y "${CONDA_PACKAGES[@]}" || {
        log_warn "${PKG_MGR} install failed; retrying with conda..."
        conda install -n "${ENV_NAME}" ${CHANNELS} -y "${CONDA_PACKAGES[@]}"
    }
fi

#===============================================================================
# INSTALL PIP PACKAGES
#===============================================================================

log_info "Installing pip packages into '${ENV_NAME}'..."
if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY RUN] Would execute: conda run -n ${ENV_NAME} python -m pip install --no-input ${PIP_PACKAGES[*]}"
else
    conda run -n "${ENV_NAME}" python -m pip install --no-input "${PIP_PACKAGES[@]}"
fi

#===============================================================================
# VERIFY
#===============================================================================

log_info "Verifying installation..."

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would verify Python imports + 'prodigy --help'"
    log_info "Done (dry-run)."
    exit 0
fi

VERIFY_PY=$(cat <<'PY'
import importlib, sys
mods = ["Bio", "freesasa", "numpy", "pandas", "scipy", "prodigy_prot"]
missing = []
for m in mods:
    try:
        importlib.import_module(m)
        print(f"  OK    {m}")
    except Exception as e:
        print(f"  FAIL  {m}: {e}")
        missing.append(m)
sys.exit(1 if missing else 0)
PY
)

if conda run -n "${ENV_NAME}" python -c "$VERIFY_PY"; then
    log_info "All Python modules importable."
else
    log_error "One or more modules failed to import. See output above."
    exit 1
fi

# PRODIGY CLI sanity check (the package installs a `prodigy` console script).
if conda run -n "${ENV_NAME}" prodigy --help >/dev/null 2>&1; then
    log_info "PRODIGY CLI ('prodigy --help') OK."
else
    log_warn "PRODIGY CLI did not respond to --help. compute_prodigy_dg.sh"
    log_warn "may still work if it invokes the Python API directly."
fi

log_info "========================================"
log_info "Environment '${ENV_NAME}' is ready."
log_info "Activate with:  conda activate ${ENV_NAME}"
log_info "Stage 14 picks it up via [conda_envs].interface_analysis and"
log_info "[conda_envs].prodigy_dg in 14_interaction_Domain_MappingCONFIG.toml."
log_info "========================================"
