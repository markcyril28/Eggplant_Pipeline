#!/bin/bash
# ============================================================================
# OpenFold3 Conda Environment Setup -- AMD MI210 Instinct (ROCm)
# ============================================================================
# Creates a dedicated 'openfold3' conda env for running OpenFold3 inference
# on the AMD MI210 Instinct GPU (gfx90a, 64 GB HBM2e) via ROCm.
#
# Source: https://github.com/aqlaboratory/openfold-3 (v0.4.1+)
# Install method: pip install openfold3 (official PyPI package)
#
# What it does:
#   1. Verifies ROCm >= 7.2 is installed on the host
#   2. Creates a fresh conda env (Python 3.13)
#   3. Installs ROCm-enabled PyTorch (rocm7.2 wheel index)
#   4. Installs openfold3 from PyPI
#   5. Runs validate-openfold3-rocm to confirm GPU access
#   6. Runs setup_openfold to download model weights (~2 GB)
#
# Requirements (host):
#   - ROCm >= 7.2  (https://rocm.docs.amd.com)
#   - AMD MI210 (gfx90a) or compatible CDNA2/CDNA3 GPU
#   - Conda or Mamba installed
#   - ~8 GB free disk space (env + weights)
#
# Usage:
#   bash setup_openfold3_rocm.sh                  # Full setup
#   bash setup_openfold3_rocm.sh --skip-weights   # Skip model weight download
#   bash setup_openfold3_rocm.sh --remove         # Remove and recreate env
#   bash setup_openfold3_rocm.sh --verify         # Verify existing install only
#
# After setup:
#   conda activate openfold3
#   run_openfold predict --query_json=<your_query.json>
#
# Invoking from another pipeline (e.g. CRISPR v3):
#   conda run -n openfold3 run_openfold predict --query_json=<query.json>
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------

ENV_NAME="openfold3"
PYTHON_VERSION="3.13"
ROCM_VERSION="7.2"
ROCM_TARGET_ARCH="gfx90a"           # MI210 Instinct
PYTORCH_WHEEL_INDEX="https://download.pytorch.org/whl/rocm${ROCM_VERSION}"

SKIP_WEIGHTS=false
MODE="create"                        # create | remove | verify

# ----------------------------------------------------------------------------
# ARGUMENT PARSING
# ----------------------------------------------------------------------------

for arg in "$@"; do
    case "$arg" in
        --skip-weights) SKIP_WEIGHTS=true ;;
        --remove)       MODE="remove" ;;
        --verify)       MODE="verify" ;;
        --help)
            sed -n '3,35p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: $0 [--skip-weights] [--remove] [--verify]" >&2
            exit 1
            ;;
    esac
done

# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
ok()   { echo "[$(date '+%H:%M:%S')] OK: $*"; }
fail() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

detect_conda() {
    if command -v mamba &>/dev/null; then
        echo "mamba"
    elif command -v conda &>/dev/null; then
        echo "conda"
    else
        fail "Neither mamba nor conda found. Install Miniconda or Mambaforge first."
    fi
}

env_exists() {
    conda env list 2>/dev/null | grep -qE "^${ENV_NAME}\s"
}

# ----------------------------------------------------------------------------
# ROCm VERIFICATION
# ----------------------------------------------------------------------------

verify_rocm() {
    log "Checking ROCm installation..."

    if ! command -v rocm-smi &>/dev/null; then
        fail "rocm-smi not found. Install ROCm ${ROCM_VERSION}: https://rocm.docs.amd.com"
    fi

    local rocm_path="${ROCM_PATH:-/opt/rocm}"
    if [[ ! -d "$rocm_path" ]]; then
        fail "ROCm directory not found at $rocm_path. Set ROCM_PATH if installed elsewhere."
    fi

    # Check version — rocm-smi --version prints e.g. "ROCm-SMI version: 7.2.0"
    local rocm_ver
    rocm_ver=$(rocm-smi --version 2>/dev/null | grep -oP '[\d]+\.[\d]+' | head -1 || echo "unknown")
    log "Detected ROCm version: ${rocm_ver}"

    # Check that MI210 (gfx90a) is visible
    if rocm-smi --showproductname 2>/dev/null | grep -qi "MI210\|Instinct"; then
        ok "AMD MI210 Instinct detected."
    else
        log "WARNING: MI210 not confirmed by rocm-smi. Proceeding anyway — verify device visibility."
    fi
}

# ----------------------------------------------------------------------------
# INSTALL
# ----------------------------------------------------------------------------

install_env() {
    local CONDA
    CONDA=$(detect_conda)
    log "Using: $CONDA"

    if [[ "$MODE" == "remove" ]] && env_exists; then
        log "Removing existing '$ENV_NAME' env..."
        $CONDA env remove -n "$ENV_NAME" -y
    fi

    if env_exists && [[ "$MODE" != "remove" ]]; then
        log "Env '$ENV_NAME' already exists. Skipping creation. Use --remove to recreate."
    else
        log "Creating conda env '$ENV_NAME' (Python ${PYTHON_VERSION})..."
        $CONDA create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    fi

    log "Installing ROCm-enabled PyTorch (rocm${ROCM_VERSION} wheel)..."
    conda run -n "$ENV_NAME" pip install \
        torch torchvision torchaudio \
        --index-url "$PYTORCH_WHEEL_INDEX"

    log "Installing openfold3..."
    conda run -n "$ENV_NAME" pip install openfold3

    # Set PYTORCH_ROCM_ARCH so torch compiles/JITs for MI210 specifically.
    # Without this, PyTorch may try to JIT for all supported archs, wasting time.
    log "Setting PYTORCH_ROCM_ARCH=${ROCM_TARGET_ARCH} in env activation..."
    local env_activate_dir
    env_activate_dir=$(conda run -n "$ENV_NAME" python3 -c \
        "import sys, pathlib; print(pathlib.Path(sys.prefix) / 'etc/conda/activate.d')")
    mkdir -p "$env_activate_dir"
    cat > "${env_activate_dir}/openfold3_rocm_env.sh" << EOF
#!/bin/bash
export ROCM_PATH="\${ROCM_PATH:-/opt/rocm}"
export HIP_PATH="\${HIP_PATH:-/opt/rocm}"
export PYTORCH_ROCM_ARCH="${ROCM_TARGET_ARCH}"
export PATH="\$ROCM_PATH/bin:\$PATH"
export LD_LIBRARY_PATH="\$ROCM_PATH/lib:\${LD_LIBRARY_PATH:-}"
EOF
    chmod +x "${env_activate_dir}/openfold3_rocm_env.sh"
    ok "ROCm env vars written to conda activation script."
}

# ----------------------------------------------------------------------------
# VALIDATE
# ----------------------------------------------------------------------------

validate_install() {
    log "Running validate-openfold3-rocm..."
    if conda run -n "$ENV_NAME" validate-openfold3-rocm; then
        ok "OpenFold3 ROCm validation passed."
    else
        fail "validate-openfold3-rocm failed. Check ROCm drivers and GPU visibility."
    fi
}

# ----------------------------------------------------------------------------
# MODEL WEIGHTS
# ----------------------------------------------------------------------------

download_weights() {
    if [[ "$SKIP_WEIGHTS" == "true" ]]; then
        log "Skipping model weight download (--skip-weights)."
        log "Run 'setup_openfold' after activating the env to download weights later."
        return
    fi

    log "Downloading model weights via setup_openfold (~2 GB)..."
    conda run -n "$ENV_NAME" setup_openfold
    ok "Model weights downloaded."
}

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

main() {
    echo "================================================================"
    echo " OpenFold3 ROCm Setup -- AMD MI210 Instinct (${ROCM_TARGET_ARCH})"
    echo " Mode: ${MODE} | Skip weights: ${SKIP_WEIGHTS}"
    echo "================================================================"

    if [[ "$MODE" == "verify" ]]; then
        validate_install
        exit 0
    fi

    verify_rocm
    install_env
    validate_install
    download_weights

    echo ""
    echo "================================================================"
    echo " Setup complete. Activate with:"
    echo "   conda activate ${ENV_NAME}"
    echo ""
    echo " Quick inference test:"
    echo "   run_openfold predict \\"
    echo "     --query_json=<path/to/query.json>"
    echo ""
    echo " From another pipeline script (no activation needed):"
    echo "   conda run -n ${ENV_NAME} run_openfold predict \\"
    echo "     --query_json=<path/to/query.json>"
    echo "================================================================"
}

main "$@"
