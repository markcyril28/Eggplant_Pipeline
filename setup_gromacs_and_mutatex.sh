#!/usr/bin/env bash
#######################################################################
# GROMACS + MutateX + Stage 14 Environment Setup Script (Auto-Detect)
#
# Automatically detects GPU vendor (NVIDIA or AMD) and builds GROMACS
# with the appropriate backend (CUDA or HIP). Also installs the
# environments and binaries needed by Stage 14
# (14_Interaction_Domain_Mapping.sh): the GROMACS env (md_equilibration,
# currently postponed), the MutaTeX env (alanine scan + FoldX wrapper),
# the gmxmmpbsa env (binding-energy mmpbsa backend), Stage 14 Python
# deps in the egg env (PRODIGY, freesasa, gemmi, Biopython), and a
# scaffold for the FoldX binary at tools/foldx/.
#
# What it does:
# 1. Auto-detects GPU vendor and architecture
# 2. Creates a conda environment with all dependencies
# 3. Installs prebuilt CUDA GROMACS from conda-forge (or builds from source)
# 4. Installs Python analysis packages
# 5. Configures environment with GPU-optimized defaults
# 6. (Optional) Sets up MutateX with FoldX integration
# 7. (Optional) Sets up gmx_MMPBSA env for Stage 14 binding-energy backend
# 8. (Optional) Installs Stage 14 Python deps into the egg env
# 9. (Optional) Scaffolds tools/foldx/ for the FoldX binary drop-zone
#
# Supported GPUs:
#   NVIDIA: RTX 5050 (Blackwell), RTX 40xx (Ada), RTX 30xx (Ampere), etc.
#   AMD:    MI210 (gfx90a), MI100 (gfx908), MI50/60 (gfx906), etc.
#
# Requirements:
# - NVIDIA: Driver 570+ and/or CUDA Toolkit 12.x
# - AMD:    ROCm 6.x or 7.x
# - Conda/Mamba installed
# - ~10GB disk space
# - For --stage14: the egg env from setup_unified_conda_env.sh (for the
#   PRODIGY/freesasa/Biopython additions to apply)
#
# Usage:
#   ./setup_gromacs_and_mutatex.sh              # Full setup (auto-detect GPU)
#   ./setup_gromacs_and_mutatex.sh --deps-only  # Only install dependencies
#   ./setup_gromacs_and_mutatex.sh --build-only # Only build GROMACS
#   ./setup_gromacs_and_mutatex.sh --verify     # Only verify installation
#   ./setup_gromacs_and_mutatex.sh --force-cuda # Force CUDA backend
#   ./setup_gromacs_and_mutatex.sh --force-hip  # Force HIP backend
#   ./setup_gromacs_and_mutatex.sh --standalone [prefix]  # Build without conda
#   ./setup_gromacs_and_mutatex.sh --mutatex    # Also set up MutateX (FoldX)
#   ./setup_gromacs_and_mutatex.sh --mutatex-only       # Only set up MutateX
#   ./setup_gromacs_and_mutatex.sh --gmxmmpbsa          # Also set up gmx_MMPBSA env
#   ./setup_gromacs_and_mutatex.sh --gmxmmpbsa-only     # Only set up gmx_MMPBSA env
#   ./setup_gromacs_and_mutatex.sh --stage14            # Full Stage 14 install
#                                                       # (GROMACS + MutaTeX +
#                                                       #  gmxmmpbsa + egg-env
#                                                       #  Python deps + FoldX
#                                                       #  scaffold)
#   ./setup_gromacs_and_mutatex.sh --stage14-only       # Same minus GROMACS build
#                                                       # (use when MD is postponed)
#
# Stage 14 software map (14_Interaction_Domain_Mapping.sh + its TOML):
#   GROMACS                  -> built here as gromacs_CUDA / gromacs_HIP env
#                               (md_equilibration; currently POSTPONED in
#                               14_Interaction_Domain_MappingCONFIG.toml)
#   gmx_MMPBSA               -> gmxmmpbsa env (binding_energy mmpbsa backend;
#                               POSTPONED with GROMACS but kept ready)
#   MutaTeX                  -> PPI env (alanine_scan)
#   FoldX                    -> tools/foldx/foldx_<date> binary (binding_energy
#                               foldx backend + alanine_scan); MANUAL download
#                               from https://foldxsuite.crg.eu/
#   PRODIGY (prodigy-prot)   -> egg env (binding_energy prodigy backend +
#                               interface_analysis)
#   freesasa, gemmi, Biopython -> egg env (interface_analysis: BSA, contacts)
#   AlphaFold3               -> external, web-only submission at
#                               https://alphafoldserver.com/ (free non-
#                               commercial tier; NOT installed here)
#
# Standalone mode (no conda required):
#   Builds GROMACS directly to a custom prefix using system tools.
#   Useful for HPC/server environments without conda.
#   Default prefix: $HOME/.local/opt/gromacs-{version}
#
# After setup:
#   conda activate gromacs_CUDA      (NVIDIA) or gromacs_HIP (AMD)
#   source <prefix>/bin/GMXRC        (standalone mode)
#   gmx --version
#
#######################################################################

# NOTE: `nounset` (`-u`) is intentionally OMITTED here. Conda's compiler
# activate/deactivate scripts (e.g. envs/.../etc/conda/deactivate.d/
# deactivate-gcc_linux-64.sh) reference CONDA_BACKUP_* variables that may be
# unbound during in-place env re-activation. Under `set -u` those references
# abort the script with messages like:
#     deactivate-gcc_linux-64.sh: line 39: CONDA_BACKUP_CONDA_BUILD_SYSROOT: unbound variable
# The `safe_conda_activate` wrapper below also relaxes `-u`, but dropping it at
# the script level is the robust fix because conda's hook function can be
# invoked from multiple entry points.
set -eo pipefail

# Script directory — resolve once here, used throughout
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------

# Conda environment names
ENV_NAME_CUDA="gromacs_CUDA"  # Name for NVIDIA/CUDA environment
ENV_NAME_HIP="gromacs_HIP"    # Name for AMD/ROCm environment

# GROMACS version
GROMACS_VERSION="2026.0"
GROMACS_URL="https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz"

# Build directory — use real disk (not /tmp tmpfs) to avoid OOM on WSL2
BUILD_DIR="${HOME}/.local/share/gromacs-build"

# Number of parallel jobs for compilation
NJOBS=$(nproc 2>/dev/null || echo 8)

# ROCm path (for AMD GPUs)
ROCM_PATH="${ROCM_PATH:-/opt/rocm}"

# CUDA path (for NVIDIA GPUs, will auto-detect)
CUDA_PATH="${CUDA_HOME:-${CUDA_PATH:-}}"

# MutateX (FoldX wrapper) settings
MUTATEX_ENV_NAME="PPI"              # Conda environment for MutateX
MUTATEX_PYTHON_VERSION="3.11"       # Python version for MutateX env
SETUP_MUTATEX=false                  # Whether to set up MutateX

# gmx_MMPBSA settings (Stage 14 binding-energy mmpbsa backend)
# Name MUST match [tools].gmx_mmpbsa_env in 14_Interaction_Domain_MappingCONFIG.toml
GMXMMPBSA_ENV_NAME="gmxmmpbsa"
GMXMMPBSA_PYTHON_VERSION="3.10"     # gmx_MMPBSA pins ParmEd, easier on 3.10
SETUP_GMXMMPBSA=true                # Whether to set up gmx_MMPBSA env

# Stage 14 (14_Interaction_Domain_Mapping.sh) Python deps in the egg env.
# Name MUST match [tools].prodigy_conda_env in 14_Interaction_Domain_MappingCONFIG.toml
STAGE14_PRODIGY_ENV="egg"
SETUP_STAGE14=false                  # Full Stage 14 install (implies mutatex + gmxmmpbsa)
SETUP_STAGE14_ONLY=true             # Stage 14 install minus GROMACS build

# Auto-detection outputs (populated by detect_gpu):
GPU_BACKEND=""         # "CUDA" or "HIP"
GPU_ARCH=""            # e.g. "120" for CUDA, "gfx90a" for HIP
GPU_NAME=""            # Human-readable GPU name
GPU_VRAM=""            # VRAM in MiB (e.g. "8192 MiB")
ENV_NAME=""            # Conda environment name (set from ENV_NAME_CUDA/HIP)
USE_CONDA_CUDA=false   # Whether to install CUDA via conda
STANDALONE_MODE=false  # Build without conda (system tools only)
INSTALL_PREFIX=""      # Custom install prefix (standalone mode)
FORCE_BACKEND=""       # Force "CUDA" or "HIP" (set via --force-cuda / --force-hip)

#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------

# -- Installation log file (set by setup_installation_logging, empty = no file logging)
INSTALL_LOG_FILE=""

strip_ansi() { sed 's/\x1b\[[0-9;]*m//g'; }

setup_installation_logging() {
    local log_dir="$SCRIPT_DIR/logs/installation_logs"
    mkdir -p "$log_dir"
    local timestamp
    timestamp=$(date '+%Y%m%d_%H%M%S')
    INSTALL_LOG_FILE="$log_dir/setup_gromacs_${timestamp}.log"

    # Redirect all stdout+stderr through tee so every line hits the log file
    exec > >(tee >(strip_ansi >> "$INSTALL_LOG_FILE")) 2>&1
}

teardown_installation_logging() {
    # Close the process-substitution tee by restoring original FDs
    exec 1>&- 2>&-          # close tee pipes
    exec 1>/dev/tty 2>&1    # restore terminal output
    if [[ -n "$INSTALL_LOG_FILE" && -f "$INSTALL_LOG_FILE" ]]; then
        echo ""
        echo "[LOG] Installation log saved to: $INSTALL_LOG_FILE"
    fi
}

log() {
    echo -e "\n\033[1;32m[$(date '+%H:%M:%S')]\033[0m $1"
}
log_warn() {
    echo -e "\n\033[1;33m[WARNING]\033[0m $1" >&2
}
log_error() {
    echo -e "\n\033[1;31m[ERROR]\033[0m $1" >&2
    exit 1
}

# Conda's activation scripts (notably binutils-feedstock's
# activate-binutils_linux-64.sh) reference variables like $ADDR2LINE / $AR / $AS
# that are unset on a fresh shell. With `set -u` (from `set -euo pipefail` at
# script top) those bare references become fatal:
#     activate-binutils_linux-64.sh: line 68: ADDR2LINE: unbound variable
# Wrap every `conda activate` (and `conda deactivate`) through this helper so
# nounset is relaxed only for the duration of the activation, then restored.
safe_conda_activate() {
    local _had_u=0
    [[ $- == *u* ]] && _had_u=1
    set +u
    conda activate "$@"
    local _rc=$?
    (( _had_u )) && set -u
    return $_rc
}

safe_conda_deactivate() {
    local _had_u=0
    [[ $- == *u* ]] && _had_u=1
    set +u
    conda deactivate
    local _rc=$?
    (( _had_u )) && set -u
    return $_rc
}

#------------------------------------------------------------------------------
# GPU AUTO-DETECTION
#------------------------------------------------------------------------------

detect_nvidia_arch() {
    # Map NVIDIA compute capability from GPU name
    # Uses nvidia-smi to get compute capability if available
    local cc=""

    # Try to get compute capability directly
    if command -v nvidia-smi &> /dev/null; then
        cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')
    fi

    # If that didn't work, infer from GPU name
    if [ -z "$cc" ] || [ "$cc" = "N/A" ]; then
        local name="$1"
        case "$name" in
            *"RTX 50"*|*"Blackwell"*)   cc="120" ;;
            *"RTX 40"*|*"Ada"*|*"L40"*) cc="89"  ;;
            *"A100"*|*"A30"*)              cc="80"  ;;
            *"RTX 30"*|*"A10"*)            cc="86"  ;;
            *"RTX 20"*|*"Turing"*|*"T4"*) cc="75" ;;
            *"GTX 16"*)                  cc="75"  ;;
            *"V100"*)                    cc="70"  ;;
            *)                           cc="86"  ;; # Safe default (Ampere)
        esac
        log_warn "Could not query compute capability directly, inferred sm_${cc} from GPU name"
    fi

    echo "$cc"
}

detect_amd_arch() {
    # Detect AMD GPU architecture from rocm-smi or rocminfo
    local arch=""

    if command -v rocminfo &> /dev/null; then
        arch=$(rocminfo 2>/dev/null | grep -Eo 'gfx[0-9]+[a-z]?' | head -1)
    fi

    if [ -z "$arch" ] && command -v rocm-smi &> /dev/null; then
        # Try to infer from device name
        local name
        name=$(rocm-smi --showproductname 2>/dev/null | grep -i "GPU" | head -1)
        case "$name" in
            *"MI210"*|*"MI250"*) arch="gfx90a" ;;
            *"MI100"*)           arch="gfx908" ;;
            *"MI50"*|*"MI60"*)   arch="gfx906" ;;
            *"7900"*)            arch="gfx1100" ;;
            *"6900"*|*"6800"*)   arch="gfx1030" ;;
            *)                   arch="gfx90a"  ;; # Safe default
        esac
        log_warn "Could not query GPU arch directly, inferred ${arch} from device name"
    fi

    if [ -z "$arch" ]; then
        arch="gfx90a"
        log_warn "Could not detect AMD GPU arch; defaulting to gfx90a. Use --force-hip to override."
    fi
    echo "$arch"
}

detect_gpu() {
    log "Auto-detecting GPU..."

    local has_nvidia=false
    local has_amd=false

    # Check for NVIDIA GPU (nvidia-smi -L is faster than bare nvidia-smi)
    if command -v nvidia-smi &> /dev/null && nvidia-smi -L &> /dev/null; then
        has_nvidia=true
    fi

    # Check for AMD GPU (ROCm)
    if [ -d "$ROCM_PATH" ] && (command -v rocm-smi &> /dev/null || [ -f "$ROCM_PATH/bin/hipcc" ]); then
        has_amd=true
    fi

    # Handle force flags
    if [ "${FORCE_BACKEND:-}" = "CUDA" ]; then
        has_nvidia=true
        has_amd=false
        log "Forced CUDA backend"
    elif [ "${FORCE_BACKEND:-}" = "HIP" ]; then
        has_nvidia=false
        has_amd=true
        log "Forced HIP backend"
    fi

    # Decide backend
    if $has_nvidia && $has_amd; then
        log_warn "Both NVIDIA and AMD GPUs detected. Defaulting to NVIDIA (CUDA)."
        log_warn "Use --force-hip to build for AMD instead."
        has_amd=false
    fi

    if $has_nvidia; then
        GPU_BACKEND="CUDA"
        # Single nvidia-smi call for name + VRAM (avoids 3 sequential invocations)
        local _nv_info
        _nv_info=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1)
        GPU_NAME=$(cut -d',' -f1 <<< "$_nv_info" | xargs)
        GPU_VRAM=$(cut -d',' -f2 <<< "$_nv_info" | xargs)
        GPU_ARCH=$(detect_nvidia_arch "$GPU_NAME")
        ENV_NAME="$ENV_NAME_CUDA"
        BUILD_DIR="${HOME}/.local/share/gromacs-build-cuda"

        log "Detected NVIDIA GPU: $GPU_NAME ($GPU_VRAM)"
        log "CUDA compute capability: sm_${GPU_ARCH}"

        # Detect CUDA toolkit
        detect_cuda_toolkit

    elif $has_amd; then
        GPU_BACKEND="HIP"
        if command -v rocm-smi &> /dev/null; then
            GPU_NAME=$(rocm-smi --showproductname 2>/dev/null | grep -i "card\|GPU" | head -1 | sed 's/.*: *//')
        fi
        GPU_NAME="${GPU_NAME:-AMD GPU}"
        GPU_ARCH=$(detect_amd_arch)
        ENV_NAME="$ENV_NAME_HIP"
        BUILD_DIR="${HOME}/.local/share/gromacs-build-hip"

        local rocm_ver
        rocm_ver=$(cat "$ROCM_PATH/.info/version" 2>/dev/null || echo "unknown")
        log "Detected AMD GPU: $GPU_NAME"
        log "HIP architecture: $GPU_ARCH"
        log "ROCm version: $rocm_ver"

        # Validate ROCm
        if [ ! -f "$ROCM_PATH/bin/hipcc" ] && [ ! -f "$ROCM_PATH/bin/amdclang++" ]; then
            log_error "ROCm HIP compiler not found. Please install ROCm."
        fi
    else
        log_error "No supported GPU detected.\n  NVIDIA: install drivers + nvidia-smi\n  AMD: install ROCm"
    fi

    echo ""
    log "========================================="
    log "  Backend:      $GPU_BACKEND"
    log "  GPU:          $GPU_NAME"
    log "  Architecture: ${GPU_BACKEND}=${GPU_ARCH}"
    log "  Environment:  $ENV_NAME"
    log "========================================="
}

detect_cuda_toolkit() {
    log "Detecting CUDA toolkit..."

    local cuda_candidates=(
        "$CUDA_PATH"
        "$HOME/cuda-12.8"
        "$HOME/cuda"
        "/usr/local/cuda"
        "/usr/local/cuda-12"
        "/usr/local/cuda-12.8"
        "/usr"
    )

    for candidate in "${cuda_candidates[@]}"; do
        if [ -n "$candidate" ] && [ -f "$candidate/bin/nvcc" ]; then
            CUDA_PATH="$candidate"
            break
        fi
    done

    if [ -z "$CUDA_PATH" ] || [ ! -f "$CUDA_PATH/bin/nvcc" ]; then
        log_warn "System CUDA toolkit not found. Will install via conda."
        USE_CONDA_CUDA=true
    else
        USE_CONDA_CUDA=false
        log "CUDA toolkit: $CUDA_PATH"
        log "nvcc: $($CUDA_PATH/bin/nvcc --version 2>&1 | grep release | awk '{print $NF}')"
    fi
}

install_cuda_runfile() {
    # Fallback CUDA installer: downloads NVIDIA's official runfile and installs
    # the toolkit to $HOME/cuda-12.8 (no sudo required).
    # Called when conda's SAT solver cannot resolve CUDA packages.
    local cuda_ver="12.8.1"
    local driver_ver="570.124.06"
    local runfile="${TMPDIR:-/tmp}/cuda_${cuda_ver}_linux.run"
    local install_dir="$HOME/cuda-12.8"
    local url="https://developer.download.nvidia.com/compute/cuda/${cuda_ver}/local_installers/cuda_${cuda_ver}_${driver_ver}_linux.run"

    if [ -f "$install_dir/bin/nvcc" ]; then
        log "CUDA toolkit already installed at $install_dir"
        CUDA_PATH="$install_dir"
        USE_CONDA_CUDA=false
        return 0
    fi

    # Download if not present or incomplete (expect ~5 GB)
    local expected_min_bytes=5000000000
    if [ ! -f "$runfile" ] || [ "$(stat -c%s "$runfile" 2>/dev/null || echo 0)" -lt "$expected_min_bytes" ]; then
        log "Downloading CUDA ${cuda_ver} runfile installer (~5 GB)..."
        rm -f "$runfile"
        wget --progress=bar:force:noscroll "$url" -O "$runfile" \
            || log_error "CUDA runfile download failed. Check network and retry."
    else
        log "CUDA runfile already downloaded: $runfile"
    fi

    chmod +x "$runfile"
    log "Installing CUDA toolkit to $install_dir (no sudo required)..."
    "$runfile" --toolkit --toolkitpath="$install_dir" --silent --override \
        || log_error "CUDA runfile installation failed."

    if [ -f "$install_dir/bin/nvcc" ]; then
        CUDA_PATH="$install_dir"
        USE_CONDA_CUDA=false
        export PATH="$install_dir/bin:$PATH"
        export LD_LIBRARY_PATH="$install_dir/lib64:${LD_LIBRARY_PATH:-}"
        log "CUDA toolkit installed: $($install_dir/bin/nvcc --version 2>&1 | grep release)"
    else
        log_error "CUDA runfile install completed but nvcc not found at $install_dir/bin/nvcc"
    fi
}

#------------------------------------------------------------------------------
# PREREQUISITE CHECKS (standalone mode)
#------------------------------------------------------------------------------

check_prerequisites() {
    log "Checking prerequisites..."

    # Common checks
    if ! command -v cmake &> /dev/null; then
        log_error "cmake not found. Please install cmake >= 3.18"
    fi
    log "CMake: $(cmake --version | head -1)"

    if ! command -v ninja &> /dev/null; then
        log_error "ninja not found. Please install ninja-build."
    fi
    log "Ninja: $(ninja --version)"

    if [ "$GPU_BACKEND" = "HIP" ]; then
        # ROCm checks
        if [ ! -d "$ROCM_PATH" ]; then
            log_error "ROCm not found at $ROCM_PATH\nPlease install ROCm: https://rocm.docs.amd.com/en/latest/"
        fi

        local rocm_ver
        rocm_ver=$(cat "$ROCM_PATH/.info/version" 2>/dev/null || echo "unknown")
        log "ROCm version: $rocm_ver"

        if ! command -v hipcc &> /dev/null && [ ! -f "$ROCM_PATH/bin/hipcc" ]; then
            log_error "hipcc not found. Is ROCm properly installed?"
        fi
        log "HIP compiler: $(hipcc --version 2>/dev/null | head -1 || echo 'available')"

        # Show available GPUs
        if command -v rocm-smi &> /dev/null; then
            log "Available GPUs:"
            rocm-smi --showproductname 2>/dev/null | grep -E "GPU|gfx" || true
        fi

    elif [ "$GPU_BACKEND" = "CUDA" ]; then
        # NVIDIA checks
        if ! command -v nvidia-smi &> /dev/null; then
            log_error "nvidia-smi not found. Please install NVIDIA drivers."
        fi
        log "GPU: $GPU_NAME ($GPU_VRAM)"

        if [ "$STANDALONE_MODE" = true ]; then
            if [ -z "$CUDA_PATH" ] || [ ! -f "$CUDA_PATH/bin/nvcc" ]; then
                log_error "CUDA toolkit not found. Install CUDA or set CUDA_HOME."
            fi
            log "nvcc: $($CUDA_PATH/bin/nvcc --version 2>&1 | grep release | awk '{print $NF}')"
        fi
    fi
}

#------------------------------------------------------------------------------
# CONDA CHECK
#------------------------------------------------------------------------------

check_conda() {
    log "Checking conda..."

    if ! command -v conda &> /dev/null; then
        log_error "Conda not found. Please install Miniconda or Anaconda."
    fi

    log "Conda: $(conda --version)"

    if ! command -v mamba &> /dev/null; then
        log "Installing mamba for faster package resolution..."
        conda install -y -c conda-forge mamba || log_warn "Could not install mamba, using conda"
    fi

    if command -v mamba &> /dev/null; then
        PKG_MGR="mamba"
    else
        PKG_MGR="conda"
    fi
    log "Package manager: $PKG_MGR"
}

#------------------------------------------------------------------------------
# INSTALL DEPENDENCIES
#------------------------------------------------------------------------------

install_dependencies() {
    log "Creating/updating conda environment: $ENV_NAME"

    if conda env list | grep -q "^${ENV_NAME} "; then
        log "Environment '$ENV_NAME' already exists."
    else
        log "Creating new environment..."
        conda create -n "$ENV_NAME" python=3.11 -y
    fi

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$ENV_NAME"

    # ---- Skip conda packages if key indicators are already present ----
    local _need_conda=false
    for _pkg in cmake ninja fftw numpy scipy pandas; do
        if ! conda list -n "$ENV_NAME" "^${_pkg}$" 2>/dev/null | grep -q "$_pkg"; then
            _need_conda=true; break
        fi
    done

    if $_need_conda; then
        log "Installing all conda dependencies (single solver pass)..."
        # vmd: required by Step 6 visualize_results for headless MP4 movie rendering
        # ffmpeg: encodes the per-frame TGAs that VMD writes into H.264 MP4
        $PKG_MGR install -y -c conda-forge \
            cmake ninja make pkg-config fftw libiconv \
            openmpi openmpi-mpicc \
            numpy scipy pandas matplotlib seaborn biopython networkx \
            pillow imageio imageio-ffmpeg scikit-learn scikit-image \
            tqdm pyyaml h5py gnuplot ffmpeg vmd \
            || {
                [ "$GPU_BACKEND" = "HIP" ] && log_error "Failed to install required dependencies."
                log_warn "Some optional packages unavailable; retrying with essential packages only..."
                $PKG_MGR install -y -c conda-forge \
                    cmake ninja make pkg-config fftw libiconv \
                    openmpi openmpi-mpicc \
                    numpy scipy pandas matplotlib seaborn biopython \
                    networkx pillow tqdm pyyaml h5py
                log_warn "VMD/ffmpeg not installed by fallback path. Install manually with:"
                log_warn "  $PKG_MGR install -n $ENV_NAME -c conda-forge vmd ffmpeg"
            }
    else
        log "Conda dependencies already installed, skipping."
    fi

    # CUDA toolkit via conda if needed — install ONLY the compiler + runtime.
    # IMPORTANT: Do NOT install cuda-libraries-dev or cuda-toolkit here!
    # Those metapackages pull hundreds of transitive deps whose SAT resolution
    # can hang mamba for 8+ hours and consume >10 GB RAM.
    # GROMACS only needs: nvcc (compiler), cudart (runtime), cccl (C++ headers),
    # and libcufft-dev (required for PME GPU acceleration).
    #
    # Strategy: try conda with a 5-minute timeout. If the solver hangs (common
    # with 200+ packages already installed), fall back to NVIDIA's official
    # runfile installer which installs to $HOME/cuda-12.8 without sudo.
    if [ "$GPU_BACKEND" = "CUDA" ] && [ "$USE_CONDA_CUDA" = true ]; then
        log "Installing minimal CUDA build dependencies via conda-forge (5 min timeout)..."
        if timeout 300 $PKG_MGR install -y -c conda-forge \
                cuda-nvcc cuda-cudart-dev cuda-cccl libcufft-dev 2>&1; then
            log "CUDA packages installed via conda."
        else
            log_warn "Conda CUDA install failed or timed out (SAT solver hang is common)."
            log "Falling back to NVIDIA runfile installer..."
            install_cuda_runfile
        fi
    fi

    # ---- Skip pip packages if key indicators are already present ----
    if python3 -c "import MDAnalysis, prolif, biotite" 2>/dev/null; then
        log "Pip analysis packages already installed, skipping."
    else
        log "Installing Python analysis packages..."
        pip install --quiet \
            MDAnalysis prolif py3Dmol nglview moviepy plotly bokeh pdb-tools biotite
    fi

    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME
    log "Dependencies installed to: $CONDA_PREFIX"
}

#------------------------------------------------------------------------------
# BUILD GROMACS
#------------------------------------------------------------------------------

download_gromacs() {
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    if [ ! -f "gromacs-${GROMACS_VERSION}.tar.gz" ]; then
        log "Downloading GROMACS ${GROMACS_VERSION}..."
        wget -q "$GROMACS_URL" -O "gromacs-${GROMACS_VERSION}.tar.gz" \
            || log_error "Download failed — check GROMACS_VERSION ($GROMACS_VERSION)"
    fi

    if [ ! -d "gromacs-${GROMACS_VERSION}" ]; then
        log "Extracting..."
        tar -xzf "gromacs-${GROMACS_VERSION}.tar.gz"
    fi
}

# Build FFTW 3.3.10 from source using the given CC/CXX compilers.
# GMX_BUILD_OWN_FFTW=ON requires Make (not Ninja), so we pre-build FFTW here
# and pass the resulting library to GROMACS via -DFFTWF_LIBRARY / -DFFTWF_INCLUDE_DIR.
# Results are cached in BUILD_DIR/fftw-install so repeated builds skip this step.
build_fftw_from_source() {
    local cc="$1"    # C compiler to use
    local cxx="$2"   # C++ compiler (unused by FFTW but kept for symmetry)
    local install_dir="$BUILD_DIR/fftw-install"

    # Validate cached build: the .a must exist AND contain the symbols GROMACS needs.
    if [ -f "$install_dir/lib/libfftw3f.a" ]; then
        if nm "$install_dir/lib/libfftw3f.a" 2>/dev/null | grep -q 'T fftwf_plan_many_dft$'; then
            log "FFTW3 already built at $install_dir (symbols OK), skipping."
            return 0
        else
            log_warn "FFTW3 cache at $install_dir is corrupt or incomplete — rebuilding."
            rm -rf "$install_dir" "$BUILD_DIR/fftw-build"
        fi
    fi

    local fftw_version="3.3.10"
    local fftw_tar="$BUILD_DIR/fftw-${fftw_version}.tar.gz"
    local fftw_url="https://fftw.org/fftw-${fftw_version}.tar.gz"

    mkdir -p "$BUILD_DIR"

    if [ ! -f "$fftw_tar" ]; then
        log "Downloading FFTW ${fftw_version}..."
        wget -q "$fftw_url" -O "$fftw_tar" \
            || log_error "FFTW download failed. Check network and retry."
    fi

    local fftw_src="$BUILD_DIR/fftw-${fftw_version}"
    if [ ! -d "$fftw_src" ]; then
        tar -xzf "$fftw_tar" -C "$BUILD_DIR"
    fi

    log "Building FFTW3 (single precision) from source..."
    mkdir -p "$BUILD_DIR/fftw-build"
    cd "$BUILD_DIR/fftw-build"

    CC="$cc" "$fftw_src/configure" \
        --prefix="$install_dir" \
        --enable-float \
        --enable-sse2 \
        --enable-avx \
        --enable-avx2 \
        --enable-shared \
        --enable-static \
        --disable-fortran \
        --quiet
    make -j"$NJOBS"
    make install

    log "FFTW3 installed to $install_dir"
    cd - > /dev/null
}

build_gromacs_cuda() {
    local conda_prefix="$1"
    local mpi_flag="$2"       # "ON" or "OFF"
    local build_suffix="$3"   # e.g. "cuda" or "cuda_mpi"

    local build_dir="$BUILD_DIR/gromacs-${GROMACS_VERSION}/build_${build_suffix}"
    rm -rf "$build_dir"
    mkdir -p "$build_dir"
    cd "$build_dir"

    # CUDA toolkit path for cmake
    local cuda_cmake_args=""
    if [ "$USE_CONDA_CUDA" = true ]; then
        # conda-forge cuda-nvcc places nvcc in bin/ but headers/libs under targets/;
        # CUDAToolkit_ROOT lets CMake's FindCUDAToolkit handle the layout correctly
        [ -f "$conda_prefix/bin/nvcc" ] && cuda_cmake_args="-DCUDAToolkit_ROOT=$conda_prefix"
    elif [ -n "$CUDA_PATH" ] && [ -d "$CUDA_PATH" ]; then
        cuda_cmake_args="-DCUDAToolkit_ROOT=$CUDA_PATH"
    fi

    # Ensure C and C++ compilers are from the same toolchain.
    # Conda may inject its own GCC for CC but leave CXX pointing at a different
    # system version, causing ABI mismatches.  Use a matched system pair.
    local compiler_args=""
    if [ "$STANDALONE_MODE" != true ]; then
        local _sys_cc _sys_cxx
        _sys_cc=$(command -v gcc 2>/dev/null || true)
        _sys_cxx=$(command -v g++ 2>/dev/null || true)
        # Prefer the system compiler pair — they are always the same GCC version
        if [ -n "$_sys_cc" ] && [ -n "$_sys_cxx" ]; then
            compiler_args="-DCMAKE_C_COMPILER=$_sys_cc -DCMAKE_CXX_COMPILER=$_sys_cxx"
            local _cc_ver _cxx_ver
            _cc_ver=$("$_sys_cc" -dumpversion 2>/dev/null || echo "?")
            _cxx_ver=$("$_sys_cxx" -dumpversion 2>/dev/null || echo "?")
            log "Using matched system compilers: gcc ${_cc_ver} / g++ ${_cxx_ver}"
        fi
    fi

    # Check for cuda_profiler_api.h — missing when CUDA toolkit is incomplete
    # (e.g. runfile install without CUPTI).  Disable profiling if absent.
    local profiler_args=""
    local _cuda_root="${CUDA_PATH:-}"
    [ "$USE_CONDA_CUDA" = true ] && _cuda_root="$conda_prefix"
    if [ -n "$_cuda_root" ]; then
        if ! find "$_cuda_root" -name "cuda_profiler_api.h" -print -quit 2>/dev/null | grep -q .; then
            log_warn "cuda_profiler_api.h not found in $_cuda_root — disabling CUDA profiling"
            profiler_args="-DGMX_CUDA_PROFILING=OFF"
        fi
    fi

    log "Configuring CMake (CUDA, MPI=$mpi_flag)..."

    # FFTW: GMX_BUILD_OWN_FFTW=ON requires Make (not Ninja), so we pre-build
    # FFTW from source ourselves and hand the result to CMake.  This is
    # necessary when system compilers are used because conda's libfftw3f was
    # compiled with conda's GCC and its symbols fail CMake link tests with the
    # system linker.  When conda compilers are active (compiler_args empty),
    # use conda's library directly to avoid the extra download+build.
    local fftw_args=""
    if [ "$STANDALONE_MODE" = true ] || [ -n "$compiler_args" ]; then
        local _fftw_cc="${_sys_cc:-gcc}"
        build_fftw_from_source "$_fftw_cc" "${_sys_cxx:-g++}"
        local _fftw_dir="$BUILD_DIR/fftw-install"
        fftw_args="-DGMX_BUILD_OWN_FFTW=OFF -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_INCLUDE_DIR=$_fftw_dir/include -DFFTWF_LIBRARY=$_fftw_dir/lib/libfftw3f.so"
    else
        fftw_args="-DGMX_BUILD_OWN_FFTW=OFF -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_INCLUDE_DIR=$conda_prefix/include -DFFTWF_LIBRARY=$conda_prefix/lib/libfftw3f.so"
    fi

    # Tell CMake where conda libraries live (HDF5, MPI, etc. still come from
    # conda even when FFTW is built from source).
    local prefix_path_args=""
    if [ "$STANDALONE_MODE" != true ]; then
        prefix_path_args="-DCMAKE_PREFIX_PATH=$conda_prefix"
    fi

    # Conda activation injects CFLAGS/CXXFLAGS/LDFLAGS that reference its
    # cross-compilation sysroot (/lib64/libm.so.6 etc.).  When building with
    # the system GCC these paths don't exist and break CMake's
    # check_library_exists() link tests.  Unsetting them lets CMake use its
    # own defaults with the system toolchain.
    #
    # Additionally, conda's own CMake (4.x) embeds the conda sysroot into
    # CMAKE_SYSTEM_LIBRARY_PATH, causing the same /lib64/ failures even after
    # unsetting env vars.  Use the system CMake (/usr/bin/cmake) to avoid this.
    local _cmake="cmake"
    if [ -n "$compiler_args" ]; then
        unset CFLAGS CXXFLAGS LDFLAGS FFLAGS DEBUG_CFLAGS DEBUG_CXXFLAGS DEBUG_FFLAGS 2>/dev/null || true
        unset CMAKE_PREFIX_PATH CMAKE_ARGS CONDA_BUILD_SYSROOT HOST CONDA_TOOLCHAIN_HOST CONDA_TOOLCHAIN_BUILD 2>/dev/null || true
        if [ -x /usr/bin/cmake ]; then
            _cmake=/usr/bin/cmake
            log "Using system CMake: $($_cmake --version | head -1)"
        fi
    fi

    # Clear stale CMake cache (previous runs may have cached bad paths).
    rm -f CMakeCache.txt 2>/dev/null || true

    $_cmake .. \
        -GNinja \
        -DCMAKE_INSTALL_PREFIX="$conda_prefix" \
        -DCMAKE_BUILD_TYPE=Release \
        $prefix_path_args \
        $compiler_args \
        $fftw_args \
        -DGMX_GPU=CUDA \
        -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH" \
        $cuda_cmake_args \
        $profiler_args \
        -DGMX_MPI="$mpi_flag" \
        -DGMX_OPENMP=ON \
        -DGMX_SIMD=AVX2_256 \
        -DREGRESSIONTEST_DOWNLOAD=OFF \
        -DGMX_BUILD_HELP=OFF \
        2>&1 | tee cmake_output.log

    log "Building (using $NJOBS cores)..."
    ninja -j${NJOBS} 2>&1 | tee build_output.log

    log "Installing..."
    ninja install 2>&1 | tee install_output.log
}

build_gromacs_hip() {
    local conda_prefix="$1"
    local mpi_flag="$2"
    local build_suffix="$3"

    local build_dir="$BUILD_DIR/gromacs-${GROMACS_VERSION}/build_${build_suffix}"
    rm -rf "$build_dir"
    mkdir -p "$build_dir"
    cd "$build_dir"

    # Use ROCm compilers
    export CC="$ROCM_PATH/bin/amdclang"
    export CXX="$ROCM_PATH/bin/amdclang++"
    export HIP_PATH="$ROCM_PATH"

    # Fallback to llvm clang if amdclang not available
    if [ ! -f "$CC" ]; then
        export CC="$ROCM_PATH/llvm/bin/clang"
        export CXX="$ROCM_PATH/llvm/bin/clang++"
    fi

    # FFTW: pre-build from source using the ROCm C compiler so the library is
    # ABI-compatible.  GMX_BUILD_OWN_FFTW=ON requires Make (not Ninja), so we
    # build FFTW ourselves and pass the result to GROMACS via -DFFTWF_LIBRARY.
    build_fftw_from_source "$CC" "$CXX"
    local _fftw_dir="$BUILD_DIR/fftw-install"
    local fftw_args="-DGMX_BUILD_OWN_FFTW=OFF -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_INCLUDE_DIR=$_fftw_dir/include -DFFTWF_LIBRARY=$_fftw_dir/lib/libfftw3f.so"

    # Tell CMake where conda libraries live (HDF5, MPI, etc.).
    local prefix_path_args=""
    if [ "$STANDALONE_MODE" != true ]; then
        prefix_path_args="-DCMAKE_PREFIX_PATH=$conda_prefix"
    fi

    # Clear conda's sysroot-based compiler/linker flags (same issue as CUDA
    # path — conda sysroot references /lib64/ which doesn't exist on Debian).
    # Also clear CMAKE_PREFIX_PATH/CONDA_BUILD_SYSROOT injected by activation:
    # conda activation appends ".../sysroot/usr" to CMAKE_PREFIX_PATH which
    # makes CMake add "-L .../sysroot/usr/lib" to link tests.  That directory
    # contains a libm.so linker script hardcoded to /lib64/libm.so.6, which
    # does not exist on Debian/Ubuntu (libm lives in /lib/x86_64-linux-gnu).
    # The result is ld.lld errors like "cannot open /lib64/libm.so.6" during
    # CMake's check_library_exists() for FFTW symbols.  We re-pass the conda
    # prefix explicitly via $prefix_path_args, so unsetting is safe.
    unset CFLAGS CXXFLAGS LDFLAGS FFLAGS DEBUG_CFLAGS DEBUG_CXXFLAGS DEBUG_FFLAGS 2>/dev/null || true
    unset CMAKE_PREFIX_PATH CMAKE_ARGS CONDA_BUILD_SYSROOT HOST CONDA_TOOLCHAIN_HOST CONDA_TOOLCHAIN_BUILD 2>/dev/null || true

    # Use system CMake to avoid conda CMake's sysroot injection.
    local _cmake="cmake"
    if [ -x /usr/bin/cmake ]; then
        _cmake=/usr/bin/cmake
        log "Using system CMake: $($_cmake --version | head -1)"
    fi

    # Clear stale CMake cache.
    rm -f CMakeCache.txt 2>/dev/null || true

    log "Configuring CMake (HIP, MPI=$mpi_flag)..."
    $_cmake .. \
        -GNinja \
        -DCMAKE_INSTALL_PREFIX="$conda_prefix" \
        -DCMAKE_BUILD_TYPE=Release \
        $prefix_path_args \
        $fftw_args \
        -DGMX_GPU=HIP \
        -DGMX_HIP_TARGET_ARCH="$GPU_ARCH" \
        -DCMAKE_HIP_COMPILER="$CXX" \
        -DCMAKE_HIP_ARCHITECTURES="$GPU_ARCH" \
        -DHIP_PATH="$ROCM_PATH" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DGMX_MPI="$mpi_flag" \
        -DGMX_OPENMP=ON \
        -DGMX_SIMD=AVX2_256 \
        -DREGRESSIONTEST_DOWNLOAD=OFF \
        -DGMX_BUILD_HELP=OFF \
        2>&1 | tee cmake_output.log

    log "Building (using $NJOBS cores)..."
    ninja -j${NJOBS} 2>&1 | tee build_output.log

    log "Installing..."
    ninja install 2>&1 | tee install_output.log
}

try_conda_prebuilt() {
    # Try installing a prebuilt CUDA GROMACS from conda-forge (much faster than
    # building from source).  Only applies to CUDA + conda mode.
    [ "$STANDALONE_MODE" = true ] && return 1
    [ "$GPU_BACKEND" != "CUDA" ] && return 1

    local prefix
    prefix=$(conda info --base)/envs/$ENV_NAME
    local build_str="nompi_cuda"

    log "Attempting prebuilt CUDA GROMACS from conda-forge (faster than source build)..."
    if timeout 300 $PKG_MGR install -n "$ENV_NAME" -y -c conda-forge \
            "gromacs=${GROMACS_VERSION}=${build_str}*" 2>&1; then
        # Verify the installed binary actually has CUDA backend
        local _gmx="$prefix/bin/gmx"
        if [ -f "$_gmx" ]; then
            local _backend
            _backend=$("$_gmx" --version 2>&1 | awk '/GPU support:/{print $NF}')
            if [ "$_backend" = "CUDA" ]; then
                log "Prebuilt CUDA GROMACS ${GROMACS_VERSION} installed successfully."
                CONDA_PREFIX="$prefix"
                return 0
            fi
            log_warn "Prebuilt package has $_backend backend instead of CUDA, will build from source."
        fi
    else
        log_warn "Prebuilt install failed or timed out, falling back to source build."
    fi
    return 1
}

build_gromacs() {
    log "Building GROMACS ${GROMACS_VERSION} with ${GPU_BACKEND} support..."

    local prefix=""

    if [ "$STANDALONE_MODE" = true ]; then
        prefix="$INSTALL_PREFIX"
        mkdir -p "$prefix"
    else
        eval "$(conda shell.bash hook)"
        safe_conda_activate "$ENV_NAME"
        prefix=$(conda info --base)/envs/$ENV_NAME
    fi

    # ---- Skip build if GROMACS is already installed with correct version AND backend ----
    local _gmx_bin=""
    [ -f "$prefix/bin/gmx" ] && _gmx_bin="$prefix/bin/gmx"
    [ -z "$_gmx_bin" ] && [ -f "$prefix/bin/gmx_mpi" ] && _gmx_bin="$prefix/bin/gmx_mpi"
    if [ -n "$_gmx_bin" ]; then
        local _installed_ver _installed_backend
        _installed_ver=$("$_gmx_bin" --version 2>&1 | sed -n 's/.*GROMACS version:[[:space:]]*\([0-9][0-9.]*\).*/\1/p' | head -1 || true)
        _installed_backend=$("$_gmx_bin" --version 2>&1 | awk '/GPU support:/{print $NF}')
        if [ "$_installed_ver" = "$GROMACS_VERSION" ] && [ "$_installed_backend" = "$GPU_BACKEND" ]; then
            log "GROMACS $GROMACS_VERSION ($GPU_BACKEND) already installed at $prefix, skipping build."
            CONDA_PREFIX="$prefix"
            return 0
        fi
        if [ "$_installed_ver" = "$GROMACS_VERSION" ] && [ "$_installed_backend" != "$GPU_BACKEND" ]; then
            log_warn "GROMACS $GROMACS_VERSION found but with $_installed_backend backend (need $GPU_BACKEND). Rebuilding..."
        fi
    fi

    # Try prebuilt conda-forge package first (much faster than source build)
    if try_conda_prebuilt; then
        return 0
    fi

    CONDA_PREFIX="$prefix"
    download_gromacs

    local build_func="build_gromacs_${GPU_BACKEND,,}"

    # Build 1: Non-MPI (gmx), then Build 2: MPI (gmx_mpi) if available
    if [ "$GPU_BACKEND" = "CUDA" ]; then
        log "=== Build 1: Non-MPI version (gmx) ==="
        $build_func "$prefix" "OFF" "cuda"

        # Build 2: MPI version (gmx_mpi) - if MPI is available
        if command -v mpirun &> /dev/null || [ -f "$prefix/bin/mpirun" ]; then
            log "=== Build 2: MPI version (gmx_mpi) ==="
            $build_func "$prefix" "ON" "cuda_mpi"
        else
            log_warn "MPI not found, skipping MPI build"
        fi
    else
        # HIP: build non-MPI first, then MPI if available
        log "=== Build 1: HIP Non-MPI version (gmx) ==="
        $build_func "$prefix" "OFF" "hip"

        if command -v mpirun &> /dev/null || [ -f "$prefix/bin/mpirun" ]; then
            log "=== Build 2: HIP + MPI version (gmx_mpi) ==="
            $build_func "$prefix" "ON" "hip_mpi"
        else
            log_warn "MPI not found, skipping HIP MPI build"
        fi
    fi

    log "GROMACS installed to: $prefix"
}

#------------------------------------------------------------------------------
# POST-INSTALL CONFIGURATION
#------------------------------------------------------------------------------

configure_environment() {
    log "Configuring environment..."

    if [ "$STANDALONE_MODE" = true ]; then
        configure_standalone_env
        return
    fi

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$ENV_NAME"

    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME

    ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
    DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"

    # ---- Skip if activation scripts already exist ----
    if [ -f "$ACTIVATE_DIR/gromacs_gpu.sh" ] && [ -f "$DEACTIVATE_DIR/gromacs_gpu.sh" ]; then
        log "Activation scripts already exist, skipping configuration."
        return 0
    fi

    mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

    if [ "$GPU_BACKEND" = "CUDA" ]; then
        configure_cuda_env "$ACTIVATE_DIR" "$DEACTIVATE_DIR"
    else
        configure_hip_env "$ACTIVATE_DIR" "$DEACTIVATE_DIR"
    fi

    log "Activation scripts created"
}

configure_cuda_env() {
    local activate_dir="$1"
    local deactivate_dir="$2"

    cat > "$activate_dir/gromacs_gpu.sh" << ACTIVATE_EOF
# GROMACS CUDA environment activation

# CUDA settings - check detected path first, then standard locations
for _cuda_candidate in "${CUDA_PATH:-}" "\$HOME/cuda-12.8" "/usr/local/cuda" "/usr/local/cuda-12"; do
    if [ -n "\$_cuda_candidate" ] && [ -d "\$_cuda_candidate" ]; then
        export CUDA_HOME="\${CUDA_HOME:-\$_cuda_candidate}"
        export PATH="\$CUDA_HOME/bin:\$PATH"
        export LD_LIBRARY_PATH="\$CUDA_HOME/lib64:\${LD_LIBRARY_PATH:-}"
        break
    fi
done

# GPU visibility
export CUDA_VISIBLE_DEVICES=0

# GROMACS settings
export GMX_MAXBACKUP=-1

# OpenMP settings - use half cores by default (thermal friendly)
TOTAL_CORES=\$(nproc 2>/dev/null || echo 8)
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-\$((TOTAL_CORES / 2))}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Source GROMACS environment
if [ -f "\$CONDA_PREFIX/bin/GMXRC.bash" ]; then
    source "\$CONDA_PREFIX/bin/GMXRC.bash"
fi

echo "GROMACS CUDA environment activated"
echo "  GPU: NVIDIA CUDA"
echo "  OMP threads: \$OMP_NUM_THREADS"
echo "  Use: gmx mdrun ... or gmx_mpi mdrun ..."
ACTIVATE_EOF

    cat > "$deactivate_dir/gromacs_gpu.sh" << 'DEACTIVATE_EOF'
unset GMX_MAXBACKUP
unset CUDA_VISIBLE_DEVICES
unset CUDA_HOME
unset OMP_NUM_THREADS
unset OMP_PLACES
unset OMP_PROC_BIND
DEACTIVATE_EOF

    chmod +x "$activate_dir/gromacs_gpu.sh" "$deactivate_dir/gromacs_gpu.sh"
}

configure_hip_env() {
    local activate_dir="$1"
    local deactivate_dir="$2"

    cat > "$activate_dir/gromacs_gpu.sh" << ACTIVATE_EOF
# GROMACS HIP environment activation

# ROCm settings
export ROCM_PATH="\${ROCM_PATH:-/opt/rocm}"
export PATH="\$ROCM_PATH/bin:\$PATH"
export LD_LIBRARY_PATH="\$ROCM_PATH/lib:\${LD_LIBRARY_PATH:-}"

# AMD GPU settings
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GPU_MAX_HW_QUEUES=8
export HIP_VISIBLE_DEVICES=0
export ROCR_VISIBLE_DEVICES=0

# DO NOT set HSA_OVERRIDE_GFX_VERSION - causes kernel mismatch!

# GROMACS settings
export GMX_MAXBACKUP=-1

# OpenMP settings - use half cores by default (thermal friendly)
TOTAL_CORES=\$(nproc 2>/dev/null || echo 8)
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-\$((TOTAL_CORES / 2))}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Source GROMACS environment
if [ -f "\$CONDA_PREFIX/bin/GMXRC.bash" ]; then
    source "\$CONDA_PREFIX/bin/GMXRC.bash"
fi

echo "GROMACS HIP environment activated"
echo "  GPU: AMD HIP"
echo "  Use: gmx_mpi or mpirun -np 1 gmx_mpi mdrun ..."
ACTIVATE_EOF

    cat > "$deactivate_dir/gromacs_gpu.sh" << 'DEACTIVATE_EOF'
unset GMX_ENABLE_DIRECT_GPU_COMM
unset GPU_MAX_HW_QUEUES
unset GMX_MAXBACKUP
unset HIP_VISIBLE_DEVICES
unset ROCR_VISIBLE_DEVICES
unset ROCM_PATH
unset OMP_NUM_THREADS
unset OMP_PLACES
unset OMP_PROC_BIND
DEACTIVATE_EOF

    chmod +x "$activate_dir/gromacs_gpu.sh" "$deactivate_dir/gromacs_gpu.sh"
}

configure_standalone_env() {
    log "Standalone mode — no conda activation scripts to create."
    log "Use 'source $INSTALL_PREFIX/bin/GMXRC' to activate GROMACS."
    create_module_file
}

#------------------------------------------------------------------------------
# MODULE FILE (standalone mode)
#------------------------------------------------------------------------------

create_module_file() {
    local module_dir="/etc/modulefiles/gromacs"
    local module_name="${GROMACS_VERSION}-${GPU_BACKEND,,}"

    # Skip if modulefiles dir doesn't exist or not writable
    if [ ! -d "/etc/modulefiles" ] || [ ! -w "/etc/modulefiles" ]; then
        log_warn "Cannot create module file (no /etc/modulefiles or not writable)"
        log "To use: source $INSTALL_PREFIX/bin/GMXRC"
        return
    fi

    log "Creating environment module file..."
    mkdir -p "$module_dir"

    cat > "$module_dir/$module_name" << EOF
#%Module1.0
proc ModulesHelp { } {
    puts stderr "GROMACS ${GROMACS_VERSION} with ${GPU_BACKEND} GPU support"
}

module-whatis "GROMACS ${GROMACS_VERSION} with ${GPU_BACKEND} support"

set GROMACS_HOME $INSTALL_PREFIX
prepend-path PATH \$GROMACS_HOME/bin
prepend-path LD_LIBRARY_PATH \$GROMACS_HOME/lib64
prepend-path MANPATH \$GROMACS_HOME/share/man
setenv GMXBIN \$GROMACS_HOME/bin
setenv GMXDATA \$GROMACS_HOME/share/gromacs
EOF

    if [ "$GPU_BACKEND" = "HIP" ]; then
        cat >> "$module_dir/$module_name" << EOF

# AMD GPU settings
setenv   ROCM_PATH          /opt/rocm
setenv   GMX_ENABLE_DIRECT_GPU_COMM 1
setenv   GPU_MAX_HW_QUEUES  8
EOF
    fi

    log "Module file created: $module_dir/$module_name"
    log "Use: module load gromacs/$module_name"
}

#------------------------------------------------------------------------------
# TEST INSTALLATION (standalone mode)
#------------------------------------------------------------------------------

test_standalone_installation() {
    log "Testing GROMACS installation..."

    # Source GROMACS
    if [ -f "$INSTALL_PREFIX/bin/GMXRC" ]; then
        source "$INSTALL_PREFIX/bin/GMXRC"
    fi

    # Check for binaries
    local gmx_bin=""
    if [ -f "$INSTALL_PREFIX/bin/gmx" ]; then
        gmx_bin="$INSTALL_PREFIX/bin/gmx"
    elif [ -f "$INSTALL_PREFIX/bin/gmx_mpi" ]; then
        gmx_bin="$INSTALL_PREFIX/bin/gmx_mpi"
    else
        log_error "No GROMACS binary found in $INSTALL_PREFIX/bin/"
    fi

    log "GROMACS version:"
    "$gmx_bin" --version 2>&1 | head -20

    log "Testing GPU detection..."
    "$gmx_bin" mdrun -version 2>&1 | grep -iE "gpu|cuda|hip|rocm" || true

    echo ""
    log "=============================================="
    log "Installation complete!"
    log "=============================================="
    echo ""
    log "To use this GROMACS installation, run:"
    log "  source $INSTALL_PREFIX/bin/GMXRC"
    echo ""
    log "Then use gmx (non-MPI) or gmx_mpi (MPI version)"
    echo ""
    log "Example mdrun with GPU:"
    if [ "$GPU_BACKEND" = "HIP" ]; then
        log "  gmx_mpi mdrun -v -deffnm md -nb gpu -pme gpu -bonded cpu"
    else
        log "  gmx mdrun -v -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu"
    fi
    echo ""
}

#------------------------------------------------------------------------------
# VERIFY INSTALLATION
#------------------------------------------------------------------------------

verify_installation() {
    if [ "$STANDALONE_MODE" = true ]; then
        test_standalone_installation
        return
    fi

    log "Verifying installation..."

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$ENV_NAME"

    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME

    # Check for gmx or gmx_mpi
    local gmx_bin=""
    if [ -f "$CONDA_PREFIX/bin/gmx" ]; then
        gmx_bin="$CONDA_PREFIX/bin/gmx"
        log "GROMACS binary: $gmx_bin"
    fi
    if [ -f "$CONDA_PREFIX/bin/gmx_mpi" ]; then
        log "GROMACS MPI binary: $CONDA_PREFIX/bin/gmx_mpi"
        [ -z "$gmx_bin" ] && gmx_bin="$CONDA_PREFIX/bin/gmx_mpi"
    fi

    if [ -z "$gmx_bin" ]; then
        log_error "No GROMACS binary found!"
    fi

    "$gmx_bin" --version 2>&1 | head -10

    # Check GPU detection
    log "Checking GPU detection..."
    "$gmx_bin" mdrun -version 2>&1 | grep -iE "gpu|cuda|hip|rocm" || true

    # Check Python imports
    log "Checking Python packages..."
    python3 -c "
import numpy, scipy, pandas, matplotlib, MDAnalysis
print('All Python packages OK')
" || log_error "Some Python packages failed to import"

    # Check video-rendering tools (Step 6 visualize_results)
    log "Checking video tools (VMD + ffmpeg)..."
    if command -v vmd &> /dev/null; then
        log "VMD: $(vmd -h 2>&1 | head -1 || echo 'available')"
    else
        log_warn "vmd not found in $ENV_NAME; Step 6 video rendering will be skipped."
        log_warn "  Install with: $PKG_MGR install -n $ENV_NAME -c conda-forge vmd"
    fi
    if command -v ffmpeg &> /dev/null; then
        log "ffmpeg: $(ffmpeg -version 2>&1 | head -1)"
    else
        log_warn "ffmpeg not found in $ENV_NAME; movie encoding will be skipped."
    fi

    log "Installation verified successfully!"

    echo ""
    echo "============================================================"
    echo " Setup Complete! (${GPU_BACKEND:-auto-detected})"
    echo "============================================================"
    echo ""
    echo "To use GROMACS:"
    echo "  conda activate $ENV_NAME"

    if [ "${GPU_BACKEND:-CUDA}" = "CUDA" ]; then
        echo ""
        echo "Run simulations with GPU:"
        echo "  gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu"
        echo ""
        echo "For large systems (limited VRAM):"
        echo "  gmx mdrun -deffnm md -nb gpu -pme cpu -bonded cpu"
        echo ""
        echo "Monitor GPU: nvidia-smi -l 5"
    else
        echo ""
        echo "Run simulations with GPU:"
        echo "  mpirun -np 1 gmx_mpi mdrun -deffnm md -nb gpu -pme gpu -bonded cpu"
        echo ""
        echo "HIP limitations:"
        echo "  - Use '-bonded cpu' (GPU bonded not implemented)"
        echo "  - Do NOT set HSA_OVERRIDE_GFX_VERSION"
    fi
    echo ""
}

#------------------------------------------------------------------------------
# MUTATEX (FoldX) SETUP
#------------------------------------------------------------------------------

install_mutatex_dependencies() {
    log "Setting up MutateX environment: $MUTATEX_ENV_NAME"

    # ---- Fast path: skip entirely if env exists and mutatex is importable ----
    if conda env list | grep -q "^${MUTATEX_ENV_NAME} "; then
        if conda run -n "$MUTATEX_ENV_NAME" bash -c 'command -v mutatex >/dev/null && python3 -c "import numpy, scipy, pandas, Bio"' 2>/dev/null; then
            log "MutateX and dependencies already installed in '$MUTATEX_ENV_NAME', skipping."
            return 0
        fi
        log "Environment '$MUTATEX_ENV_NAME' exists but mutatex missing/broken; will complete install."
    else
        log "Creating new environment: $MUTATEX_ENV_NAME"
        conda create -n "$MUTATEX_ENV_NAME" python="$MUTATEX_PYTHON_VERSION" -y
    fi

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$MUTATEX_ENV_NAME"

    # ---- All conda packages in one solver pass ----
    log "Installing conda dependencies..."
    $PKG_MGR install -y -c conda-forge \
        numpy scipy pandas matplotlib seaborn biopython networkx pillow tqdm pyyaml

    # ---- All pip packages in one invocation ----
    log "Installing MutateX and analysis tools..."
    pip install --quiet mutatex pdb-tools biotite

    log "MutateX dependencies installed to: $(conda info --base)/envs/$MUTATEX_ENV_NAME"
}

verify_mutatex() {
    log "Verifying MutateX installation..."

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$MUTATEX_ENV_NAME"

    # Check mutatex command
    if ! command -v mutatex &> /dev/null; then
        log_error "mutatex command not found in $MUTATEX_ENV_NAME environment"
    fi
    log "MutateX: $(mutatex --version 2>&1 || echo 'installed (no --version flag)')"

    # Check FoldX binary
    local foldx_binary
    foldx_binary=$(ls "$SCRIPT_DIR"/modules/PPI/foldx_* 2>/dev/null | head -1)
    if [ -n "$foldx_binary" ] && [ -f "$foldx_binary" ]; then
        log "FoldX binary: $foldx_binary (found)"
        if [ ! -x "$foldx_binary" ]; then
            log_warn "FoldX binary is not executable. Fixing permissions..."
            chmod +x "$foldx_binary"
        fi
    else
        log_warn "No FoldX binary found in $SCRIPT_DIR/modules/PPI/"
        log_warn "Download FoldX from https://foldxsuite.crg.eu/ and place it there."
    fi

    # Check rotabase
    local rotabase="$SCRIPT_DIR/modules/PPI/rotabase.txt"
    if [ -f "$rotabase" ]; then
        log "Rotabase: $rotabase (found)"
    else
        log_warn "rotabase.txt not found at $rotabase"
    fi

    # Check Python imports
    log "Checking Python packages..."
    python3 -c "
import numpy, scipy, pandas, matplotlib, seaborn
from Bio import PDB
print('All MutateX Python packages OK')
" || log_error "Some Python packages failed to import"

    log "MutateX installation verified successfully!"

    echo ""
    echo "============================================================"
    echo " MutateX Setup Complete!"
    echo "============================================================"
    echo ""
    echo "To use MutateX:"
    echo "  conda activate $MUTATEX_ENV_NAME"
    echo "  ./k_mutatex_repo_version.sh"
    echo ""
    echo "Configuration: config/PPI/mutatex/a_mutatex.toml"
    echo ""
}

#------------------------------------------------------------------------------
# GMX_MMPBSA SETUP (Stage 14 binding-energy mmpbsa backend)
#------------------------------------------------------------------------------
# Stage 14 reads [tools].gmx_mmpbsa_env from the TOML and invokes gmx_MMPBSA
# via `conda run -n <env>` inside modules/14_special_pipeline/compute_mmpbsa.sh.
# That env is currently UNUSED while md_equilibration is postponed in the
# config, but we install it ahead of time so re-enabling MD is a one-line TOML
# flip rather than a tooling step.
install_gmxmmpbsa_dependencies() {
    log "Setting up gmx_MMPBSA environment: $GMXMMPBSA_ENV_NAME"

    # ---- Fast path: skip entirely if env exists and gmx_MMPBSA is importable ----
    if conda env list | grep -q "^${GMXMMPBSA_ENV_NAME} "; then
        local _gmx_bin
        _gmx_bin="$(conda info --base)/envs/${GMXMMPBSA_ENV_NAME}/bin/gmx_MMPBSA"
        if [ -x "$_gmx_bin" ] && \
           conda run -n "$GMXMMPBSA_ENV_NAME" python3 -c "import GMXMMPBSA" 2>/dev/null; then
            log "gmx_MMPBSA already installed in '$GMXMMPBSA_ENV_NAME', skipping."
            return 0
        fi
        log "Environment '$GMXMMPBSA_ENV_NAME' exists but gmx_MMPBSA missing/broken; removing and recreating."
        conda env remove -n "$GMXMMPBSA_ENV_NAME" -y 2>/dev/null || true
    fi

    # gmx_mmpbsa pins pandas=1.2.2 which requires python<3.10.  Creating the
    # env with python=3.10 first, then installing gmx_mmpbsa, results in an
    # unsolvable conflict (the python=3.10 pin file blocks the downgrade).
    # Solution: create the env *with* gmx_mmpbsa in a single solver pass and
    # let conda pick a compatible python (typically 3.9).
    log "Creating env '$GMXMMPBSA_ENV_NAME' with gmx_MMPBSA + AmberTools (single solver pass)..."
    log "  NOTE: solving ambertools+gmx_mmpbsa typically takes 2-10 minutes; libmamba is silent during solve."
    # --override-channels avoids pulling defaults/anaconda.com (the commercial
    # warnings above) and speeds up the solve.  -v gives the user feedback that
    # the solver is actually running.
    if ! $PKG_MGR create -n "$GMXMMPBSA_ENV_NAME" -y -v \
            --override-channels -c conda-forge -c bioconda \
            --channel-priority flexible \
            "gmx_mmpbsa>=1.6" "ambertools>=23" parmed mpi4py \
            numpy scipy pandas matplotlib h5py 2>&1; then
        log_warn "Combined create+install failed. Trying minimal spec set..."
        # Fallback: minimal spec without secondary packages that may conflict.
        conda env remove -n "$GMXMMPBSA_ENV_NAME" -y 2>/dev/null || true
        $PKG_MGR create -n "$GMXMMPBSA_ENV_NAME" -y -v \
            --override-channels -c conda-forge -c bioconda \
            --channel-priority flexible \
            "gmx_mmpbsa" "ambertools" \
            || log_error "gmx_MMPBSA install failed. Manual: $PKG_MGR create -n $GMXMMPBSA_ENV_NAME --override-channels -c conda-forge -c bioconda gmx_mmpbsa ambertools"
    fi

    log "gmx_MMPBSA installed to: $(conda info --base)/envs/$GMXMMPBSA_ENV_NAME"
}

verify_gmxmmpbsa() {
    log "Verifying gmx_MMPBSA installation..."

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$GMXMMPBSA_ENV_NAME"

    if ! command -v gmx_MMPBSA &> /dev/null; then
        log_warn "gmx_MMPBSA CLI not found in $GMXMMPBSA_ENV_NAME."
        log_warn "  Reinstall with: $PKG_MGR install -n $GMXMMPBSA_ENV_NAME -c conda-forge -c bioconda gmx_mmpbsa"
        return 1
    fi
    log "gmx_MMPBSA: $(gmx_MMPBSA --version 2>&1 | head -1)"

    # AmberTools sanity check (tleap is needed for parameterisation)
    if command -v tleap &> /dev/null; then
        log "AmberTools tleap: $(tleap -help 2>&1 | head -1 || echo 'available')"
    else
        log_warn "tleap not found - AmberTools may be incomplete."
    fi

    log "gmx_MMPBSA installation verified."

    echo ""
    echo "============================================================"
    echo " gmx_MMPBSA Setup Complete!"
    echo "============================================================"
    echo ""
    echo "To use gmx_MMPBSA (Stage 14 binding_energy mmpbsa backend):"
    echo "  conda activate $GMXMMPBSA_ENV_NAME"
    echo "  gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ci index.ndx ..."
    echo ""
    echo "Pipeline integration: 14_Interaction_Domain_MappingCONFIG.toml"
    echo "  [tools].gmx_mmpbsa_env  = \"$GMXMMPBSA_ENV_NAME\""
    echo "  [binding_energy].backends includes \"mmpbsa\"  (currently gated"
    echo "                                                  off while GROMACS"
    echo "                                                  is postponed)"
    echo ""
}

#------------------------------------------------------------------------------
# STAGE 14 PYTHON DEPS (egg env)
#------------------------------------------------------------------------------
# Stage 14 calls interface_analysis (analyze_interface.py), binding_energy
# prodigy backend (compute_prodigy_dg.sh), and comparative_report
# (compare_variants.py) all via `conda run -n egg python3 ...`. Those scripts
# need PRODIGY (prodigy-prot), freesasa, gemmi (CIF parser), and Biopython.
# The egg env itself is built by setup_unified_conda_env.sh; here we only
# augment it with the Stage 14 additions if it already exists.
install_stage14_python_deps() {
    local target_env="${1:-$STAGE14_PRODIGY_ENV}"
    log "Augmenting '$target_env' env with Stage 14 Python deps..."

    if ! conda env list | grep -q "^${target_env} "; then
        log_warn "Env '$target_env' not found. Stage 14 expects [tools].prodigy_conda_env = \"$target_env\"."
        log_warn "  Create the egg env first: bash setup_unified_conda_env.sh"
        log_warn "  Then re-run: bash setup_gromacs_and_mutatex.sh --stage14"
        return 1
    fi

    # ---- Fast path: skip entirely if all imports succeed (no env activation needed) ----
    if conda run -n "$target_env" python3 -c "import freesasa, gemmi, prodigy_prot; from Bio import PDB" 2>/dev/null; then
        log "Stage 14 Python deps already present in '$target_env', skipping."
        return 0
    fi

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$target_env"

    log "Installing freesasa, gemmi, biopython via conda-forge..."
    $PKG_MGR install -y -c conda-forge biopython freesasa gemmi \
        || log_warn "conda install partially failed; will try pip fallback for missing packages."

    # PRODIGY (prodigy-prot) is PyPI-only. Keep --quiet so the install log
    # stays readable.
    if ! python3 -c "import prodigy_prot" 2>/dev/null; then
        log "Installing prodigy-prot (PyPI)..."
        pip install --quiet "prodigy-prot>=2.1" \
            || log_warn "prodigy-prot install failed. Install manually: pip install prodigy-prot"
    fi

    log "Stage 14 Python deps installed in '$target_env'."
}

verify_stage14_python_deps() {
    local target_env="${1:-$STAGE14_PRODIGY_ENV}"
    log "Verifying Stage 14 Python deps in '$target_env'..."

    if ! conda env list | grep -q "^${target_env} "; then
        log_warn "Env '$target_env' not found - skipping verification."
        return 1
    fi

    eval "$(conda shell.bash hook)"
    safe_conda_activate "$target_env"

    python3 -c "
import freesasa
import gemmi
from Bio import PDB
import prodigy_prot
print(f'  freesasa     {freesasa.__version__ if hasattr(freesasa, \"__version__\") else \"installed\"}')
print(f'  gemmi        {gemmi.__version__}')
print(f'  biopython    {PDB.__doc__.splitlines()[0] if PDB.__doc__ else \"installed\"}')
print(f'  prodigy_prot installed')
print('Stage 14 Python deps OK')
" || { log_warn "Some Stage 14 Python imports failed - re-run install_stage14_python_deps."; return 1; }
}

#------------------------------------------------------------------------------
# FOLDX SCAFFOLD (Stage 14: tools/foldx/foldx_<date>)
#------------------------------------------------------------------------------
# FoldX is closed-source; we can only scaffold the drop-zone and remind the
# user. Stage 14's manual step [M4] (printed by `show_manual = true`) covers
# the same ground; this just creates the directory so the path the TOML
# points at actually exists.
scaffold_foldx_dir() {
    local foldx_dir="$SCRIPT_DIR/tools/foldx"
    log "Scaffolding FoldX drop-zone: $foldx_dir"
    mkdir -p "$foldx_dir"

    # Skip if a foldx binary is already present
    if ls "$foldx_dir"/foldx_* 2>/dev/null | grep -qv README; then
        local _bin
        _bin=$(ls "$foldx_dir"/foldx_* 2>/dev/null | grep -v README | head -1)
        log "FoldX binary already present: $_bin"
        [ ! -x "$_bin" ] && { log_warn "Not executable, fixing..."; chmod +x "$_bin"; }
        return 0
    fi

    if [ ! -f "$foldx_dir/README.txt" ]; then
        cat > "$foldx_dir/README.txt" << 'FOLDX_README'
FoldX binary drop-zone (Stage 14)
==================================

Stage 14 (14_Interaction_Domain_Mapping.sh) uses FoldX for:
  - binding_energy.foldx backend (compute_foldx_dg.sh)
  - alanine_scan operation (alanine_scan.sh)

FoldX is closed-source; this script cannot download it for you.

Steps:
  1. Request a free academic license at https://foldxsuite.crg.eu/
  2. Download the Linux binary (e.g. foldx_20251231)
  3. Place it in this directory:  tools/foldx/
  4. chmod +x foldx_<date>
  5. Update the TOML:
        [tools]
        foldx_binary = "tools/foldx/foldx_20251231"
     in 14_Interaction_Domain_MappingCONFIG.toml

Also drop rotabase.txt here if your FoldX version needs it
(older versions did; newer 5.x packages embed it).
FOLDX_README
    fi

    log_warn "FoldX binary not present. See $foldx_dir/README.txt for instructions."
    log_warn "Without FoldX, Stage 14 will skip the foldx binding-energy backend AND the alanine_scan operation."
}

#------------------------------------------------------------------------------
# STAGE 14 STATUS BANNER
#------------------------------------------------------------------------------
print_stage14_status() {
    local foldx_dir="$SCRIPT_DIR/tools/foldx"
    local foldx_bin
    foldx_bin=$(ls "$foldx_dir"/foldx_* 2>/dev/null | grep -v README | head -1 || true)
    local foldx_state="MISSING (manual download)"
    [ -n "$foldx_bin" ] && [ -x "$foldx_bin" ] && foldx_state="OK ($foldx_bin)"

    echo ""
    echo "============================================================"
    echo " Stage 14 Setup Status (14_Interaction_Domain_Mapping.sh)"
    echo "============================================================"
    printf "  GROMACS env (%s)    : %s\n" "$ENV_NAME" \
        "$(conda env list | grep -q "^${ENV_NAME} " && echo "OK" || echo "NOT INSTALLED")"
    printf "  MutaTeX env (%s)               : %s\n" "$MUTATEX_ENV_NAME" \
        "$(conda env list | grep -q "^${MUTATEX_ENV_NAME} " && echo "OK" || echo "NOT INSTALLED")"
    printf "  gmx_MMPBSA env (%s)     : %s\n" "$GMXMMPBSA_ENV_NAME" \
        "$(conda env list | grep -q "^${GMXMMPBSA_ENV_NAME} " && echo "OK" || echo "NOT INSTALLED")"
    printf "  Stage 14 Python deps in '%s' : " "$STAGE14_PRODIGY_ENV"
    if conda env list | grep -q "^${STAGE14_PRODIGY_ENV} "; then
        if conda run -n "$STAGE14_PRODIGY_ENV" python3 -c "import freesasa, gemmi, prodigy_prot" 2>/dev/null; then
            echo "OK"
        else
            echo "PARTIAL (re-run install_stage14_python_deps)"
        fi
    else
        echo "egg env not found (run setup_unified_conda_env.sh)"
    fi
    printf "  FoldX binary (tools/foldx/)  : %s\n" "$foldx_state"
    echo ""
    echo "Manual next steps (set [run].show_manual = true in the TOML for the full panel):"
    echo "  [M3] Submit AF3 complexes at https://alphafoldserver.com/"
    echo "  [M4] Download FoldX from https://foldxsuite.crg.eu/ -> tools/foldx/"
    echo "  [M5] Confirm \`gmx --version\` reports 'GPU support: CUDA' (when MD is re-enabled)"
    echo "  [M6] (optional) CHARMM-GUI membrane prep if [md_equilibration].membrane = true"
    echo ""
    echo "Note: 14_Interaction_Domain_MappingCONFIG.toml currently has md_equilibration"
    echo "      commented out and the mmpbsa backend disabled. The GROMACS + gmxmmpbsa"
    echo "      envs are installed ahead of time so re-enabling MD is a TOML flip."
    echo ""
}

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------

main() {
    setup_installation_logging
    trap teardown_installation_logging EXIT

    echo ""
    echo "============================================================"
    echo " GROMACS GPU Environment Setup (Auto-Detect)"
    echo "============================================================"
    echo ""

    # Single-pass argument parsing (fixes prev-variable scoping bug in 3-loop approach)
    local action=""
    local _args=("$@")
    local _i=0
    while [[ $_i -lt ${#_args[@]} ]]; do
        local _arg="${_args[$_i]}"
        case "$_arg" in
            --force-cuda)    FORCE_BACKEND="CUDA"   ;;
            --force-hip)     FORCE_BACKEND="HIP"    ;;
            --mutatex)       SETUP_MUTATEX=true     ;;
            --gmxmmpbsa)     SETUP_GMXMMPBSA=true   ;;
            --stage14|--stage-14)
                # Full Stage 14 surface: implies GROMACS + MutaTeX + gmxmmpbsa
                # + egg-env Python deps + FoldX scaffold. Re-enables MD-side
                # envs ahead of time even though md_equilibration is currently
                # commented out in 14_Interaction_Domain_MappingCONFIG.toml.
                SETUP_STAGE14=true
                SETUP_MUTATEX=true
                SETUP_GMXMMPBSA=true
                ;;
            --standalone)
                STANDALONE_MODE=true
                local _next="${_args[$((_i+1))]:-}"
                if [[ -n "$_next" && "$_next" != --* ]]; then
                    INSTALL_PREFIX="$_next"
                    _i=$(( _i + 1 ))
                fi
                ;;
            --deps-only|--build-only|--verify|--mutatex-only|--gmxmmpbsa-only|--stage14-only|--stage-14-only)
                action="$_arg"
                ;;
            *)
                [[ "$_arg" == --* ]] && log_warn "Unknown flag: $_arg"
                ;;
        esac
        _i=$(( _i + 1 ))
    done

    # Handle --mutatex-only (standalone mutatex setup, no GROMACS)
    if [[ "${action:-}" == "--mutatex-only" ]]; then
        check_conda
        install_mutatex_dependencies
        verify_mutatex
        return
    fi

    # Handle --gmxmmpbsa-only (standalone gmx_MMPBSA env setup, no GROMACS)
    if [[ "${action:-}" == "--gmxmmpbsa-only" ]]; then
        check_conda
        install_gmxmmpbsa_dependencies
        verify_gmxmmpbsa
        return
    fi

    # Handle --stage14-only: install all Stage 14 tooling EXCEPT the GROMACS
    # build. Use this when MD is postponed (md_equilibration commented out in
    # the TOML) but you want the other envs (MutaTeX, gmxmmpbsa, PRODIGY/
    # freesasa in egg) and the FoldX scaffold ready to go.
    if [[ "${action:-}" == "--stage14-only" || "${action:-}" == "--stage-14-only" ]]; then
        check_conda
        install_mutatex_dependencies
        verify_mutatex
        install_gmxmmpbsa_dependencies
        verify_gmxmmpbsa
        install_stage14_python_deps "$STAGE14_PRODIGY_ENV" || true
        verify_stage14_python_deps "$STAGE14_PRODIGY_ENV"  || true
        scaffold_foldx_dir
        # Detect GPU just for the status banner ENV_NAME; tolerate failure
        # (no GPU is fine for the python/foldx side).
        detect_gpu 2>/dev/null || ENV_NAME="${ENV_NAME:-(not built)}"
        print_stage14_status
        return
    fi

    case "${action:-}" in
        --deps-only)
            detect_gpu
            if [ "$STANDALONE_MODE" = true ]; then
                check_prerequisites
                log "Standalone mode: no conda deps to install. Prerequisites OK."
            else
                check_conda
                install_dependencies
                if [ "$SETUP_MUTATEX" = true ]; then
                    install_mutatex_dependencies
                fi
                if [ "$SETUP_GMXMMPBSA" = true ]; then
                    install_gmxmmpbsa_dependencies
                fi
                if [ "$SETUP_STAGE14" = true ]; then
                    install_stage14_python_deps "$STAGE14_PRODIGY_ENV" || true
                    scaffold_foldx_dir
                fi
            fi
            ;;
        --build-only)
            detect_gpu
            if [ "$STANDALONE_MODE" = true ]; then
                check_prerequisites
                [ -z "$INSTALL_PREFIX" ] && INSTALL_PREFIX="${HOME}/.local/opt/gromacs-${GROMACS_VERSION}"
            else
                check_conda
            fi
            build_gromacs
            configure_environment
            verify_installation
            ;;
        --verify)
            detect_gpu
            if [ "$STANDALONE_MODE" = true ]; then
                [ -z "$INSTALL_PREFIX" ] && INSTALL_PREFIX="${HOME}/.local/opt/gromacs-${GROMACS_VERSION}"
            fi
            verify_installation
            if [ "$SETUP_MUTATEX" = true ]; then
                verify_mutatex
            fi
            if [ "$SETUP_GMXMMPBSA" = true ]; then
                verify_gmxmmpbsa
            fi
            if [ "$SETUP_STAGE14" = true ]; then
                verify_stage14_python_deps "$STAGE14_PRODIGY_ENV" || true
                print_stage14_status
            fi
            ;;
        *)
            detect_gpu
            if [ "$STANDALONE_MODE" = true ]; then
                [ -z "$INSTALL_PREFIX" ] && INSTALL_PREFIX="${HOME}/.local/opt/gromacs-${GROMACS_VERSION}"
                check_prerequisites

                echo ""
                log "Standalone build to: $INSTALL_PREFIX"
                read -r -n 1 -p "Continue with build? [y/N] " REPLY || true
                echo ""
                [[ ! $REPLY =~ ^[Yy]$ ]] && { log "Aborted"; exit 0; }

                build_gromacs
                configure_environment
                verify_installation
            else
                check_conda
                install_dependencies
                build_gromacs
                configure_environment
                verify_installation
                if [ "$SETUP_MUTATEX" = true ]; then
                    install_mutatex_dependencies
                    verify_mutatex
                fi
                if [ "$SETUP_GMXMMPBSA" = true ]; then
                    install_gmxmmpbsa_dependencies
                    verify_gmxmmpbsa
                fi
                if [ "$SETUP_STAGE14" = true ]; then
                    install_stage14_python_deps "$STAGE14_PRODIGY_ENV" || true
                    verify_stage14_python_deps "$STAGE14_PRODIGY_ENV"  || true
                    scaffold_foldx_dir
                    print_stage14_status
                fi
            fi
            ;;
    esac
}

main "$@"
