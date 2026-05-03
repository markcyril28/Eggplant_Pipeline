#!/usr/bin/env bash
set -euo pipefail
#######################################################################
# GROMACS Common Functions and Configuration
#
# Source this file in other GROMACS scripts:
#   source "$(dirname "$0")/modules/gromacs_common.sh"
#
# Provides:
# - Common environment variables and paths
# - GPU configuration for NVIDIA CUDA
# - Logging functions
# - Structure preparation functions
# - GROMACS wrapper functions
#######################################################################

#------------------------------------------------------------------------------
# PATH CONFIGURATION
#------------------------------------------------------------------------------

# Detect script directory (where this file lives)
if [[ -z "${GROMACS_COMMON_DIR:-}" ]]; then
    GROMACS_COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Project root (parent of modules)
if [[ -z "${PROJECT_ROOT:-}" ]]; then
    PROJECT_ROOT="$(dirname "$GROMACS_COMMON_DIR")"
fi

# Modules directory
MODULES_DIR="${GROMACS_COMMON_DIR}/gromacs_utils"

# Add modules to Python path for CLI access
export PYTHONPATH="${GROMACS_COMMON_DIR}:${PYTHONPATH:-}"

#------------------------------------------------------------------------------
# DEFAULT GROMACS CONFIGURATION
#------------------------------------------------------------------------------

# Fix MPI library conflicts (use conda's OpenMPI, not system)
# Prefer the gromacs_CUDA environment for GROMACS binaries
GROMACS_CUDA_ENV="${HOME}/miniconda3/envs/gromacs_CUDA"

# Query GPU backend of a given gmx binary (CUDA, OpenCL, SYCL, disabled, or "")
_gmx_gpu_backend() {
    "$1" --version 2>/dev/null | awk '/GPU support:/{print $NF}'
}

# Determine the best GROMACS binary — prefer CUDA backend for GPU offload.
# Search order: gromacs_CUDA env → current conda env → system gmx.
# Among candidates, CUDA wins over OpenCL (OpenCL is broken on WSL2).
if [[ -z "${GMX_BIN:-}" ]]; then
    _best_gmx=""
    _best_priority=0  # 0=none, 1=no-gpu, 2=OpenCL, 3=CUDA

    _candidates=()
    [[ -x "${GROMACS_CUDA_ENV}/bin/gmx" ]] && _candidates+=("${GROMACS_CUDA_ENV}/bin/gmx")
    [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/gmx" ]] && _candidates+=("${CONDA_PREFIX}/bin/gmx")
    command -v gmx &>/dev/null && _candidates+=("$(command -v gmx)")

    for _cand in "${_candidates[@]}"; do
        _backend=$(_gmx_gpu_backend "$_cand")
        case "$_backend" in
            CUDA)  _pri=3 ;;
            OpenCL|SYCL) _pri=2 ;;
            *)     _pri=1 ;;
        esac
        if (( _pri > _best_priority )); then
            _best_gmx="$_cand"
            _best_priority=$_pri
        fi
        # CUDA is the best possible — stop early
        (( _best_priority == 3 )) && break
    done

    if [[ -n "$_best_gmx" ]]; then
        GMX_BIN="$_best_gmx"
        GMX_CONDA_PATH="$(dirname "$(dirname "$_best_gmx")")"
    else
        GMX_CONDA_PATH="${CONDA_PREFIX:-${GROMACS_CUDA_ENV}}"
    fi
    unset _best_gmx _best_priority _candidates _cand _backend _pri
else
    # GMX_BIN already set externally — derive conda path from it
    GMX_CONDA_PATH="$(dirname "$(dirname "${GMX_BIN}")" 2>/dev/null)" || GMX_CONDA_PATH="${GROMACS_CUDA_ENV}"
fi

# CRITICAL: Prepend conda lib path to avoid system MPI conflicts
# This fixes: "undefined symbol: ompi_mpi_errors_throw_exceptions"
export LD_LIBRARY_PATH="${GMX_CONDA_PATH}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Use non-MPI gmx for single-node runs (avoids MPI library conflicts).
# GMX_BIN is set by the candidate loop above; this block is only a fallback.
if [[ -z "${GMX_BIN:-}" ]]; then
    if [[ -x "${GMX_CONDA_PATH}/bin/gmx" ]]; then
        GMX_BIN="${GMX_CONDA_PATH}/bin/gmx"
    elif command -v gmx &>/dev/null; then
        GMX_BIN="$(command -v gmx)"
    else
        GMX_BIN="${GMX_CONDA_PATH}/bin/gmx"  # Will fail later with clear error
    fi
fi

# Force field and water model
FORCEFIELD="${FORCEFIELD:-amber99sb-ildn}"
WATERMODEL="${WATERMODEL:-tip3p}"

# Box and solvation
BOX_DISTANCE="${BOX_DISTANCE:-1.5}"
BOX_TYPE="${BOX_TYPE:-dodecahedron}"
ION_CONCENTRATION="${ION_CONCENTRATION:-0.15}"

# Simulation parameters
EM_STEPS="${EM_STEPS:-50000}"
NVT_STEPS="${NVT_STEPS:-100000}"
NPT_STEPS="${NPT_STEPS:-100000}"
MD_STEPS="${MD_STEPS:-5000000}"

#------------------------------------------------------------------------------
# GPU CONFIGURATION (NVIDIA CUDA)
#------------------------------------------------------------------------------

# Number of threads — auto-detect available cores; callers may override via NTHREADS env var
NTHREADS="${NTHREADS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)}"
# Synchronize OMP_NUM_THREADS with NTHREADS to prevent mismatch
# (conda env may set OMP_NUM_THREADS to a different value than -ntomp)
export OMP_NUM_THREADS="$NTHREADS"

# GPU flags for different simulation types
# EM: Can't use -pme gpu or -update gpu with steep integrator; -bonded gpu OK on NVIDIA
GPU_EM_FLAGS="${GPU_EM_FLAGS:--nb gpu -pme cpu -bonded gpu}"
# MD: Maximum GPU offload — all forces on GPU
GPU_MD_FLAGS="${GPU_MD_FLAGS:--nb gpu -pme gpu -bonded gpu -update gpu}"

# Environment variables for NVIDIA CUDA GPU performance
setup_nvidia_gpu_env() {
    export GMX_ENABLE_DIRECT_GPU_COMM=1
    export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
}

# Call setup by default
setup_nvidia_gpu_env

#------------------------------------------------------------------------------
# LOGGING FUNCTIONS (skipped when logging_utils.sh is already sourced)
#------------------------------------------------------------------------------

if [[ "${LOGGING_UTILS_SOURCED:-}" != "true" ]]; then

# Simple log with timestamp
log() {
    echo "[$(date '+%H:%M:%S')] $1"
}

log_info() {
    echo "[$(date '+%H:%M:%S')] INFO: $1"
}

log_warn() {
    echo "[$(date '+%H:%M:%S')] WARN: $1" >&2
}

log_error() {
    echo "[$(date '+%H:%M:%S')] ERROR: $1" >&2
}

log_section() {
    echo ""
    echo "================================================================"
    echo " $1"
    echo "================================================================"
}

fi  # LOGGING_UTILS_SOURCED guard

#------------------------------------------------------------------------------
# REQUIREMENT CHECKS
#------------------------------------------------------------------------------

check_gromacs() {
    if [[ ! -x "$GMX_BIN" ]]; then
        log_error "GROMACS not found at: $GMX_BIN"
        log_error "Please activate the gromacs_CUDA conda environment: conda activate gromacs_CUDA"
        return 1
    fi
    
    # Verify GROMACS can actually run (check for library issues)
    local version_output
    version_output=$("$GMX_BIN" --version 2>&1) || {
        log_error "GROMACS failed to run: $version_output"
        log_error "This may be a library path issue. Try: conda activate gromacs_CUDA"
        return 1
    }
    
    local gpu_backend
    gpu_backend=$(_gmx_gpu_backend "$GMX_BIN")
    log_info "GROMACS: $(echo "$version_output" | head -1) [GPU: ${gpu_backend:-none}]"
    return 0
}

check_file() {
    local file="$1"
    local desc="${2:-File}"
    if [[ ! -f "$file" ]]; then
        log_error "$desc not found: $file"
        return 1
    fi
    return 0
}

check_directory() {
    local dir="$1"
    local desc="${2:-Directory}"
    if [[ ! -d "$dir" ]]; then
        log_error "$desc not found: $dir"
        return 1
    fi
    return 0
}

check_python_modules() {
    if [[ ! -d "$MODULES_DIR" ]]; then
        log_error "Python modules not found at: $MODULES_DIR"
        return 1
    fi
    return 0
}

#------------------------------------------------------------------------------
# STRUCTURE PREPARATION
#------------------------------------------------------------------------------

# Prepare structure using Python CLI (handles CIF/PDB, cleaning)
prepare_structure() {
    local input="$1"
    local output="${2:-clean.pdb}"
    
    python3 -m gromacs_utils.cli prepare-structure "$input" -o "$output"
}

# Quick PDB cleaning (bash-only, for simple cases)
clean_pdb_quick() {
    local input="$1"
    local output="$2"
    
    grep -E "^(ATOM|TER|END)" "$input" | \
        sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > "$output"
}

# Convert CIF to PDB using Python CLI
convert_cif_to_pdb() {
    local cif_file="$1"
    local pdb_file="$2"
    
    python3 -m gromacs_utils.cli prepare-structure "$cif_file" -o "$pdb_file"
}

#------------------------------------------------------------------------------
# MDP FILE GENERATION
#------------------------------------------------------------------------------

# Generate MDP files using Python CLI
generate_mdp() {
    local mdp_type="$1"
    local output="$2"
    shift 2
    
    python3 -m gromacs_utils.cli generate-mdp "$mdp_type" -o "$output" "$@"
}

# Generate all MDP files to a directory
generate_all_mdp() {
    local output_dir="$1"
    local em_steps="${2:-$EM_STEPS}"
    local nvt_steps="${3:-$NVT_STEPS}"
    local npt_steps="${4:-$NPT_STEPS}"
    local md_steps="${5:-$MD_STEPS}"
    
    python3 -m gromacs_utils.cli generate-mdp all -o "$output_dir" \
        --em-steps "$em_steps" \
        --nvt-steps "$nvt_steps" \
        --npt-steps "$npt_steps" \
        --md-steps "$md_steps"
}

#------------------------------------------------------------------------------
# GROMACS WRAPPER FUNCTIONS
#------------------------------------------------------------------------------

# Detect whether a usable GPU is present (cached for the session)
_detect_gpu_available() {
    # Return cached result if already checked
    if [[ -n "${_GPU_AVAILABLE_CACHED:-}" ]]; then
        [[ "$_GPU_AVAILABLE_CACHED" == "true" ]]; return $?
    fi
    _GPU_AVAILABLE_CACHED="false"

    # 1. Check GPU hardware
    local hw_found=false
    if command -v nvidia-smi &>/dev/null; then
        local count
        count=$(nvidia-smi --query-gpu=count --format=csv,noheader 2>/dev/null | head -1) || count=0
        [[ "${count:-0}" -gt 0 ]] && hw_found=true
    fi
    if [[ "$hw_found" != "true" ]] && command -v rocm-smi &>/dev/null; then
        rocm-smi --showproductname 2>/dev/null | grep -qi "GPU" && hw_found=true
    fi
    if [[ "$hw_found" != "true" ]]; then
        log_warn "No GPU hardware detected"
        return 1
    fi

    # 2. Verify GROMACS GPU backend runtime is functional
    local gpu_type
    gpu_type=$(_gmx_gpu_backend "${GMX_BIN:-gmx}")
    case "$gpu_type" in
        OpenCL)
            # OpenCL ICD may exist but the vendor library (libnvidia-opencl.so)
            # is absent in WSL2.  Check ldconfig for the actual implementation.
            if ! ldconfig -p 2>/dev/null | grep -q "libnvidia-opencl"; then
                log_warn "GROMACS compiled with OpenCL but libnvidia-opencl.so not found — forcing CPU mode"
                log_warn "Fix: install CUDA-compiled GROMACS:  mamba install -n gromacs_CUDA -c conda-forge gromacs=2026.0=nompi_cuda*"
                log_warn "Then re-run, or set GMX_BIN to the CUDA binary:  export GMX_BIN=\$HOME/miniconda3/envs/gromacs_CUDA/bin/gmx"
                return 1
            fi
            ;;
        CUDA)
            if ! ldconfig -p 2>/dev/null | grep -q "libcuda\.so"; then
                log_warn "GROMACS compiled with CUDA but libcuda.so not found — forcing CPU mode"
                log_warn "Ensure the NVIDIA driver is installed and /usr/lib/wsl/lib is in LD_LIBRARY_PATH"
                return 1
            fi
            ;;
        disabled|none|"")
            log_warn "GROMACS at ${GMX_BIN:-gmx} built without GPU support — forcing CPU mode"
            log_warn "Fix: install CUDA-compiled GROMACS:  mamba install -n gromacs_CUDA -c conda-forge gromacs=2026.0=nompi_cuda*"
            return 1
            ;;
    esac

    _GPU_AVAILABLE_CACHED="true"
    log_info "GPU backend: $gpu_type — GPU acceleration enabled"
    return 0
}

# Run mdrun with GPU, fall back to CPU if needed
run_mdrun() {
    local name="$1"
    local gpu_flags="$2"
    local log_dir="${3:-.}"

    mkdir -p "$log_dir"

    # Sync OMP_NUM_THREADS with current NTHREADS to avoid GROMACS mismatch error
    export OMP_NUM_THREADS="$NTHREADS"
    
    # Determine if GPU should be attempted:
    # 1. Honour explicit USE_GPU=false (set by callers e.g. _check_gpu)
    # 2. Otherwise auto-detect GPU hardware
    local try_gpu=true
    if [[ "${USE_GPU:-}" == "false" ]]; then
        try_gpu=false
    elif ! _detect_gpu_available; then
        try_gpu=false
    fi
    
    # `-pin off` is mandatory for parallel-dispatch callers: the default
    # `-pin auto` makes every concurrent mdrun pin to cores 0..NTHREADS-1,
    # so N parallel jobs all collide on the same cores and stall.
    # OS scheduler handles balancing; manual -pinoffset would also work but
    # would require slot-index tracking that the dispatcher does not expose.
    if [[ "$try_gpu" == true && -n "${GPU_ID:-}" ]]; then
        local gpu_id_flag="-gpu_id $GPU_ID"
        # GROMACS 2025+ requires -ntmpi when using GPU with OpenMP threads
        # gpu_flags must NOT be quoted — it contains multiple separate flags
        if $GMX_BIN mdrun -v -deffnm "$name" -ntmpi 1 -ntomp "$NTHREADS" -pin off $gpu_id_flag $gpu_flags 2>&1 | tee "${log_dir}/mdrun_${name}.log"; then
            return 0
        else
            log_warn "GPU failed, using CPU..."
        fi
    fi

    # CPU-only run
    $GMX_BIN mdrun -v -deffnm "$name" -ntmpi 1 -ntomp "$NTHREADS" -pin off 2>&1 | tee "${log_dir}/mdrun_${name}_cpu.log"
}

# Run energy minimization
run_em() {
    local name="${1:-em}"
    local log_dir="${2:-.}"
    
    run_mdrun "$name" "$GPU_EM_FLAGS" "$log_dir"
}

# Run MD simulation (NVT, NPT, or production)
run_md() {
    local name="$1"
    local log_dir="${2:-.}"
    
    run_mdrun "$name" "$GPU_MD_FLAGS" "$log_dir"
}

#------------------------------------------------------------------------------
# CHAIN HANDLING
#------------------------------------------------------------------------------

# Create chain index using Python CLI
create_chain_index() {
    local pdb="$1"
    local gro="$2"
    local index="$3"
    
    python3 -m gromacs_utils.cli chain-index --pdb "$pdb" --gro "$gro" --index "$index"
}

# Get chain info
get_chain_info() {
    local pdb="$1"
    local outdir="${2:-.}"
    
    python3 -m gromacs_utils.cli chain-info "$pdb" --outdir "$outdir"
}

#------------------------------------------------------------------------------
# ANALYSIS HELPERS
#------------------------------------------------------------------------------

# Extract metrics from GROMACS outputs
extract_metrics() {
    local workdir="$1"
    local output="${2:-metrics.txt}"
    
    python3 -m gromacs_utils.cli extract-metrics --workdir "$workdir" -o "$output"
}

# Generate plots from analysis
generate_plots() {
    local analysis_dir="$1"
    local plots_dir="$2"
    
    python3 -m gromacs_utils.cli generate-plots md -i "$analysis_dir" -o "$plots_dir"
}

# Generate gnuplot scripts and run if available
run_gnuplot_scripts() {
    local plots_dir="$1"
    
    if command -v gnuplot &> /dev/null; then
        for gp in "$plots_dir"/*.gp; do
            [[ -f "$gp" ]] && gnuplot "$gp" 2>/dev/null || true
        done
    fi
}

#------------------------------------------------------------------------------
# OUTPUT ORGANIZATION
#------------------------------------------------------------------------------

# Create standard output structure
setup_output_dirs() {
    local base_dir="$1"
    
    mkdir -p "${base_dir}"/{logs,structures,analysis,trajectories,statistics,visualization,plots}
}

# Setup output structure using Python CLI
setup_output_structure() {
    local workdir="$1"
    
    python3 -m gromacs_utils.cli setup-output --workdir "$workdir"
}

#------------------------------------------------------------------------------
# VISUALIZATION
#------------------------------------------------------------------------------

# Generate visualization scripts
generate_visualization() {
    local workdir="$1"
    local viz_type="${2:-all}"
    
    python3 -m gromacs_utils.visualization_generator --workdir "$workdir" --type "$viz_type"
}

# Run batch analysis
run_batch_analysis() {
    local workdir="$1"
    local with_plots="${2:-true}"
    
    if [[ "$with_plots" == "true" ]]; then
        python3 -m gromacs_utils.cli batch-analysis --workdir "$workdir" --plots
    else
        python3 -m gromacs_utils.cli batch-analysis --workdir "$workdir"
    fi
}

#------------------------------------------------------------------------------
# MODULE INFO
#------------------------------------------------------------------------------

gromacs_common_info() {
    echo "GROMACS Common Configuration"
    echo "============================="
    echo "Project Root: $PROJECT_ROOT"
    echo "Modules Dir:  $MODULES_DIR"
    echo "GMX Binary:   $GMX_BIN"
    echo "Force Field:  $FORCEFIELD"
    echo "Water Model:  $WATERMODEL"
    echo "GPU EM Flags: $GPU_EM_FLAGS"
    echo "GPU MD Flags: $GPU_MD_FLAGS"
    echo "Threads:      $NTHREADS"
}

# Export for use in subshells
export GMX_BIN FORCEFIELD WATERMODEL BOX_DISTANCE BOX_TYPE ION_CONCENTRATION
export EM_STEPS NVT_STEPS NPT_STEPS MD_STEPS
export NTHREADS GPU_EM_FLAGS GPU_MD_FLAGS
export PROJECT_ROOT MODULES_DIR
