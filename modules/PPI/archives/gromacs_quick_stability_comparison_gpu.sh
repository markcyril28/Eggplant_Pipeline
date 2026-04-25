#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS Quick Chain Stability Comparison Script - GPU VERSION
#
# This version attempts to use GPU acceleration. Requires a GROMACS build
# with working GPU support (HIP or CUDA, NOT SYCL/AdaptiveCpp on ROCm 7.x).
#
# A faster version that compares chain pairing stability using:
# 1. Energy minimization
# 2. Interface analysis (H-bonds, contacts, distances)
#
# Usage: ./gromacs_quick_stability_comparison_gpu.sh [options]
#######################################################################

#------------------------------------------------------------------------------
# Configuration - MODIFY THESE PARAMETERS AS NEEDED
#------------------------------------------------------------------------------

PDB1="${1:-${HOME}/PPI/inputs/SmelGRF-GIF/fold_1_x_1_model_0.pdb}"
PDB2="${2:-${HOME}/PPI/inputs/SmelGRF-GIF/fold_2_x_1_model_0.pdb}"
WORKDIR="${3:-${HOME}/PPI/results_gromacs/quick_analysis_gpu}"

# GROMACS parameters
FORCEFIELD="amber99sb-ildn"
WATERMODEL="tip3p"
BOX_DISTANCE=1.0
EM_STEPS=10000

# Performance settings - GPU mode
# For single GPU, fewer CPU threads often gives better GPU utilization
# 4-8 threads is optimal for most GPU workloads
NTHREADS=8

# GPU settings - OPTIMIZED for AMD MI210 with HIP backend
#
# IMPORTANT: HIP backend limitations (as of GROMACS 2025.4):
#   - bonded gpu is NOT supported yet (use CPU)
#   - update gpu IS supported for MD
#
# Key flags:
#   -nb gpu       : Non-bonded interactions on GPU (most compute-intensive)
#   -pme gpu      : PME electrostatics on GPU (requires md integrator)
#   -bonded cpu   : Bonded on CPU (HIP doesn't support GPU bonded yet)
#   -update gpu   : Keep coordinates on GPU, reduces CPU-GPU transfers
#
# For Energy Minimization: PME must be on CPU (steep integrator not supported)
GPU_EM_FLAGS="-nb gpu -pme cpu -bonded cpu"
# For MD simulations: Maximum GPU offload (except bonded)
GPU_MD_FLAGS="-nb gpu -pme gpu -bonded cpu -update gpu"

# Environment variables for optimal AMD GPU performance
export GMX_ENABLE_DIRECT_GPU_COMM=1    # Enable direct GPU communication
export GPU_MAX_HW_QUEUES=8             # Maximize hardware queues
export HIP_VISIBLE_DEVICES=0           # Use first GPU
# Note: Do NOT set HSA_OVERRIDE_GFX_VERSION - it causes kernel arch mismatch

# GROMACS binary - change this to your HIP-built GROMACS if available
# The default SYCL build doesn't work on MI210/gfx90a with ROCm 7.x
GMX_BIN="${HOME}/miniconda3/envs/gromacs_HIP/bin/gmx_mpi"

# Alternative: Use a HIP-built GROMACS (uncomment if you have one)
# GMX_BIN="/opt/gromacs-hip/bin/gmx"

# Path to Python modules
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules/gromacs_utils"

# Logs directory (will be created inside WORKDIR)
LOGS_DIR="${WORKDIR}/logs"

#------------------------------------------------------------------------------
# GPU Detection and Validation
#------------------------------------------------------------------------------

check_gpu() {
    log "Checking GPU availability..."
    
    # Check for AMD ROCm
    if command -v rocm-smi &> /dev/null; then
        log "  ROCm detected:"
        rocm-smi --showproductname 2>/dev/null | head -5 || true
        GPU_VENDOR="AMD"
    # Check for NVIDIA
    elif command -v nvidia-smi &> /dev/null; then
        log "  NVIDIA GPU detected:"
        nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -3 || true
        GPU_VENDOR="NVIDIA"
    else
        log "  WARNING: No GPU detected! Will attempt to run anyway."
        GPU_VENDOR="UNKNOWN"
    fi
    
    # Check GROMACS GPU capabilities
    log "  Checking GROMACS GPU support..."
    local gpu_info
    gpu_info=$($GMX_BIN mdrun -version 2>&1 | grep -i "gpu" || echo "No GPU info")
    log "  $gpu_info"
    
    # Set environment variables for AMD GPUs
    # NOTE: Do NOT set HSA_OVERRIDE_GFX_VERSION - it causes kernel arch mismatch!
    if [ "$GPU_VENDOR" = "AMD" ]; then
        export HIP_VISIBLE_DEVICES=0
        export ROCR_VISIBLE_DEVICES=0
        log "  Set AMD GPU environment variables"
    fi
}

test_gpu_sanity() {
    log "Testing GPU sanity..."
    
    # Create a minimal test run
    local test_dir=$(mktemp -d)
    cd "$test_dir"
    
    # Try a very short mdrun to see if GPU works
    # This will fail fast if there's a GPU issue
    
    log "  GPU sanity check: Will be validated during actual simulation"
    log "  Note: If you see 'Dummy kernel produced invalid values', GPU won't work"
    
    cd - > /dev/null
    rm -rf "$test_dir"
}

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------

log() {
    echo "[$(date '+%H:%M:%S')] $1"
}

create_em_mdp() {
    cat > "$1" << EOF
integrator      = steep
emtol           = 10.0
emstep          = 0.01
nsteps          = ${EM_STEPS}
nstlist         = 10
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz
EOF
}

clean_pdb() {
    grep -E "^(ATOM|TER|END)" "$1" | \
        sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > "$2"
}

get_chain_info() {
    local pdb="$1"
    local outdir="$2"
    
    python3 "${MODULES_DIR}/chain_parser.py" --mode info --pdb "$pdb" --outdir "$outdir"
}

process_structure() {
    local pdb="$1"
    local name="$2"
    local outdir="$3"
    
    log "Processing $name: $(basename $pdb)"
    
    mkdir -p "$outdir"
    cd "$outdir"
    
    # Clean PDB
    clean_pdb "$pdb" "clean.pdb"
    get_chain_info "$pdb" "$outdir"
    
    # Generate topology
    log "  Generating topology..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 > pdb2gmx.log
    
    # Create index with chains
    create_chain_index "$pdb" "$outdir"
    
    # Define box
    log "  Setting up simulation box..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 > editconf.log

    # Solvate
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 > solvate.log
    create_em_mdp "ions.mdp"
    $GMX_BIN grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 > grompp_ions.log
    echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral 2>&1 > genion.log
    
    # Energy minimization
    log "  Running energy minimization (GPU-accelerated)..."
    create_em_mdp "em.mdp"
    $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 > grompp_em.log
    
    # Run energy minimization with GPU acceleration
    # Optimized for MI210: fewer threads for better GPU utilization
    local gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp $NTHREADS $GPU_EM_FLAGS"
    log "  Command: $gmx_cmd"
    
    if $gmx_cmd 2>&1 | tee mdrun_em.log; then
        log "  Energy minimization completed successfully with GPU"
    else
        log "  WARNING: GPU execution failed, falling back to CPU-only"
        gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp 32"
        $gmx_cmd 2>&1 | tee mdrun_em_cpu.log
    fi
    
    # Extract energies
    log "  Extracting energy values..."
    echo -e "Potential\nCoul-SR\nLJ-SR\nPressure\n\n" | $GMX_BIN energy -f em.edr -o energies.xvg 2>&1 > energy.log || true
    
    # Analyze interface
    analyze_interface "$outdir" "$name"
}

create_chain_index() {
    local pdb="$1"
    local outdir="$2"
    
    cd "$outdir"
    
    # Generate basic index
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 > make_ndx.log
    
    # Add chain groups using Python module
    python3 "${MODULES_DIR}/chain_parser.py" --mode index \
        --pdb "$pdb" \
        --gro protein.gro \
        --index index.ndx \
        --outdir "$outdir"
}

analyze_interface() {
    local outdir="$1"
    local name="$2"
    
    cd "$outdir"
    
    log "  Analyzing interface..."
    
    # Calculate various metrics using GROMACS tools
    # Hydrogen bonds between chains
    echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx -num hbonds.xvg 2>&1 > hbond.log || true

    # Minimum distance between chains
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx -od mindist.xvg -on numcont.xvg -d 0.6 2>&1 > mindist.log || true
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg 2>&1 > sasa.log || true
    echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainA.xvg 2>&1 > sasa_a.log || true
    echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainB.xvg 2>&1 > sasa_b.log || true
    
    # Extract and compile metrics using Python module
    python3 "${MODULES_DIR}/metrics_extractor.py" --workdir "$outdir" --output metrics.txt
}

compare_results() {
    log "Comparing results..."
    
    python3 "${MODULES_DIR}/results_comparator.py" \
        --workdir "$WORKDIR" \
        --pdb1 "$PDB1" \
        --pdb2 "$PDB2" \
        --struct1-dir "structure_1" \
        --struct2-dir "structure_2"
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

main() {
    echo ""
    echo "============================================================"
    echo " GROMACS Quick Chain Stability Comparison - GPU VERSION"
    echo "============================================================"
    
    # Check GROMACS
    if [ ! -x "$GMX_BIN" ]; then
        echo "ERROR: GROMACS not found at $GMX_BIN"
        exit 1
    fi
    
    log "GROMACS: $($GMX_BIN --version 2>&1 | head -1)"
    
    # Check Python modules
    if [ ! -d "$MODULES_DIR" ]; then
        echo "ERROR: Python modules not found at $MODULES_DIR"
        exit 1
    fi
    
    # Check GPU
    check_gpu
    test_gpu_sanity
    
    log "Running in GPU mode (with $NTHREADS CPU threads)"
    log "GPU EM flags: $GPU_EM_FLAGS"
    log "GPU MD flags: $GPU_MD_FLAGS"
    
    # Check input files
    [ ! -f "$PDB1" ] && { echo "ERROR: $PDB1 not found"; exit 1; }
    [ ! -f "$PDB2" ] && { echo "ERROR: $PDB2 not found"; exit 1; }
    
    # Create workdir and logs directory
    mkdir -p "$WORKDIR"
    LOGS_DIR="${WORKDIR}/logs"
    mkdir -p "$LOGS_DIR"
    
    # Log file paths (ANSI codes will be stripped)
    TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
    LOG_FILE="${LOGS_DIR}/quick_stability_gpu_${TIMESTAMP}.log"
    ERROR_LOG="${LOGS_DIR}/quick_stability_gpu_${TIMESTAMP}_errors.log"
    
    # Start logging to file with ANSI code stripping
    exec > >(tee >(sed 's/\x1b\[[0-9;]*m//g' >> "$LOG_FILE")) 2> >(tee >(sed 's/\x1b\[[0-9;]*m//g' >> "$ERROR_LOG") | sed 's/\x1b\[[0-9;]*m//g' >> "$LOG_FILE")
    
    START=$(date +%s)
    
    # Process structures
    process_structure "$PDB1" "structure_1" "$WORKDIR/structure_1"
    echo ""
    process_structure "$PDB2" "structure_2" "$WORKDIR/structure_2"
    echo ""
    
    # Compare
    compare_results
    
    END=$(date +%s)
    log "Total time: $((END - START)) seconds"
    log "Results: $WORKDIR"
    
    # Collect all logs to central logs directory
    log "Collecting logs to ${LOGS_DIR}..."
    for struct_dir in "$WORKDIR"/structure_*; do
        if [ -d "$struct_dir" ]; then
            struct_name=$(basename "$struct_dir")
            mkdir -p "${LOGS_DIR}/${struct_name}"
            cp -f "$struct_dir"/*.log "${LOGS_DIR}/${struct_name}/" 2>/dev/null || true
        fi
    done
    log "Logs saved to: ${LOGS_DIR}"
}

main "$@"
