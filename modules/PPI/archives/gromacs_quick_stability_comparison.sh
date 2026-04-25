#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS Quick Chain Stability Comparison Script
#
# A faster version that compares chain pairing stability using:
# 1. Energy minimization
# 2. Interface analysis (H-bonds, contacts, distances)
#
# This is faster than full MD but provides useful comparative metrics.
#
# Usage: ./quick_stability_comparison.sh [options]
#######################################################################

#------------------------------------------------------------------------------
# Configuration - MODIFY THESE PARAMETERS AS NEEDED
#------------------------------------------------------------------------------

PDB1="${1:-${HOME}/PPI/inputs/SmelGRF-GIF/fold_1_x_1_model_0.pdb}"
PDB2="${2:-${HOME}/PPI/inputs/SmelGRF-GIF/fold_2_x_1_model_0.pdb}"
WORKDIR="${3:-${HOME}/PPI/results_gromacs/quick_analysis}"

# GROMACS parameters
FORCEFIELD="amber99sb-ildn"
WATERMODEL="tip3p"
BOX_DISTANCE=1.0
EM_STEPS=10000

# Performance settings
# Using 32 threads for CPU-only mode (GPU disabled due to AdaptiveCpp SYCL bug)
# Note: GPU acceleration is disabled because AdaptiveCpp/hipSYCL produces incorrect
# kernel results on AMD MI210 (gfx90a) with ROCm 7.x - the sanity check kernel
# returns all zeros instead of expected values.
NTHREADS=64

# GROMACS binary (use gmx_mpi from conda)
GMX_BIN="${HOME}/miniconda3/envs/gromacs/bin/gmx_mpi"

# Path to Python modules
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules/gromacs_utils"

# Logs directory (will be created inside WORKDIR)
LOGS_DIR="${WORKDIR}/logs"

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
    log "  Running energy minimization..."
    create_em_mdp "em.mdp"
    $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 > grompp_em.log
    
    # Run energy minimization (CPU-only mode)
    # Note: GPU acceleration via SYCL/AdaptiveCpp is disabled due to a kernel bug
    # where the sanity check kernel produces all zeros instead of expected values.
    # This affects gfx90a (MI210) with ROCm 7.x and AdaptiveCpp develop branch.
    local gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp $NTHREADS"
    $gmx_cmd 2>&1 | tee mdrun_em.log
    
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
    echo " GROMACS Quick Chain Stability Comparison"
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
    
    # GPU is disabled due to AdaptiveCpp SYCL bug - using CPU-only mode
    log "Running in CPU-only mode ($NTHREADS threads)"
    log "Note: GPU disabled - AdaptiveCpp/SYCL kernel bug on MI210/gfx90a"
    
    # Check input files
    [ ! -f "$PDB1" ] && { echo "ERROR: $PDB1 not found"; exit 1; }
    [ ! -f "$PDB2" ] && { echo "ERROR: $PDB2 not found"; exit 1; }
    
    # Create workdir and logs directory
    mkdir -p "$WORKDIR"
    LOGS_DIR="${WORKDIR}/logs"
    mkdir -p "$LOGS_DIR"
    
    # Log file paths (ANSI codes will be stripped)
    TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
    LOG_FILE="${LOGS_DIR}/quick_stability_${TIMESTAMP}.log"
    ERROR_LOG="${LOGS_DIR}/quick_stability_${TIMESTAMP}_errors.log"
    
    # Start logging to file with ANSI code stripping
    # stdout -> main log, stderr -> both main log and error log
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
