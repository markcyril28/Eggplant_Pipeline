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
WORKDIR="${3:-${HOME}/PPI/results_gromacs}"

# GROMACS parameters
FORCEFIELD="amber99sb-ildn"
WATERMODEL="tip3p"
BOX_DISTANCE=1.0
EM_STEPS=10000

# Performance settings
# Optimal threads for 128-thread server with GPU (16 for GPU, 32 for CPU-only)
NTHREADS=16
GPU_ID="auto"  # Auto-detect GPU

# Path to Python modules
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules/gromacs_utils"

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
    echo "1" | gmx pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 > pdb2gmx.log
    
    # Create index with chains
    create_chain_index "$pdb" "$outdir"
    
    # Define box
    log "  Setting up simulation box..."
    gmx editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 > editconf.log
    
    # Solvate
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 > solvate.log
    
    # Add ions
    create_em_mdp "ions.mdp"
    gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 > grompp_ions.log
    echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral 2>&1 > genion.log
    
    # Energy minimization
    log "  Running energy minimization..."
    create_em_mdp "em.mdp"
    gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 > grompp_em.log
    
    # Run with GPU if available
    local gmx_cmd="gmx mdrun -v -deffnm em -nt $NTHREADS"
    if [ "$GPU_ID" != "" ] && gmx -version 2>&1 | grep -q "GPU support.*enabled"; then
        if [ "$GPU_ID" = "auto" ]; then
            gmx_cmd="$gmx_cmd -nb gpu -pme gpu"
        else
            gmx_cmd="$gmx_cmd -gpu_id $GPU_ID -nb gpu -pme gpu"
        fi
    fi
    $gmx_cmd 2>&1 | tee mdrun_em.log
    
    # Extract energies
    log "  Extracting energy values..."
    echo -e "Potential\nCoul-SR\nLJ-SR\nPressure\n\n" | gmx energy -f em.edr -o energies.xvg 2>&1 > energy.log || true
    
    # Analyze interface
    analyze_interface "$outdir" "$name"
}

create_chain_index() {
    local pdb="$1"
    local outdir="$2"
    
    cd "$outdir"
    
    # Generate basic index
    echo "q" | gmx make_ndx -f protein.gro -o index.ndx 2>&1 > make_ndx.log
    
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
    echo -e "ChainA\nChainB" | gmx hbond -s em.tpr -f em.gro -n index.ndx -num hbonds.xvg 2>&1 > hbond.log || true
    
    # Minimum distance between chains
    echo -e "ChainA\nChainB" | gmx mindist -s em.tpr -f em.gro -n index.ndx -od mindist.xvg -on numcont.xvg -d 0.6 2>&1 > mindist.log || true
    
    # SASA calculations
    echo "Protein" | gmx sasa -s em.tpr -f em.gro -o sasa.xvg 2>&1 > sasa.log || true
    echo "ChainA" | gmx sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainA.xvg 2>&1 > sasa_a.log || true
    echo "ChainB" | gmx sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainB.xvg 2>&1 > sasa_b.log || true
    
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
    if ! command -v gmx &> /dev/null; then
        echo "ERROR: gmx not found. Please source GROMACS first."
        exit 1
    fi
    
    log "GROMACS: $(gmx --version 2>&1 | head -1)"
    
    # Check Python modules
    if [ ! -d "$MODULES_DIR" ]; then
        echo "ERROR: Python modules not found at $MODULES_DIR"
        exit 1
    fi
    
    # Detect GPU and adjust threads
    if gmx -version 2>&1 | grep -q "GPU support.*enabled"; then
        log "GPU support detected - using $NTHREADS threads with GPU acceleration"
    else
        log "No GPU detected - increasing to 32 CPU threads"
        NTHREADS=32
    fi
    
    # Check input files
    [ ! -f "$PDB1" ] && { echo "ERROR: $PDB1 not found"; exit 1; }
    [ ! -f "$PDB2" ] && { echo "ERROR: $PDB2 not found"; exit 1; }
    
    # Create workdir
    mkdir -p "$WORKDIR"
    
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
}

main "$@"
