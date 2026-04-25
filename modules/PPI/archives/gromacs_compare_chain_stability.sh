#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS Protein-Protein Interaction Stability Comparison Script
#
# This script compares the stability of chain pairings (protein-protein
# interactions) between two PDB structures using GROMACS.
#
# It performs:
# 1. Structure preparation and topology generation
# 2. Energy minimization
# 3. Short equilibration (NVT and NPT)
# 4. Short production MD
# 5. Interaction energy calculation between chains
# 6. Comparative analysis
#
# Usage: ./compare_chain_stability.sh [options]
#
# Requirements: GROMACS 2020+ installed and sourced
#######################################################################

#------------------------------------------------------------------------------
# Configuration - Modify these paths as needed
#------------------------------------------------------------------------------

# Input PDB files to compare
PDB1="${HOME}/PPI/inputs/SmelGRF-GIF/fold_1_x_1_model_0.pdb"
PDB2="${HOME}/PPI/inputs/SmelGRF-GIF/fold_2_x_1_model_0.pdb"

# Output directory
WORKDIR="${HOME}/PPI/results_gromacs/full_md"

# Force field selection (choose one: amber99sb-ildn, charmm27, oplsaa, gromos54a7)
FORCEFIELD="amber99sb-ildn"

# Water model
WATERMODEL="tip3p"

# Simulation parameters
BOX_DISTANCE=1.2          # nm from protein to box edge
ION_CONCENTRATION=0.15    # M NaCl
EM_STEPS=5000             # Energy minimization steps
NVT_STEPS=50000           # NVT equilibration (100 ps with 2fs timestep)
NPT_STEPS=50000           # NPT equilibration (100 ps with 2fs timestep)
MD_STEPS=250000           # Production MD (500 ps with 2fs timestep)

# Number of CPU threads for CPU-only mode
# Note: GPU acceleration is disabled because AdaptiveCpp/hipSYCL produces incorrect
# kernel results on AMD MI210 (gfx90a) with ROCm 7.x - the sanity check kernel
# returns all zeros instead of expected values.
NTHREADS=32

# GROMACS binary (use gmx_mpi from conda)
GMX_BIN="${HOME}/miniconda3/envs/gromacs/bin/gmx_mpi"

# Path to Python modules
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules/gromacs_utils"

# Logs directory (will be created inside WORKDIR)
LOGS_DIR="${WORKDIR}/logs"

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

log_section() {
    echo ""
    echo "=============================================================================="
    echo " $1"
    echo "=============================================================================="
}

check_gromacs() {
    if [ ! -x "$GMX_BIN" ]; then
        log_error "GROMACS not found at $GMX_BIN"
        log_error "Please verify the GMX_BIN path in the configuration"
        exit 1
    fi
    log_info "GROMACS version: $($GMX_BIN --version 2>&1 | head -1)"
    
    # GPU is disabled due to AdaptiveCpp SYCL bug - using CPU-only mode
    log_info "Running in CPU-only mode ($NTHREADS threads)"
    log_info "Note: GPU disabled - AdaptiveCpp/SYCL kernel bug on MI210/gfx90a"
}

check_pdb() {
    local pdb="$1"
    if [ ! -f "$pdb" ]; then
        log_error "PDB file not found: $pdb"
        exit 1
    fi
    
    # Check for chain identifiers
    local chains
    chains=$(grep "^ATOM" "$pdb" | cut -c22 | sort -u | tr '\n' ' ')
    log_info "Chains in $pdb: $chains"
}

clean_pdb() {
    local input="$1"
    local output="$2"
    
    # Remove HETATM (except common cofactors), ANISOU, and fix common issues
    grep -E "^(ATOM|TER|END)" "$input" | \
        sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > "$output"
}

#------------------------------------------------------------------------------
# MDP File Generation
#------------------------------------------------------------------------------

create_em_mdp() {
    local output="$1"
    cat > "$output" << 'EOF'
; Energy Minimization Parameters
integrator      = steep         ; Steepest descent minimization
emtol           = 100.0         ; Stop when max force < 100 kJ/mol/nm
emstep          = 0.01          ; Energy step size
nsteps          = 5000          ; Maximum number of steps

; Neighbor searching
nstlist         = 10            ; Frequency to update neighbor list
cutoff-scheme   = Verlet
ns_type         = grid          ; Method to determine neighbor list

; Electrostatics
coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.0           ; Short-range electrostatic cutoff (nm)

; VdW
vdwtype         = Cut-off
rvdw            = 1.0           ; Short-range van der Waals cutoff (nm)

; Periodic boundary conditions
pbc             = xyz           ; 3D PBC
EOF
    sed -i "s/nsteps          = 5000/nsteps          = $EM_STEPS/" "$output"
}

create_nvt_mdp() {
    local output="$1"
    cat > "$output" << 'EOF'
; NVT Equilibration Parameters
define          = -DPOSRES      ; Position restrain the protein

integrator      = md            ; Leap-frog integrator
nsteps          = 50000         ; 100 ps
dt              = 0.002         ; 2 fs timestep

; Output control
nstxout         = 500           ; Save coordinates every 1 ps
nstvout         = 500           ; Save velocities every 1 ps
nstenergy       = 500           ; Save energies every 1 ps
nstlog          = 500           ; Update log file every 1 ps
nstxout-compressed = 500        ; Save compressed coordinates every 1 ps

; Bond parameters
continuation    = no            ; First dynamics run
constraint_algorithm = lincs    ; Holonomic constraints
constraints     = h-bonds       ; H-bonds constrained
lincs_iter      = 1             ; Accuracy of LINCS
lincs_order     = 4             ; Also related to accuracy

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20            ; 20 fs, largely irrelevant with Verlet
rcoulomb        = 1.0           ; Short-range electrostatic cutoff (nm)
rvdw            = 1.0           ; Short-range van der Waals cutoff (nm)

; Electrostatics
coulombtype     = PME           ; Particle Mesh Ewald
pme_order       = 4             ; Cubic interpolation
fourierspacing  = 0.12          ; Grid spacing for FFT

; Temperature coupling
tcoupl          = V-rescale     ; Modified Berendsen thermostat
tc-grps         = Protein Non-Protein ; Two coupling groups
tau_t           = 0.1     0.1   ; Time constant (ps)
ref_t           = 300     300   ; Reference temperature (K)

; Pressure coupling - off for NVT
pcoupl          = no

; Periodic boundary conditions
pbc             = xyz           ; 3D PBC

; Velocity generation
gen_vel         = yes           ; Generate velocities from Maxwell distribution
gen_temp        = 300           ; Temperature for velocity generation
gen_seed        = -1            ; Random seed
EOF
    sed -i "s/nsteps          = 50000/nsteps          = $NVT_STEPS/" "$output"
}

create_npt_mdp() {
    local output="$1"
    cat > "$output" << 'EOF'
; NPT Equilibration Parameters
define          = -DPOSRES      ; Position restrain the protein

integrator      = md            ; Leap-frog integrator
nsteps          = 50000         ; 100 ps
dt              = 0.002         ; 2 fs timestep

; Output control
nstxout         = 500           ; Save coordinates every 1 ps
nstvout         = 500           ; Save velocities every 1 ps
nstenergy       = 500           ; Save energies every 1 ps
nstlog          = 500           ; Update log file every 1 ps
nstxout-compressed = 500        ; Save compressed coordinates

; Bond parameters
continuation    = yes           ; Continuing after NVT
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = 300     300

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0           ; Time constant (ps)
ref_p           = 1.0           ; Reference pressure (bar)
compressibility = 4.5e-5        ; Isothermal compressibility of water

; Periodic boundary conditions
pbc             = xyz

; Velocity generation - no, continuing from NVT
gen_vel         = no
EOF
    sed -i "s/nsteps          = 50000/nsteps          = $NPT_STEPS/" "$output"
}

create_md_mdp() {
    local output="$1"
    cat > "$output" << 'EOF'
; Production MD Parameters
integrator      = md            ; Leap-frog integrator
nsteps          = 250000        ; 500 ps
dt              = 0.002         ; 2 fs timestep

; Output control
nstxout         = 0             ; Don't save .trr coordinates
nstvout         = 0             ; Don't save .trr velocities
nstfout         = 0             ; Don't save forces
nstenergy       = 500           ; Save energies every 1 ps
nstlog          = 500           ; Update log file every 1 ps
nstxout-compressed = 1000       ; Save .xtc coordinates every 2 ps

; Bond parameters
continuation    = yes           ; Continuing after NPT
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = 300     300

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; Periodic boundary conditions
pbc             = xyz

; Velocity generation - no
gen_vel         = no
EOF
    sed -i "s/nsteps          = 250000/nsteps          = $MD_STEPS/" "$output"
}

create_ie_mdp() {
    # MDP for rerunning trajectory to calculate interaction energies
    local output="$1"
    cat > "$output" << 'EOF'
; Interaction Energy Calculation (rerun)
integrator      = md
nsteps          = 0
dt              = 0.002

; Energy output - detailed
nstenergy       = 1
nstlog          = 1

; Energy groups for interaction energy calculation
energygrps      = ChainA ChainB

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12

; PBC
pbc             = xyz
EOF
}

#------------------------------------------------------------------------------
# Processing Functions
#------------------------------------------------------------------------------

process_structure() {
    local pdb_input="$1"
    local name="$2"
    local outdir="$3"
    
    log_section "Processing: $name"
    
    mkdir -p "$outdir"
    cd "$outdir"
    
    # Clean PDB
    log_info "Cleaning PDB file..."
    clean_pdb "$pdb_input" "clean.pdb"
    
    # Generate topology
    log_info "Generating topology with pdb2gmx..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 | tee pdb2gmx.log
    
    # Create index file with chain selections
    log_info "Creating index file with chain groups..."
    create_chain_index "$outdir"
    
    # Define box
    log_info "Defining simulation box..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 | tee editconf.log
    
    # Solvate
    log_info "Solvating the system..."
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tee solvate.log
    
    # Add ions
    log_info "Adding ions..."
    create_em_mdp "ions.mdp"
    $GMX_BIN grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 | tee grompp_ions.log
    echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc $ION_CONCENTRATION 2>&1 | tee genion.log
    
    # Energy minimization
    log_info "Running energy minimization..."
    create_em_mdp "em.mdp"
    $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 | tee grompp_em.log
    
    # Run energy minimization (CPU-only mode)
    # Note: GPU acceleration via SYCL/AdaptiveCpp is disabled due to a kernel bug
    # where the sanity check kernel produces all zeros instead of expected values.
    local gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp $NTHREADS"
    $gmx_cmd 2>&1 | tee mdrun_em.log
    
    # Extract EM potential energy
    echo "Potential" | $GMX_BIN energy -f em.edr -o em_potential.xvg 2>&1 | tee energy_em.log
    
    # NVT equilibration
    log_info "Running NVT equilibration..."
    create_nvt_mdp "nvt.mdp"
    $GMX_BIN grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 2>&1 | tee grompp_nvt.log
    
    # Run NVT (CPU-only mode due to AdaptiveCpp SYCL bug)
    gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm nvt -ntomp $NTHREADS"
    $gmx_cmd 2>&1 | tee mdrun_nvt.log
    
    # NPT equilibration  
    log_info "Running NPT equilibration..."
    create_npt_mdp "npt.mdp"
    $GMX_BIN grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr 2>&1 | tee grompp_npt.log
    
    # Run NPT (CPU-only mode due to AdaptiveCpp SYCL bug)
    gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm npt -ntomp $NTHREADS"
    $gmx_cmd 2>&1 | tee mdrun_npt.log
    
    # Production MD
    log_info "Running production MD..."
    create_md_mdp "md.mdp"
    $GMX_BIN grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr 2>&1 | tee grompp_md.log
    
    # Run production MD (CPU-only mode due to AdaptiveCpp SYCL bug)
    gmx_cmd="mpirun -np 1 $GMX_BIN mdrun -v -deffnm md -ntomp $NTHREADS"
    $gmx_cmd 2>&1 | tee mdrun_md.log
    
    log_info "MD simulation completed for $name"
}

create_chain_index() {
    local outdir="$1"
    cd "$outdir"
    
    # Get chain information from the original PDB
    # Create a custom index file with Chain A and Chain B groups
    
    # First, generate basic index
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 | tee make_ndx.log
    
    # Create chain-specific index groups using Python module
    python3 "${MODULES_DIR}/chain_indexer.py" \
        --gro protein.gro \
        --pdb clean.pdb \
        --index index.ndx
}

calculate_interaction_energy() {
    local outdir="$1"
    local name="$2"
    
    log_section "Calculating Interaction Energy: $name"
    
    cd "$outdir"
    
    # Create MDP for energy rerun
    create_ie_mdp "ie.mdp"
    
    # Prepare TPR with energy groups
    # We need to modify the topology to include energy group definitions
    $GMX_BIN grompp -f ie.mdp -c md.gro -p topol.top -n index.ndx -o ie.tpr -maxwarn 5 2>&1 | tee grompp_ie.log
    
    # Rerun trajectory to calculate energies
    log_info "Rerunning trajectory for energy calculation..."
    $GMX_BIN mdrun -rerun md.xtc -s ie.tpr -e ie.edr -g ie.log 2>&1 | tee mdrun_ie.log
    
    # Extract interaction energies
    log_info "Extracting interaction energies..."
    
    # Get Coulomb and LJ energies between chains
    echo -e "Coul-SR:ChainA-ChainB\nLJ-SR:ChainA-ChainB\n\n" | \
        $GMX_BIN energy -f ie.edr -o interaction_energy.xvg 2>&1 | tee energy_ie.log || true
    
    # Alternative: Calculate total energies if energy groups fail
    if [ ! -s "interaction_energy.xvg" ]; then
        log_info "Energy groups not available, using alternative analysis..."
        
        # Calculate RMSD for stability assessment
        echo -e "Backbone\nBackbone" | $GMX_BIN rms -s md.tpr -f md.xtc -o rmsd.xvg 2>&1 | tee rmsd.log
        
        # Calculate radius of gyration
        echo "Protein" | $GMX_BIN gyrate -s md.tpr -f md.xtc -o gyrate.xvg 2>&1 | tee gyrate.log
        
        # Calculate RMSF
        echo "Backbone" | $GMX_BIN rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res 2>&1 | tee rmsf.log
        
        # Calculate hydrogen bonds between chains
        echo -e "ChainA\nChainB" | $GMX_BIN hbond -s md.tpr -f md.xtc -n index.ndx -num hbonds.xvg 2>&1 | tee hbond.log || true
        
        # Calculate contact map / minimum distance between chains
        echo -e "ChainA\nChainB" | $GMX_BIN mindist -s md.tpr -f md.xtc -n index.ndx -od mindist.xvg 2>&1 | tee mindist.log || true
    fi
    
    # Calculate binding energy using MM-PBSA/GBSA approach (simplified)
    calculate_binding_metrics "$outdir" "$name"
}

calculate_binding_metrics() {
    local outdir="$1"
    local name="$2"
    
    cd "$outdir"
    
    log_info "Calculating binding metrics for $name..."
    
    # Extract final frame for analysis
    echo "System" | $GMX_BIN trjconv -s md.tpr -f md.xtc -o final_frame.gro -dump 0 -pbc mol 2>&1 || true
    
    # Calculate SASA (Solvent Accessible Surface Area)
    echo "Protein" | $GMX_BIN sasa -s md.tpr -f md.xtc -o sasa.xvg -or resarea.xvg 2>&1 | tee sasa.log || true
    
    # Analyze stability metrics using Python module
    python3 "${MODULES_DIR}/stability_analyzer.py" \
        --workdir "$outdir" \
        --output stability_metrics.txt
}

compare_results() {
    log_section "Comparing Stability Results"
    
    python3 "${MODULES_DIR}/md_results_comparator.py" \
        --workdir "$WORKDIR" \
        --pdb1 "$PDB1" \
        --pdb2 "$PDB2" \
        --struct1-dir "structure_1" \
        --struct2-dir "structure_2"
}

#------------------------------------------------------------------------------
# Main Execution
#------------------------------------------------------------------------------

main() {
    log_section "GROMACS Chain Stability Comparison"
    
    # Verify requirements
    check_gromacs
    check_pdb "$PDB1"
    check_pdb "$PDB2"
    
    # Check Python modules
    if [ ! -d "$MODULES_DIR" ]; then
        log_error "Python modules not found at $MODULES_DIR"
        exit 1
    fi
    
    # Create working directory and logs directory
    mkdir -p "$WORKDIR"
    LOGS_DIR="${WORKDIR}/logs"
    mkdir -p "$LOGS_DIR"
    cd "$WORKDIR"
    
    # Log file paths (ANSI codes will be stripped)
    TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
    LOG_FILE="${LOGS_DIR}/chain_stability_${TIMESTAMP}.log"
    ERROR_LOG="${LOGS_DIR}/chain_stability_${TIMESTAMP}_errors.log"
    
    # Start logging to file with ANSI code stripping
    # stdout -> main log, stderr -> both main log and error log
    exec > >(tee >(sed 's/\x1b\[[0-9;]*m//g' >> "$LOG_FILE")) 2> >(tee >(sed 's/\x1b\[[0-9;]*m//g' >> "$ERROR_LOG") | sed 's/\x1b\[[0-9;]*m//g' >> "$LOG_FILE")
    
    # Record start time
    START_TIME=$(date +%s)
    
    # Process first structure
    process_structure "$PDB1" "structure_1" "$WORKDIR/structure_1"
    calculate_interaction_energy "$WORKDIR/structure_1" "structure_1"
    
    # Process second structure
    process_structure "$PDB2" "structure_2" "$WORKDIR/structure_2"
    calculate_interaction_energy "$WORKDIR/structure_2" "structure_2"
    
    # Compare results
    compare_results
    
    # Record end time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    log_section "Comparison Complete"
    log_info "Total time: $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"
    log_info "Results saved to: $WORKDIR"
    
    # Collect all logs to central logs directory
    log_info "Collecting logs to ${LOGS_DIR}..."
    for struct_dir in "$WORKDIR"/structure_*; do
        if [ -d "$struct_dir" ]; then
            struct_name=$(basename "$struct_dir")
            mkdir -p "${LOGS_DIR}/${struct_name}"
            cp -f "$struct_dir"/*.log "${LOGS_DIR}/${struct_name}/" 2>/dev/null || true
        fi
    done
    log_info "Logs saved to: ${LOGS_DIR}"
}

# Run main function
main "$@"
