#!/usr/bin/env bash
set -euo pipefail

#######################################################################
# GROMACS Unified Pipeline
# Runs all GROMACS PPI analysis steps from a single script.
# Steps are configured in config/PPI/gromacs/gromacs_pipeline.toml.
#
# Usage:
#   ./run_gromacs_pipeline.sh                          # Run all enabled steps
#   ./run_gromacs_pipeline.sh --steps quick_stability,interface_analysis
#   ./run_gromacs_pipeline.sh --steps compare_chain_stability --mode sim-only
#   ./run_gromacs_pipeline.sh --dry-run                # Show what would run
#   ./run_gromacs_pipeline.sh --list                   # List available steps
#
# Output structure (auto-generated):
#   <output_base>/
#   └── run_YYYYMMDD_HHMMSS/          (or flat if auto_timestamp_dir=false)
#       ├── 1_quick_stability/
#       │   └── <structure_names>/
#       ├── 2_compare_chain_stability/
#       │   └── <structure_names>/
#       ├── 3_interface_analysis/
#       │   └── <structure_name>/
#       ├── 4_batch_comparison/
#       │   └── <dataset_name>/
#       └── pipeline_summary.txt
#######################################################################

#------------------------------------------------------------------------------
# STEPS TO RUN  —  comment/uncomment to toggle individual steps
# (ignored when --steps is passed on the command line)
#------------------------------------------------------------------------------

INLINE_STEPS=(
    #"quick_stability"
    "compare_chain_stability"
    "interface_analysis"
    #"batch_comparison"
    "production_md"
    "visualize_results"
)

#------------------------------------------------------------------------------
# GENE GROUP  —  structures are loaded from config/<GENE_GROUP>/
#------------------------------------------------------------------------------

GENE_GROUP="DMP-HAP2"

#------------------------------------------------------------------------------
# STRUCTURES TO ANALYZE  -  override config/<GENE_GROUP>/h_protein_structure_analysis_<GENE_GROUP>.toml
# Paths are relative to INPUT_BASE (set in common.toml)
# Leave empty to use [ppi_structures].active from the gene-group config.
#------------------------------------------------------------------------------

INLINE_STRUCTURES=(
)

#------------------------------------------------------------------------------
# ARGUMENT PARSING
#------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Logging setup
PROJECT_ROOT="$SCRIPT_DIR"
mkdir -p "$PROJECT_ROOT/logs"
source "${SCRIPT_DIR}/modules/logging/logging_utils.sh"

# Defaults
CLI_STEPS=""
CLI_MODE=""
DRY_RUN=false
LIST_ONLY=false

show_help() {
    sed -n '3,26p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --steps)      [[ $# -gt 1 ]] || { echo "Option $1 requires an argument" >&2; exit 1; }
                      CLI_STEPS="$2";  shift 2 ;;
        --mode|-m)    [[ $# -gt 1 ]] || { echo "Option $1 requires an argument" >&2; exit 1; }
                      CLI_MODE="$2";   shift 2 ;;
        --dry-run)    DRY_RUN=true;    shift   ;;
        --list)       LIST_ONLY=true;  shift   ;;
        -h|--help)    show_help              ;;
        *)            echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------

source "${SCRIPT_DIR}/modules/PPI/config_parser.sh"

CONFIG_DIR="${SCRIPT_DIR}/config/PPI"

# Load combined pipeline config
load_config "${SCRIPT_DIR}/10_run_gromacs_pipelineCONFIG.toml"
INPUT_BASE=$(toml_get "paths.input_base" "III_RESULT/DMP-HAP2/08_Protein_Structure/GPE001970_SMEL5/AlphaFold3_Results")
OUTPUT_BASE=$(toml_get "paths.results.gromacs" "III_RESULT/DMP-HAP2/11_PPI")

# Resolve relative paths to absolute (prevents breakage after pushd)
[[ "$INPUT_BASE" != /* ]] && INPUT_BASE="${SCRIPT_DIR}/${INPUT_BASE}"
[[ "$OUTPUT_BASE" != /* ]] && OUTPUT_BASE="${SCRIPT_DIR}/${OUTPUT_BASE}"

STOP_ON_ERROR=$(toml_get "pipeline.stop_on_error" "false")
AUTO_TIMESTAMP=$(toml_get "pipeline.auto_timestamp_dir" "true")

# Build steps list: CLI > inline list > config
STEPS=()
if [[ -n "$CLI_STEPS" ]]; then
    IFS=',' read -ra STEPS <<< "$CLI_STEPS"
elif [[ ${#INLINE_STEPS[@]} -gt 0 ]]; then
    STEPS=("${INLINE_STEPS[@]}")
else
    while IFS= read -r step; do
        [[ -n "$step" ]] && STEPS+=("$step")
    done < <(toml_get_array "pipeline.steps")
fi

# Build shared structures list: inline list > gene-group config > pipeline config
SHARED_STRUCTURES=()
if [[ ${#INLINE_STRUCTURES[@]} -gt 0 ]]; then
    for rel_path in "${INLINE_STRUCTURES[@]}"; do
        SHARED_STRUCTURES+=("${INPUT_BASE}/${rel_path}")
    done
else
    # Load structures from gene-group config (config/<GENE_GROUP>/h_protein_structure_analysis_<GENE_GROUP>.toml)
    GENE_STRUCTURES_CONFIG="${SCRIPT_DIR}/config/${GENE_GROUP}/h_protein_structure_analysis_${GENE_GROUP}.toml"
    if [[ -f "$GENE_STRUCTURES_CONFIG" ]]; then
        load_config "$GENE_STRUCTURES_CONFIG"
        while IFS= read -r rel_path; do
            [[ -n "$rel_path" ]] && SHARED_STRUCTURES+=("${INPUT_BASE}/${rel_path}")
        done < <(toml_get_array "ppi_structures.active")
    fi
    # Fallback to pipeline config if gene-group config had no structures
    if [[ ${#SHARED_STRUCTURES[@]} -eq 0 ]]; then
        load_config "${SCRIPT_DIR}/10_run_gromacs_pipelineCONFIG.toml"
        while IFS= read -r rel_path; do
            [[ -n "$rel_path" ]] && SHARED_STRUCTURES+=("${INPUT_BASE}/${rel_path}")
        done < <(toml_get_array "structures.active")
    fi
    # Restore pipeline config as active
    load_config "${SCRIPT_DIR}/10_run_gromacs_pipelineCONFIG.toml"
fi

# Source common GROMACS functions
source "${SCRIPT_DIR}/modules/PPI/gromacs_common.sh"

# Override gromacs_common.sh simple log functions with logging_utils.sh versions
# so output is captured in the central logs/ directory
log()         { log_info "$@"; }
log_section() { log_step "$@"; }

setup_logging "${GENE_GROUP}_gromacs_pipeline"
trap 'wait; teardown_logging' EXIT

#------------------------------------------------------------------------------
# LIST MODE
#------------------------------------------------------------------------------

if [[ "$LIST_ONLY" == true ]]; then
    echo "Available pipeline steps:"
    echo "  quick_stability           - Energy minimization comparison (fast)"
    echo "  compare_chain_stability   - Full MD chain stability comparison (GPU)"
    echo "  interface_analysis        - Detailed PPI interface analysis"
    echo "  batch_comparison          - Batch analysis of all structures in dataset folders"
    echo "  production_md             - Full MD simulation with maximum GPU offload"
    echo "  visualize_results         - (Re)generate plots and viz scripts from prior outputs"
    echo ""
    echo "Enabled in config:"
    for step in "${STEPS[@]}"; do
        local_enabled=$(toml_get "step.${step}.enabled" "true")
        echo "  [$([ "$local_enabled" == "true" ] && echo "x" || echo " ")] $step"
    done
    exit 0
fi

#------------------------------------------------------------------------------
# OUTPUT DIRECTORY SETUP
#------------------------------------------------------------------------------

if [[ "$AUTO_TIMESTAMP" == "true" ]]; then
    RUN_DIR="${OUTPUT_BASE}/run_$(date '+%Y%m%d_%H%M%S')"
else
    RUN_DIR="${OUTPUT_BASE}"
fi

# Step number mapping for output folder prefixes
declare -A STEP_PREFIX=(
    [quick_stability]="1_quick_stability"
    [compare_chain_stability]="2_compare_chain_stability"
    [interface_analysis]="3_interface_analysis"
    [batch_comparison]="4_batch_comparison"
    [production_md]="5_production_md"
    [visualize_results]="6_visualize_results"
)

#------------------------------------------------------------------------------
# HELPER: Load per-step structures or fall back to shared list
#------------------------------------------------------------------------------

load_step_structures() {
    local step_name="$1"
    local -n _out_arr=$2  # nameref to output array

    _out_arr=()
    local has_override=false

    while IFS= read -r rel_path; do
        if [[ -n "$rel_path" ]]; then
            _out_arr+=("${INPUT_BASE}/${rel_path}")
            has_override=true
        fi
    done < <(toml_get_array "step.${step_name}.structures.active")

    if [[ "$has_override" != true ]]; then
        _out_arr=("${SHARED_STRUCTURES[@]}")
    fi
}

#------------------------------------------------------------------------------
# HELPER: Generate short name from PDB path
#------------------------------------------------------------------------------

shorten_pdb_name() {
    local pdb="$1"
    basename "$pdb" .pdb \
        | sed 's/_model_0.*//; s/_20[0-9][0-9]_[0-9][0-9]_[0-9][0-9].*//; s/^fold_//; s/^[a-zA-Z0-9]*_and_//'
}

# Build a compact workdir name from an array of shortened structure names.
# When >3 structures, abbreviates to: first_vs_last_Nstructs
make_workdir_name() {
    local -n _mwn_names=$1
    local sep="${2:-_vs_}"
    if [[ ${#_mwn_names[@]} -le 3 ]]; then
        local joined="${_mwn_names[0]}"
        for ((k=1; k<${#_mwn_names[@]}; k++)); do
            joined="${joined}${sep}${_mwn_names[$k]}"
        done
        echo "$joined"
    else
        echo "${_mwn_names[0]}${sep}${_mwn_names[${#_mwn_names[@]}-1]}_${#_mwn_names[@]}structs"
    fi
}

# Write a MANIFEST.txt mapping structure_N → original PDB path
write_manifest() {
    local outdir="$1"; shift
    local -a pdb_list=("$@")
    {
        echo "# Structure Manifest — maps structure_N to original PDB"
        echo "# Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        for i in "${!pdb_list[@]}"; do
            echo "structure_$((i+1))  $(basename "${pdb_list[$i]}")"
        done
    } > "${outdir}/MANIFEST.txt"
}

#------------------------------------------------------------------------------
# STEP 1: QUICK STABILITY
#------------------------------------------------------------------------------

run_quick_stability() {
    local step_dir="$1"

    # Load settings
    local OVERRIDE_EXISTING=$(toml_get "step.quick_stability.override_existing" "false")
    local BOX_DISTANCE=$(toml_get "step.quick_stability.box_distance" "1.0")
    local EM_STEPS=$(toml_get "step.quick_stability.em_steps" "10000")
    local MAX_THREADS=$(toml_get "step.quick_stability.max_threads" "0")
    local GPU_ID=$(toml_get "step.quick_stability.gpu_id" "0")
    NTHREADS=$(( MAX_THREADS > 0 ? MAX_THREADS : $(nproc) ))

    # Load structures
    local PDB_LIST=()
    load_step_structures "quick_stability" PDB_LIST

    if [[ ${#PDB_LIST[@]} -lt 2 ]]; then
        log_error "quick_stability requires at least 2 structures"
        return 1
    fi

    # Generate output dir name from structure names
    local names=()
    for pdb in "${PDB_LIST[@]}"; do
        names+=("$(shorten_pdb_name "$pdb")")
    done
    local joined
    joined=$(make_workdir_name names "_")
    local OUTPUT_DIR="${step_dir}/${joined}"

    log "Structures: ${#PDB_LIST[@]}"
    for i in "${!PDB_LIST[@]}"; do
        log "  Structure $((i+1)): $(basename "${PDB_LIST[$i]}")"
    done
    log "Output: $OUTPUT_DIR"
    echo ""

    mkdir -p "$OUTPUT_DIR"
    write_manifest "$OUTPUT_DIR" "${PDB_LIST[@]}"
    check_gromacs || return 1
    check_python_modules || return 1

    # Process each structure
    for i in "${!PDB_LIST[@]}"; do
        local struct_num=$((i+1))
        local struct_name="structure_${struct_num}"
        local pdb="${PDB_LIST[$i]}"
        local outdir="$OUTPUT_DIR/$struct_name"

        log "Processing $struct_name: $(basename "$pdb")"

        if [[ "$OVERRIDE_EXISTING" != "true" && -f "$outdir/metrics.txt" && -f "$outdir/em.gro" ]]; then
            log "  ✓ Already processed, skipping"
            continue
        fi

        [[ "$OVERRIDE_EXISTING" == "true" && -d "$outdir" ]] && rm -rf "$outdir"

        mkdir -p "$outdir/logs"
        cd "$outdir"

        if ! python3 -m gromacs_utils.cli prepare-structure "$pdb" -o clean.pdb; then
            log_error "Failed to prepare structure: $pdb"
            continue
        fi
        get_chain_info "$pdb" "$outdir"

        log "  Generating topology..."
        echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > logs/pdb2gmx.log 2>&1
        echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx > logs/make_ndx.log 2>&1
        create_chain_index clean.pdb protein.gro index.ndx

        log "  Setting up simulation box..."
        $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron > logs/editconf.log 2>&1
        $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > logs/solvate.log 2>&1

        python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 10.0

        $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 > logs/grompp_ions.log 2>&1
        echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral > logs/genion.log 2>&1

        log "  Running energy minimization (GPU)..."
        $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr > logs/grompp_em.log 2>&1
        run_em em logs

        if [[ ! -f em.gro ]]; then
            log_error "EM failed — em.gro not produced for $struct_name"
            continue
        fi

        log "  Extracting metrics..."
        echo -e "Potential\nCoul-SR\nLJ-SR\nPressure\n\n" | $GMX_BIN energy -f em.edr -o energies.xvg > logs/energy.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx -num hbonds.xvg > logs/hbond.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx -od mindist.xvg -on numcont.xvg -d 0.6 > logs/mindist.log 2>&1 || true
        echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg > logs/sasa.log 2>&1 || true
        echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainA.xvg > logs/sasa_a.log 2>&1 || true
        echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainB.xvg > logs/sasa_b.log 2>&1 || true

        extract_metrics "$outdir" metrics.txt
        generate_visualization "$outdir" interface
        echo ""
    done

    # Compare results
    log "Comparing results..."
    local summary_file="$OUTPUT_DIR/structures_summary.txt"
    {
        echo "=============================================="
        echo "QUICK STABILITY COMPARISON SUMMARY"
        echo "=============================================="
        echo ""
        echo "Structures analyzed: ${#PDB_LIST[@]}"
        echo "Generated: $(date)"
        echo ""

        for i in "${!PDB_LIST[@]}"; do
            local struct_num=$((i+1))
            local struct_dir="$OUTPUT_DIR/structure_${struct_num}"
            local pdb_name=$(basename "${PDB_LIST[$i]}" .pdb)

            echo "--- Structure $struct_num: $pdb_name ---"
            if [[ -f "$struct_dir/metrics.txt" ]]; then
                cat "$struct_dir/metrics.txt" | sed 's/^/  /'
            else
                echo "  (metrics unavailable)"
            fi
            echo ""
        done
    } > "$summary_file"
    log "Summary: $summary_file"

    # Pairwise comparisons
    local num_structs=${#PDB_LIST[@]}
    if [[ $num_structs -ge 2 ]]; then
        log "Running pairwise comparisons..."
        for ((i=0; i<num_structs-1; i++)); do
            for ((j=i+1; j<num_structs; j++)); do
                local s1="structure_$((i+1))" s2="structure_$((j+1))"
                if [[ -d "$OUTPUT_DIR/$s1" && -d "$OUTPUT_DIR/$s2" ]]; then
                    python3 -m gromacs_utils.results_comparator \
                        --workdir "$OUTPUT_DIR" \
                        --pdb1 "${PDB_LIST[$i]}" --pdb2 "${PDB_LIST[$j]}" \
                        --struct1-dir "$s1" --struct2-dir "$s2" \
                        --output "comparison_$((i+1))_vs_$((j+1)).txt" 2>&1 || true
                fi
            done
        done

        if [[ $num_structs -gt 2 ]]; then
            local pdb_args=() struct_args=()
            for i in "${!PDB_LIST[@]}"; do
                pdb_args+=("${PDB_LIST[$i]}")
                struct_args+=("structure_$((i+1))")
            done
            python3 -m gromacs_utils.results_comparator \
                --workdir "$OUTPUT_DIR" --multi \
                --pdb-list "${pdb_args[@]}" --struct-dirs "${struct_args[@]}" 2>&1 || true
        fi
    fi

    cat "$summary_file" 2>/dev/null || true
}

#------------------------------------------------------------------------------
# STEP 2: COMPARE CHAIN STABILITY (GPU)
#------------------------------------------------------------------------------

run_compare_chain_stability() {
    local step_dir="$1"

    # Load settings
    local OVERRIDE_EXISTING=$(toml_get "step.compare_chain_stability.override_existing" "false")
    local MODE="${CLI_MODE:-$(toml_get "step.compare_chain_stability.mode" "full")}"
    local BOX_DISTANCE=$(toml_get "step.compare_chain_stability.box_distance" "1.5")
    local ION_CONCENTRATION=$(toml_get "step.compare_chain_stability.ion_concentration" "0.15")
    local EM_STEPS=$(toml_get "step.compare_chain_stability.em_steps" "50000")
    local NVT_STEPS=$(toml_get "step.compare_chain_stability.nvt_steps" "100000")
    local NPT_STEPS=$(toml_get "step.compare_chain_stability.npt_steps" "100000")
    local MD_STEPS=$(toml_get "step.compare_chain_stability.md_steps" "5000000")
    local GPU_ID=$(toml_get "step.compare_chain_stability.gpu_id" "0")
    local MAX_THREADS=$(toml_get "step.compare_chain_stability.max_threads" "0")
    local MDRUN_OPTION=$(toml_get "step.compare_chain_stability.mdrun_option" "Option_A")
    NTHREADS=$(( MAX_THREADS > 0 ? MAX_THREADS : $(nproc) ))

    # GPU flag selection
    case "$MDRUN_OPTION" in
        Option_A)
            local GPU_MD_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu"
            local GPU_EM_FLAGS="-nb gpu -pme cpu -bonded gpu"
            ;;
        Option_B)
            local GPU_MD_FLAGS="-nb gpu -pme gpu -bonded cpu -update gpu"
            local GPU_EM_FLAGS="-nb gpu -pme cpu -bonded cpu"  # Conservative: bonded on CPU
            ;;
        *)
            log_error "Unknown MDRUN_OPTION '$MDRUN_OPTION'"
            return 1
            ;;
    esac

    # Load structures
    local PDB_LIST=()
    load_step_structures "compare_chain_stability" PDB_LIST

    if [[ ${#PDB_LIST[@]} -lt 2 ]]; then
        log_error "compare_chain_stability requires at least 2 structures"
        return 1
    fi

    # Generate workdir name
    local names=()
    for pdb in "${PDB_LIST[@]}"; do
        names+=("$(shorten_pdb_name "$pdb")")
    done
    local joined
    joined=$(make_workdir_name names)
    local WORKDIR="${step_dir}/${joined}"

    log "Mode: $MODE"
    log "Structures: ${#PDB_LIST[@]}"
    for i in "${!PDB_LIST[@]}"; do
        log "  Structure $((i+1)): $(basename "${PDB_LIST[$i]}")"
    done
    log "Output: $WORKDIR"

    mkdir -p "$WORKDIR/logs"
    write_manifest "$WORKDIR" "${PDB_LIST[@]}"
    cd "$WORKDIR"

    # GPU detection — use centralized detection from gromacs_common.sh
    USE_GPU=true
    export USE_GPU
    local GPU_INFO=""

    _check_gpu() {
        if _detect_gpu_available; then
            GPU_INFO=$( (nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || rocm-smi --showproductname 2>/dev/null) | head -1 )
            log "  Found GPU: ${GPU_INFO:-detected}"
        else
            log_warn "GPU not usable by GROMACS, will use CPU mode"
            USE_GPU=false
            export USE_GPU
        fi
        return 0
    }

    _build_mdrun_cmd() {
        local stage="$1" deffnm="$2"
        local cmd="$GMX_BIN mdrun -v -deffnm $deffnm -ntmpi 1 -ntomp $NTHREADS"
        if [[ "$USE_GPU" == true ]]; then
            case "$stage" in
                em)       cmd="$cmd -gpu_id $GPU_ID $GPU_EM_FLAGS" ;;
                nvt|npt|md) cmd="$cmd -gpu_id $GPU_ID $GPU_MD_FLAGS" ;;
            esac
        fi
        echo "$cmd"
    }

    # Run mdrun with GPU→CPU fallback (replaces raw $gmx_cmd || true)
    _try_mdrun() {
        local stage="$1" deffnm="$2" log_dir="$3"
        export OMP_NUM_THREADS="$NTHREADS"
        local gmx_cmd
        gmx_cmd=$(_build_mdrun_cmd "$stage" "$deffnm")
        if $gmx_cmd 2>&1 | tee "${log_dir}/mdrun_${deffnm}.log"; then
            return 0
        fi
        if [[ "$USE_GPU" == true ]]; then
            log_warn "GPU mdrun failed for $deffnm, retrying CPU-only..."
            $GMX_BIN mdrun -v -deffnm "$deffnm" -ntmpi 1 -ntomp "$NTHREADS" 2>&1 \
                | tee "${log_dir}/mdrun_${deffnm}_cpu.log"
            return $?
        fi
        return 1
    }

    _process_structure() {
        local pdb_input="$1" name="$2" outdir="$3"

        log_section "Processing: $name"

        if [[ -f "$outdir/md.gro" && -f "$outdir/md.edr" ]]; then
            log "  ✓ Simulation complete, skipping"
            return 0
        fi

        setup_output_dirs "$outdir"
        pushd "$outdir" > /dev/null

        if ! python3 -m gromacs_utils.cli prepare-structure "$pdb_input" -o clean.pdb; then
            log_error "Failed to prepare structure: $pdb_input"
            popd > /dev/null
            return 1
        fi

        log "Generating topology..."
        echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 | tee logs/pdb2gmx.log

        log "Creating index file..."
        echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 | tee logs/make_ndx.log
        python3 -m gromacs_utils.cli chain-index --pdb clean.pdb --gro protein.gro --index index.ndx

        log "Setting up simulation box..."
        $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 | tee logs/editconf.log
        $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tee logs/solvate.log

        log "Generating MDP files..."
        python3 -m gromacs_utils.cli generate-mdp all -o . \
            --em-steps $EM_STEPS --nvt-steps $NVT_STEPS --npt-steps $NPT_STEPS --md-steps $MD_STEPS

        log "Adding ions..."
        $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 | tee logs/grompp_ions.log
        echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc $ION_CONCENTRATION 2>&1 | tee logs/genion.log

        log "Running energy minimization..."
        $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 | tee logs/grompp_em.log
        _try_mdrun "em" "em" "logs"
        if [[ ! -f em.gro ]]; then
            log_error "EM failed — em.gro not produced for $name"
            popd > /dev/null
            return 1
        fi
        echo "Potential" | $GMX_BIN energy -f em.edr -o em_potential.xvg 2>&1 | tee logs/energy_em.log

        log "Running NVT equilibration..."
        $GMX_BIN grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 2>&1 | tee logs/grompp_nvt.log
        _try_mdrun "nvt" "nvt" "logs"
        if [[ ! -f nvt.gro ]]; then
            log_error "NVT equilibration failed — nvt.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        log "Running NPT equilibration..."
        $GMX_BIN grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr 2>&1 | tee logs/grompp_npt.log
        _try_mdrun "npt" "npt" "logs"
        if [[ ! -f npt.gro ]]; then
            log_error "NPT equilibration failed — npt.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        log "Running production MD..."
        $GMX_BIN grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr 2>&1 | tee logs/grompp_md.log
        _try_mdrun "md" "md" "logs"
        if [[ ! -f md.gro ]]; then
            log_error "Production MD failed — md.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        popd > /dev/null
        log "✓ MD simulation completed for $name"
    }

    _calculate_interaction_energy() {
        local outdir="$1" name="$2"
        [[ -f "$outdir/ie.edr" || -f "$outdir/interaction_energy.xvg" ]] && { log "  ✓ IE already calculated"; return 0; }
        [[ ! -f "$outdir/md.gro" || ! -f "$outdir/md.xtc" ]] && { log_warn "  No MD data, skipping IE"; return 1; }

        log "Calculating Interaction Energy: $name"
        pushd "$outdir" > /dev/null
        [[ -f ie.mdp ]] || python3 -m gromacs_utils.cli generate-mdp ie -o ie.mdp
        $GMX_BIN grompp -f ie.mdp -c md.gro -p topol.top -n index.ndx -o ie.tpr 2>&1 | tee logs/grompp_ie.log || true
        # IE rerun: use GPU for non-bonded forces (nb); PME/bonded stay CPU (rerun limitation)
        local ie_gpu_flags=""
        if [[ "$USE_GPU" == true ]]; then
            ie_gpu_flags="-gpu_id $GPU_ID -nb gpu"
        fi
        $GMX_BIN mdrun -s ie.tpr -rerun md.xtc -e ie.edr $ie_gpu_flags 2>&1 | tee logs/mdrun_ie.log || true
        echo -e "Coul-SR:ChainA-ChainB\nLJ-SR:ChainA-ChainB\n\n" | $GMX_BIN energy -f ie.edr -o interaction_energy.xvg 2>&1 | tee logs/energy_ie.log || true
        popd > /dev/null
    }

    _generate_trajectory_outputs() {
        local outdir="$1" name="$2"
        [[ -f "$outdir/statistics/md_statistics.json" ]] && { log "  ✓ Trajectory outputs exist"; return 0; }
        [[ ! -f "$outdir/md.tpr" ]] && { log_warn "  No MD data, skipping trajectory outputs"; return 1; }

        log "Generating trajectory outputs: $name"
        pushd "$outdir" > /dev/null
        echo "Protein Protein" | $GMX_BIN trjconv -s md.tpr -f md.xtc -o trajectories/md_center.xtc -center -pbc mol 2>&1 | tee logs/trjconv.log || true
        echo "Protein" | $GMX_BIN trjconv -s md.tpr -f md.gro -o structures/final_structure.pdb >> logs/trjconv.log 2>&1 || true
        echo "Backbone Backbone" | $GMX_BIN rms -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsd.xvg -tu ps > logs/rmsd.log 2>&1 || true
        echo "Backbone" | $GMX_BIN rmsf -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsf.xvg -res > logs/rmsf.log 2>&1 || true
        echo "Protein" | $GMX_BIN gyrate -s md.tpr -f trajectories/md_center.xtc -o analysis/gyrate.xvg > logs/gyrate.log 2>&1 || true
        echo "Protein" | $GMX_BIN sasa -s md.tpr -f trajectories/md_center.xtc -o analysis/sasa.xvg > logs/sasa.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN hbond -s md.tpr -f trajectories/md_center.xtc -n index.ndx -num analysis/hbonds.xvg > logs/hbond.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN mindist -s md.tpr -f trajectories/md_center.xtc -n index.ndx -od analysis/mindist.xvg > logs/mindist.log 2>&1 || true
        python3 -m gromacs_utils.md_statistics --workdir "$outdir"
        popd > /dev/null
    }

    _generate_structure_visualization() {
        local outdir="$1" name="$2"
        [[ -f "$outdir/visualization/visualize_trajectory.pml" || -f "$outdir/visualization/visualize_interface.pml" ]] && { log "  ✓ Visualization exists"; return 0; }
        [[ ! -d "$outdir/analysis" ]] && { log_warn "  No analysis dir, skipping viz"; return 1; }

        log "Generating visualization for $name..."
        pushd "$outdir" > /dev/null
        generate_visualization "$outdir" all
        run_gnuplot_scripts plots
        python3 -m gromacs_utils.cli generate-plots md -i "$outdir/analysis" -o "$outdir/plots" 2>&1 || true
        popd > /dev/null
    }

    _compare_chain_results() {
        log_section "Comparing Stability Results"
        local num_structs=${#PDB_LIST[@]}

        local summary_file="$WORKDIR/structures_summary.txt"
        {
            echo "=============================================="
            echo "MULTI-STRUCTURE STABILITY COMPARISON SUMMARY"
            echo "=============================================="
            echo ""
            echo "Structures analyzed: $num_structs"
            echo "Generated: $(date)"
            echo ""

            for i in "${!PDB_LIST[@]}"; do
                local struct_dir="$WORKDIR/structure_$((i+1))"
                local pdb_name=$(basename "${PDB_LIST[$i]}" .pdb)
                echo "--- Structure $((i+1)): $pdb_name ---"
                if [[ -f "$struct_dir/statistics/md_statistics.json" ]]; then
                    python3 -c "
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
for key, val in data.items():
    if isinstance(val, (int, float)):
        print(f'  {key}: {val:.4f}' if isinstance(val, float) else f'  {key}: {val}')
    elif isinstance(val, dict):
        for k, v in val.items():
            if isinstance(v, (int, float)):
                print(f'  {key}.{k}: {v:.4f}' if isinstance(v, float) else f'  {key}.{k}: {v}')
" "$struct_dir/statistics/md_statistics.json" 2>/dev/null || echo "  (statistics unavailable)"
                else
                    echo "  (no statistics file found)"
                fi
                echo ""
            done
        } > "$summary_file"
        log "Summary: $summary_file"

        if [[ $num_structs -ge 2 ]]; then
            for ((i=0; i<num_structs-1; i++)); do
                for ((j=i+1; j<num_structs; j++)); do
                    local s1="structure_$((i+1))" s2="structure_$((j+1))"
                    if [[ -d "$WORKDIR/$s1" && -d "$WORKDIR/$s2" ]]; then
                        python3 -m gromacs_utils.md_results_comparator \
                            --workdir "$WORKDIR" \
                            --pdb1 "${PDB_LIST[$i]}" --pdb2 "${PDB_LIST[$j]}" \
                            --struct1-dir "$s1" --struct2-dir "$s2" \
                            --output "comparison_$((i+1))_vs_$((j+1)).txt" 2>&1 || true
                    fi
                done
            done

            if [[ $num_structs -gt 2 ]]; then
                local pdb_args="" struct_args=""
                for i in "${!PDB_LIST[@]}"; do
                    pdb_args="$pdb_args ${PDB_LIST[$i]}"
                    struct_args="$struct_args structure_$((i+1))"
                done
                python3 -m gromacs_utils.md_results_comparator \
                    --workdir "$WORKDIR" --multi \
                    --pdb-list $pdb_args --struct-dirs $struct_args 2>&1 || true
            fi
        fi
    }

    # --- Execute based on mode ---
    check_gromacs || return 1
    check_python_modules || return 1
    local STRUCT_FAIL_COUNT=0

    case $MODE in
        full)
            _check_gpu
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local failed=false
                _process_structure "${PDB_LIST[$i]}" "$sn" "$WORKDIR/$sn" || failed=true
                _calculate_interaction_energy "$WORKDIR/$sn" "$sn" || failed=true
                _generate_trajectory_outputs "$WORKDIR/$sn" "$sn" || failed=true
                _generate_structure_visualization "$WORKDIR/$sn" "$sn" || failed=true
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            _compare_chain_results
            ;;
        sim-only)
            _check_gpu
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local failed=false
                _process_structure "${PDB_LIST[$i]}" "$sn" "$WORKDIR/$sn" || failed=true
                _calculate_interaction_energy "$WORKDIR/$sn" "$sn" || failed=true
                _generate_trajectory_outputs "$WORKDIR/$sn" "$sn" || failed=true
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            _compare_chain_results
            ;;
        viz-only)
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                [[ ! -d "$WORKDIR/$sn" ]] && { log_error "No data for $sn"; return 1; }
                _generate_structure_visualization "$WORKDIR/$sn" "$sn"
            done
            _compare_chain_results
            ;;
        compare-only)
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local sd="$WORKDIR/$sn"
                if [[ ! -d "$sd" ]]; then
                    STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
                    log_warn "No simulation data for $sn, skipping"
                    continue
                fi
                local failed=false
                _calculate_interaction_energy "$sd" "$sn" || failed=true
                _generate_trajectory_outputs "$sd" "$sn" || failed=true
                _generate_structure_visualization "$sd" "$sn" || failed=true
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            _compare_chain_results
            ;;
    esac

    if [[ $STRUCT_FAIL_COUNT -gt 0 ]]; then
        log_error "compare_chain_stability completed with failures: $STRUCT_FAIL_COUNT/${#PDB_LIST[@]} structure(s) incomplete"
        return 1
    fi

    return 0
}

#------------------------------------------------------------------------------
# STEP 3: INTERFACE ANALYSIS
#------------------------------------------------------------------------------

run_interface_analysis() {
    local step_dir="$1"

    local OVERRIDE_EXISTING=$(toml_get "step.interface_analysis.override_existing" "false")
    local BOX_DISTANCE=$(toml_get "step.interface_analysis.box_distance" "1.0")
    local EM_STEPS=$(toml_get "step.interface_analysis.em_steps" "5000")
    local GPU_ID=$(toml_get "step.interface_analysis.gpu_id" "0")
    local MAX_THREADS=$(toml_get "step.interface_analysis.max_threads" "0")
    NTHREADS=$(( MAX_THREADS > 0 ? MAX_THREADS : $(nproc) ))

    local PDB_FILES=()
    load_step_structures "interface_analysis" PDB_FILES

    if [[ ${#PDB_FILES[@]} -eq 0 ]]; then
        log_error "No structures defined for interface_analysis"
        return 1
    fi

    check_gromacs || return 1
    check_python_modules || return 1

    log "Found ${#PDB_FILES[@]} PDB files to process:"
    for pdb in "${PDB_FILES[@]}"; do
        log "  - $(basename "$pdb")"
    done
    echo ""

    local PROCESSED=0

    for PDB_INPUT in "${PDB_FILES[@]}"; do
        local STRUCT_NAME
        STRUCT_NAME=$(shorten_pdb_name "$PDB_INPUT")
        local WORKDIR="${step_dir}/${STRUCT_NAME}"

        log_section "Processing: $STRUCT_NAME"
        log "Output: $WORKDIR"

        if [[ "$OVERRIDE_EXISTING" != "true" && -f "$WORKDIR/interface_metrics.json" ]]; then
            log "  ✓ Already processed, skipping"
            PROCESSED=$((PROCESSED + 1))
            continue
        fi

        [[ "$OVERRIDE_EXISTING" == "true" && -d "$WORKDIR" ]] && rm -rf "$WORKDIR"

        mkdir -p "$WORKDIR"
        pushd "$WORKDIR" > /dev/null
        setup_output_structure "$WORKDIR"

        local START_T=$(date +%s)

        # Prepare structure
        log "Preparing structure..."
        if ! python3 -m gromacs_utils.cli prepare-structure "$PDB_INPUT" -o clean.pdb; then
            log_error "Failed to prepare structure: $PDB_INPUT"
            popd > /dev/null
            continue
        fi
        log "  Generating topology..."
        echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > logs/pdb2gmx.log 2>&1
        echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx > logs/make_ndx.log 2>&1
        create_chain_index clean.pdb protein.gro index.ndx

        log "  Solvating..."
        $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron > logs/editconf.log 2>&1
        $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > logs/solvate.log 2>&1

        python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 100.0

        log "  Running energy minimization..."
        $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 2 > logs/grompp.log 2>&1
        run_em em logs

        if [[ ! -f em.gro ]]; then
            log_error "EM failed — em.gro not produced for $STRUCT_NAME"
            popd > /dev/null
            continue
        fi

        # GROMACS analysis
        log "Running GROMACS analysis tools..."
        echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx -od analysis/mindist.xvg -on analysis/numcont.xvg -d 0.6 > logs/mindist.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx -num analysis/hbonds.xvg > logs/hbond.log 2>&1 || true
        echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o analysis/sasa_complex.xvg > logs/sasa_complex.log 2>&1 || true
        echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o analysis/sasa_chainA.xvg > logs/sasa_A.log 2>&1 || true
        echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o analysis/sasa_chainB.xvg > logs/sasa_B.log 2>&1 || true
        echo -e "Potential\nCoul-SR\nLJ-SR\n\n" | $GMX_BIN energy -f em.edr -o analysis/energies.xvg > logs/energy.log 2>&1 || true

        # Contact map and interface report
        log "Generating contact map..."
        python3 -m gromacs_utils.contact_map --workdir "$WORKDIR" --gro em.gro --pdb clean.pdb
        log "Generating interface report..."
        python3 -m gromacs_utils.interface_analyzer --workdir "$WORKDIR" --name "$STRUCT_NAME" --output "$WORKDIR"

        generate_visualization "$WORKDIR" interface
        echo "Protein" | $GMX_BIN trjconv -s em.tpr -f em.gro -o structures/structure_minimized.pdb > logs/trjconv.log 2>&1 || true

        log "Generating plots..."
        run_gnuplot_scripts plots
        python3 -m gromacs_utils.cli generate-plots contact -i analysis/contact_map.txt -o plots/contact_map_heatmap.png 2>&1 || true
        python3 -c "
from gromacs_utils.plotting import plot_focused_contact_map
plot_focused_contact_map('analysis/interface_residues.txt', 'plots/contact_map_focused.png', max_pairs=50)
" 2>&1 || true

        popd > /dev/null

        local END_T=$(date +%s)
        log "  ✓ Completed in $((END_T - START_T)) seconds"
        PROCESSED=$((PROCESSED + 1))
        echo ""
    done

    log_section "Interface Analysis Complete"
    log "Processed: $PROCESSED / ${#PDB_FILES[@]}"
}

#------------------------------------------------------------------------------
# STEP 4: BATCH COMPARISON
#------------------------------------------------------------------------------

run_batch_comparison() {
    local step_dir="$1"

    local OVERRIDE_EXISTING=$(toml_get "step.batch_comparison.override_existing" "false")
    local EM_STEPS=$(toml_get "step.batch_comparison.em_steps" "2000")
    local BOX_DISTANCE=$(toml_get "step.batch_comparison.box_distance" "0.8")
    local BOX_TYPE=$(toml_get "step.batch_comparison.box_type" "cubic")
    local GPU_ID=$(toml_get "step.batch_comparison.gpu_id" "0")
    local MAX_THREADS=$(toml_get "step.batch_comparison.max_threads" "0")
    NTHREADS=$(( MAX_THREADS > 0 ? MAX_THREADS : $(nproc) ))

    local DATASETS=()
    while IFS= read -r ds; do
        [[ -n "$ds" ]] && DATASETS+=("$ds")
    done < <(toml_get_array "step.batch_comparison.datasets.active")

    if [[ ${#DATASETS[@]} -eq 0 ]]; then
        log_error "No datasets enabled for batch_comparison"
        return 1
    fi

    check_gromacs || return 1
    check_python_modules || return 1

    log "Datasets: ${DATASETS[*]}"
    echo ""

    for dataset_name in "${DATASETS[@]}"; do
        local input_dir="${INPUT_BASE}/${dataset_name}"
        local output_dir="${step_dir}/${dataset_name}"

        log_section "Processing Dataset: $dataset_name"
        log "Input:  $input_dir"
        log "Output: $output_dir"

        if [[ ! -d "$input_dir" ]]; then
            log_warn "Input directory does not exist: $input_dir"
            continue
        fi

        mkdir -p "$output_dir"
        cd "$output_dir"

        mapfile -t STRUCTURES < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.pdb" -o -name "*.cif" \) | sort)

        if [[ ${#STRUCTURES[@]} -eq 0 ]]; then
            log "No structures found in $input_dir"
            continue
        fi

        log "Found ${#STRUCTURES[@]} structures"
        echo ""

        local ds_start=$(date +%s)

        for struct_file in "${STRUCTURES[@]}"; do
            local struct_name=$(basename "$struct_file" | sed 's/\.[^.]*$//')
            local struct_dir="$output_dir/$struct_name"

            log "Analyzing: $struct_name"

            if [[ "$OVERRIDE_EXISTING" != "true" && -f "$struct_dir/metrics.json" && -f "$struct_dir/em.gro" ]]; then
                log "  ✓ Already processed, skipping"
                echo "OK" > "$struct_dir/status.txt"
                continue
            fi

            [[ "$OVERRIDE_EXISTING" == "true" && -d "$struct_dir" ]] && rm -rf "$struct_dir"

            mkdir -p "$struct_dir"
            cd "$struct_dir"

            if ! python3 -m gromacs_utils.cli prepare-structure "$struct_file" -o clean.pdb; then
                log_error "Failed to prepare structure: $struct_file"
                echo "ERROR" > status.txt
                continue
            fi

            echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > pdb2gmx.log 2>&1 || {
                echo "ERROR" > status.txt
                log "  ✗ Failed"
                echo ""
                continue
            }

            echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx > make_ndx.log 2>&1
            $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt $BOX_TYPE > editconf.log 2>&1
            $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > solvate.log 2>&1

            python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 500.0

            $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 3 > grompp.log 2>&1
            run_em em .

            if [[ ! -f em.gro ]]; then
                log_error "EM failed — em.gro not produced for $struct_name"
                echo "ERROR" > status.txt
                continue
            fi

            echo -e "Potential\n\n" | $GMX_BIN energy -f em.edr -o energy.xvg > energy.log 2>&1 || true
            echo "Protein" | $GMX_BIN gyrate -s em.tpr -f em.gro -o gyrate.xvg > gyrate.log 2>&1 || true
            echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg > sasa.log 2>&1 || true

            python3 -m gromacs_utils.cli extract-metrics --workdir . -o metrics.json
            echo "OK" > status.txt
            log "  ✓ Complete"
            echo ""
        done

        # Generate comparison report
        log "Generating comparison report for $dataset_name..."
        cd "$output_dir"
        python3 -m gromacs_utils.cli batch-analysis --workdir "$output_dir" --plots
        run_gnuplot_scripts "$output_dir/plots"

        local ds_end=$(date +%s)
        log_section "$dataset_name Complete"
        log "Time: $((ds_end - ds_start))s | Structures: ${#STRUCTURES[@]}"
        log "Output: $output_dir"
        head -15 "$output_dir/comparison_report.txt" 2>/dev/null || true
        echo ""
    done
}

#------------------------------------------------------------------------------
# STEP 5: PRODUCTION MD (Full GPU Offload)
#
# gmx mdrun -deffnm md -v -maxh <max_hours> -ntmpi 1
#           -nb gpu -pme gpu -bonded gpu -update gpu
#
# Flags:
#   -deffnm md    : all I/O files use "md" prefix (md.tpr -> md.xtc, md.edr, md.log)
#   -v            : verbose – prints step number, time, performance to terminal
#   -maxh         : maximum wall-clock runtime in hours; checkpoints and stops at limit
#   -ntmpi 1      : single thread-MPI rank (one rank per GPU)
#   -nb gpu       : non-bonded interactions (Coulomb + van der Waals) on GPU
#   -pme gpu      : PME electrostatics (long-range Coulomb via Particle Mesh Ewald) on GPU
#   -bonded gpu   : bonded interactions (bonds, angles, dihedrals) on GPU
#   -update gpu   : coordinate update and constraints (LINCS/SETTLE) on GPU
#------------------------------------------------------------------------------

run_production_md() {
    local step_dir="$1"

    # Load settings
    local OVERRIDE_EXISTING=$(toml_get "step.production_md.override_existing" "false")
    local MODE="${CLI_MODE:-$(toml_get "step.production_md.mode" "full")}"
    local BOX_DISTANCE=$(toml_get "step.production_md.box_distance" "1.2")
    local ION_CONCENTRATION=$(toml_get "step.production_md.ion_concentration" "0.15")
    local EM_STEPS=$(toml_get "step.production_md.em_steps" "5000")
    local NVT_STEPS=$(toml_get "step.production_md.nvt_steps" "50000")
    local NPT_STEPS=$(toml_get "step.production_md.npt_steps" "50000")
    local MD_STEPS=$(toml_get "step.production_md.md_steps" "250000")
    local GPU_ID=$(toml_get "step.production_md.gpu_id" "0")
    local MAX_THREADS=$(toml_get "step.production_md.max_threads" "0")
    local MAX_HOURS=$(toml_get "step.production_md.max_hours" "5000")
    NTHREADS=$(( MAX_THREADS > 0 ? MAX_THREADS : $(nproc) ))

    # Full GPU offload flags (all four components on GPU)
    local GPU_FULL_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu"
    # EM: -pme cpu required (steep integrator limitation); -bonded gpu OK
    local GPU_EM_FLAGS="-nb gpu -pme cpu -bonded gpu"

    # Load structures
    local PDB_LIST=()
    load_step_structures "production_md" PDB_LIST

    if [[ ${#PDB_LIST[@]} -eq 0 ]]; then
        log_error "production_md requires at least 1 structure"
        return 1
    fi

    # Generate workdir name
    local names=()
    for pdb in "${PDB_LIST[@]}"; do
        names+=("$(shorten_pdb_name "$pdb")")
    done
    local joined
    joined=$(make_workdir_name names)
    local WORKDIR="${step_dir}/${joined}"

    log "Mode: $MODE"
    log "Max wall-clock hours: $MAX_HOURS"
    log "GPU flags (MD): $GPU_FULL_FLAGS"
    log "Structures: ${#PDB_LIST[@]}"
    for i in "${!PDB_LIST[@]}"; do
        log "  Structure $((i+1)): $(basename "${PDB_LIST[$i]}")"
    done
    log "Output: $WORKDIR"

    mkdir -p "$WORKDIR/logs"
    write_manifest "$WORKDIR" "${PDB_LIST[@]}"
    cd "$WORKDIR"

    # GPU detection — use centralized detection from gromacs_common.sh
    USE_GPU=true
    export USE_GPU
    local GPU_INFO=""

    _check_gpu_prod() {
        if _detect_gpu_available; then
            GPU_INFO=$( (nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || rocm-smi --showproductname 2>/dev/null) | head -1 )
            log "  Found GPU: ${GPU_INFO:-detected}"
        else
            log_warn "GPU not usable by GROMACS, will use CPU mode"
            USE_GPU=false
            export USE_GPU
        fi
        return 0
    }

    _build_mdrun_cmd_prod() {
        local stage="$1" deffnm="$2"
        local cmd="$GMX_BIN mdrun -deffnm $deffnm -v -ntmpi 1 -ntomp $NTHREADS"
        if [[ "$USE_GPU" == true ]]; then
            case "$stage" in
                em)
                    cmd="$cmd -gpu_id $GPU_ID $GPU_EM_FLAGS"
                    ;;
                nvt|npt)
                    cmd="$cmd -gpu_id $GPU_ID $GPU_FULL_FLAGS"
                    ;;
                md)
                    # Production MD: full GPU offload + wall-clock limit
                    cmd="$cmd -maxh $MAX_HOURS -gpu_id $GPU_ID $GPU_FULL_FLAGS"
                    ;;
            esac
        fi
        echo "$cmd"
    }

    # Run mdrun with GPU→CPU fallback for production MD
    _try_mdrun_prod() {
        local stage="$1" deffnm="$2" log_dir="$3"
        export OMP_NUM_THREADS="$NTHREADS"
        local gmx_cmd
        gmx_cmd=$(_build_mdrun_cmd_prod "$stage" "$deffnm")
        log "  CMD: $gmx_cmd"
        if $gmx_cmd 2>&1 | tee "${log_dir}/mdrun_${deffnm}.log"; then
            return 0
        fi
        if [[ "$USE_GPU" == true ]]; then
            log_warn "GPU mdrun failed for $deffnm, retrying CPU-only..."
            local cpu_cmd="$GMX_BIN mdrun -deffnm $deffnm -v -ntmpi 1 -ntomp $NTHREADS"
            [[ "$stage" == "md" ]] && cpu_cmd="$cpu_cmd -maxh $MAX_HOURS"
            $cpu_cmd 2>&1 | tee "${log_dir}/mdrun_${deffnm}_cpu.log"
            return $?
        fi
        return 1
    }

    _process_structure_prod() {
        local pdb_input="$1" name="$2" outdir="$3"

        log_section "Processing: $name"

        if [[ -f "$outdir/md.gro" && -f "$outdir/md.edr" && "$OVERRIDE_EXISTING" != "true" ]]; then
            log "  ✓ Simulation complete, skipping"
            return 0
        fi

        [[ "$OVERRIDE_EXISTING" == "true" && -d "$outdir" ]] && rm -rf "$outdir"

        setup_output_dirs "$outdir"
        pushd "$outdir" > /dev/null

        if ! python3 -m gromacs_utils.cli prepare-structure "$pdb_input" -o clean.pdb; then
            log_error "Failed to prepare structure: $pdb_input"
            popd > /dev/null
            return 1
        fi

        log "Generating topology..."
        echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 | tee logs/pdb2gmx.log

        log "Creating index file..."
        echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 | tee logs/make_ndx.log
        python3 -m gromacs_utils.cli chain-index --pdb clean.pdb --gro protein.gro --index index.ndx

        log "Setting up simulation box..."
        $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 | tee logs/editconf.log
        $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tee logs/solvate.log

        log "Generating MDP files..."
        python3 -m gromacs_utils.cli generate-mdp all -o . \
            --em-steps $EM_STEPS --nvt-steps $NVT_STEPS --npt-steps $NPT_STEPS --md-steps $MD_STEPS

        log "Adding ions..."
        $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 | tee logs/grompp_ions.log
        echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc $ION_CONCENTRATION 2>&1 | tee logs/genion.log

        log "Running energy minimization..."
        $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 | tee logs/grompp_em.log
        _try_mdrun_prod "em" "em" "logs"
        if [[ ! -f em.gro ]]; then
            log_error "EM failed — em.gro not produced for $name"
            popd > /dev/null
            return 1
        fi
        echo "Potential" | $GMX_BIN energy -f em.edr -o em_potential.xvg 2>&1 | tee logs/energy_em.log

        log "Running NVT equilibration..."
        $GMX_BIN grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 2>&1 | tee logs/grompp_nvt.log
        _try_mdrun_prod "nvt" "nvt" "logs"
        if [[ ! -f nvt.gro ]]; then
            log_error "NVT equilibration failed — nvt.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        log "Running NPT equilibration..."
        $GMX_BIN grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr 2>&1 | tee logs/grompp_npt.log
        _try_mdrun_prod "npt" "npt" "logs"
        if [[ ! -f npt.gro ]]; then
            log_error "NPT equilibration failed — npt.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        log "Running production MD (full GPU offload, maxh=${MAX_HOURS}h)..."
        $GMX_BIN grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr 2>&1 | tee logs/grompp_md.log
        _try_mdrun_prod "md" "md" "logs"
        if [[ ! -f md.gro ]]; then
            log_error "Production MD failed — md.gro not produced for $name"
            popd > /dev/null
            return 1
        fi

        popd > /dev/null
        log "✓ Production MD completed for $name"
    }

    _generate_prod_trajectory_outputs() {
        local outdir="$1" name="$2"
        [[ -f "$outdir/statistics/md_statistics.json" ]] && { log "  ✓ Trajectory outputs exist"; return 0; }
        [[ ! -f "$outdir/md.tpr" ]] && { log_warn "  No MD data, skipping trajectory outputs"; return 1; }

        log "Generating trajectory outputs: $name"
        pushd "$outdir" > /dev/null
        echo "Protein Protein" | $GMX_BIN trjconv -s md.tpr -f md.xtc -o trajectories/md_center.xtc -center -pbc mol 2>&1 | tee logs/trjconv.log || true
        echo "Protein" | $GMX_BIN trjconv -s md.tpr -f md.gro -o structures/final_structure.pdb >> logs/trjconv.log 2>&1 || true
        echo "Backbone Backbone" | $GMX_BIN rms -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsd.xvg -tu ps > logs/rmsd.log 2>&1 || true
        echo "Backbone" | $GMX_BIN rmsf -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsf.xvg -res > logs/rmsf.log 2>&1 || true
        echo "Protein" | $GMX_BIN gyrate -s md.tpr -f trajectories/md_center.xtc -o analysis/gyrate.xvg > logs/gyrate.log 2>&1 || true
        echo "Protein" | $GMX_BIN sasa -s md.tpr -f trajectories/md_center.xtc -o analysis/sasa.xvg > logs/sasa.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN hbond -s md.tpr -f trajectories/md_center.xtc -n index.ndx -num analysis/hbonds.xvg > logs/hbond.log 2>&1 || true
        echo -e "ChainA\nChainB" | $GMX_BIN mindist -s md.tpr -f trajectories/md_center.xtc -n index.ndx -od analysis/mindist.xvg > logs/mindist.log 2>&1 || true
        python3 -m gromacs_utils.md_statistics --workdir "$outdir"
        popd > /dev/null
    }

    _generate_prod_visualization() {
        local outdir="$1" name="$2"
        [[ -f "$outdir/visualization/visualize_trajectory.pml" || -f "$outdir/visualization/visualize_interface.pml" ]] && { log "  ✓ Visualization exists"; return 0; }
        [[ ! -d "$outdir/analysis" ]] && { log_warn "  No analysis dir, skipping viz"; return 1; }

        log "Generating visualization for $name..."
        pushd "$outdir" > /dev/null
        generate_visualization "$outdir" all
        run_gnuplot_scripts plots
        python3 -m gromacs_utils.cli generate-plots md -i "$outdir/analysis" -o "$outdir/plots" 2>&1 || true
        popd > /dev/null
    }

    # --- Execute based on mode ---
    check_gromacs || return 1
    check_python_modules || return 1
    local STRUCT_FAIL_COUNT=0

    case $MODE in
        full)
            _check_gpu_prod
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local failed=false
                _process_structure_prod "${PDB_LIST[$i]}" "$sn" "$WORKDIR/$sn" || failed=true
                _generate_prod_trajectory_outputs "$WORKDIR/$sn" "$sn" || failed=true
                _generate_prod_visualization "$WORKDIR/$sn" "$sn" || failed=true
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            ;;
        sim-only)
            _check_gpu_prod
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local failed=false
                _process_structure_prod "${PDB_LIST[$i]}" "$sn" "$WORKDIR/$sn" || failed=true
                _generate_prod_trajectory_outputs "$WORKDIR/$sn" "$sn" || failed=true
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            ;;
        md-only)
            # Only run the production MD step (assumes em/nvt/npt already done)
            _check_gpu_prod
            for i in "${!PDB_LIST[@]}"; do
                local sn="structure_$((i+1))"
                local outdir="$WORKDIR/$sn"
                local failed=false
                if [[ ! -f "$outdir/npt.gro" || ! -f "$outdir/npt.cpt" ]]; then
                    log_error "No NPT checkpoint found for $sn — run full mode first"
                    STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
                    continue
                fi
                pushd "$outdir" > /dev/null
                log "Running production MD only for $sn (full GPU offload, maxh=${MAX_HOURS}h)..."
                $GMX_BIN grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr 2>&1 | tee logs/grompp_md.log
                _try_mdrun_prod "md" "md" "logs"
                if [[ ! -f md.gro ]]; then
                    log_error "Production MD failed — md.gro not produced for $sn"
                    failed=true
                    popd > /dev/null
                else
                    popd > /dev/null
                    _generate_prod_trajectory_outputs "$outdir" "$sn" || failed=true
                fi
                [[ "$failed" == true ]] && STRUCT_FAIL_COUNT=$((STRUCT_FAIL_COUNT + 1))
            done
            ;;
    esac

    # Summary
    local summary_file="$WORKDIR/production_md_summary.txt"
    {
        echo "=============================================="
        echo "PRODUCTION MD SIMULATION SUMMARY"
        echo "=============================================="
        echo ""
        echo "Structures analyzed: ${#PDB_LIST[@]}"
        echo "Failed structures:   $STRUCT_FAIL_COUNT"
        echo "Mode: $MODE"
        echo "Max wall-clock hours: $MAX_HOURS"
        echo "GPU flags: $GPU_FULL_FLAGS"
        echo "Generated: $(date)"
        echo ""

        for i in "${!PDB_LIST[@]}"; do
            local struct_dir="$WORKDIR/structure_$((i+1))"
            local pdb_name=$(basename "${PDB_LIST[$i]}" .pdb)
            echo "--- Structure $((i+1)): $pdb_name ---"
            if [[ -f "$struct_dir/statistics/md_statistics.json" ]]; then
                python3 -c "
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
for key, val in data.items():
    if isinstance(val, (int, float)):
        print(f'  {key}: {val:.4f}' if isinstance(val, float) else f'  {key}: {val}')
    elif isinstance(val, dict):
        for k, v in val.items():
            if isinstance(v, (int, float)):
                print(f'  {key}.{k}: {v:.4f}' if isinstance(v, float) else f'  {key}.{k}: {v}')
" "$struct_dir/statistics/md_statistics.json" 2>/dev/null || echo "  (statistics unavailable)"
            elif [[ -f "$struct_dir/md.gro" ]]; then
                echo "  MD completed (no statistics generated yet)"
            else
                echo "  (simulation not completed)"
            fi
            echo ""
        done
    } > "$summary_file"
    log "Summary: $summary_file"

    if [[ $STRUCT_FAIL_COUNT -gt 0 ]]; then
        log_error "production_md completed with failures: $STRUCT_FAIL_COUNT/${#PDB_LIST[@]} structure(s) incomplete"
        return 1
    fi

    return 0
}

#------------------------------------------------------------------------------
# STEP 6: VISUALIZE RESULTS
#
# Post-processing step — scans existing simulation outputs (steps 1-5) and
# (re)generates matplotlib plots, PyMOL/VMD/ChimeraX scripts, gnuplot renders,
# contact maps, and statistics.  No GROMACS runs or GPU required.
#------------------------------------------------------------------------------

run_visualize_results() {
    local step_dir="$1"

    # Load settings
    local OVERRIDE_EXISTING=$(toml_get "step.visualize_results.override_existing" "false")
    local DPI=$(toml_get "step.visualize_results.dpi" "150")
    local PLOT_FORMAT=$(toml_get "step.visualize_results.plot_format" "png")
    local MAX_CONTACT_PAIRS=$(toml_get "step.visualize_results.max_contact_pairs" "50")
    local GEN_PYMOL=$(toml_get "step.visualize_results.generate_pymol" "true")
    local GEN_VMD=$(toml_get "step.visualize_results.generate_vmd" "true")
    local GEN_CHIMERAX=$(toml_get "step.visualize_results.generate_chimerax" "true")

    check_python_modules || return 1

    # Determine which source steps to process
    local SOURCE_STEPS=()
    while IFS= read -r s; do
        [[ -n "$s" ]] && SOURCE_STEPS+=("$s")
    done < <(toml_get_array "step.visualize_results.source_steps.active")

    # Default: scan all known step prefixes that exist in the run directory
    if [[ ${#SOURCE_STEPS[@]} -eq 0 ]]; then
        for name in quick_stability compare_chain_stability interface_analysis batch_comparison production_md; do
            local prefix="${STEP_PREFIX[$name]:-$name}"
            [[ -d "$RUN_DIR/$prefix" ]] && SOURCE_STEPS+=("$name")
        done
    fi

    if [[ ${#SOURCE_STEPS[@]} -eq 0 ]]; then
        log_warn "No prior step outputs found in $RUN_DIR — nothing to visualize"
        return 0
    fi

    log_section "Visualize Results"
    log "Source steps: ${SOURCE_STEPS[*]}"
    log "Settings: dpi=$DPI format=$PLOT_FORMAT override=$OVERRIDE_EXISTING"
    log "Generators: PyMOL=$GEN_PYMOL VMD=$GEN_VMD ChimeraX=$GEN_CHIMERAX"
    echo ""

    local TOTAL_VIZ=0

    # Build viz-type flags based on config
    local VIZ_TYPES=""
    [[ "$GEN_PYMOL" == "true" || "$GEN_VMD" == "true" || "$GEN_CHIMERAX" == "true" ]] && VIZ_TYPES="all"

    for source_step in "${SOURCE_STEPS[@]}"; do
        local prefix="${STEP_PREFIX[$source_step]:-$source_step}"
        local source_dir="$RUN_DIR/$prefix"
        [[ ! -d "$source_dir" ]] && continue

        log "--- Processing outputs from: $source_step ($prefix/) ---"

        # Find structure directories (structure_1, structure_2, ...) or flat workdirs
        local struct_dirs=()
        while IFS= read -r d; do
            struct_dirs+=("$d")
        done < <(find "$source_dir" -maxdepth 1 -type d -name 'structure_*' 2>/dev/null | sort)

        # Also check if source_dir itself is a workdir (batch_comparison layout)
        if [[ ${#struct_dirs[@]} -eq 0 ]]; then
            # Check for analysis/ or em.gro directly in source_dir or its subdirs
            while IFS= read -r d; do
                struct_dirs+=("$d")
            done < <(find "$source_dir" -maxdepth 2 -type d -name 'analysis' -exec dirname {} \; 2>/dev/null | sort -u)
        fi

        if [[ ${#struct_dirs[@]} -eq 0 ]]; then
            log "  No structure workdirs found, skipping"
            continue
        fi

        for workdir in "${struct_dirs[@]}"; do
            local wname=$(basename "$workdir")

            # Skip if viz already exists and override is off
            if [[ "$OVERRIDE_EXISTING" != "true" ]]; then
                local has_plots=false has_viz=false
                [[ -d "$workdir/plots" ]] && [[ $(find "$workdir/plots" -name "*.${PLOT_FORMAT}" 2>/dev/null | head -1) ]] && has_plots=true
                [[ -d "$workdir/visualization" ]] && [[ $(find "$workdir/visualization" -name '*.pml' -o -name '*.vmd' -o -name '*.cxc' 2>/dev/null | head -1) ]] && has_viz=true
                if [[ "$has_plots" == true && "$has_viz" == true ]]; then
                    log "  ✓ $wname — visualization already exists, skipping"
                    TOTAL_VIZ=$((TOTAL_VIZ + 1))
                    continue
                fi
            fi

            log "  Generating visualization for $wname..."
            mkdir -p "$workdir/plots" "$workdir/visualization"

            pushd "$workdir" > /dev/null

            # 3D visualization scripts (PyMOL, VMD, ChimeraX)
            if [[ -n "$VIZ_TYPES" ]]; then
                generate_visualization "$workdir" "$VIZ_TYPES" 2>&1 || log_warn "  Visualization script generation had warnings"
            fi

            # Gnuplot scripts
            run_gnuplot_scripts plots 2>&1 || true

            # MD analysis plots (RMSD, RMSF, Rg, SASA, H-bonds, min-distance)
            if [[ -d "$workdir/analysis" ]]; then
                python3 -m gromacs_utils.cli generate-plots md \
                    -i "$workdir/analysis" -o "$workdir/plots" \
                    --dpi "$DPI" 2>&1 || log_warn "  MD plot generation had warnings"
            fi

            # Contact maps (from interface_analysis outputs)
            if [[ -f "$workdir/analysis/contact_map.txt" ]]; then
                python3 -m gromacs_utils.cli generate-plots contact \
                    -i "analysis/contact_map.txt" \
                    -o "plots/contact_map_heatmap.${PLOT_FORMAT}" \
                    --dpi "$DPI" 2>&1 || true
            fi

            # Focused contact map
            if [[ -f "$workdir/analysis/interface_residues.txt" ]]; then
                python3 -c "
from gromacs_utils.plotting import plot_focused_contact_map
plot_focused_contact_map('analysis/interface_residues.txt', 'plots/contact_map_focused.${PLOT_FORMAT}', max_pairs=${MAX_CONTACT_PAIRS}, dpi=${DPI})
" 2>&1 || true
            fi

            popd > /dev/null
            TOTAL_VIZ=$((TOTAL_VIZ + 1))
        done

        # Batch comparison plots (step 4 generates a batch_results directory)
        if [[ "$source_step" == "batch_comparison" ]]; then
            local batch_results="$source_dir"
            if [[ -f "$batch_results/batch_summary.csv" || -f "$batch_results/statistics/batch_metrics.json" ]]; then
                log "  Generating batch comparison plots..."
                python3 -m gromacs_utils.cli generate-plots batch \
                    -i "$batch_results" -o "$batch_results/plots" \
                    --dpi "$DPI" 2>&1 || log_warn "  Batch plot generation had warnings"
            fi
        fi

        echo ""
    done

    # Copy consolidated plots to the step_dir for easy access
    local consolidated_dir="$step_dir/all_plots"
    mkdir -p "$consolidated_dir"
    for source_step in "${SOURCE_STEPS[@]}"; do
        local prefix="${STEP_PREFIX[$source_step]:-$source_step}"
        local source_dir="$RUN_DIR/$prefix"
        find "$source_dir" -path '*/plots/*.'"$PLOT_FORMAT" -exec cp -n {} "$consolidated_dir/" \; 2>/dev/null || true
    done
    local plot_count=$(find "$consolidated_dir" -name "*.${PLOT_FORMAT}" 2>/dev/null | wc -l)

    log_section "Visualization Complete"
    log "Workdirs processed: $TOTAL_VIZ"
    log "Consolidated plots: $consolidated_dir/ ($plot_count files)"
}

#------------------------------------------------------------------------------
# PIPELINE RUNNER
#------------------------------------------------------------------------------

# Step name -> function dispatch
declare -A STEP_FUNC=(
    [quick_stability]=run_quick_stability
    [compare_chain_stability]=run_compare_chain_stability
    [interface_analysis]=run_interface_analysis
    [batch_comparison]=run_batch_comparison
    [production_md]=run_production_md
    [visualize_results]=run_visualize_results
)

run_pipeline() {
    log_section "GROMACS PPI Analysis Pipeline"
    log "Run directory: $RUN_DIR"
    log "Steps: ${STEPS[*]}"
    log "GROMACS: $($GMX_BIN --version 2>&1 | head -1)"
    echo ""

    if [[ "$DRY_RUN" == true ]]; then
        echo "=== DRY RUN — No simulations will be executed ==="
        echo ""
        for step in "${STEPS[@]}"; do
            local enabled=$(toml_get "step.${step}.enabled" "true")
            local prefix="${STEP_PREFIX[$step]:-$step}"
            echo "  Step: $step"
            echo "    Enabled:    $enabled"
            echo "    Output dir: $RUN_DIR/$prefix"
            echo ""
        done
        echo "Shared structures:"
        for s in "${SHARED_STRUCTURES[@]}"; do
            echo "  - $(basename "$s")"
        done
        exit 0
    fi

    mkdir -p "$RUN_DIR"

    local PIPELINE_START=$(date +%s)
    local STEP_RESULTS=()
    local TOTAL_STEPS=${#STEPS[@]}
    local CURRENT_STEP=0

    for step in "${STEPS[@]}"; do
        CURRENT_STEP=$((CURRENT_STEP + 1))

        # Check if step is enabled
        local enabled=$(toml_get "step.${step}.enabled" "true")
        if [[ "$enabled" != "true" ]]; then
            log "[${CURRENT_STEP}/${TOTAL_STEPS}] SKIP: $step (disabled in config)"
            STEP_RESULTS+=("$step: SKIPPED")
            continue
        fi

        # Check if step function exists
        local func="${STEP_FUNC[$step]:-}"
        if [[ -z "$func" ]]; then
            log_error "Unknown step: $step"
            STEP_RESULTS+=("$step: UNKNOWN")
            [[ "$STOP_ON_ERROR" == "true" ]] && break
            continue
        fi

        local prefix="${STEP_PREFIX[$step]:-$step}"
        local step_dir="$RUN_DIR/$prefix"
        mkdir -p "$step_dir"

        log_section "[${CURRENT_STEP}/${TOTAL_STEPS}] Running: $step"
        local step_start=$(date +%s)

        if $func "$step_dir"; then
            local step_end=$(date +%s)
            local elapsed=$((step_end - step_start))
            log "✓ $step completed in ${elapsed}s"
            STEP_RESULTS+=("$step: OK (${elapsed}s)")
        else
            local step_end=$(date +%s)
            local elapsed=$((step_end - step_start))
            log_error "✗ $step failed after ${elapsed}s"
            STEP_RESULTS+=("$step: FAILED (${elapsed}s)")
            [[ "$STOP_ON_ERROR" == "true" ]] && break
        fi

        echo ""
    done

    local PIPELINE_END=$(date +%s)
    local TOTAL_ELAPSED=$((PIPELINE_END - PIPELINE_START))

    # Write pipeline summary
    local summary="$RUN_DIR/pipeline_summary.txt"
    {
        echo "=============================================="
        echo "GROMACS PPI PIPELINE SUMMARY"
        echo "=============================================="
        echo ""
        echo "Run directory: $RUN_DIR"
        echo "Started:       $(date -d @$PIPELINE_START '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date)"
        echo "Finished:      $(date)"
        echo "Total time:    $((TOTAL_ELAPSED / 3600))h $((TOTAL_ELAPSED % 3600 / 60))m $((TOTAL_ELAPSED % 60))s"
        echo ""
        echo "Step Results:"
        for result in "${STEP_RESULTS[@]}"; do
            echo "  - $result"
        done
        echo ""
        echo "Output structure:"
        ls -1d "$RUN_DIR"/*/ 2>/dev/null | while read -r d; do
            echo "  $(basename "$d")/"
        done
    } > "$summary"

    log_section "Pipeline Complete"
    log "Total time: $((TOTAL_ELAPSED / 3600))h $((TOTAL_ELAPSED % 3600 / 60))m $((TOTAL_ELAPSED % 60))s"
    log "Results: $RUN_DIR"
    echo ""
    log "Step Results:"
    for result in "${STEP_RESULTS[@]}"; do
        log "  - $result"
    done
    echo ""
    log "Summary: $summary"
}

run_pipeline

wait
