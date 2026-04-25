#!/bin/bash
set -euo pipefail

# PyMOL Automation Script for SmelGRF-GIF and SmelGIF-SWI2 Visualization
# Processes AlphaFold3 models and generates visualizations

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_DIR="${SCRIPT_DIR}"
MODULES_DIR="${SCRIPT_DIR}/modules"
REFERENCE_FASTA="${SCRIPT_DIR}/0_Reference_FASTA_file.fasta"

# Model input files (toggle by commenting/uncommenting)
MODEL_INPUTS=(
    #"fold_1_x_1/fold_1_x_1_model_0.cif"
    #"fold_test_recruitment_of_smelgif_to_swi/#fold_test_recruitment_of_smelgif_to_swi_model_0.cif"
    #"fold_test_recruit_of_smelgrf05_smelgif_to_swi/#fold_test_recruit_of_smelgrf05_smelgif_to_swi_model_0.cif"
    "fold_set_2_smelgrf08_140_smelgif11_070_smelswi2_complex_three_protein_interaction_v2/fold_set_2_smelgrf08_140_smelgif11_070_smelswi2_complex_three_protein_interaction_v2_model_0.cif"
)

# Configuration toggles
CLEAR_LOGS=false  # Toggle: Clear logs or not
CREATE_IMAGES=true  # Default: session file only (set to true for images)

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --clear-logs)
            CLEAR_LOGS=true
            shift
            ;;
        --create-images)
            CREATE_IMAGES=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--clear-logs] [--create-images]"
            exit 1
            ;;
    esac
done

echo "=== PyMOL Visualization Pipeline ==="
echo "Workspace: ${WORKSPACE_DIR}"
echo "Clear logs: ${CLEAR_LOGS}"
echo "Create images: ${CREATE_IMAGES}"
echo ""

# Create logs_protein directory
LOGS_DIR="${WORKSPACE_DIR}/logs_protein"
mkdir -p "${LOGS_DIR}"

# Clear logs if requested
if [ "${CLEAR_LOGS}" = true ]; then
    echo "Clearing previous logs..."
    rm -f "${LOGS_DIR}"/*.log
fi

# Check for reference FASTA
if [ ! -f "${REFERENCE_FASTA}" ]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FASTA}"
    exit 1
fi

# Process model inputs specified in MODEL_INPUTS array
for model_input in "${MODEL_INPUTS[@]}"; do
    cif_file="${WORKSPACE_DIR}/${model_input}"
    
    if [ ! -f "${cif_file}" ]; then
        echo "WARNING: CIF file not found: ${cif_file}"
        echo "  Skipping..."
        continue
    fi
    
    folder_dir="$(dirname "${cif_file}")"
    folder_name="$(basename "${folder_dir}")"
    
    echo "Processing folder: ${folder_name}"
    echo "  CIF file: ${cif_file}"
    
    # Create output directory inside input folder
    output_dir="${folder_dir}/b_run_PyMOL_for_SmelGRF-GIF_and_SmelGIF-SWI2_output"
    mkdir -p "${output_dir}"
    
    # Generate Python script for this model
    python_script="${output_dir}/pymol_script.py"
    
    cat > "${python_script}" << 'PYTHON_EOF'
import os
import sys

# Add modules directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
# script_dir is: .../fold_X/output_dir/
# We need to go up 2 levels to workspace root, then into modules
workspace_root = os.path.dirname(os.path.dirname(script_dir))
modules_dir = os.path.join(workspace_root, 'modules')
sys.path.insert(0, modules_dir)

from pymol import cmd
import pymol_config_protein as config
import pymol_utils_protein as utils

# Configuration
cif_file = sys.argv[1]
output_dir = sys.argv[2]
reference_fasta = sys.argv[3]
create_images = sys.argv[4] == 'true'

model_name = os.path.splitext(os.path.basename(cif_file))[0]

# Initialize PyMOL
cmd.reinitialize()
utils.apply_visual_settings(config.VISUAL_SETTINGS)

# Load model
utils.load_cif_model(cif_file, model_name)

# Parse reference sequences
fasta_sequences = utils.parse_fasta_sequences(reference_fasta)

# Identify chains
chain_mapping = utils.identify_chains(cif_file, fasta_sequences)

# Process each bubble size variant
for bubble_idx, bubble_size in enumerate(config.BUBBLE_SIZE_VARIANTS, 1):
    print(f"\n=== Processing bubble size variant {bubble_idx}: {bubble_size} ===")
    
    # Clear previous styles
    cmd.reinitialize()
    utils.apply_visual_settings(config.VISUAL_SETTINGS)
    utils.load_cif_model(cif_file, model_name)
    
    # Apply styles to each chain
    for chain_id, protein_id in chain_mapping.items():
        color_hex = config.PROTEIN_COLORS.get(protein_id, '#CCCCCC')
        utils.apply_cartoon_bubble_style(model_name, chain_id, color_hex, bubble_size)
    
    # Orient structure
    cmd.orient()
    cmd.zoom('all', buffer=2)
    
    # Apply cartoon settings
    cmd.set('cartoon_sampling', config.RENDER_SETTINGS['cartoon_sampling'])
    cmd.set('cartoon_fancy_helices', config.RENDER_SETTINGS['cartoon_fancy_helices'])
    
    # Save session with dynamic name
    session_name = f"{model_name}_bubble_{bubble_idx}.pse"
    session_path = os.path.join(output_dir, session_name)
    utils.save_session(session_path)
    
    # Render images if requested
    if create_images:
        for bg_name, bg_config in config.BACKGROUND_MODES.items():
            image_name = f"{model_name}_bubble_{bubble_idx}_{bg_name}.{bg_config['format']}"
            image_path = os.path.join(output_dir, image_name)
            utils.render_image(
                image_path,
                bg_config,
                config.RENDER_SETTINGS['width'],
                config.RENDER_SETTINGS['height'],
                config.RENDER_SETTINGS['dpi'],
                config.RENDER_SETTINGS['ray_trace_width'],
                config.RENDER_SETTINGS['ray_trace_height']
            )
    
    print(f"Completed bubble size variant {bubble_idx}")

print(f"\n=== All variants processed ===")
cmd.quit()
PYTHON_EOF
    
    # Execute PyMOL script using conda environment
    echo "  Running PyMOL visualization..."
    log_file="${LOGS_DIR}/${folder_name}_pymol_execution.log"
    
    # Use conda run to execute with pymol environment
    conda run -n pymol --no-capture-output python "${python_script}" \
        "${cif_file}" \
        "${output_dir}" \
        "${REFERENCE_FASTA}" \
        "${CREATE_IMAGES}" \
        > "${log_file}" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Successfully processed ${folder_name}"
    else
        echo "  ✗ Error processing ${folder_name} (check log: ${log_file})"
        cat "${log_file}"
    fi
    
    echo ""
done

echo "=== Pipeline complete ==="
