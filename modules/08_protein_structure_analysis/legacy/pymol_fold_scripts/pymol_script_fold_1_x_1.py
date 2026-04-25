#!/usr/bin/env python3
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
