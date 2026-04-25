#!/usr/bin/env python3
"""
PyMOL Visualization for SmelGRF-GIF with Residue Labels and Leader Lines
Target: fold_1_x_1 (SmelGRF-GIF complex)

This script creates detailed residue-focused visualizations with:
  - Leader lines connecting labels to residues
  - Detailed residue highlighting (GLN, LEU, VAL, PRO, TYR, ILE in GRF; LEU, ILE in GIF)
  - Both full structure and helix-only PSE sessions
  - Custom bubble surfaces with adjustable sizing
  - Automatic view centering on all structures

Usage:
  pymol -c visualize_gif_grf_residue.py
"""

print("\n" + "="*60, flush=True)
print("Loading SmelGRF-GIF Residue Visualization Pipeline...", flush=True)
print("="*60 + "\n", flush=True)

from pymol import cmd
import os
import sys
import logging
from datetime import datetime

# Add script directory to Python path
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.path.join(os.getcwd(), 'modules')

if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

modules_dir = os.path.join(os.getcwd(), 'modules')
if modules_dir not in sys.path and os.path.exists(modules_dir):
    sys.path.insert(0, modules_dir)

print(f"Python path includes: {script_dir}", flush=True)
print(f"Current directory: {os.getcwd()}", flush=True)

from pymol_config_gif_grf_residue import (
    COLORS, COLORS_HIGHLIGHTED_STICK, COLORS_HIGHLIGHTED_BUBBLE,
    RENDER_CONFIG, SURFACE_CONFIG, CARTOON_CONFIG,
    VISUALIZATION_STYLE, SEQUENCES, PROTEIN_TYPES, SHOW_LABELS, OUTPUT_FORMATS,
    HIGHLIGHT_HYDROPHOBIC, HYDROPHOBIC_RESIDUES_BY_PROTEIN,
    STICK_RADIUS, LABEL_RESN_ONLY, USE_LEADER_LINES,
    CREATE_BUBBLE_SURFACE, BUBBLE_SURFACE_QUALITY, BUBBLE_SURFACE_RADIUS,
    BUBBLE_SURFACE_TRANSPARENCY, BUBBLE_SURFACE_TYPE, BUBBLE_SURFACE_SOLVENT,
    TARGET_FOLD, CREATE_HELIX_ONLY_PSE, HELIX_ONLY_PSE_SUFFIX, RENDER_IMAGES
)
from pymol_utils import (
    hex_to_rgb, identify_chains_from_structure, label_chains,
    highlight_specific_residues
)

print("✓ All modules imported successfully\n", flush=True)


def setup_logging():
    """Setup logging to both console and file."""
    logs_dir = "logs"
    os.makedirs(logs_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(logs_dir, f"pymol_gif_grf_residue_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stderr)
        ],
        force=True
    )
    
    return log_filename


def apply_leader_lines():
    """Enable leader lines for labels."""
    cmd.set('label_connector', 1)
    cmd.set('label_connector_mode', 1)  # Mode 1: leader lines
    cmd.set('label_connector_width', 2.0)
    cmd.set('label_connector_color', 'white')
    cmd.set('label_size', 14)
    cmd.set('label_position', (2, 2, 2))
    cmd.set('label_outline_color', 'black')
    logging.info("  ✓ Leader lines enabled with outline")


def label_hydrophobic_residues(chain_mapping):
    """Apply labels to hydrophobic residues with leader lines."""
    logging.info("Applying residue labels with leader lines...")
    
    for chain_id, protein_name in chain_mapping.items():
        protein_type = PROTEIN_TYPES.get(protein_name)
        if not protein_type or protein_type not in HYDROPHOBIC_RESIDUES_BY_PROTEIN:
            continue
        
        residues = HYDROPHOBIC_RESIDUES_BY_PROTEIN[protein_type]
        resn_list = '+'.join(residues)
        selection = f"(chain {chain_id} and resn {resn_list} and name CA)"
        
        if LABEL_RESN_ONLY:
            cmd.label(selection, '"%s" % resn')
        else:
            cmd.label(selection, '"%s%s" % (resn, resi)')
        
        logging.info(f"  ✓ Labeled {protein_type} ({chain_id}): {', '.join(residues)}")


def visualize_fold():
    """Main visualization function for fold_1_x_1."""
    fold_dir = TARGET_FOLD
    model_file = os.path.join(fold_dir, f"{TARGET_FOLD}_model_0.cif")
    
    if not os.path.exists(model_file):
        logging.error(f"Model file not found: {model_file}")
        return False
    
    # Output directory inside the fold directory as per instructions
    output_dir = os.path.join(fold_dir, "OUTPUT_RESIDUE")
    os.makedirs(output_dir, exist_ok=True)
    
    logging.info(f"Processing: {TARGET_FOLD}")
    logging.info(f"Model: {model_file}")
    logging.info(f"Output: {output_dir}")
    
    # Load model
    cmd.load(model_file, TARGET_FOLD)
    logging.info("✓ Model loaded")
    
    # Identify chains
    chain_mapping = identify_chains_from_structure(TARGET_FOLD, SEQUENCES)
    if not chain_mapping:
        logging.error("Failed to identify chains")
        return False
    
    logging.info(f"✓ Chains identified: {chain_mapping}")
    
    # Apply cartoon settings
    for key, value in CARTOON_CONFIG.items():
        cmd.set(key, value)
    cmd.set('ray_shadows', 0)
    
    # Color chains
    for chain_id, protein_name in chain_mapping.items():
        protein_type = PROTEIN_TYPES[protein_name]
        color_hex = COLORS[protein_type]
        color_rgb = hex_to_rgb(color_hex)
        color_name = f"{protein_type}_base"
        cmd.set_color(color_name, color_rgb)
        cmd.color(color_name, f"chain {chain_id}")
        logging.info(f"✓ Colored {protein_type} (chain {chain_id}): {color_hex}")
    
    # Highlight hydrophobic residues
    if HIGHLIGHT_HYDROPHOBIC:
        for chain_id, protein_name in chain_mapping.items():
            protein_type = PROTEIN_TYPES[protein_name]
            if protein_type in HYDROPHOBIC_RESIDUES_BY_PROTEIN:
                residues = HYDROPHOBIC_RESIDUES_BY_PROTEIN[protein_type]
                
                # Stick color (darkest)
                stick_hex = COLORS_HIGHLIGHTED_STICK[protein_type]
                stick_rgb = hex_to_rgb(stick_hex)
                stick_color = f"{protein_type}_stick"
                cmd.set_color(stick_color, stick_rgb)
                
                # Bubble color (medium-dark)
                bubble_hex = COLORS_HIGHLIGHTED_BUBBLE[protein_type]
                bubble_rgb = hex_to_rgb(bubble_hex)
                bubble_color = f"{protein_type}_bubble"
                cmd.set_color(bubble_color, bubble_rgb)
                
                # Select and show sticks
                resn_list = '+'.join(residues)
                selection_name = f"{protein_type}_hydrophobic"
                cmd.select(selection_name, f"chain {chain_id} and resn {resn_list}")
                cmd.show('sticks', selection_name)
                cmd.set('stick_radius', STICK_RADIUS, selection_name)
                cmd.color(stick_color, selection_name)
                
                # Create bubble surface
                if CREATE_BUBBLE_SURFACE:
                    bubble_name = f"{selection_name}_bubble"
                    cmd.create(bubble_name, selection_name)
                    cmd.show('surface', bubble_name)
                    cmd.set('surface_quality', BUBBLE_SURFACE_QUALITY, bubble_name)
                    cmd.set('surface_type', BUBBLE_SURFACE_TYPE, bubble_name)  # Solid surface
                    cmd.set('surface_solvent', BUBBLE_SURFACE_SOLVENT, bubble_name)
                    cmd.set('surface_ramp_above_mode', 1, bubble_name)
                    cmd.set('solvent_radius', BUBBLE_SURFACE_RADIUS, bubble_name)
                    cmd.set('transparency', BUBBLE_SURFACE_TRANSPARENCY, bubble_name)
                    cmd.color(bubble_color, bubble_name)
                    logging.info(f"✓ Created smooth bubble surface for {protein_type}")
                
                logging.info(f"✓ Highlighted {protein_type}: {', '.join(residues)}")
    
    # Apply labels with leader lines
    if SHOW_LABELS:
        label_hydrophobic_residues(chain_mapping)
        if USE_LEADER_LINES:
            apply_leader_lines()
    
    # Center view on all atoms
    cmd.zoom('all', buffer=5.0)
    
    # Save full PSE session
    pse_file = os.path.join(output_dir, f"{TARGET_FOLD}_model_0_session_residue.pse")
    cmd.save(pse_file)
    logging.info(f"✓ Saved full session: {pse_file}")
    
    # Create helix-only version
    if CREATE_HELIX_ONLY_PSE:
        logging.info("Creating helix-only version...")
        cmd.hide('everything')
        cmd.show('cartoon', 'ss h')  # Show only helices
        
        # Show sticks and bubbles ONLY for hydrophobic residues in helix regions
        for chain_id, protein_name in chain_mapping.items():
            protein_type = PROTEIN_TYPES[protein_name]
            if protein_type in HYDROPHOBIC_RESIDUES_BY_PROTEIN:
                residues = HYDROPHOBIC_RESIDUES_BY_PROTEIN[protein_type]
                resn_list = '+'.join(residues)
                helix_hydro = f"{protein_type}_helix_hydro"
                cmd.select(helix_hydro, f"chain {chain_id} and resn {resn_list} and ss h")
                
                if cmd.count_atoms(helix_hydro) > 0:
                    cmd.show('sticks', helix_hydro)
                    
                    # Apply labels only to helix residues
                    if SHOW_LABELS:
                        if LABEL_RESN_ONLY:
                            cmd.label(f"{helix_hydro} and name CA", '"%s" % resn')
                        else:
                            cmd.label(f"{helix_hydro} and name CA", '"%s%s" % (resn, resi)')
                    
                    if CREATE_BUBBLE_SURFACE:
                        bubble_name = f"{helix_hydro}_bubble"
                        cmd.create(bubble_name, helix_hydro)
                        cmd.show('surface', bubble_name)
                        cmd.set('surface_quality', BUBBLE_SURFACE_QUALITY, bubble_name)
                        cmd.set('surface_type', BUBBLE_SURFACE_TYPE, bubble_name)
                        cmd.set('surface_solvent', BUBBLE_SURFACE_SOLVENT, bubble_name)
                        cmd.set('transparency', BUBBLE_SURFACE_TRANSPARENCY, bubble_name)
        
        cmd.zoom('all', buffer=5.0)
        pse_helix_file = os.path.join(output_dir, f"{TARGET_FOLD}_model_0_session_residue{HELIX_ONLY_PSE_SUFFIX}.pse")
        cmd.save(pse_helix_file)
        logging.info(f"✓ Saved helix-only session: {pse_helix_file}")
        
        # Restore full view for rendering
        cmd.show('cartoon')
        for chain_id, protein_name in chain_mapping.items():
            protein_type = PROTEIN_TYPES[protein_name]
            if protein_type in HYDROPHOBIC_RESIDUES_BY_PROTEIN:
                selection_name = f"{protein_type}_hydrophobic"
                cmd.show('sticks', selection_name)
    
    # Render images
    if RENDER_IMAGES:
        render_output(output_dir, TARGET_FOLD)
    
    cmd.delete('all')
    return True


def render_output(output_dir, fold_name):
    """Render images in multiple backgrounds with ray tracing."""
    logging.info("Rendering images with ray tracing...")
    
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_shadows', 0)
    
    for bg_name, bg_info in OUTPUT_FORMATS.items():
        cmd.bg_color(bg_info['bg'])
        
        # Set transparency for PNG format
        if bg_info['ext'] == 'png' and bg_name == 'transparent':
            cmd.set('ray_opaque_background', 0)
        else:
            cmd.set('ray_opaque_background', 1)
        
        output_file = os.path.join(
            output_dir,
            f"{fold_name}_model_0_render_residue_{bg_name.capitalize()}_BG.{bg_info['ext']}"
        )
        
        # Apply ray-trace with higher resolution for softer edges
        cmd.ray(RENDER_CONFIG['ray_width'], RENDER_CONFIG['ray_height'])
        cmd.png(output_file, dpi=RENDER_CONFIG['dpi'])
        
        logging.info(f"  ✓ Rendered {bg_name}: {output_file}")


def main():
    """Main execution."""
    log_file = setup_logging()
    logging.info("="*60)
    logging.info("SmelGRF-GIF Residue Visualization with Leader Lines")
    logging.info("="*60)
    logging.info(f"Log file: {log_file}")
    logging.info(f"Target fold: {TARGET_FOLD}")
    logging.info(f"Bubble size multiplier: {BUBBLE_SURFACE_RADIUS}")
    logging.info(f"Leader lines: {'Enabled' if USE_LEADER_LINES else 'Disabled'}")
    logging.info("")
    
    success = visualize_fold()
    
    if success:
        logging.info("")
        logging.info("="*60)
        logging.info("✓ Visualization completed successfully!")
        logging.info("="*60)
    else:
        logging.error("")
        logging.error("="*60)
        logging.error("✗ Visualization failed")
        logging.error("="*60)
        sys.exit(1)


if __name__ == '__main__' or __name__ == 'pymol':
    main()
