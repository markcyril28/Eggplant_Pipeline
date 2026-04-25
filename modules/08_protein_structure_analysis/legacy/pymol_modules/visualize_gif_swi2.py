#!/usr/bin/env python3
"""
PyMOL Visualization for SmelGIF-SWI2 Interaction
Target: fold_test_recruitment_of_smelgif_to_swi (SmelGIF-SWI2 complex)

This script visualizes the SmelGIF-SWI2 protein interaction with:
  - Chain identification via sequence alignment
  - Interface residue highlighting (LEU/ILE in SmelGIF)
  - Custom bubble surfaces around key residues
  - Multiple output formats (Black/White/Transparent)

Usage:
  pymol -c visualize_gif_swi2.py
  pymol -c visualize_gif_swi2.py -- --skip
"""

print("\n" + "="*60, flush=True)
print("Loading SmelGIF-SWI2 Visualization Pipeline...", flush=True)
print("="*60 + "\n", flush=True)

from pymol import cmd
import os
import sys
import logging
from datetime import datetime

# Add script directory to Python path
# Try multiple methods to find the modules directory
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # If __file__ is not defined, use current working directory
    script_dir = os.path.join(os.getcwd(), 'modules')

# Add both the script directory and current working directory to path
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Also add the modules subdirectory if we're in the parent directory
modules_dir = os.path.join(os.getcwd(), 'modules')
if modules_dir not in sys.path and os.path.exists(modules_dir):
    sys.path.insert(0, modules_dir)

print(f"Python path includes: {script_dir}", flush=True)
print(f"Current directory: {os.getcwd()}", flush=True)

# Import configuration and utilities from same directory
from pymol_config_gif_swi2 import (
    COLORS, COLORS_HIGHLIGHTED_STICK, COLORS_HIGHLIGHTED_BUBBLE,
    RENDER_CONFIG, SURFACE_CONFIG, CARTOON_CONFIG,
    VISUALIZATION_STYLE, SEQUENCES, PROTEIN_TYPES, OUTPUT_FORMATS,
    HIGHLIGHT_HYDROPHOBIC, HYDROPHOBIC_RESIDUES_BY_PROTEIN,
    STICK_RADIUS,
    CREATE_BUBBLE_SURFACE, BUBBLE_SURFACE_QUALITY, BUBBLE_SURFACE_RADIUS,
    BUBBLE_SURFACE_TRANSPARENCY, BUBBLE_VERSION,
    TARGET_FOLD, RENDER_IMAGES
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
    log_filename = os.path.join(logs_dir, f"pymol_gif_swi2_{timestamp}.log")
    
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


def check_outputs_exist(fold_name, output_dir):
    """Check if all expected output files exist."""
    expected_files = [
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Black_BG.jpeg",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_White_BG.jpeg",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Transparent.png",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_session.pse"
    ]
    
    for filename in expected_files:
        if not os.path.exists(os.path.join(output_dir, filename)):
            return False
    
    return True


def delete_outputs(fold_name, output_dir):
    """Delete existing output files for fresh run."""
    if not os.path.exists(output_dir):
        return
    
    logging.info(f"  → Deleting existing outputs in {output_dir}...")
    
    output_files = [
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Black_BG.jpeg",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_White_BG.jpeg",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Transparent.png",
        f"{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_session.pse"
    ]
    
    deleted_count = 0
    for filename in output_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
                deleted_count += 1
            except Exception as e:
                logging.warning(f"  ⚠ Could not delete {filename}: {e}")
    
    if deleted_count > 0:
        logging.info(f"  ✓ Deleted {deleted_count} previous output file(s)")


def render_output(object_name, chain_identities, bg_type, bg_color, output_path):
    """Render output image with specified background."""
    try:
        cmd.bg_color(bg_color)
        
        fmt_info = OUTPUT_FORMATS[bg_type]
        
        if fmt_info['ext'] == 'png' and bg_type == 'transparent':
            cmd.set('ray_opaque_background', 0)
        
        cmd.png(output_path,
                width=RENDER_CONFIG['width'],
                height=RENDER_CONFIG['height'],
                dpi=RENDER_CONFIG['dpi'],
                ray=RENDER_CONFIG['ray'])
        
        cmd.set('ray_opaque_background', 1)
        
        logging.info(f"  ✓ Rendered: {os.path.basename(output_path)}")
        
    except Exception as e:
        logging.error(f"  ✗ Error rendering {os.path.basename(output_path)}: {e}")
        raise


def visualize_fold(fold_name, skip_existing=False):
    """Visualize the SmelGIF-SWI2 interaction fold."""
    logging.info("="*60)
    logging.info(f"Processing: {fold_name} (SmelGIF-SWI2)")
    logging.info("="*60)
    
    fold_dir = fold_name
    cif_file = f"{fold_dir}/{fold_name}_model_0.cif"
    output_dir = f"{fold_dir}/output"
    
    os.makedirs(output_dir, exist_ok=True)
    
    if skip_existing and check_outputs_exist(fold_name, output_dir):
        logging.info("✓ Outputs already exist. Skipping...")
        return True
    
    if not skip_existing:
        delete_outputs(fold_name, output_dir)
    
    if not os.path.exists(cif_file):
        logging.error(f"CIF file not found: {cif_file}")
        return False
    
    # Clear PyMOL session
    cmd.reinitialize()
    logging.info("✓ PyMOL session reinitialized")
    
    # Load structure
    logging.info(f"Loading CIF file: {cif_file}")
    try:
        cmd.load(cif_file, fold_name)
        logging.info("✓ Successfully loaded structure")
    except Exception as e:
        logging.error(f"Error loading structure: {e}")
        return False
    
    # Identify chains
    logging.info("Identifying chains...")
    chain_map = identify_chains_from_structure(fold_name, SEQUENCES)
    
    if not chain_map:
        logging.error("Could not identify chains")
        return False
    
    # Build chain identities
    chain_identities = {}
    for chain, protein_name in chain_map.items():
        chain_identities[chain] = {
            'name': protein_name,
            'type': PROTEIN_TYPES[protein_name],
            'color': COLORS[PROTEIN_TYPES[protein_name]]
        }
        logging.info(f"  Chain {chain}: {protein_name} ({PROTEIN_TYPES[protein_name]})")
    
    # Apply colors
    logging.info("Applying colors...")
    for chain, info in chain_identities.items():
        color_name = f"{info['type']}_color_{chain}_{fold_name}"
        rgb_values = hex_to_rgb(info['color'])
        cmd.set_color(color_name, rgb_values)
        cmd.color(color_name, f"{fold_name} and chain {chain}")
        logging.info(f"  Chain {chain}: {info['name']} → {color_name} {info['color']}")
    logging.info("✓ Colors applied")
    
    # Set visualization style
    logging.info("Configuring visualization...")
    cmd.hide('everything', fold_name)
    
    if VISUALIZATION_STYLE == 'cartoon':
        cmd.show('cartoon', fold_name)
        for setting_name, setting_value in CARTOON_CONFIG.items():
            cmd.set(setting_name, setting_value)
        logging.info("✓ Cartoon representation configured")
    elif VISUALIZATION_STYLE == 'surface':
        cmd.show('surface', fold_name)
        cmd.set('surface_quality', SURFACE_CONFIG['surface_quality'])
        cmd.set('solvent_radius', SURFACE_CONFIG['solvent_radius'])
        for chain, info in chain_identities.items():
            color_name = f"{info['type']}_color_{chain}_{fold_name}"
            cmd.set('surface_color', color_name, f"{fold_name} and chain {chain}")
        logging.info("✓ Surface representation configured")
    
    cmd.orient(fold_name)
    cmd.zoom(fold_name, complete=1)
    
    # Set rendering quality (shadows off for clarity)
    cmd.set('ray_shadows', 0)
    cmd.set('antialias', 2)
    logging.info("✓ Rendering settings applied (shadows off)")
    
    # Highlight interface residues
    if HIGHLIGHT_HYDROPHOBIC:
        logging.info("Highlighting interface residues...")
        highlight_specific_residues(
            fold_name, HYDROPHOBIC_RESIDUES_BY_PROTEIN, chain_identities,
            STICK_RADIUS, False,  # No labels
            CREATE_BUBBLE_SURFACE, BUBBLE_SURFACE_QUALITY, BUBBLE_SURFACE_RADIUS,
            BUBBLE_SURFACE_TRANSPARENCY, COLORS_HIGHLIGHTED_STICK, COLORS_HIGHLIGHTED_BUBBLE
        )
    
    # Generate outputs
    logging.info("Generating outputs...")
    
    try:
        if RENDER_IMAGES:
            outputs = {
                'black': f"{output_dir}/{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Black_BG.jpeg",
                'white': f"{output_dir}/{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_White_BG.jpeg",
                'transparent': f"{output_dir}/{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_render_Transparent.png"
            }
            
            for bg_type, output_path in outputs.items():
                fmt_info = OUTPUT_FORMATS[bg_type]
                render_output(fold_name, chain_identities, bg_type, fmt_info['bg'], output_path)
        else:
            logging.info("Image rendering disabled (RENDER_IMAGES=False)")
        
        # Save session
        session_file = f"{output_dir}/{fold_name}_model_0_bubble_v{BUBBLE_VERSION}_session.pse"
        cmd.save(session_file)
        logging.info(f"  ✓ Session saved: {os.path.basename(session_file)}")
    except Exception as e:
        logging.error(f"Error generating outputs: {e}")
        return False
    
    # Summary
    logging.info("="*60)
    logging.info(f"Visualization Summary for {fold_name}")
    logging.info("="*60)
    for chain, info in chain_identities.items():
        logging.info(f"Chain {chain}: {info['name']}")
        logging.info(f"  Type: {info['type']} | Color: {info['color']}")
    logging.info(f"Outputs saved to: {output_dir}/")
    if RENDER_IMAGES:
        logging.info("  • Black Background (JPEG)")
        logging.info("  • White Background (JPEG)")
        logging.info("  • Transparent (PNG)")
    logging.info("  • PyMOL Session (.pse)")
    logging.info("="*60)
    
    return True


def main():
    """Main execution function."""
    print("\n" + "="*60, flush=True)
    print("SmelGIF-SWI2 Visualization Pipeline", flush=True)
    print("="*60 + "\n", flush=True)
    
    log_file = setup_logging()
    print(f"✓ Logging initialized: {log_file}\n", flush=True)
    
    skip_existing = False
    if len(sys.argv) > 1:
        if sys.argv[1].lower() in ['--skip', '-s', 'skip']:
            skip_existing = True
    
    mode_text = 'SKIP EXISTING' if skip_existing else 'FRESH RUN'
    print(f"Mode: {mode_text}\n", flush=True)
    
    logging.info("="*60)
    logging.info("SmelGIF-SWI2 Visualization Pipeline")
    logging.info("="*60)
    logging.info(f"Mode: {mode_text}")
    logging.info(f"Log file: {log_file}")
    logging.info(f"Target fold: {TARGET_FOLD}")
    logging.info("="*60)
    
    # Process target fold
    result = visualize_fold(TARGET_FOLD, skip_existing)
    
    if result:
        print("\n" + "="*60, flush=True)
        print("[OK] SmelGIF-SWI2 visualization complete!", flush=True)
        print(f"Log saved to: {log_file}", flush=True)
        print("="*60 + "\n", flush=True)
        
        logging.info("="*60)
        logging.info("[OK] SmelGIF-SWI2 visualization complete!")
        logging.info(f"Log saved to: {log_file}")
        logging.info("="*60)
    else:
        print("\n" + "="*60, flush=True)
        print("[ERROR] Visualization failed", flush=True)
        print("="*60 + "\n", flush=True)
        
        logging.error("="*60)
        logging.error("[ERROR] Visualization failed")
        logging.error("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()
else:
    main()
