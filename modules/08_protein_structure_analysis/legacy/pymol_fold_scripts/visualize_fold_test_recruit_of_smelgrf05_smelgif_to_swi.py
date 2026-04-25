#!/usr/bin/env python3
"""
PyMOL visualization script for fold_test_recruit_of_smelgrf05_smelgif_to_swi_model_0.cif
Identifies chains by comparing with AlphaFold3_Inputs.fasta
Colors: GRF (#00BCD4), GIF (#673AB7), SWI2 (#4CAF50)
Outputs: Black BG (JPEG), White BG (JPEG), Transparent (PNG)
"""

from pymol import cmd
import os
import sys

# Add parent directory to path to import configuration and utilities
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pymol_config import (COLORS, RENDER_CONFIG, SURFACE_CONFIG, SEQUENCES, 
                          PROTEIN_TYPES, SHOW_LABELS, OUTPUT_FORMATS)
from pymol_utils import (hex_to_rgb, identify_chains_from_json, 
                        identify_chains_from_structure, label_chains)

# Define fold name and paths
FOLD_NAME = 'fold_test_recruit_of_smelgrf05_smelgif_to_swi'
FOLD_DIR = f'{FOLD_NAME}'
CIF_FILE = f'{FOLD_DIR}/{FOLD_NAME}_model_0.cif'
JSON_FILE = f'{FOLD_DIR}/{FOLD_NAME}_job_request.json'
OUTPUT_DIR = f'{FOLD_DIR}/output'

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Import configuration
sequences = SEQUENCES
protein_types = PROTEIN_TYPES
colors = COLORS

print(f"\n{'='*60}")
print(f"PyMOL Visualization Pipeline for {FOLD_NAME}")
print(f"{'='*60}\n")

# Clear everything and start fresh
cmd.reinitialize()
print("✓ PyMOL session reinitialized")

# Load structure
print(f"Loading CIF file: {CIF_FILE}")
if not os.path.exists(CIF_FILE):
    print(f"ERROR: CIF file not found: {CIF_FILE}")
    sys.exit(1)

cmd.load(CIF_FILE, FOLD_NAME)
print(f"✓ Successfully loaded structure\n")

# Identify chains - try JSON method first (faster), then structure alignment
print("Identifying chains...")
chain_map = identify_chains_from_json(JSON_FILE, sequences)

if chain_map:
    print("✓ Chains identified from JSON file (fast method)")
else:
    print("⚠ JSON method failed, using structure alignment (slower method)")
    chain_map = identify_chains_from_structure(FOLD_NAME, sequences)

if not chain_map:
    print("ERROR: Could not identify chains")
    sys.exit(1)

# Build chain identities dictionary
chain_identities = {}
for chain, protein_name in chain_map.items():
    chain_identities[chain] = {
        'name': protein_name,
        'type': protein_types[protein_name],
        'color': colors[protein_types[protein_name]]
    }
    print(f"  Chain {chain}: {protein_name} ({protein_types[protein_name]})")

print()

# Apply colors to chains
print("Applying colors...")
for chain, info in chain_identities.items():
    color_name = f"{info['type']}_color_{chain}"
    cmd.set_color(color_name, hex_to_rgb(info['color']))
    cmd.color(color_name, f'{FOLD_NAME} and chain {chain}')
print("✓ Colors applied\n")

# Set visualization style
print("Configuring visualization...")
cmd.hide('everything', FOLD_NAME)
cmd.show('surface', FOLD_NAME)
cmd.set('surface_quality', SURFACE_CONFIG['surface_quality'])
cmd.set('solvent_radius', SURFACE_CONFIG['solvent_radius'])
print("✓ Surface representation configured\n")

# Orient the view
cmd.orient(FOLD_NAME)
cmd.zoom(FOLD_NAME, complete=1)

# Function to render outputs
def render_output(bg_type, bg_color, output_path):
    """Render output with specific background"""
    cmd.bg_color(bg_color)
    
    # Add labels if enabled
    if SHOW_LABELS:
        label_chains(FOLD_NAME, chain_identities)
    
    # Get output format info
    fmt_info = OUTPUT_FORMATS[bg_type]
    
    # Render based on format
    if fmt_info['ext'] == 'png':
        # PNG format supports transparency
        if bg_type == 'transparent':
            cmd.set('ray_opaque_background', 0)
        cmd.png(output_path, 
                width=RENDER_CONFIG['width'], 
                height=RENDER_CONFIG['height'], 
                dpi=RENDER_CONFIG['dpi'], 
                ray=RENDER_CONFIG['ray'])
        # Reset transparency
        cmd.set('ray_opaque_background', 1)
    else:
        # JPEG format
        cmd.png(output_path, 
                width=RENDER_CONFIG['width'], 
                height=RENDER_CONFIG['height'], 
                dpi=RENDER_CONFIG['dpi'], 
                ray=RENDER_CONFIG['ray'])
    
    # Remove labels for next render
    if SHOW_LABELS:
        cmd.label('all', '')
    
    print(f"✓ Rendered: {output_path}")

# Generate all outputs
print("Generating outputs...")
outputs = {
    'black': f'{OUTPUT_DIR}/{FOLD_NAME}_model_0_render_Black_BG.jpeg',
    'white': f'{OUTPUT_DIR}/{FOLD_NAME}_model_0_render_White_BG.jpeg',
    'transparent': f'{OUTPUT_DIR}/{FOLD_NAME}_model_0_render_Transparent.png'
}

for bg_type, output_path in outputs.items():
    fmt_info = OUTPUT_FORMATS[bg_type]
    render_output(bg_type, fmt_info['bg'], output_path)

# Save PyMOL session
session_file = f'{OUTPUT_DIR}/{FOLD_NAME}_model_0_session.pse'
cmd.save(session_file)
print(f"✓ Session saved: {session_file}\n")

# Print summary
print(f"{'='*60}")
print(f"Visualization Summary for {FOLD_NAME}")
print(f"{'='*60}")
for chain, info in chain_identities.items():
    print(f"Chain {chain}: {info['name']}")
    print(f"  Type: {info['type']} | Color: {info['color']}")
print(f"\nOutputs saved to: {OUTPUT_DIR}/")
print(f"  • Black Background (JPEG)")
print(f"  • White Background (JPEG)")
print(f"  • Transparent (PNG)")
print(f"  • PyMOL Session (.pse)")
print(f"{'='*60}\n")
