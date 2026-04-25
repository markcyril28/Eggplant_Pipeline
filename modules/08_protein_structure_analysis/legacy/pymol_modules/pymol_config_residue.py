#!/usr/bin/env python3
"""
PyMOL Configuration Module for Residue Highlighting
Central configuration for protein visualization with hydrophobic residue highlighting

Proofread: Configuration matches instruction specifications:
- Color schemes (cyan, purple, green)
- Residue highlighting (GLN,LEU,VAL,PRO,TYR,ILE for GRF; LEU,ILE for GIF)
- Render settings (1920x1080, 300 DPI, ray tracing 2000x1500)
- Style toggles (cartoon with adjustable thickness)
- Label modes (LABEL_RESN_ONLY toggle)
"""

# Target genes and their color schemes
PROTEIN_COLORS = {
    'SmelGRF_SMEL4.1_05g020970.1.01': '#00BCD4',  # cyan
    'SmelGIF_SMEL4.1_11g026070.1.01': '#673AB7',  # purple
    'SmelSWI2_SMEL4.1_07g019160.1.01': '#4CAF50'  # green
}

# Short names for display
PROTEIN_SHORT_NAMES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'SmelGRF',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'SmelGIF',
    'SmelSWI2_SMEL4.1_07g019160.1.01': 'SmelSWI2'
}

# Hydrophobic residues to highlight by protein
HYDROPHOBIC_RESIDUES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': ['GLN', 'LEU', 'VAL', 'PRO', 'TYR', 'ILE'],
    'SmelGIF_SMEL4.1_11g026070.1.01': ['LEU', 'ILE']
}

# Single-letter to three-letter amino acid mapping
AA_CODE_MAP = {
    'Q': 'GLN', 'L': 'LEU', 'V': 'VAL', 'P': 'PRO',
    'Y': 'TYR', 'I': 'ILE', 'A': 'ALA', 'C': 'CYS',
    'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
    'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'N': 'ASN',
    'R': 'ARG', 'S': 'SER', 'T': 'THR', 'W': 'TRP'
}

# Three-letter to single-letter amino acid mapping (reverse)
AA_CODE_MAP_REVERSE = {v: k for k, v in AA_CODE_MAP.items()}

# Label configuration
LABEL_RESN_ONLY = True  # TRUE: show only residue name, FALSE: show residue name + number

# Background modes
BACKGROUND_MODES = {
    'black': {'bg_color': 'black', 'format': 'jpeg'},
    'white': {'bg_color': 'white', 'format': 'jpeg'},
    'transparent': {'bg_color': 'white', 'format': 'png'}
}

# Render settings
RENDER_SETTINGS = {
    'width': 1920,
    'height': 1080,
    'dpi': 300,
    'ray_trace_mode': 1,
    'ray_width': 2000,
    'ray_height': 1500,
    'cartoon_sampling': 12,
    'cartoon_fancy_helices': 1
}

# Style settings
STYLE_SETTINGS = {
    'stick_radius': 0.2,
    'stick_color_factor': 0.6,  # Darker than base color (stick more visible)
    'bubble_transparency': 0.35,  # Less dark than stick
    'bubble_color_factor': 0.75,  # Bubble less dark than stick
    'surface_quality': 2,
    'surface_solvent': 0
}

# PyMOL visual settings
VISUAL_SETTINGS = {
    'antialias': 2,
    'ambient': 0.4,
    'specular': 0.5,
    'shininess': 10,
    'depth_cue': 0,
    'ray_shadows': 0
}
