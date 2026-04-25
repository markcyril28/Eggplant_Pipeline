#!/usr/bin/env python3
"""
PyMOL Configuration Module
Central configuration for protein visualization
"""

# Target gene identifiers
TARGET_GENES = {
    'SmelGRF': 'SmelGRF_SMEL4.1_05g020970.1.01',
    'SmelGIF': 'SmelGIF_SMEL4.1_11g026070.1.01',
    'SmelSWI2': 'SmelSWI2_SMEL4.1_07g019160.1.01'
}

# Color scheme for target genes
PROTEIN_COLORS = {
    'SmelGRF_SMEL4.1_05g020970.1.01': '#00BCD4',  # cyan
    'SmelGIF_SMEL4.1_11g026070.1.01': '#673AB7',  # purple
    'SmelSWI2_SMEL4.1_07g019160.1.01': '#4CAF50'  # green
}

# Model input files (toggle by commenting/uncommenting)
MODEL_INPUTS = [
    'fold_1_x_1/fold_1_x_1_model_0.cif',
    # 'fold_test_recruitment_of_smelgif_to_swi/fold_test_recruitment_of_smelgif_to_swi_model_0.cif',
    # 'fold_test_recruit_of_smelgrf05_smelgif_to_swi/fold_test_recruit_of_smelgrf05_smelgif_to_swi_model_0.cif',
]

# BUBBLE SIZE VARIANTS - Adjust this value to control bubble size
# Location: BUBBLE_SIZE_VARIANTS (Easy to locate for adjustments)
BUBBLE_SIZE_VARIANTS = [1.0, 1.5, 2.0, 2.5]  # 4 versions with increasing size

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
    'ray_trace_width': 2000,
    'ray_trace_height': 1500,
    'cartoon_sampling': 12,  # Cartoon smoothing
    'cartoon_fancy_helices': 1
}

# Style modes (config-toggleable)
STYLE_MODES = {
    'cartoon': {
        'enabled': True,
        'cartoon_thickness': 1.5,
        'cartoon_loop_radius': 0.3
    },
    'surface': {
        'enabled': True,
        'surface_quality': 2,
        'transparency': 0.15  # Less dark than stick
    }
}

# FASTA reference file
FASTA_REFERENCE = '0_Reference_FASTA_file.fasta'

# PyMOL visual settings
VISUAL_SETTINGS = {
    'antialias': 2,
    'ambient': 0.4,
    'specular': 0.5,
    'shininess': 10,
    'depth_cue': 0,
    'ray_shadows': 0
}
