#!/usr/bin/env python3
"""
PyMOL Configuration for SmelGRF-GIF Residue-Focused Visualization
Target: fold_1_x_1 (SmelGRF-GIF complex with detailed residue labels)

Configuration for residue highlighting with leader_lines labels
"""

# ============================================================
# COLOR SCHEME
# ============================================================
COLORS = {
    'GRF': '#00BCD4',   # Cyan
    'GIF': '#673AB7'    # Purple
}

COLORS_HIGHLIGHTED_STICK = {
    'GRF': '#007A8C',   # Darker cyan for sticks (more visible)
    'GIF': '#38006B'    # Darker purple for sticks (more visible)
}

COLORS_HIGHLIGHTED_BUBBLE = {
    'GRF': '#009AB3',   # Medium-dark cyan for bubbles (less dark than stick)
    'GIF': '#4A0E7A'    # Medium-dark purple for bubbles (less dark than stick)
}

# ============================================================
# RENDER CONFIGURATION
# ============================================================
RENDER_CONFIG = {
    'width': 1920,
    'height': 1080,
    'dpi': 300,
    'ray': 1,
    'ray_width': 2000,
    'ray_height': 1500
}

# ============================================================
# BACKGROUND COLORS AND OUTPUT FORMATS
# ============================================================
OUTPUT_FORMATS = {
    'black': {'ext': 'png', 'bg': 'black'},
    'white': {'ext': 'png', 'bg': 'white'},
    'transparent': {'ext': 'png', 'bg': 'white'}
}

# ============================================================
# VISUALIZATION STYLE
# ============================================================
VISUALIZATION_STYLE = 'cartoon'

# ============================================================
# SURFACE QUALITY AND SOLVENT ACCESSIBILITY
# ============================================================
SURFACE_CONFIG = {
    'surface_quality': 2,
    'solvent_radius': 3.0,
    'thickness': 1.0
}

# ============================================================
# CARTOON REPRESENTATION SETTINGS
# ============================================================
CARTOON_CONFIG = {
    'cartoon_sampling': 7,
    'cartoon_loop_radius': 0.2,
    'cartoon_rect_length': 1.4,
    'cartoon_rect_width': 0.4,
    'cartoon_oval_length': 1.2,
    'cartoon_oval_width': 0.25,
    'cartoon_tube_radius': 0.4,
    'cartoon_tube_cap': 1,
    'cartoon_cylindrical_helices': 0,
    'cartoon_helix_radius': 1.6,
    'cartoon_ring_finder': 3,
    'cartoon_ring_mode': 3,
    'cartoon_ring_width': 0.22,
    'cartoon_ring_radius': 0.45,
    'cartoon_ladder_mode': 1,
    'cartoon_ladder_radius': 0.25,
    'cartoon_nucleic_acid_mode': 0,
    'cartoon_fancy_helices': 1,
    'cartoon_fancy_sheets': 1,
    'cartoon_smooth_loops': 1,
    'cartoon_flat_sheets': 1,
    'cartoon_side_chain_helper': 0,
    'cartoon_transparency': 0.0,
    'cartoon_discrete_colors': 0
}

# ============================================================
# LABEL SETTINGS WITH LEADER LINES
# ============================================================
SHOW_LABELS = True
LABEL_RESN_ONLY = False  # Show both residue name and number
USE_LEADER_LINES = True  # Enable leader lines for labels

# ============================================================
# AMINO ACID CODE MAPPING
# ============================================================
AMINO_ACID_CODES = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

AMINO_ACID_CODES_REVERSE = {v: k for k, v in AMINO_ACID_CODES.items()}

# ============================================================
# HYDROPHOBIC RESIDUE HIGHLIGHTING
# ============================================================
HIGHLIGHT_HYDROPHOBIC = True

HYDROPHOBIC_RESIDUES_BY_PROTEIN = {
    'GRF': ['GLN', 'LEU', 'VAL', 'PRO', 'TYR', 'ILE'],
    'GIF': ['LEU', 'ILE']
}

STICK_RADIUS = 0.2

# ============================================================
# BUBBLE SURFACE SETTINGS
# ============================================================
BUBBLE_SIZE_MULTIPLIER = 3.0  # <-- EASILY ADJUSTABLE BUBBLE SIZE

CREATE_BUBBLE_SURFACE = True
BUBBLE_SURFACE_QUALITY = 2
BUBBLE_SURFACE_RADIUS = BUBBLE_SIZE_MULTIPLIER
BUBBLE_SURFACE_TRANSPARENCY = 0.35  # Less dark than stick
BUBBLE_RAY_TRACE = True
BUBBLE_RAY_WIDTH = 2000
BUBBLE_RAY_HEIGHT = 1500
BUBBLE_BLUR = True
BUBBLE_SURFACE_TYPE = 0  # Solid surface
BUBBLE_SURFACE_SOLVENT = 0  # No solvent exclusion

# ============================================================
# TARGET FOLD CONFIGURATION
# ============================================================
TARGET_FOLD = 'fold_2_x_1'

# ============================================================
# ALPHA-HELIX ONLY PSE CONFIGURATION
# ============================================================
CREATE_HELIX_ONLY_PSE = True
HELIX_ONLY_PSE_SUFFIX = '_helix_only'

# ============================================================
# OUTPUT CONTROL
# ============================================================
RENDER_IMAGES = True

# ============================================================
# PROTEIN SEQUENCES
# ============================================================
SEQUENCES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'MSTTAETVAEGYRTPFTAVQWQELEHQAMIYKYLVAGVPVPADLVVPIRRSFEPISARFCHHPSLGYYSYYGKKFDPEPGRCRRTDGKKWRCAKDAYPDSKYCERHMHRGRNRSRKHVESQSTAPALLTSVSHNTTGSSKTSGNFQRSSSGSFQNTPLYSAANSEGPSYGSATTKMQTEPTTYAIDFKGYFHGMNSDEQNFSFEASAGTRSLGMGSNTDSMWCLMPLQLPSNPMVKPKKDSQLPDSSQPIRMPNPFEPMNDATISGQQHQHCFFSSDIGSPGTVKQEQRSMRPFFDEWPTTKESWSNLDDDGSNKNNFCTPQLSISIPMTPPDFSSRSSCSPNGVSGAALSRQILISTSRWNEPWPRMSKLPLVPPALLTSVSHNTTGSSKTSGNFQRSSSGSFQNTPLYSAANSEGPSYGSATTKMQTEPTTYAIDFKGYFHGMNSDEQNFSFEASAGTRSLGMGSNTDSMWCLMPLQLPSNPMVKPKKDSQLPDSSQPIRMPNPFEPMNDATISGQQHQHCFFSSDIGSPGTVKQEQRSMRPFFDEWPTTKESWSNLDDDGSNKNNFCTPQLSISIPMTPPDFSSRSSCSPNGELTSSSSSPFIQLKSHSIATQQ',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'MQQHLMQMQPMMAAYYPTNVTTDHIQQYLDENKSLILKIVESQNSGKLSECAENQARLQRNLMYLAAIADSQPQPSSMHSQFSSGGMMQPGTHNYLQQQQQQVQQMATQQLMAARSSSMLYGQQQQQPQLSPFQQGLHGSQLGMSSGSGGSTGLHHMLQSESSPHGGGFSHDFVRANKQDIGSSMSAEGRGGNSGGDGGENLYLKASED'
}

# ============================================================
# PROTEIN TYPE MAPPING
# ============================================================
PROTEIN_TYPES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'GRF',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'GIF'
}

SHORT_NAMES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'SmelGRF',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'SmelGIF'
}
