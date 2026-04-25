#!/usr/bin/env python3
"""
PyMOL Configuration for SmelGIF-SWI2 Visualization
Target: fold_test_recruitment_of_smelgif_to_swi (SmelGIF-SWI2 complex)

This configuration file contains all settings for:
  - Color schemes for protein types
  - Render settings (resolution, DPI, ray tracing)
  - Combined cartoon and bubble visualization
  - Background colors and output formats
  - Style mode toggles (cartoon/surface)
  - Residue highlighting for interface analysis
"""

# ============================================================
# COLOR SCHEME
# ============================================================
COLORS = {
    'GRF': '#00BCD4',   # Cyan
    'GIF': '#673AB7',   # Purple
    'SWI2': '#4CAF50'   # Green
}

COLORS_HIGHLIGHTED_STICK = {
    'GRF': '#008BA3',
    'GIF': '#4A148C',
    'SWI2': '#2E7D32'
}

COLORS_HIGHLIGHTED_BUBBLE = {
    'GRF': '#00A5BB',
    'GIF': '#5E35B1',
    'SWI2': '#388E3C'
}

# ============================================================
# RENDER CONFIGURATION
# ============================================================
RENDER_CONFIG = {
    'width': 1920,
    'height': 1080,
    'dpi': 300,
    'ray': 1
}

# ============================================================
# BACKGROUND COLORS AND OUTPUT FORMATS
# ============================================================
OUTPUT_FORMATS = {
    'black': {'ext': 'jpeg', 'bg': 'black'},
    'white': {'ext': 'jpeg', 'bg': 'white'},
    'transparent': {'ext': 'png', 'bg': 'white'}
}

# ============================================================
# VISUALIZATION STYLE (CONFIG-TOGGLEABLE)
# ============================================================
# Options: 'cartoon' | 'surface'
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
# CARTOON REPRESENTATION SETTINGS (ADJUSTABLE THICKNESS)
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

# Interface residues highlighting (by residue name, all occurrences)
HYDROPHOBIC_RESIDUES_BY_PROTEIN = {
    'GIF': ['LEU', 'ILE']
}

STICK_RADIUS = 0.2

# ============================================================
# BUBBLE SURFACE SETTINGS (COMBINED WITH CARTOON)
# ============================================================
# Creates smooth, rounded "blobs" around residues
# Cartoon is more visible, bubble is less dark than stick
# 
# ADJUST THIS VALUE TO CHANGE BUBBLE SIZE ⬇️
# 4 versions with increasing bubble sizes:
#   Version 1: 2.0 (smallest, tighter bubbles)
#   Version 2: 3.0 (small-medium bubbles)
#   Version 3: 4.0 (medium-large bubbles)
#   Version 4: 5.0 (largest bubbles)
BUBBLE_SIZE_MULTIPLIER = 3.0  # <-- EASILY ADJUSTABLE BUBBLE SIZE (default: version 2)

# Bubble version for output filename (automatically set based on BUBBLE_SIZE_MULTIPLIER)
BUBBLE_VERSION = 2  # 1=small, 2=medium, 3=large, 4=xlarge

CREATE_BUBBLE_SURFACE = True
BUBBLE_SURFACE_QUALITY = 2  # 0-4, higher = smoother
BUBBLE_SURFACE_RADIUS = BUBBLE_SIZE_MULTIPLIER  # Uses the multiplier value
BUBBLE_SURFACE_TRANSPARENCY = 0.5  # Less dark than stick

# Ray-tracing to soften edges and create smooth blobs
BUBBLE_RAY_TRACE = True
BUBBLE_RAY_WIDTH = 2000  # ray 2000,1500 for softening edges
BUBBLE_RAY_HEIGHT = 1500
BUBBLE_BLUR = True  # Slightly blur for smooth appearance

# ============================================================
# TARGET FOLD CONFIGURATION
# ============================================================
TARGET_FOLD = 'fold_test_recruitment_of_smelgif_to_swi'

# ============================================================
# OUTPUT CONTROL
# ============================================================
# Default: Session file only (images disabled)
RENDER_IMAGES = False

# ============================================================
# PROTEIN SEQUENCES
# ============================================================
SEQUENCES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'MSTTAETVAEGYRTPFTAVQWQELEHQAMIYKYLVAGVPVPADLVVPIRRSFEPISARFCHHPSLGYYSYYGKKFDPEPGRCRRTDGKKWRCAKDAYPDSKYCERHMHRGRNRSRKHVESQSTAPALLTSVSHNTTGSSKTSGNFQRSSSGSFQNTPLYSAANSEGPSYGSATTKMQTEPTTYAIDFKGYFHGMNSDEQNFSFEASAGTRSLGMGSNTDSMWCLMPLQLPSNPMVKPKKDSQLPDSSQPIRMPNPFEPMNDATISGQQHQHCFFSSDIGSPGTVKQEQRSMRPFFDEWPTTKESWSNLDDDGSNKNNFCTPQLSISIPMTPPDFSSRSSCSPNGVSGAALSRQILISTSRWNEPWPRMSKLPLVPPALLTSVSHNTTGSSKTSGNFQRSSSGSFQNTPLYSAANSEGPSYGSATTKMQTEPTTYAIDFKGYFHGMNSDEQNFSFEASAGTRSLGMGSNTDSMWCLMPLQLPSNPMVKPKKDSQLPDSSQPIRMPNPFEPMNDATISGQQHQHCFFSSDIGSPGTVKQEQRSMRPFFDEWPTTKESWSNLDDDGSNKNNFCTPQLSISIPMTPPDFSSRSSCSPNGELTSSSSSPFIQLKSHSIATQQ',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'MQQHLMQMQPMMAAYYPTNVTTDHIQQYLDENKSLILKIVESQNSGKLSECAENQARLQRNLMYLAAIADSQPQPSSMHSQFSSGGMMQPGTHNYLQQQQQQVQQMATQQLMAARSSSMLYGQQQQQPQLSPFQQGLHGSQLGMSSGSGGSTGLHHMLQSESSPHGGGFSHDFVRANKQDIGSSMSAEGRGGNSGGDGGENLYLKASED',
    'SmelSWI2_SMEL4.1_07g019160.1.01': 'MSLNALKETLKPCTNQSFSQSSSTSYNFDTKSVNPRKPPKSSLSQQLLRLEDHTSVLQNQPQTLKKQNHFDLKRKYEKSEEEEEEEEEEEEKGIGFGRPKLDSLLLDQAGPYEPLVLSSPNEKPLVQVPASINCRLMEHQREGVKFLYSLYQNNHGGVLGDDMGLGKTIQSIAFLAAVYCKYGDLPESSVSKERRRTMGPVLIICPSSLIHNWENEFSKWATFSICIYHGSNRDLMIDKLEARGVEIFITSFDTYRIHGRILSDIQWEIVIIDEAHRLKNEKSKLYEACLAIKTQKRYGLTGTIMQNRLMELFNLFDWVIPGCLGTREHFREFYEEPLKHGQRSSAPDRFVRVADERKQHLVSVLRKYLLRRTKEETIGHLMLGKEDNVVFCAMSELQKRVYQRMLLLPDVQCLINKDVPCSCGSPLKQVECCGRTAPDGVIWPYLHRDNPDGCDHCPFCLVLPCLVKLQQISNHLELIKPNPRDDPDKQRRDAEFAAVVFGKDVDLVGGNTQNKSFLGLSNVEHCGKMRALEKLMSSWVLQSDKILLFSYSVRMLDILEKFIIRKGYGFSRLDGSTPTGLRQSLVDDFNSSPSKQVFLLSTKAGGLGLNLVSANRVVIFDPNWNPAHDLQAQDRSFRFGQKRHVIVFRLLAAGSLEELVYTRQVYKQQLSNIAVSGNMEKRYFEGVQVENHFPFFVFKDLTSIFVSFFQLEMITYSCCNFQDSKEFQGELFGICNLFRDLSDKLFTSEIIELHKKNGKEDDGTHSKQDLNVLGMNFVPEKEITTESFVGAESSKHKEEECKAVAPVLEDLGIVYAHRYEDIVNLGLAKIKEKKEQTMHLDYPPRQPKFSTIGKRKSNTITGKESVGTVNPITIRKKSQYGLLARSMGMDVVQFSKWLLSASPAEHEKVLKDYCKRKEKIPNG'
}

# ============================================================
# PROTEIN TYPE MAPPING
# ============================================================
PROTEIN_TYPES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'GRF',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'GIF',
    'SmelSWI2_SMEL4.1_07g019160.1.01': 'SWI2'
}

SHORT_NAMES = {
    'SmelGRF_SMEL4.1_05g020970.1.01': 'SmelGRF',
    'SmelGIF_SMEL4.1_11g026070.1.01': 'SmelGIF',
    'SmelSWI2_SMEL4.1_07g019160.1.01': 'SmelSWI2'
}
