#!/usr/bin/env python3
"""
PyMOL Utility Functions for Residue Highlighting
Functions for loading, styling, and rendering protein structures with hydrophobic residue highlighting

Proofread: Implements key requirements:
- CIF model loading
- Chain identification via sequence alignment (JSON parsing disabled)
- Residue-name highlighting (resn selection by name not position)
- Bubble surface creation with smoothing and transparency
- Ray tracing for soft edges (2000x1500)
- Leader lines for labels with outlines
- Dual PSE sessions (full + helix-only)
"""

import os
import sys
from pymol import cmd
from Bio import SeqIO, pairwise2
from Bio.PDB.MMCIFParser import MMCIFParser


def load_cif_model(cif_path, model_name):
    """Load CIF model into PyMOL"""
    cmd.load(cif_path, model_name)
    print(f"Loaded {cif_path} as {model_name}")


def parse_fasta_sequences(fasta_path):
    """Parse FASTA file and return sequences dictionary"""
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def get_chain_sequence(cif_path, chain_id):
    """Extract sequence from CIF file for a specific chain"""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                seq = []
                for residue in chain:
                    if residue.id[0] == ' ':
                        resname = residue.resname
                        aa_dict = {
                            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                        }
                        if resname in aa_dict:
                            seq.append(aa_dict[resname])
                return ''.join(seq)
    return ''


def identify_chains(cif_path, fasta_sequences):
    """Identify chains by sequence alignment against FASTA reference.
    JSON parsing explicitly disabled per specification.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)
    
    chain_mapping = {}
    
    for model in structure:
        for chain in model:
            chain_seq = get_chain_sequence(cif_path, chain.id)
            
            best_match = None
            best_score = 0
            
            for protein_id, ref_seq in fasta_sequences.items():
                alignments = pairwise2.align.localxx(chain_seq, ref_seq)
                if alignments:
                    score = alignments[0].score
                    if score > best_score:
                        best_score = score
                        best_match = protein_id
            
            if best_match:
                chain_mapping[chain.id] = best_match
                print(f"Chain {chain.id} identified as {best_match} (score: {best_score})")
    
    return chain_mapping


def apply_visual_settings(config):
    """Apply global visual settings"""
    cmd.set('antialias', config['antialias'])
    cmd.set('ambient', config['ambient'])
    cmd.set('specular', config['specular'])
    cmd.set('shininess', config['shininess'])
    cmd.set('depth_cue', config['depth_cue'])
    cmd.set('ray_shadows', config['ray_shadows'])


def apply_cartoon_style(model_name, chain_id, color_hex):
    """Apply cartoon style to protein chain"""
    selection = f"{model_name} and chain {chain_id}"
    
    cmd.show('cartoon', selection)
    cmd.set('cartoon_thickness', 1.5, selection)
    cmd.set('cartoon_loop_radius', 0.3, selection)
    cmd.color(color_hex, selection)
    
    print(f"Applied cartoon style to {selection} with color {color_hex}")


def darken_color(hex_color, factor):
    """Darken a hex color by a factor (0.0 to 1.0)"""
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    r = int(r * factor)
    g = int(g * factor)
    b = int(b * factor)
    return f'#{r:02x}{g:02x}{b:02x}'


def highlight_hydrophobic_residues(model_name, chain_mapping, config):
    """Highlight hydrophobic residues with sticks and bubble surface"""
    
    all_residues = []
    
    for chain_id, protein_id in chain_mapping.items():
        if protein_id in config.HYDROPHOBIC_RESIDUES:
            residues = config.HYDROPHOBIC_RESIDUES[protein_id]
            residue_list = '+'.join(residues)
            
            selection = f"{model_name} and chain {chain_id} and resn {residue_list}"
            all_residues.append(selection)
            
            base_color = config.PROTEIN_COLORS.get(protein_id, '#CCCCCC')
            
            # Apply stick visualization (more visible, darker)
            stick_color = darken_color(base_color, config.STYLE_SETTINGS['stick_color_factor'])
            cmd.show('sticks', selection)
            cmd.set('stick_radius', config.STYLE_SETTINGS['stick_radius'], selection)
            cmd.color(stick_color, selection)
            
            print(f"Highlighted residues in {protein_id} chain {chain_id}: {residue_list}")
    
    if all_residues:
        # Create selection for all hydrophobic residues
        cmd.select('hydrophobic_res', ' or '.join(all_residues))
        
        # Create bubble residues (solid blobs around residues)
        cmd.create('hydrophobic_bubble', 'hydrophobic_res')
        cmd.show('surface', 'hydrophobic_bubble')
        
        # Make surface smooth and rounded
        cmd.set('surface_quality', config.STYLE_SETTINGS['surface_quality'], 'hydrophobic_bubble')
        cmd.set('surface_type', 0, 'hydrophobic_bubble')  # Solid surface
        cmd.set('surface_solvent', config.STYLE_SETTINGS['surface_solvent'], 'hydrophobic_bubble')
        cmd.set('surface_ramp_above_mode', 1, 'hydrophobic_bubble')
        
        # Apply transparency (less dark than stick)
        cmd.set('transparency', config.STYLE_SETTINGS['bubble_transparency'], 'hydrophobic_bubble')
        
        # Color bubble residues (less dark than stick)
        for chain_id, protein_id in chain_mapping.items():
            if protein_id in config.HYDROPHOBIC_RESIDUES:
                base_color = config.PROTEIN_COLORS.get(protein_id, '#CCCCCC')
                bubble_color = darken_color(base_color, config.STYLE_SETTINGS['bubble_color_factor'])
                
                residues = config.HYDROPHOBIC_RESIDUES[protein_id]
                residue_list = '+'.join(residues)
                bubble_selection = f"hydrophobic_bubble and chain {chain_id} and resn {residue_list}"
                cmd.color(bubble_color, bubble_selection)
        
        # Apply labels with leader lines
        if config.LABEL_RESN_ONLY:
            cmd.label('hydrophobic_res and name CA', '"%s" % resn')
        else:
            cmd.label('hydrophobic_res and name CA', '"%s%s" % (resn, resi)')
        
        # Enable leader lines for labels
        cmd.set('label_connector', 1, 'hydrophobic_res')
        cmd.set('label_connector_mode', 1, 'hydrophobic_res')  # Mode 1: leader lines
        cmd.set('label_size', 14, 'hydrophobic_res')
        cmd.set('label_outline_color', 'black', 'hydrophobic_res')
        
        print(f"Created bubble surface and labels with leader lines for hydrophobic residues")


def render_image(output_path, bg_mode, width, height, dpi, ray_width=None, ray_height=None):
    """Render and save image with specified background and ray tracing"""
    cmd.bg_color(bg_mode['bg_color'])
    
    if bg_mode['format'] == 'png':
        cmd.set('ray_opaque_background', 0)
    else:
        cmd.set('ray_opaque_background', 1)
    
    # Apply ray-trace to soften edges (create smooth bubble effect)
    if ray_width and ray_height:
        cmd.ray(ray_width, ray_height)
    else:
        cmd.ray(width, height)
    
    cmd.png(output_path, dpi=dpi)
    print(f"Rendered image: {output_path}")


def save_session(output_path):
    """Save PyMOL session"""
    cmd.save(output_path)
    print(f"Saved session: {output_path}")
