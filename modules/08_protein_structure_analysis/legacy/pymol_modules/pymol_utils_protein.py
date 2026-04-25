#!/usr/bin/env python3
"""
PyMOL Utility Functions
Functions for loading, styling, and rendering protein structures
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
                    if residue.id[0] == ' ':  # Standard residue
                        resname = residue.resname
                        # Convert 3-letter to 1-letter code
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
    """Identify chains by sequence alignment against FASTA reference"""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_path)
    
    chain_mapping = {}
    
    for model in structure:
        for chain in model:
            chain_seq = get_chain_sequence(cif_path, chain.id)
            
            best_match = None
            best_score = 0
            
            for protein_id, ref_seq in fasta_sequences.items():
                # Perform pairwise alignment
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


def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple (0-1 range)"""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def apply_cartoon_bubble_style(model_name, chain_id, color_hex, bubble_size=1.0):
    """Apply combined cartoon and bubble style with cartoon more visible than bubble"""
    selection = f"{model_name} and chain {chain_id}"
    
    # Convert hex to RGB and set color
    rgb = hex_to_rgb(color_hex)
    color_name = f"custom_color_{chain_id}"
    cmd.set_color(color_name, rgb)
    
    # Cartoon representation (more visible)
    cmd.show('cartoon', selection)
    cmd.set('cartoon_rect_length', 1.5, model_name)
    cmd.set('cartoon_oval_length', 1.5, model_name)
    cmd.set('cartoon_loop_radius', 0.3, model_name)
    cmd.color(color_name, selection)
    
    # Bubble representation (smooth rounded blobs, less dark than stick)
    cmd.show('surface', selection)
    cmd.set('surface_quality', 2, model_name)  # Smooth and rounded
    cmd.set('transparency', 0.15, selection)  # Less dark than stick
    cmd.set('solvent_radius', 1.4 * bubble_size, model_name)  # BUBBLE SIZE CONTROL
    
    print(f"Applied cartoon+bubble style to {selection} with color {color_hex}, bubble_size={bubble_size}")


def render_image(output_path, bg_mode, width, height, dpi, ray_width=2000, ray_height=1500):
    """Render and save image with specified background and ray tracing"""
    cmd.bg_color(bg_mode['bg_color'])
    
    if bg_mode['format'] == 'png':
        cmd.set('ray_opaque_background', 0)
    else:
        cmd.set('ray_opaque_background', 1)
    
    # Apply ray trace to soften edges
    cmd.ray(ray_width, ray_height)
    cmd.png(output_path, dpi=dpi)
    print(f"Rendered image: {output_path}")


def save_session(output_path):
    """Save PyMOL session"""
    cmd.save(output_path)
    print(f"Saved session: {output_path}")
