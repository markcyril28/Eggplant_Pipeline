#!/usr/bin/env python3
"""
PyMOL Utility Functions for SmelGIF-SWI2 Visualization
Utility functions for loading models, identifying chains, and applying visualizations
"""

from pymol import cmd
from Bio.Align import PairwiseAligner
import os
import logging


def convert_aa_codes(codes, reverse=False):
    """Convert amino acid codes between single-letter and three-letter formats."""
    from pymol_config_gif_swi2 import AMINO_ACID_CODES, AMINO_ACID_CODES_REVERSE
    
    mapping = AMINO_ACID_CODES_REVERSE if reverse else AMINO_ACID_CODES
    
    if isinstance(codes, str):
        return mapping.get(codes.upper(), codes)
    elif isinstance(codes, list):
        return [mapping.get(code.upper(), code) for code in codes]
    else:
        return codes


def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple (0-1 range) for PyMOL."""
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return [r/255.0, g/255.0, b/255.0]


def get_chain_sequence(chain_id):
    """Extract sequence from a chain in the loaded structure."""
    try:
        seq = cmd.get_fastastr(f'chain {chain_id}')
        seq_lines = seq.strip().split('\n')
        if len(seq_lines) > 1:
            return ''.join(seq_lines[1:])
        return ''
    except Exception as e:
        logging.warning(f"Error extracting sequence from chain {chain_id}: {e}")
        return ''


def identify_chain(chain_seq, fasta_sequences):
    """Identify chain by alignment with FASTA sequences."""
    if not chain_seq:
        return None
        
    best_match = None
    best_score = 0
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    try:
        for name, ref_seq in fasta_sequences.items():
            score = aligner.score(chain_seq, ref_seq)
            if score > best_score:
                best_score = score
                best_match = name
    except Exception as e:
        logging.warning(f"Error during sequence alignment: {e}")
        return None
    
    return best_match


def identify_chains_from_structure(object_name, fasta_sequences):
    """Identify chains by extracting sequences from structure and aligning."""
    chain_identities = {}
    
    try:
        stored_chains = []
        cmd.iterate(f'{object_name}', 'stored_chains.append(chain)', space={'stored_chains': stored_chains})
        unique_chains = list(set(stored_chains))
        
        for chain in unique_chains:
            if not chain:
                continue
                
            chain_seq = get_chain_sequence(chain)
            if chain_seq:
                match = identify_chain(chain_seq, fasta_sequences)
                if match:
                    chain_identities[chain] = match
                    logging.debug(f"Chain {chain} identified from structure: {match}")
                else:
                    logging.warning(f"Chain {chain} could not be matched to any reference sequence")
            else:
                logging.warning(f"Could not extract sequence from chain {chain}")
                
    except Exception as e:
        logging.error(f"Error identifying chains from structure: {e}")
        return {}
    
    return chain_identities


def highlight_specific_residues(object_name, residues_by_protein, chain_identities, 
                               stick_radius=0.2, label_resn_only=True,
                               create_bubble=False, bubble_quality=2, bubble_radius=3.0,
                               bubble_transparency=0.5, colors_highlighted_stick=None, 
                               colors_highlighted_bubble=None):
    """Highlight specific residues by name in specific protein chains."""
    try:
        all_selections = []
        
        for protein_type, residue_spec in residues_by_protein.items():
            target_chains = [chain for chain, info in chain_identities.items() 
                           if info.get('type') == protein_type]
            
            if not target_chains:
                logging.warning(f"  No chains found for protein type '{protein_type}'")
                continue
            
            logging.info(f"  Highlighting residues in {protein_type} chains: {', '.join(target_chains)}")
            chain_str = '+'.join(target_chains)
            
            converted_residues = []
            for resname in residue_spec:
                if len(resname) == 1:
                    resname = convert_aa_codes(resname)
                converted_residues.append(resname.upper())
            
            resn_str = '+'.join(converted_residues)
            selection_name = f'res_{protein_type}_all'
            selection_cmd = f'{object_name} and chain {chain_str} and resn {resn_str}'
            
            cmd.select(selection_name, selection_cmd)
            
            count = cmd.count_atoms(selection_name)
            if count == 0:
                logging.warning(f"  ⚠ No atoms found for {', '.join(converted_residues)} in {protein_type}")
                continue
            
            cmd.show('sticks', selection_name)
            cmd.set('stick_radius', stick_radius, selection_name)
            
            if colors_highlighted_stick:
                color_hex = colors_highlighted_stick.get(protein_type)
                if color_hex:
                    stick_color_name = f'{protein_type}_highlight_all'
                    rgb_values = hex_to_rgb(color_hex)
                    cmd.set_color(stick_color_name, rgb_values)
                    cmd.color(stick_color_name, selection_name)
            
            if label_resn_only:
                cmd.label(f'{selection_name} and name CA', '"%s" % resn')
            else:
                cmd.label(f'{selection_name} and name CA', '"%s%s" % (resn, resi)')
            
            if create_bubble:
                bubble_name = f'bubble_{protein_type}_all'
                cmd.create(bubble_name, selection_name)
                cmd.show('surface', bubble_name)
                cmd.set('surface_quality', bubble_quality, bubble_name)
                cmd.set('solvent_radius', bubble_radius, bubble_name)
                cmd.set('transparency', bubble_transparency, bubble_name)
                
                if colors_highlighted_bubble:
                    color_hex = colors_highlighted_bubble.get(protein_type)
                    if color_hex:
                        bubble_color_name = f'{protein_type}_bubble_all'
                        rgb_values = hex_to_rgb(color_hex)
                        cmd.set_color(bubble_color_name, rgb_values)
                        cmd.set('surface_color', bubble_color_name, bubble_name)
                
                all_selections.append(bubble_name)
                logging.info(f"    ✓ Created bubble surface for all {', '.join(converted_residues)}")
            
            all_selections.append(selection_name)
            logging.info(f"    ✓ Highlighted all {', '.join(converted_residues)} ({count} atoms)")
        
        cmd.deselect()
        return all_selections
        
    except Exception as e:
        logging.error(f"Error highlighting specific residues: {e}")
        return []


def center_on_residue(object_name, chain_identities, protein_type, resname, position, zoom_buffer=5.0):
    """Center camera view on a specific residue."""
    try:
        target_chains = [chain for chain, info in chain_identities.items() 
                       if info.get('type') == protein_type]
        
        if not target_chains:
            logging.warning(f"  Cannot center: No chains found for protein type '{protein_type}'")
            return False
        
        if len(resname) == 1:
            resname = convert_aa_codes(resname)
        
        chain_str = '+'.join(target_chains)
        center_selection = f'{object_name} and chain {chain_str} and resn {resname} and resi {position}'
        
        count = cmd.count_atoms(center_selection)
        if count == 0:
            logging.warning(f"  Cannot center: Residue {resname} {position} not found in {protein_type}")
            return False
        
        cmd.center(center_selection)
        cmd.zoom(center_selection, buffer=zoom_buffer)
        
        logging.info(f"  ✓ Centered view on {protein_type} {resname} {position}")
        return True
        
    except Exception as e:
        logging.error(f"Error centering on residue: {e}")
        return False
