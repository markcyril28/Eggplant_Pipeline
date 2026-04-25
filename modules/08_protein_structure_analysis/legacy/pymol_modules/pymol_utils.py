#!/usr/bin/env python3
"""
PyMOL Utility Functions
Shared utility functions for all fold visualizations

This module provides:
  - Hex to RGB color conversion for PyMOL
  - Structure-based chain identification (primary method)
  - Chain labeling utilities
  - Residue highlighting with bubble surfaces
  - Ray-traced rendering for smooth visualization
"""

from pymol import cmd
from Bio.Align import PairwiseAligner
# import json  # JSON method commented out
import os
import logging


def convert_aa_codes(codes, reverse=False):
    """
    Convert amino acid codes between single-letter and three-letter formats.
    
    Args:
        codes (str or list): Single code(s) or list of codes to convert
        reverse (bool): If True, convert three-letter to single-letter; 
                       if False (default), convert single-letter to three-letter
        
    Returns:
        str or list: Converted code(s) in the same format as input
        
    Examples:
        >>> convert_aa_codes('L')  # Returns 'LEU'
        >>> convert_aa_codes(['L', 'I', 'V'])  # Returns ['LEU', 'ILE', 'VAL']
        >>> convert_aa_codes('LEU', reverse=True)  # Returns 'L'
    """
    # Import here to avoid circular dependency
    from pymol_config import AMINO_ACID_CODES, AMINO_ACID_CODES_REVERSE
    
    mapping = AMINO_ACID_CODES_REVERSE if reverse else AMINO_ACID_CODES
    
    if isinstance(codes, str):
        return mapping.get(codes.upper(), codes)
    elif isinstance(codes, list):
        return [mapping.get(code.upper(), code) for code in codes]
    else:
        return codes


def hex_to_rgb(hex_color):
    """
    Convert hex color to RGB tuple (0-1 range) for PyMOL.
    
    Args:
        hex_color (str): Hex color string (e.g., '#00BCD4' or '00BCD4')
        
    Returns:
        list: RGB values as [R, G, B] in 0-1 range
        
    Example:
        >>> hex_to_rgb('#00BCD4')
        [0.0, 0.7372549019607844, 0.8313725490196079]
    """
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return [r/255.0, g/255.0, b/255.0]


# ============================================================
# JSON-BASED CHAIN IDENTIFICATION (COMMENTED OUT)
# ============================================================
# def load_job_request_json(json_file):
#     """
#     Load and parse the AlphaFold job request JSON file.
#     
#     Args:
#         json_file (str): Path to the job request JSON file
#         
#     Returns:
#         dict or None: Parsed JSON data, or None if file doesn't exist or has errors
#     """
#     if not os.path.exists(json_file):
#         logging.debug(f"Job request JSON not found: {json_file}")
#         return None
#     
#     try:
#         with open(json_file, 'r') as f:
#             data = json.load(f)
#         return data
#     except json.JSONDecodeError as e:
#         logging.warning(f"Invalid JSON format in {json_file}: {e}")
#         return None
#     except Exception as e:
#         logging.warning(f"Error loading JSON file {json_file}: {e}")
#         return None


# def identify_chains_from_json(json_file, fasta_sequences):
#     """
#     Identify chains from AlphaFold job request JSON file (OPTIONAL FAST METHOD).
#     This is approximately 10x faster than structure alignment method.
#     Falls back to structure alignment if JSON is unavailable.
#     
#     Args:
#         json_file (str): Path to AlphaFold job request JSON file
#         fasta_sequences (dict): Dictionary mapping protein names to sequences
#         
#     Returns:
#         dict or None: Dictionary mapping chain IDs (A, B, C, etc.) to protein names,
#                      or None if JSON is missing/invalid or chains cannot be identified
#     """
#     data = load_job_request_json(json_file)
#     if not data or len(data) == 0:
#         return None
#     
#     # Extract sequences from JSON (typically in data[0]['sequences'])
#     try:
#         job_sequences = data[0].get('sequences', [])
#     except (KeyError, IndexError, TypeError) as e:
#         logging.warning(f"Unable to extract sequences from JSON: {e}")
#         return None
#     
#     chain_identities = {}
#     
#     # Map each sequence to a chain (A, B, C, ...)
#     chain_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
#     
#     for idx, seq_entry in enumerate(job_sequences):
#         if idx >= len(chain_letters):
#             logging.warning(f"More than {len(chain_letters)} chains detected, only processing first {len(chain_letters)}")
#             break
#             
#         chain = chain_letters[idx]
#         
#         # Extract protein sequence from JSON
#         try:
#             if 'proteinChain' in seq_entry:
#                 json_seq = seq_entry['proteinChain'].get('sequence', '')
#                 
#                 # Match with reference sequences
#                 for name, ref_seq in fasta_sequences.items():
#                     if json_seq == ref_seq:
#                         chain_identities[chain] = name
#                         logging.debug(f"Chain {chain} identified from JSON: {name}")
#                         break
#         except (KeyError, TypeError) as e:
#             logging.warning(f"Error processing sequence entry {idx}: {e}")
#             continue
#     
#     return chain_identities if chain_identities else None


# ============================================================
# STRUCTURE-BASED CHAIN IDENTIFICATION (PRIMARY METHOD)
# ============================================================


def get_chain_sequence(chain_id):
    """
    Extract sequence from a chain in the loaded structure.
    
    Args:
        chain_id (str): Chain identifier (e.g., 'A', 'B', 'C')
        
    Returns:
        str: Amino acid sequence of the chain (empty string if extraction fails)
    """
    try:
        seq = cmd.get_fastastr(f'chain {chain_id}')
        # Remove header line (first line starting with '>')
        seq_lines = seq.strip().split('\n')
        if len(seq_lines) > 1:
            return ''.join(seq_lines[1:])
        return ''
    except Exception as e:
        logging.warning(f"Error extracting sequence from chain {chain_id}: {e}")
        return ''


def identify_chain(chain_seq, fasta_sequences):
    """
    Identify chain by alignment with FASTA sequences (primary identification method).
    
    Args:
        chain_seq (str): Amino acid sequence from the structure
        fasta_sequences (dict): Dictionary mapping protein names to reference sequences
        
    Returns:
        str or None: Best matching protein name, or None if no good match found
    """
    if not chain_seq:
        return None
        
    best_match = None
    best_score = 0
    
    # Create aligner with global alignment mode
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    try:
        for name, ref_seq in fasta_sequences.items():
            # Use global alignment
            score = aligner.score(chain_seq, ref_seq)
            if score > best_score:
                best_score = score
                best_match = name
    except Exception as e:
        logging.warning(f"Error during sequence alignment: {e}")
        return None
    
    return best_match


def identify_chains_from_structure(object_name, fasta_sequences):
    """
    Identify chains by extracting sequences from structure and aligning (PRIMARY METHOD).
    Uses sequence alignment with reference FASTA to identify chains.
    
    Args:
        object_name (str): Name of the loaded PyMOL object
        fasta_sequences (dict): Dictionary mapping protein names to sequences
        
    Returns:
        dict: Dictionary mapping chain IDs to protein names (empty dict if no chains identified)
    """
    chain_identities = {}
    
    try:
        # Get all chains in the object
        stored_chains = []
        cmd.iterate(f'{object_name}', 'stored_chains.append(chain)', space={'stored_chains': stored_chains})
        unique_chains = list(set(stored_chains))
        
        for chain in unique_chains:
            if not chain:  # Skip empty chain identifiers
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


def label_chains(object_name, chain_identities):
    """
    Add labels to chains with full protein names.
    
    Args:
        object_name (str): Name of the loaded PyMOL object
        chain_identities (dict): Dictionary with chain info including 'name' key
        
    Note:
        Labels are added to the first CA atom (residue 1) of each chain.
        Label properties (size, color) are set globally.
    """
    try:
        for chain, info in chain_identities.items():
            # Get the center of mass for the chain
            selection = f'{object_name} and chain {chain}'
            
            # Create label with full protein name
            label_text = info['name']
            
            # Add label at chain center (first CA atom)
            cmd.label(f'{selection} and name CA and resi 1', f'"{label_text}"')
        
        # Set label properties globally (not per selection to avoid warnings)
        cmd.set('label_size', 20)
        cmd.set('label_color', 'white')
        
    except Exception as e:
        logging.warning(f"Error adding labels to chains: {e}")


def highlight_hydrophobic_residues(object_name, residue_names, stick_radius=0.2, label_resn_only=True, 
                                   chain_identities=None, protein_filter=None):
    """
    Highlight hydrophobic residues by showing them as sticks and optionally labeling them.
    DEPRECATED: Use highlight_specific_residues for position-specific highlighting.
    
    Args:
        object_name (str): Name of the loaded PyMOL object
        residue_names (list): List of residue names to highlight
        stick_radius (float): Radius of stick representation (default: 0.2)
        label_resn_only (bool): If True, show only residue name; if False, show name+number
        chain_identities (dict): Dictionary mapping chains to protein info
        protein_filter (str): Only highlight residues in chains of this protein type
    """
    try:
        converted_residues = []
        for resn in residue_names:
            if len(resn) == 1:
                converted = convert_aa_codes(resn)
                converted_residues.append(converted)
                logging.debug(f"  Converted {resn} -> {converted}")
            else:
                converted_residues.append(resn.upper())
        
        resn_str = '+'.join(converted_residues)
        
        if protein_filter and chain_identities:
            target_chains = [chain for chain, info in chain_identities.items() 
                           if info.get('type') == protein_filter]
            
            if not target_chains:
                logging.warning(f"  No chains found for protein type '{protein_filter}'")
                return
            
            chain_str = '+'.join(target_chains)
            selection_name = 'hydrophobic_res'
            cmd.select(selection_name, f'{object_name} and chain {chain_str} and resn {resn_str}')
            logging.info(f"  Highlighting residues in {protein_filter} chains: {', '.join(target_chains)}")
        else:
            selection_name = 'hydrophobic_res'
            cmd.select(selection_name, f'{object_name} and resn {resn_str}')
        
        cmd.show('sticks', selection_name)
        cmd.set('stick_radius', stick_radius, selection_name)
        
        if label_resn_only:
            cmd.label(f'{selection_name} and name CA', '"%s" % resn')
        else:
            cmd.label(f'{selection_name} and name CA', '"%s%s" % (resn, resi)')
        
        count = cmd.count_atoms(selection_name)
        filter_msg = f" (filtered to {protein_filter} chains)" if protein_filter else ""
        logging.info(f"  ✓ Highlighted {count} atoms from {len(converted_residues)} residue types: {', '.join(converted_residues)}{filter_msg}")
        
        cmd.deselect()
        
    except Exception as e:
        logging.warning(f"Error highlighting hydrophobic residues: {e}")


def highlight_specific_residues(object_name, residues_by_protein, chain_identities, 
                               stick_radius=0.2, label_resn_only=True,
                               create_bubble=False, bubble_quality=2, bubble_radius=3.0,
                               bubble_transparency=0.5, colors_highlighted_stick=None, 
                               colors_highlighted_bubble=None):
    """
    Highlight specific residues by position or by name in specific protein chains.
    Creates combined cartoon and bubble visualization where:
      - Cartoon is more visible than bubble
      - Bubble (surface) is less dark than stick
      - Smooth, rounded blobs around residues
    
    Args:
        object_name (str): Name of the loaded PyMOL object
        residues_by_protein (dict): Dict mapping protein types to residue specifications:
                                   - List of strings for name-only: ['GLN', 'LEU']
                                   - List of tuples for position-specific: [('GLN', 20), ('LEU', 24)]
        chain_identities (dict): Dictionary mapping chains to protein info (with 'type' key)
        stick_radius (float): Radius of stick representation
        label_resn_only (bool): If True, show only residue name; if False, show name+number
        create_bubble (bool): Create smooth surface around selected residues
        bubble_quality (int): Surface quality for bubble (0-4, higher = smoother)
        bubble_radius (float): Solvent radius for bubble surface (controls size)
        bubble_transparency (float): Transparency of bubble surface (0.0=opaque, 1.0=transparent)
        colors_highlighted_stick (dict): Darker colors for highlighted residue sticks
        colors_highlighted_bubble (dict): Medium-dark colors for highlighted residue bubbles
        
    Returns:
        list: List of selection names created for highlighted residues
    """
    try:
        all_selections = []
        
        for protein_type, residue_spec in residues_by_protein.items():
            # Find chains of this protein type
            target_chains = [chain for chain, info in chain_identities.items() 
                           if info.get('type') == protein_type]
            
            if not target_chains:
                logging.warning(f"  No chains found for protein type '{protein_type}'")
                continue
            
            logging.info(f"  Highlighting residues in {protein_type} chains: {', '.join(target_chains)}")
            chain_str = '+'.join(target_chains)
            
            # Check if residue_spec contains tuples (position-specific) or strings (name-only)
            if residue_spec and isinstance(residue_spec[0], tuple):
                # Position-specific mode
                for resname, position in residue_spec:
                    # Convert single-letter code if needed
                    if len(resname) == 1:
                        resname = convert_aa_codes(resname)
                    
                    # Create selection for this specific residue
                    selection_name = f'res_{protein_type}_{resname}{position}'
                    selection_cmd = f'{object_name} and chain {chain_str} and resn {resname} and resi {position}'
                    
                    cmd.select(selection_name, selection_cmd)
                    
                    # Check if selection found anything
                    count = cmd.count_atoms(selection_name)
                    if count == 0:
                        logging.warning(f"  ⚠ No atoms found for {resname} {position} in {protein_type}")
                        continue
                    
                    # Show as sticks
                    cmd.show('sticks', selection_name)
                    cmd.set('stick_radius', stick_radius, selection_name)
                    
                    # Apply darker color to sticks if provided
                    if colors_highlighted_stick:
                        color_hex = colors_highlighted_stick.get(protein_type)
                        if color_hex:
                            stick_color_name = f'{protein_type}_highlight_{resname}{position}'
                            rgb_values = hex_to_rgb(color_hex)
                            cmd.set_color(stick_color_name, rgb_values)
                            cmd.color(stick_color_name, selection_name)
                    
                    # Label the residue
                    if label_resn_only:
                        cmd.label(f'{selection_name} and name CA', '"%s" % resn')
                    else:
                        cmd.label(f'{selection_name} and name CA', '"%s%s" % (resn, resi)')
                    
                    # Create bubble surface if requested
                    if create_bubble:
                        bubble_name = f'bubble_{protein_type}_{resname}{position}'
                        cmd.create(bubble_name, selection_name)
                        cmd.show('surface', bubble_name)
                        cmd.set('surface_quality', bubble_quality, bubble_name)
                        cmd.set('solvent_radius', bubble_radius, bubble_name)
                        cmd.set('transparency', bubble_transparency, bubble_name)
                        
                        # Apply intermediate bubble color (lighter than sticks, darker than standard)
                        if colors_highlighted_bubble:
                            color_hex = colors_highlighted_bubble.get(protein_type)
                            if color_hex:
                                bubble_color_name = f'{protein_type}_bubble_{resname}{position}'
                                rgb_values = hex_to_rgb(color_hex)
                                cmd.set_color(bubble_color_name, rgb_values)
                                cmd.set('surface_color', bubble_color_name, bubble_name)
                        else:
                            # Fallback to chain color
                            chain_info = chain_identities.get(target_chains[0])
                            if chain_info:
                                color_name = f"{chain_info['type']}_color"
                                cmd.set('surface_color', color_name, bubble_name)
                        
                        all_selections.append(bubble_name)
                        logging.info(f"    ✓ Created bubble surface for {resname} {position}")
                    
                    all_selections.append(selection_name)
                    logging.info(f"    ✓ Highlighted {resname} {position} ({count} atoms)")
            
            else:
                # Name-only mode (all residues of specified types)
                converted_residues = []
                for resname in residue_spec:
                    # Convert single-letter code if needed
                    if len(resname) == 1:
                        resname = convert_aa_codes(resname)
                    converted_residues.append(resname.upper())
                
                # Build selection for all residues of these types
                resn_str = '+'.join(converted_residues)
                selection_name = f'res_{protein_type}_all'
                selection_cmd = f'{object_name} and chain {chain_str} and resn {resn_str}'
                
                cmd.select(selection_name, selection_cmd)
                
                # Check if selection found anything
                count = cmd.count_atoms(selection_name)
                if count == 0:
                    logging.warning(f"  ⚠ No atoms found for {', '.join(converted_residues)} in {protein_type}")
                    continue
                
                # Show as sticks
                cmd.show('sticks', selection_name)
                cmd.set('stick_radius', stick_radius, selection_name)
                
                # Apply darker color to sticks if provided
                if colors_highlighted_stick:
                    color_hex = colors_highlighted_stick.get(protein_type)
                    if color_hex:
                        stick_color_name = f'{protein_type}_highlight_all'
                        rgb_values = hex_to_rgb(color_hex)
                        cmd.set_color(stick_color_name, rgb_values)
                        cmd.color(stick_color_name, selection_name)
                
                # Label the residues
                if label_resn_only:
                    cmd.label(f'{selection_name} and name CA', '"%s" % resn')
                else:
                    cmd.label(f'{selection_name} and name CA', '"%s%s" % (resn, resi)')
                
                # Create bubble surface if requested
                if create_bubble:
                    bubble_name = f'bubble_{protein_type}_all'
                    cmd.create(bubble_name, selection_name)
                    cmd.show('surface', bubble_name)
                    cmd.set('surface_quality', bubble_quality, bubble_name)
                    cmd.set('solvent_radius', bubble_radius, bubble_name)
                    cmd.set('transparency', bubble_transparency, bubble_name)
                    
                    # Apply intermediate bubble color (lighter than sticks, darker than standard)
                    if colors_highlighted_bubble:
                        color_hex = colors_highlighted_bubble.get(protein_type)
                        if color_hex:
                            bubble_color_name = f'{protein_type}_bubble_all'
                            rgb_values = hex_to_rgb(color_hex)
                            cmd.set_color(bubble_color_name, rgb_values)
                            cmd.set('surface_color', bubble_color_name, bubble_name)
                    else:
                        # Fallback to chain color
                        chain_info = chain_identities.get(target_chains[0])
                        if chain_info:
                            color_name = f"{chain_info['type']}_color"
                            cmd.set('surface_color', color_name, bubble_name)
                    
                    all_selections.append(bubble_name)
                    logging.info(f"    ✓ Created bubble surface for all {', '.join(converted_residues)}")
                
                all_selections.append(selection_name)
                logging.info(f"    ✓ Highlighted all {', '.join(converted_residues)} ({count} atoms)")
        
        # Deselect to clean up visual
        cmd.deselect()
        
        return all_selections
        
    except Exception as e:
        logging.error(f"Error highlighting specific residues: {e}")
        return []


def center_on_residue(object_name, chain_identities, protein_type, resname, position, zoom_buffer=5.0):
    """
    Center camera view on a specific residue.
    
    Args:
        object_name (str): Name of the loaded PyMOL object
        chain_identities (dict): Dictionary mapping chains to protein info
        protein_type (str): Protein type to find (e.g., 'GRF')
        resname (str): Residue name (e.g., 'LEU')
        position (int): Residue position number
        zoom_buffer (float): Buffer distance for zoom (Angstroms)
    """
    try:
        # Find chains of this protein type
        target_chains = [chain for chain, info in chain_identities.items() 
                       if info.get('type') == protein_type]
        
        if not target_chains:
            logging.warning(f"  Cannot center: No chains found for protein type '{protein_type}'")
            return False
        
        # Convert single-letter code if needed
        if len(resname) == 1:
            resname = convert_aa_codes(resname)
        
        # Create selection for the center residue
        chain_str = '+'.join(target_chains)
        center_selection = f'{object_name} and chain {chain_str} and resn {resname} and resi {position}'
        
        # Check if residue exists
        count = cmd.count_atoms(center_selection)
        if count == 0:
            logging.warning(f"  Cannot center: Residue {resname} {position} not found in {protein_type}")
            return False
        
        # Center the view on this residue
        cmd.center(center_selection)
        cmd.zoom(center_selection, buffer=zoom_buffer)
        
        logging.info(f"  ✓ Centered view on {protein_type} {resname} {position}")
        return True
        
    except Exception as e:
        logging.error(f"Error centering on residue: {e}")
        return False
