#!/usr/bin/env python3
"""
Contact Map Generator Module for GROMACS PPI Analysis

Generates inter-chain contact maps from structure files.
"""

import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np


def get_residue_coords(gro_file: Path) -> Dict[int, Dict]:
    """
    Extract CA coordinates per residue from GRO file.
    
    Args:
        gro_file: Path to GRO file
        
    Returns:
        Dictionary mapping residue number to {'name': str, 'coords': np.array}
    """
    residues = {}
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header (2 lines) and box vector (last line)
    for line in lines[2:-1]:
        if len(line) >= 44:
            try:
                resnum = int(line[0:5].strip())
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                
                # Use CA atoms for residue position
                if atomname == 'CA':
                    x = float(line[20:28])
                    y = float(line[28:36])
                    z = float(line[36:44])
                    residues[resnum] = {
                        'name': resname, 
                        'coords': np.array([x, y, z])
                    }
            except (ValueError, IndexError):
                continue
    
    return residues


def get_chain_residues(pdb_file: Path) -> Tuple[List[int], List[int]]:
    """
    Get chain A and chain B residue lists from PDB file.
    
    Because GROMACS renumbers residues continuously, we need to 
    return the continuous residue numbers, not the per-chain ones.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        Tuple of (chain_a_residues, chain_b_residues) with continuous numbering
    """
    chain_a_res = []
    chain_b_res = []
    continuous_resnum = 0
    last_resnum = None
    last_chain = None
    
    if not pdb_file.exists():
        return chain_a_res, chain_b_res
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    pdb_resnum = int(line[22:26].strip())
                    chain = line[21].strip() or 'A'
                    
                    # Track when residue changes (new residue or chain change)
                    if pdb_resnum != last_resnum or chain != last_chain:
                        continuous_resnum += 1
                        last_resnum = pdb_resnum
                        last_chain = chain
                        
                        if chain == 'A':
                            chain_a_res.append(continuous_resnum)
                        else:
                            chain_b_res.append(continuous_resnum)
                except (ValueError, IndexError):
                    continue
    
    return chain_a_res, chain_b_res


def calculate_contact_map(residues: Dict[int, Dict], 
                          chain_a_res: List[int],
                          chain_b_res: List[int],
                          cutoff: float = 0.8) -> Tuple[np.ndarray, List[Tuple]]:
    """
    Calculate contact map between two chains.
    
    Args:
        residues: Dictionary of residue coordinates
        chain_a_res: List of residue numbers for chain A
        chain_b_res: List of residue numbers for chain B
        cutoff: Distance cutoff in nm
        
    Returns:
        Tuple of (contact_map_array, interface_pairs_list)
    """
    contact_map = np.zeros((len(chain_a_res), len(chain_b_res)))
    interface_pairs = []
    
    for i, res_a in enumerate(chain_a_res):
        for j, res_b in enumerate(chain_b_res):
            if res_a in residues and res_b in residues:
                dist = np.linalg.norm(
                    residues[res_a]['coords'] - residues[res_b]['coords']
                )
                if dist < cutoff:
                    contact_map[i, j] = 1
                    interface_pairs.append((
                        res_a, 
                        residues[res_a]['name'],
                        res_b, 
                        residues[res_b]['name'], 
                        dist
                    ))
    
    return contact_map, interface_pairs


def generate_contact_map(workdir: Path, 
                        gro_file: str = "em.gro",
                        pdb_file: str = "clean.pdb",
                        cutoff: float = 0.8) -> Dict:
    """
    Generate contact map and related outputs.
    
    Args:
        workdir: Working directory
        gro_file: Name of GRO structure file
        pdb_file: Name of PDB file for chain info
        cutoff: Contact distance cutoff in nm
        
    Returns:
        Dictionary with contact map info
    """
    workdir = Path(workdir)
    
    # Get residue coordinates
    residues = get_residue_coords(workdir / gro_file)
    
    if not residues:
        print("No residues found")
        return {'error': 'No residues found'}
    
    # Get chain residue lists with continuous numbering
    chain_a_res, chain_b_res = get_chain_residues(workdir / pdb_file)
    
    # If chain detection failed, split by residue number (middle split)
    if not chain_a_res or not chain_b_res:
        all_res = sorted(residues.keys())
        mid = len(all_res) // 2
        chain_a_res = all_res[:mid]
        chain_b_res = all_res[mid:]
        print(f"Using fallback chain split: A={len(chain_a_res)}, B={len(chain_b_res)}")
    
    print(f"Chain A residues: {len(chain_a_res)}")
    print(f"Chain B residues: {len(chain_b_res)}")
    
    # Calculate contact map
    contact_map, interface_pairs = calculate_contact_map(
        residues, chain_a_res, chain_b_res, cutoff
    )
    
    # Sort interface pairs by distance
    interface_pairs.sort(key=lambda x: x[4])
    
    # Create output directory
    analysis_dir = workdir / 'analysis'
    analysis_dir.mkdir(exist_ok=True)
    
    # Save contact map as text
    np.savetxt(analysis_dir / 'contact_map.txt', contact_map, fmt='%d')
    
    # Save interface residue pairs
    with open(analysis_dir / 'interface_residues.txt', 'w') as f:
        f.write(f"Interface Residue Pairs (distance < {cutoff} nm)\n")
        f.write("=" * 50 + "\n")
        f.write(f"{'ChainA':<15} {'ChainB':<15} {'Distance (nm)':<15}\n")
        f.write("-" * 50 + "\n")
        for res_a, name_a, res_b, name_b, dist in interface_pairs[:30]:
            f.write(f"{name_a}{res_a:<10} {name_b}{res_b:<10} {dist:.3f}\n")
        f.write(f"\nTotal interface residue pairs: {len(interface_pairs)}\n")
    
    # Save plotting data as CSV
    with open(analysis_dir / 'contact_map_plot.csv', 'w') as f:
        f.write('chainA_res,chainB_res,distance_nm,contact\n')
        for i, res_a in enumerate(chain_a_res):
            for j, res_b in enumerate(chain_b_res):
                if res_a in residues and res_b in residues:
                    dist = np.linalg.norm(
                        residues[res_a]['coords'] - residues[res_b]['coords']
                    )
                    contact = 1 if dist < cutoff else 0
                    f.write(f'{res_a},{res_b},{dist:.3f},{contact}\n')
    
    result = {
        'chain_a_residues': len(chain_a_res),
        'chain_b_residues': len(chain_b_res),
        'total_contacts': len(interface_pairs),
        'interface_pairs': interface_pairs[:30],  # Top 30
        'files': {
            'contact_map': str(analysis_dir / 'contact_map.txt'),
            'interface_residues': str(analysis_dir / 'interface_residues.txt'),
            'plot_data': str(analysis_dir / 'contact_map_plot.csv')
        }
    }
    
    print(f"Contact map saved: {analysis_dir / 'contact_map.txt'}")
    print(f"Interface residues saved: {analysis_dir / 'interface_residues.txt'}")
    print(f"Found {len(interface_pairs)} interface residue pairs")
    
    return result


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate contact map")
    parser.add_argument("--workdir", required=True, help="Working directory")
    parser.add_argument("--gro", default="em.gro", help="GRO file name")
    parser.add_argument("--pdb", default="clean.pdb", help="PDB file name")
    parser.add_argument("--cutoff", type=float, default=0.8, help="Contact cutoff (nm)")
    args = parser.parse_args()
    
    result = generate_contact_map(
        Path(args.workdir),
        args.gro,
        args.pdb,
        args.cutoff
    )
    
    print(f"\nChain A residues: {result['chain_a_residues']}")
    print(f"Chain B residues: {result['chain_b_residues']}")
    print(f"Total contacts: {result['total_contacts']}")
