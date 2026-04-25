#!/usr/bin/env python3
"""
Chain Indexer Module for GROMACS MD Simulations

Creates chain-based index groups for GROMACS simulations by parsing
PDB and GRO files to determine chain boundaries.
"""

import os
import argparse
from typing import Dict, List, Set, Tuple


def parse_gro_for_chains(
    gro_file: str, 
    pdb_file: str
) -> Tuple[List[int], List[int], List[int]]:
    """
    Parse GRO file and create chain-based index groups.
    Uses original PDB to determine chain boundaries by atom count.
    
    Note: GRO files lose chain information and renumber residues sequentially,
    so we use cumulative atom counts from the PDB to map atoms to chains.
    
    Args:
        gro_file: Path to the GRO file
        pdb_file: Path to the original PDB file
        
    Returns:
        Tuple of (chain_a_atoms, chain_b_atoms, non_protein_atoms) lists
    """
    # Count atoms per chain from PDB (preserving order)
    chain_atom_counts = {}
    chain_order = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip() or 'A'
                if chain not in chain_atom_counts:
                    chain_atom_counts[chain] = 0
                    chain_order.append(chain)
                chain_atom_counts[chain] += 1
    
    for chain in chain_order:
        print(f"Chain {chain}: {chain_atom_counts[chain]} atoms in PDB")
    
    # Read GRO file and collect protein atom indices and non-protein atoms
    protein_atoms: List[int] = []
    non_protein_atoms: List[int] = []
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # Non-protein residue names (water, ions, etc.)
    non_protein_resnames = {'SOL', 'NA', 'CL', 'WAT', 'HOH', 'Na+', 'Cl-', 'K', 'MG', 'CA', 'ZN', 'K+', 'Mg2+', 'Ca2+', 'Zn2+'}
    
    for i, line in enumerate(lines[2:-1], start=1):  # Skip header and box line
        if len(line) >= 15:
            try:
                resname = line[5:10].strip()
                # Separate protein from non-protein atoms
                if resname in non_protein_resnames:
                    non_protein_atoms.append(i)
                else:
                    protein_atoms.append(i)
            except (ValueError, IndexError):
                continue
    
    # Assign atoms to chains based on cumulative counts
    chain_a_atoms: List[int] = []
    chain_b_atoms: List[int] = []
    
    atom_idx = 0
    for chain in chain_order:
        count = chain_atom_counts[chain]
        for _ in range(count):
            if atom_idx < len(protein_atoms):
                if chain == 'A' or (len(chain_order) == 1):
                    chain_a_atoms.append(protein_atoms[atom_idx])
                elif chain == 'B':
                    chain_b_atoms.append(protein_atoms[atom_idx])
                # For other chains (C, D, etc.), add to chain B for simplicity
                else:
                    chain_b_atoms.append(protein_atoms[atom_idx])
                atom_idx += 1
    
    return chain_a_atoms, chain_b_atoms, non_protein_atoms


def create_chain_index(
    gro_file: str,
    pdb_file: str,
    index_file: str
) -> Tuple[int, int]:
    """
    Create chain-specific index groups and append to existing index file.
    Also creates a Non-Protein group for temperature coupling.
    
    Args:
        gro_file: Path to the GRO file
        pdb_file: Path to the original PDB file
        index_file: Path to the index.ndx file to append to
        
    Returns:
        Tuple of (chain_a_count, chain_b_count) atom counts
    """
    if not os.path.exists(pdb_file):
        print(f"Warning: PDB file not found: {pdb_file}")
        return 0, 0
        
    if not os.path.exists(gro_file):
        print(f"Warning: GRO file not found: {gro_file}")
        return 0, 0
    
    chain_a, chain_b, non_protein = parse_gro_for_chains(gro_file, pdb_file)
    
    # Append to existing index file
    with open(index_file, "a") as f:
        f.write("\n[ ChainA ]\n")
        for i, atom in enumerate(chain_a):
            f.write(f"{atom:6d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        f.write("\n")
        
        f.write("\n[ ChainB ]\n")
        for i, atom in enumerate(chain_b):
            f.write(f"{atom:6d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        f.write("\n")
        
        # Add Non-Protein group (water + ions) for temperature coupling
        f.write("\n[ Non-Protein ]\n")
        for i, atom in enumerate(non_protein):
            f.write(f"{atom:6d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        f.write("\n")
    
    print(f"Added ChainA ({len(chain_a)} atoms), ChainB ({len(chain_b)} atoms), and Non-Protein ({len(non_protein)} atoms) to {index_file}")
    
    return len(chain_a), len(chain_b)


def main():
    """Command-line interface for chain indexing."""
    parser = argparse.ArgumentParser(
        description='Create chain-based index groups for GROMACS'
    )
    parser.add_argument('--gro', required=True,
                       help='Path to GRO file')
    parser.add_argument('--pdb', required=True,
                       help='Path to original PDB file')
    parser.add_argument('--index', required=True,
                       help='Path to index.ndx file to append to')
    
    args = parser.parse_args()
    
    create_chain_index(args.gro, args.pdb, args.index)


if __name__ == '__main__':
    main()
