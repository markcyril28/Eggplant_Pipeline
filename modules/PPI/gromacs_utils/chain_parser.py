#!/usr/bin/env python3
"""
Chain Parser Module for GROMACS Stability Analysis

Parses PDB files to extract chain information and create index groups.
"""

import os
import sys
import argparse


def get_chain_info(pdb_file: str, outdir: str) -> dict:
    """
    Parse chain information from a PDB file.
    
    Args:
        pdb_file: Path to the PDB file
        outdir: Output directory for chain_info.txt
        
    Returns:
        Dictionary with chain atoms and residues info
    """
    chain_atoms = {}
    chain_residues = {}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip() or 'A'
                resnum = int(line[22:26].strip())
                
                if chain not in chain_atoms:
                    chain_atoms[chain] = 0
                    chain_residues[chain] = set()
                
                chain_atoms[chain] += 1
                chain_residues[chain].add(resnum)

    # Write chain info
    with open(os.path.join(outdir, 'chain_info.txt'), 'w') as f:
        for chain in sorted(chain_atoms.keys()):
            res = sorted(chain_residues[chain])
            f.write(f"Chain {chain}: {chain_atoms[chain]} atoms, residues {min(res)}-{max(res)} ({len(res)} total)\n")
            print(f"  Chain {chain}: {chain_atoms[chain]} atoms, residues {min(res)}-{max(res)}")
    
    return {'atoms': chain_atoms, 'residues': chain_residues}


def parse_chains_for_index(pdb_file: str, gro_file: str, index_file: str) -> dict:
    """
    Parse chains from PDB and GRO files to create index groups.
    
    The GRO file loses chain information, so we use atom counts from the PDB
    to determine which atoms belong to which chain.
    
    Args:
        pdb_file: Path to the PDB file
        gro_file: Path to the GRO file  
        index_file: Path to the index.ndx file to append to
        
    Returns:
        Dictionary mapping chain IDs to atom indices
    """
    # Count atoms per chain from PDB (in order of appearance)
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
    
    # Count protein atoms in GRO (excluding solvent/ions)
    gro_protein_atoms = []
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines[2:-1], start=1):
        if len(line) >= 15:
            try:
                resname = line[5:10].strip()
                # Skip solvent and ions
                if resname in ['SOL', 'NA', 'CL', 'WAT', 'HOH', 'K', 'MG', 'CA', 'ZN']:
                    continue
                gro_protein_atoms.append(i)
            except (ValueError, IndexError):
                continue
    
    # Assign atoms to chains based on cumulative counts
    chain_atoms = {c: [] for c in chain_order}
    atom_idx = 0
    
    for chain in chain_order:
        count = chain_atom_counts[chain]
        for _ in range(count):
            if atom_idx < len(gro_protein_atoms):
                chain_atoms[chain].append(gro_protein_atoms[atom_idx])
                atom_idx += 1
    
    # Append chain groups to index file
    with open(index_file, "a") as f:
        for chain in chain_order:
            atoms = chain_atoms[chain]
            f.write(f"\n[ Chain{chain} ]\n")
            for i, atom in enumerate(atoms):
                f.write(f"{atom:6d}")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n")
    
    for chain in chain_order:
        print(f"  Added Chain{chain}: {len(chain_atoms[chain])} atoms")
    
    return chain_atoms


def main():
    """Command-line interface for chain parsing."""
    parser = argparse.ArgumentParser(
        description='Parse chain information from PDB files'
    )
    parser.add_argument('--pdb', required=True, help='Path to PDB file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--gro', help='Path to GRO file (for index creation)')
    parser.add_argument('--index', help='Path to index.ndx file to append to')
    parser.add_argument('--mode', choices=['info', 'index'], default='info',
                       help='Mode: info (chain info) or index (create index groups)')
    
    args = parser.parse_args()
    
    if args.mode == 'info':
        get_chain_info(args.pdb, args.outdir)
    elif args.mode == 'index':
        if not args.gro or not args.index:
            parser.error("--gro and --index required for index mode")
        parse_chains_for_index(args.pdb, args.gro, args.index)


if __name__ == '__main__':
    main()
