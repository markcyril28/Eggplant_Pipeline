#!/usr/bin/env python3
"""
Structure Converter Module for GROMACS PPI Analysis

Utilities for converting between molecular structure formats (CIF, PDB, GRO).
"""

import os
import sys
from pathlib import Path
from typing import Optional, Union


def convert_cif_to_pdb(cif_file: Union[str, Path], 
                        output_file: Optional[Union[str, Path]] = None,
                        quiet: bool = True) -> Optional[Path]:
    """
    Convert a CIF file to PDB format using BioPython.
    
    Args:
        cif_file: Path to input CIF file
        output_file: Optional output PDB file path (defaults to same name with .pdb)
        quiet: Suppress BioPython warnings
    
    Returns:
        Path to output PDB file or None if failed
    """
    try:
        from Bio.PDB import MMCIFParser, PDBIO
    except ImportError:
        print("Error: BioPython not installed. Install with: pip install biopython")
        return None
    
    cif_path = Path(cif_file)
    if not cif_path.exists():
        print(f"Error: File not found: {cif_file}")
        return None
    
    # Determine output path
    if output_file:
        pdb_path = Path(output_file)
    else:
        pdb_path = cif_path.with_suffix('.pdb')
    
    # Ensure output directory exists
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Parse and convert
    try:
        parser = MMCIFParser(QUIET=quiet)
        structure = parser.get_structure(cif_path.stem, str(cif_path))
        
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(str(pdb_path))
        
        return pdb_path
    except Exception as e:
        print(f"Error converting {cif_file}: {e}")
        return None


def clean_pdb_for_gromacs(input_file: Union[str, Path],
                          output_file: Union[str, Path]) -> Path:
    """
    Clean a PDB file for GROMACS compatibility.
    
    - Keeps only ATOM, TER, END records
    - Converts non-standard histidine names (HSD, HSE, HSP) to HIS
    
    Args:
        input_file: Path to input PDB file
        output_file: Path to output cleaned PDB file
    
    Returns:
        Path to cleaned PDB file
    """
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            # Only keep ATOM, TER, END records
            if line.startswith(('ATOM', 'TER', 'END')):
                # Fix non-standard histidine names
                line = line.replace('HSD', 'HIS')
                line = line.replace('HSE', 'HIS')
                line = line.replace('HSP', 'HIS')
                f_out.write(line)
    
    return output_path


def prepare_structure(input_file: Union[str, Path],
                      output_file: Union[str, Path]) -> Path:
    """
    Prepare a structure file for GROMACS.
    
    Handles both CIF and PDB input files, converts and cleans as needed.
    
    Args:
        input_file: Path to input structure (CIF or PDB)
        output_file: Path to output cleaned PDB file
    
    Returns:
        Path to prepared PDB file
    """
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    # Handle CIF files
    if input_path.suffix.lower() == '.cif':
        temp_pdb = output_path.with_suffix('.temp.pdb')
        pdb_result = convert_cif_to_pdb(input_path, temp_pdb)
        if pdb_result is None:
            raise RuntimeError(f"Failed to convert CIF: {input_path}")
        clean_pdb_for_gromacs(temp_pdb, output_path)
        temp_pdb.unlink()  # Remove temporary file
    else:
        # Assume PDB, just clean it
        clean_pdb_for_gromacs(input_path, output_path)
    
    return output_path


def get_structure_info(pdb_file: Union[str, Path]) -> dict:
    """
    Get basic information about a PDB structure.
    
    Args:
        pdb_file: Path to PDB file
    
    Returns:
        Dictionary with structure info (chains, residues, atoms)
    """
    pdb_path = Path(pdb_file)
    info = {
        'chains': set(),
        'residues': 0,
        'atoms': 0,
        'chain_residues': {}
    }
    
    current_chain = None
    residues_seen = set()
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = line[22:26].strip()
                resname = line[17:20].strip()
                
                info['atoms'] += 1
                info['chains'].add(chain)
                
                res_key = f"{chain}_{resnum}"
                if res_key not in residues_seen:
                    residues_seen.add(res_key)
                    info['residues'] += 1
                    if chain not in info['chain_residues']:
                        info['chain_residues'][chain] = 0
                    info['chain_residues'][chain] += 1
    
    info['chains'] = sorted(list(info['chains']))
    return info


# Command-line interface
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert and prepare structure files')
    parser.add_argument('input', help='Input structure file (CIF or PDB)')
    parser.add_argument('-o', '--output', help='Output PDB file')
    parser.add_argument('--info', action='store_true', help='Print structure info only')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if args.info:
        if not input_path.exists():
            print(f"Error: File not found: {args.input}")
            sys.exit(1)
        info = get_structure_info(input_path)
        print(f"Structure: {input_path.name}")
        print(f"  Chains: {', '.join(info['chains'])}")
        print(f"  Residues: {info['residues']}")
        print(f"  Atoms: {info['atoms']}")
        for chain, count in info['chain_residues'].items():
            print(f"  Chain {chain}: {count} residues")
    else:
        output_path = Path(args.output) if args.output else input_path.with_suffix('.clean.pdb')
        result = prepare_structure(input_path, output_path)
        print(f"Prepared: {result}")
