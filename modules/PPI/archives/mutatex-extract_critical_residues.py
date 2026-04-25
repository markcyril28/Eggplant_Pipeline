#!/usr/bin/env python3
"""
Extract Critical Residues from MutateX Results
Creates FASTA files and position lists for MSA analysis
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def parse_mutatex_results(results_dir):
    """Parse MutateX final_averages CSV files."""
    results = {}
    
    # Look for interface DDG results
    interface_dir = Path(results_dir) / "results" / "interface_ddgs" / "final_averages"
    if interface_dir.exists():
        for csv_file in interface_dir.glob("*.csv"):
            chain_pair = csv_file.stem
            df = pd.read_csv(csv_file)
            results[f"interface_{chain_pair}"] = df
    
    # Look for mutation DDG results
    mutate_dir = Path(results_dir) / "results" / "mutate_ddgs" / "final_averages"
    if mutate_dir.exists():
        for csv_file in mutate_dir.glob("*.csv"):
            chain = csv_file.stem
            df = pd.read_csv(csv_file)
            results[f"stability_{chain}"] = df
    
    return results


def extract_sequence_from_pdb(pdb_file, chain_id):
    """Extract sequence from PDB file for a specific chain."""
    from Bio.PDB import PDBParser
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Standard 3-letter to 1-letter code
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    sequence = []
    residue_numbers = []
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residue
                        res_name = residue.resname
                        if res_name in three_to_one:
                            sequence.append(three_to_one[res_name])
                            residue_numbers.append(residue.id[1])
    
    return ''.join(sequence), residue_numbers


def identify_critical_residues(df, threshold_destabilizing=2.0, threshold_stabilizing=-2.0, top_n=20):
    """
    Identify critical residues based on DDG thresholds.
    
    Parameters:
    - threshold_destabilizing: DDG > this value (default 2.0 kcal/mol)
    - threshold_stabilizing: DDG < this value (default -2.0 kcal/mol)
    - top_n: Number of top residues to include
    """
    critical = {
        'destabilizing': [],
        'stabilizing': [],
        'interface': [],
        'hotspot': []
    }
    
    # Get average DDG per position
    if 'position' in df.columns and 'average_total_energy' in df.columns:
        position_avg = df.groupby('position')['average_total_energy'].mean()
        
        # Destabilizing mutations (positive DDG)
        destabilizing = position_avg[position_avg > threshold_destabilizing].sort_values(ascending=False)
        critical['destabilizing'] = destabilizing.head(top_n).index.tolist()
        
        # Stabilizing mutations (negative DDG)
        stabilizing = position_avg[position_avg < threshold_stabilizing].sort_values()
        critical['stabilizing'] = stabilizing.head(top_n).index.tolist()
        
        # Hotspots (high absolute DDG change)
        hotspots = position_avg.abs().sort_values(ascending=False)
        critical['hotspot'] = hotspots.head(top_n).index.tolist()
    
    return critical


def create_fasta_files(results_dir, pdb_file, output_dir, thresholds):
    """Create FASTA files with critical residues annotated."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse MutateX results
    results = parse_mutatex_results(results_dir)
    
    if not results:
        print(f"[ERROR] No MutateX results found in {results_dir}")
        return
    
    all_critical_positions = {}
    
    # Process each result type
    for result_type, df in results.items():
        print(f"\n[INFO] Processing {result_type}")
        
        # Extract chain ID
        if 'interface' in result_type:
            # Format: interface_chainA_chainB
            chains = result_type.replace('interface_', '').split('_')
            chain_id = chains[0] if chains else 'A'
        else:
            # Format: stability_chainA
            chain_id = result_type.replace('stability_', '')
        
        # Identify critical residues
        critical = identify_critical_residues(
            df, 
            threshold_destabilizing=thresholds['destabilizing'],
            threshold_stabilizing=thresholds['stabilizing'],
            top_n=thresholds['top_n']
        )
        
        if chain_id not in all_critical_positions:
            all_critical_positions[chain_id] = {
                'destabilizing': set(),
                'stabilizing': set(),
                'hotspot': set(),
                'interface': set()
            }
        
        for category, positions in critical.items():
            all_critical_positions[chain_id][category].update(positions)
            print(f"  {category.capitalize()}: {len(positions)} positions")
    
    # Extract sequences and create FASTA files
    for chain_id, critical_sets in all_critical_positions.items():
        try:
            sequence, residue_numbers = extract_sequence_from_pdb(pdb_file, chain_id)
            
            if not sequence:
                print(f"[WARN] No sequence found for chain {chain_id}")
                continue
            
            # Create position mapping (residue number -> sequence index)
            pos_to_idx = {num: idx for idx, num in enumerate(residue_numbers)}
            
            # Full sequence FASTA
            record = SeqRecord(
                Seq(sequence),
                id=f"{Path(pdb_file).stem}_chain{chain_id}",
                description=f"Full sequence | Chain {chain_id}"
            )
            
            fasta_file = output_dir / f"chain_{chain_id}_full_sequence.fasta"
            SeqIO.write(record, fasta_file, "fasta")
            print(f"\n[SAVED] {fasta_file}")
            
            # Create annotated position files for each category
            for category, positions in critical_sets.items():
                if not positions:
                    continue
                
                # Position list file (for MSA highlighting)
                pos_file = output_dir / f"chain_{chain_id}_{category}_positions.txt"
                sorted_positions = sorted(positions)
                with open(pos_file, 'w') as f:
                    f.write(f"# Critical {category} residues for chain {chain_id}\n")
                    f.write(f"# Total: {len(sorted_positions)} positions\n")
                    f.write("# Format: Position\tResidue\tIndex_in_sequence\n")
                    for pos in sorted_positions:
                        if pos in pos_to_idx:
                            idx = pos_to_idx[pos]
                            residue = sequence[idx] if idx < len(sequence) else 'X'
                            f.write(f"{pos}\t{residue}\t{idx}\n")
                
                print(f"[SAVED] {pos_file}")
                
                # Create highlighted sequence (lowercase = normal, UPPERCASE = critical)
                highlighted_seq = list(sequence.lower())
                for pos in positions:
                    if pos in pos_to_idx:
                        idx = pos_to_idx[pos]
                        if idx < len(highlighted_seq):
                            highlighted_seq[idx] = highlighted_seq[idx].upper()
                
                highlighted_record = SeqRecord(
                    Seq(''.join(highlighted_seq)),
                    id=f"{Path(pdb_file).stem}_chain{chain_id}_{category}",
                    description=f"{category.capitalize()} residues highlighted (UPPERCASE) | Chain {chain_id}"
                )
                
                highlight_file = output_dir / f"chain_{chain_id}_{category}_highlighted.fasta"
                SeqIO.write(highlighted_record, highlight_file, "fasta")
                print(f"[SAVED] {highlight_file}")
            
            # Create combined critical residues file
            all_critical = set()
            for positions in critical_sets.values():
                all_critical.update(positions)
            
            if all_critical:
                combined_file = output_dir / f"chain_{chain_id}_all_critical_positions.txt"
                with open(combined_file, 'w') as f:
                    f.write(f"# All critical residues for chain {chain_id}\n")
                    f.write(f"# Total: {len(all_critical)} unique positions\n")
                    f.write("# Format: Position\tResidue\tCategories\n")
                    
                    for pos in sorted(all_critical):
                        if pos in pos_to_idx:
                            idx = pos_to_idx[pos]
                            residue = sequence[idx] if idx < len(sequence) else 'X'
                            categories = [cat for cat, positions in critical_sets.items() if pos in positions]
                            f.write(f"{pos}\t{residue}\t{','.join(categories)}\n")
                
                print(f"[SAVED] {combined_file}")
        
        except Exception as e:
            print(f"[ERROR] Failed to process chain {chain_id}: {e}")
            continue
    
    print(f"\n[COMPLETE] All files saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract critical residues from MutateX results for MSA analysis"
    )
    parser.add_argument(
        "results_dir",
        help="Directory containing MutateX results"
    )
    parser.add_argument(
        "pdb_file",
        help="Original PDB file used for analysis"
    )
    parser.add_argument(
        "-o", "--output",
        default="critical_residues",
        help="Output directory (default: critical_residues)"
    )
    parser.add_argument(
        "--destabilizing",
        type=float,
        default=2.0,
        help="DDG threshold for destabilizing mutations (default: 2.0 kcal/mol)"
    )
    parser.add_argument(
        "--stabilizing",
        type=float,
        default=-2.0,
        help="DDG threshold for stabilizing mutations (default: -2.0 kcal/mol)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top residues to extract per category (default: 20)"
    )
    
    args = parser.parse_args()
    
    thresholds = {
        'destabilizing': args.destabilizing,
        'stabilizing': args.stabilizing,
        'top_n': args.top_n
    }
    
    print(f"[INFO] Extracting critical residues from {args.results_dir}")
    print(f"[INFO] Using PDB file: {args.pdb_file}")
    print(f"[INFO] Thresholds: destabilizing >{thresholds['destabilizing']}, "
          f"stabilizing <{thresholds['stabilizing']}, top N={thresholds['top_n']}")
    
    create_fasta_files(args.results_dir, args.pdb_file, args.output, thresholds)


if __name__ == "__main__":
    main()
