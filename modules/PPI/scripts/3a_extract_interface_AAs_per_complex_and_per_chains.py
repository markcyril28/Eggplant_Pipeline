#!/usr/bin/env python3
"""
Interface Amino Acid Extraction Script

This script extracts amino acids from the interface regions of protein complexes
and generates various output formats including:
- Raw FASTA files (continuous sequence)
- FASTA with dashes (showing gaps/positions)
- Detailed annotation files with positions and amino acid names
- Separate files for general and critical interface residues

Author: PPI Analysis Pipeline
Date: 2026-02-04
"""

import os
import sys
import re
import json
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional

# ============================================================================
# AMINO ACID CONVERSION TABLES
# ============================================================================

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Non-standard amino acids
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XLE': 'J',
    'XAA': 'X', 'UNK': 'X'
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

ONE_TO_NAME = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
    'C': 'Cysteine', 'Q': 'Glutamine', 'E': 'Glutamic Acid', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    'U': 'Selenocysteine', 'O': 'Pyrrolysine', 'B': 'Asparagine/Aspartic Acid',
    'Z': 'Glutamine/Glutamic Acid', 'J': 'Leucine/Isoleucine', 'X': 'Unknown'
}

# ============================================================================
# DATA CLASSES
# ============================================================================

class InterfaceResidue:
    """Represents a single interface residue."""
    def __init__(self, resname: str, resnum: int, chain: str, 
                 partner_resname: str = None, partner_resnum: int = None,
                 distance: float = None):
        self.resname = resname.upper()
        self.resnum = resnum
        self.chain = chain
        self.partner_resname = partner_resname
        self.partner_resnum = partner_resnum
        self.distance = distance
        self.one_letter = THREE_TO_ONE.get(self.resname, 'X')
    
    def __repr__(self):
        return f"{self.resname}{self.resnum}({self.chain})"
    
    def __hash__(self):
        return hash((self.resname, self.resnum, self.chain))
    
    def __eq__(self, other):
        return (self.resname, self.resnum, self.chain) == (other.resname, other.resnum, other.chain)


class InterfaceData:
    """Container for all interface data of a complex."""
    def __init__(self, complex_name: str):
        self.complex_name = complex_name
        self.chain_a_residues: Dict[int, InterfaceResidue] = {}
        self.chain_b_residues: Dict[int, InterfaceResidue] = {}
        self.residue_pairs: List[Tuple[InterfaceResidue, InterfaceResidue, float]] = []
        self.min_distance = float('inf')
        self.max_distance = 0.0
        
    def add_pair(self, res_a: InterfaceResidue, res_b: InterfaceResidue, distance: float):
        self.residue_pairs.append((res_a, res_b, distance))
        self.chain_a_residues[res_a.resnum] = res_a
        self.chain_b_residues[res_b.resnum] = res_b
        self.min_distance = min(self.min_distance, distance)
        self.max_distance = max(self.max_distance, distance)
    
    def get_critical_residues(self, threshold_percentile: float = 25.0) -> Tuple[Dict, Dict]:
        """Get critical residues based on distance threshold (closest contacts)."""
        if not self.residue_pairs:
            return {}, {}
        
        # Calculate distance threshold (bottom percentile = closest)
        distances = sorted([p[2] for p in self.residue_pairs])
        threshold_idx = max(1, int(len(distances) * threshold_percentile / 100))
        threshold_distance = distances[threshold_idx - 1]
        
        critical_a = {}
        critical_b = {}
        
        for res_a, res_b, dist in self.residue_pairs:
            if dist <= threshold_distance:
                critical_a[res_a.resnum] = res_a
                critical_b[res_b.resnum] = res_b
        
        return critical_a, critical_b


# ============================================================================
# PARSING FUNCTIONS
# ============================================================================

def parse_residue_string(res_str: str) -> Tuple[str, int]:
    """Parse residue string like 'ALA63' into ('ALA', 63)."""
    res_str = res_str.strip()
    match = re.match(r'([A-Z]{3})(\d+)', res_str)
    if match:
        return match.group(1), int(match.group(2))
    return None, None


def parse_interface_residues_file(filepath: Path) -> InterfaceData:
    """Parse interface_residues.txt file and extract residue pairs."""
    complex_name = filepath.parent.parent.name
    data = InterfaceData(complex_name)
    
    if not filepath.exists():
        print(f"  Warning: {filepath} not found")
        return data
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    in_data_section = False
    for line in lines:
        line = line.strip()
        
        # Skip header lines
        if line.startswith('Interface Residue Pairs') or line.startswith('=') or line.startswith('-'):
            continue
        if 'ChainA' in line and 'ChainB' in line:
            in_data_section = True
            continue
        if line.startswith('Total interface'):
            break
        if not line or not in_data_section:
            continue
        
        # Parse data line: ALA63         TYR357        0.382
        parts = line.split()
        if len(parts) >= 3:
            resname_a, resnum_a = parse_residue_string(parts[0])
            resname_b, resnum_b = parse_residue_string(parts[1])
            try:
                distance = float(parts[2])
            except ValueError:
                continue
            
            if resname_a and resname_b:
                res_a = InterfaceResidue(resname_a, resnum_a, 'A', 
                                         resname_b, resnum_b, distance)
                res_b = InterfaceResidue(resname_b, resnum_b, 'B',
                                         resname_a, resnum_a, distance)
                data.add_pair(res_a, res_b, distance)
    
    return data


def parse_protein_gro(filepath: Path) -> Dict[str, Dict[int, str]]:
    """Parse protein.gro file to get full sequence for each chain."""
    chains = {'A': {}, 'B': {}}
    
    if not filepath.exists():
        return chains
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        return chains
    
    # Skip header and count
    current_resnum = 0
    current_resname = None
    residue_count = 0
    chain_a_start = None
    chain_b_start = None
    
    for line in lines[2:-1]:  # Skip header, count, and box line
        if len(line) < 20:
            continue
        
        # GRO format: residue number (5), residue name (5), atom name (5), atom number (5)
        try:
            resnum = int(line[:5].strip())
            resname = line[5:10].strip()
            
            # Skip non-amino acid entries
            if resname not in THREE_TO_ONE:
                continue
            
            if resnum != current_resnum:
                current_resnum = resnum
                current_resname = resname
                residue_count += 1
                
                # Detect chain boundary (residue number reset or large gap)
                if chain_a_start is None:
                    chain_a_start = resnum
                    chains['A'][resnum] = resname
                elif chain_b_start is None and resnum < current_resnum - 10:
                    chain_b_start = resnum
                    chains['B'][resnum] = resname
                elif chain_b_start is not None:
                    chains['B'][resnum] = resname
                else:
                    chains['A'][resnum] = resname
        except (ValueError, IndexError):
            continue
    
    return chains


# ============================================================================
# OUTPUT GENERATION FUNCTIONS
# ============================================================================

def generate_fasta_raw(residues: Dict[int, InterfaceResidue], 
                       sequence_name: str, 
                       description: str = "") -> str:
    """Generate raw FASTA format (continuous sequence)."""
    if not residues:
        return f">{sequence_name} {description}\n"
    
    sorted_residues = sorted(residues.values(), key=lambda x: x.resnum)
    sequence = ''.join([r.one_letter for r in sorted_residues])
    
    # Format sequence in lines of 60 characters
    formatted_seq = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
    
    return f">{sequence_name} {description}\n{formatted_seq}\n"


def generate_fasta_with_dashes(residues: Dict[int, InterfaceResidue],
                                sequence_name: str,
                                description: str = "",
                                max_resnum: int = None) -> str:
    """Generate FASTA with dashes showing gaps/positions."""
    if not residues:
        return f">{sequence_name} {description}\n"
    
    if max_resnum is None:
        max_resnum = max(residues.keys()) if residues else 0
    
    sequence = []
    for i in range(1, max_resnum + 1):
        if i in residues:
            sequence.append(residues[i].one_letter)
        else:
            sequence.append('-')
    
    seq_str = ''.join(sequence)
    # Format sequence in lines of 60 characters
    formatted_seq = '\n'.join([seq_str[i:i+60] for i in range(0, len(seq_str), 60)])
    
    return f">{sequence_name} {description}\n{formatted_seq}\n"


def generate_annotation_file(residues: Dict[int, InterfaceResidue],
                              sequence_name: str,
                              description: str = "") -> str:
    """Generate detailed annotation with position, 3-letter, 1-letter, and full name."""
    if not residues:
        return f"# {sequence_name}\n# {description}\n# No residues found\n"
    
    lines = [
        f"# {sequence_name}",
        f"# {description}",
        f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# Total residues: {len(residues)}",
        "#",
        "# Position  3-Letter  1-Letter  Full Name            Partner    Distance(nm)",
        "#" + "=" * 80
    ]
    
    sorted_residues = sorted(residues.values(), key=lambda x: x.resnum)
    
    for res in sorted_residues:
        full_name = ONE_TO_NAME.get(res.one_letter, 'Unknown')
        partner_str = f"{res.partner_resname}{res.partner_resnum}" if res.partner_resname else "-"
        distance_str = f"{res.distance:.3f}" if res.distance else "-"
        
        # Format with proper column alignment
        line = f"{res.resnum:>8}  {res.resname:<8}  {res.one_letter:<8}  {full_name:<18}  {partner_str:<9}  {distance_str:>11}"
        lines.append(line)
    
    return '\n'.join(lines) + '\n'


def generate_annotation_csv(residues: Dict[int, InterfaceResidue],
                            complex_name: str,
                            chain: str,
                            residue_type: str) -> str:
    """Generate CSV format annotation file."""
    import csv
    from io import StringIO
    
    output = StringIO()
    writer = csv.writer(output)
    
    # Write header
    writer.writerow([
        'Complex', 'Chain', 'Type', 'Position', '3-Letter', '1-Letter', 
        'Full Name', 'Partner Residue', 'Partner Position', 'Distance (nm)'
    ])
    
    # Write data rows
    sorted_residues = sorted(residues.values(), key=lambda x: x.resnum)
    for res in sorted_residues:
        full_name = ONE_TO_NAME.get(res.one_letter, 'Unknown')
        partner_resname = res.partner_resname if res.partner_resname else ''
        partner_resnum = res.partner_resnum if res.partner_resnum else ''
        distance = f"{res.distance:.3f}" if res.distance else ''
        
        writer.writerow([
            complex_name, chain, residue_type, res.resnum, res.resname, 
            res.one_letter, full_name, partner_resname, partner_resnum, distance
        ])
    
    return output.getvalue()


def generate_summary_json(interface_data: InterfaceData,
                          critical_a: Dict[int, InterfaceResidue],
                          critical_b: Dict[int, InterfaceResidue]) -> dict:
    """Generate summary data as JSON-compatible dict."""
    return {
        "complex_name": interface_data.complex_name,
        "generated": datetime.now().isoformat(),
        "statistics": {
            "total_residue_pairs": len(interface_data.residue_pairs),
            "chainA_interface_residues": len(interface_data.chain_a_residues),
            "chainB_interface_residues": len(interface_data.chain_b_residues),
            "chainA_critical_residues": len(critical_a),
            "chainB_critical_residues": len(critical_b),
            "min_distance_nm": interface_data.min_distance if interface_data.min_distance != float('inf') else None,
            "max_distance_nm": interface_data.max_distance if interface_data.max_distance != 0 else None
        },
        "chainA": {
            "all_interface": {
                "residue_count": len(interface_data.chain_a_residues),
                "sequence": ''.join([interface_data.chain_a_residues[k].one_letter 
                                    for k in sorted(interface_data.chain_a_residues.keys())]),
                "positions": sorted(interface_data.chain_a_residues.keys())
            },
            "critical": {
                "residue_count": len(critical_a),
                "sequence": ''.join([critical_a[k].one_letter for k in sorted(critical_a.keys())]),
                "positions": sorted(critical_a.keys())
            }
        },
        "chainB": {
            "all_interface": {
                "residue_count": len(interface_data.chain_b_residues),
                "sequence": ''.join([interface_data.chain_b_residues[k].one_letter 
                                    for k in sorted(interface_data.chain_b_residues.keys())]),
                "positions": sorted(interface_data.chain_b_residues.keys())
            },
            "critical": {
                "residue_count": len(critical_b),
                "sequence": ''.join([critical_b[k].one_letter for k in sorted(critical_b.keys())]),
                "positions": sorted(critical_b.keys())
            }
        }
    }


# ============================================================================
# MAIN PROCESSING FUNCTIONS
# ============================================================================

def process_complex(complex_dir: Path, output_base: Path = None) -> Optional[dict]:
    """Process a single protein complex directory."""
    complex_name = complex_dir.name
    print(f"\nProcessing: {complex_name}")
    
    # Find interface_residues.txt
    interface_file = complex_dir / "analysis" / "interface_residues.txt"
    if not interface_file.exists():
        print(f"  Skipping: No interface_residues.txt found")
        return None
    
    # Parse interface data
    interface_data = parse_interface_residues_file(interface_file)
    
    if not interface_data.residue_pairs:
        print(f"  Warning: No interface residue pairs found")
        return None
    
    print(f"  Found {len(interface_data.residue_pairs)} residue pairs")
    print(f"  Chain A: {len(interface_data.chain_a_residues)} interface residues")
    print(f"  Chain B: {len(interface_data.chain_b_residues)} interface residues")
    
    # Get critical residues (closest 25% by distance)
    critical_a, critical_b = interface_data.get_critical_residues(threshold_percentile=25.0)
    print(f"  Critical Chain A: {len(critical_a)} residues")
    print(f"  Critical Chain B: {len(critical_b)} residues")
    
    # Determine output directory
    if output_base is None:
        output_dir = complex_dir / "AAs_interface_extraction"
    else:
        output_dir = output_base / complex_name
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories for organization
    fasta_dir = output_dir / "fasta"
    annotation_dir = output_dir / "annotations"
    csv_dir = output_dir / "csv"
    fasta_dir.mkdir(exist_ok=True)
    annotation_dir.mkdir(exist_ok=True)
    csv_dir.mkdir(exist_ok=True)
    
    # Get max residue numbers for proper alignment
    max_resnum_a = max(interface_data.chain_a_residues.keys()) if interface_data.chain_a_residues else 0
    max_resnum_b = max(interface_data.chain_b_residues.keys()) if interface_data.chain_b_residues else 0
    
    # ============ Generate Chain A outputs ============
    # All interface residues - raw FASTA
    with open(fasta_dir / f"{complex_name}_chainA_interface.fasta", 'w') as f:
        f.write(generate_fasta_raw(
            interface_data.chain_a_residues,
            f"{complex_name}_chainA_interface",
            f"All interface residues from Chain A ({len(interface_data.chain_a_residues)} residues)"
        ))
    
    # All interface residues - FASTA with dashes
    with open(fasta_dir / f"{complex_name}_chainA_interface_aligned.fasta", 'w') as f:
        f.write(generate_fasta_with_dashes(
            interface_data.chain_a_residues,
            f"{complex_name}_chainA_interface_aligned",
            f"Interface residues with gaps ({len(interface_data.chain_a_residues)} residues, {max_resnum_a} positions)",
            max_resnum_a
        ))
    
    # Critical residues - raw FASTA
    with open(fasta_dir / f"{complex_name}_chainA_critical.fasta", 'w') as f:
        f.write(generate_fasta_raw(
            critical_a,
            f"{complex_name}_chainA_critical",
            f"Critical interface residues from Chain A ({len(critical_a)} residues)"
        ))
    
    # Critical residues - FASTA with dashes
    with open(fasta_dir / f"{complex_name}_chainA_critical_aligned.fasta", 'w') as f:
        f.write(generate_fasta_with_dashes(
            critical_a,
            f"{complex_name}_chainA_critical_aligned",
            f"Critical residues with gaps ({len(critical_a)} residues, {max_resnum_a} positions)",
            max_resnum_a
        ))
    
    # Annotation files
    with open(annotation_dir / f"{complex_name}_chainA_interface_annotated.txt", 'w') as f:
        f.write(generate_annotation_file(
            interface_data.chain_a_residues,
            f"{complex_name} Chain A - All Interface Residues",
            f"Total: {len(interface_data.chain_a_residues)} residues"
        ))
    
    with open(annotation_dir / f"{complex_name}_chainA_critical_annotated.txt", 'w') as f:
        f.write(generate_annotation_file(
            critical_a,
            f"{complex_name} Chain A - Critical Interface Residues",
            f"Critical: {len(critical_a)} residues (closest 25% contacts)"
        ))
    
    # CSV files for Chain A
    with open(csv_dir / f"{complex_name}_chainA_interface.csv", 'w') as f:
        f.write(generate_annotation_csv(
            interface_data.chain_a_residues, complex_name, 'A', 'Interface'
        ))
    
    with open(csv_dir / f"{complex_name}_chainA_critical.csv", 'w') as f:
        f.write(generate_annotation_csv(
            critical_a, complex_name, 'A', 'Critical'
        ))
    
    # ============ Generate Chain B outputs ============
    # All interface residues - raw FASTA
    with open(fasta_dir / f"{complex_name}_chainB_interface.fasta", 'w') as f:
        f.write(generate_fasta_raw(
            interface_data.chain_b_residues,
            f"{complex_name}_chainB_interface",
            f"All interface residues from Chain B ({len(interface_data.chain_b_residues)} residues)"
        ))
    
    # All interface residues - FASTA with dashes
    with open(fasta_dir / f"{complex_name}_chainB_interface_aligned.fasta", 'w') as f:
        f.write(generate_fasta_with_dashes(
            interface_data.chain_b_residues,
            f"{complex_name}_chainB_interface_aligned",
            f"Interface residues with gaps ({len(interface_data.chain_b_residues)} residues, {max_resnum_b} positions)",
            max_resnum_b
        ))
    
    # Critical residues - raw FASTA
    with open(fasta_dir / f"{complex_name}_chainB_critical.fasta", 'w') as f:
        f.write(generate_fasta_raw(
            critical_b,
            f"{complex_name}_chainB_critical",
            f"Critical interface residues from Chain B ({len(critical_b)} residues)"
        ))
    
    # Critical residues - FASTA with dashes
    with open(fasta_dir / f"{complex_name}_chainB_critical_aligned.fasta", 'w') as f:
        f.write(generate_fasta_with_dashes(
            critical_b,
            f"{complex_name}_chainB_critical_aligned",
            f"Critical residues with gaps ({len(critical_b)} residues, {max_resnum_b} positions)",
            max_resnum_b
        ))
    
    # Annotation files
    with open(annotation_dir / f"{complex_name}_chainB_interface_annotated.txt", 'w') as f:
        f.write(generate_annotation_file(
            interface_data.chain_b_residues,
            f"{complex_name} Chain B - All Interface Residues",
            f"Total: {len(interface_data.chain_b_residues)} residues"
        ))
    
    with open(annotation_dir / f"{complex_name}_chainB_critical_annotated.txt", 'w') as f:
        f.write(generate_annotation_file(
            critical_b,
            f"{complex_name} Chain B - Critical Interface Residues",
            f"Critical: {len(critical_b)} residues (closest 25% contacts)"
        ))
    
    # CSV files for Chain B
    with open(csv_dir / f"{complex_name}_chainB_interface.csv", 'w') as f:
        f.write(generate_annotation_csv(
            interface_data.chain_b_residues, complex_name, 'B', 'Interface'
        ))
    
    with open(csv_dir / f"{complex_name}_chainB_critical.csv", 'w') as f:
        f.write(generate_annotation_csv(
            critical_b, complex_name, 'B', 'Critical'
        ))
    
    # ============ Generate combined outputs ============
    # Combined FASTA (all chains, all interface)
    combined_fasta = (
        generate_fasta_raw(
            interface_data.chain_a_residues,
            f"{complex_name}_chainA",
            "Chain A interface residues"
        ) +
        generate_fasta_raw(
            interface_data.chain_b_residues,
            f"{complex_name}_chainB",
            "Chain B interface residues"
        )
    )
    with open(fasta_dir / f"{complex_name}_all_interface.fasta", 'w') as f:
        f.write(combined_fasta)
    
    # Combined FASTA (all chains, critical only)
    combined_critical = (
        generate_fasta_raw(critical_a, f"{complex_name}_chainA_critical", "Chain A critical") +
        generate_fasta_raw(critical_b, f"{complex_name}_chainB_critical", "Chain B critical")
    )
    with open(fasta_dir / f"{complex_name}_all_critical.fasta", 'w') as f:
        f.write(combined_critical)
    
    # Summary JSON
    summary = generate_summary_json(interface_data, critical_a, critical_b)
    with open(output_dir / f"{complex_name}_interface_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Human-readable summary
    with open(output_dir / f"{complex_name}_interface_summary.txt", 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"INTERFACE AMINO ACID EXTRACTION SUMMARY\n")
        f.write(f"Complex: {complex_name}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total residue pairs: {len(interface_data.residue_pairs)}\n")
        f.write(f"Distance range: {interface_data.min_distance:.3f} - {interface_data.max_distance:.3f} nm\n\n")
        
        f.write("CHAIN A (Interface)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total interface residues: {len(interface_data.chain_a_residues)}\n")
        f.write(f"Critical residues: {len(critical_a)}\n")
        f.write(f"Positions: {sorted(interface_data.chain_a_residues.keys())}\n")
        seq_a = ''.join([interface_data.chain_a_residues[k].one_letter 
                         for k in sorted(interface_data.chain_a_residues.keys())])
        f.write(f"Sequence: {seq_a}\n\n")
        
        f.write("CHAIN B (Interface)\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total interface residues: {len(interface_data.chain_b_residues)}\n")
        f.write(f"Critical residues: {len(critical_b)}\n")
        f.write(f"Positions: {sorted(interface_data.chain_b_residues.keys())}\n")
        seq_b = ''.join([interface_data.chain_b_residues[k].one_letter 
                         for k in sorted(interface_data.chain_b_residues.keys())])
        f.write(f"Sequence: {seq_b}\n\n")
        
        f.write("OUTPUT FILES\n")
        f.write("-" * 40 + "\n")
        f.write("fasta/\n")
        f.write("  - *_interface.fasta      : Raw FASTA of all interface residues\n")
        f.write("  - *_interface_aligned.fasta : FASTA with dashes for gaps\n")
        f.write("  - *_critical.fasta       : Raw FASTA of critical residues\n")
        f.write("  - *_critical_aligned.fasta : Critical residues with gaps\n")
        f.write("  - *_all_interface.fasta  : Combined chains interface\n")
        f.write("  - *_all_critical.fasta   : Combined chains critical\n")
        f.write("\nannotations/\n")
        f.write("  - *_interface_annotated.txt : Detailed residue information\n")
        f.write("  - *_critical_annotated.txt  : Critical residue details\n")
        f.write("\ncsv/\n")
        f.write("  - *_interface.csv        : CSV format of interface residues\n")
        f.write("  - *_critical.csv         : CSV format of critical residues\n")
    
    print(f"  Output saved to: {output_dir}")
    
    return summary


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Extract interface amino acids from protein complexes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                              # Process all complexes in default directory
  %(prog)s -i /path/to/analysis         # Specify input directory
  %(prog)s -o /path/to/output           # Specify output directory
  %(prog)s -c SmelDMP01.730_SmelHAP2    # Process specific complex
  %(prog)s --critical-threshold 20      # Use 20%% for critical residue cutoff
        """
    )
    
    parser.add_argument('-i', '--input-dir', type=Path,
                        default=Path.home() / 'PPI/RESULTS/results_gromacs/3_interface_analysis',
                        help='Input directory containing complex folders')
    parser.add_argument('-o', '--output-dir', type=Path,
                        default=None,
                        help='Output directory (default: creates "AAs_interface_extraction" in each complex folder)')
    parser.add_argument('-c', '--complex', type=str,
                        default=None,
                        help='Process only specific complex by name')
    parser.add_argument('--critical-threshold', type=float,
                        default=25.0,
                        help='Percentile threshold for critical residues (default: 25.0)')
    parser.add_argument('--save-combined', action='store_true',
                        help='Also save combined summary of all complexes')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not args.input_dir.exists():
        print(f"Error: Input directory not found: {args.input_dir}")
        sys.exit(1)
    
    # Find complex directories
    if args.complex:
        complex_dirs = [args.input_dir / args.complex]
        if not complex_dirs[0].exists():
            print(f"Error: Complex not found: {args.complex}")
            sys.exit(1)
    else:
        complex_dirs = [d for d in args.input_dir.iterdir() 
                       if d.is_dir() and not d.name.startswith('.')]
    
    if not complex_dirs:
        print(f"No complex directories found in: {args.input_dir}")
        sys.exit(1)
    
    print("=" * 80)
    print("INTERFACE AMINO ACID EXTRACTION")
    print("=" * 80)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir or 'per-complex folders'}")
    print(f"Critical threshold: {args.critical_threshold}%")
    print(f"Complexes to process: {len(complex_dirs)}")
    
    # Process each complex
    all_summaries = []
    successful = 0
    failed = 0
    
    for complex_dir in sorted(complex_dirs):
        try:
            summary = process_complex(complex_dir, args.output_dir)
            if summary:
                all_summaries.append(summary)
                successful += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  Error processing {complex_dir.name}: {str(e)}")
            failed += 1
    
    # Generate combined summary if requested
    if args.save_combined and all_summaries:
        combined_output = args.output_dir if args.output_dir else args.input_dir
        combined_output.mkdir(parents=True, exist_ok=True)
        
        combined_summary = {
            "generated": datetime.now().isoformat(),
            "total_complexes": len(all_summaries),
            "complexes": all_summaries
        }
        
        with open(combined_output / "all_complexes_interface_summary.json", 'w') as f:
            json.dump(combined_summary, f, indent=2)
        
        print(f"\nCombined summary saved to: {combined_output / 'all_complexes_interface_summary.json'}")
    
    # Final summary
    print("\n" + "=" * 80)
    print("PROCESSING COMPLETE")
    print("=" * 80)
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total: {len(complex_dirs)}")


if __name__ == "__main__":
    main()
