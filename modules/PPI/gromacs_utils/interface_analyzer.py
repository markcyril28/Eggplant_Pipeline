#!/usr/bin/env python3
"""
Interface Analyzer Module for GROMACS PPI Analysis

Analyzes protein-protein interfaces: BSA, H-bonds, contacts, residue pairs.
Generates comprehensive reports with visualizations.
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import numpy as np

from .xvg_parser import parse_xvg, parse_xvg_last


# Metric explanations for detailed reporting
METRIC_EXPLANATIONS = {
    'bsa': {
        'name': 'Buried Surface Area (BSA)',
        'unit': 'nm²',
        'better': 'higher',
        'explanation': '''Buried Surface Area measures the surface area that becomes buried (inaccessible 
to solvent) when two proteins bind. It is calculated as:
  BSA = (SASA_chainA + SASA_chainB - SASA_complex) / 2

Reference values for protein-protein interactions:
  • Strong binding:    > 8 nm² (> 800 Å²)
  • Moderate binding:  4-8 nm² (400-800 Å²)
  • Weak binding:      < 4 nm² (< 400 Å²)

Larger BSA generally correlates with higher binding affinity, as more atoms 
participate in the interaction and more desolvation occurs.''',
    },
    'hbonds': {
        'name': 'Interface Hydrogen Bonds',
        'unit': 'count',
        'better': 'higher',
        'explanation': '''Hydrogen bonds at the interface provide specificity and contribute to binding 
affinity. Each H-bond typically contributes 1-5 kJ/mol to binding energy.

Reference values:
  • Strong network:   ≥ 10 H-bonds
  • Moderate network: 5-9 H-bonds  
  • Weak network:     < 5 H-bonds

More H-bonds indicate more specific recognition between the proteins and 
typically correlate with tighter, more stable binding.''',
    },
    'contacts': {
        'name': 'Interface Contacts',
        'unit': 'count',
        'better': 'higher',
        'explanation': '''Atomic contacts between chains (typically within 0.6 nm cutoff) indicate the 
extent of the binding interface. More contacts generally indicate:
  • Larger binding interface
  • Better shape complementarity
  • More van der Waals interactions

Reference values:
  • Extensive interface: ≥ 50 contacts
  • Moderate interface:  20-49 contacts
  • Small interface:     < 20 contacts''',
    },
    'min_distance': {
        'name': 'Minimum Inter-chain Distance',
        'unit': 'nm',
        'better': 'lower',
        'explanation': '''The closest approach between any two atoms of different chains. Lower values 
indicate tighter packing at the interface:
  • Very tight:   < 0.25 nm (typical H-bond distance)
  • Normal:       0.25-0.35 nm
  • Loose:        > 0.35 nm

This metric helps identify if chains are in close contact and forming 
direct interactions vs. being separated.''',
    },
    'potential_energy': {
        'name': 'Potential Energy',
        'unit': 'kJ/mol',
        'better': 'lower',
        'explanation': '''Total potential energy of the system after energy minimization. Lower (more 
negative) values indicate a more thermodynamically favorable conformation.

This includes bonded interactions (bonds, angles, dihedrals) and non-bonded 
interactions (electrostatics, van der Waals). Comparing between structures 
indicates relative stability of the complex conformation.''',
    },
    'coulomb_sr': {
        'name': 'Coulomb Interaction (Short-range)',
        'unit': 'kJ/mol',
        'better': 'lower',
        'explanation': '''Short-range electrostatic (Coulomb) interaction energy between chains. More 
negative values indicate stronger favorable electrostatic complementarity:
  • Strong attraction:  < -200 kJ/mol
  • Moderate:          -200 to -50 kJ/mol
  • Weak/repulsive:    > -50 kJ/mol

Complementary charged residues (e.g., Lys-Glu salt bridges) contribute 
significantly to this energy component.''',
    },
    'lj_sr': {
        'name': 'Lennard-Jones Interaction (Short-range)',
        'unit': 'kJ/mol',
        'better': 'lower',
        'explanation': '''Short-range Lennard-Jones (van der Waals) interaction energy. More negative 
values indicate better shape complementarity at the interface:
  • Good packing:     < -100 kJ/mol
  • Moderate:         -100 to -30 kJ/mol
  • Poor packing:     > -30 kJ/mol

This reflects how well the molecular surfaces "fit" together - hydrophobic 
core packing and aromatic stacking contribute to this term.''',
    },
}


@dataclass
class InterfaceMetrics:
    """Container for interface analysis metrics."""
    bsa: float = 0.0              # Buried Surface Area (nm²)
    hbonds: int = 0               # Number of hydrogen bonds
    contacts: int = 0             # Number of contacts
    min_distance: float = 0.0     # Minimum inter-chain distance (nm)
    potential_energy: float = 0.0 # Potential energy (kJ/mol)
    coulomb_sr: float = 0.0       # Short-range Coulomb (kJ/mol)
    lj_sr: float = 0.0            # Short-range LJ (kJ/mol)
    interface_residues_a: List[int] = None
    interface_residues_b: List[int] = None
    quality_score: int = 0
    quality_assessment: str = ""
    
    def __post_init__(self):
        if self.interface_residues_a is None:
            self.interface_residues_a = []
        if self.interface_residues_b is None:
            self.interface_residues_b = []


def calculate_bsa(sasa_complex: float, sasa_a: float, sasa_b: float) -> float:
    """
    Calculate Buried Surface Area.
    
    BSA = (SASA_A + SASA_B - SASA_complex) / 2
    
    Args:
        sasa_complex: SASA of the complex
        sasa_a: SASA of chain A alone
        sasa_b: SASA of chain B alone
        
    Returns:
        Buried surface area in nm²
    """
    return (sasa_a + sasa_b - sasa_complex) / 2


def assess_binding_quality(metrics: InterfaceMetrics, 
                          bsa_good: float = 8.0,
                          bsa_moderate: float = 4.0,
                          hbond_good: int = 10,
                          hbond_moderate: int = 5,
                          contact_good: int = 50,
                          contact_moderate: int = 20) -> Tuple[int, str, List[str]]:
    """
    Assess binding quality based on interface metrics.
    
    Args:
        metrics: InterfaceMetrics object
        bsa_good/moderate: BSA thresholds
        hbond_good/moderate: H-bond thresholds
        contact_good/moderate: Contact thresholds
        
    Returns:
        Tuple of (score, assessment_string, list_of_assessments)
    """
    score = 0
    assessments = []
    
    # BSA assessment
    if metrics.bsa > bsa_good:
        score += 2
        assessments.append("✓ Large binding interface")
    elif metrics.bsa > bsa_moderate:
        score += 1
        assessments.append("○ Moderate binding interface")
    else:
        assessments.append("✗ Small binding interface")
    
    # H-bond assessment
    if metrics.hbonds >= hbond_good:
        score += 2
        assessments.append("✓ Strong H-bond network")
    elif metrics.hbonds >= hbond_moderate:
        score += 1
        assessments.append("○ Moderate H-bond network")
    else:
        assessments.append("✗ Few H-bonds")
    
    # Contact assessment
    if metrics.contacts >= contact_good:
        score += 2
        assessments.append("✓ Extensive contacts")
    elif metrics.contacts >= contact_moderate:
        score += 1
        assessments.append("○ Moderate contacts")
    else:
        assessments.append("✗ Limited contacts")
    
    # Overall assessment
    if score >= 5:
        assessment = "STRONG binding predicted"
    elif score >= 3:
        assessment = "MODERATE binding predicted"
    else:
        assessment = "WEAK binding predicted"
    
    return score, assessment, assessments


def extract_interface_metrics(workdir: Path) -> InterfaceMetrics:
    """
    Extract interface metrics from GROMACS output files.
    
    Args:
        workdir: Working directory containing GROMACS outputs
        
    Returns:
        InterfaceMetrics object
    """
    metrics = InterfaceMetrics()
    
    # Parse SASA files for BSA
    sasa_complex = parse_xvg_last(workdir / 'sasa_complex.xvg')
    sasa_a = parse_xvg_last(workdir / 'sasa_chainA.xvg')
    sasa_b = parse_xvg_last(workdir / 'sasa_chainB.xvg')
    
    if all([sasa_complex, sasa_a, sasa_b]):
        complex_val = sasa_complex[0] if len(sasa_complex) == 1 else sasa_complex[1] if len(sasa_complex) > 1 else 0
        a_val = sasa_a[0] if len(sasa_a) == 1 else sasa_a[1] if len(sasa_a) > 1 else 0
        b_val = sasa_b[0] if len(sasa_b) == 1 else sasa_b[1] if len(sasa_b) > 1 else 0
        metrics.bsa = calculate_bsa(complex_val, a_val, b_val)
    
    # Parse H-bonds
    hbond_data = parse_xvg_last(workdir / 'hbonds.xvg')
    if hbond_data:
        metrics.hbonds = int(hbond_data[0] if len(hbond_data) == 1 else hbond_data[1] if len(hbond_data) > 1 else 0)
    
    # Parse contacts
    contact_data = parse_xvg_last(workdir / 'numcont.xvg')
    if contact_data:
        metrics.contacts = int(contact_data[0] if len(contact_data) == 1 else contact_data[1] if len(contact_data) > 1 else 0)
    
    # Parse minimum distance
    mindist_data = parse_xvg_last(workdir / 'mindist.xvg')
    if mindist_data:
        metrics.min_distance = mindist_data[0] if len(mindist_data) == 1 else mindist_data[1] if len(mindist_data) > 1 else 0
    
    # Parse energies
    energy_data = parse_xvg_last(workdir / 'energies.xvg')
    if energy_data and len(energy_data) >= 1:
        metrics.potential_energy = energy_data[0]
        if len(energy_data) >= 2:
            metrics.coulomb_sr = energy_data[1]
        if len(energy_data) >= 3:
            metrics.lj_sr = energy_data[2]
    
    # Assess quality
    score, assessment, _ = assess_binding_quality(metrics)
    metrics.quality_score = score
    metrics.quality_assessment = assessment
    
    return metrics


def generate_interface_report(metrics: InterfaceMetrics, 
                             output_file: Path,
                             structure_name: str = "Unknown") -> None:
    """
    Generate a comprehensive text report of interface analysis with detailed explanations.
    
    Args:
        metrics: InterfaceMetrics object
        output_file: Path to output file
        structure_name: Name of the structure
    """
    score, assessment, assessments = assess_binding_quality(metrics)
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" PROTEIN-PROTEIN INTERFACE ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"\n  Structure: {structure_name}\n")
        f.write(f"  Report generated by GROMACS Interface Analyzer\n")
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Summary box
        f.write("┌" + "─" * 78 + "┐\n")
        f.write("│  QUICK SUMMARY" + " " * 63 + "│\n")
        f.write("├" + "─" * 78 + "┤\n")
        f.write(f"│  Binding Assessment: {assessment:<55} │\n")
        f.write(f"│  Quality Score: {score}/6" + " " * 58 + "│\n")
        f.write(f"│  Buried Surface Area: {metrics.bsa:.2f} nm²" + " " * (52 - len(f"{metrics.bsa:.2f}")) + "│\n")
        f.write(f"│  Interface H-bonds: {metrics.hbonds}" + " " * (56 - len(str(metrics.hbonds))) + "│\n")
        f.write(f"│  Interface Contacts: {metrics.contacts}" + " " * (55 - len(str(metrics.contacts))) + "│\n")
        f.write("└" + "─" * 78 + "┘\n\n")
        
        # Detailed metrics with explanations
        f.write("-" * 80 + "\n")
        f.write(" DETAILED METRIC ANALYSIS\n")
        f.write("-" * 80 + "\n\n")
        
        # BSA
        bsa_info = METRIC_EXPLANATIONS['bsa']
        f.write(f"▸ {bsa_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.bsa:.3f} {bsa_info['unit']}\n\n")
        for line in bsa_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # H-bonds
        hb_info = METRIC_EXPLANATIONS['hbonds']
        f.write(f"▸ {hb_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.hbonds} {hb_info['unit']}\n\n")
        for line in hb_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # Contacts
        ct_info = METRIC_EXPLANATIONS['contacts']
        f.write(f"▸ {ct_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.contacts} {ct_info['unit']}\n\n")
        for line in ct_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # Minimum distance
        md_info = METRIC_EXPLANATIONS['min_distance']
        f.write(f"▸ {md_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.min_distance:.4f} {md_info['unit']}\n\n")
        for line in md_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # Energy section
        f.write("-" * 80 + "\n")
        f.write(" ENERGY ANALYSIS\n")
        f.write("-" * 80 + "\n\n")
        
        # Potential energy
        pe_info = METRIC_EXPLANATIONS['potential_energy']
        f.write(f"▸ {pe_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.potential_energy:.2f} {pe_info['unit']}\n\n")
        for line in pe_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # Coulomb
        coul_info = METRIC_EXPLANATIONS['coulomb_sr']
        f.write(f"▸ {coul_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.coulomb_sr:.2f} {coul_info['unit']}\n\n")
        for line in coul_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # LJ
        lj_info = METRIC_EXPLANATIONS['lj_sr']
        f.write(f"▸ {lj_info['name']}\n")
        f.write(f"  {'─' * 70}\n")
        f.write(f"  Value: {metrics.lj_sr:.2f} {lj_info['unit']}\n\n")
        for line in lj_info['explanation'].split('\n'):
            f.write(f"  {line}\n")
        f.write("\n")
        
        # Quality assessment section
        f.write("-" * 80 + "\n")
        f.write(" BINDING QUALITY ASSESSMENT\n")
        f.write("-" * 80 + "\n\n")
        
        for a in assessments:
            f.write(f"  {a}\n")
        
        f.write(f"\n  Overall Assessment: {assessment}\n")
        f.write(f"  Quality Score: {score}/6\n\n")
        
        # Scoring explanation
        f.write("  Scoring methodology:\n")
        f.write("  • BSA > 8 nm²: +2 points, BSA > 4 nm²: +1 point\n")
        f.write("  • H-bonds ≥ 10: +2 points, H-bonds ≥ 5: +1 point\n")
        f.write("  • Contacts ≥ 50: +2 points, Contacts ≥ 20: +1 point\n")
        f.write("  • Score 5-6: STRONG binding, 3-4: MODERATE, 0-2: WEAK\n\n")
        
        # Methodology section
        f.write("=" * 80 + "\n")
        f.write(" METHODOLOGY & INTERPRETATION\n")
        f.write("=" * 80 + "\n")
        f.write("""
  This analysis evaluates the protein-protein interface using several 
  complementary metrics. Each provides different insights:

  STRUCTURAL METRICS:
  • BSA quantifies the extent of the binding interface
  • H-bonds indicate specificity of the interaction
  • Contacts reflect overall interface size and packing

  ENERGETIC METRICS:
  • Coulomb energy shows electrostatic complementarity
  • LJ energy indicates shape/hydrophobic complementarity
  • Lower (more negative) interaction energies favor binding

  CAVEATS:
  • These are based on energy-minimized structures
  • For binding affinity prediction, consider:
    - MD simulations to assess stability over time
    - MM-PBSA/MM-GBSA free energy calculations
    - Experimental validation (SPR, ITC, etc.)
  • Values should be compared between alternative structures
    rather than used as absolute predictors

""")
        f.write("=" * 80 + "\n")


def generate_interface_plots(metrics: InterfaceMetrics,
                            output_dir: Path,
                            structure_name: str = "Unknown") -> List[str]:
    """
    Generate visualization plots for interface analysis.
    
    Args:
        metrics: InterfaceMetrics object
        output_dir: Directory to save plots
        structure_name: Name of the structure
        
    Returns:
        List of generated plot paths
    """
    generated_plots = []
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available, skipping plots")
        return generated_plots
    
    short_name = structure_name[:30] + "..." if len(structure_name) > 30 else structure_name
    
    # 1. Interface metrics bar chart
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Interface Analysis: {short_name}', fontsize=14, fontweight='bold')
    
    # BSA with reference lines
    ax = axes[0, 0]
    bar = ax.bar(['BSA'], [metrics.bsa], color='steelblue', width=0.5)
    ax.axhline(y=8, color='green', linestyle='--', alpha=0.7, label='Strong (>8)')
    ax.axhline(y=4, color='orange', linestyle='--', alpha=0.7, label='Moderate (>4)')
    ax.set_ylabel('nm²')
    ax.set_title('Buried Surface Area')
    ax.legend(loc='upper right', fontsize=8)
    ax.bar_label(bar, fmt='%.2f', fontsize=10)
    
    # H-bonds with reference lines
    ax = axes[0, 1]
    bar = ax.bar(['H-bonds'], [metrics.hbonds], color='coral', width=0.5)
    ax.axhline(y=10, color='green', linestyle='--', alpha=0.7, label='Strong (≥10)')
    ax.axhline(y=5, color='orange', linestyle='--', alpha=0.7, label='Moderate (≥5)')
    ax.set_ylabel('Count')
    ax.set_title('Interface Hydrogen Bonds')
    ax.legend(loc='upper right', fontsize=8)
    ax.bar_label(bar, fontsize=10)
    
    # Contacts with reference lines
    ax = axes[1, 0]
    bar = ax.bar(['Contacts'], [metrics.contacts], color='mediumseagreen', width=0.5)
    ax.axhline(y=50, color='green', linestyle='--', alpha=0.7, label='Extensive (≥50)')
    ax.axhline(y=20, color='orange', linestyle='--', alpha=0.7, label='Moderate (≥20)')
    ax.set_ylabel('Count')
    ax.set_title('Interface Contacts')
    ax.legend(loc='upper right', fontsize=8)
    ax.bar_label(bar, fontsize=10)
    
    # Energy components
    ax = axes[1, 1]
    energies = ['Potential', 'Coulomb', 'LJ']
    values = [metrics.potential_energy/1000, metrics.coulomb_sr, metrics.lj_sr]
    colors = ['steelblue' if v < 0 else 'coral' for v in [metrics.coulomb_sr, metrics.lj_sr, metrics.lj_sr]]
    colors[0] = 'gray'
    bars = ax.bar(energies, values, color=['gray', 'steelblue', 'mediumseagreen'], width=0.6)
    ax.set_ylabel('kJ/mol (Potential: ×10³)')
    ax.set_title('Energy Components')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02*abs(bar.get_height()),
                f'{val:.1f}', ha='center', fontsize=9)
    
    plt.tight_layout()
    metrics_plot = plots_dir / "interface_metrics.png"
    plt.savefig(metrics_plot, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(str(metrics_plot))
    
    # 2. Quality score gauge-style plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    score = metrics.quality_score
    max_score = 6
    
    # Create stacked horizontal bar
    categories = ['BSA', 'H-bonds', 'Contacts']
    
    # Determine individual scores
    bsa_score = 2 if metrics.bsa > 8 else (1 if metrics.bsa > 4 else 0)
    hb_score = 2 if metrics.hbonds >= 10 else (1 if metrics.hbonds >= 5 else 0)
    ct_score = 2 if metrics.contacts >= 50 else (1 if metrics.contacts >= 20 else 0)
    
    scores = [bsa_score, hb_score, ct_score]
    colors = ['steelblue', 'coral', 'mediumseagreen']
    
    for i, (cat, s, c) in enumerate(zip(categories, scores, colors)):
        ax.barh(cat, s, color=c, alpha=0.8)
        ax.barh(cat, 2 - s, left=s, color='lightgray', alpha=0.3)
        ax.text(s/2, i, f'{s}/2', ha='center', va='center', fontweight='bold', color='white' if s > 0 else 'gray')
    
    ax.set_xlim(0, 2)
    ax.set_xlabel('Score (max 2 per category)')
    ax.set_title(f'Quality Score Breakdown: {score}/{max_score} ({metrics.quality_assessment})', fontsize=12)
    
    plt.tight_layout()
    quality_plot = plots_dir / "quality_score.png"
    plt.savefig(quality_plot, dpi=150)
    plt.close()
    generated_plots.append(str(quality_plot))
    
    return generated_plots


def save_metrics_json(metrics: InterfaceMetrics, 
                      output_file: Path,
                      structure_name: str = "Unknown") -> None:
    """
    Save interface metrics as JSON.
    
    Args:
        metrics: InterfaceMetrics object
        output_file: Path to output file
        structure_name: Name of the structure
    """
    data = {
        'structure': structure_name,
        'interface': {
            'buried_surface_area_nm2': metrics.bsa,
            'hydrogen_bonds': metrics.hbonds,
            'contacts': metrics.contacts,
            'min_distance_nm': metrics.min_distance,
            'interface_residues_chainA': metrics.interface_residues_a,
            'interface_residues_chainB': metrics.interface_residues_b,
        },
        'energy': {
            'potential_kJ_mol': metrics.potential_energy,
            'coulomb_sr_kJ_mol': metrics.coulomb_sr,
            'lj_sr_kJ_mol': metrics.lj_sr,
        },
        'quality': {
            'score': metrics.quality_score,
            'assessment': metrics.quality_assessment,
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze PPI interface")
    parser.add_argument("--workdir", required=True, help="Working directory with GROMACS outputs")
    parser.add_argument("--name", default="Unknown", help="Structure name")
    parser.add_argument("--output", help="Output directory (default: workdir)")
    args = parser.parse_args()
    
    workdir = Path(args.workdir)
    output_dir = Path(args.output) if args.output else workdir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    metrics = extract_interface_metrics(workdir)
    generate_interface_report(metrics, output_dir / "interface_analysis.txt", args.name)
    save_metrics_json(metrics, output_dir / "interface_metrics.json", args.name)
    
    # Generate plots
    print("\n  Generating interface plots...")
    plots = generate_interface_plots(metrics, output_dir, args.name)
    if plots:
        print(f"  Generated {len(plots)} plots in {output_dir / 'plots'}/")
    
    print(f"\n  Interface analysis complete:")
    print(f"    BSA: {metrics.bsa:.2f} nm²")
    print(f"    H-bonds: {metrics.hbonds}")
    print(f"    Contacts: {metrics.contacts}")
    print(f"    Quality: {metrics.quality_assessment}")
    print(f"\n  Report saved to: {output_dir / 'interface_analysis.txt'}")
