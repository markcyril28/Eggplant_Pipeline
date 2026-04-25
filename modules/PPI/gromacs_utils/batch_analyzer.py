#!/usr/bin/env python3
"""
Batch Analyzer Module for GROMACS PPI Analysis

Analyzes multiple structures and generates comprehensive comparison reports
with detailed explanations and visualizations.
"""

import os
import json
import csv
import glob
from pathlib import Path
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field

from .xvg_parser import parse_xvg_last


# Metric explanations for detailed reporting
BATCH_METRIC_EXPLANATIONS = {
    'potential': {
        'name': 'Potential Energy',
        'unit': 'kJ/mol',
        'better': 'lower',
        'explanation': '''Total potential energy of the system after energy minimization.
Lower (more negative) values indicate a more thermodynamically stable conformation.
This is the primary metric for ranking structures - the structure with the lowest
potential energy represents the most stable predicted complex.''',
    },
    'gyration': {
        'name': 'Radius of Gyration',
        'unit': 'nm',
        'better': 'context',
        'explanation': '''Measure of overall protein complex compactness.
Smaller values indicate a more compact structure. For comparing similar complexes,
consistent Rg values suggest structural similarity. Large variations may indicate
different binding modes or structural instability.''',
    },
    'sasa': {
        'name': 'Solvent Accessible Surface Area',
        'unit': 'nm²',
        'better': 'context',
        'explanation': '''Total surface area of the complex accessible to solvent.
This is an informational metric - the buried surface area (BSA) at the interface
is more relevant for binding assessment. SASA depends on overall protein size.''',
    },
    'hbonds': {
        'name': 'Interface Hydrogen Bonds',
        'unit': 'count',
        'better': 'higher',
        'explanation': '''Number of hydrogen bonds between the protein chains.
More H-bonds indicate stronger, more specific protein-protein interactions.
Reference: ≥10 H-bonds suggests strong binding specificity.''',
    },
    'contacts': {
        'name': 'Interface Contacts',
        'unit': 'count',
        'better': 'higher',
        'explanation': '''Number of atomic contacts between chains (within cutoff).
More contacts indicate a larger, more extensive binding interface.
Reference: ≥50 contacts suggests a substantial interface.''',
    },
    'min_distance': {
        'name': 'Minimum Inter-chain Distance',
        'unit': 'nm',
        'better': 'lower',
        'explanation': '''Closest approach between atoms of different chains.
Lower values indicate tighter packing at the interface.
Values around 0.25-0.35 nm are typical for well-packed interfaces.''',
    },
}


@dataclass
class StructureMetrics:
    """Metrics for a single structure."""
    name: str = ""
    potential: Optional[float] = None
    gyration: Optional[float] = None
    sasa: Optional[float] = None
    hbonds: Optional[int] = None
    contacts: Optional[int] = None
    min_distance: Optional[float] = None
    bsa: Optional[float] = None
    status: str = "OK"


def collect_structure_metrics(struct_dir: Path) -> StructureMetrics:
    """
    Collect metrics from a structure directory.
    
    Args:
        struct_dir: Path to structure directory
        
    Returns:
        StructureMetrics object
    """
    metrics = StructureMetrics(name=struct_dir.name)
    
    # Check for existing metrics.json
    metrics_file = struct_dir / 'metrics.json'
    if metrics_file.exists():
        try:
            with open(metrics_file, 'r') as f:
                content = f.read().strip()
                if not content:
                    print(f"  Warning: Empty metrics.json in {struct_dir.name}, parsing XVG files instead")
                else:
                    data = json.loads(content)
                    metrics.potential = data.get('potential')
                    metrics.gyration = data.get('gyration')
                    metrics.sasa = data.get('sasa')
                    metrics.hbonds = data.get('hbonds')
                    metrics.contacts = data.get('contacts')
                    metrics.min_distance = data.get('min_distance')
                    metrics.bsa = data.get('bsa')
                    metrics.status = data.get('status', 'OK')
                    return metrics
        except json.JSONDecodeError as e:
            print(f"  Warning: Invalid JSON in {metrics_file}: {e}, parsing XVG files instead")
        except Exception as e:
            print(f"  Warning: Error reading {metrics_file}: {e}, parsing XVG files instead")
    
    # Parse individual XVG files
    energy = parse_xvg_last(struct_dir / 'energy.xvg')
    if energy:
        metrics.potential = energy[-1] if energy else None
    
    gyrate = parse_xvg_last(struct_dir / 'gyrate.xvg')
    if gyrate:
        metrics.gyration = gyrate[-1] if len(gyrate) > 1 else gyrate[0]
    
    sasa = parse_xvg_last(struct_dir / 'sasa.xvg')
    if sasa:
        metrics.sasa = sasa[-1] if len(sasa) > 1 else sasa[0]
    
    hbonds = parse_xvg_last(struct_dir / 'hbonds.xvg')
    if hbonds:
        metrics.hbonds = int(hbonds[-1] if len(hbonds) > 1 else hbonds[0])
    
    numcont = parse_xvg_last(struct_dir / 'numcont.xvg')
    if numcont:
        metrics.contacts = int(numcont[-1] if len(numcont) > 1 else numcont[0])
    
    mindist = parse_xvg_last(struct_dir / 'mindist.xvg')
    if mindist:
        metrics.min_distance = mindist[-1] if len(mindist) > 1 else mindist[0]
    
    # Check status
    status_file = struct_dir / 'status.txt'
    if status_file.exists():
        metrics.status = status_file.read_text().strip()
    
    return metrics


def analyze_all_structures(output_dir: Path) -> List[StructureMetrics]:
    """
    Analyze all structures in output directory.
    
    Args:
        output_dir: Directory containing structure subdirectories
        
    Returns:
        List of StructureMetrics sorted by potential energy
    """
    results = []
    
    # Find all structure directories
    for struct_dir in sorted(output_dir.iterdir()):
        if struct_dir.is_dir() and not struct_dir.name.startswith(('.', '_')):
            # Skip special directories
            if struct_dir.name in ['logs', 'plots', 'visualization', 'statistics', 'analysis']:
                continue
            
            metrics = collect_structure_metrics(struct_dir)
            if metrics.status != "ERROR":
                results.append(metrics)
    
    # Sort by potential energy (more negative = more stable)
    results.sort(key=lambda x: x.potential if x.potential is not None else float('inf'))
    
    return results


def generate_comparison_report(results: List[StructureMetrics], 
                               output_dir: Path) -> Dict[str, Path]:
    """
    Generate comprehensive comparison report with detailed explanations.
    
    Args:
        results: List of StructureMetrics
        output_dir: Output directory
        
    Returns:
        Dictionary mapping file type to path
    """
    files = {}
    
    # Detailed text report
    report_file = output_dir / 'comparison_report.txt'
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" BATCH PPI STRUCTURE COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  Total structures analyzed: {len(results)}\n")
        f.write(f"  Output directory: {output_dir}\n")
        f.write("\n")
        
        # Winner box
        if results:
            best = results[0]
            f.write("┌" + "─" * 78 + "┐\n")
            f.write("│  MOST STABLE STRUCTURE" + " " * 55 + "│\n")
            f.write("├" + "─" * 78 + "┤\n")
            best_name = best.name[:60] if len(best.name) > 60 else best.name
            padding = 76 - len(best_name)
            f.write(f"│  {best_name}" + " " * padding + "│\n")
            if best.potential is not None:
                pot_str = f"Potential Energy: {best.potential:.1f} kJ/mol"
                padding2 = 76 - len(pot_str)
                f.write(f"│  {pot_str}" + " " * padding2 + "│\n")
            f.write("└" + "─" * 78 + "┘\n\n")
        
        # Ranking table
        f.write("-" * 80 + "\n")
        f.write(" STRUCTURE RANKINGS (by Potential Energy)\n")
        f.write("-" * 80 + "\n\n")
        
        f.write(f"{'Rank':<6}{'Structure':<35}{'Potential':<14}{'Rg':<10}{'H-bonds':<10}{'Contacts':<10}\n")
        f.write("─" * 85 + "\n")
        
        for i, r in enumerate(results, 1):
            pot = f"{r.potential:.0f}" if r.potential is not None else "N/A"
            rg = f"{r.gyration:.2f}" if r.gyration is not None else "N/A"
            hb = f"{r.hbonds}" if r.hbonds is not None else "N/A"
            ct = f"{r.contacts}" if r.contacts is not None else "N/A"
            name = r.name[:33] if len(r.name) > 33 else r.name
            
            # Highlight best
            marker = "★" if i == 1 else " "
            f.write(f"{marker}{i:<5}{name:<35}{pot:<14}{rg:<10}{hb:<10}{ct:<10}\n")
        
        f.write("\n")
        
        # Energy difference analysis
        if len(results) >= 2 and results[0].potential is not None and results[-1].potential is not None:
            best_pot = results[0].potential
            worst_pot = results[-1].potential
            diff = worst_pot - best_pot
            f.write("-" * 80 + "\n")
            f.write(" ENERGY ANALYSIS\n")
            f.write("-" * 80 + "\n\n")
            f.write(f"  Best structure:   {results[0].name}\n")
            f.write(f"    Potential:      {best_pot:.1f} kJ/mol\n\n")
            f.write(f"  Worst structure:  {results[-1].name}\n")
            f.write(f"    Potential:      {worst_pot:.1f} kJ/mol\n\n")
            f.write(f"  Energy difference: {diff:.1f} kJ/mol\n")
            f.write(f"  (More negative energy = more stable)\n\n")
            
            # Top 3 comparison if available
            if len(results) >= 3:
                f.write("  Top 3 structures:\n")
                for i, r in enumerate(results[:3], 1):
                    if r.potential is not None:
                        diff_from_best = r.potential - best_pot
                        f.write(f"    {i}. {r.name[:40]}\n")
                        f.write(f"       Potential: {r.potential:.1f} kJ/mol")
                        if i > 1:
                            f.write(f" (+{diff_from_best:.1f} vs best)")
                        f.write("\n")
        
        # Detailed metric explanations
        f.write("\n" + "=" * 80 + "\n")
        f.write(" METRIC EXPLANATIONS\n")
        f.write("=" * 80 + "\n\n")
        
        for key, info in BATCH_METRIC_EXPLANATIONS.items():
            f.write(f"▸ {info['name']} ({info['unit']})\n")
            f.write(f"  {'─' * 70}\n")
            f.write(f"  Better when: {info['better']}\n\n")
            for line in info['explanation'].split('\n'):
                f.write(f"  {line}\n")
            f.write("\n")
        
        # Methodology section
        f.write("=" * 80 + "\n")
        f.write(" METHODOLOGY & INTERPRETATION\n")
        f.write("=" * 80 + "\n")
        f.write("""
  This batch analysis compares protein-protein interaction (PPI) stability 
  across multiple AlphaFold-predicted structures. The ranking is based on:

  PRIMARY CRITERION: Potential Energy
  - Lower (more negative) potential energy indicates more thermodynamic stability
  - The structure with the lowest energy represents the most favorable complex

  SECONDARY CRITERIA: Interface Quality
  - H-bonds: More hydrogen bonds suggest stronger, more specific binding
  - Contacts: More contacts indicate a larger binding interface
  - Rg: Similar values across structures suggest similar overall conformations

  IMPORTANT CAVEATS:
  1. Energy minimization provides a static snapshot - MD simulations give 
     better insight into dynamic stability
  2. Large energy differences (>100 kJ/mol) are significant; small differences
     may be within noise
  3. Consider running full MD simulations on top candidates
  4. Experimental validation is essential for confirming predictions

  RECOMMENDATIONS:
  - Focus on the top 3-5 structures for further analysis
  - Run full MD simulations (Script 2) on promising candidates
  - Perform interface analysis (Script 3) for detailed binding characterization
  - Consider MutateX analysis for stability mutations

""")
        f.write("=" * 80 + "\n")
    
    files['report'] = report_file
    
    # CSV data
    csv_file = output_dir / 'comparison_data.csv'
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank', 'structure', 'potential_kJ_mol', 'gyration_nm', 
                        'sasa_nm2', 'hbonds', 'contacts', 'min_distance_nm', 'status'])
        for i, r in enumerate(results, 1):
            writer.writerow([
                i, r.name,
                r.potential if r.potential is not None else '',
                r.gyration if r.gyration is not None else '',
                r.sasa if r.sasa is not None else '',
                r.hbonds if r.hbonds is not None else '',
                r.contacts if r.contacts is not None else '',
                r.min_distance if r.min_distance is not None else '',
                r.status
            ])
    files['csv'] = csv_file
    
    # JSON data
    json_file = output_dir / 'comparison_data.json'
    with open(json_file, 'w') as f:
        json.dump({
            'structures': [
                {
                    'rank': i,
                    'name': r.name,
                    'potential_kJ_mol': r.potential,
                    'gyration_nm': r.gyration,
                    'sasa_nm2': r.sasa,
                    'hbonds': r.hbonds,
                    'contacts': r.contacts,
                    'min_distance_nm': r.min_distance,
                    'status': r.status
                }
                for i, r in enumerate(results, 1)
            ],
            'best_structure': results[0].name if results else None,
            'total_analyzed': len(results)
        }, f, indent=2)
    files['json'] = json_file
    
    return files


def generate_plotting_scripts(output_dir: Path) -> Dict[str, Path]:
    """
    Generate plotting scripts for comparison data.
    
    Args:
        output_dir: Output directory
        
    Returns:
        Dictionary mapping script type to path
    """
    files = {}
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Gnuplot script
    gnuplot_script = '''# Bar chart of potential energies
set terminal pngcairo size 1200,600 enhanced font 'Arial,10'
set output 'comparison_barchart.png'

set title 'Structure Comparison - Potential Energy'
set ylabel 'Potential Energy (kJ/mol)'
set style fill solid 0.8
set boxwidth 0.8
set xtics rotate by -45

set datafile separator ','
plot '../comparison_data.csv' skip 1 using 0:3:xtic(2) with boxes lc rgb 'blue' notitle
'''
    gp_file = plots_dir / 'plot_comparison.gp'
    with open(gp_file, 'w') as f:
        f.write(gnuplot_script)
    files['gnuplot'] = gp_file
    
    # Python plotting script
    python_script = '''#!/usr/bin/env python3
"""Generate comparison plots."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Read data
df = pd.read_csv("comparison_data.csv")

# Set style
sns.set_style("whitegrid")

# 1. Bar plot of potential energies
fig, ax = plt.subplots(figsize=(12, max(6, len(df) * 0.4)))
colors = ['green' if i == 0 else 'steelblue' for i in range(len(df))]
bars = ax.barh(df['structure'], df['potential_kJ_mol'], color=colors)
ax.set_xlabel('Potential Energy (kJ/mol)')
ax.set_title('Structure Comparison - Potential Energy\\n(Lower = More Stable, Green = Best)')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('plots/comparison_barplot.png', dpi=150)
plt.close()

# 2. Scatter: Rg vs SASA (if available)
if 'gyration_nm' in df.columns and 'sasa_nm2' in df.columns:
    df_valid = df.dropna(subset=['gyration_nm', 'sasa_nm2'])
    if len(df_valid) > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df_valid['gyration_nm'], df_valid['sasa_nm2'], 
                             c=df_valid['potential_kJ_mol'], cmap='RdYlGn_r', 
                             s=100, edgecolors='black')
        for i, row in df_valid.iterrows():
            ax.annotate(row['structure'][:15], 
                        (row['gyration_nm'], row['sasa_nm2']),
                        fontsize=8, alpha=0.7)
        plt.colorbar(scatter, label='Potential Energy (kJ/mol)')
        ax.set_xlabel('Radius of Gyration (nm)')
        ax.set_ylabel('SASA (nm²)')
        ax.set_title('Structure Compactness vs Surface Area')
        plt.tight_layout()
        plt.savefig('plots/rg_vs_sasa_scatter.png', dpi=150)
        plt.close()

print("Plots saved to plots/ directory")
'''
    py_file = plots_dir / 'generate_plots.py'
    with open(py_file, 'w') as f:
        f.write(python_script)
    os.chmod(py_file, 0o755)
    files['python'] = py_file
    
    # R script
    r_script = '''# R script for PPI comparison analysis
# Run with: Rscript analyze_comparison.R

library(ggplot2)

setwd(dirname(dirname(sys.frame(1)$ofile)))

# Read data
data <- read.csv("comparison_data.csv")

# Bar plot of potential energies
p1 <- ggplot(data, aes(x=reorder(structure, potential_kJ_mol), y=potential_kJ_mol)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Structure Comparison", x="Structure", y="Potential Energy (kJ/mol)") +
  theme_minimal()

ggsave("plots/comparison_barplot_r.png", p1, width=10, height=max(6, nrow(data)*0.4))

print("Plots saved to plots/ directory")
'''
    r_file = plots_dir / 'analyze_comparison.R'
    with open(r_file, 'w') as f:
        f.write(r_script)
    files['r'] = r_file
    
    return files


def generate_batch_plots(results: List[StructureMetrics], output_dir: Path) -> List[str]:
    """
    Generate comparison plots directly using matplotlib.
    
    Args:
        results: List of StructureMetrics sorted by potential energy
        output_dir: Output directory
        
    Returns:
        List of generated plot file paths
    """
    generated_plots = []
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available, skipping direct plots")
        return generated_plots
    
    if not results:
        return generated_plots
    
    # Prepare data
    names = [r.name[:20] + "..." if len(r.name) > 20 else r.name for r in results]
    potentials = [r.potential if r.potential is not None else 0 for r in results]
    hbonds = [r.hbonds if r.hbonds is not None else 0 for r in results]
    contacts = [r.contacts if r.contacts is not None else 0 for r in results]
    
    # 1. Horizontal bar chart of potential energies (ranked)
    fig, ax = plt.subplots(figsize=(12, max(6, len(results) * 0.4)))
    
    colors = ['green' if i == 0 else 'steelblue' for i in range(len(results))]
    y_pos = np.arange(len(results))
    
    bars = ax.barh(y_pos, potentials, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.invert_yaxis()  # Best at top
    ax.set_xlabel('Potential Energy (kJ/mol)', fontsize=12)
    ax.set_title('Structure Ranking by Potential Energy\n(Lower = More Stable, Green = Best)', 
                 fontsize=14, fontweight='bold')
    
    # Add value labels
    for bar, val in zip(bars, potentials):
        ax.text(bar.get_width() + 0.01 * abs(min(potentials)), 
                bar.get_y() + bar.get_height()/2,
                f'{val:.0f}', va='center', fontsize=9)
    
    plt.tight_layout()
    ranking_plot = plots_dir / 'energy_ranking.png'
    plt.savefig(ranking_plot, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(str(ranking_plot))
    
    # 2. Multi-metric comparison (top structures)
    n_show = min(10, len(results))
    if n_show >= 2:
        fig, axes = plt.subplots(1, 3, figsize=(15, max(5, n_show * 0.5)))
        
        top_names = names[:n_show]
        top_potentials = potentials[:n_show]
        top_hbonds = hbonds[:n_show]
        top_contacts = contacts[:n_show]
        
        y_pos = np.arange(n_show)
        
        # Potential energy
        ax = axes[0]
        colors = ['green' if i == 0 else 'steelblue' for i in range(n_show)]
        ax.barh(y_pos, top_potentials, color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_names)
        ax.invert_yaxis()
        ax.set_xlabel('kJ/mol')
        ax.set_title('Potential Energy\n(lower = better)')
        
        # H-bonds
        ax = axes[1]
        colors = ['green' if h == max(top_hbonds) else 'coral' for h in top_hbonds]
        ax.barh(y_pos, top_hbonds, color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_names)
        ax.invert_yaxis()
        ax.set_xlabel('count')
        ax.set_title('H-bonds\n(higher = better)')
        
        # Contacts
        ax = axes[2]
        colors = ['green' if c == max(top_contacts) else 'mediumseagreen' for c in top_contacts]
        ax.barh(y_pos, top_contacts, color=colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_names)
        ax.invert_yaxis()
        ax.set_xlabel('count')
        ax.set_title('Contacts\n(higher = better)')
        
        plt.suptitle(f'Top {n_show} Structures - Multi-metric Comparison', fontsize=14, fontweight='bold')
        plt.tight_layout()
        multi_plot = plots_dir / 'multi_metric_comparison.png'
        plt.savefig(multi_plot, dpi=150, bbox_inches='tight')
        plt.close()
        generated_plots.append(str(multi_plot))
    
    # 3. Summary stats
    if len(results) >= 3:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Show top 3 as podium-style
        positions = [1, 0, 2]  # Gold center, Silver left, Bronze right
        heights = [0.8, 1.0, 0.6]
        bar_colors = ['silver', 'gold', '#CD7F32']
        
        top3 = results[:3]
        for i, (pos, height, color, r) in enumerate(zip(positions, heights, bar_colors, top3)):
            ax.bar(pos, height, color=color, edgecolor='black', linewidth=2, width=0.8)
            name = r.name[:15] + "..." if len(r.name) > 15 else r.name
            ax.text(pos, height + 0.05, f"#{i+1}", ha='center', fontsize=14, fontweight='bold')
            ax.text(pos, height/2, name, ha='center', va='center', fontsize=9, rotation=90)
            if r.potential is not None:
                ax.text(pos, -0.1, f"{r.potential:.0f}\nkJ/mol", ha='center', fontsize=9)
        
        ax.set_xlim(-0.6, 2.6)
        ax.set_ylim(-0.3, 1.3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Top 3 Most Stable Structures', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        plt.tight_layout()
        podium_plot = plots_dir / 'top3_podium.png'
        plt.savefig(podium_plot, dpi=150)
        plt.close()
        generated_plots.append(str(podium_plot))
    
    # 4. Radar/Spider plot for multi-metric comparison
    if len(results) >= 2:
        # Collect metrics for radar plot
        radar_metrics = ['potential', 'gyration', 'sasa', 'hbonds', 'contacts', 'min_distance']
        available_metrics = []
        
        for metric in radar_metrics:
            values = [getattr(r, metric, None) for r in results]
            if any(v is not None for v in values):
                available_metrics.append(metric)
        
        if len(available_metrics) >= 3:
            # Prepare normalized data
            normalized_data = []
            labels = []
            
            for metric in available_metrics:
                info = BATCH_METRIC_EXPLANATIONS.get(metric, {'name': metric, 'better': 'context'})
                labels.append(info['name'])
                
                values = [getattr(r, metric, None) for r in results]
                # Replace None with median
                valid_values = [v for v in values if v is not None]
                if not valid_values:
                    continue
                    
                median_val = sorted(valid_values)[len(valid_values)//2]
                values = [v if v is not None else median_val for v in values]
                
                min_val = min(values)
                max_val = max(values)
                range_val = max_val - min_val if max_val != min_val else 1
                
                better = info.get('better', 'context')
                
                # Normalize to 0-1, where 1 is "better"
                norm_values = []
                for v in values:
                    normalized = (v - min_val) / range_val
                    # Invert if lower is better
                    if better == 'lower':
                        normalized = 1 - normalized
                    norm_values.append(normalized)
                
                normalized_data.append(norm_values)
            
            if normalized_data:
                # Transpose for per-structure data
                per_structure_data = list(zip(*normalized_data))
                
                # Generate colors
                colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
                
                # Limit to top structures for readability
                n_show = min(10, len(results))
                
                # Radar Plot
                fig, ax = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection='polar'))
                
                angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
                angles += angles[:1]  # Complete the loop
                
                for i in range(n_show):
                    data = per_structure_data[i]
                    data_closed = list(data) + [data[0]]  # Close the polygon
                    name = names[i]
                    ax.plot(angles, data_closed, 'o-', linewidth=2, label=f"S{i+1}: {name}", color=colors[i])
                    ax.fill(angles, data_closed, alpha=0.1, color=colors[i])
                
                ax.set_xticks(angles[:-1])
                ax.set_xticklabels(labels, size=9)
                ax.set_title(f'Batch Stability Comparison - Radar Plot\n({n_show} structures, outer = better)', 
                             pad=20, fontsize=14, fontweight='bold')
                ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.0), fontsize=8)
                ax.set_ylim(0, 1)
                
                plt.tight_layout()
                radar_plot = plots_dir / 'batch_radar.png'
                plt.savefig(radar_plot, dpi=150, bbox_inches='tight')
                plt.close()
                generated_plots.append(str(radar_plot))
                
                # 5. Combined Ranking Stacked Bar Chart
                fig, ax = plt.subplots(figsize=(14, max(8, n_show * 0.5)))
                
                # Calculate combined score for each structure
                scores = [sum(data) for data in per_structure_data[:n_show]]
                
                # Sort by score
                sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
                sorted_names = [f"S{sorted_indices[i]+1}: {names[sorted_indices[i]]}" for i in range(len(sorted_indices))]
                sorted_scores = [scores[i] for i in sorted_indices]
                
                # Create stacked bars
                bar_width = 0.6
                y_pos = np.arange(len(sorted_names))
                
                bottom = np.zeros(len(sorted_names))
                metric_colors = plt.cm.Set2(np.linspace(0, 1, len(labels)))
                
                for j, (label, mc) in enumerate(zip(labels, metric_colors)):
                    metric_values = [per_structure_data[sorted_indices[i]][j] for i in range(len(sorted_names))]
                    ax.barh(y_pos, metric_values, bar_width, left=bottom, label=label, color=mc, edgecolor='white')
                    bottom += metric_values
                
                ax.set_yticks(y_pos)
                ax.set_yticklabels(sorted_names)
                ax.invert_yaxis()  # Best at top
                ax.set_xlabel('Combined Normalized Score (higher = better)', fontsize=12)
                ax.set_title('Batch Structure Ranking\n(Stacked contribution from each metric)', 
                             fontsize=14, fontweight='bold')
                ax.legend(loc='lower right', fontsize=9)
                
                # Add total score labels
                for i, score in enumerate(sorted_scores):
                    ax.text(score + 0.05, i, f'{score:.2f}', va='center', fontsize=10, fontweight='bold')
                
                # Mark the winner
                ax.annotate('★ BEST', xy=(sorted_scores[0], 0), fontsize=12, fontweight='bold',
                            color='green', va='center', ha='left',
                            xytext=(sorted_scores[0] + 0.3, 0))
                
                plt.tight_layout()
                ranking_plot = plots_dir / 'batch_combined_ranking.png'
                plt.savefig(ranking_plot, dpi=150, bbox_inches='tight')
                plt.close()
                generated_plots.append(str(ranking_plot))
    
    return generated_plots


def run_batch_analysis(output_dir: Path) -> Dict[str, Any]:
    """
    Run complete batch analysis with detailed reports and plots.
    
    Args:
        output_dir: Directory containing analyzed structures
        
    Returns:
        Dictionary with results and file paths
    """
    output_dir = Path(output_dir)
    
    # Collect results
    results = analyze_all_structures(output_dir)
    
    if not results:
        print("No valid structures found for comparison")
        return {'error': 'No valid structures'}
    
    # Generate reports
    report_files = generate_comparison_report(results, output_dir)
    
    # Generate plotting scripts (for manual use)
    plot_files = generate_plotting_scripts(output_dir)
    
    # Generate plots directly
    print("\n  Generating comparison plots...")
    generated_plots = generate_batch_plots(results, output_dir)
    if generated_plots:
        print(f"  Generated {len(generated_plots)} plots in {output_dir / 'plots'}/")
    
    print(f"\n  Batch comparison complete: {len(results)} structures analyzed")
    print(f"  Report: {report_files['report']}")
    print(f"  CSV: {report_files['csv']}")
    print(f"  JSON: {report_files['json']}")
    
    if results:
        print(f"\n  ★ Best structure: {results[0].name}")
        if results[0].potential is not None:
            print(f"    Potential: {results[0].potential:.1f} kJ/mol")
    
    return {
        'structures': len(results),
        'best': results[0].name if results else None,
        'files': {**report_files, **plot_files},
        'plots': generated_plots
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Batch structure comparison")
    parser.add_argument("--output-dir", required=True, help="Output directory with analyzed structures")
    args = parser.parse_args()
    
    result = run_batch_analysis(Path(args.output_dir))
