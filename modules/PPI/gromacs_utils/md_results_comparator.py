#!/usr/bin/env python3
"""
MD Results Comparator Module for GROMACS Stability Analysis

Compares full MD simulation stability metrics between two structures
and generates comprehensive comparison reports with visualizations.
Saves comparison data in CSV format for later analysis.
"""

import os
import json
import csv
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple, Optional, Any, List


# Metric definitions: (display_name, unit, better_when, explanation)
METRIC_DESCRIPTIONS = {
    'rmsd_mean': ('RMSD (mean)', 'nm', 'lower',
        'Root Mean Square Deviation from initial structure. Lower values indicate '
        'the protein complex maintains its initial conformation better throughout simulation, '
        'suggesting structural stability.'),
    'rmsd_std': ('RMSD (std dev)', 'nm', 'lower',
        'Standard deviation of RMSD. Lower values indicate consistent stability with '
        'fewer structural fluctuations during simulation.'),
    'rmsd_final': ('RMSD (final)', 'nm', 'lower',
        'Final RMSD value at end of simulation. Lower values indicate the complex '
        'remained close to its initial structure.'),
    'rg_mean': ('Radius of Gyration (mean)', 'nm', 'similar',
        'Measure of overall protein compactness. Stable values indicate '
        'the complex maintains its overall shape during simulation.'),
    'rg_std': ('Radius of Gyration (std)', 'nm', 'lower',
        'Fluctuation in radius of gyration. Lower values indicate the complex '
        'does not expand/contract significantly during simulation.'),
    'hbonds_mean': ('H-bonds (mean)', 'count', 'higher',
        'Average number of hydrogen bonds between chains. More H-bonds indicate '
        'stronger, more specific protein-protein interactions at the interface.'),
    'hbonds_std': ('H-bonds (std)', 'count', 'lower',
        'Fluctuation in H-bond count. Lower values indicate stable, '
        'persistent hydrogen bonding throughout simulation.'),
    'mindist_mean': ('Min Distance (mean)', 'nm', 'lower',
        'Average minimum distance between chains. Lower values indicate '
        'tighter packing at the protein-protein interface.'),
    'mindist_std': ('Min Distance (std)', 'nm', 'lower',
        'Fluctuation in minimum distance. Lower values indicate consistent, '
        'stable interface contact throughout simulation.'),
    'sasa_mean': ('SASA (mean)', 'nm²', 'context',
        'Solvent Accessible Surface Area. Informational metric - the buried '
        'surface area (difference from unbound) is more relevant for binding.'),
    'sasa_std': ('SASA (std)', 'nm²', 'lower',
        'SASA fluctuation. Lower values indicate the exposed surface area '
        'remains consistent during simulation.'),
    'coul_mean': ('Coulomb Energy (mean)', 'kJ/mol', 'lower',
        'Electrostatic interaction energy between chains. More negative values '
        'indicate stronger favorable electrostatic interactions.'),
    'lj_mean': ('LJ Energy (mean)', 'kJ/mol', 'lower',
        'Lennard-Jones (van der Waals) interaction energy. More negative values '
        'indicate better shape complementarity at the interface.'),
    'total_ie_mean': ('Total Interaction Energy', 'kJ/mol', 'lower',
        'Sum of Coulomb and LJ energies. More negative values indicate '
        'stronger overall binding affinity between the protein chains.'),
    'potential': ('Potential Energy', 'kJ/mol', 'lower',
        'Total potential energy of the system. Lower values indicate a more '
        'thermodynamically stable complex conformation.'),
    'contacts': ('Interface Contacts', 'count', 'higher',
        'Number of atomic contacts at the interface. More contacts indicate '
        'a larger, more extensive binding interface.'),
    'min_distance': ('Minimum Distance', 'nm', 'lower',
        'Closest approach between chains. Lower values indicate tighter binding.'),
}


def read_metrics(filepath: str) -> Dict[str, float]:
    """Read stability metrics from a metrics file."""
    metrics = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                if ':' in line and not line.startswith('='):
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        try:
                            metrics[parts[0].strip()] = float(parts[1].strip())
                        except ValueError:
                            pass
    return metrics


def read_json_metrics(filepath: str) -> Dict[str, Any]:
    """Read metrics from JSON file."""
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return {}


def collect_all_metrics(workdir: str, struct_dir: str) -> Dict[str, float]:
    """Collect metrics from all available sources in a structure directory."""
    struct_path = os.path.join(workdir, struct_dir)
    metrics = {}
    
    # Try stability_metrics.txt first (MD analysis)
    stability_file = os.path.join(struct_path, "stability_metrics.txt")
    if os.path.exists(stability_file):
        metrics.update(read_metrics(stability_file))
    
    # Try metrics.json (quick stability)
    json_file = os.path.join(struct_path, "metrics.json")
    if os.path.exists(json_file):
        data = read_json_metrics(json_file)
        for k, v in data.items():
            if isinstance(v, (int, float)) and v is not None:
                metrics[k] = float(v)
    
    # Try statistics/md_statistics.json
    stats_file = os.path.join(struct_path, "statistics", "md_statistics.json")
    if os.path.exists(stats_file):
        data = read_json_metrics(stats_file)
        for k, v in data.items():
            if isinstance(v, (int, float)) and v is not None:
                metrics[k] = float(v)
            elif isinstance(v, dict):
                for sub_k, sub_v in v.items():
                    if isinstance(sub_v, (int, float)) and sub_v is not None:
                        metrics[f"{k}_{sub_k}"] = float(sub_v)
    
    return metrics


def generate_md_comparison_plots(
    workdir: str,
    m1: Dict[str, float],
    m2: Dict[str, float],
    name1: str,
    name2: str
) -> List[str]:
    """Generate comparison plots for MD simulation results."""
    generated_plots = []
    plots_dir = os.path.join(workdir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available, skipping plots")
        return generated_plots
    
    # Shorten names
    short1 = name1[:20] + "..." if len(name1) > 20 else name1
    short2 = name2[:20] + "..." if len(name2) > 20 else name2
    
    # Get common metrics
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    common_keys = [k for k in all_keys 
                   if k in m1 and k in m2 
                   and m1[k] is not None and m2[k] is not None]
    
    if not common_keys:
        return generated_plots
    
    # 1. Multi-panel bar chart comparison
    n_metrics = len(common_keys)
    fig, axes = plt.subplots(min(n_metrics, 8), 1, figsize=(12, min(n_metrics, 8) * 2.5))
    if n_metrics == 1:
        axes = [axes]
    
    for ax, key in zip(axes, common_keys[:8]):
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        display_name = info[0]
        unit = info[1]
        better = info[2]
        
        v1, v2 = m1[key], m2[key]
        
        # Determine colors
        if better == 'lower':
            colors = ['green' if v1 < v2 else 'coral', 
                      'green' if v2 < v1 else 'coral']
        elif better == 'higher':
            colors = ['green' if v1 > v2 else 'coral',
                      'green' if v2 > v1 else 'coral']
        else:
            colors = ['steelblue', 'steelblue']
        
        bars = ax.barh(['S1', 'S2'], [v1, v2], color=colors)
        ax.set_xlabel(f'{display_name} ({unit})' if unit else display_name)
        
        for bar, val in zip(bars, [v1, v2]):
            ax.text(bar.get_width() + 0.02 * max(abs(v1), abs(v2)), 
                    bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', va='center', fontsize=9)
    
    plt.tight_layout()
    metrics_plot = os.path.join(plots_dir, "md_comparison_metrics.png")
    plt.savefig(metrics_plot, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(metrics_plot)
    
    # 2. Score summary plot
    scores = {'Structure 1': 0, 'Structure 2': 0}
    for key in common_keys:
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        better = info[2]
        v1, v2 = m1[key], m2[key]
        
        if better == 'lower' and v1 < v2:
            scores['Structure 1'] += 1
        elif better == 'lower' and v2 < v1:
            scores['Structure 2'] += 1
        elif better == 'higher' and v1 > v2:
            scores['Structure 1'] += 1
        elif better == 'higher' and v2 > v1:
            scores['Structure 2'] += 1
    
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['green' if scores['Structure 1'] > scores['Structure 2'] else 'coral',
              'green' if scores['Structure 2'] > scores['Structure 1'] else 'coral']
    
    bars = ax.bar([f'Structure 1\n{short1}', f'Structure 2\n{short2}'],
                  [scores['Structure 1'], scores['Structure 2']],
                  color=colors, edgecolor='black', linewidth=2)
    
    ax.set_ylabel('Favorable Metrics Count', fontsize=12)
    ax.set_title('MD Simulation Stability Score Comparison', fontsize=14)
    
    for bar, score in zip(bars, [scores['Structure 1'], scores['Structure 2']]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(score), ha='center', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    summary_plot = os.path.join(plots_dir, "md_comparison_summary.png")
    plt.savefig(summary_plot, dpi=150)
    plt.close()
    generated_plots.append(summary_plot)
    
    # Save plot data as CSV in stats folder
    stats_dir = os.path.join(workdir, "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    # Save the data used for bar chart plot
    plot_data_csv = os.path.join(stats_dir, "plot_comparison_data.csv")
    with open(plot_data_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric_key', 'display_name', 'unit', 'better_when',
                        'structure_1_value', 'structure_2_value', 'winner'])
        for key in common_keys:
            info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
            v1, v2 = m1[key], m2[key]
            better = info[2]
            
            if better == 'lower':
                winner = 'S1' if v1 < v2 else 'S2' if v2 < v1 else 'tie'
            elif better == 'higher':
                winner = 'S1' if v1 > v2 else 'S2' if v2 > v1 else 'tie'
            else:
                winner = 'n/a'
            
            writer.writerow([key, info[0], info[1], better,
                           f"{v1:.6f}", f"{v2:.6f}", winner])
    
    # Save score summary
    scores_csv = os.path.join(stats_dir, "scores_summary.csv")
    with open(scores_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['structure', 'name', 'favorable_count'])
        writer.writerow(['Structure_1', name1, scores['Structure 1']])
        writer.writerow(['Structure_2', name2, scores['Structure 2']])
    
    return generated_plots


def generate_multi_structure_md_plots(
    workdir: str,
    pdb_list: List[str],
    struct_dirs: List[str] = None
) -> List[str]:
    """
    Generate combined radar plot and ranking for multiple MD simulation structures.
    
    Args:
        workdir: Working directory containing structure subdirectories
        pdb_list: List of PDB file paths
        struct_dirs: Optional list of structure directory names (default: structure_1, structure_2, ...)
        
    Returns:
        List of generated plot file paths
    """
    generated_plots = []
    plots_dir = os.path.join(workdir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    if struct_dirs is None:
        struct_dirs = [f"structure_{i+1}" for i in range(len(pdb_list))]
    
    # Collect metrics from all structures
    all_metrics = []
    names = []
    
    for i, (pdb, struct_dir) in enumerate(zip(pdb_list, struct_dirs)):
        metrics = collect_all_metrics(workdir, struct_dir)
        if metrics:
            all_metrics.append(metrics)
            # Create short name from PDB
            basename = os.path.basename(pdb).replace('.pdb', '')
            short_name = basename[:20] + "..." if len(basename) > 20 else basename
            names.append(f"S{i+1}: {short_name}")
        else:
            print(f"  Warning: No metrics found for {struct_dir}")
    
    if len(all_metrics) < 2:
        print("  Not enough structures with metrics for combined plot")
        return generated_plots
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("  matplotlib not available, skipping combined plots")
        return generated_plots
    
    # Metrics to include in radar plot (MD-specific)
    radar_metrics = ['rmsd_mean', 'hbonds_mean', 'mindist_mean', 'total_ie_mean', 'rg_mean']
    
    # Find metrics available in all structures
    available_metrics = []
    for key in radar_metrics:
        if all(key in m and m[key] is not None for m in all_metrics):
            available_metrics.append(key)
    
    # Fall back to other common metrics if MD-specific ones not available
    if len(available_metrics) < 3:
        fallback_metrics = ['potential', 'hbonds', 'contacts', 'min_distance', 'sasa']
        for key in fallback_metrics:
            if key not in available_metrics:
                if all(key in m and m[key] is not None for m in all_metrics):
                    available_metrics.append(key)
    
    if len(available_metrics) < 3:
        print(f"  Not enough common metrics for radar plot (found: {available_metrics})")
        return generated_plots
    
    # Normalize metrics for radar plot
    normalized_data = []
    labels = []
    
    for key in available_metrics:
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        labels.append(info[0])
        
        values = [m[key] for m in all_metrics]
        min_val = min(values)
        max_val = max(values)
        range_val = max_val - min_val if max_val != min_val else 1
        
        better = info[2]
        
        # Normalize to 0-1, where 1 is "better"
        norm_values = []
        for v in values:
            normalized = (v - min_val) / range_val
            # Invert if lower is better
            if better == 'lower':
                normalized = 1 - normalized
            norm_values.append(normalized)
        
        normalized_data.append(norm_values)
    
    # Transpose for per-structure data
    per_structure_data = list(zip(*normalized_data))
    
    # Generate colors
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_metrics)))
    
    # 1. Combined Radar Plot
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection='polar'))
    
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    angles += angles[:1]  # Complete the loop
    
    for i, (data, name, color) in enumerate(zip(per_structure_data, names, colors)):
        data_closed = list(data) + [data[0]]  # Close the polygon
        ax.plot(angles, data_closed, 'o-', linewidth=2, label=name, color=color)
        ax.fill(angles, data_closed, alpha=0.15, color=color)
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, size=9)
    ax.set_title(f'Multi-Structure MD Stability Comparison\n({len(all_metrics)} structures, outer = better)', 
                 pad=20, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.0), fontsize=9)
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    radar_path = os.path.join(plots_dir, "md_combined_radar.png")
    plt.savefig(radar_path, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(radar_path)
    print(f"  Generated: {radar_path}")
    
    # 2. Stacked bar chart for ranking
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Calculate "score" for each structure (sum of normalized values)
    scores = [sum(data) for data in per_structure_data]
    
    # Sort structures by score
    sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    sorted_names = [names[i] for i in sorted_indices]
    sorted_scores = [scores[i] for i in sorted_indices]
    sorted_colors = [colors[i] for i in sorted_indices]
    
    # Create stacked bars showing contribution from each metric
    bar_width = 0.6
    y_pos = np.arange(len(sorted_names))
    
    # Stack the metrics
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
    ax.set_title('MD Simulation Structure Ranking\n(Stacked contribution from each metric)', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    
    # Add total score labels
    for i, (name, score) in enumerate(zip(sorted_names, sorted_scores)):
        ax.text(score + 0.05, i, f'{score:.2f}', va='center', fontsize=10, fontweight='bold')
    
    # Mark the winner
    ax.annotate('★ BEST', xy=(sorted_scores[0], 0), fontsize=12, fontweight='bold',
                color='green', va='center', ha='left',
                xytext=(sorted_scores[0] + 0.3, 0))
    
    plt.tight_layout()
    ranking_path = os.path.join(plots_dir, "md_combined_ranking.png")
    plt.savefig(ranking_path, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(ranking_path)
    print(f"  Generated: {ranking_path}")
    
    # 3. Individual metric comparison bars
    n_metrics = len(available_metrics)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(14, 3 * n_metrics))
    if n_metrics == 1:
        axes = [axes]
    
    for ax, key in zip(axes, available_metrics):
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        display_name = info[0]
        unit = info[1]
        better = info[2]
        
        values = [m[key] for m in all_metrics]
        
        # Determine best value
        if better == 'lower':
            best_idx = values.index(min(values))
        elif better == 'higher':
            best_idx = values.index(max(values))
        else:
            best_idx = -1
        
        bar_colors = ['green' if i == best_idx else colors[i] for i in range(len(values))]
        
        bars = ax.barh(range(len(names)), values, color=bar_colors)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.set_xlabel(f'{display_name} ({unit})' if unit else display_name)
        ax.set_title(f'{display_name} - {"Lower" if better == "lower" else "Higher" if better == "higher" else "Context"} is better')
        
        # Add value labels
        for bar, val in zip(bars, values):
            ax.text(bar.get_width() + 0.01 * max(abs(v) for v in values), 
                    bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', va='center', fontsize=9)
    
    plt.tight_layout()
    metrics_path = os.path.join(plots_dir, "md_combined_metrics.png")
    plt.savefig(metrics_path, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(metrics_path)
    print(f"  Generated: {metrics_path}")
    
    # Generate summary report
    summary_path = os.path.join(workdir, "md_comparison_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" MULTI-STRUCTURE MD STABILITY COMPARISON SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  Structures compared: {len(all_metrics)}\n")
        f.write(f"  Metrics analyzed: {', '.join(labels)}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write(" OVERALL RANKING (by combined normalized score)\n")
        f.write("-" * 80 + "\n\n")
        
        for rank, (idx, score) in enumerate(zip(sorted_indices, sorted_scores), 1):
            marker = "★" if rank == 1 else " "
            f.write(f"  {marker} Rank {rank}: {names[idx]}\n")
            f.write(f"       Combined Score: {score:.3f}\n")
            f.write(f"       PDB: {os.path.basename(pdb_list[idx])}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write(" WINNER ANALYSIS\n")
        f.write("-" * 80 + "\n\n")
        
        winner_idx = sorted_indices[0]
        f.write(f"  ╔══════════════════════════════════════════════════════════════════════╗\n")
        f.write(f"  ║  WINNER: {names[winner_idx]:<60}║\n")
        f.write(f"  ╚══════════════════════════════════════════════════════════════════════╝\n\n")
        
        f.write(f"  The winning structure shows the best combined performance across\n")
        f.write(f"  {len(available_metrics)} stability metrics from MD simulation analysis.\n\n")
        
        f.write("  Metric breakdown for winner:\n")
        for j, key in enumerate(available_metrics):
            info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
            val = all_metrics[winner_idx][key]
            f.write(f"    • {info[0]}: {val:.4f} {info[1]}\n")
        
        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write(" INTERPRETATION GUIDE\n")
        f.write("-" * 80 + "\n")
        f.write("""
  This ranking evaluates protein-protein interaction stability based on
  molecular dynamics simulation results. Higher combined scores indicate:

  • More stable protein complex during simulation (lower RMSD)
  • Stronger interfacial hydrogen bonding (higher H-bonds)  
  • Tighter interface packing (lower minimum distance)
  • More favorable interaction energies (more negative IE)
  • Consistent complex shape (stable radius of gyration)

  NOTE: For definitive binding affinity predictions, consider:
  - Longer MD simulations (100+ ns)
  - MM-PBSA/MM-GBSA free energy calculations
  - Multiple independent simulation replicas
  - Experimental validation (ITC, SPR, etc.)
""")
        f.write("=" * 80 + "\n")
    
    print(f"  Generated: {summary_path}")
    
    # ============================================================
    # Save comparison data as CSV files in 'stats' folder
    # ============================================================
    stats_dir = os.path.join(workdir, "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    # 1. Combined ranking CSV
    ranking_csv = os.path.join(stats_dir, "combined_ranking.csv")
    with open(ranking_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank', 'structure_index', 'structure_name', 'pdb_file', 
                        'combined_score', 'normalized_scores'])
        for rank, (idx, score) in enumerate(zip(sorted_indices, sorted_scores), 1):
            norm_scores = ';'.join([f"{labels[j]}={per_structure_data[idx][j]:.4f}" 
                                   for j in range(len(labels))])
            writer.writerow([
                rank,
                idx + 1,
                names[idx],
                os.path.basename(pdb_list[idx]),
                f"{score:.6f}",
                norm_scores
            ])
    print(f"  Generated: {ranking_csv}")
    
    # 2. Raw metrics for all structures CSV
    raw_metrics_csv = os.path.join(stats_dir, "all_structures_metrics.csv")
    with open(raw_metrics_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        # Get all unique metrics across all structures
        all_metric_keys = sorted(set(k for m in all_metrics for k in m.keys()))
        header = ['structure_index', 'structure_name', 'pdb_file'] + all_metric_keys
        writer.writerow(header)
        
        for i, (metrics, name, pdb) in enumerate(zip(all_metrics, names, pdb_list)):
            row = [i + 1, name, os.path.basename(pdb)]
            for key in all_metric_keys:
                val = metrics.get(key)
                if val is not None:
                    row.append(f"{val:.6f}" if isinstance(val, float) else val)
                else:
                    row.append('')
            writer.writerow(row)
    print(f"  Generated: {raw_metrics_csv}")
    
    # 3. Normalized metrics CSV (for radar plot data)
    normalized_csv = os.path.join(stats_dir, "normalized_metrics.csv")
    with open(normalized_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['structure_index', 'structure_name'] + labels
        writer.writerow(header)
        
        for i, (data, name) in enumerate(zip(per_structure_data, names)):
            row = [i + 1, name] + [f"{v:.6f}" for v in data]
            writer.writerow(row)
    print(f"  Generated: {normalized_csv}")
    
    # 4. Metric descriptions CSV (for reference)
    descriptions_csv = os.path.join(stats_dir, "metric_descriptions.csv")
    with open(descriptions_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric_key', 'display_name', 'unit', 'better_when', 'explanation'])
        for key in available_metrics:
            info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
            writer.writerow([key, info[0], info[1], info[2], info[3] if len(info) > 3 else ''])
    print(f"  Generated: {descriptions_csv}")
    
    # 5. Comparison metadata CSV
    metadata_csv = os.path.join(stats_dir, "comparison_metadata.csv")
    with open(metadata_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['key', 'value'])
        writer.writerow(['generated_at', datetime.now().isoformat()])
        writer.writerow(['num_structures', len(all_metrics)])
        writer.writerow(['num_metrics', len(available_metrics)])
        writer.writerow(['metrics_analyzed', ','.join(available_metrics)])
        writer.writerow(['winner_index', sorted_indices[0] + 1])
        writer.writerow(['winner_name', names[sorted_indices[0]]])
        writer.writerow(['winner_pdb', os.path.basename(pdb_list[sorted_indices[0]])])
        writer.writerow(['winner_score', f"{sorted_scores[0]:.6f}"])
    print(f"  Generated: {metadata_csv}")
    
    # Also save as JSON for programmatic access
    ranking_json = os.path.join(workdir, "combined_ranking.json")
    ranking_data = {
        'generated_at': datetime.now().isoformat(),
        'num_structures': len(all_metrics),
        'metrics_analyzed': available_metrics,
        'winner': {
            'rank': 1,
            'index': sorted_indices[0] + 1,
            'name': names[sorted_indices[0]],
            'pdb': os.path.basename(pdb_list[sorted_indices[0]]),
            'score': sorted_scores[0]
        },
        'ranking': [
            {
                'rank': rank,
                'index': idx + 1,
                'name': names[idx],
                'pdb': os.path.basename(pdb_list[idx]),
                'combined_score': score,
                'raw_metrics': all_metrics[idx],
                'normalized_scores': {labels[j]: per_structure_data[idx][j] 
                                     for j in range(len(labels))}
            }
            for rank, (idx, score) in enumerate(zip(sorted_indices, sorted_scores), 1)
        ]
    }
    with open(ranking_json, 'w') as f:
        json.dump(ranking_data, f, indent=2)
    print(f"  Generated: {ranking_json}")
    
    return generated_plots


def compare_md_results(
    workdir: str,
    pdb1_path: str,
    pdb2_path: str,
    struct1_dir: str = "structure_1",
    struct2_dir: str = "structure_2",
    output_file: str = "md_comparison_report.txt"
) -> Tuple[Dict[str, int], Optional[str]]:
    """
    Compare MD simulation stability results between two structures.
    
    Returns:
        Tuple of (scores dict, winner string or None)
    """
    # Collect all metrics from both structures
    m1 = collect_all_metrics(workdir, struct1_dir)
    m2 = collect_all_metrics(workdir, struct2_dir)
    
    pdb1_name = os.path.basename(pdb1_path)
    pdb2_name = os.path.basename(pdb2_path)
    
    # Build detailed report
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append(" MD SIMULATION STABILITY COMPARISON")
    report_lines.append("=" * 80)
    report_lines.append("")
    report_lines.append(f"Structure 1: {pdb1_name}")
    report_lines.append(f"  Path: {pdb1_path}")
    report_lines.append(f"Structure 2: {pdb2_name}")
    report_lines.append(f"  Path: {pdb2_path}")
    report_lines.append("")
    report_lines.append("-" * 80)
    report_lines.append(" METRIC-BY-METRIC ANALYSIS")
    report_lines.append("-" * 80)
    
    scores = {"s1": 0, "s2": 0}
    detailed_comparisons = []
    
    # Get all metrics
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    
    for key in all_keys:
        v1 = m1.get(key)
        v2 = m2.get(key)
        
        if v1 is None and v2 is None:
            continue
        
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', f'Metric: {key}'))
        display_name = info[0]
        unit = info[1]
        better = info[2]
        explanation = info[3] if len(info) > 3 else ''
        
        report_lines.append("")
        report_lines.append(f"▸ {display_name}")
        report_lines.append(f"  {'─' * 70}")
        
        v1_str = f"{v1:.4f} {unit}" if v1 is not None else "N/A"
        v2_str = f"{v2:.4f} {unit}" if v2 is not None else "N/A"
        
        report_lines.append(f"  Structure 1: {v1_str}")
        report_lines.append(f"  Structure 2: {v2_str}")
        
        winner = ""
        winner_explanation = ""
        
        if v1 is not None and v2 is not None:
            diff = abs(v1 - v2)
            pct_diff = (diff / max(abs(v1), abs(v2), 0.0001)) * 100
            
            if better == 'lower':
                if v1 < v2:
                    winner = "Structure 1"
                    scores['s1'] += 1
                    winner_explanation = f"Structure 1 is {pct_diff:.1f}% lower (more favorable)"
                elif v2 < v1:
                    winner = "Structure 2"
                    scores['s2'] += 1
                    winner_explanation = f"Structure 2 is {pct_diff:.1f}% lower (more favorable)"
            elif better == 'higher':
                if v1 > v2:
                    winner = "Structure 1"
                    scores['s1'] += 1
                    winner_explanation = f"Structure 1 is {pct_diff:.1f}% higher (more favorable)"
                elif v2 > v1:
                    winner = "Structure 2"
                    scores['s2'] += 1
                    winner_explanation = f"Structure 2 is {pct_diff:.1f}% higher (more favorable)"
            elif better == 'similar':
                winner_explanation = f"Similar values preferred (diff: {pct_diff:.1f}%)"
            else:
                winner_explanation = "Informational metric"
        
        if winner:
            report_lines.append(f"  ★ Winner: {winner} — {winner_explanation}")
            detailed_comparisons.append((display_name, winner, winner_explanation))
        else:
            report_lines.append(f"  ○ {winner_explanation}")
        
        if explanation:
            report_lines.append(f"  ")
            report_lines.append(f"  Interpretation: {explanation}")
    
    # Conclusion
    report_lines.append("")
    report_lines.append("=" * 80)
    report_lines.append(" OVERALL CONCLUSION")
    report_lines.append("=" * 80)
    report_lines.append("")
    report_lines.append(f"  Favorable metrics for Structure 1: {scores['s1']}")
    report_lines.append(f"  Favorable metrics for Structure 2: {scores['s2']}")
    report_lines.append("")
    
    winner_str = None
    if scores['s1'] > scores['s2']:
        winner_str = "Structure 1"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  WINNER: Structure 1 ({pdb1_name[:50]})  ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 1 demonstrates MORE STABLE protein-protein interaction")
        report_lines.append(f"  based on {scores['s1']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    elif scores['s2'] > scores['s1']:
        winner_str = "Structure 2"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  WINNER: Structure 2 ({pdb2_name[:50]})  ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 2 demonstrates MORE STABLE protein-protein interaction")
        report_lines.append(f"  based on {scores['s2']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    else:
        winner_str = "Tie"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  RESULT: Both structures show SIMILAR stability                      ║")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
    
    # Key factors
    if detailed_comparisons:
        report_lines.append("")
        report_lines.append("  Key factors:")
        for metric_name, winner, explanation in detailed_comparisons[:5]:
            report_lines.append(f"    • {metric_name}: {winner} — {explanation}")
    
    # Methodology
    report_lines.append("")
    report_lines.append("-" * 80)
    report_lines.append(" METHODOLOGY & INTERPRETATION")
    report_lines.append("-" * 80)
    report_lines.append("""
  This comparison evaluates protein-protein interaction stability based on
  molecular dynamics simulation results. Key stability indicators include:

  1. RMSD (Root Mean Square Deviation)
     Measures structural deviation from the starting structure.
     Lower RMSD indicates the complex maintains its fold during simulation.

  2. HYDROGEN BONDS
     More H-bonds at the interface indicate stronger, specific interactions.
     Stable H-bond count throughout simulation suggests persistent binding.

  3. INTERACTION ENERGY (Coulomb + LJ)
     More negative values indicate stronger binding between chains.
     This is a key indicator of binding affinity.

  4. MINIMUM DISTANCE
     Lower values indicate tighter packing at the interface.
     Stable minimum distance suggests consistent binding.

  5. RADIUS OF GYRATION
     Stable values indicate the complex maintains its shape.
     Large fluctuations may indicate structural instability.

  NOTE: For publication-quality binding affinity predictions, consider:
  - Longer MD simulations (100+ ns)
  - MM-PBSA or MM-GBSA free energy calculations
  - Multiple independent replicas
  - Experimental validation
""")
    report_lines.append("=" * 80)
    
    # Print and save
    report_text = "\n".join(report_lines)
    print(report_text)
    
    report_file = os.path.join(workdir, output_file)
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\n  Report saved to: {report_file}")
    
    # Generate plots
    print("\n  Generating comparison plots...")
    plots = generate_md_comparison_plots(workdir, m1, m2, pdb1_name, pdb2_name)
    if plots:
        print(f"  Generated {len(plots)} plots in {os.path.join(workdir, 'plots')}/")
    
    # ============================================================
    # Save pairwise comparison data as CSV files in 'stats' folder
    # ============================================================
    stats_dir = os.path.join(workdir, "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    # 1. Pairwise comparison metrics CSV
    pairwise_csv = os.path.join(stats_dir, "pairwise_comparison.csv")
    with open(pairwise_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric_key', 'display_name', 'unit', 'better_when',
                        'structure_1_value', 'structure_2_value', 'difference', 
                        'pct_difference', 'winner'])
        
        for key in all_keys:
            v1 = m1.get(key)
            v2 = m2.get(key)
            
            if v1 is None or v2 is None:
                continue
            
            info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
            diff = v1 - v2
            pct_diff = (abs(diff) / max(abs(v1), abs(v2), 0.0001)) * 100
            
            # Determine winner
            better = info[2]
            winner = ''
            if better == 'lower':
                winner = 'structure_1' if v1 < v2 else 'structure_2' if v2 < v1 else 'tie'
            elif better == 'higher':
                winner = 'structure_1' if v1 > v2 else 'structure_2' if v2 > v1 else 'tie'
            else:
                winner = 'n/a'
            
            writer.writerow([
                key, info[0], info[1], better,
                f"{v1:.6f}", f"{v2:.6f}", f"{diff:.6f}", f"{pct_diff:.2f}",
                winner
            ])
    print(f"  Generated: {pairwise_csv}")
    
    # 2. Structure metrics CSV (both structures side by side)
    structures_csv = os.path.join(stats_dir, "structures_metrics.csv")
    with open(structures_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['metric'] + [f'structure_1_{pdb1_name}', f'structure_2_{pdb2_name}']
        writer.writerow(header)
        
        for key in sorted(all_keys):
            v1 = m1.get(key)
            v2 = m2.get(key)
            row = [key]
            row.append(f"{v1:.6f}" if v1 is not None else '')
            row.append(f"{v2:.6f}" if v2 is not None else '')
            writer.writerow(row)
    print(f"  Generated: {structures_csv}")
    
    # 3. Comparison summary CSV
    summary_csv = os.path.join(stats_dir, "comparison_summary.csv")
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['key', 'value'])
        writer.writerow(['generated_at', datetime.now().isoformat()])
        writer.writerow(['structure_1_pdb', pdb1_name])
        writer.writerow(['structure_2_pdb', pdb2_name])
        writer.writerow(['structure_1_path', pdb1_path])
        writer.writerow(['structure_2_path', pdb2_path])
        writer.writerow(['structure_1_favorable_metrics', scores['s1']])
        writer.writerow(['structure_2_favorable_metrics', scores['s2']])
        writer.writerow(['winner', winner_str])
        writer.writerow(['total_comparable_metrics', scores['s1'] + scores['s2']])
    print(f"  Generated: {summary_csv}")
    
    # 4. Detailed comparisons CSV (key factors)
    if detailed_comparisons:
        factors_csv = os.path.join(stats_dir, "key_factors.csv")
        with open(factors_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['metric_name', 'winner', 'explanation'])
            for metric_name, winner, explanation in detailed_comparisons:
                writer.writerow([metric_name, winner, explanation])
        print(f"  Generated: {factors_csv}")
    
    return scores, winner_str


def main():
    """Command-line interface for MD results comparison."""
    parser = argparse.ArgumentParser(
        description='Compare MD simulation stability results between structures'
    )
    parser.add_argument('--workdir', required=True,
                       help='Base working directory')
    parser.add_argument('--multi', action='store_true',
                       help='Multi-structure comparison mode (generates radar plot and ranking)')
    parser.add_argument('--pdb-list', nargs='+',
                       help='List of PDB file paths (for --multi mode)')
    parser.add_argument('--struct-dirs', nargs='+',
                       help='List of structure directory names (for --multi mode)')
    parser.add_argument('--pdb1',
                       help='Path to first PDB file (pairwise mode)')
    parser.add_argument('--pdb2',
                       help='Path to second PDB file (pairwise mode)')
    parser.add_argument('--struct1-dir', default='structure_1',
                       help='Subdirectory for structure 1')
    parser.add_argument('--struct2-dir', default='structure_2',
                       help='Subdirectory for structure 2')
    parser.add_argument('--output', '-o', default='md_comparison_report.txt',
                       help='Output report filename')
    
    args = parser.parse_args()
    
    if args.multi:
        # Multi-structure comparison mode
        if not args.pdb_list:
            parser.error("--pdb-list is required for --multi mode")
        
        pdb_list = args.pdb_list
        struct_dirs = args.struct_dirs if args.struct_dirs else [f"structure_{i+1}" for i in range(len(pdb_list))]
        
        print(f"  Multi-structure MD comparison: {len(pdb_list)} structures")
        plots = generate_multi_structure_md_plots(args.workdir, pdb_list, struct_dirs)
        if plots:
            print(f"  Generated {len(plots)} plots")
    else:
        # Pairwise comparison mode
        if not args.pdb1 or not args.pdb2:
            parser.error("--pdb1 and --pdb2 are required for pairwise comparison")
        
        compare_md_results(
            args.workdir,
            args.pdb1,
            args.pdb2,
            args.struct1_dir,
            args.struct2_dir,
            args.output
        )


if __name__ == '__main__':
    main()
