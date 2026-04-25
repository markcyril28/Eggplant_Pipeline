#!/usr/bin/env python3
"""
Results Comparator Module for GROMACS Stability Analysis

Compares stability metrics between two structures and generates reports
with detailed explanations and visualizations.
"""

import os
import json
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Any


# Metric definitions: (display_name, better_when, explanation, unit)
METRIC_INFO = {
    'potential': ('Potential Energy', 'lower', 
        'Total potential energy of the system. Lower values indicate a more thermodynamically '
        'stable conformation. The protein complex has found a lower energy state, suggesting '
        'more favorable atomic interactions and bond geometries.', 'kJ/mol'),
    'potential_energy': ('Potential Energy', 'lower',
        'Total potential energy of the system. Lower values indicate a more thermodynamically '
        'stable conformation.', 'kJ/mol'),
    'hbonds': ('Hydrogen Bonds', 'higher',
        'Number of hydrogen bonds formed between the two protein chains at the interface. '
        'More H-bonds indicate stronger, more specific intermolecular interactions. H-bonds '
        'are critical for protein-protein recognition and binding specificity.', 'count'),
    'min_distance': ('Minimum Distance', 'lower',
        'Closest distance between atoms of the two chains. Lower values indicate tighter '
        'packing at the interface. Very low values (<0.15 nm) suggest intimate contact, '
        'while values >0.3 nm may indicate weaker association.', 'nm'),
    'contacts': ('Interface Contacts', 'higher',
        'Number of atomic contacts within 0.6 nm between the two chains. More contacts '
        'indicate a larger interaction interface and potentially stronger binding. '
        'This correlates with binding affinity in many protein complexes.', 'count'),
    'contacts_0.6nm': ('Interface Contacts (<0.6nm)', 'higher',
        'Number of atomic contacts within 0.6 nm between the two chains.', 'count'),
    'sasa': ('Solvent Accessible Surface Area', 'context',
        'Total surface area of the complex accessible to solvent. This is informational - '
        'the buried surface area (reduction upon binding) is more relevant for binding strength.', 'nm²'),
    'buried_surface_area': ('Buried Surface Area', 'higher',
        'Surface area buried upon complex formation. Larger buried surface indicates '
        'a more extensive binding interface, often correlating with tighter binding.', 'nm²'),
    'sasa_total': ('Total SASA', 'context',
        'Total solvent accessible surface area of the complex.', 'nm²'),
    'sasa_chainA': ('SASA Chain A', 'context', 'Solvent accessible surface area of chain A.', 'nm²'),
    'sasa_chainB': ('SASA Chain B', 'context', 'Solvent accessible surface area of chain B.', 'nm²'),
    'gyration': ('Radius of Gyration', 'context',
        'Measure of the overall compactness of the structure. Lower values indicate '
        'a more compact structure.', 'nm'),
}


def read_metrics_txt(filepath: str) -> Dict[str, float]:
    """Read metrics from a metrics.txt file."""
    metrics = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                if ':' in line:
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        try:
                            metrics[parts[0].strip()] = float(parts[1].strip())
                        except ValueError:
                            pass
    return metrics


def read_metrics_json(filepath: str) -> Dict[str, Any]:
    """Read metrics from a metrics.json file."""
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return {}


def read_metrics(workdir: str, struct_dir: str) -> Dict[str, float]:
    """Read metrics from either JSON or TXT file."""
    json_path = os.path.join(workdir, struct_dir, "metrics.json")
    txt_path = os.path.join(workdir, struct_dir, "metrics.txt")
    
    # Prefer JSON
    if os.path.exists(json_path):
        data = read_metrics_json(json_path)
        # Flatten if needed
        metrics = {}
        for k, v in data.items():
            if isinstance(v, (int, float)) and v is not None:
                metrics[k] = float(v)
        return metrics
    
    return read_metrics_txt(txt_path)


def generate_comparison_plots(
    workdir: str,
    m1: Dict[str, float],
    m2: Dict[str, float],
    name1: str,
    name2: str
) -> List[str]:
    """
    Generate comparison plots for the two structures.
    
    Returns list of generated plot paths.
    """
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
    
    # Shorten names for display
    short1 = name1[:25] + "..." if len(name1) > 25 else name1
    short2 = name2[:25] + "..." if len(name2) > 25 else name2
    
    # 1. Bar chart comparison of all metrics
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    
    # Filter to common metrics with values
    common_keys = [k for k in all_keys 
                   if k in m1 and k in m2 
                   and m1[k] is not None and m2[k] is not None]
    
    if common_keys:
        fig, axes = plt.subplots(len(common_keys), 1, figsize=(12, 3 * len(common_keys)))
        if len(common_keys) == 1:
            axes = [axes]
        
        for ax, key in zip(axes, common_keys):
            info = METRIC_INFO.get(key, (key, 'context', '', ''))
            display_name = info[0]
            better = info[1]
            unit = info[3] if len(info) > 3 else ''
            
            v1, v2 = m1[key], m2[key]
            
            # Determine colors based on which is better
            if better == 'lower':
                colors = ['green' if v1 < v2 else 'coral', 
                          'green' if v2 < v1 else 'coral']
            elif better == 'higher':
                colors = ['green' if v1 > v2 else 'coral',
                          'green' if v2 > v1 else 'coral']
            else:
                colors = ['steelblue', 'steelblue']
            
            bars = ax.barh(['Structure 1', 'Structure 2'], [v1, v2], color=colors)
            ax.set_xlabel(f'{display_name} ({unit})' if unit else display_name)
            ax.set_title(f'{display_name} Comparison')
            
            # Add value labels
            for bar, val in zip(bars, [v1, v2]):
                ax.text(bar.get_width() + 0.01 * max(abs(v1), abs(v2)), 
                        bar.get_y() + bar.get_height()/2,
                        f'{val:.2f}', va='center', fontsize=10)
            
            # Add legend for color meaning
            if better != 'context':
                ax.annotate(f'Green = Better ({better} is better)', 
                            xy=(0.98, 0.02), xycoords='axes fraction',
                            ha='right', fontsize=8, color='gray')
        
        plt.tight_layout()
        bar_plot_path = os.path.join(plots_dir, "comparison_metrics.png")
        plt.savefig(bar_plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        generated_plots.append(bar_plot_path)
        print(f"  Generated: {bar_plot_path}")
    
    # 2. Radar/Spider plot for normalized metrics
    # Normalize metrics to 0-1 scale for comparison
    comparison_metrics = ['potential', 'hbonds', 'contacts', 'min_distance', 'sasa']
    available = [k for k in comparison_metrics if k in m1 and k in m2]
    
    if len(available) >= 3:
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        
        # Normalize values (handle direction - some lower is better)
        normalized_1 = []
        normalized_2 = []
        labels = []
        
        for key in available:
            v1, v2 = m1[key], m2[key]
            info = METRIC_INFO.get(key, (key, 'context', '', ''))
            better = info[1]
            
            # Normalize to 0-1 where 1 is "better"
            max_val = max(abs(v1), abs(v2)) if max(abs(v1), abs(v2)) > 0 else 1
            
            if better == 'lower':
                # For "lower is better", invert so higher normalized = better
                n1 = 1 - (v1 / (max_val * 2))
                n2 = 1 - (v2 / (max_val * 2))
            else:
                n1 = v1 / (max_val * 1.2)
                n2 = v2 / (max_val * 1.2)
            
            normalized_1.append(max(0, min(1, n1)))
            normalized_2.append(max(0, min(1, n2)))
            labels.append(info[0])
        
        # Complete the loop
        angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
        normalized_1 += normalized_1[:1]
        normalized_2 += normalized_2[:1]
        angles += angles[:1]
        
        ax.plot(angles, normalized_1, 'o-', linewidth=2, label=f'Structure 1: {short1}', color='blue')
        ax.fill(angles, normalized_1, alpha=0.25, color='blue')
        ax.plot(angles, normalized_2, 'o-', linewidth=2, label=f'Structure 2: {short2}', color='red')
        ax.fill(angles, normalized_2, alpha=0.25, color='red')
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels)
        ax.set_title('Normalized Stability Metrics\n(Outer = Better)', pad=20)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
        
        plt.tight_layout()
        radar_plot_path = os.path.join(plots_dir, "comparison_radar.png")
        plt.savefig(radar_plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        generated_plots.append(radar_plot_path)
        print(f"  Generated: {radar_plot_path}")
    
    # 3. Winner summary plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    scores = {'Structure 1': 0, 'Structure 2': 0}
    metric_winners = []
    
    for key in common_keys:
        info = METRIC_INFO.get(key, (key, 'context', '', ''))
        better = info[1]
        v1, v2 = m1[key], m2[key]
        
        if better == 'lower':
            if v1 < v2:
                scores['Structure 1'] += 1
                metric_winners.append((info[0], 'S1', v1, v2))
            elif v2 < v1:
                scores['Structure 2'] += 1
                metric_winners.append((info[0], 'S2', v1, v2))
        elif better == 'higher':
            if v1 > v2:
                scores['Structure 1'] += 1
                metric_winners.append((info[0], 'S1', v1, v2))
            elif v2 > v1:
                scores['Structure 2'] += 1
                metric_winners.append((info[0], 'S2', v1, v2))
    
    colors = ['green' if scores['Structure 1'] > scores['Structure 2'] else 'coral',
              'green' if scores['Structure 2'] > scores['Structure 1'] else 'coral']
    
    bars = ax.bar(['Structure 1\n' + short1, 'Structure 2\n' + short2], 
                  [scores['Structure 1'], scores['Structure 2']], 
                  color=colors, edgecolor='black', linewidth=2)
    
    ax.set_ylabel('Number of Favorable Metrics', fontsize=12)
    ax.set_title('Stability Score Comparison\n(Higher = More Favorable Metrics)', fontsize=14)
    
    for bar, score in zip(bars, [scores['Structure 1'], scores['Structure 2']]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(score), ha='center', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    summary_plot_path = os.path.join(plots_dir, "comparison_summary.png")
    plt.savefig(summary_plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(summary_plot_path)
    print(f"  Generated: {summary_plot_path}")
    
    return generated_plots


def compare_stability_results(
    workdir: str,
    pdb1_path: str,
    pdb2_path: str,
    struct1_dir: str = "structure_1",
    struct2_dir: str = "structure_2",
    output_file: str = "comparison_summary.txt"
) -> Tuple[Dict[str, int], Optional[str]]:
    """
    Compare stability results between two structures with detailed explanations.
    
    Args:
        workdir: Base working directory
        pdb1_path: Path to first PDB file
        pdb2_path: Path to second PDB file
        struct1_dir: Subdirectory name for structure 1
        struct2_dir: Subdirectory name for structure 2
        output_file: Name of output comparison file
        
    Returns:
        Tuple of (scores dict, winner string or None)
    """
    # Read metrics
    m1 = read_metrics(workdir, struct1_dir)
    m2 = read_metrics(workdir, struct2_dir)
    
    pdb1_name = os.path.basename(pdb1_path)
    pdb2_name = os.path.basename(pdb2_path)
    
    # Build detailed report
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append(" PROTEIN-PROTEIN INTERACTION STABILITY COMPARISON")
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
    
    scores = {'s1': 0, 's2': 0}
    detailed_comparisons = []
    
    # Compare each metric
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    
    for key in all_keys:
        v1 = m1.get(key)
        v2 = m2.get(key)
        
        if v1 is None and v2 is None:
            continue
            
        info = METRIC_INFO.get(key, (key, 'context', f'Metric: {key}', ''))
        display_name = info[0]
        better = info[1]
        explanation = info[2] if len(info) > 2 else ''
        unit = info[3] if len(info) > 3 else ''
        
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
                else:
                    winner_explanation = "Both structures are equal"
            elif better == 'higher':
                if v1 > v2:
                    winner = "Structure 1"
                    scores['s1'] += 1
                    winner_explanation = f"Structure 1 is {pct_diff:.1f}% higher (more favorable)"
                elif v2 > v1:
                    winner = "Structure 2"
                    scores['s2'] += 1
                    winner_explanation = f"Structure 2 is {pct_diff:.1f}% higher (more favorable)"
                else:
                    winner_explanation = "Both structures are equal"
            else:
                winner_explanation = "Informational metric (no preference)"
        
        if winner:
            report_lines.append(f"  ★ Winner: {winner} — {winner_explanation}")
            detailed_comparisons.append((display_name, winner, winner_explanation))
        else:
            report_lines.append(f"  ○ {winner_explanation}")
        
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
        report_lines.append(f"  ║  WINNER: Structure 1 ({pdb1_name})                                   ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 1 demonstrates MORE STABLE protein-protein interaction based on")
        report_lines.append(f"  {scores['s1']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    elif scores['s2'] > scores['s1']:
        winner_str = "Structure 2"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  WINNER: Structure 2 ({pdb2_name})                                   ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 2 demonstrates MORE STABLE protein-protein interaction based on")
        report_lines.append(f"  {scores['s2']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    else:
        winner_str = "Tie"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  RESULT: Both structures show SIMILAR stability                      ║")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
    
    # Add specific winning metrics
    if detailed_comparisons:
        report_lines.append("")
        report_lines.append("  Key factors:")
        for metric_name, winner, explanation in detailed_comparisons:
            report_lines.append(f"    • {metric_name}: {winner} — {explanation}")
    
    # Add methodology note
    report_lines.append("")
    report_lines.append("-" * 80)
    report_lines.append(" METHODOLOGY & INTERPRETATION GUIDE")
    report_lines.append("-" * 80)
    report_lines.append("""
  This comparison evaluates protein-protein interaction (PPI) stability based on
  energy minimization results. Key indicators for a stable PPI include:

  1. POTENTIAL ENERGY (Lower = Better)
     The total potential energy reflects the thermodynamic stability of the complex.
     A lower value indicates more favorable atomic interactions (bonds, angles,
     electrostatics, van der Waals forces).

  2. HYDROGEN BONDS (Higher = Better)
     H-bonds at the interface provide specificity and strength to the interaction.
     More H-bonds typically correlate with tighter, more specific binding.

  3. INTERFACE CONTACTS (Higher = Better)
     The number of atomic contacts (<0.6 nm) at the interface indicates the extent
     of the binding surface. More contacts often correlate with higher affinity.

  4. MINIMUM DISTANCE (Lower = Better)
     The closest approach between chains indicates how tightly they pack together.
     Very short distances (<0.15 nm) suggest intimate contact.

  5. BURIED SURFACE AREA (Higher = Better)
     The surface area hidden from solvent upon binding. Larger values indicate
     a more extensive interface, often correlating with stronger binding.

  NOTE: This is a quick stability assessment based on energy minimization only.
  For more accurate binding predictions, consider:
  - Full MD simulation with trajectory analysis
  - Free energy perturbation or MM-PBSA calculations
  - Experimental validation (ITC, SPR, etc.)
""")
    report_lines.append("=" * 80)
    
    # Print to console
    report_text = "\n".join(report_lines)
    print(report_text)
    
    # Save detailed report
    report_file = os.path.join(workdir, output_file)
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\n  Report saved to: {report_file}")
    
    # Generate plots
    print("\n  Generating comparison plots...")
    generated_plots = generate_comparison_plots(workdir, m1, m2, pdb1_name, pdb2_name)
    
    if generated_plots:
        print(f"  Generated {len(generated_plots)} plots in {os.path.join(workdir, 'plots')}/")
    
    return scores, winner_str


def generate_multi_structure_radar(
    workdir: str,
    pdb_list: List[str],
    struct_dirs: List[str] = None
) -> List[str]:
    """
    Generate a combined radar plot comparing all structures.
    
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
        metrics = read_metrics(workdir, struct_dir)
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
    
    # Metrics to include in radar plot
    radar_metrics = ['potential', 'hbonds', 'contacts', 'min_distance', 'sasa']
    
    # Find metrics available in all structures
    available_metrics = []
    for key in radar_metrics:
        if all(key in m and m[key] is not None for m in all_metrics):
            available_metrics.append(key)
    
    if len(available_metrics) < 3:
        print(f"  Not enough common metrics for radar plot (found: {available_metrics})")
        return generated_plots
    
    # Normalize metrics for radar plot
    # For each metric, find min/max across all structures
    normalized_data = []
    labels = []
    
    for key in available_metrics:
        info = METRIC_INFO.get(key, (key, 'context', '', ''))
        labels.append(info[0])
        
        values = [m[key] for m in all_metrics]
        min_val = min(values)
        max_val = max(values)
        range_val = max_val - min_val if max_val != min_val else 1
        
        better = info[1]
        
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
    
    # Generate colors for structures
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
    ax.set_xticklabels(labels, size=10)
    ax.set_title(f'Multi-Structure Stability Comparison\n({len(all_metrics)} structures, outer = better)', 
                 pad=20, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.0), fontsize=9)
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    radar_path = os.path.join(plots_dir, "combined_radar.png")
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
    ax.set_title('Structure Ranking by Stability Metrics\n(Stacked contribution from each metric)', 
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
    ranking_path = os.path.join(plots_dir, "combined_ranking.png")
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
        info = METRIC_INFO.get(key, (key, 'context', '', ''))
        display_name = info[0]
        better = info[1]
        unit = info[3] if len(info) > 3 else ''
        
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
                    f'{val:.2f}', va='center', fontsize=9)
    
    plt.tight_layout()
    metrics_path = os.path.join(plots_dir, "combined_metrics.png")
    plt.savefig(metrics_path, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(metrics_path)
    print(f"  Generated: {metrics_path}")
    
    # Generate summary report
    summary_path = os.path.join(workdir, "comparison_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(" MULTI-STRUCTURE STABILITY COMPARISON SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  Structures compared: {len(all_metrics)}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write(" OVERALL RANKING (by combined normalized score)\n")
        f.write("-" * 80 + "\n\n")
        
        for rank, (idx, score) in enumerate(zip(sorted_indices, sorted_scores), 1):
            marker = "★" if rank == 1 else " "
            f.write(f"  {marker} Rank {rank}: {names[idx]}\n")
            f.write(f"       Combined Score: {score:.3f}\n")
            f.write(f"       PDB: {os.path.basename(pdb_list[idx])}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write(" METRIC-BY-METRIC BEST PERFORMERS\n")
        f.write("-" * 80 + "\n\n")
        
        for key in available_metrics:
            info = METRIC_INFO.get(key, (key, 'context', '', ''))
            display_name = info[0]
            better = info[1]
            unit = info[3] if len(info) > 3 else ''
            
            values = [m[key] for m in all_metrics]
            
            if better == 'lower':
                best_idx = values.index(min(values))
                best_val = min(values)
            elif better == 'higher':
                best_idx = values.index(max(values))
                best_val = max(values)
            else:
                best_idx = 0
                best_val = values[0]
            
            f.write(f"  {display_name}:\n")
            f.write(f"    Best: {names[best_idx]} = {best_val:.3f} {unit}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write(" PLOTS GENERATED\n")
        f.write("=" * 80 + "\n\n")
        for plot in generated_plots:
            f.write(f"  - {os.path.basename(plot)}\n")
        f.write("\n")
    
    print(f"  Summary saved to: {summary_path}")
    
    return generated_plots


def main():
    """Command-line interface for results comparison."""
    parser = argparse.ArgumentParser(
        description='Compare stability results between structures'
    )
    parser.add_argument('--workdir', required=True,
                       help='Base working directory')
    parser.add_argument('--pdb1',
                       help='Path to first PDB file (for pairwise comparison)')
    parser.add_argument('--pdb2',
                       help='Path to second PDB file (for pairwise comparison)')
    parser.add_argument('--struct1-dir', default='structure_1',
                       help='Subdirectory for structure 1')
    parser.add_argument('--struct2-dir', default='structure_2',
                       help='Subdirectory for structure 2')
    parser.add_argument('--output', default='comparison_summary.txt',
                       help='Output filename for comparison report')
    # Multi-structure mode arguments
    parser.add_argument('--multi', action='store_true',
                       help='Enable multi-structure comparison mode')
    parser.add_argument('--pdb-list', nargs='+',
                       help='List of PDB files for multi-structure comparison')
    parser.add_argument('--struct-dirs', nargs='+',
                       help='List of structure directories (default: structure_1, structure_2, ...)')
    
    args = parser.parse_args()
    
    if args.multi or args.pdb_list:
        # Multi-structure comparison mode
        if not args.pdb_list:
            # Auto-detect structures in workdir
            import glob
            struct_dirs = sorted([d for d in os.listdir(args.workdir) 
                                  if os.path.isdir(os.path.join(args.workdir, d)) 
                                  and d.startswith('structure_')])
            if len(struct_dirs) < 2:
                print("Error: Not enough structure directories found")
                return
            # Create placeholder PDB list
            pdb_list = [f"structure_{i+1}.pdb" for i in range(len(struct_dirs))]
        else:
            pdb_list = args.pdb_list
            struct_dirs = args.struct_dirs if args.struct_dirs else None
        
        print(f"\n  Running multi-structure comparison ({len(pdb_list)} structures)...")
        generate_multi_structure_radar(args.workdir, pdb_list, struct_dirs)
    else:
        # Pairwise comparison mode (original behavior)
        if not args.pdb1 or not args.pdb2:
            parser.error("--pdb1 and --pdb2 are required for pairwise comparison")
        
        compare_stability_results(
            args.workdir,
            args.pdb1,
            args.pdb2,
            args.struct1_dir,
            args.struct2_dir,
            args.output
        )


if __name__ == '__main__':
    main()
