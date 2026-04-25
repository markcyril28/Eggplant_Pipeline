#!/usr/bin/env python3
"""
Plotting Module for GROMACS PPI Analysis

Provides centralized matplotlib-based plotting functions for analysis outputs.
"""

import os
from pathlib import Path
from typing import Optional, Union, List, Tuple, Dict, Any
import warnings


def check_matplotlib():
    """Check if matplotlib is available."""
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        return True
    except ImportError:
        return False


def plot_comparison_barplot(data: Dict[str, float],
                            output_file: Union[str, Path],
                            title: str = "Structure Comparison",
                            xlabel: str = "Potential Energy (kJ/mol)",
                            highlight_best: bool = True,
                            figsize: Tuple[int, int] = (12, 6),
                            dpi: int = 150) -> Optional[Path]:
    """
    Create a horizontal bar plot for structure comparison.
    
    Args:
        data: Dictionary mapping structure names to values
        output_file: Path to save the plot
        title: Plot title
        xlabel: X-axis label
        highlight_best: Whether to highlight the best (first) bar
        figsize: Figure size tuple
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        print("matplotlib not available, skipping plot")
        return None
    
    import matplotlib.pyplot as plt
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    names = list(data.keys())
    values = list(data.values())
    
    fig, ax = plt.subplots(figsize=(figsize[0], max(figsize[1], len(names) * 0.4)))
    
    colors = ['green' if i == 0 and highlight_best else 'steelblue' 
              for i in range(len(names))]
    
    ax.barh(names, values, color=colors)
    ax.set_xlabel(xlabel)
    ax.set_title(f"{title} (Lower = More Stable, Green = Best)" if highlight_best else title)
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    
    return output_path


def plot_scatter_comparison(df_or_dict: Any,
                            x_col: str,
                            y_col: str,
                            color_col: str,
                            output_file: Union[str, Path],
                            xlabel: str = "",
                            ylabel: str = "",
                            colorbar_label: str = "",
                            title: str = "",
                            figsize: Tuple[int, int] = (10, 8),
                            dpi: int = 150) -> Optional[Path]:
    """
    Create a scatter plot with color-coded points.
    
    Args:
        df_or_dict: DataFrame or dict with data
        x_col, y_col, color_col: Column names for x, y, and color
        output_file: Path to save the plot
        xlabel, ylabel: Axis labels
        colorbar_label: Label for color bar
        title: Plot title
        figsize: Figure size
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    import matplotlib.pyplot as plt
    
    try:
        import pandas as pd
        if isinstance(df_or_dict, dict):
            df = pd.DataFrame(df_or_dict)
        else:
            df = df_or_dict
        
        df_valid = df.dropna(subset=[x_col, y_col, color_col])
        if len(df_valid) < 2:
            return None
    except Exception:
        return None
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(df_valid[x_col], df_valid[y_col],
                         c=df_valid[color_col], cmap='RdYlGn_r',
                         s=100, edgecolors='black')
    plt.colorbar(scatter, label=colorbar_label or color_col)
    ax.set_xlabel(xlabel or x_col)
    ax.set_ylabel(ylabel or y_col)
    ax.set_title(title)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    
    return output_path


def plot_timeseries(times: List[float],
                    values: List[float],
                    output_file: Union[str, Path],
                    xlabel: str = "Time (ps)",
                    ylabel: str = "Value",
                    title: str = "",
                    color: str = 'blue',
                    figsize: Tuple[int, int] = (10, 6),
                    dpi: int = 150) -> Optional[Path]:
    """
    Create a simple timeseries plot.
    
    Args:
        times: List of time values
        values: List of y values
        output_file: Path to save the plot
        xlabel, ylabel: Axis labels
        title: Plot title
        color: Line color
        figsize: Figure size
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    import matplotlib.pyplot as plt
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(times, values, color=color, linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    
    return output_path


def plot_rmsd(xvg_file: Union[str, Path],
              output_file: Union[str, Path],
              title: str = "RMSD",
              dpi: int = 150) -> Optional[Path]:
    """
    Plot RMSD from XVG file.
    
    Args:
        xvg_file: Path to RMSD XVG file
        output_file: Path to save the plot
        title: Plot title
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    from .xvg_parser import parse_xvg
    
    data, meta = parse_xvg(Path(xvg_file))
    if not data:
        return None
    
    times = [d[0] for d in data]
    values = [d[1] for d in data]
    
    return plot_timeseries(
        times, values, output_file,
        xlabel=meta.get('xlabel', 'Time (ps)'),
        ylabel=meta.get('ylabel', 'RMSD (nm)'),
        title=title,
        color='blue',
        dpi=dpi
    )


def plot_rmsf(xvg_file: Union[str, Path],
              output_file: Union[str, Path],
              title: str = "RMSF per Residue",
              dpi: int = 150) -> Optional[Path]:
    """
    Plot RMSF from XVG file.
    
    Args:
        xvg_file: Path to RMSF XVG file
        output_file: Path to save the plot
        title: Plot title
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    import matplotlib.pyplot as plt
    from .xvg_parser import parse_xvg
    
    data, meta = parse_xvg(Path(xvg_file))
    if not data:
        return None
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    residues = [d[0] for d in data]
    values = [d[1] for d in data]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(residues, values, color='steelblue', width=0.8)
    ax.set_xlabel('Residue')
    ax.set_ylabel(meta.get('ylabel', 'RMSF (nm)'))
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    
    return output_path


def plot_gyration(xvg_file: Union[str, Path],
                  output_file: Union[str, Path],
                  title: str = "Radius of Gyration",
                  dpi: int = 150) -> Optional[Path]:
    """
    Plot radius of gyration from XVG file.
    
    Args:
        xvg_file: Path to gyration XVG file
        output_file: Path to save the plot
        title: Plot title
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    from .xvg_parser import parse_xvg
    
    data, meta = parse_xvg(Path(xvg_file))
    if not data:
        return None
    
    times = [d[0] for d in data]
    values = [d[1] for d in data]
    
    return plot_timeseries(
        times, values, output_file,
        xlabel=meta.get('xlabel', 'Time (ps)'),
        ylabel='Radius of Gyration (nm)',
        title=title,
        color='green',
        dpi=dpi
    )


def plot_hbonds(xvg_file: Union[str, Path],
                output_file: Union[str, Path],
                title: str = "Hydrogen Bonds",
                dpi: int = 150) -> Optional[Path]:
    """
    Plot hydrogen bonds from XVG file.
    
    Args:
        xvg_file: Path to hbonds XVG file
        output_file: Path to save the plot
        title: Plot title
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    from .xvg_parser import parse_xvg
    
    data, meta = parse_xvg(Path(xvg_file))
    if not data:
        return None
    
    times = [d[0] for d in data]
    values = [d[1] for d in data if len(d) > 1]
    
    if not values:
        return None
    
    return plot_timeseries(
        times[:len(values)], values, output_file,
        xlabel=meta.get('xlabel', 'Time (ps)'),
        ylabel='Number of H-bonds',
        title=title,
        color='purple',
        dpi=dpi
    )


def plot_focused_contact_map(interface_residues_file: Union[str, Path],
                              output_file: Union[str, Path],
                              title: str = "Interface Contact Map (Critical Residues)",
                              cmap: str = "YlOrRd",
                              figsize: Tuple[int, int] = (12, 10),
                              dpi: int = 150,
                              max_pairs: int = 50) -> Optional[Path]:
    """
    Plot a focused contact map showing only residues with interactions.
    
    Args:
        interface_residues_file: Path to interface_residues.txt file
        output_file: Path to save the plot
        title: Plot title
        cmap: Colormap name
        figsize: Figure size
        dpi: Output resolution
        max_pairs: Maximum number of pairs to show
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    interface_file = Path(interface_residues_file)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if not interface_file.exists():
        return None
    
    # Parse interface residues file
    pairs = []
    with open(interface_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip headers
            if not line or line.startswith(('Interface', '=', '-', 'Chain', 'Total')):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    res_a = parts[0]  # e.g., "ASN178"
                    res_b = parts[1]  # e.g., "THR179"
                    dist = float(parts[2])
                    pairs.append((res_a, res_b, dist))
                except (ValueError, IndexError):
                    continue
    
    if not pairs:
        return None
    
    # Limit to top pairs by distance
    pairs = pairs[:max_pairs]
    
    # Get unique residues
    chain_a_res = list(dict.fromkeys([p[0] for p in pairs]))
    chain_b_res = list(dict.fromkeys([p[1] for p in pairs]))
    
    # Create contact matrix
    matrix = np.zeros((len(chain_a_res), len(chain_b_res)))
    for res_a, res_b, dist in pairs:
        if res_a in chain_a_res and res_b in chain_b_res:
            i = chain_a_res.index(res_a)
            j = chain_b_res.index(res_b)
            # Use inverse distance as intensity (closer = stronger)
            matrix[i, j] = 1.0 / dist if dist > 0 else 1.0
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest')
    cbar = plt.colorbar(im, label='Contact Strength (1/distance)')
    
    # Set tick labels
    ax.set_xticks(range(len(chain_b_res)))
    ax.set_xticklabels(chain_b_res, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(chain_a_res)))
    ax.set_yticklabels(chain_a_res, fontsize=8)
    
    ax.set_xlabel('Chain B Residue', fontsize=12)
    ax.set_ylabel('Chain A Residue', fontsize=12)
    ax.set_title(title, fontsize=14)
    
    # Add grid
    ax.set_xticks(np.arange(-0.5, len(chain_b_res), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(chain_a_res), 1), minor=True)
    ax.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_path


def plot_contact_map(contact_matrix: Any,
                     output_file: Union[str, Path],
                     title: str = "Inter-chain Contact Map",
                     xlabel: str = "Chain B Residue",
                     ylabel: str = "Chain A Residue",
                     cmap: str = "Reds",
                     figsize: Tuple[int, int] = (10, 8),
                     dpi: int = 150) -> Optional[Path]:
    """
    Plot a contact map heatmap.
    
    Args:
        contact_matrix: 2D numpy array or path to contact map file
        output_file: Path to save the plot
        title, xlabel, ylabel: Labels
        cmap: Colormap name
        figsize: Figure size
        dpi: Output resolution
    
    Returns:
        Path to saved file or None if failed
    """
    if not check_matplotlib():
        return None
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Load data if path provided
    if isinstance(contact_matrix, (str, Path)):
        try:
            data = np.loadtxt(contact_matrix)
        except Exception:
            return None
    else:
        data = contact_matrix
    
    if data.size == 0:
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data, cmap=cmap, aspect='auto')
    plt.colorbar(im, label='Contact')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    
    return output_path


def plot_all_md_analysis(analysis_dir: Union[str, Path],
                         plots_dir: Union[str, Path],
                         dpi: int = 150) -> Dict[str, Path]:
    """
    Generate all standard MD analysis plots.
    
    Args:
        analysis_dir: Directory containing XVG analysis files
        plots_dir: Directory to save plots
        dpi: Output resolution
    
    Returns:
        Dictionary mapping plot names to file paths
    """
    analysis_path = Path(analysis_dir)
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)
    
    generated = {}
    
    # RMSD
    rmsd_file = analysis_path / 'rmsd.xvg'
    if rmsd_file.exists():
        result = plot_rmsd(rmsd_file, plots_path / 'rmsd.png', dpi=dpi)
        if result:
            generated['rmsd'] = result
    
    # RMSF
    rmsf_file = analysis_path / 'rmsf.xvg'
    if rmsf_file.exists():
        result = plot_rmsf(rmsf_file, plots_path / 'rmsf.png', dpi=dpi)
        if result:
            generated['rmsf'] = result
    
    # Gyration
    gyrate_file = analysis_path / 'gyrate.xvg'
    if gyrate_file.exists():
        result = plot_gyration(gyrate_file, plots_path / 'gyrate.png', dpi=dpi)
        if result:
            generated['gyrate'] = result
    
    # H-bonds
    hbonds_file = analysis_path / 'hbonds.xvg'
    if hbonds_file.exists():
        result = plot_hbonds(hbonds_file, plots_path / 'hbonds.png', dpi=dpi)
        if result:
            generated['hbonds'] = result
    
    # Contact map
    contact_file = analysis_path / 'contact_map.txt'
    if contact_file.exists():
        result = plot_contact_map(contact_file, plots_path / 'contact_map.png', dpi=dpi)
        if result:
            generated['contact_map'] = result
    
    # Also export plot data as CSV in stats folder
    stats_dir = analysis_path.parent / 'stats'
    if stats_dir.exists() or analysis_path.parent.exists():
        stats_dir.mkdir(parents=True, exist_ok=True)
        export_plot_data_csv(analysis_path, stats_dir)
    
    return generated


def export_plot_data_csv(analysis_dir: Union[str, Path],
                         stats_dir: Union[str, Path]) -> Dict[str, Path]:
    """
    Export all XVG plot data as CSV files for later visualization.
    
    Args:
        analysis_dir: Directory containing XVG analysis files
        stats_dir: Directory to save CSV files
        
    Returns:
        Dictionary mapping plot names to CSV file paths
    """
    import csv
    from .xvg_parser import parse_xvg
    
    analysis_path = Path(analysis_dir)
    stats_path = Path(stats_dir)
    stats_path.mkdir(parents=True, exist_ok=True)
    
    exported = {}
    
    # Define XVG files to export with their column headers
    xvg_files = {
        'rmsd': {'file': 'rmsd.xvg', 'headers': ['time_ps', 'rmsd_nm']},
        'rmsf': {'file': 'rmsf.xvg', 'headers': ['residue', 'rmsf_nm']},
        'gyrate': {'file': 'gyrate.xvg', 'headers': ['time_ps', 'rg_nm', 'rgx_nm', 'rgy_nm', 'rgz_nm']},
        'sasa': {'file': 'sasa.xvg', 'headers': ['time_ps', 'sasa_nm2']},
        'hbonds': {'file': 'hbonds.xvg', 'headers': ['time_ps', 'num_hbonds']},
        'mindist': {'file': 'mindist.xvg', 'headers': ['time_ps', 'min_distance_nm']},
        'numcont': {'file': 'numcont.xvg', 'headers': ['time_ps', 'num_contacts']},
    }
    
    for name, info in xvg_files.items():
        xvg_file = analysis_path / info['file']
        if not xvg_file.exists():
            continue
        
        csv_file = stats_path / f"{name}_data.csv"
        
        try:
            data, metadata = parse_xvg(xvg_file)
            if not data:
                continue
            
            with open(csv_file, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # Use provided headers or generate from metadata
                headers = info['headers'][:len(data[0])] if len(info['headers']) >= len(data[0]) else \
                          info['headers'] + [f'value_{i}' for i in range(len(info['headers']), len(data[0]))]
                writer.writerow(headers)
                
                for row in data:
                    formatted_row = [f"{v:.6f}" if isinstance(v, float) else v for v in row]
                    writer.writerow(formatted_row)
            
            exported[name] = csv_file
        except Exception as e:
            print(f"  Warning: Could not export {name}: {e}")
    
    return exported


def generate_batch_plots(csv_file: Union[str, Path],
                         plots_dir: Union[str, Path],
                         dpi: int = 150) -> Dict[str, Path]:
    """
    Generate plots for batch comparison results.
    
    Args:
        csv_file: Path to comparison CSV file
        plots_dir: Directory to save plots
        dpi: Output resolution
    
    Returns:
        Dictionary mapping plot names to file paths
    """
    if not check_matplotlib():
        return {}
    
    import matplotlib.pyplot as plt
    
    try:
        import pandas as pd
    except ImportError:
        print("pandas not available, skipping batch plots")
        return {}
    
    csv_path = Path(csv_file)
    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)
    
    generated = {}
    
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return generated
    
    # Bar plot of potential energies
    if 'structure' in df.columns and 'potential_kJ_mol' in df.columns:
        data = dict(zip(df['structure'], df['potential_kJ_mol']))
        result = plot_comparison_barplot(
            data, plots_path / 'comparison_barplot.png',
            title='Structure Comparison - Potential Energy',
            xlabel='Potential Energy (kJ/mol)',
            dpi=dpi
        )
        if result:
            generated['barplot'] = result
    
    # Scatter plot: Rg vs SASA
    if all(col in df.columns for col in ['gyration_nm', 'sasa_nm2', 'potential_kJ_mol']):
        result = plot_scatter_comparison(
            df, 'gyration_nm', 'sasa_nm2', 'potential_kJ_mol',
            plots_path / 'rg_vs_sasa.png',
            xlabel='Radius of Gyration (nm)',
            ylabel='SASA (nm²)',
            colorbar_label='Potential Energy (kJ/mol)',
            title='Compactness vs Surface Area',
            dpi=dpi
        )
        if result:
            generated['rg_vs_sasa'] = result
    
    return generated


# Command-line interface
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate GROMACS analysis plots')
    parser.add_argument('type', choices=['rmsd', 'rmsf', 'gyrate', 'hbonds', 'contact', 'batch', 'all'],
                        help='Type of plot to generate')
    parser.add_argument('-i', '--input', required=True, help='Input file or directory')
    parser.add_argument('-o', '--output', required=True, help='Output file or directory')
    parser.add_argument('--dpi', type=int, default=150, help='Output DPI')
    
    args = parser.parse_args()
    
    if args.type == 'all':
        results = plot_all_md_analysis(args.input, args.output, args.dpi)
        for name, path in results.items():
            print(f"Created: {path}")
    elif args.type == 'batch':
        results = generate_batch_plots(args.input, args.output, args.dpi)
        for name, path in results.items():
            print(f"Created: {path}")
    else:
        funcs = {
            'rmsd': plot_rmsd,
            'rmsf': plot_rmsf,
            'gyrate': plot_gyration,
            'hbonds': plot_hbonds,
            'contact': plot_contact_map
        }
        result = funcs[args.type](args.input, args.output, dpi=args.dpi)
        if result:
            print(f"Created: {result}")
        else:
            print("Failed to create plot")
