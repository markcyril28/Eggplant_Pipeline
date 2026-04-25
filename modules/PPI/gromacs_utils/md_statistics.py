#!/usr/bin/env python3
"""
MD Statistics Module for GROMACS Analysis

Extracts and calculates statistics from MD simulation outputs.
Generates comprehensive CSV data files in 'stats' folder for later visualization.
"""

import os
import json
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
import numpy as np

from .xvg_parser import parse_xvg, calculate_statistics, xvg_to_csv


@dataclass
class MDStatistics:
    """Container for MD simulation statistics."""
    
    structure: str = ""
    simulation_type: str = "MD"
    timestep_fs: float = 2.0
    temperature_K: float = 300.0
    duration_ps: float = 0.0
    
    # Structure metrics
    rmsd: Dict[str, float] = field(default_factory=dict)
    rmsf: Dict[str, float] = field(default_factory=dict)
    gyration: Dict[str, float] = field(default_factory=dict)
    sasa: Dict[str, float] = field(default_factory=dict)
    
    # Interface metrics
    hbonds: Dict[str, float] = field(default_factory=dict)
    min_distance: Dict[str, float] = field(default_factory=dict)
    contacts: Dict[str, float] = field(default_factory=dict)
    
    # Energy metrics
    interaction_energy: Dict[str, float] = field(default_factory=dict)
    potential_energy: Dict[str, float] = field(default_factory=dict)
    
    # Flexible residues
    flexible_residues: List[int] = field(default_factory=list)
    
    # Raw timeseries (for plotting and CSV export)
    timeseries: Dict[str, Dict] = field(default_factory=dict)
    
    # RMSF per residue data
    rmsf_per_residue: List[Tuple[int, float]] = field(default_factory=list)


def extract_md_statistics(workdir: Path) -> MDStatistics:
    """
    Extract MD statistics from GROMACS output files.
    
    Args:
        workdir: Working directory containing XVG files
        
    Returns:
        MDStatistics object
    """
    stats = MDStatistics(structure=workdir.name)
    workdir = Path(workdir)
    
    # Check for analysis subfolder (used by compare script)
    analysis_dir = workdir / 'analysis'
    if not analysis_dir.exists():
        analysis_dir = workdir  # Fallback to workdir itself
    
    # Parse RMSD
    rmsd_file = analysis_dir / 'rmsd.xvg'
    if not rmsd_file.exists():
        rmsd_file = workdir / 'rmsd.xvg'
    rmsd_data, _ = parse_xvg(rmsd_file)
    if rmsd_data:
        times = [d[0] for d in rmsd_data]
        values = [d[1] for d in rmsd_data]
        stats.rmsd = calculate_statistics(values)
        stats.timeseries['rmsd'] = {'time_ps': times, 'value_nm': values}
        stats.duration_ps = times[-1] if times else 0
    
    # Parse radius of gyration
    gyrate_file = analysis_dir / 'gyrate.xvg'
    if not gyrate_file.exists():
        gyrate_file = workdir / 'gyrate.xvg'
    gyrate_data, _ = parse_xvg(gyrate_file)
    if gyrate_data:
        times = [d[0] for d in gyrate_data]
        values = [d[1] for d in gyrate_data]
        stats.gyration = calculate_statistics(values)
        stats.timeseries['gyration'] = {'time_ps': times, 'value_nm': values}
    
    # Parse RMSF
    rmsf_file = analysis_dir / 'rmsf.xvg'
    if not rmsf_file.exists():
        rmsf_file = workdir / 'rmsf.xvg'
    rmsf_data, _ = parse_xvg(rmsf_file)
    if rmsf_data:
        residues = [int(d[0]) for d in rmsf_data]
        values = [d[1] for d in rmsf_data]
        stats.rmsf = calculate_statistics(values)
        stats.rmsf_per_residue = list(zip(residues, values))
        
        # Find flexible regions (top 10% RMSF)
        if values:
            threshold = np.percentile(values, 90)
            stats.flexible_residues = [r for r, v in zip(residues, values) if v > threshold][:20]
    
    # Parse SASA
    sasa_file = analysis_dir / 'sasa.xvg'
    if not sasa_file.exists():
        sasa_file = workdir / 'sasa.xvg'
    sasa_data, _ = parse_xvg(sasa_file)
    if sasa_data:
        times = [d[0] for d in sasa_data]
        values = [d[1] for d in sasa_data if len(d) > 1]
        if values:
            stats.sasa = calculate_statistics(values)
            stats.timeseries['sasa'] = {'time_ps': times[:len(values)], 'value_nm2': values}
    
    # Parse H-bonds
    hbonds_file = analysis_dir / 'hbonds.xvg'
    if not hbonds_file.exists():
        hbonds_file = workdir / 'hbonds.xvg'
    hbonds_data, _ = parse_xvg(hbonds_file)
    if hbonds_data:
        times = [d[0] for d in hbonds_data]
        values = [d[1] for d in hbonds_data if len(d) > 1]
        if values:
            stats.hbonds = calculate_statistics(values)
            stats.timeseries['hbonds'] = {'time_ps': times[:len(values)], 'count': values}
    
    # Parse minimum distance
    mindist_file = analysis_dir / 'mindist.xvg'
    if not mindist_file.exists():
        mindist_file = workdir / 'mindist.xvg'
    mindist_data, _ = parse_xvg(mindist_file)
    if mindist_data:
        times = [d[0] for d in mindist_data]
        values = [d[1] for d in mindist_data if len(d) > 1]
        if values:
            stats.min_distance = calculate_statistics(values)
            stats.timeseries['mindist'] = {'time_ps': times[:len(values)], 'value_nm': values}
    
    # Parse contacts
    contacts_file = analysis_dir / 'numcont.xvg'
    if not contacts_file.exists():
        contacts_file = workdir / 'numcont.xvg'
    contacts_data, _ = parse_xvg(contacts_file)
    if contacts_data:
        times = [d[0] for d in contacts_data]
        values = [d[1] for d in contacts_data if len(d) > 1]
        if values:
            stats.contacts = calculate_statistics(values)
            stats.timeseries['contacts'] = {'time_ps': times[:len(values)], 'count': values}
    
    # Parse interaction energy (if available)
    ie_file = workdir / 'interaction_energy.xvg'
    if ie_file.exists():
        ie_data, ie_meta = parse_xvg(ie_file)
        if ie_data:
            legends = ie_meta.get('legends', [])
            coul_vals = []
            lj_vals = []
            for i, leg in enumerate(legends):
                col_idx = i + 1
                if 'Coul' in leg:
                    coul_vals = [d[col_idx] for d in ie_data if len(d) > col_idx]
                elif 'LJ' in leg:
                    lj_vals = [d[col_idx] for d in ie_data if len(d) > col_idx]
            
            if coul_vals:
                stats.interaction_energy['coulomb'] = calculate_statistics(coul_vals)
            if lj_vals:
                stats.interaction_energy['lj'] = calculate_statistics(lj_vals)
            if coul_vals and lj_vals:
                total = [c + l for c, l in zip(coul_vals, lj_vals)]
                stats.interaction_energy['total'] = calculate_statistics(total)
    
    # Parse potential energy from em_potential.xvg
    em_file = workdir / 'em_potential.xvg'
    if em_file.exists():
        em_data, _ = parse_xvg(em_file)
        if em_data:
            values = [d[1] for d in em_data if len(d) > 1]
            if values:
                stats.potential_energy = calculate_statistics(values)
    
    return stats


def save_statistics_json(stats: MDStatistics, output_file: Path,
                         include_timeseries: bool = False) -> None:
    """
    Save MD statistics to JSON file.
    
    Args:
        stats: MDStatistics object
        output_file: Path to output file
        include_timeseries: Whether to include raw timeseries data
    """
    data = {
        'structure': stats.structure,
        'generated_at': datetime.now().isoformat(),
        'simulation': {
            'type': stats.simulation_type,
            'timestep_fs': stats.timestep_fs,
            'temperature_K': stats.temperature_K,
            'duration_ps': stats.duration_ps
        },
        'structure_metrics': {
            'rmsd_nm': stats.rmsd,
            'rmsf_nm': stats.rmsf,
            'gyration_nm': stats.gyration,
            'sasa_nm2': stats.sasa,
            'flexible_residues': stats.flexible_residues
        },
        'interface': {
            'hbonds': stats.hbonds,
            'min_distance_nm': stats.min_distance,
            'contacts': stats.contacts
        },
        'energy': {
            'interaction_energy': stats.interaction_energy,
            'potential_energy': stats.potential_energy
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    # Save timeseries separately if requested
    if include_timeseries and stats.timeseries:
        ts_file = output_file.parent / 'timeseries_data.json'
        with open(ts_file, 'w') as f:
            json.dump(stats.timeseries, f)


def save_statistics_csv(stats: MDStatistics, output_file: Path) -> None:
    """
    Save MD statistics summary to CSV file.
    
    Args:
        stats: MDStatistics object
        output_file: Path to output file
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'mean', 'std', 'min', 'max', 'median', 'unit', 'category'])
        
        # Structure metrics
        for name, data, unit in [
            ('rmsd', stats.rmsd, 'nm'),
            ('rmsf', stats.rmsf, 'nm'),
            ('gyration', stats.gyration, 'nm'),
            ('sasa', stats.sasa, 'nm²')
        ]:
            if data and 'mean' in data:
                writer.writerow([
                    name,
                    f"{data['mean']:.6f}",
                    f"{data['std']:.6f}",
                    f"{data['min']:.6f}",
                    f"{data['max']:.6f}",
                    f"{data.get('median', data['mean']):.6f}",
                    unit,
                    'structure'
                ])
        
        # Interface metrics
        for name, data, unit in [
            ('hbonds', stats.hbonds, 'count'),
            ('min_distance', stats.min_distance, 'nm'),
            ('contacts', stats.contacts, 'count')
        ]:
            if data and 'mean' in data:
                writer.writerow([
                    name,
                    f"{data['mean']:.6f}",
                    f"{data['std']:.6f}",
                    f"{data['min']:.6f}",
                    f"{data['max']:.6f}",
                    f"{data.get('median', data['mean']):.6f}",
                    unit,
                    'interface'
                ])
        
        # Energy metrics - interaction energy
        if stats.interaction_energy:
            for name, data in stats.interaction_energy.items():
                if data and 'mean' in data:
                    writer.writerow([
                        f'ie_{name}',
                        f"{data['mean']:.6f}",
                        f"{data['std']:.6f}",
                        f"{data['min']:.6f}",
                        f"{data['max']:.6f}",
                        f"{data.get('median', data['mean']):.6f}",
                        'kJ/mol',
                        'energy'
                    ])
        
        # Potential energy
        if stats.potential_energy and 'mean' in stats.potential_energy:
            writer.writerow([
                'potential_energy',
                f"{stats.potential_energy['mean']:.6f}",
                f"{stats.potential_energy['std']:.6f}",
                f"{stats.potential_energy['min']:.6f}",
                f"{stats.potential_energy['max']:.6f}",
                f"{stats.potential_energy.get('median', stats.potential_energy['mean']):.6f}",
                'kJ/mol',
                'energy'
            ])


def save_timeseries_csv(stats: MDStatistics, output_dir: Path) -> Dict[str, Path]:
    """
    Save all timeseries data as individual CSV files.
    
    Args:
        stats: MDStatistics object
        output_dir: Directory to save CSV files
        
    Returns:
        Dictionary mapping timeseries name to CSV file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = {}
    
    for name, data in stats.timeseries.items():
        if not data:
            continue
        
        csv_file = output_dir / f'{name}_timeseries.csv'
        
        # Get column names and values
        columns = list(data.keys())
        if not columns:
            continue
        
        # Get the length from the first column
        first_col = columns[0]
        n_rows = len(data[first_col])
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            
            for i in range(n_rows):
                row = []
                for col in columns:
                    vals = data[col]
                    if i < len(vals):
                        val = vals[i]
                        if isinstance(val, float):
                            row.append(f"{val:.6f}")
                        else:
                            row.append(val)
                    else:
                        row.append('')
                writer.writerow(row)
        
        saved_files[name] = csv_file
    
    return saved_files


def save_rmsf_per_residue_csv(stats: MDStatistics, output_file: Path) -> Optional[Path]:
    """
    Save RMSF per residue data to CSV.
    
    Args:
        stats: MDStatistics object
        output_file: Path to output CSV file
        
    Returns:
        Path to saved file or None if no data
    """
    if not stats.rmsf_per_residue:
        return None
    
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['residue_number', 'rmsf_nm', 'is_flexible'])
        
        # Calculate threshold for flexible residues
        values = [v for _, v in stats.rmsf_per_residue]
        threshold = np.percentile(values, 90) if values else 0
        
        for residue, rmsf in stats.rmsf_per_residue:
            is_flexible = 'yes' if rmsf > threshold else 'no'
            writer.writerow([residue, f"{rmsf:.6f}", is_flexible])
    
    return output_file


def convert_analysis_xvg_to_csv(workdir: Path, stats_dir: Path) -> Dict[str, Path]:
    """
    Convert all XVG files in analysis folder to CSV format.
    
    Args:
        workdir: Working directory containing XVG files
        stats_dir: Output directory for CSV files
        
    Returns:
        Dictionary mapping XVG file name to CSV path
    """
    workdir = Path(workdir)
    stats_dir = Path(stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    converted = {}
    
    # Check both analysis subfolder and workdir
    analysis_dir = workdir / 'analysis'
    search_dirs = [analysis_dir, workdir] if analysis_dir.exists() else [workdir]
    
    xvg_files = set()
    for search_dir in search_dirs:
        xvg_files.update(search_dir.glob('*.xvg'))
    
    for xvg_file in xvg_files:
        csv_file = stats_dir / f'{xvg_file.stem}.csv'
        
        # Skip if already converted
        if csv_file.exists():
            converted[xvg_file.stem] = csv_file
            continue
        
        try:
            xvg_to_csv(xvg_file, csv_file)
            converted[xvg_file.stem] = csv_file
        except Exception as e:
            print(f"  Warning: Could not convert {xvg_file.name}: {e}")
    
    return converted


def generate_md_statistics(workdir: Path) -> Dict[str, Path]:
    """
    Generate all MD statistics files including comprehensive CSV data.
    
    Creates:
    - statistics/md_statistics.json - Summary statistics
    - statistics/md_statistics.csv - Summary statistics in CSV
    - statistics/timeseries_data.json - Raw timeseries as JSON
    - stats/ - Folder with all raw data in CSV format:
        - stats/summary_statistics.csv - Copy of summary stats
        - stats/rmsf_per_residue.csv - RMSF for each residue
        - stats/*_timeseries.csv - Timeseries for each metric
        - stats/*.csv - Converted XVG files
    
    Args:
        workdir: Working directory
        
    Returns:
        Dictionary mapping file type to path
    """
    workdir = Path(workdir)
    stats_dir = workdir / 'statistics'
    stats_dir.mkdir(exist_ok=True)
    
    # Create 'stats' folder for comprehensive CSV data
    csv_stats_dir = workdir / 'stats'
    csv_stats_dir.mkdir(exist_ok=True)
    
    # Extract statistics
    stats = extract_md_statistics(workdir)
    
    # Save files
    files = {}
    
    # 1. Main JSON file
    json_file = stats_dir / 'md_statistics.json'
    save_statistics_json(stats, json_file, include_timeseries=True)
    files['json'] = json_file
    print(f"Created: {json_file}")
    
    # 2. Summary CSV in statistics folder
    csv_file = stats_dir / 'md_statistics.csv'
    save_statistics_csv(stats, csv_file)
    files['csv'] = csv_file
    print(f"Created: {csv_file}")
    
    # 3. Copy summary CSV to stats folder
    summary_csv = csv_stats_dir / 'summary_statistics.csv'
    save_statistics_csv(stats, summary_csv)
    files['stats_summary'] = summary_csv
    print(f"Created: {summary_csv}")
    
    # 4. Save timeseries data as individual CSVs
    ts_files = save_timeseries_csv(stats, csv_stats_dir)
    for name, path in ts_files.items():
        files[f'ts_{name}'] = path
        print(f"Created: {path}")
    
    # 5. Save RMSF per residue
    rmsf_file = save_rmsf_per_residue_csv(stats, csv_stats_dir / 'rmsf_per_residue.csv')
    if rmsf_file:
        files['rmsf_per_residue'] = rmsf_file
        print(f"Created: {rmsf_file}")
    
    # 6. Convert all XVG files to CSV
    xvg_csvs = convert_analysis_xvg_to_csv(workdir, csv_stats_dir)
    for name, path in xvg_csvs.items():
        files[f'xvg_{name}'] = path
        print(f"Created: {path}")
    
    # 7. Save timeseries JSON
    if stats.timeseries:
        ts_json = stats_dir / 'timeseries_data.json'
        files['timeseries_json'] = ts_json
        print(f"Created: {ts_json}")
    
    # 8. Create a metadata file with generation info
    meta_file = csv_stats_dir / 'metadata.csv'
    with open(meta_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['key', 'value'])
        writer.writerow(['structure', stats.structure])
        writer.writerow(['generated_at', datetime.now().isoformat()])
        writer.writerow(['simulation_type', stats.simulation_type])
        writer.writerow(['duration_ps', stats.duration_ps])
        writer.writerow(['temperature_K', stats.temperature_K])
        writer.writerow(['num_flexible_residues', len(stats.flexible_residues)])
        writer.writerow(['flexible_residues', ','.join(map(str, stats.flexible_residues))])
    files['metadata'] = meta_file
    print(f"Created: {meta_file}")
    
    return files


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate MD statistics")
    parser.add_argument("--workdir", required=True, help="Working directory")
    args = parser.parse_args()
    
    files = generate_md_statistics(Path(args.workdir))
    
    print("\nStatistics files generated:")
    for name, path in files.items():
        print(f"  {name}: {path}")
