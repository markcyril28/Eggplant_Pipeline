#!/usr/bin/env python3
"""
XVG Parser Module for GROMACS Analysis

Utilities for parsing GROMACS XVG output files.
"""

import os
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import numpy as np


def parse_xvg(filename: Path) -> Tuple[List[List[float]], Dict[str, str]]:
    """
    Parse XVG file and return all data points with metadata.
    
    Args:
        filename: Path to the XVG file
        
    Returns:
        Tuple of (data_list, metadata_dict)
        data_list: List of [time, value1, value2, ...] for each row
        metadata_dict: Contains 'title', 'xlabel', 'ylabel'
    """
    data = []
    metadata = {'title': '', 'xlabel': '', 'ylabel': '', 'legends': []}
    
    if not Path(filename).exists():
        return data, metadata
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@ title'):
                metadata['title'] = line.split('"')[1] if '"' in line else ''
            elif line.startswith('@ xaxis label'):
                metadata['xlabel'] = line.split('"')[1] if '"' in line else ''
            elif line.startswith('@ yaxis label'):
                metadata['ylabel'] = line.split('"')[1] if '"' in line else ''
            elif line.startswith('@ s') and 'legend' in line:
                legend = line.split('"')[1] if '"' in line else ''
                metadata['legends'].append(legend)
            elif not line.startswith(('#', '@')) and line:
                parts = line.split()
                if parts:
                    try:
                        data.append([float(x) for x in parts])
                    except ValueError:
                        continue
    
    return data, metadata


def parse_xvg_last(filename: Path) -> Optional[List[float]]:
    """
    Read the last data value from an XVG file.
    
    Args:
        filename: Path to the XVG file
        
    Returns:
        List of float values from the last line, or None if not found
    """
    if not Path(filename).exists():
        return None
    
    last_val = None
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('@'):
                parts = line.split()
                if len(parts) >= 1:
                    try:
                        last_val = [float(x) for x in parts]
                    except ValueError:
                        pass
    return last_val


def parse_xvg_column(filename: Path, column: int = 1) -> Tuple[List[float], List[float]]:
    """
    Extract a specific column from XVG file.
    
    Args:
        filename: Path to the XVG file
        column: Column index (0=time, 1=first data column, etc.)
        
    Returns:
        Tuple of (times, values)
    """
    data, _ = parse_xvg(filename)
    if not data:
        return [], []
    
    times = [row[0] for row in data if len(row) > column]
    values = [row[column] for row in data if len(row) > column]
    return times, values


def calculate_statistics(values: List[float]) -> Dict[str, float]:
    """
    Calculate basic statistics for a list of values.
    
    Args:
        values: List of numeric values
        
    Returns:
        Dictionary with mean, std, min, max, median
    """
    if not values:
        return {'mean': 0, 'std': 0, 'min': 0, 'max': 0, 'median': 0}
    
    arr = np.array(values)
    return {
        'mean': float(np.mean(arr)),
        'std': float(np.std(arr)),
        'min': float(np.min(arr)),
        'max': float(np.max(arr)),
        'median': float(np.median(arr))
    }


def parse_energy_xvg(filename: Path) -> Dict[str, float]:
    """
    Parse energy XVG file and return named energy components.
    
    Args:
        filename: Path to energy XVG file
        
    Returns:
        Dictionary mapping energy names to values
    """
    data, metadata = parse_xvg(filename)
    
    if not data:
        return {}
    
    # Use the last row of data
    last_row = data[-1]
    
    # Map legends to values
    energies = {}
    legends = metadata.get('legends', [])
    
    for i, legend in enumerate(legends):
        col_idx = i + 1  # Skip time column
        if col_idx < len(last_row):
            # Clean up legend name
            key = legend.lower().replace(' ', '_').replace('-', '_')
            energies[key] = last_row[col_idx]
    
    # Also include unnamed columns
    for i in range(1, len(last_row)):
        if i > len(legends):
            energies[f'value_{i}'] = last_row[i]
    
    return energies


def xvg_to_csv(xvg_file: Path, csv_file: Path, 
               headers: Optional[List[str]] = None) -> None:
    """
    Convert XVG file to CSV format.
    
    Args:
        xvg_file: Input XVG file path
        csv_file: Output CSV file path
        headers: Optional column headers
    """
    data, metadata = parse_xvg(xvg_file)
    
    if not data:
        return
    
    # Generate headers if not provided
    if headers is None:
        headers = ['time']
        if metadata['legends']:
            headers.extend(metadata['legends'])
        else:
            headers.extend([f'value_{i}' for i in range(1, len(data[0]))])
    
    with open(csv_file, 'w') as f:
        f.write(','.join(headers) + '\n')
        for row in data:
            f.write(','.join(str(v) for v in row) + '\n')


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Parse GROMACS XVG files")
    parser.add_argument("xvg_file", help="XVG file to parse")
    parser.add_argument("--to-csv", help="Convert to CSV file")
    parser.add_argument("--stats", action="store_true", help="Print statistics")
    args = parser.parse_args()
    
    data, metadata = parse_xvg(args.xvg_file)
    
    print(f"Title: {metadata['title']}")
    print(f"X-axis: {metadata['xlabel']}")
    print(f"Y-axis: {metadata['ylabel']}")
    print(f"Legends: {metadata['legends']}")
    print(f"Data rows: {len(data)}")
    
    if data:
        print(f"First row: {data[0]}")
        print(f"Last row: {data[-1]}")
    
    if args.stats and data:
        for i in range(1, len(data[0])):
            values = [row[i] for row in data if len(row) > i]
            stats = calculate_statistics(values)
            label = metadata['legends'][i-1] if i-1 < len(metadata['legends']) else f'Column {i}'
            print(f"\n{label}:")
            for k, v in stats.items():
                print(f"  {k}: {v:.4f}")
    
    if args.to_csv:
        xvg_to_csv(args.xvg_file, args.to_csv)
        print(f"\nConverted to: {args.to_csv}")
