#!/usr/bin/env python3
"""
Metrics Extractor Module for GROMACS Stability Analysis

Extracts and compiles metrics from GROMACS output files (XVG, logs).
"""

import os
import re
import argparse
from typing import Optional, List, Dict, Any


def read_xvg_last(filename: str) -> Optional[List[float]]:
    """
    Read the last data value from an XVG file.
    
    Args:
        filename: Path to the XVG file
        
    Returns:
        List of float values from the last line, or None if not found
    """
    if not os.path.exists(filename):
        return None
    
    last_val = None
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#') and not line.startswith('@'):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        last_val = [float(x) for x in parts[1:]]
                    except ValueError:
                        pass
    return last_val


def extract_energy_from_log(log_file: str = 'mdrun_em.log') -> Dict[str, float]:
    """
    Extract final energy values from GROMACS mdrun log.
    
    Args:
        log_file: Path to the mdrun log file
        
    Returns:
        Dictionary with energy values
    """
    if not os.path.exists(log_file):
        return {}
    
    energies = {}
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Look for final potential energy
    match = re.search(r'Potential Energy\s+=\s+([-\d.e+]+)', content)
    if match:
        energies['potential_energy'] = float(match.group(1))
    
    return energies


def extract_metrics(workdir: str, output_file: str = 'metrics.txt') -> Dict[str, Any]:
    """
    Extract and compile all interface metrics from GROMACS outputs.
    
    Args:
        workdir: Working directory containing GROMACS output files
        output_file: Name of output metrics file
        
    Returns:
        Dictionary containing all extracted metrics
    """
    original_dir = os.getcwd()
    os.chdir(workdir)
    
    try:
        metrics = {}
        
        # Energy from minimization
        metrics.update(extract_energy_from_log())
        
        # Read XVG files
        hbonds = read_xvg_last('hbonds.xvg')
        if hbonds:
            metrics['hbonds'] = hbonds[0]
        
        mindist = read_xvg_last('mindist.xvg')
        if mindist:
            metrics['min_distance'] = mindist[0]
        
        numcont = read_xvg_last('numcont.xvg')
        if numcont:
            metrics['contacts_0.6nm'] = numcont[0]
        
        sasa_total = read_xvg_last('sasa.xvg')
        if sasa_total:
            metrics['sasa_total'] = sasa_total[0]
        
        sasa_a = read_xvg_last('sasa_chainA.xvg')
        if sasa_a:
            metrics['sasa_chainA'] = sasa_a[0]
        
        sasa_b = read_xvg_last('sasa_chainB.xvg')
        if sasa_b:
            metrics['sasa_chainB'] = sasa_b[0]
        
        # Calculate buried surface area (interface area)
        if sasa_a and sasa_b and sasa_total:
            bsa = (sasa_a[0] + sasa_b[0]) - sasa_total[0]
            metrics['buried_surface_area'] = bsa
        
        # Save metrics
        with open(output_file, 'w') as f:
            for key, value in sorted(metrics.items()):
                f.write(f"{key}: {value}\n")
                print(f"    {key}: {value:.2f}")
        
        print(f"  Metrics saved to {output_file}")
        
        return metrics
        
    finally:
        os.chdir(original_dir)


def main():
    """Command-line interface for metrics extraction."""
    parser = argparse.ArgumentParser(
        description='Extract metrics from GROMACS output files'
    )
    parser.add_argument('--workdir', required=True, 
                       help='Working directory with GROMACS outputs')
    parser.add_argument('--output', default='metrics.txt',
                       help='Output metrics file name')
    
    args = parser.parse_args()
    
    metrics = extract_metrics(args.workdir, args.output)
    
    return metrics


if __name__ == '__main__':
    main()
