#!/usr/bin/env python3
"""
Stability Analyzer Module for GROMACS MD Simulations

Analyzes stability metrics from GROMACS MD output files including
RMSD, radius of gyration, hydrogen bonds, SASA, and interaction energies.
"""

import os
import argparse
import numpy as np
from typing import Dict, Optional, Tuple, Any


def read_xvg(filename: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Read GROMACS XVG file and return time and data arrays.
    
    Args:
        filename: Path to the XVG file
        
    Returns:
        Tuple of (times array, values array) or (None, None) if file not found
    """
    if not os.path.exists(filename):
        return None, None
    
    times = []
    values = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    times.append(float(parts[0]))
                    values.append([float(x) for x in parts[1:]])
                except ValueError:
                    continue
    
    if not times:
        return None, None
        
    return np.array(times), np.array(values)


def analyze_stability(workdir: str, output_file: str = 'stability_metrics.txt') -> Dict[str, float]:
    """
    Analyze stability metrics from GROMACS output files.
    
    Args:
        workdir: Working directory containing GROMACS output files
        output_file: Name of output metrics file
        
    Returns:
        Dictionary containing all extracted stability metrics
    """
    original_dir = os.getcwd()
    os.chdir(workdir)
    
    try:
        results = {}
        
        # RMSD analysis
        times, rmsd = read_xvg('rmsd.xvg')
        if rmsd is not None and len(rmsd) > 0:
            rmsd_vals = rmsd[:, 0] if len(rmsd.shape) > 1 else rmsd
            results['rmsd_mean'] = float(np.mean(rmsd_vals))
            results['rmsd_std'] = float(np.std(rmsd_vals))
            results['rmsd_final'] = float(rmsd_vals[-1]) if len(rmsd_vals) > 0 else 0.0
        
        # Radius of gyration
        times, rg = read_xvg('gyrate.xvg')
        if rg is not None and len(rg) > 0:
            rg_vals = rg[:, 0] if len(rg.shape) > 1 else rg
            results['rg_mean'] = float(np.mean(rg_vals))
            results['rg_std'] = float(np.std(rg_vals))
        
        # Hydrogen bonds
        times, hb = read_xvg('hbonds.xvg')
        if hb is not None and len(hb) > 0:
            hb_vals = hb[:, 0] if len(hb.shape) > 1 else hb
            results['hbonds_mean'] = float(np.mean(hb_vals))
            results['hbonds_std'] = float(np.std(hb_vals))
        
        # Minimum distance
        times, mindist = read_xvg('mindist.xvg')
        if mindist is not None and len(mindist) > 0:
            mindist_vals = mindist[:, 0] if len(mindist.shape) > 1 else mindist
            results['mindist_mean'] = float(np.mean(mindist_vals))
            results['mindist_std'] = float(np.std(mindist_vals))
        
        # SASA
        times, sasa = read_xvg('sasa.xvg')
        if sasa is not None and len(sasa) > 0:
            sasa_vals = sasa[:, 0] if len(sasa.shape) > 1 else sasa
            results['sasa_mean'] = float(np.mean(sasa_vals))
            results['sasa_std'] = float(np.std(sasa_vals))
        
        # Interaction energy (if available)
        times, ie = read_xvg('interaction_energy.xvg')
        if ie is not None and len(ie) > 0:
            if len(ie.shape) > 1:
                coul = ie[:, 0]
                lj = ie[:, 1] if ie.shape[1] > 1 else np.zeros_like(coul)
            else:
                coul = ie
                lj = np.zeros_like(coul)
            
            results['coul_mean'] = float(np.mean(coul))
            results['lj_mean'] = float(np.mean(lj))
            results['total_ie_mean'] = float(np.mean(coul + lj))
        
        # Save results
        with open(output_file, 'w') as f:
            f.write("Stability Metrics Summary\n")
            f.write("=" * 50 + "\n\n")
            
            for key, value in sorted(results.items()):
                f.write(f"{key}: {value:.4f}\n")
        
        print(f"Stability metrics saved to {output_file}")
        return results
        
    finally:
        os.chdir(original_dir)


def main():
    """Command-line interface for stability analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze stability metrics from GROMACS MD output files'
    )
    parser.add_argument('--workdir', required=True,
                       help='Working directory with GROMACS outputs')
    parser.add_argument('--output', default='stability_metrics.txt',
                       help='Output metrics file name')
    
    args = parser.parse_args()
    
    analyze_stability(args.workdir, args.output)


if __name__ == '__main__':
    main()
