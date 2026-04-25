#!/usr/bin/env python3
"""
GROMACS CLI Tools - Command-line interface for gromacs_utils modules.

This module provides CLI wrappers for common operations, eliminating
the need for inline Python in bash scripts.

Usage from bash:
    python3 -m gromacs_utils.cli <command> [args]
    
Or via the wrapper script:
    ${MODULES_DIR}/gromacs_cli.py <command> [args]
"""

import sys
import os
import json
import argparse
from pathlib import Path


def cmd_prepare_structure(args):
    """Prepare a structure file for GROMACS (convert CIF, clean PDB)."""
    from .structure_converter import prepare_structure
    result = prepare_structure(args.input, args.output)
    print(f"Prepared: {result}")
    return 0


def cmd_generate_mdp(args):
    """Generate MDP file(s) for GROMACS simulations."""
    from .mdp_generator import (
        MDPConfig, generate_em_mdp, generate_nvt_mdp, 
        generate_npt_mdp, generate_md_mdp, generate_ie_mdp,
        generate_all_mdp_files
    )
    
    config = MDPConfig(
        em_steps=args.em_steps,
        em_tolerance=args.em_tolerance,
        nvt_steps=args.nvt_steps,
        npt_steps=args.npt_steps,
        md_steps=args.md_steps,
        temperature=args.temperature
    )
    
    if args.type == 'all':
        files = generate_all_mdp_files(args.output, config)
        for name, path in files.items():
            print(f"Generated: {path}")
    else:
        generators = {
            'em': generate_em_mdp,
            'nvt': generate_nvt_mdp,
            'npt': generate_npt_mdp,
            'md': generate_md_mdp,
            'ie': generate_ie_mdp
        }
        generators[args.type](args.output, config)
        print(f"Generated: {args.output}")
    
    return 0


def cmd_extract_metrics(args):
    """Extract and save metrics from XVG files."""
    from .xvg_parser import parse_xvg_last
    
    workdir = Path(args.workdir)
    
    metrics = {
        'potential': None,
        'gyration': None,
        'sasa': None,
        'hbonds': None,
        'contacts': None,
        'min_distance': None,
        'status': 'OK'
    }
    
    # Parse energy
    energy_file = workdir / 'energy.xvg'
    if not energy_file.exists():
        energy_file = workdir / 'energies.xvg'
    if energy_file.exists():
        energy = parse_xvg_last(energy_file)
        if energy:
            metrics['potential'] = energy[-1] if energy else None
    
    # Parse gyration
    gyrate_file = workdir / 'gyrate.xvg'
    if gyrate_file.exists():
        gyrate = parse_xvg_last(gyrate_file)
        if gyrate:
            metrics['gyration'] = gyrate[-1] if len(gyrate) > 1 else gyrate[0]
    
    # Parse SASA
    sasa_file = workdir / 'sasa.xvg'
    if sasa_file.exists():
        sasa = parse_xvg_last(sasa_file)
        if sasa:
            metrics['sasa'] = sasa[-1] if len(sasa) > 1 else sasa[0]
    
    # Parse H-bonds
    hbonds_file = workdir / 'hbonds.xvg'
    if hbonds_file.exists():
        hbonds = parse_xvg_last(hbonds_file)
        if hbonds:
            metrics['hbonds'] = int(hbonds[-1] if len(hbonds) > 1 else hbonds[0])
    
    # Parse contacts
    numcont_file = workdir / 'numcont.xvg'
    if numcont_file.exists():
        numcont = parse_xvg_last(numcont_file)
        if numcont:
            metrics['contacts'] = int(numcont[-1] if len(numcont) > 1 else numcont[0])
    
    # Parse min distance
    mindist_file = workdir / 'mindist.xvg'
    if mindist_file.exists():
        mindist = parse_xvg_last(mindist_file)
        if mindist:
            metrics['min_distance'] = mindist[-1] if len(mindist) > 1 else mindist[0]
    
    # Determine output file (default to metrics.json for JSON)
    output_file = workdir / args.output if args.output else workdir / 'metrics.json'
    
    # Check if user wants JSON or text format
    if str(output_file).endswith('.json'):
        # Save JSON
        with open(output_file, 'w') as f:
            json.dump(metrics, f, indent=2)
        print(f"Saved JSON metrics to: {output_file}")
    else:
        # Save text summary
        with open(output_file, 'w') as f:
            f.write("Structure Metrics\n")
            f.write("=" * 40 + "\n")
            for key, value in metrics.items():
                if value is not None and key != 'status':
                    if isinstance(value, float):
                        f.write(f"{key}: {value:.4f}\n")
                    else:
                        f.write(f"{key}: {value}\n")
        print(f"Saved text metrics to: {output_file}")
        
        # Also save JSON version alongside text
        json_file = output_file.parent / 'metrics.json'
        with open(json_file, 'w') as f:
            json.dump(metrics, f, indent=2)
    
    # Print summary
    if metrics['potential'] is not None:
        print(f"Potential: {metrics['potential']:.1f} kJ/mol")
    
    return 0


def cmd_run_batch_analysis(args):
    """Run batch analysis on structure directories."""
    from .batch_analyzer import run_batch_analysis
    from .plotting import generate_batch_plots
    
    output_dir = Path(args.workdir)
    
    # Run analysis
    run_batch_analysis(output_dir)
    
    # Generate plots if requested
    if args.plots:
        plots_dir = output_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        csv_file = output_dir / 'comparison_data.csv'
        if csv_file.exists():
            generate_batch_plots(csv_file, plots_dir)
    
    print(f"Analysis complete: {output_dir}")
    return 0


def cmd_generate_plots(args):
    """Generate plots from analysis files."""
    from .plotting import plot_all_md_analysis, generate_batch_plots, plot_contact_map
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if args.type == 'contact':
        # Single contact map plot
        result = plot_contact_map(input_path, output_path, dpi=args.dpi)
        if result:
            print(f"Created: {result}")
        results = {}
    else:
        output_path.mkdir(parents=True, exist_ok=True)
        
        if args.type == 'md':
            results = plot_all_md_analysis(input_path, output_path, dpi=args.dpi)
        elif args.type == 'batch':
            results = generate_batch_plots(input_path, output_path, dpi=args.dpi)
        else:
            results = {}
        
        for name, path in results.items():
            print(f"Created: {path}")
    
    return 0


def cmd_setup_output(args):
    """Create standard output directory structure."""
    from .output_organizer import OutputOrganizer
    
    organizer = OutputOrganizer(Path(args.workdir))
    dirs = organizer.create_structure()
    
    for name, path in dirs.items():
        print(f"Created: {path}")
    
    return 0


def cmd_chain_index(args):
    """Create chain index file for GROMACS."""
    from .chain_parser import parse_chains_for_index
    from .chain_indexer import create_chain_index
    
    # Use the chain_indexer module
    create_chain_index(
        gro_file=args.gro,
        pdb_file=args.pdb,
        index_file=args.index
    )
    
    print(f"Updated index: {args.index}")
    return 0


def cmd_chain_info(args):
    """Get chain information from a PDB file."""
    from .chain_parser import get_chain_info
    
    outdir = args.outdir if args.outdir else os.path.dirname(args.pdb) or '.'
    
    info = get_chain_info(args.pdb, outdir)
    
    # Print summary
    print(f"Chains found: {len(info.get('atoms', {}))}")
    for chain, count in info.get('atoms', {}).items():
        residues = info.get('residues', {}).get(chain, set())
        print(f"  Chain {chain}: {count} atoms, {len(residues)} residues")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='GROMACS CLI Tools',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # prepare-structure
    p = subparsers.add_parser('prepare-structure', help='Prepare structure for GROMACS')
    p.add_argument('input', help='Input structure file (PDB or CIF)')
    p.add_argument('-o', '--output', default='clean.pdb', help='Output PDB file')
    p.set_defaults(func=cmd_prepare_structure)
    
    # generate-mdp
    p = subparsers.add_parser('generate-mdp', help='Generate MDP files')
    p.add_argument('type', choices=['em', 'nvt', 'npt', 'md', 'ie', 'all'],
                   help='Type of MDP file')
    p.add_argument('-o', '--output', default='.', help='Output file or directory')
    p.add_argument('--em-steps', type=int, default=5000)
    p.add_argument('--em-tolerance', type=float, default=100.0)
    p.add_argument('--nvt-steps', type=int, default=50000)
    p.add_argument('--npt-steps', type=int, default=50000)
    p.add_argument('--md-steps', type=int, default=250000)
    p.add_argument('--temperature', type=float, default=300.0)
    p.set_defaults(func=cmd_generate_mdp)
    
    # extract-metrics
    p = subparsers.add_parser('extract-metrics', help='Extract metrics from XVG files')
    p.add_argument('--workdir', default='.', help='Working directory with XVG files')
    p.add_argument('-o', '--output', default='metrics.txt', help='Output text file')
    p.set_defaults(func=cmd_extract_metrics)
    
    # batch-analysis
    p = subparsers.add_parser('batch-analysis', help='Run batch structure analysis')
    p.add_argument('--workdir', default='.', help='Directory with structure subdirs')
    p.add_argument('--plots', action='store_true', help='Generate plots')
    p.set_defaults(func=cmd_run_batch_analysis)
    
    # generate-plots
    p = subparsers.add_parser('generate-plots', help='Generate analysis plots')
    p.add_argument('type', choices=['md', 'batch', 'contact'], help='Type of plots')
    p.add_argument('-i', '--input', required=True, help='Input directory or file')
    p.add_argument('-o', '--output', required=True, help='Output directory or file')
    p.add_argument('--dpi', type=int, default=150, help='Plot DPI')
    p.set_defaults(func=cmd_generate_plots)
    
    # setup-output
    p = subparsers.add_parser('setup-output', help='Create output directory structure')
    p.add_argument('--workdir', default='.', help='Base directory')
    p.set_defaults(func=cmd_setup_output)
    
    # chain-index
    p = subparsers.add_parser('chain-index', help='Create/update chain index')
    p.add_argument('--pdb', required=True, help='Input PDB file')
    p.add_argument('--gro', required=True, help='Input GRO file')
    p.add_argument('--index', required=True, help='Index file to update')
    p.set_defaults(func=cmd_chain_index)
    
    # chain-info
    p = subparsers.add_parser('chain-info', help='Get chain information')
    p.add_argument('pdb', help='Input PDB file')
    p.add_argument('--outdir', help='Output directory for chain_info.json')
    p.set_defaults(func=cmd_chain_info)
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return 1
    
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
