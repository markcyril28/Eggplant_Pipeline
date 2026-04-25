#!/usr/bin/env python3
"""
Extract critical residues from MutateX results.
Identifies destabilizing and stabilizing mutations based on ΔΔG thresholds.

Usage:
    python extract_critical_residues.py <results_dir> <pdb_file> [options]
    
Example:
    python extract_critical_residues.py RESULTS/results_mutatex/SmelDMP10.200_SmelHAP2 inputs/SmelDMP/SmelDMP10.200_SmelHAP2.pdb
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from mutatex_parser import parse_selfmutation_dat_file


def find_and_parse_data(results_dir):
    """Find and parse all MutateX output data."""
    interface_dfs = []
    folding_dfs = []
    
    # Check for selfmutation results structure
    selfmut_dir = os.path.join(results_dir, "results_selfmutation")
    
    if os.path.isdir(selfmut_dir):
        # Interface DDGs
        interface_ddg_dir = os.path.join(selfmut_dir, "interface_ddgs")
        if os.path.isdir(interface_ddg_dir):
            for model_dir in os.listdir(interface_ddg_dir):
                model_path = os.path.join(interface_ddg_dir, model_dir)
                if os.path.isdir(model_path):
                    for f in os.listdir(model_path):
                        if f.endswith('.dat'):
                            df = parse_selfmutation_dat_file(os.path.join(model_path, f), 'interface')
                            if not df.empty:
                                interface_dfs.append(df)
        
        # Folding/Mutation DDGs
        mutation_ddg_dir = os.path.join(selfmut_dir, "mutation_ddgs")
        if os.path.isdir(mutation_ddg_dir):
            for model_dir in os.listdir(mutation_ddg_dir):
                model_path = os.path.join(mutation_ddg_dir, model_dir)
                if os.path.isdir(model_path):
                    for f in os.listdir(model_path):
                        if f.endswith('.dat'):
                            df = parse_selfmutation_dat_file(os.path.join(model_path, f), 'folding')
                            if not df.empty:
                                folding_dfs.append(df)
    
    # Also check for full mutation scan results
    interface_dir = os.path.join(results_dir, "results", "interface_ddgs", "final_averages")
    if os.path.isdir(interface_dir):
        for chain_file in os.listdir(interface_dir):
            if chain_file.endswith('.txt') or chain_file.endswith('.dat'):
                # Parse full mutation format if needed
                pass
    
    interface_df = pd.concat(interface_dfs, ignore_index=True) if interface_dfs else pd.DataFrame()
    folding_df = pd.concat(folding_dfs, ignore_index=True) if folding_dfs else pd.DataFrame()
    
    return interface_df, folding_df


def extract_critical_residues(df, destab_threshold, stab_threshold, top_n, data_type):
    """Extract critical residues based on thresholds."""
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    # Destabilizing residues (positive ΔΔG above threshold)
    destabilizing = df[df['ddg'] > destab_threshold].copy()
    destabilizing = destabilizing.sort_values('ddg', ascending=False)
    
    # Stabilizing residues (negative ΔΔG below threshold)
    stabilizing = df[df['ddg'] < stab_threshold].copy()
    stabilizing = stabilizing.sort_values('ddg', ascending=True)
    
    # Top N most significant (by absolute value)
    df_copy = df.copy()
    df_copy['abs_ddg'] = df_copy['ddg'].abs()
    top_significant = df_copy.nlargest(top_n, 'abs_ddg').drop(columns=['abs_ddg'])
    
    return destabilizing, stabilizing, top_significant


def main():
    parser = argparse.ArgumentParser(
        description="Extract critical residues from MutateX results"
    )
    parser.add_argument(
        "results_dir",
        help="Directory containing MutateX results"
    )
    parser.add_argument(
        "pdb_file",
        help="Path to the original PDB file"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output directory (default: <results_dir>/critical_residues)"
    )
    parser.add_argument(
        "--destabilizing",
        type=float,
        default=2.0,
        help="DDG threshold for destabilizing mutations (default: 2.0 kcal/mol)"
    )
    parser.add_argument(
        "--stabilizing",
        type=float,
        default=-2.0,
        help="DDG threshold for stabilizing mutations (default: -2.0 kcal/mol)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top residues to extract (default: 20)"
    )
    
    args = parser.parse_args()
    
    results_dir = args.results_dir
    pdb_file = args.pdb_file
    output_dir = args.output or os.path.join(results_dir, "critical_residues")
    
    if not os.path.isdir(results_dir):
        print(f"[ERROR] Results directory not found: {results_dir}")
        sys.exit(1)
    
    if not os.path.isfile(pdb_file):
        print(f"[WARNING] PDB file not found: {pdb_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse data
    print(f"[INFO] Parsing MutateX results from {results_dir}")
    interface_df, folding_df = find_and_parse_data(results_dir)
    
    if interface_df.empty and folding_df.empty:
        print(f"[WARNING] No data found in {results_dir}")
        sys.exit(0)
    
    pdb_name = Path(pdb_file).stem
    
    for data_type, df in [('interface', interface_df), ('folding', folding_df)]:
        if df.empty:
            continue
        print(f"[INFO] Processing {len(df)} {data_type} residues")
        destab, stab, top = extract_critical_residues(
            df, args.destabilizing, args.stabilizing, args.top_n, data_type
        )
        
        if not destab.empty:
            destab_file = os.path.join(output_dir, f"{pdb_name}_{data_type}_destabilizing.csv")
            destab.to_csv(destab_file, index=False)
            print(f"  [SAVED] {len(destab)} destabilizing residues -> {destab_file}")
        
        if not stab.empty:
            stab_file = os.path.join(output_dir, f"{pdb_name}_{data_type}_stabilizing.csv")
            stab.to_csv(stab_file, index=False)
            print(f"  [SAVED] {len(stab)} stabilizing residues -> {stab_file}")
        
        top_file = os.path.join(output_dir, f"{pdb_name}_{data_type}_top{args.top_n}.csv")
        top.to_csv(top_file, index=False)
        print(f"  [SAVED] Top {len(top)} significant residues -> {top_file}")
    
    print(f"[SUCCESS] Critical residues extracted to {output_dir}")


if __name__ == "__main__":
    main()
