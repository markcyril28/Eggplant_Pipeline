#!/usr/bin/env python3
"""Utility to explore PlantCARE data structure."""

import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict

def analyze_raw_data(raw_dir, genes):
    """Analyze raw PlantCARE data."""
    print("=== PlantCARE Data Analysis ===\n")
    
    all_functions = set()
    all_groups = set()
    all_motifs = set()
    group_func_map = defaultdict(set)
    
    for gene in genes:
        tab_file = Path(raw_dir) / f"{gene}.tab"
        if not tab_file.exists():
            print(f"Warning: {gene}.tab not found")
            continue
            
        df = pd.read_csv(tab_file, sep='\t', header=None,
                        names=['ID', 'Motif', 'Seq', 'Pos', 'Len', 'Strand', 'Species', 'Function', 'MajorGroup'])
        
        for _, row in df.iterrows():
            func = str(row['Function']).strip()
            group = str(row['MajorGroup']).strip()
            motif = str(row['Motif']).strip()
            
            if func and func != 'nan':
                all_functions.add(func)
            if group and group != 'nan':
                all_groups.add(group)
            if motif:
                all_motifs.add(motif)
            if func and group and func != 'nan' and group != 'nan':
                group_func_map[group].add(func)
    
    print(f"Genes analyzed: {len(genes)}")
    print(f"Unique Motifs: {len(all_motifs)}")
    print(f"Unique Functions: {len(all_functions)}")
    print(f"Unique Major Groups: {len(all_groups)}\n")
    
    print("Major Groups:")
    for i, group in enumerate(sorted(all_groups), 1):
        func_count = len(group_func_map[group])
        print(f"  {i}. {group} ({func_count} functions)")
    
    print("\nFunctions by Major Group:")
    for group in sorted(all_groups):
        print(f"\n{group}:")
        for func in sorted(group_func_map[group]):
            print(f"  - {func}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_dir', default='2_PlantCARE_Results/Raw_Results')
    parser.add_argument('--genes', nargs='+', 
                       default=['SmelGRF05_970', 'SmelGRF08_140', 'SmelGIF11_070', 'SmelGIF11_650', 'SmelGIF11_790'])
    args = parser.parse_args()
    
    analyze_raw_data(args.raw_dir, args.genes)

if __name__ == '__main__':
    main()
