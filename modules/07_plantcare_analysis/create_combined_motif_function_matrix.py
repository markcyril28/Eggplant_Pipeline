#!/usr/bin/env python3
"""Create combined motif and function matrix for PlantCARE data."""

import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict

def load_function_mapping(function_matrix_path):
    """Load function names from matrix (short names)."""
    df = pd.read_csv(function_matrix_path, sep='\t', index_col=0)
    return list(df.columns)

def load_raw_data(raw_dir, genes):
    """Load raw PlantCARE tab files."""
    data = {}
    for gene in genes:
        tab_file = Path(raw_dir) / f"{gene}.tab"
        if tab_file.exists():
            df = pd.read_csv(tab_file, sep='\t', header=None, 
                           names=['ID', 'Motif', 'Seq', 'Pos', 'Len', 'Strand', 'Species', 'Function'])
            data[gene] = df
    return data

def create_function_to_motif_mapping(raw_data):
    """Map functions to their motifs from raw data."""
    func_motif_map = defaultdict(set)
    for gene, df in raw_data.items():
        for _, row in df.iterrows():
            func = str(row['Function']).strip()
            motif = str(row['Motif']).strip()
            if func and func != 'nan' and motif:
                func_motif_map[func].add(motif)
    return {k: sorted(list(v)) for k, v in func_motif_map.items()}

def count_motifs_by_function_and_gene(raw_data, function_list):
    """Count motifs grouped by function for each gene."""
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    
    for gene, df in raw_data.items():
        for _, row in df.iterrows():
            func = str(row['Function']).strip()
            motif = str(row['Motif']).strip()
            
            if func in function_list and motif:
                result[gene][func][motif] += 1
    
    return result

def create_combined_matrix(raw_data, function_list, genes):
    """Create matrix with motifs grouped by function."""
    func_motif_map = create_function_to_motif_mapping(raw_data)
    counts = count_motifs_by_function_and_gene(raw_data, function_list)
    
    # Build column list: Function label (no data column) | Motif1 | Motif2 | ...
    columns = []
    col_metadata = []
    
    for func in function_list:
        if func in func_motif_map:
            motifs = func_motif_map[func]
            # Add function as label only (no data column)
            for motif in motifs:
                columns.append(motif)
                col_metadata.append({'type': 'motif', 'name': motif, 'function': func})
    
    # Build data matrix
    matrix_data = []
    for gene in genes:
        row = []
        for col_meta in col_metadata:
            # Only motif counts (no function total column)
            func = col_meta['function']
            motif = col_meta['name']
            count = counts[gene][func][motif] if gene in counts and func in counts[gene] else 0
            row.append(count)
        matrix_data.append(row)
    
    df = pd.DataFrame(matrix_data, index=genes, columns=columns)
    return df, col_metadata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_dir', required=True)
    parser.add_argument('--function_matrix', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--genes', nargs='+', required=True)
    parser.add_argument('--functions', nargs='*', default=[])
    parser.add_argument('--option', default='A')
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    raw_data = load_raw_data(args.raw_dir, args.genes)
    
    # Get all unique functions from raw data
    all_functions = set()
    for gene, df in raw_data.items():
        for func in df['Function'].dropna().unique():
            func_str = str(func).strip()
            if func_str and func_str != 'nan':
                all_functions.add(func_str)
    
    # Filter functions based on option
    if args.option == 'C' and args.functions:
        function_list = [f for f in args.functions if f in all_functions]
        print(f"Option C: Selected {len(function_list)} out of {len(args.functions)} requested functions")
        if len(function_list) < len(args.functions):
            missing = [f for f in args.functions if f not in all_functions]
            print(f"Warning: {len(missing)} functions not found in data")
    else:
        function_list = sorted(list(all_functions))
        print(f"Option {args.option}: Using all {len(function_list)} functions from data")
    
    # Create combined matrix
    matrix, metadata = create_combined_matrix(raw_data, function_list, args.genes)
    
    # Get function to motif mapping for metadata
    func_motif_map = create_function_to_motif_mapping(raw_data)
    
    # Save outputs
    matrix.to_csv(output_dir / 'combined_matrix.tsv', sep='\t')
    
    # Save metadata (includes function grouping info)
    with open(output_dir / 'column_metadata.txt', 'w') as f:
        for i, meta in enumerate(metadata):
            f.write(f"{i}\t{meta['name']}\tmotif\t{meta['function']}\n")
    
    # Save function list for R script
    with open(output_dir / 'function_groups.txt', 'w') as f:
        for func in function_list:
            if func in func_motif_map:
                motif_count = len(func_motif_map[func])
                f.write(f"{func}\t{motif_count}\n")
    
    print(f"Combined matrix created: {matrix.shape}")
    print(f"Genes: {len(args.genes)}, Functions: {len(function_list)}, Columns: {len(metadata)}")

if __name__ == '__main__':
    main()
