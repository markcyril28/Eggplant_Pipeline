#!/usr/bin/env python3
"""Create matrix with functions organized by major groups."""

import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict

def load_raw_data(raw_dir, genes):
    """Load raw PlantCARE tab files."""
    data = {}
    for gene in genes:
        tab_file = Path(raw_dir) / f"{gene}.tab"
        if tab_file.exists():
            df = pd.read_csv(tab_file, sep='\t', header=None, 
                           names=['ID', 'Motif', 'Seq', 'Pos', 'Len', 'Strand', 'Species', 'Function', 'MajorGroup'])
            data[gene] = df
    return data

def get_function_groups(raw_data):
    """Extract unique functions organized by major groups."""
    group_func_map = defaultdict(set)
    func_group_map = {}
    
    for gene, df in raw_data.items():
        for _, row in df.iterrows():
            func = str(row['Function']).strip()
            group = str(row['MajorGroup']).strip()
            
            if func and func != 'nan' and group and group != 'nan':
                group_func_map[group].add(func)
                func_group_map[func] = group
    
    # Sort groups and functions
    sorted_groups = sorted(group_func_map.keys())
    ordered_functions = []
    for group in sorted_groups:
        funcs = sorted(list(group_func_map[group]))
        ordered_functions.extend(funcs)
    
    return ordered_functions, func_group_map, sorted_groups

def count_functions_by_gene(raw_data, function_list):
    """Count function occurrences for each gene."""
    result = defaultdict(lambda: defaultdict(int))
    
    for gene, df in raw_data.items():
        for _, row in df.iterrows():
            func = str(row['Function']).strip()
            if func in function_list:
                result[gene][func] += 1
    
    return result

def create_matrix(raw_data, genes, selected_groups=None):
    """Create matrix with functions organized by major groups."""
    all_functions, func_group_map, all_groups = get_function_groups(raw_data)
    
    # Filter by selected groups
    if selected_groups:
        # Preserve order from selected_groups and organize functions within each group
        filtered_functions = []
        for group in selected_groups:
            group_funcs = sorted([f for f in all_functions if func_group_map.get(f) == group])
            filtered_functions.extend(group_funcs)

        if not filtered_functions:
            print("Warning: no functions matched selected major groups.")
            print(f"Selected groups ({len(selected_groups)}): {', '.join(selected_groups)}")
            print(f"Available groups ({len(all_groups)}): {', '.join(all_groups)}")
    else:
        filtered_functions = all_functions
    
    # Count occurrences
    counts = count_functions_by_gene(raw_data, filtered_functions)
    
    # Build matrix
    matrix_data = []
    for gene in genes:
        row = [counts[gene][func] for func in filtered_functions]
        matrix_data.append(row)
    
    df = pd.DataFrame(matrix_data, index=genes, columns=filtered_functions)
    
    # Build metadata for R script preserving order
    metadata = []
    for func in filtered_functions:
        group = func_group_map.get(func, 'Unknown')
        metadata.append({'function': func, 'group': group})
    
    return df, metadata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--genes', nargs='+', required=True)
    parser.add_argument('--functions', nargs='*', default=[])
    parser.add_argument('--groups', nargs='*', default=[])
    parser.add_argument('--option', default='A')
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    raw_data = load_raw_data(args.raw_dir, args.genes)
    
    # Determine functions to include
    selected_groups = args.groups if args.option == 'C' and args.groups else None
    
    # Create matrix
    matrix, metadata = create_matrix(raw_data, args.genes, selected_groups)
    
    # Save matrix
    matrix.to_csv(output_dir / 'function_groups_matrix.tsv', sep='\t')
    
    # Save metadata
    with open(output_dir / 'function_metadata.txt', 'w') as f:
        for i, meta in enumerate(metadata):
            f.write(f"{i}\t{meta['function']}\t{meta['group']}\n")
    
    # Save group info - preserve order
    groups_dict = {}
    for meta in metadata:
        group = meta['group']
        if group not in groups_dict:
            groups_dict[group] = 0
        groups_dict[group] += 1
    
    with open(output_dir / 'group_info.txt', 'w') as f:
        for group, count in groups_dict.items():
            f.write(f"{group}\t{count}\n")
    
    print(f"Matrix: {matrix.shape[0]} genes x {matrix.shape[1]} functions")
    print(f"Major groups: {len(groups_dict)}")

if __name__ == '__main__':
    main()
