#!/usr/bin/env python3
"""
PlantCARE to Matrix Converter
Converts PlantCARE tab-delimited output into matrices of CAREs.
Version 1: Matrix by Motif Name
Version 2: Matrix by Function Description
"""

import pandas as pd
import argparse
import os

def create_count_matrix(df, group_by_col):
    """
    Create a matrix with counts of each motif/function per sequence.
    """
    count_matrix = df.groupby(['Sequence_ID', group_by_col]).size().unstack(fill_value=0)
    return count_matrix

def process_heatmap_files(input_dir, output_dir, gene_order=None):
    """
    Reads all heatmap files and creates two versions:
    Version 1: Matrix by Motif Name
    Version 2: Matrix by Function Description
    """
    all_dfs = []
    for filename in os.listdir(input_dir):
        if filename.endswith('_heatmap.tab'):
            file_path = os.path.join(input_dir, filename)
            df = pd.read_csv(file_path, sep='\t', header=0)
            all_dfs.append(df)

    if not all_dfs:
        print(f"No heatmap files found in '{input_dir}'.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # Clean data
    combined_df = combined_df[combined_df['Motif_Name'].notna()]
    combined_df = combined_df[combined_df['Motif_Name'] != '']

    # Version 1: Matrix by Motif Name
    count_matrix_v1 = create_count_matrix(combined_df, 'Motif_Name')
    if gene_order:
        missing_genes = set(gene_order) - set(count_matrix_v1.index)
        for gene in missing_genes:
            count_matrix_v1.loc[gene] = 0
        count_matrix_v1 = count_matrix_v1.reindex(gene_order)
    
    output_v1 = os.path.join(output_dir, 'plantcare_heatmap_matrix_v1_motif.tsv')
    count_matrix_v1.to_csv(output_v1, sep='\t')
    print(f"Version 1 (Motif Name) matrix saved to '{output_v1}'")

    # Version 2: Matrix by Function Description
    combined_df_v2 = combined_df[combined_df['Function'].notna()]
    combined_df_v2 = combined_df_v2[combined_df_v2['Function'] != '']
    
    if len(combined_df_v2) > 0:
        count_matrix_v2 = create_count_matrix(combined_df_v2, 'Function')
        if gene_order:
            missing_genes = set(gene_order) - set(count_matrix_v2.index)
            for gene in missing_genes:
                count_matrix_v2.loc[gene] = 0
            count_matrix_v2 = count_matrix_v2.reindex(gene_order)
        
        output_v2 = os.path.join(output_dir, 'plantcare_heatmap_matrix_v2_function.tsv')
        count_matrix_v2.to_csv(output_v2, sep='\t')
        print(f"Version 2 (Function) matrix saved to '{output_v2}'")
    else:
        print("Warning: No functions found in data. Version 2 matrix not created.")

def process_tbtools_files(input_dir, output_dir):
    """
    Reads all TBTools files and creates a matrix.
    """
    all_dfs = []
    for filename in os.listdir(input_dir):
        if filename.endswith('_tbtools.tab'):
            file_path = os.path.join(input_dir, filename)
            columns=['Sequence_ID', 'Motif_Name', 'Start', 'End', 'Strand']
            df = pd.read_csv(file_path, sep='\t', header=None)
            df.columns = columns
            all_dfs.append(df)

    if not all_dfs:
        print(f"No TBTools files found in '{input_dir}'.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df = combined_df[combined_df['Motif_Name'].notna()]
    combined_df = combined_df[combined_df['Motif_Name'] != '']

    count_matrix = create_count_matrix(combined_df, 'Motif_Name')
    
    output_path = os.path.join(output_dir, 'plantcare_tbtools_matrix.tsv')
    count_matrix.to_csv(output_path, sep='\t')
    print(f"TBTools matrix saved to '{output_path}'")

def main():
    """Main function to run the matrix generation."""
    parser = argparse.ArgumentParser(description='Generate matrices from post-processed PlantCARE files.')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory with post-processed tab files.')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for matrices.')
    parser.add_argument('--genes', default=None, help='Comma-separated gene order for matrix rows (optional; defaults to sorted).')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Gene order is derived from data (sorted alphabetically); override via --genes
    gene_order = None
    if hasattr(args, 'genes') and args.genes:
        gene_order = [g.strip() for g in args.genes.split(',')]

    # Process heatmap files (creates both Version 1 and Version 2)
    process_heatmap_files(args.input_dir, args.output_dir, gene_order=gene_order)

    # Process TBTools files
    process_tbtools_files(args.input_dir, args.output_dir)

if __name__ == '__main__':
    main()