#!/usr/bin/env python3
"""
Matrix Filtering Utility
Filters rows (genes) and columns (motifs/functions) from PlantCARE matrices
"""

import pandas as pd
import argparse
import sys

def filter_matrix(input_file, output_file, genes=None, columns=None):
    """
    Filter matrix by rows (genes) and/or columns (motifs/functions)
    """
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    # Filter genes (rows)
    if genes:
        gene_list = [g.strip() for g in genes.split(',')]
        available_genes = [g for g in gene_list if g in df.index]
        if not available_genes:
            print(f"Error: None of the specified genes found in matrix")
            sys.exit(1)
        df = df.loc[available_genes]
        print(f"Filtered to {len(available_genes)} gene(s)")
    
    # Filter columns (motifs/functions)
    if columns:
        col_list = [c.strip() for c in columns.split(',')]
        available_cols = [c for c in col_list if c in df.columns]
        if not available_cols:
            print(f"Error: None of the specified columns found in matrix")
            sys.exit(1)
        df = df[available_cols]
        print(f"Filtered to {len(available_cols)} column(s)")
    
    df.to_csv(output_file, sep='\t')
    print(f"Filtered matrix saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Filter PlantCARE matrix by genes and/or columns.')
    parser.add_argument('-i', '--input', required=True, help='Input matrix file (TSV format)')
    parser.add_argument('-o', '--output', required=True, help='Output filtered matrix file')
    parser.add_argument('--genes', help='Comma-separated list of genes to keep')
    parser.add_argument('--columns', help='Comma-separated list of columns to keep')
    
    args = parser.parse_args()
    
    if not args.genes and not args.columns:
        print("Warning: No filters specified. Output will be identical to input.")
    
    filter_matrix(args.input, args.output, args.genes, args.columns)

if __name__ == '__main__':
    main()
