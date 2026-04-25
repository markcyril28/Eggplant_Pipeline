#!/usr/bin/env python3
"""
PlantCARE to Matrix Converter
Converts PlantCARE tab-delimited output into a matrix of CAREs.

Input: Tab-delimited file with columns:
    [0] Sequence_ID, [1] Motif_Name, [2] Motif_Sequence, [3] Position, 
    [4] Length, [5] Strand, [6] Organism, [7] Function

Output: Multiple matrix formats for analysis
"""

import pandas as pd
import numpy as np
import sys
import os
from collections import defaultdict
import argparse


def parse_plantcare_file(input_file):
    """
    Parse PlantCARE tab-delimited file.
    
    Args:
        input_file: Path to input tab file
        
    Returns:
        DataFrame with parsed data
    """
    # Define column names
    columns = ['Sequence_ID', 'Motif_Name', 'Motif_Sequence', 'Position', 
               'Length', 'Strand', 'Organism', 'Function']
    
    # Read the file
    df = pd.read_csv(input_file, sep='\t', names=columns, header=None)
    
    # Clean the data
    # Remove rows with placeholder data
    df = df[df['Motif_Sequence'] != 'motif_sequence']
    df = df[df['Organism'] != 'organism']
    
    # Remove empty motif names
    df = df[df['Motif_Name'].notna()]
    df = df[df['Motif_Name'] != '']
    
    # Convert position to integer
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
    df['Length'] = pd.to_numeric(df['Length'], errors='coerce')
    
    # Fill NaN in Function column with empty string
    df['Function'] = df['Function'].fillna('')
    
    return df


def create_count_matrix(df):
    """
    Create a matrix with counts of each motif per sequence.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with sequences as rows and motifs as columns
    """
    # Count occurrences of each motif per sequence
    count_matrix = df.groupby(['Sequence_ID', 'Motif_Name']).size().unstack(fill_value=0)
    
    return count_matrix


def create_position_matrix(df):
    """
    Create a matrix with positions of each motif.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with motif positions (comma-separated if multiple)
    """
    # Group positions by sequence and motif
    position_groups = df.groupby(['Sequence_ID', 'Motif_Name'])['Position'].apply(
        lambda x: ','.join(map(str, sorted(x)))
    ).unstack(fill_value='')
    
    return position_groups


def create_strand_matrix(df):
    """
    Create a matrix showing strand distribution (+/-) for each motif.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with strand information
    """
    def strand_summary(strands):
        plus = (strands == '+').sum()
        minus = (strands == '-').sum()
        return f"+:{plus}/-:{minus}"
    
    strand_matrix = df.groupby(['Sequence_ID', 'Motif_Name'])['Strand'].apply(
        strand_summary
    ).unstack(fill_value='')
    
    return strand_matrix


def create_detailed_matrix(df):
    """
    Create a detailed matrix with all information per motif.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with detailed motif information
    """
    detailed_list = []
    
    for seq_id in df['Sequence_ID'].unique():
        seq_data = df[df['Sequence_ID'] == seq_id]
        
        for motif_name in seq_data['Motif_Name'].unique():
            motif_data = seq_data[seq_data['Motif_Name'] == motif_name]
            
            # Get representative function (most common non-empty)
            functions = motif_data['Function'].value_counts()
            function = functions.index[0] if len(functions) > 0 and functions.index[0] != '' else ''
            
            detailed_list.append({
                'Sequence_ID': seq_id,
                'Motif_Name': motif_name,
                'Count': len(motif_data),
                'Positions': ','.join(map(str, sorted(motif_data['Position']))),
                'Strands': ','.join(motif_data['Strand']),
                'Sequences': ','.join(motif_data['Motif_Sequence'].unique()),
                'Function': function
            })
    
    return pd.DataFrame(detailed_list)


def create_functional_category_matrix(df):
    """
    Create a matrix grouped by functional categories.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with functional categories
    """
    # Define functional categories based on keywords
    def categorize_function(function_text):
        function_lower = str(function_text).lower()
        
        if any(keyword in function_lower for keyword in ['light', 'box 4', 'ae-box', 'g-box']):
            return 'Light Responsiveness'
        elif any(keyword in function_lower for keyword in ['drought', 'mbs']):
            return 'Drought Response'
        elif any(keyword in function_lower for keyword in ['abscisic', 'abre', 'aba']):
            return 'ABA Response'
        elif any(keyword in function_lower for keyword in ['anaerobic', 'are']):
            return 'Anaerobic Response'
        elif any(keyword in function_lower for keyword in ['stress', 'stre']):
            return 'Stress Response'
        elif any(keyword in function_lower for keyword in ['promoter', 'tata', 'caat']):
            return 'Core Promoter'
        elif any(keyword in function_lower for keyword in ['myb', 'myc']):
            return 'Transcription Factor Binding'
        else:
            return 'Other/Unknown'
    
    # Add category column
    df_cat = df.copy()
    df_cat['Category'] = df_cat['Function'].apply(categorize_function)
    
    # Count by category
    category_matrix = df_cat.groupby(['Sequence_ID', 'Category']).size().unstack(fill_value=0)
    
    return category_matrix


def create_summary_statistics(df):
    """
    Create summary statistics for the data.
    
    Args:
        df: DataFrame with parsed PlantCARE data
        
    Returns:
        DataFrame with summary statistics
    """
    summary_data = []
    
    for seq_id in df['Sequence_ID'].unique():
        seq_data = df[df['Sequence_ID'] == seq_id]
        
        summary_data.append({
            'Sequence_ID': seq_id,
            'Total_Motifs': len(seq_data),
            'Unique_Motif_Types': seq_data['Motif_Name'].nunique(),
            'Plus_Strand': (seq_data['Strand'] == '+').sum(),
            'Minus_Strand': (seq_data['Strand'] == '-').sum(),
            'Min_Position': seq_data['Position'].min(),
            'Max_Position': seq_data['Position'].max(),
            'Avg_Motif_Length': seq_data['Length'].mean()
        })
    
    return pd.DataFrame(summary_data)


def main():
    """Main function to run the pipeline."""
    parser = argparse.ArgumentParser(
        description='Convert PlantCARE tab file to matrix format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python plantCARE_to_matrix.py input.tab
  
  # Specify output directory
  python plantCARE_to_matrix.py input.tab -o output_dir/
  
  # Specify output prefix
  python plantCARE_to_matrix.py input.tab -p my_analysis
        """
    )
    
    parser.add_argument('input_file', help='Input PlantCARE tab file')
    parser.add_argument('-o', '--output_dir', default='.', 
                        help='Output directory (default: current directory)')
    parser.add_argument('-p', '--prefix', default='plantCARE_matrix',
                        help='Output file prefix (default: plantCARE_matrix)')
    parser.add_argument('--no-header', action='store_true',
                        help='Do not include headers in output files')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found!")
        sys.exit(1)
    
    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse input file
    print(f"Parsing input file: {args.input_file}")
    df = parse_plantcare_file(args.input_file)
    print(f"Found {len(df)} valid motif entries")
    print(f"Unique motif types: {df['Motif_Name'].nunique()}")
    
    # Generate output file paths
    output_base = os.path.join(args.output_dir, args.prefix)
    
    # Create and save count matrix
    print("\nGenerating count matrix...")
    count_matrix = create_count_matrix(df)
    count_output = f"{output_base}_count_matrix.tsv"
    count_matrix.to_csv(count_output, sep='\t')
    print(f"  Saved: {count_output}")
    
    # Create and save position matrix
    print("Generating position matrix...")
    position_matrix = create_position_matrix(df)
    position_output = f"{output_base}_position_matrix.tsv"
    position_matrix.to_csv(position_output, sep='\t')
    print(f"  Saved: {position_output}")
    
    # Create and save strand matrix
    print("Generating strand distribution matrix...")
    strand_matrix = create_strand_matrix(df)
    strand_output = f"{output_base}_strand_matrix.tsv"
    strand_matrix.to_csv(strand_output, sep='\t')
    print(f"  Saved: {strand_output}")
    
    # Create and save detailed matrix
    print("Generating detailed matrix...")
    detailed_matrix = create_detailed_matrix(df)
    detailed_output = f"{output_base}_detailed.tsv"
    detailed_matrix.to_csv(detailed_output, sep='\t', index=False)
    print(f"  Saved: {detailed_output}")
    
    # Create and save functional category matrix
    print("Generating functional category matrix...")
    category_matrix = create_functional_category_matrix(df)
    category_output = f"{output_base}_functional_categories.tsv"
    category_matrix.to_csv(category_output, sep='\t')
    print(f"  Saved: {category_output}")
    
    # Create and save summary statistics
    print("Generating summary statistics...")
    summary_stats = create_summary_statistics(df)
    summary_output = f"{output_base}_summary.tsv"
    summary_stats.to_csv(summary_output, sep='\t', index=False)
    print(f"  Saved: {summary_output}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total motif entries processed: {len(df)}")
    print(f"Unique sequences: {df['Sequence_ID'].nunique()}")
    print(f"Unique motif types: {df['Motif_Name'].nunique()}")
    print("\nTop 10 most frequent motifs:")
    print(df['Motif_Name'].value_counts().head(10))
    print("\n" + "="*60)
    print("Pipeline completed successfully!")
    print("="*60)


if __name__ == '__main__':
    main()
