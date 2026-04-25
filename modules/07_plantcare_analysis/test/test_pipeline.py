#!/usr/bin/env python3
"""
Quick test script for PlantCARE to Matrix pipeline
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plantCARE_to_matrix import (
    parse_plantcare_file,
    create_count_matrix,
    create_position_matrix,
    create_detailed_matrix,
    create_functional_category_matrix,
    create_summary_statistics
)

def test_pipeline():
    """Test the pipeline with the sample file."""
    
    input_file = "plantCARE_output_PlantCARE_5913.tab"
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found!")
        return False
    
    print("=" * 60)
    print("Testing PlantCARE to Matrix Pipeline")
    print("=" * 60)
    
    # Parse file
    print("\n1. Parsing input file...")
    df = parse_plantcare_file(input_file)
    print(f"   ✓ Parsed {len(df)} valid entries")
    print(f"   ✓ Found {df['Sequence_ID'].nunique()} unique sequence(s)")
    print(f"   ✓ Found {df['Motif_Name'].nunique()} unique motif types")
    
    # Create count matrix
    print("\n2. Creating count matrix...")
    count_matrix = create_count_matrix(df)
    print(f"   ✓ Matrix shape: {count_matrix.shape}")
    print(f"   ✓ Preview:")
    print(count_matrix.head().to_string())
    
    # Create position matrix
    print("\n3. Creating position matrix...")
    position_matrix = create_position_matrix(df)
    print(f"   ✓ Matrix shape: {position_matrix.shape}")
    
    # Create detailed matrix
    print("\n4. Creating detailed matrix...")
    detailed_matrix = create_detailed_matrix(df)
    print(f"   ✓ Detailed entries: {len(detailed_matrix)}")
    print(f"   ✓ Sample entries:")
    print(detailed_matrix.head().to_string())
    
    # Create functional category matrix
    print("\n5. Creating functional category matrix...")
    category_matrix = create_functional_category_matrix(df)
    print(f"   ✓ Matrix shape: {category_matrix.shape}")
    print(f"   ✓ Categories found:")
    print(category_matrix.T.to_string())
    
    # Create summary statistics
    print("\n6. Creating summary statistics...")
    summary = create_summary_statistics(df)
    print(f"   ✓ Summary:")
    print(summary.to_string(index=False))
    
    # Top motifs
    print("\n7. Top 10 most frequent motifs:")
    top_motifs = df['Motif_Name'].value_counts().head(10)
    for motif, count in top_motifs.items():
        print(f"   • {motif}: {count}")
    
    print("\n" + "=" * 60)
    print("✓ All tests passed successfully!")
    print("=" * 60)
    
    return True

if __name__ == '__main__':
    success = test_pipeline()
    sys.exit(0 if success else 1)
