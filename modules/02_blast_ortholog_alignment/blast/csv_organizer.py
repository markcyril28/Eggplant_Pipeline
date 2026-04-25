#!/usr/bin/env python3
"""
Organize BLAST-like CSV by sorting first on Subject ID,
then on E-value (ascending, smaller values first).
"""

import pandas as pd
import argparse
import os

def organize_csv(input_csv, output_csv=None, create_raw=True):
    # Read CSV
    df = pd.read_csv(input_csv)

    # Sort by Subject ID (asc), then by E-value (asc)
    df_sorted = df.sort_values(by=["Subject ID", "E-value"], ascending=[True, True])

    # Generate output filenames if not provided
    if output_csv is None:
        # Get the directory and base filename
        input_dir = os.path.dirname(input_csv)
        input_basename = os.path.basename(input_csv)
        # Remove .csv extension
        name_without_ext = os.path.splitext(input_basename)[0]
        
        # Create organized filename
        organized_filename = f"{name_without_ext}_organized.csv"
        output_csv = os.path.join(input_dir, organized_filename)
        
        # Always create raw version when using automatic naming
        create_raw = True

    # Create raw version if requested (when using automatic naming)
    if create_raw:
        input_dir = os.path.dirname(input_csv)
        input_basename = os.path.basename(input_csv)
        name_without_ext = os.path.splitext(input_basename)[0]
        raw_filename = f"{name_without_ext}_raw.csv"
        raw_csv = os.path.join(input_dir, raw_filename)
        
        # Save raw version (copy of original)
        df.to_csv(raw_csv, index=False)
        print(f"Raw CSV saved to: {raw_csv}")

    # Save sorted CSV
    df_sorted.to_csv(output_csv, index=False)
    print(f"Organized CSV saved to: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organize a CSV by Subject ID and E-value.")
    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("output_csv", nargs='?', default=None, help="Path to the output organized CSV file (optional - will auto-generate if not provided)")

    args = parser.parse_args()
    organize_csv(args.input_csv, args.output_csv)
