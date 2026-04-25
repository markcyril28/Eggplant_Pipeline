#!/usr/bin/env python3
"""
Filter CRISPR gRNA scoring results by on-target score threshold.

Reads raw CSV results (seq_id, sgRNA_id, Score, Sequence, strand, pos, %GC),
filters by score >= threshold, assigns a Tier column (High >= 0.7, Moderate >= 0.5),
and writes per-gene and combined filtered CSVs.

Usage:
    python3 filter_scores.py --input-dir <raw_dir> --output-dir <out_dir> --threshold 0.5
    python3 filter_scores.py --input-dir <raw_dir> --output-dir <out_dir> --threshold 0.7
"""

import argparse
import csv
import os
import sys


def assign_tier(score: float) -> str:
    """Assign quality tier based on on-target score."""
    if score >= 0.7:
        return "High"
    elif score >= 0.5:
        return "Moderate"
    else:
        return "Low"


def filter_csv(input_path: str, output_path: str, threshold: float) -> int:
    """Filter a single CSV by score threshold. Returns number of rows kept."""
    kept = []
    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            score = float(row["Score"])
            if score >= threshold:
                row["Tier"] = assign_tier(score)
                kept.append(row)
    if not kept:
        return 0
    fieldnames = list(kept[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(kept)
    return len(kept)


def main():
    parser = argparse.ArgumentParser(description="Filter CRISPR scores by threshold")
    parser.add_argument("--input-dir", required=True, help="Directory with raw per-gene CSVs")
    parser.add_argument("--output-dir", required=True, help="Output directory for filtered CSVs")
    parser.add_argument("--threshold", type=float, required=True, help="Minimum score threshold")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    combined_rows = []
    gene_files = sorted(
        f for f in os.listdir(args.input_dir)
        if f.endswith(".csv") and f != "CRISPR_results_combined.csv"
    )

    if not gene_files:
        print(f"WARNING: No per-gene CSV files found in {args.input_dir}", file=sys.stderr)
        sys.exit(0)

    total_kept = 0
    for fname in gene_files:
        gene_name = fname.replace(".csv", "")
        in_path = os.path.join(args.input_dir, fname)
        out_name = f"{gene_name}_filtered_{args.threshold}.csv"
        out_path = os.path.join(args.output_dir, out_name)

        kept = filter_csv(in_path, out_path, args.threshold)
        total_kept += kept

        # Collect rows for combined output
        if kept > 0:
            with open(out_path, newline="") as fh:
                reader = csv.DictReader(fh)
                combined_rows.extend(list(reader))

    # Write combined filtered CSV
    if combined_rows:
        combined_path = os.path.join(
            args.output_dir, f"CRISPR_results_combined_filtered_{args.threshold}.csv"
        )
        fieldnames = list(combined_rows[0].keys())
        with open(combined_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(combined_rows)

    print(f"Filtered {total_kept} guides (score >= {args.threshold}) -> {args.output_dir}")


if __name__ == "__main__":
    main()
