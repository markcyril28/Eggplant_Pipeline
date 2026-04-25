#!/usr/bin/env python3
"""
CRISPR Off-Target Analysis — Summary Report Generator

Reads filtered scoring CSVs, BLAST off-target results, and Cas-OFFinder output,
then produces:
  - guide_summary.csv        (per-guide consolidated table)
  - CRISPR_summary_report.txt (human-readable text report)

Usage:
    python3 generate_report.py --crispr-dir <dir> --genome <name> \
        [--score-thresholds 0.5 0.7] [--blast-dir <dir>] [--casoff-dir <dir>]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path


# -- helpers ----------------------------------------------------------------

def read_csv_rows(path: str) -> list[dict]:
    """Read a CSV file and return list of row dicts."""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def count_blast_offtargets(blast_results_dir: str, guide_name: str) -> dict:
    """Count BLAST off-target hits for a given guide name (sgRNA_id) across all result dirs.
    Reads tabular blastn output (format 6) and counts rows where qseqid == guide_name.
    Filtered counts come from per-guide files written by off_target_blast.sh (awk split).
    """
    counts = {"blast_total_hits": 0, "blast_filtered_hits": 0}
    if not blast_results_dir or not os.path.isdir(blast_results_dir):
        return counts

    results_dir = os.path.join(blast_results_dir, "results")
    if not os.path.isdir(results_dir):
        return counts

    for grna_dir_name in os.listdir(results_dir):
        grna_dir = os.path.join(results_dir, grna_dir_name)
        if not os.path.isdir(grna_dir):
            continue

        # Count total hits: scan aggregate tabular file, match qseqid (column 0) to guide_name
        all_file = os.path.join(grna_dir, f"{grna_dir_name}_all.txt")
        if os.path.isfile(all_file):
            with open(all_file) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if cols and cols[0] == guide_name:
                        counts["blast_total_hits"] += 1

        # Count filtered hits: per-guide file written by awk split in off_target_blast.sh
        per_guide_file = os.path.join(grna_dir, f"{guide_name}_filtered.txt")
        if os.path.isfile(per_guide_file):
            with open(per_guide_file) as fh:
                counts["blast_filtered_hits"] += sum(1 for line in fh if line.strip())

    return counts


def count_casoff_offtargets(casoff_dir: str, guide_seq: str) -> dict:
    """Count Cas-OFFinder hits for a guide sequence."""
    counts = {"casoff_total_hits": 0, "casoff_0mm": 0, "casoff_1mm": 0,
              "casoff_2mm": 0, "casoff_3mm": 0, "casoff_4mm": 0}
    if not casoff_dir or not os.path.isdir(casoff_dir):
        return counts

    # Cas-OFFinder output format: guide_seq chrom position DNA_seq strand mismatches
    output_dir = os.path.join(casoff_dir, "output")
    if os.path.isdir(output_dir):
        search_dir = output_dir
    elif os.path.isdir(os.path.join(casoff_dir, "cas_offinder", "output")):
        search_dir = os.path.join(casoff_dir, "cas_offinder", "output")
    else:
        return counts

    seq_clean = guide_seq.upper().replace(" ", "")
    for fname in os.listdir(search_dir):
        if not fname.endswith("_output.txt"):
            continue
        fpath = os.path.join(search_dir, fname)
        with open(fpath) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                if cols[0].upper() == seq_clean:
                    counts["casoff_total_hits"] += 1
                    mm = int(cols[5])
                    key = f"casoff_{mm}mm"
                    if key in counts:
                        counts[key] += 1

    return counts


# -- main -------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate CRISPR summary report")
    parser.add_argument("--crispr-dir", required=True, help="Root CRISPR output dir for one genome")
    parser.add_argument("--genome", required=True, help="Genome name (for report header)")
    parser.add_argument("--output-dir", required=True, help="Where to write report files")
    parser.add_argument("--score-thresholds", nargs="+", type=float, default=[0.5, 0.7])
    parser.add_argument("--blast-dir", default="", help="BLAST off-target results dir")
    parser.add_argument("--casoff-dir", default="", help="Cas-OFFinder results dir")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # -- Collect raw scoring data ---------------------------------------
    raw_dir = os.path.join(args.crispr_dir, "01_Raw_Scoring_Results_from_CRISPR-P_V2_0")
    if not os.path.isdir(raw_dir):
        print(f"ERROR: Raw scoring results not found: {raw_dir}", file=sys.stderr)
        sys.exit(1)

    # Read combined raw CSV
    combined_csv = os.path.join(raw_dir, "CRISPR_results_combined.csv")
    if os.path.isfile(combined_csv):
        all_guides = read_csv_rows(combined_csv)
    else:
        # Merge per-gene files
        all_guides = []
        for fname in sorted(os.listdir(raw_dir)):
            if fname.endswith(".csv"):
                all_guides.extend(read_csv_rows(os.path.join(raw_dir, fname)))

    if not all_guides:
        print("WARNING: No guide data found.", file=sys.stderr)
        sys.exit(0)

    # -- Parse Cas-OFFinder results into per-sequence lookup ------------
    casoff_lookup: dict[str, dict] = {}
    if args.casoff_dir and os.path.isdir(args.casoff_dir):
        for seq_row in all_guides:
            seq = seq_row.get("Sequence", "").upper().strip()
            if seq and seq not in casoff_lookup:
                casoff_lookup[seq] = count_casoff_offtargets(args.casoff_dir, seq)

    # -- Parse BLAST results into per-guide-name lookup -----------------
    blast_lookup: dict[str, dict] = {}
    if args.blast_dir and os.path.isdir(args.blast_dir):
        for seq_row in all_guides:
            gname = seq_row.get("sgRNA_id", "").strip()
            if gname and gname not in blast_lookup:
                blast_lookup[gname] = count_blast_offtargets(args.blast_dir, gname)

    # -- Build per-guide summary ----------------------------------------
    genes = defaultdict(list)
    guide_rows = []
    for row in all_guides:
        seq_id = row.get("seq_id", "unknown").split()[0]
        score = float(row.get("Score", 0))
        seq = row.get("Sequence", "").upper().strip()
        gc = row.get("%GC", "")
        strand = row.get("strand", "")
        sgRNA_id = row.get("sgRNA_id", "")

        tier = "High" if score >= 0.7 else ("Moderate" if score >= 0.5 else "Low")

        casoff = casoff_lookup.get(seq, {})
        blast = blast_lookup.get(sgRNA_id, {})

        guide_entry = {
            "gene": seq_id,
            "sgRNA_id": sgRNA_id,
            "Sequence": row.get("Sequence", ""),
            "Score": f"{score:.4f}",
            "Tier": tier,
            "strand": strand,
            "pos": row.get("pos", ""),
            "%GC": gc,
            "blast_total": str(blast.get("blast_total_hits", "")),
            "blast_filtered": str(blast.get("blast_filtered_hits", "")),
            "casoff_total": str(casoff.get("casoff_total_hits", "")),
            "casoff_0mm": str(casoff.get("casoff_0mm", "")),
            "casoff_1mm": str(casoff.get("casoff_1mm", "")),
            "casoff_2mm": str(casoff.get("casoff_2mm", "")),
            "casoff_3mm": str(casoff.get("casoff_3mm", "")),
            "casoff_4mm": str(casoff.get("casoff_4mm", "")),
        }
        guide_rows.append(guide_entry)
        genes[seq_id].append(guide_entry)

    # -- Write guide_summary.csv ----------------------------------------
    summary_csv = os.path.join(args.output_dir, "guide_summary.csv")
    fieldnames = list(guide_rows[0].keys())
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(guide_rows)
    print(f"  guide_summary.csv: {len(guide_rows)} guides -> {summary_csv}")

    # -- Write human-readable text report -------------------------------
    report_path = os.path.join(args.output_dir, "CRISPR_summary_report.txt")
    with open(report_path, "w", encoding="utf-8") as rpt:
        rpt.write("=" * 78 + "\n")
        rpt.write("  CRISPR Off-Target Analysis — Summary Report\n")
        rpt.write("=" * 78 + "\n\n")
        rpt.write(f"  Genome:         {args.genome}\n")
        rpt.write(f"  Total guides:   {len(all_guides)}\n")
        rpt.write(f"  Genes analyzed: {len(genes)}\n")

        # Tier breakdown
        tier_counts = defaultdict(int)
        for g in guide_rows:
            tier_counts[g["Tier"]] += 1
        rpt.write(f"\n  Tier breakdown:\n")
        for tier in ["High", "Moderate", "Low"]:
            count = tier_counts.get(tier, 0)
            pct = count / len(guide_rows) * 100 if guide_rows else 0
            rpt.write(f"    {tier:10s}  {count:3d}  ({pct:5.1f}%)\n")

        # Score statistics
        scores = [float(g["Score"]) for g in guide_rows]
        rpt.write(f"\n  Score statistics:\n")
        rpt.write(f"    Min:    {min(scores):.4f}\n")
        rpt.write(f"    Max:    {max(scores):.4f}\n")
        rpt.write(f"    Mean:   {sum(scores)/len(scores):.4f}\n")
        median = sorted(scores)[len(scores) // 2]
        rpt.write(f"    Median: {median:.4f}\n")

        # Filter summary for each threshold
        for t in args.score_thresholds:
            passed = sum(1 for s in scores if s >= t)
            rpt.write(f"\n  Score >= {t}: {passed}/{len(scores)} guides pass\n")

        # Per-gene breakdown
        rpt.write(f"\n{'-' * 78}\n")
        rpt.write(f"  Per-Gene Breakdown\n")
        rpt.write(f"{'-' * 78}\n\n")

        for gene in sorted(genes.keys()):
            gene_guides = genes[gene]
            gene_scores = [float(g["Score"]) for g in gene_guides]
            high = sum(1 for g in gene_guides if g["Tier"] == "High")
            mod = sum(1 for g in gene_guides if g["Tier"] == "Moderate")
            low = sum(1 for g in gene_guides if g["Tier"] == "Low")

            rpt.write(f"  {gene}\n")
            rpt.write(f"    Total guides: {len(gene_guides)}  "
                       f"(High: {high}, Moderate: {mod}, Low: {low})\n")
            rpt.write(f"    Score range:  {min(gene_scores):.4f} - {max(gene_scores):.4f}\n")

            # Top 3 guides
            top = sorted(gene_guides, key=lambda g: float(g["Score"]), reverse=True)[:3]
            rpt.write(f"    Top guides:\n")
            for i, g in enumerate(top, 1):
                casoff_str = ""
                if g["casoff_total"] and g["casoff_total"] != "0":
                    casoff_str = f"  off-targets: {g['casoff_total']}"
                rpt.write(f"      {i}. {g['sgRNA_id']:12s} score={g['Score']}  "
                           f"{g['Sequence']}  {g['strand']}  GC={g['%GC']}{casoff_str}\n")
            rpt.write("\n")

        # Cas-OFFinder summary (if data exists)
        has_casoff = any(g["casoff_total"] not in ("", "0") for g in guide_rows)
        if has_casoff:
            rpt.write(f"{'-' * 78}\n")
            rpt.write(f"  Cas-OFFinder Off-Target Summary\n")
            rpt.write(f"{'-' * 78}\n\n")
            rpt.write(f"  {'Guide':<14s} {'Score':>7s}  {'0mm':>4s} {'1mm':>4s} "
                       f"{'2mm':>4s} {'3mm':>4s} {'4mm':>4s} {'Total':>6s}\n")
            rpt.write(f"  {'-'*14} {'-'*7}  {'-'*4} {'-'*4} {'-'*4} {'-'*4} {'-'*4} {'-'*6}\n")
            for g in sorted(guide_rows, key=lambda x: float(x["Score"]), reverse=True):
                if g["casoff_total"] in ("", "0"):
                    continue
                rpt.write(f"  {g['sgRNA_id']:<14s} {g['Score']:>7s}  "
                           f"{g['casoff_0mm']:>4s} {g['casoff_1mm']:>4s} "
                           f"{g['casoff_2mm']:>4s} {g['casoff_3mm']:>4s} "
                           f"{g['casoff_4mm']:>4s} {g['casoff_total']:>6s}\n")
            rpt.write("\n")

        rpt.write("=" * 78 + "\n")
        rpt.write("  End of report\n")
        rpt.write("=" * 78 + "\n")

    print(f"  CRISPR_summary_report.txt -> {report_path}")


if __name__ == "__main__":
    main()
