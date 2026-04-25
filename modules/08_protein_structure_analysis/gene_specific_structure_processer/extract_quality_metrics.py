#!/usr/bin/env python3
"""Op 6: Extract quality metrics from AlphaFold3 and SWISS-MODEL outputs.

Parses pLDDT (AF3 PDB B-factors), ptm/ranking_score (AF3 JSON), QMEAN local
(SWISS PDB B-factors), and GMQE (SWISS PDB REMARK 3) for each gene.

Outputs:
  Comparison_Results/quality_metrics.csv          — one row per gene
  Comparison_Results/{gene}/per_residue_metrics.csv — per-residue pLDDT & QMEAN

Usage (standalone):
    python3 extract_quality_metrics.py --run-dir /path/to/genome_dir

Usage (orchestrated):
    Called by process_structures.sh
"""

import argparse
import csv
import json
import os
import re
from pathlib import Path

# Matches raw AlphaFold3 timestamped extract folders: YYYY_MM_DD_HH_MM[_SS]
_AF3_TIMESTAMP_RE = re.compile(r"^\d{4}_\d{2}_\d{2}_\d{2}_\d{2}")


# ── Gene-name mapping ───────────────────────────────────────────────────────

def af3_to_swiss_name(af3_name: str) -> str:
    """Convert AF3 gene name to SWISS-MODEL folder name.

    smel5_01g008730_1  →  SMEL5_01g008730.1
    """
    idx = af3_name.rfind("_")
    if idx == -1:
        return af3_name
    base = af3_name[:idx]
    version = af3_name[idx + 1 :]
    first_under = base.find("_")
    if first_under == -1:
        return base.upper() + "." + version
    prefix = base[:first_under].upper()
    return prefix + base[first_under:] + "." + version


# ── PDB B-factor extraction ────────────────────────────────────────────────

def extract_ca_bfactors(pdb_path: Path) -> list[tuple[int, str, float]]:
    """Extract per-residue B-factors from CA atoms.

    Returns list of (residue_number, residue_name, bfactor).
    """
    residues = []
    with open(pdb_path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            resi = int(line[22:26].strip())
            resn = line[17:20].strip()
            bfactor = float(line[60:66].strip())
            residues.append((resi, resn, bfactor))
    return residues


def parse_swiss_gmqe(pdb_path: Path) -> float | None:
    """Extract GMQE from SWISS-MODEL PDB REMARK 3."""
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("REMARK   3  GMQE"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    val = parts[-1]
                    if val != "NA":
                        return float(val)
    return None


def parse_swiss_sid(pdb_path: Path) -> float | None:
    """Extract sequence identity (SID) from SWISS-MODEL PDB REMARK 3."""
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("REMARK   3  SID"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    val = parts[-1]
                    if val != "NA":
                        return float(val)
    return None


# ── AF3 JSON confidence ────────────────────────────────────────────────────

def find_af3_json(run_dir: Path, gene_name: str, suffix: str) -> Path | None:
    """Find AF3 JSON file matching gene name under run_dir (including archive)."""
    pattern = re.compile(re.escape(gene_name) + r"_" + suffix.replace(".", r"\.") + "$")
    for root, _dirs, files in os.walk(str(run_dir)):
        for f in files:
            if pattern.search(f):
                return Path(root) / f
    # Fallback: broader search
    for root, _dirs, files in os.walk(str(run_dir)):
        for f in files:
            if gene_name in f and suffix in f:
                return Path(root) / f
    return None


def parse_af3_summary(json_path: Path) -> dict:
    """Parse AF3 summary_confidences JSON."""
    with open(json_path) as fh:
        data = json.load(fh)
    return {
        "ptm": data.get("ptm"),
        "ranking_score": data.get("ranking_score"),
        "fraction_disordered": data.get("fraction_disordered"),
    }


# ── Main ────────────────────────────────────────────────────────────────────

def extract_metrics(run_dir: Path, overwrite: bool = False) -> int:
    """Extract and tabulate quality metrics for all genes.  Returns gene count."""
    af3_dir = run_dir / "AlphaFold3_Results"
    swiss_dir = run_dir / "SWISS_Results"
    out_dir = run_dir / "Comparison_Results"

    if not af3_dir.is_dir():
        print(f"  AlphaFold3_Results not found in {run_dir.name}")
        return 0
    if not swiss_dir.is_dir():
        print(f"  SWISS_Results not found in {run_dir.name}")
        return 0

    summary_csv = out_dir / "quality_metrics.csv"
    if summary_csv.exists() and not overwrite:
        print(f"  quality_metrics.csv already exists — skipping (overwrite={overwrite})")
        return 0

    os.makedirs(out_dir, exist_ok=True)

    # Collect gene directories from AF3
    gene_dirs = sorted(
        d for d in af3_dir.iterdir()
        if d.is_dir() and d.name not in ("data", "__pycache__")
        and not _AF3_TIMESTAMP_RE.match(d.name)
    )

    summary_rows = []
    processed = 0

    for gene_dir in gene_dirs:
        gene_name = gene_dir.name
        swiss_name = af3_to_swiss_name(gene_name)

        # AF3 PDB
        af3_pdb = gene_dir / f"{gene_name}.pdb"
        if not af3_pdb.exists():
            print(f"  SKIP {gene_name}: AF3 PDB not found")
            continue

        # SWISS PDB
        swiss_gene_dir = swiss_dir / swiss_name
        swiss_pdbs = sorted(swiss_gene_dir.glob("*_model_*.pdb")) if swiss_gene_dir.is_dir() else []
        if not swiss_pdbs:
            print(f"  SKIP {gene_name}: SWISS PDB not found for {swiss_name}")
            continue
        swiss_pdb = swiss_pdbs[0]

        print(f"  Extracting metrics: {gene_name}")

        # Per-residue B-factors
        af3_residues = extract_ca_bfactors(af3_pdb)
        swiss_residues = extract_ca_bfactors(swiss_pdb)

        af3_plddt_map = {resi: bf for resi, _, bf in af3_residues}
        swiss_qmean_map = {resi: bf for resi, _, bf in swiss_residues}

        mean_plddt = sum(bf for _, _, bf in af3_residues) / len(af3_residues) if af3_residues else None
        mean_qmean = sum(bf for _, _, bf in swiss_residues) / len(swiss_residues) if swiss_residues else None

        # AF3 summary confidences JSON
        af3_summary = {}
        json_path = find_af3_json(run_dir, gene_name, "summary_confidences_0.json")
        if json_path:
            af3_summary = parse_af3_summary(json_path)

        # SWISS GMQE and SID
        swiss_gmqe = parse_swiss_gmqe(swiss_pdb)
        swiss_sid = parse_swiss_sid(swiss_pdb)

        # Per-residue CSV
        gene_out = out_dir / gene_name
        os.makedirs(gene_out, exist_ok=True)
        per_res_csv = gene_out / "per_residue_metrics.csv"

        all_resi = sorted(set(list(af3_plddt_map.keys()) + list(swiss_qmean_map.keys())))
        with open(per_res_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["residue", "af3_pLDDT", "swiss_QMEAN_local"])
            for resi in all_resi:
                plddt_val = af3_plddt_map.get(resi, "")
                qmean_val = swiss_qmean_map.get(resi, "")
                writer.writerow([resi, plddt_val, qmean_val])

        # Summary row
        summary_rows.append({
            "gene": gene_name,
            "swiss_gene": swiss_name,
            "af3_residues": len(af3_residues),
            "swiss_residues": len(swiss_residues),
            "af3_mean_pLDDT": f"{mean_plddt:.2f}" if mean_plddt else "",
            "af3_ptm": af3_summary.get("ptm", ""),
            "af3_ranking_score": af3_summary.get("ranking_score", ""),
            "af3_fraction_disordered": af3_summary.get("fraction_disordered", ""),
            "swiss_GMQE": swiss_gmqe if swiss_gmqe is not None else "",
            "swiss_SID": swiss_sid if swiss_sid is not None else "",
            "swiss_mean_QMEAN_local": f"{mean_qmean:.4f}" if mean_qmean else "",
        })
        processed += 1

    # Write summary CSV
    if summary_rows:
        with open(summary_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=summary_rows[0].keys())
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"  Wrote {summary_csv.relative_to(run_dir)}")

    return processed


def main():
    parser = argparse.ArgumentParser(description="Extract AF3/SWISS quality metrics")
    parser.add_argument("--run-dir", required=True, help="Genome run directory")
    parser.add_argument("--overwrite", default="false", help="Overwrite existing outputs")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    overwrite = args.overwrite.lower() in ("true", "1", "yes")

    n = extract_metrics(run_dir, overwrite)
    print(f"  Extracted metrics for {n} gene(s)")


if __name__ == "__main__":
    main()
