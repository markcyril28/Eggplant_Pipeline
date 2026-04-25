#!/usr/bin/env python3
"""Op 7: Structural alignment between AlphaFold3 and SWISS-MODEL PDBs via PyMOL.

Superimposes each matched AF3/SWISS gene pair, records RMSD and per-residue
CA–CA distances.

Outputs per gene:
  Comparison_Results/{gene}/structural_alignment.csv   — per-residue distances
  Comparison_Results/{gene}/superposed_af3.pdb         — AF3 after superposition
  Comparison_Results/{gene}/superposed_swiss.pdb       — SWISS after superposition

Usage:
    python3 structural_alignment.py --run-dir /path/to/genome_dir
"""

import argparse
import csv
import math
import os
import re
from pathlib import Path

# Matches raw AlphaFold3 timestamped extract folders: YYYY_MM_DD_HH_MM[_SS]
_AF3_TIMESTAMP_RE = re.compile(r"^\d{4}_\d{2}_\d{2}_\d{2}_\d{2}")


def af3_to_swiss_name(af3_name: str) -> str:
    """Convert AF3 gene name to SWISS-MODEL folder name."""
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


def _init_pymol():
    """Start PyMOL in headless / quiet mode."""
    import pymol  # noqa: E402
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd  # noqa: E402
    return cmd


def align_gene(cmd, af3_pdb: Path, swiss_pdb: Path, out_dir: Path) -> dict:
    """Align one AF3/SWISS pair. Returns dict with RMSD and per-residue data."""
    cmd.reinitialize()

    cmd.load(str(af3_pdb), "af3")
    cmd.load(str(swiss_pdb), "swiss")

    # Superpose SWISS onto AF3 (super is robust to sequence differences)
    try:
        result = cmd.super("swiss", "af3")
        rmsd = result[0]
        n_aligned = result[1]
    except Exception as e:
        print(f"    PyMOL super failed: {e}")
        return {"error": str(e)}

    # Save superposed structures
    os.makedirs(out_dir, exist_ok=True)
    cmd.save(str(out_dir / "superposed_af3.pdb"), "af3")
    cmd.save(str(out_dir / "superposed_swiss.pdb"), "swiss")

    # Per-residue CA distances
    af3_cas = {}
    swiss_cas = {}

    af3_model = cmd.get_model("af3 and name CA")
    for atom in af3_model.atom:
        af3_cas[int(atom.resi)] = (atom.coord[0], atom.coord[1], atom.coord[2])

    swiss_model = cmd.get_model("swiss and name CA")
    for atom in swiss_model.atom:
        swiss_cas[int(atom.resi)] = (atom.coord[0], atom.coord[1], atom.coord[2])

    distances = []
    common_resi = sorted(set(af3_cas.keys()) & set(swiss_cas.keys()))
    for resi in common_resi:
        a = af3_cas[resi]
        s = swiss_cas[resi]
        dist = math.sqrt(sum((a[i] - s[i]) ** 2 for i in range(3)))
        distances.append((resi, dist))

    # Write per-residue alignment CSV
    csv_path = out_dir / "structural_alignment.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["residue", "ca_distance_angstrom"])
        for resi, dist in distances:
            writer.writerow([resi, f"{dist:.4f}"])

    return {
        "rmsd": rmsd,
        "n_aligned": n_aligned,
        "n_common_ca": len(common_resi),
        "mean_ca_dist": sum(d for _, d in distances) / len(distances) if distances else 0,
        "max_ca_dist": max(d for _, d in distances) if distances else 0,
    }


def align_all(run_dir: Path, overwrite: bool = False) -> int:
    """Align all matched gene pairs. Returns count of processed genes."""
    af3_dir = run_dir / "AlphaFold3_Results"
    swiss_dir = run_dir / "SWISS_Results"
    out_dir = run_dir / "Comparison_Results"

    if not af3_dir.is_dir() or not swiss_dir.is_dir():
        print(f"  Missing AF3 or SWISS directory in {run_dir.name}")
        return 0

    # Init PyMOL once
    cmd = _init_pymol()

    gene_dirs = sorted(
        d for d in af3_dir.iterdir()
        if d.is_dir() and d.name not in ("data", "__pycache__")
        and not _AF3_TIMESTAMP_RE.match(d.name)
    )

    # Collect alignment summary for the RMSD CSV
    summary_rows = []
    processed = 0

    for gene_dir in gene_dirs:
        gene_name = gene_dir.name
        swiss_name = af3_to_swiss_name(gene_name)

        af3_pdb = gene_dir / f"{gene_name}.pdb"
        swiss_gene_dir = swiss_dir / swiss_name
        swiss_pdbs = sorted(swiss_gene_dir.glob("*_model_*.pdb")) if swiss_gene_dir.is_dir() else []

        if not af3_pdb.exists() or not swiss_pdbs:
            print(f"  {gene_name}: no SWISS match for '{swiss_name}' — skipping")
            continue

        gene_out = out_dir / gene_name
        marker = gene_out / "structural_alignment.csv"
        if marker.exists() and not overwrite:
            print(f"  {gene_name}: alignment exists — skipping")
            summary_rows.append({"gene": gene_name, "status": "skipped"})
            processed += 1
            continue

        print(f"  Aligning: {gene_name}")
        result = align_gene(cmd, af3_pdb, swiss_pdbs[0], gene_out)

        if "error" in result:
            summary_rows.append({"gene": gene_name, "status": f"error: {result['error']}"})
        else:
            summary_rows.append({
                "gene": gene_name,
                "rmsd": f"{result['rmsd']:.4f}",
                "n_aligned": result["n_aligned"],
                "n_common_ca": result["n_common_ca"],
                "mean_ca_dist": f"{result['mean_ca_dist']:.4f}",
                "max_ca_dist": f"{result['max_ca_dist']:.4f}",
                "status": "ok",
            })
            processed += 1

    # Write alignment summary CSV
    if summary_rows:
        os.makedirs(out_dir, exist_ok=True)
        summary_csv = out_dir / "alignment_summary.csv"
        fieldnames = ["gene", "rmsd", "n_aligned", "n_common_ca", "mean_ca_dist", "max_ca_dist", "status"]
        with open(summary_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"  Wrote {summary_csv.relative_to(run_dir)}")

    return processed


def main():
    parser = argparse.ArgumentParser(description="Structural alignment AF3 vs SWISS")
    parser.add_argument("--run-dir", required=True, help="Genome run directory")
    parser.add_argument("--overwrite", default="false", help="Overwrite existing outputs")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    overwrite = args.overwrite.lower() in ("true", "1", "yes")

    n = align_all(run_dir, overwrite)
    print(f"  Aligned {n} gene(s)")


if __name__ == "__main__":
    main()
