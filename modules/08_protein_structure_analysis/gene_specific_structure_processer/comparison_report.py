#!/usr/bin/env python3
"""Op 9: Aggregate comparison metrics and generate summary report.

Reads quality_metrics.csv and alignment_summary.csv produced by Ops 6-7,
merges them into a single report CSV, and generates matplotlib figures.

Outputs:
  Comparison_Results/comparison_report.csv
  Comparison_Results/figures/rmsd_per_gene.jpg
  Comparison_Results/figures/plddt_vs_gmqe.jpg
  Comparison_Results/figures/residue_comparison_{gene}.jpg (per gene)

Usage:
    python3 comparison_report.py --run-dir /path/to/genome_dir
"""

import argparse
import csv
import os
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def build_report(run_dir: Path, overwrite: bool = False) -> int:
    comp_dir = run_dir / "Comparison_Results"
    quality_csv = comp_dir / "quality_metrics.csv"
    alignment_csv = comp_dir / "alignment_summary.csv"
    report_csv = comp_dir / "comparison_report.csv"

    if report_csv.exists() and not overwrite:
        print(f"  comparison_report.csv exists — skipping (overwrite={overwrite})")
        return 0

    if not quality_csv.exists():
        print("  quality_metrics.csv not found — run Op 6 first")
        return 0

    # Load quality metrics
    quality_rows = {}
    with open(quality_csv) as fh:
        for row in csv.DictReader(fh):
            quality_rows[row["gene"]] = row

    # Load alignment summary (optional)
    align_rows = {}
    if alignment_csv.exists():
        with open(alignment_csv) as fh:
            for row in csv.DictReader(fh):
                align_rows[row["gene"]] = row

    # Merge into report
    report = []
    for gene in sorted(quality_rows.keys()):
        q = quality_rows[gene]
        a = align_rows.get(gene, {})
        report.append({
            "gene": gene,
            "swiss_gene": q.get("swiss_gene", ""),
            "af3_residues": q.get("af3_residues", ""),
            "swiss_residues": q.get("swiss_residues", ""),
            "af3_mean_pLDDT": q.get("af3_mean_pLDDT", ""),
            "af3_ptm": q.get("af3_ptm", ""),
            "af3_ranking_score": q.get("af3_ranking_score", ""),
            "af3_fraction_disordered": q.get("af3_fraction_disordered", ""),
            "swiss_GMQE": q.get("swiss_GMQE", ""),
            "swiss_SID": q.get("swiss_SID", ""),
            "swiss_mean_QMEAN_local": q.get("swiss_mean_QMEAN_local", ""),
            "rmsd": a.get("rmsd", ""),
            "n_aligned": a.get("n_aligned", ""),
            "mean_ca_dist": a.get("mean_ca_dist", ""),
            "max_ca_dist": a.get("max_ca_dist", ""),
        })

    if not report:
        print("  No genes in quality_metrics.csv — nothing to report")
        return 0

    os.makedirs(comp_dir, exist_ok=True)
    with open(report_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=report[0].keys())
        writer.writeheader()
        writer.writerows(report)
    print(f"  Wrote {report_csv.relative_to(run_dir)}")

    # Figures
    if not HAS_MPL:
        print("  matplotlib not available — skipping figures")
        return len(report)

    fig_dir = comp_dir / "figures"
    os.makedirs(fig_dir, exist_ok=True)

    genes = [r["gene"] for r in report]
    short_labels = [g.replace("smel5_", "") for g in genes]

    # ── Figure 1: RMSD per gene ──────────────────────────────────────────
    rmsd_genes, rmsd_labels, rmsds = [], [], []
    for r, lbl in zip(report, short_labels):
        try:
            val = float(r["rmsd"])
            rmsd_genes.append(r["gene"])
            rmsd_labels.append(lbl)
            rmsds.append(val)
        except (ValueError, KeyError):
            pass  # skip genes with no alignment data

    if rmsds:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(rmsd_labels, rmsds, color="#4a90d9", edgecolor="black", linewidth=0.5)
        ax.set_ylabel("RMSD (Å)")
        ax.set_title("AF3 vs SWISS-MODEL — Global RMSD per Gene")
        ax.set_xlabel("Gene")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.tight_layout()
        fig_path = fig_dir / "rmsd_per_gene.jpg"
        fig.savefig(str(fig_path), dpi=300, format="jpeg")
        plt.close(fig)
        print(f"  Wrote {fig_path.relative_to(run_dir)}")

    # ── Figure 2: pLDDT vs GMQE scatter ─────────────────────────────────
    plddts, gmqes, labels = [], [], []
    for r in report:
        try:
            p = float(r["af3_mean_pLDDT"])
            g = float(r["swiss_GMQE"])
            plddts.append(p)
            gmqes.append(g)
            labels.append(r["gene"].replace("smel5_", ""))
        except (ValueError, KeyError):
            pass

    if plddts:
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(gmqes, plddts, s=60, color="#e74c3c", edgecolors="black", linewidth=0.5, zorder=3)
        for i, lbl in enumerate(labels):
            ax.annotate(lbl, (gmqes[i], plddts[i]), fontsize=7,
                        xytext=(5, 5), textcoords="offset points")
        ax.set_xlabel("SWISS-MODEL GMQE")
        ax.set_ylabel("AlphaFold3 Mean pLDDT")
        ax.set_title("Model Confidence Comparison")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fig_path = fig_dir / "plddt_vs_gmqe.jpg"
        fig.savefig(str(fig_path), dpi=300, format="jpeg")
        plt.close(fig)
        print(f"  Wrote {fig_path.relative_to(run_dir)}")

    # ── Figure 3: Per-residue comparison plots (one per gene) ────────────
    for r in report:
        gene = r["gene"]
        res_csv = comp_dir / gene / "per_residue_metrics.csv"
        align_csv = comp_dir / gene / "structural_alignment.csv"
        if not res_csv.exists():
            continue

        residues, af3_vals, swiss_vals = [], [], []
        with open(res_csv) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                resi = int(row["residue"])
                af3_v = float(row["af3_pLDDT"]) if row["af3_pLDDT"] else None
                sw_v = float(row["swiss_QMEAN_local"]) if row["swiss_QMEAN_local"] else None
                if af3_v is not None and sw_v is not None:
                    residues.append(resi)
                    af3_vals.append(af3_v)
                    swiss_vals.append(sw_v * 100)  # scale to 0-100 for comparison

        ca_dists = {}
        if align_csv.exists():
            with open(align_csv) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    ca_dists[int(row["residue"])] = float(row["ca_distance_angstrom"])

        if not residues:
            continue

        n_panels = 3 if ca_dists else 2
        fig, axes = plt.subplots(n_panels, 1, figsize=(12, 3 * n_panels), sharex=True)
        if n_panels == 2:
            axes = list(axes)

        # Panel 1: pLDDT
        axes[0].plot(residues, af3_vals, color="#0053d6", linewidth=0.8, label="AF3 pLDDT")
        axes[0].axhline(70, color="grey", linestyle="--", linewidth=0.5, alpha=0.5)
        axes[0].set_ylabel("pLDDT (0-100)")
        axes[0].set_title(f"{gene} — Per-Residue Confidence Comparison")
        axes[0].legend(fontsize=8)
        axes[0].set_ylim(0, 100)

        # Panel 2: QMEAN local (scaled to 0-100)
        axes[1].plot(residues, swiss_vals, color="#e74c3c", linewidth=0.8, label="SWISS QMEAN×100")
        axes[1].axhline(70, color="grey", linestyle="--", linewidth=0.5, alpha=0.5)
        axes[1].set_ylabel("QMEAN local ×100")
        axes[1].legend(fontsize=8)
        axes[1].set_ylim(0, 100)

        # Panel 3: CA distance
        if n_panels == 3:
            dist_resi = sorted(ca_dists.keys())
            dist_vals = [ca_dists[r] for r in dist_resi]
            axes[2].fill_between(dist_resi, dist_vals, color="#f39c12", alpha=0.6)
            axes[2].plot(dist_resi, dist_vals, color="#e67e22", linewidth=0.6)
            axes[2].set_ylabel("CA Distance (Å)")
            axes[2].set_xlabel("Residue Number")

        plt.tight_layout()
        fig_path = fig_dir / f"residue_comparison_{gene}.jpg"
        fig.savefig(str(fig_path), dpi=200, format="jpeg")
        plt.close(fig)
        print(f"  Wrote {fig_path.relative_to(run_dir)}")

    return len(report)


def main():
    parser = argparse.ArgumentParser(description="Comparison summary report")
    parser.add_argument("--run-dir", required=True, help="Genome run directory")
    parser.add_argument("--overwrite", default="false", help="Overwrite existing outputs")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    overwrite = args.overwrite.lower() in ("true", "1", "yes")

    n = build_report(run_dir, overwrite)
    print(f"  Report covers {n} gene(s)")


if __name__ == "__main__":
    main()
