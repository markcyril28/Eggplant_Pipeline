#!/usr/bin/env python3
"""
Module: compare_tools.py
Stage [C] — CRISPOR vs CRISPR-P v2.0 side-by-side comparison.

Reads the final ranked-guide tables from the CRISPOR arm and the
CRISPR-P v2.0 arm (both produced by stage [8] in their respective
output trees) and generates:

  comparison_metrics.tsv       — per-guide metric table (merged, labelled)
  comparison_summary.tsv       — scalar summary statistics
  fig_score_scatter.png        — composite KO score scatter (CRISPOR vs CRISPR-P)
  fig_metric_boxplot.png       — side-by-side box plots for key metrics
  fig_venn_topN.png            — Venn diagram of top-N guide overlap
  fig_rank_correlation.png     — Spearman rank correlation heatmap

Usage:
    python3 compare_tools.py \\
        --crispor-ranked   <crispor_arm/08_Ranked_Guides/ranked_guides.tsv>  \\
        --crisprp-ranked   <crispr_p_arm/08_Ranked_Guides/ranked_guides.tsv> \\
        --outdir           <comparison/>                                      \\
        --metrics          doench2016OnTarget crisprOn_score best_fs_frac cfd_sum composite_ko_score \\
        --jaccard-top-n    10                                                 \\
        --dpi              600                                                \\
        --format           png
"""

import argparse
import csv
import math
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="CRISPOR vs CRISPR-P v2.0 comparison")
    p.add_argument("--crispor-ranked",  required=True)
    p.add_argument("--crisprp-ranked",  required=True)
    p.add_argument("--outdir",          required=True)
    p.add_argument("--metrics",         nargs="+",
                   default=["doench2016OnTarget", "crisprOn_score",
                             "best_fs_frac", "cfd_sum", "composite_ko_score"])
    p.add_argument("--jaccard-top-n",   type=int,   default=10)
    p.add_argument("--dpi",             type=int,   default=600)
    p.add_argument("--format",          default="png", choices=["png", "svg", "pdf"])
    p.add_argument("--overwrite",       action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


# ─── IO helpers ───────────────────────────────────────────────────────────────

def load_table(path: str) -> list[dict]:
    p = Path(path)
    if not p.exists():
        print(f"[compare] File not found: {path}", file=sys.stderr)
        return []
    delim = "," if path.endswith(".csv") else "\t"
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def guide_id(row: dict) -> str:
    for col in ("guideId", "guide_id", "name", "ID", "targetSeq"):
        if row.get(col):
            return str(row[col]).strip()
    return ""


def safe_float(val, default=float("nan")) -> float:
    if val in (None, "", "NA", "N/A"):
        return default
    try:
        f = float(val)
        return f if math.isfinite(f) else default
    except (ValueError, TypeError):
        return default


# ─── Summary statistics ───────────────────────────────────────────────────────

def summary_stats(rows: list[dict], metrics: list[str], label: str) -> list[dict]:
    out = []
    for m in metrics:
        vals = [safe_float(r.get(m)) for r in rows]
        vals = [v for v in vals if not math.isnan(v)]
        if not vals:
            out.append({"tool": label, "metric": m,
                        "n": 0, "mean": "NA", "median": "NA",
                        "min": "NA", "max": "NA", "stdev": "NA"})
            continue
        n    = len(vals)
        mean = sum(vals) / n
        srt  = sorted(vals)
        med  = srt[n // 2] if n % 2 else (srt[n // 2 - 1] + srt[n // 2]) / 2
        var  = sum((v - mean) ** 2 for v in vals) / n
        out.append({"tool": label, "metric": m,
                    "n": n,
                    "mean":   f"{mean:.4f}",
                    "median": f"{med:.4f}",
                    "min":    f"{srt[0]:.4f}",
                    "max":    f"{srt[-1]:.4f}",
                    "stdev":  f"{var**0.5:.4f}"})
    return out


def jaccard(set_a: set, set_b: set) -> float:
    if not set_a and not set_b:
        return 1.0
    return len(set_a & set_b) / len(set_a | set_b)


# ─── Plots ────────────────────────────────────────────────────────────────────

def make_plots(crispor_rows: list[dict], crisprp_rows: list[dict],
               metrics: list[str], top_n: int,
               outdir: Path, dpi: int, fmt: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[compare] matplotlib not installed — skipping plots.", file=sys.stderr)
        return

    # ── 1. Composite KO score scatter ────────────────────────────────────────
    crispor_scores = {guide_id(r): safe_float(r.get("composite_ko_score"))
                      for r in crispor_rows}
    crisprp_scores = {guide_id(r): safe_float(r.get("composite_ko_score"))
                      for r in crisprp_rows}
    common = set(crispor_scores) & set(crisprp_scores)
    if common:
        # Filter to guides where BOTH scores are non-NaN, building paired lists
        valid = [g for g in common
                 if not math.isnan(crispor_scores[g]) and not math.isnan(crisprp_scores[g])]
        xs = [crispor_scores[g] for g in valid]
        ys = [crisprp_scores[g] for g in valid]
        if xs:
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(xs, ys, alpha=0.6, edgecolors="none", s=25)
            lim = max(max(xs, default=1), max(ys, default=1)) * 1.05
            ax.plot([0, lim], [0, lim], "r--", lw=0.8, label="y = x")
            ax.set_xlabel("CRISPOR composite KO score")
            ax.set_ylabel("CRISPR-P v2.0 composite KO score")
            ax.set_title("Composite KO score: CRISPOR vs CRISPR-P v2.0")
            ax.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(outdir / f"fig_score_scatter.{fmt}", dpi=dpi)
            plt.close()

    # ── 2. Side-by-side box plots ─────────────────────────────────────────────
    plot_metrics = [m for m in metrics if m != "composite_ko_score"]
    if plot_metrics:
        n_m = len(plot_metrics)
        fig, axes = plt.subplots(1, n_m, figsize=(3 * n_m, 4), sharey=False)
        if n_m == 1:
            axes = [axes]
        for ax, m in zip(axes, plot_metrics):
            c_vals = [safe_float(r.get(m)) for r in crispor_rows]
            p_vals = [safe_float(r.get(m)) for r in crisprp_rows]
            c_vals = [v for v in c_vals if not math.isnan(v)]
            p_vals = [v for v in p_vals if not math.isnan(v)]
            ax.boxplot([c_vals, p_vals], labels=["CRISPOR", "CRISPR-P\nv2.0"],
                       patch_artist=True,
                       boxprops=dict(facecolor="#90CAF9"),
                       medianprops=dict(color="red"))
            ax.set_title(m.replace("_", "\n"), fontsize=8)
        plt.suptitle("Metric distributions: CRISPOR vs CRISPR-P v2.0", fontsize=9)
        plt.tight_layout()
        plt.savefig(outdir / f"fig_metric_boxplot.{fmt}", dpi=dpi)
        plt.close()

    # ── 3. Venn diagram of top-N guide overlap ────────────────────────────────
    top_crispor = {guide_id(r) for r in crispor_rows[:top_n]}
    top_crisprp = {guide_id(r) for r in crisprp_rows[:top_n]}
    only_c  = len(top_crispor - top_crisprp)
    only_p  = len(top_crisprp - top_crispor)
    both    = len(top_crispor & top_crisprp)
    jac     = jaccard(top_crispor, top_crisprp)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    try:
        from matplotlib_venn import venn2
        venn2(subsets=(only_c, only_p, both),
              set_labels=(f"CRISPOR\ntop-{top_n}", f"CRISPR-P v2.0\ntop-{top_n}"),
              ax=ax)
        ax.set_title(f"Top-{top_n} guide overlap  (Jaccard = {jac:.2f})")
    except ImportError:
        # Fallback: text-only summary
        ax.axis("off")
        ax.text(0.5, 0.5,
                f"CRISPOR only: {only_c}\n"
                f"CRISPR-P v2.0 only: {only_p}\n"
                f"Shared: {both}\n"
                f"Jaccard: {jac:.2f}",
                ha="center", va="center", fontsize=11,
                transform=ax.transAxes)
        ax.set_title(f"Top-{top_n} guide overlap")
    plt.tight_layout()
    plt.savefig(outdir / f"fig_venn_top{top_n}.{fmt}", dpi=dpi)
    plt.close()

    # ── 4. Spearman rank correlation heatmap ──────────────────────────────────
    try:
        from scipy.stats import spearmanr
        import numpy as np

        all_rows = ([(r, "CRISPOR") for r in crispor_rows] +
                    [(r, "CRISPR-P") for r in crisprp_rows])
        corr_metrics = [m for m in metrics
                        if any(safe_float(r.get(m)) == safe_float(r.get(m))
                               for r, _ in all_rows)]
        if len(corr_metrics) >= 2:
            mat = np.array([[safe_float(r.get(m), 0.0)
                             for m in corr_metrics]
                            for r, _ in all_rows])
            n_m2 = len(corr_metrics)
            corr_matrix = np.zeros((n_m2, n_m2))
            for i in range(n_m2):
                for j in range(n_m2):
                    rho, _ = spearmanr(mat[:, i], mat[:, j])
                    corr_matrix[i, j] = rho

            fig, ax = plt.subplots(figsize=(max(4, n_m2), max(4, n_m2)))
            im = ax.imshow(corr_matrix, vmin=-1, vmax=1, cmap="RdBu_r")
            ax.set_xticks(range(n_m2))
            ax.set_yticks(range(n_m2))
            labels = [m.replace("_", "\n") for m in corr_metrics]
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
            ax.set_yticklabels(labels, fontsize=8)
            for i in range(n_m2):
                for j in range(n_m2):
                    ax.text(j, i, f"{corr_matrix[i,j]:.2f}",
                            ha="center", va="center", fontsize=7)
            plt.colorbar(im, ax=ax, label="Spearman ρ")
            ax.set_title("Metric rank correlations (both tools combined)")
            plt.tight_layout()
            plt.savefig(outdir / f"fig_rank_correlation.{fmt}", dpi=dpi)
            plt.close()
    except ImportError:
        print("[compare] scipy/numpy not available — skipping correlation heatmap.",
              file=sys.stderr)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_path = outdir / "comparison_summary.tsv"
    metrics_path = outdir / "comparison_metrics.tsv"

    if not args.overwrite and summary_path.exists():
        print(f"[compare] Skipping (overwrite=false): {summary_path}")
        return

    crispor_rows = load_table(args.crispor_ranked)
    crisprp_rows = load_table(args.crisprp_ranked)

    if not crispor_rows:
        print("[compare] CRISPOR ranked table is empty — aborting.", file=sys.stderr)
        return
    if not crisprp_rows:
        print("[compare] CRISPR-P v2.0 ranked table is empty — results may be incomplete.",
              file=sys.stderr)

    # ── Merged per-guide table ────────────────────────────────────────────────
    all_fields = (list(next(iter(crispor_rows), {}).keys()) +
                  list(next(iter(crisprp_rows), {}).keys()))
    # Deduplicate while preserving order
    seen: set = set()
    merged_fields = []
    for f in all_fields:
        if f not in seen:
            seen.add(f)
            merged_fields.append(f)
    merged_fields = ["tool_source"] + merged_fields

    merged: list[dict] = []
    for r in crispor_rows:
        merged.append({"tool_source": "CRISPOR", **r})
    for r in crisprp_rows:
        merged.append({"tool_source": "CRISPR-P_v2.0", **r})

    with open(metrics_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=merged_fields,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    # ── Summary stats ─────────────────────────────────────────────────────────
    stats = (summary_stats(crispor_rows, args.metrics, "CRISPOR") +
             summary_stats(crisprp_rows, args.metrics, "CRISPR-P_v2.0"))

    top_c = {guide_id(r) for r in crispor_rows[:args.jaccard_top_n]}
    top_p = {guide_id(r) for r in crisprp_rows[:args.jaccard_top_n]}
    jac   = jaccard(top_c, top_p)

    stats.append({"tool": "comparison",
                  "metric": f"jaccard_top{args.jaccard_top_n}",
                  "n":      len(top_c | top_p),
                  "mean":   f"{jac:.4f}",
                  "median": "NA", "min": "NA", "max": "NA", "stdev": "NA"})

    stat_fields = ["tool", "metric", "n", "mean", "median", "min", "max", "stdev"]
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=stat_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(stats)

    # ── Figures ───────────────────────────────────────────────────────────────
    make_plots(crispor_rows, crisprp_rows, args.metrics,
               args.jaccard_top_n, outdir, args.dpi, args.format)

    print(f"[compare] CRISPOR guides: {len(crispor_rows)}, "
          f"CRISPR-P v2.0 guides: {len(crisprp_rows)}")
    print(f"[compare] Top-{args.jaccard_top_n} Jaccard overlap: {jac:.3f}")
    print(f"[compare] Outputs -> {outdir}")


if __name__ == "__main__":
    main()
