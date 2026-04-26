#!/usr/bin/env python3
"""
CRISPR Off-Target Analysis: Publication-Quality Visualizations

Reads guide_summary.csv (from generate_report.py) and produces:
  1. guide_comparison.png         : score vs off-target count scatter + sequence table
  2. guide_comparison_scatter.png : scatter only
  3. guide_comparison_table.png   : sequence table only

Usage:
    python3 plot_crispr_results.py --summary-csv <guide_summary.csv> --output-dir <dir> \
        [--dpi 300] [--format png] [--y-max-cap 100]
"""

import argparse
import csv
import os
import sys

try:
    import matplotlib
    matplotlib.use("Agg")  # headless rendering
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ── Colourblind-safe palette (matches pipeline motif convention) ───────────
TIER_COLORS = {
    "High": "#2ca02c",      # green
    "Moderate": "#ff7f0e",   # orange
    "Low": "#d62728",        # red
}
GENE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
               "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]


def read_summary(path: str) -> list[dict]:
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


# ── Plot: Score vs off-target count scatter ───────────────────────────────

def plot_guide_comparison(guides: list[dict], output_dir: str, dpi: int, fmt: str,
                          y_max_cap: float = 0.0,
                          score_strict: float = 0.7,
                          score_moderate: float = 0.5,
                          offtarget_strict: int = 0,
                          offtarget_moderate: int = 10):
    """Scatter plot: on-target score (x) vs Cas-OFFinder off-target count (y).

    Saves three files:
      guide_comparison.{fmt}         : combined scatter + sequence table
      guide_comparison_scatter.{fmt} : scatter only
      guide_comparison_table.{fmt}   : sequence table only
    """
    data = []
    for g in guides:
        casoff = g.get("casoff_total", "")
        if casoff and casoff != "0":
            data.append({
                "score":      float(g["Score"]),
                "offtargets": int(casoff),
                "gene":       g["gene"].split()[0],
                "tier":       g["Tier"],
                "label":      g["sgRNA_id"],
                "sequence":   g.get("Sequence", ""),
                "gc":         g.get("%GC", ""),
                "strand":     g.get("strand", ""),
            })

    if not data:
        return  # no off-target data available

    SCORE_HIGH = score_strict
    SCORE_MOD  = score_moderate

    # Off-target threshold: explicit override if > 0, else dynamic
    # (25th percentile of observed off-target counts, floored at 5)
    offtargets_sorted = sorted(d["offtargets"] for d in data)
    if offtarget_strict and offtarget_strict > 0:
        offtarget_thresh = int(offtarget_strict)
    else:
        offtarget_thresh = max(5, offtargets_sorted[len(offtargets_sorted) // 4])
    raw_max = max(d["offtargets"] for d in data)
    ymax_auto = raw_max * 1.15
    # Capped variant: clamp to y_max_cap when configured (0 means no cap → same as auto)
    if y_max_cap > 0:
        over_cap = sum(1 for d in data if d["offtargets"] > y_max_cap)
        if over_cap:
            print(f"  Warning: {over_cap} guide(s) exceed y_max_cap={y_max_cap:g} (observed max={raw_max})")
        ymax_capped = float(y_max_cap)
    else:
        ymax_capped = ymax_auto
    # Default ymax (used by combined fig + table layout) is the capped variant
    ymax = ymax_capped

    # Candidate pool for the comparison tables: ALL guides, grouped by gene
    # and sorted by on-target score descending. Inclusion in the strict and
    # moderate tables is controlled solely by the score/off-target thresholds
    # declared in the TOML config (no hidden top-N cap).
    from collections import defaultdict as _dd
    _per_gene = _dd(list)
    for d in sorted(data, key=lambda d: d["score"], reverse=True):
        _per_gene[d["gene"]].append(d)
    best = []
    for _gene in sorted(_per_gene.keys()):
        best.extend(_per_gene[_gene])
    has_table = bool(best)

    gene_set       = sorted(set(d["gene"] for d in data))
    gene_color_map = {g: GENE_COLORS[i % len(GENE_COLORS)] for i, g in enumerate(gene_set)}

    # ── Shared helpers ────────────────────────────────────────────────────
    def _lighten(hex_color: str, alpha: float = 0.18) -> tuple:
        import matplotlib.colors as mc
        r, g, b = mc.to_rgb(hex_color)
        return (1 - alpha + alpha * r,
                1 - alpha + alpha * g,
                1 - alpha + alpha * b)

    def _draw_scatter(fig, ax, ymax_local=None):
        if ymax_local is None:
            ymax_local = ymax
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0, ymax_local)

        # Ideal-quadrant shading
        ax.axhspan(0, offtarget_thresh, xmin=SCORE_HIGH, xmax=1.0,
                   color="#2ca02c", alpha=0.08, zorder=0, label="_nolegend_")

        # Threshold lines
        ax.axvline(x=SCORE_HIGH, color="#2ca02c", linestyle="--", linewidth=1.3,
                   alpha=0.80, zorder=1)
        ax.axvline(x=SCORE_MOD,  color="#ff7f0e", linestyle=":",  linewidth=1.1,
                   alpha=0.70, zorder=1)
        ax.axhline(y=offtarget_thresh, color="#1f77b4", linestyle="--", linewidth=1.1,
                   alpha=0.70, zorder=1)

        ax.text(SCORE_HIGH + 0.012, ymax_local * 0.96,
                f"Score ≥ {SCORE_HIGH}\n(High)",
                fontsize=7.5, color="#2ca02c", va="top", alpha=0.95)
        ax.text(SCORE_MOD + 0.012, ymax_local * 0.96,
                f"Score ≥ {SCORE_MOD}\n(Moderate)",
                fontsize=7.5, color="#ff7f0e", va="top", alpha=0.90)
        ax.text(0.012, offtarget_thresh + ymax_local * 0.012,
                f"≤ {offtarget_thresh} off-targets",
                fontsize=7.5, color="#1f77b4", va="bottom", alpha=0.90)

        # Directional arrows
        ax.annotate("", xy=(0.97, -0.055), xytext=(0.75, -0.055),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=1.3),
                    annotation_clip=False)
        ax.text(0.86, -0.072, "More on-target specific →",
                transform=ax.transAxes, fontsize=8, color="#555555",
                ha="center", va="top", style="italic")

        ax.annotate("", xy=(-0.085, 0.04), xytext=(-0.085, 0.26),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=1.3),
                    annotation_clip=False)
        ax.text(-0.10, 0.15, "Fewer off-targets ↓",
                transform=ax.transAxes, fontsize=8, color="#555555",
                ha="center", va="center", style="italic", rotation=90)

        # Ideal-region label
        ideal_x = (SCORE_HIGH + 1.0) / 2.0
        ideal_y = offtarget_thresh * 0.45
        ax.text(ideal_x, ideal_y, "Ideal\nguides", fontsize=8.5,
                color="#2ca02c", ha="center", va="center", alpha=0.75,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.35", fc="white",
                          ec="#2ca02c", alpha=0.45, lw=0.9))

        # Scatter points (clamp y to ymax_local so over-cap guides render at top edge)
        for d in data:
            c      = gene_color_map[d["gene"]]
            marker = "^" if d["tier"] == "High" else ("o" if d["tier"] == "Moderate" else "s")
            y_plot = min(d["offtargets"], ymax_local)
            ax.scatter(d["score"], y_plot, c=c, marker=marker, s=70,
                       edgecolors="black", linewidth=0.5, zorder=4)
            lbl = d["label"]
            if d["offtargets"] > ymax_local:
                lbl += f" (>{int(ymax_local)})"
            ax.annotate(lbl, (d["score"], y_plot),
                        fontsize=6, alpha=0.75, xytext=(3, 3),
                        textcoords="offset points")

        ax.set_xlabel("On-Target Score", fontsize=11)
        ax.set_ylabel("Off-Target Hits (Cas-OFFinder)", fontsize=11)
        ax.set_title("Guide Selection: Score vs Off-Target Count",
                     fontsize=13, fontweight="bold")

        # Legend (upper left, avoids the 0.7 threshold line)
        gene_patches = [mpatches.Patch(color=c, label=g) for g, c in gene_color_map.items()]
        tier_handles = [
            plt.scatter([], [], marker="^", c="gray", s=50, label="High",
                        edgecolors="black", linewidth=0.4),
            plt.scatter([], [], marker="o", c="gray", s=50, label="Moderate",
                        edgecolors="black", linewidth=0.4),
            plt.scatter([], [], marker="s", c="gray", s=50, label="Low",
                        edgecolors="black", linewidth=0.4),
        ]
        thresh_handles = [
            plt.Line2D([0], [0], color="#2ca02c", linestyle="--", lw=1.3,
                       label=f"High score ≥ {SCORE_HIGH}"),
            plt.Line2D([0], [0], color="#ff7f0e", linestyle=":",  lw=1.1,
                       label=f"Moderate score ≥ {SCORE_MOD}"),
            plt.Line2D([0], [0], color="#1f77b4", linestyle="--", lw=1.1,
                       label=f"Off-target limit ({offtarget_thresh})"),
        ]
        ax.legend(handles=gene_patches + tier_handles + thresh_handles,
                  fontsize=7, loc="upper left", ncol=2, framealpha=0.92)

    # Threshold-tier palette (matches scatter line colors)
    STRICT_COLOR   = "#2ca02c"   # green: passes strict thresholds
    MODERATE_COLOR = "#ff7f0e"   # orange: passes moderate but not strict

    def _row_tier_color(row_d):
        """Return hex color for a row based on which threshold tier it satisfies.

        Strict = score >= SCORE_HIGH AND off-targets <= offtarget_thresh.
        Otherwise the row is rendered with the moderate-tier color.
        """
        passes_strict = (row_d["score"] >= SCORE_HIGH
                         and row_d["offtargets"] <= offtarget_thresh)
        return STRICT_COLOR if passes_strict else MODERATE_COLOR

    def _draw_table(fig, ax_t, variant="strict"):
        # variant: "strict" -> green header/title; "moderate" -> orange.
        header_color = STRICT_COLOR if variant == "strict" else MODERATE_COLOR
        title = ("Guides Passing Strict Thresholds  (ranked by on-target score)"
                 if variant == "strict"
                 else "Guides Passing Moderate Thresholds  (ranked by on-target score)")

        col_labels = ["Gene", "Guide", "Score", "Off-targets", "Strand", "%GC",
                      "Sequence  (5'→3')"]
        table_rows = [
            [
                d["gene"].split(".")[0],
                d["label"],
                f"{d['score']:.4f}",
                str(d["offtargets"]),
                d["strand"],
                d["gc"],
                d["sequence"],
            ]
            for d in best
        ]

        ax_t.axis("off")
        ax_t.set_xlim(0, 1)
        ax_t.set_ylim(0, 1)
        tbl = ax_t.table(
            cellText=table_rows,
            colLabels=col_labels,
            cellLoc="left",
            loc="upper center",
            bbox=[0.0, 0.0, 1.0, 1.0],   # fill the axes; no internal whitespace
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8.5)
        tbl.auto_set_column_width(list(range(len(col_labels))))

        # Header row: no fill, bold text only
        for j in range(len(col_labels)):
            cell = tbl[0, j]
            cell.set_facecolor("white")
            cell.set_text_props(fontweight="bold")

        # Data rows: leave most cells unfilled; tint only the Score (col 2)
        # and Off-targets (col 3) cells by which threshold tier the row passes.
        SCORE_COL_IDX     = 2
        OFFTARGET_COL_IDX = 3
        for i, row_d in enumerate(best, start=1):
            tier_hex  = _row_tier_color(row_d)
            cell_tint = _lighten(tier_hex, alpha=0.22)
            for j in range(len(col_labels)):
                tbl[i, j].set_facecolor("white")
            tbl[i, SCORE_COL_IDX].set_facecolor(cell_tint)
            tbl[i, OFFTARGET_COL_IDX].set_facecolor(cell_tint)
            tbl[i, 6].set_text_props(family="monospace")

        # Gene-group demarcation: thicken the border between rows whose Gene
        # column differs from the previous row (rows are clustered by gene).
        n_cols = len(col_labels)
        for i in range(1, len(best)):
            if best[i]["gene"] != best[i - 1]["gene"]:
                # tbl[i, *] is the last guide of the previous gene; tbl[i+1, *]
                # is the first of the new gene. Thicken the shared edge.
                for j in range(n_cols):
                    tbl[i,     j].set_linewidth(1.8)
                    tbl[i + 1, j].set_linewidth(1.8)
                    tbl[i,     j].set_edgecolor("#222222")
                    tbl[i + 1, j].set_edgecolor("#222222")

        # Title above the table (axes coords; va=bottom places it above y=1.0)
        ax_t.text(0.5, 1.13, title,
                  transform=ax_t.transAxes,
                  ha="center", va="bottom", fontsize=10,
                  fontweight="bold", color=header_color, clip_on=False)

        # Tier legend: explains the Score / Off-targets cell tints
        legend_handles = [
            mpatches.Patch(
                facecolor=_lighten(STRICT_COLOR, alpha=0.22),
                edgecolor="#666666", linewidth=0.6,
                label=f"Strict: score ≥ {SCORE_HIGH:g}, "
                      f"off-targets ≤ {offtarget_thresh:g}",
            ),
            mpatches.Patch(
                facecolor=_lighten(MODERATE_COLOR, alpha=0.22),
                edgecolor="#666666", linewidth=0.6,
                label=f"Moderate: score ≥ {SCORE_MOD:g}, "
                      f"off-targets ≤ {int(offtarget_moderate):g}",
            ),
        ]
        ax_t.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=2,
            fontsize=7.5,
            frameon=False,
            handlelength=1.4,
            handleheight=1.0,
            columnspacing=1.4,
        )

    # ── Output directory layout ───────────────────────────────────────────
    #   <output_dir>/
    #     guide_comparison.{fmt}              (combined scatter + table)
    #     scatter/
    #       guide_comparison_scatter_capped.{fmt}     (y-axis clamped to y_max_cap)
    #       guide_comparison_scatter_uncapped.{fmt}   (auto-scaled to data range)
    #     tables/
    #       guide_comparison_table.{png,csv,tsv}            (strict)
    #       guide_comparison_table_moderate.{png,csv,tsv}   (moderate)
    scatter_dir = os.path.join(output_dir, "scatter")
    tables_dir  = os.path.join(output_dir, "tables")
    os.makedirs(scatter_dir, exist_ok=True)
    os.makedirs(tables_dir,  exist_ok=True)

    # ── Combined (scatter + table) ─────────────────────────────────────────
    # Embed only the strict-passing shortlist so the figure doesn't blow up.
    _all_candidates = list(best)
    strict_passing = [d for d in _all_candidates
                      if d["score"] >= SCORE_HIGH and d["offtargets"] <= offtarget_thresh]
    fig_c, ax_c = plt.subplots(figsize=(11, 11))
    fig_c.subplots_adjust(bottom=0.50, top=0.94, left=0.10, right=0.95)
    _draw_scatter(fig_c, ax_c, ymax_capped)
    if strict_passing:
        best = strict_passing
        n_rows = len(best)
        ax_t_c = fig_c.add_axes([0.06, 0.10, 0.90, min(0.06 + n_rows * 0.055, 0.26)])
        _draw_table(fig_c, ax_t_c, variant="strict")
        best = _all_candidates
    path_c = os.path.join(output_dir, f"guide_comparison.{fmt}")
    fig_c.savefig(path_c, dpi=dpi, bbox_inches="tight")
    plt.close(fig_c)
    print(f"  guide_comparison.{fmt} -> {path_c}")

    # ── Scatter only: two y-axis variants ──────────────────────────────────
    scatter_variants = [
        ("capped",   ymax_capped, "y-axis clamped to y_max_cap"),
        ("uncapped", ymax_auto,   "y-axis auto-scaled to data range"),
    ]
    for v_label, v_ymax, v_desc in scatter_variants:
        fig_s, ax_s = plt.subplots(figsize=(11, 7))
        fig_s.subplots_adjust(bottom=0.12, top=0.94, left=0.12, right=0.95)
        _draw_scatter(fig_s, ax_s, v_ymax)
        path_s = os.path.join(scatter_dir,
                              f"guide_comparison_scatter_{v_label}.{fmt}")
        fig_s.savefig(path_s, dpi=dpi, bbox_inches="tight")
        plt.close(fig_s)
        print(f"  scatter/guide_comparison_scatter_{v_label}.{fmt} -> {path_s} "
              f"({v_desc}, ymax={v_ymax:.1f})")

    # ── Table only ─────────────────────────────────────────────────────────
    # Two filtered table variants:
    #   strict:   score >= SCORE_HIGH AND off-targets <= offtarget_thresh
    #             (the green "ideal" region)
    #   moderate: score >= SCORE_MOD  AND off-targets <= offtarget_moderate
    #             (a more graceful inclusion that matches the orange "Moderate"
    #             line on the scatter; surfaces near-miss candidates without
    #             flooding the table)
    table_variants = [
        {
            "suffix":      "",
            "label":       "strict",
            "score_thr":   SCORE_HIGH,
            "offt_thr":    offtarget_thresh,
        },
        {
            "suffix":      "_moderate",
            "label":       "moderate",
            "score_thr":   SCORE_MOD,
            "offt_thr":    int(offtarget_moderate),
        },
    ]

    header = ["Gene", "Guide", "Score", "Off-targets", "Strand", "%GC",
              "Sequence_5to3", "Tier",
              "Score_threshold", "Offtarget_threshold"]

    # Cache the full candidate list so we can restore it between variants
    best_for_table = list(best)

    for variant in table_variants:
        v_suffix     = variant["suffix"]
        v_label      = variant["label"]
        v_score_thr  = variant["score_thr"]
        v_offt_thr   = variant["offt_thr"]

        passing = [d for d in best_for_table
                   if d["score"] >= v_score_thr and d["offtargets"] <= v_offt_thr]

        rows = [
            [
                d["gene"].split(".")[0],
                d["label"],
                f"{d['score']:.4f}",
                d["offtargets"],
                d["strand"],
                d["gc"],
                d["sequence"],
                d["tier"],
                v_score_thr,
                v_offt_thr,
            ]
            for d in passing
        ]
        for ext, delim in (("csv", ","), ("tsv", "\t")):
            out_path = os.path.join(tables_dir, f"guide_comparison_table{v_suffix}.{ext}")
            with open(out_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter=delim)
                writer.writerow(header)
                writer.writerows(rows)
            print(f"  tables/guide_comparison_table{v_suffix}.{ext} -> {out_path} "
                  f"({len(passing)} guide(s) passed {v_label} thresholds: "
                  f"score>={v_score_thr}, off-targets<={v_offt_thr})")

        if passing:
            # _draw_table renders from the outer-scope `best`, so swap it in
            best = passing
            n_rows     = len(best)
            fig_height = max(3.0, 1.2 + n_rows * 0.42)
            fig_t, ax_t = plt.subplots(figsize=(11, fig_height))
            fig_t.subplots_adjust(top=0.82, bottom=0.01, left=0.005, right=0.995)
            _draw_table(fig_t, ax_t, variant=v_label)
            path_t = os.path.join(tables_dir, f"guide_comparison_table{v_suffix}.{fmt}")
            fig_t.savefig(path_t, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
            plt.close(fig_t)
            print(f"  tables/guide_comparison_table{v_suffix}.{fmt} -> {path_t} "
                  f"({len(passing)} of {len(best_for_table)} guides passed "
                  f"{v_label} thresholds)")
            best = best_for_table
        else:
            print(f"  tables/guide_comparison_table{v_suffix}.{fmt} -> skipped "
                  f"(no guides crossed {v_label} thresholds: "
                  f"score>={v_score_thr} AND off-targets<={v_offt_thr})")



# ── main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Plot CRISPR analysis results")
    parser.add_argument("--summary-csv", required=True, help="guide_summary.csv from generate_report.py")
    parser.add_argument("--output-dir", required=True, help="Directory for plot output")
    parser.add_argument("--dpi", type=int, default=300, help="Plot resolution (default 300)")
    parser.add_argument("--format", default="png", choices=["png", "svg", "pdf"],
                        help="Output format (default png)")
    parser.add_argument("--y-max-cap", type=float, default=0.0,
                        help="Upper bound for off-target y-axis (0 = auto-scale)")
    parser.add_argument("--score-strict", type=float, default=0.7,
                        help="Strict on-target score threshold (default 0.7)")
    parser.add_argument("--score-moderate", type=float, default=0.5,
                        help="Moderate on-target score threshold (default 0.5)")
    parser.add_argument("--offtarget-strict", type=int, default=0,
                        help="Strict off-target cap. 0 = auto (default 0)")
    parser.add_argument("--offtarget-moderate", type=int, default=10,
                        help="Moderate off-target cap (default 10)")
    args = parser.parse_args()

    if not HAS_MPL:
        print("ERROR: matplotlib not installed. Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.summary_csv):
        print(f"ERROR: Summary CSV not found: {args.summary_csv}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    guides = read_summary(args.summary_csv)

    if not guides:
        print("WARNING: No guide data in summary CSV.", file=sys.stderr)
        sys.exit(0)

    print(f"Generating plots for {len(guides)} guides...")

    plot_guide_comparison(guides, args.output_dir, args.dpi, args.format, args.y_max_cap,
                          score_strict=args.score_strict,
                          score_moderate=args.score_moderate,
                          offtarget_strict=args.offtarget_strict,
                          offtarget_moderate=args.offtarget_moderate)

    print(f"Done. Plots saved to {args.output_dir}")


if __name__ == "__main__":
    main()
