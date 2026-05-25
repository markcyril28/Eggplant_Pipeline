#!/usr/bin/env python3
"""
Compare AlphaFold 3 confidence metrics across experimental conditions for the
trimeric HAP2 + SmelDMPv5_10.610 complex.

Reads `*_summary_confidences_{0..4}.json` files under each condition subfolder
of AlphaFold3_Experiment/ and produces a single multi-panel figure:

  (a) Top-model confidence (ipTM, pTM, ranking_score) per condition.
  (b) Per-rank ranking_score (5 models per condition) — variance across seeds.
  (c) HAP2 (A/B/C) x DMP (D) interface ipTM for the top model per condition.
  (d) 4x4 chain-pair ipTM heatmaps for the top model of each condition.

Chain assignment:  A/B/C = HAP2 trimer,  D = DMP monomer.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

DEFAULT_EXP_DIR = Path(
    "III_RESULT/DMP-HAP2/08_Protein_Structure/"
    "GPE001970_SMEL5/AlphaFold3_Experiment"
)

CONDITION_LABELS = {
    "no_template":    "No template",
    "seed_20":        "Seed 20",
    "seed_100":       "Seed 100",
    "swiss_template": "SWISS template",
}

CONDITION_ORDER = ["no_template", "seed_20", "seed_100", "swiss_template"]

# Heatmap tick labels (compact) and the biological identity legend.
CHAIN_TICK_LABELS = ["A", "B", "C", "D"]
CHAIN_IDENTITY = [
    ("A", "HAP2 chain 1"),
    ("B", "HAP2 chain 2"),
    ("C", "HAP2 chain 3"),
    ("D", "SmelDMPv5_10.610 (DMP)"),
]

CONDITION_COLORS = {
    "no_template":    "#5C87BF",  # steel blue
    "seed_20":        "#66CCF2",  # cyan
    "seed_100":       "#8059C2",  # periwinkle
    "swiss_template": "#D4A621",  # golden amber
}

# Global font size scaling for readability.
plt.rcParams.update({
    "font.size":         11,
    "axes.titlesize":    12,
    "axes.labelsize":    11,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "legend.fontsize":   9,
})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exp-dir", default=None,
                   help="AlphaFold3_Experiment directory (auto-resolved if omitted)")
    p.add_argument("--output-dir", default="Figures",
                   help="Output subdirectory (default: Figures)")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def find_condition_dir(exp_dir: Path, key: str) -> Path | None:
    """Return the folder whose name ends with `_{key}` (no trailing .zip)."""
    suffix = f"_{key}"
    for sub in exp_dir.iterdir():
        if sub.is_dir() and sub.name.endswith(suffix):
            return sub
    return None


def load_summaries(cond_dir: Path) -> list[dict]:
    """Load all summary_confidences_*.json in rank order (0..4)."""
    files = sorted(cond_dir.glob("*_summary_confidences_*.json"),
                   key=lambda p: int(re.search(r"_(\d+)\.json$", p.name).group(1)))
    return [json.loads(f.read_text()) for f in files]


def panel_a_top_metrics(ax, data):
    """Grouped bar chart: ipTM / pTM / ranking_score for top-ranked model."""
    metrics = ["ipTM\n(interface)", "pTM\n(overall)", "Ranking score\n(AF3 composite)"]
    keys = ["iptm", "ptm", "ranking_score"]
    n_cond = len(data)
    x = np.arange(len(metrics))
    width = 0.8 / n_cond

    for i, (cond, summaries) in enumerate(data.items()):
        if not summaries:
            continue
        top = summaries[0]
        vals = [top.get(k, np.nan) for k in keys]
        offset = (i - (n_cond - 1) / 2) * width
        ax.bar(x + offset, vals, width,
               label=CONDITION_LABELS[cond],
               color=CONDITION_COLORS[cond],
               edgecolor="black", linewidth=0.6)
        for xi, v in zip(x + offset, vals):
            ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score (0 to 1)")
    ax.set_title("(a) Top-model confidence by condition",
                 loc="left", fontweight="bold")
    ax.axhline(0.5, linestyle="--", color="gray", linewidth=0.8, alpha=0.6)
    ax.text(ax.get_xlim()[1] * 0.99, 0.51, "0.50 reference",
            ha="right", va="bottom", fontsize=7.5, color="gray", style="italic")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.9, ncol=2)


def panel_b_per_rank(ax, data):
    """Line plot: ranking_score across the 5 model ranks for each condition."""
    for cond, summaries in data.items():
        if not summaries:
            continue
        scores = [s.get("ranking_score", np.nan) for s in summaries]
        ranks = np.arange(len(scores))
        ax.plot(ranks, scores, marker="o", linewidth=2.0, markersize=7,
                label=CONDITION_LABELS[cond],
                color=CONDITION_COLORS[cond],
                markeredgecolor="black", markeredgewidth=0.5)
    ax.set_xlabel("Model rank (0 = top model)")
    ax.set_ylabel("Ranking score")
    ax.set_xticks(range(5))
    ax.set_xticklabels(["0\n(top)", "1", "2", "3", "4"])
    ax.set_title("(b) Ranking score across the 5 sampled models",
                 loc="left", fontweight="bold")
    ax.grid(linestyle=":", alpha=0.4)
    ax.legend(loc="lower left", framealpha=0.9, ncol=2)


def panel_c_interface_iptm(ax, data):
    """Per-chain HAP2 (A/B/C) x DMP (D) interface ipTM for the top model."""
    hap2_chains = ["A", "B", "C"]
    x = np.arange(len(hap2_chains))
    n_cond = len(data)
    width = 0.8 / n_cond

    max_v = 0.0
    for i, (cond, summaries) in enumerate(data.items()):
        if not summaries:
            continue
        top = summaries[0]
        pair = np.array(top["chain_pair_iptm"])
        # Row D (index 3), columns A/B/C (0/1/2). Use mean of (D,X) and (X,D).
        d_vs = [(pair[3, j] + pair[j, 3]) / 2.0 for j in range(3)]
        max_v = max(max_v, max(d_vs))
        offset = (i - (n_cond - 1) / 2) * width
        ax.bar(x + offset, d_vs, width,
               label=CONDITION_LABELS[cond],
               color=CONDITION_COLORS[cond],
               edgecolor="black", linewidth=0.6)
        for xi, v in zip(x + offset, d_vs):
            ax.text(xi, v + 0.006, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=8.5, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"DMP (D) x HAP2-{c}" for c in hap2_chains])
    ax.set_ylabel("Pair ipTM (interface)")
    ax.set_ylim(0, max(0.40, max_v * 1.25))
    ax.set_title("(c) HAP2-DMP interface ipTM by HAP2 chain (top model)",
                 loc="left", fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.9, ncol=2)


def panel_d_heatmaps(axes, data):
    """4 small chain_pair_iptm heatmaps for each condition's top model.

    Adds:
      - biological identity in tick labels (A: HAP2, ..., D: DMP)
      - a red box outlining the HAP2-DMP interface row/column (chain D)
      - a thin separator showing the HAP2-trimer block vs DMP
    """
    vmax = 1.0
    im = None
    for ax, (cond, summaries) in zip(axes, data.items()):
        if not summaries:
            ax.set_visible(False)
            continue
        top = summaries[0]
        pair = np.array(top["chain_pair_iptm"])
        im = ax.imshow(pair, vmin=0, vmax=vmax, cmap="viridis", aspect="equal")
        ax.set_xticks(range(4))
        ax.set_yticks(range(4))
        ax.set_xticklabels(
            ["A\nHAP2", "B\nHAP2", "C\nHAP2", "D\nDMP"], fontsize=8.5)
        ax.set_yticklabels(
            ["A: HAP2", "B: HAP2", "C: HAP2", "D: DMP"], fontsize=8.5)
        ax.set_title(CONDITION_LABELS[cond], fontsize=10.5, fontweight="bold")
        for i in range(4):
            for j in range(4):
                v = pair[i, j]
                color = "white" if v < 0.55 else "black"
                weight = "bold" if (i == 3 or j == 3) and i != j else "normal"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color=color, fontsize=8, fontweight=weight)
        # Separator line between HAP2 trimer block (A/B/C) and DMP (D).
        ax.axhline(2.5, color="white", linewidth=1.2, alpha=0.85)
        ax.axvline(2.5, color="white", linewidth=1.2, alpha=0.85)
        # Highlight the HAP2-DMP interface row (D) and column (D) with a red border.
        ax.add_patch(Rectangle((-0.5, 2.5), 4, 1, fill=False,
                               edgecolor="#D32F2F", linewidth=1.4))
        ax.add_patch(Rectangle((2.5, -0.5), 1, 4, fill=False,
                               edgecolor="#D32F2F", linewidth=1.4))
    return im


def main():
    args = parse_args()

    if args.exp_dir:
        exp_dir = Path(args.exp_dir).resolve()
    else:
        script_dir = Path(__file__).resolve().parent
        pipeline_dir = script_dir.parent.parent
        exp_dir = pipeline_dir / DEFAULT_EXP_DIR

    if not exp_dir.is_dir():
        print(f"ERROR: experiment dir not found: {exp_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir = exp_dir / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Experiment dir: {exp_dir}")
    print(f"Output dir:     {out_dir}\n")

    data: dict[str, list[dict]] = {}
    for key in CONDITION_ORDER:
        cond_dir = find_condition_dir(exp_dir, key)
        if cond_dir is None:
            print(f"  [warn] no folder found for condition '{key}'", file=sys.stderr)
            continue
        summaries = load_summaries(cond_dir)
        if not summaries:
            print(f"  [warn] no summary_confidences_*.json under {cond_dir.name}",
                  file=sys.stderr)
            continue
        # AF3 server output is already rank-ordered (model_0 = top); confirm by sort.
        summaries.sort(key=lambda s: -s.get("ranking_score", -np.inf))
        data[key] = summaries
        top = summaries[0]
        print(f"  {key:15s}  models={len(summaries)}  top ipTM={top['iptm']:.2f}  "
              f"pTM={top['ptm']:.2f}  rank={top['ranking_score']:.2f}")

    if not data:
        print("\nNo data loaded.", file=sys.stderr)
        sys.exit(1)

    # ── Build figure ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 13), constrained_layout=False)
    gs = GridSpec(3, 4, figure=fig,
                  height_ratios=[1.0, 1.0, 1.15],
                  hspace=0.95, wspace=0.40,
                  left=0.07, right=0.94, top=0.90, bottom=0.10)

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_c = fig.add_subplot(gs[1, 0:4])

    panel_a_top_metrics(ax_a, data)
    panel_b_per_rank(ax_b, data)
    panel_c_interface_iptm(ax_c, data)

    # Panel d: 4 heatmaps along the bottom row.
    ax_d_axes = [fig.add_subplot(gs[2, i]) for i in range(4)]
    im = panel_d_heatmaps(ax_d_axes, data)

    # Determine the actual top edge of the heatmap row so the panel d header
    # sits cleanly above it (not on top of the subplot titles).
    fig.canvas.draw()
    d_top = max(ax.get_position().y1 for ax in ax_d_axes)
    header_y    = min(0.40, d_top + 0.055)
    subhead_y   = header_y - 0.020

    # Panel d title + explicit chain key (above the heatmaps).
    fig.text(0.07, header_y,
             "(d) Chain-pair ipTM matrices  -  top model per condition",
             fontsize=12.5, fontweight="bold")
    fig.text(0.07, subhead_y,
             "Chains  A, B, C  =  HAP2 trimer    |    Chain  D  =  SmelDMPv5_10.610 (DMP)"
             "    -    red box = HAP2-DMP interface",
             fontsize=10, style="italic", color="#444")

    # Shared colorbar for panel d heatmaps, aligned to the heatmap row.
    d_bottom = min(ax.get_position().y0 for ax in ax_d_axes)
    cbar_ax = fig.add_axes([0.955, d_bottom, 0.013, d_top - d_bottom])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Pair ipTM (0 = no confidence, 1 = high)", fontsize=10)

    fig.suptitle(
        "AlphaFold 3 condition comparison:\n"
        "trimeric HAP2 (chains A+B+C) + SmelDMPv5_10.610 (chain D)",
        fontsize=15, fontweight="bold", y=0.97,
    )

    # Footnote: metric definitions.
    fig.text(0.07, 0.022,
             "ipTM = interface predicted TM-score (cross-chain confidence)   |   "
             "pTM = predicted TM-score (overall fold)   |   "
             "Ranking score = AF3 composite used to rank the 5 sampled models.",
             fontsize=9, color="#333")

    # Save
    stem = "AF3_experiment_condition_comparison"
    out_png = out_dir / f"{stem}.png"
    out_jpg = out_dir / f"{stem}.jpg"
    out_svg = out_dir / f"{stem}.svg"
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out_jpg, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out_svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"\nSaved: {out_png.name}")
    print(f"Saved: {out_jpg.name}")
    print(f"Saved: {out_svg.name}")


if __name__ == "__main__":
    main()
