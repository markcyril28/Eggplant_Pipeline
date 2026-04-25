#!/usr/bin/env python3
"""
Module: 09_comparison_scatter.py
Stage [9] — Guide comparison scatter (mirrors v1 plot_crispr_results.py).

Aggregates all per-gene stage-8 ranked TSVs into scatter plots of on-target
efficacy (x) vs off-target burden (y), colour-coded per gene and marker-coded
per composite KO-score tier. Supports two modes:

    • Single-variant (default) — one figure set using --x-score-column.
    • Multi-variant — pass --score-variants "col:high:mod:max:sfx|…" to emit
      one figure set per scoring column. This is the mode the orchestrator
      uses; spec comes from [crispr_v2.comparison_scatter].score_variants.

Layout (per figure set):
    • Ideal quadrant (high score, low off-targets) shaded green.
    • Threshold lines: configurable per variant for the X axis; y-axis limit
      at the 25th-percentile guide count (floor 5).
    • Tier markers on composite_ko_score: ▲ High, ● Moderate, ■ Low.
    • Directional arrows: "better →" on X, "↑ worse" on Y.
    • Top-N guides per gene by composite KO score get a text label.
    • Optional y-axis clamp via --y-max-cap (over-cap guides render at the
      top edge, label suffixed with "(>cap)").
    • Optional companion table of the top-N ranked guide per gene.

Inputs:
    A directory (conventionally
    `III_RESULT/{GROUP}/09_CRISPR_v2/{genome}/{arm}/08_Ranking_Composite/per_gene/`)
    containing `*.ranked.tsv` (current stage-8 naming) or
    `*_ranked_guides.tsv` (legacy). Aggregate `ranked_guides.tsv` at the
    stage root is NOT consumed — pass the per-gene subdirectory.

Outputs:
    Multi-variant mode groups each variant's 3 files into its own subdir
    (named after the variant suffix), so the output tree reads in pipeline
    order:

        09_Guide_Scatter/
          ├── 1b_moreno_mateos/
          │     ├── guide_comparison.{fmt}
          │     ├── guide_comparison_scatter.{fmt}
          │     └── guide_comparison_table.{fmt}
          ├── 1d_oof_bae/...
          ├── 2a_crispron/...
          └── 8b_composite_ko/...

    Single-variant mode (no suffix) writes the 3 files directly into
    --outdir, preserving the original flat layout for standalone callers.

Usage:
    python3 09_comparison_scatter.py \\
        --input-dir <08_Ranking_Composite/per_gene/> \\
        --outdir    <09_Guide_Scatter/> \\
        --top-n     3 \\
        --y-max-cap 2000 \\
        --score-variants "Moreno-Mateos-Score:40:30:100:moreno_mateos|…" \\
        --dpi       600 \\
        --format    png
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ── Colour palette (matches v1 and the pipeline motif convention) ───────────
GENE_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
]

# Default tier cut-offs on composite_ko_score (0–1 by construction in
# 08_rank_guides) — controls the point MARKER shape (▲ / ● / ■).
TIER_HIGH_DEFAULT = 0.7
TIER_MOD_DEFAULT  = 0.5

# Default X-axis scoring — tuned to Moreno-Mateos-Score (0–100 scale) because
# that is the active filter column in i_crispr_v2CONFIG.toml. Override via
# --x-score-column / --x-axis-high / --x-axis-moderate / --x-axis-max.
# The Moreno-Mateos thresholds (40 strong, 30 moderate) match the cut-offs
# documented in the TOML comment table (see `score_columns` near the top of
# i_crispr_v2CONFIG.toml).
X_COLUMN_DEFAULT       = "Moreno-Mateos-Score"
X_AXIS_HIGH_DEFAULT    = 40.0
X_AXIS_MOD_DEFAULT     = 30.0
X_AXIS_MAX_DEFAULT     = 100.0


def safe_float(v, default=0.0) -> float:
    if v in (None, "", "NA", "N/A"):
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def parse_args():
    p = argparse.ArgumentParser(description="Guide comparison scatter (stage 9)")
    p.add_argument("--input-dir", required=True,
                   help="Directory containing stage-8 ranked TSVs "
                        "(*.ranked.tsv or *_ranked_guides.tsv)")
    p.add_argument("--outdir",    required=True)
    p.add_argument("--top-n",     type=int, default=3,
                   help="Top-N guides per gene to label + include in table")
    p.add_argument("--y-max-cap", type=float, default=0.0,
                   help="Upper bound for the off-target y-axis. 0 = auto "
                        "(max observed × 1.15). Guides above the cap are "
                        "plotted at the top edge with a warning.")
    p.add_argument("--x-score-column", default=X_COLUMN_DEFAULT,
                   help=f"TSV column used for the X-axis. Default: "
                        f"{X_COLUMN_DEFAULT!r} (active stage-1 filter). "
                        "Use 'c_ontarget' for the normalised 0-1 form.")
    p.add_argument("--x-axis-high",     type=float, default=X_AXIS_HIGH_DEFAULT,
                   help="High-efficacy threshold (green dashed line).")
    p.add_argument("--x-axis-moderate", type=float, default=X_AXIS_MOD_DEFAULT,
                   help="Moderate-efficacy threshold (orange dotted line).")
    p.add_argument("--x-axis-max",      type=float, default=X_AXIS_MAX_DEFAULT,
                   help="X-axis upper bound. 0 = auto (max observed × 1.05).")
    p.add_argument("--tier-high",       type=float, default=TIER_HIGH_DEFAULT,
                   help="Composite KO score ≥ this uses the High (▲) marker.")
    p.add_argument("--tier-moderate",   type=float, default=TIER_MOD_DEFAULT,
                   help="Composite KO score ≥ this uses the Moderate (●) marker.")
    p.add_argument("--dpi",       type=int, default=600)
    p.add_argument("--format",    default="png", choices=["png", "svg", "pdf"])
    p.add_argument("--overwrite", dest="overwrite", action="store_true",  default=True)
    p.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    # Multi-variant: emit one scatter set per column, each with its own
    # scale-appropriate thresholds. Pipe-separated spec list. Per-variant
    # fields separated by ':'. Fields (in order, trailing fields optional):
    #   column : high : moderate : axis_max : file_suffix
    # When --score-variants is non-empty, it OVERRIDES the single-column
    # --x-score-column / --x-axis-* flags and produces one figure set per
    # variant (files suffixed with variant[:4] to avoid name collisions).
    # Example:
    #   "Moreno-Mateos-Score:40:30:100:moreno_mateos|
    #    Doench '16-Score:60:40:100:doench16|
    #    composite_ko_score:0.7:0.5:1.0:composite"
    p.add_argument("--score-variants", default="",
                   help="Pipe-separated list of scorer plot specs; each "
                        "colon-separated as column:high:moderate:max:suffix. "
                        "If given, emits one figure set per variant and the "
                        "single --x-score-column flag is ignored.")
    return p.parse_args()


def _slug(name: str) -> str:
    """Filesystem-safe slug for a column name (keeps alphanumerics, '_', '-')."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_").lower()


def parse_variants(spec: str) -> list[dict]:
    """Parse the --score-variants spec string into a list of variant dicts.

    Accepts `column[:high[:moderate[:axis_max[:suffix]]]]` per variant,
    pipe-separated. Missing trailing fields fall back to Moreno-Mateos-style
    defaults (40 / 30 / 100) or a column-name-derived suffix via `_slug()`.
    """
    out: list[dict] = []
    for chunk in spec.split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = [p.strip() for p in chunk.split(":")]
        col = parts[0]
        if not col:
            continue
        high = float(parts[1]) if len(parts) > 1 and parts[1] else X_AXIS_HIGH_DEFAULT
        mod  = float(parts[2]) if len(parts) > 2 and parts[2] else X_AXIS_MOD_DEFAULT
        xmax = float(parts[3]) if len(parts) > 3 and parts[3] else X_AXIS_MAX_DEFAULT
        sfx  = parts[4] if len(parts) > 4 and parts[4] else _slug(col)
        out.append({"column": col, "high": high, "mod": mod,
                    "xmax": xmax, "suffix": sfx})
    return out


def get_gene_stem(path: Path) -> str:
    """Derive per-gene stem from the ranked TSV filename.

    Stage-8 writes `{stem}.ranked.tsv` (e.g. SMEL5_01g008730.ranked.tsv).
    Legacy `{stem}_ranked_guides.tsv` is also accepted for older outputs.
    """
    name = path.name
    for suf in (".ranked.tsv", "_ranked_guides.tsv"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return path.stem


def load_guides(input_dir: Path, x_score_column: str) -> list[dict]:
    """Read all ranked TSVs; each row carries the configurable X-score under
    the neutral key `x_score` plus the fallback component `c_ontarget`."""
    rows: list[dict] = []
    # Current stage-8 naming uses dot separators (*.ranked.tsv); legacy
    # underscore form (*_ranked_guides.tsv) is kept for backwards compat.
    patterns = ["*.ranked.tsv", "*_ranked_guides.tsv"]
    seen: set[Path] = set()
    tsvs: list[Path] = []
    for pat in patterns:
        for p in input_dir.glob(pat):
            if p not in seen:
                seen.add(p)
                tsvs.append(p)
    missing_col_reported = False
    for tsv in sorted(tsvs):
        gene = get_gene_stem(tsv)
        with open(tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if (reader.fieldnames is not None
                    and x_score_column not in reader.fieldnames
                    and not missing_col_reported):
                _log(f"[09_scatter] X-axis column {x_score_column!r} missing "
                     f"from {tsv.name}; falling back to c_ontarget × 100 "
                     "for affected rows.", level="WARN")
                missing_col_reported = True
            for r in reader:
                raw_x = r.get(x_score_column)
                # NotEnoughFlankSeq / NA / empty → fall back to c_ontarget
                # scaled to the same 0-100 range as most CRISPOR columns.
                if raw_x in (None, "", "NA", "N/A", "NotEnoughFlankSeq"):
                    x_val = safe_float(r.get("c_ontarget")) * 100.0
                else:
                    try:
                        x_val = float(raw_x)
                    except (ValueError, TypeError):
                        x_val = safe_float(r.get("c_ontarget")) * 100.0
                rows.append({
                    "gene":       gene,
                    "guide_id":   r.get("guideId") or r.get("guide_id") or r.get("name") or "",
                    "x_score":    x_val,
                    "c_ontarget": safe_float(r.get("c_ontarget")),
                    "offtargets": int(safe_float(r.get("offtargetCount"), 0)),
                    "cfd_sum":    safe_float(r.get("cfd_sum")),
                    "composite":  safe_float(r.get("composite_ko_score")),
                    "sequence":   r.get("targetSeq") or r.get("guideSeq") or "",
                })
    return rows


def classify_tier(composite: float, tier_high: float, tier_mod: float) -> str:
    if composite >= tier_high:
        return "High"
    if composite >= tier_mod:
        return "Moderate"
    return "Low"


def render_scatter(rows: list[dict], outdir: Path, top_n: int,
                   dpi: int, fmt: str,
                   x_score_column: str,
                   x_axis_high: float,
                   x_axis_mod: float,
                   x_axis_max: float,
                   tier_high: float,
                   tier_mod: float,
                   y_max_cap: float = 0.0,
                   file_suffix: str = "") -> None:
    if not rows:
        _log("[09_scatter] No guides in stage-8 outputs — nothing to plot.", level="WARN")
        return

    # Dynamic off-target threshold (25th percentile, floored at 5) — same rule
    # as v1 so the line adapts to whatever the CRISPOR off-target counts are.
    sorted_offs = sorted(r["offtargets"] for r in rows)
    off_thresh  = max(5, sorted_offs[len(sorted_offs) // 4])
    raw_max     = max(sorted_offs) if sorted_offs else 0
    ymax_auto   = raw_max * 1.15 if raw_max > 0 else 1.0
    # Apply the configurable cap: 0 means "auto"; otherwise clamp at cap even
    # if observed values exceed it. Over-cap guides are plotted at the cap
    # line (done via matplotlib clipping) and flagged in the log.
    if y_max_cap and y_max_cap > 0:
        ymax = float(y_max_cap)
        over = sum(1 for r in rows if r["offtargets"] > ymax)
        if over:
            _log(f"[09_scatter] {over} guide(s) exceed y_max_cap={ymax:g} "
                 f"(observed max={raw_max}); they will render at the top edge.",
                 level="WARN")
    else:
        ymax = ymax_auto

    # Per-gene colour map
    genes = sorted({r["gene"] for r in rows})
    gene_color = {g: GENE_COLORS[i % len(GENE_COLORS)] for i, g in enumerate(genes)}

    # Top-N guides per gene by composite KO score (for labels + table)
    per_gene: dict[str, list[dict]] = {}
    for r in rows:
        per_gene.setdefault(r["gene"], []).append(r)
    best = []
    for g in genes:
        best.extend(sorted(per_gene[g], key=lambda x: -x["composite"])[:top_n])
    best_ids = {(b["gene"], b["guide_id"]) for b in best}

    # Compute X-axis bounds. Auto-widens to negative lower bound when any
    # observed score is < 0 (e.g. Doench-RuleSet3-Score spans −3 to +3);
    # otherwise anchors at 0. x_axis_max = 0 means auto-scale upper bound.
    x_values = [r["x_score"] for r in rows]
    raw_xmax = max(x_values, default=1.0)
    raw_xmin = min(x_values, default=0.0)
    xlim_lo  = min(0.0, raw_xmin * 1.05) if raw_xmin < 0 else 0.0
    if x_axis_max and x_axis_max > 0:
        xmax = float(x_axis_max)
    else:
        xmax = raw_xmax * 1.05 if raw_xmax > 0 else 1.0
    x_span = xmax - xlim_lo if xmax > xlim_lo else 1.0

    def _draw_scatter(ax):
        ax.set_xlim(xlim_lo, xmax)
        ax.set_ylim(0, ymax)

        # Ideal-quadrant shading (high score, low off-targets). axhspan's
        # xmin/xmax are AXIS FRACTIONS, so normalise against the full x-span
        # (which may start below 0 for RS3-style scorers).
        xmin_ideal = max(0.0, min(1.0, (x_axis_high - xlim_lo) / x_span))
        ax.axhspan(0, off_thresh, xmin=xmin_ideal, xmax=1.0,
                   color="#2ca02c", alpha=0.08, zorder=0)

        # Threshold lines
        ax.axvline(x_axis_high, color="#2ca02c", linestyle="--", lw=1.3, alpha=0.8)
        ax.axvline(x_axis_mod,  color="#ff7f0e", linestyle=":",  lw=1.1, alpha=0.7)
        ax.axhline(off_thresh,  color="#1f77b4", linestyle="--", lw=1.1, alpha=0.7)

        # Threshold annotations (x-offset scales with axis span, anchored at
        # the left edge for the off-target label — works for negative xlim_lo).
        x_tick = x_span * 0.012
        ax.text(x_axis_high + x_tick, ymax * 0.96,
                f"≥ {x_axis_high:g}\n(High)",
                fontsize=7.5, color="#2ca02c", va="top")
        ax.text(x_axis_mod + x_tick, ymax * 0.96,
                f"≥ {x_axis_mod:g}\n(Moderate)",
                fontsize=7.5, color="#ff7f0e", va="top")
        ax.text(xlim_lo + x_tick, off_thresh + ymax * 0.012,
                f"≤ {off_thresh} off-targets",
                fontsize=7.5, color="#1f77b4", va="bottom")

        # Directional arrows
        ax.annotate("", xy=(0.97, -0.055), xytext=(0.75, -0.055),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="#555", lw=1.3),
                    annotation_clip=False)
        ax.text(0.86, -0.09, "better →", fontsize=8, color="#555",
                ha="center", transform=ax.transAxes)
        ax.annotate("", xy=(-0.055, 0.97), xytext=(-0.055, 0.75),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="#555", lw=1.3),
                    annotation_clip=False)
        ax.text(-0.085, 0.86, "↑ worse", fontsize=8, color="#555",
                rotation=90, va="center", transform=ax.transAxes)

        # Points — clamp y to ymax so over-cap guides render at the top edge
        # rather than being cut off silently. Labels use the clamped position
        # but still show the guide id so the reader can look up the real value.
        tier_marker = {"High": "^", "Moderate": "o", "Low": "s"}
        for r in rows:
            tier = classify_tier(r["composite"], tier_high, tier_mod)
            y_plot = min(r["offtargets"], ymax)
            ax.scatter(r["x_score"], y_plot,
                       c=gene_color[r["gene"]], marker=tier_marker[tier],
                       s=70, edgecolors="black", linewidth=0.5, zorder=4)
            if (r["gene"], r["guide_id"]) in best_ids and r["guide_id"]:
                lbl = r["guide_id"]
                if r["offtargets"] > ymax:
                    lbl += f" (>{int(ymax)})"
                ax.annotate(lbl,
                            (r["x_score"], y_plot),
                            fontsize=6, alpha=0.75, xytext=(3, 3),
                            textcoords="offset points")

        ax.set_xlabel(f"On-Target Score ({x_score_column})", fontsize=11)
        ax.set_ylabel("Off-Target Hits (CRISPOR offtargetCount)", fontsize=11)
        ax.set_title(
            f"Guide Selection: {x_score_column} vs Off-Target Count",
            fontsize=13, fontweight="bold")

        gene_patches = [mpatches.Patch(color=c, label=g) for g, c in gene_color.items()]
        tier_handles = [
            plt.scatter([], [], marker="^", c="gray", s=50,
                        label=f"High (KO≥{tier_high:g})",
                        edgecolors="black", linewidth=0.4),
            plt.scatter([], [], marker="o", c="gray", s=50,
                        label=f"Moderate (KO≥{tier_mod:g})",
                        edgecolors="black", linewidth=0.4),
            plt.scatter([], [], marker="s", c="gray", s=50, label="Low",
                        edgecolors="black", linewidth=0.4),
        ]
        thresh_handles = [
            plt.Line2D([0], [0], color="#2ca02c", linestyle="--", lw=1.3,
                       label=f"High score ≥ {x_axis_high:g}"),
            plt.Line2D([0], [0], color="#ff7f0e", linestyle=":",  lw=1.1,
                       label=f"Moderate score ≥ {x_axis_mod:g}"),
            plt.Line2D([0], [0], color="#1f77b4", linestyle="--", lw=1.1,
                       label=f"Off-target limit ({off_thresh})"),
        ]
        ax.legend(handles=gene_patches + tier_handles + thresh_handles,
                  fontsize=7, loc="upper left", ncol=2, framealpha=0.92)

    def _draw_table(ax_t):
        col_labels = ["Gene", "Guide", x_score_column, "Off-targets",
                      "cfd_sum", "KO score", "Sequence (5'→3')"]
        table_rows = [
            [b["gene"], b["guide_id"],
             f"{b['x_score']:.3f}", str(b["offtargets"]),
             f"{b['cfd_sum']:.3f}", f"{b['composite']:.3f}",
             b["sequence"]]
            for b in best
        ]
        ax_t.axis("off")
        if not table_rows:
            return
        tbl = ax_t.table(cellText=table_rows, colLabels=col_labels,
                         loc="center", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7.5)
        tbl.scale(1.0, 1.25)

    # Each variant gets its own subdirectory so the three files (combined,
    # scatter-only, table-only) stay grouped. The subdir name is the variant
    # suffix, which is stage-ordered (e.g. 1b_moreno_mateos) so a plain file
    # listing of outdir reads in pipeline order. Single-variant mode (no
    # suffix) writes directly into outdir for back-compat with standalone use.
    variant_dir = outdir / file_suffix if file_suffix else outdir
    variant_dir.mkdir(parents=True, exist_ok=True)

    # Combined (scatter + table)
    fig_c, ax_c = plt.subplots(figsize=(11, 11))
    fig_c.subplots_adjust(bottom=0.50, top=0.94, left=0.10, right=0.95)
    _draw_scatter(ax_c)
    n_rows = len(best)
    if n_rows:
        ax_t_c = fig_c.add_axes([0.06, 0.10, 0.90,
                                 min(0.06 + n_rows * 0.055, 0.26)])
        _draw_table(ax_t_c)
    path_c = variant_dir / f"guide_comparison.{fmt}"
    fig_c.savefig(path_c, dpi=dpi, bbox_inches="tight")
    plt.close(fig_c)
    _log(f"[09_scatter] Combined plot -> {path_c}")

    # Scatter only
    fig_s, ax_s = plt.subplots(figsize=(11, 7))
    fig_s.subplots_adjust(bottom=0.12, top=0.94, left=0.12, right=0.95)
    _draw_scatter(ax_s)
    path_s = variant_dir / f"guide_comparison_scatter.{fmt}"
    fig_s.savefig(path_s, dpi=dpi, bbox_inches="tight")
    plt.close(fig_s)
    _log(f"[09_scatter] Scatter-only -> {path_s}")

    # Table only
    if n_rows:
        fig_height = max(3.0, 1.2 + n_rows * 0.42)
        fig_t, ax_t = plt.subplots(figsize=(11, fig_height))
        fig_t.subplots_adjust(top=0.82, bottom=0.04, left=0.02, right=0.98)
        _draw_table(ax_t)
        path_t = variant_dir / f"guide_comparison_table.{fmt}"
        fig_t.savefig(path_t, dpi=dpi, bbox_inches="tight")
        plt.close(fig_t)
        _log(f"[09_scatter] Table-only -> {path_t}")


def main():
    args = parse_args()
    if not HAS_MPL:
        _log("[09_scatter] matplotlib not installed — cannot render.", level="ERROR")
        sys.exit(1)

    input_dir = Path(args.input_dir)
    outdir    = Path(args.outdir)

    if not input_dir.is_dir():
        _log(f"[09_scatter] Input directory not found: {input_dir}", level="ERROR")
        sys.exit(1)

    # Build the list of variants to render. If --score-variants is non-empty,
    # one figure set per variant; otherwise fall back to the single-column
    # --x-score-column / --x-axis-* flags (renders one "default" set with no
    # filename suffix, preserving back-compat with earlier callers).
    variants = parse_variants(args.score_variants)
    if not variants:
        variants = [{
            "column": args.x_score_column,
            "high":   args.x_axis_high,
            "mod":    args.x_axis_moderate,
            "xmax":   args.x_axis_max,
            "suffix": "",
        }]
        _log(f"[09_scatter] Single-variant mode: {variants[0]['column']!r} "
             f"(high={variants[0]['high']:g}, mod={variants[0]['mod']:g})")
    else:
        _log(f"[09_scatter] Multi-variant mode: {len(variants)} scorer(s) — "
             + ", ".join(v["column"] for v in variants))

    # Short-circuit overwrite check: only skip when ALL expected combined
    # files already exist. Each variant owns a subdirectory under outdir
    # (outdir/{suffix}/guide_comparison.{fmt}); single-variant mode uses
    # outdir directly. Partial runs (e.g. after adding a new variant) still
    # render the missing ones instead of returning early.
    if not args.overwrite:
        def _combined_path(v: dict) -> Path:
            vdir = outdir / v["suffix"] if v["suffix"] else outdir
            return vdir / f"guide_comparison.{args.format}"
        missing = [v for v in variants if not _combined_path(v).exists()]
        if not missing:
            _log(f"[09_scatter] Skipping (overwrite=false): "
                 f"all {len(variants)} variant output(s) already present in {outdir}")
            return
        variants = missing

    for v in variants:
        rows = load_guides(input_dir, v["column"])
        if not rows:
            continue
        _log(f"[09_scatter] Rendering variant {v['column']!r} "
             f"(high={v['high']:g}, mod={v['mod']:g}, xmax={v['xmax']:g})")
        render_scatter(
            rows, outdir, args.top_n, args.dpi, args.format,
            x_score_column=v["column"],
            x_axis_high=v["high"],
            x_axis_mod=v["mod"],
            x_axis_max=v["xmax"],
            tier_high=args.tier_high,
            tier_mod=args.tier_moderate,
            y_max_cap=args.y_max_cap,
            file_suffix=v["suffix"],
        )


if __name__ == "__main__":
    main()
