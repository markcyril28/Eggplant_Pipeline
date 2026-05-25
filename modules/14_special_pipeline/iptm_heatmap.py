#!/usr/bin/env python3
"""
AF3 confidence heatmaps for Stage 14 (Interaction Domain Mapping).

Reads each pair slot under <complex-dir>/<pair>/ for the rank-1 AlphaFold3
model's summary_confidences_0.json (file naming preserved from the AF3
web download zip), extracts the overall iptm / ptm / ranking_score, and
renders one HAP2-variant x DMP-variant heatmap per metric. The caller is
expected to pass an already-per-stoichiometry complex root (e.g.
$EXP_OUT_DIR/<stoich>/02_Complexes); --stoich is now name-only and used
only to suffix the output file names.

Metrics rendered:
    ipTM            interface predicted TM-score (interface accuracy; primary)
    pTM             overall predicted TM-score (sanity check on fold quality)
    ranking_score   AF3 ranking score = 0.8 * ipTM + 0.2 * pTM (AF3 default rank)

Outputs (per stoichiometry label):
    <output-dir>/iptm_heatmap_<stoich>.{pdf|png|svg}
    <output-dir>/ptm_heatmap_<stoich>.{pdf|png|svg}
    <output-dir>/ranking_score_heatmap_<stoich>.{pdf|png|svg}
    <output-dir>/iptm_table_<stoich>.tsv   (all three metrics + source JSON)

Independent of interface_analysis / binding_energy: only needs the AF3
JSON files. Run any time after prepare_complexes has produced summary
JSONs in the per-pair slots.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# CRISPR-Local scores for the guide-anchored DMP truncation ladder.
# Source: 14_Interaction_Domain_MappingCONFIG.toml lines 233-239 (dmp_variants).
# Embedded as a second tick-label line so the heatmap is self-contained.
CRISPR_LOCAL_SCORES: dict[str, float] = {
    "delGuide4":  0.1172,
    "delGuide16": 0.4235,
    "delGuide17": 0.8725,
    "delGuide20": 0.4861,
    "delGuide37": 0.5374,
    "delGuide46": 0.4683,
    "delGuide50": 0.4767,
}


def format_dmp_tick(name: str) -> str:
    score = CRISPR_LOCAL_SCORES.get(name)
    if score is None:
        return name
    return f"{name}\n(CL {score:.2f})"


def find_summary_json(slot_dir: Path, model_idx: int) -> Path | None:
    cands = sorted(slot_dir.glob(f"*_summary_confidences_{model_idx}.json"))
    return cands[0] if cands else None


def parse_pair_label(label: str) -> tuple[str, str]:
    if "__" not in label:
        return label, "WT"
    h, d = label.split("__", 1)
    return h, d


def collect_metrics(
    complex_root: Path,
    stoich: str,
    pair_labels: list[str],
    model_idx: int,
) -> pd.DataFrame:
    rows = []
    for pl in pair_labels:
        slot = complex_root / pl
        h, d = parse_pair_label(pl)
        js = find_summary_json(slot, model_idx)
        if js is None:
            print(f"[WARN] no summary_confidences_{model_idx}.json in {slot}", file=sys.stderr)
            rows.append((pl, h, d, np.nan, np.nan, np.nan, ""))
            continue
        try:
            data = json.loads(js.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"[WARN] could not parse {js}: {exc}", file=sys.stderr)
            rows.append((pl, h, d, np.nan, np.nan, np.nan, str(js.name)))
            continue
        iptm = float(data.get("iptm", np.nan))
        ptm = float(data.get("ptm", np.nan))
        rs_raw = data.get("ranking_score")
        if rs_raw is None and not (np.isnan(iptm) or np.isnan(ptm)):
            rs = 0.8 * iptm + 0.2 * ptm
        else:
            rs = float(rs_raw) if rs_raw is not None else np.nan
        rows.append((pl, h, d, iptm, ptm, rs, js.name))
    return pd.DataFrame(
        rows,
        columns=["pair", "hap2", "dmp", "iptm", "ptm", "ranking_score", "source_json"],
    )


def render_heatmap(
    df: pd.DataFrame,
    metric_col: str,
    metric_label: str,
    title: str,
    hap2_vars: list[str],
    dmp_vars: list[str],
    out_path: Path,
    cmap_name: str,
    vmin: float,
    vmax: float,
    strong: float | None = None,
    borderline: float | None = None,
    alone_controls: list[tuple[str, dict]] | None = None,
) -> None:
    """Render a single HAP2 x DMP heatmap for one metric column.

    strong / borderline (optional): if both provided, cells >= strong get **,
    cells >= borderline get *, and threshold lines are drawn on the colorbar.
    Pass None to omit threshold annotations (recommended for pTM / ranking_score
    where the ipTM 0.60/0.80 thresholds do not transfer cleanly).
    """
    mat = np.full((len(hap2_vars), len(dmp_vars)), np.nan)
    for _, row in df.iterrows():
        if row.hap2 in hap2_vars and row.dmp in dmp_vars:
            mat[hap2_vars.index(row.hap2), dmp_vars.index(row.dmp)] = row[metric_col]

    # Fixed cell size so main-heatmap squares and control-strip squares
    # render identical on disk; height/width derived from cell counts.
    cell_in = 0.55
    margin_w = 2.5   # y-labels + colorbar + tight side padding
    margin_h = 1.6   # title + x-labels (kept small; bbox_inches="tight" crops the rest)
    fig_w = len(dmp_vars) * cell_in + margin_w
    # Vertical room reserved for the manually-placed controls strip below:
    # 2.6 in = control-strip title + 1 cell + x-labels + comfortable bottom pad.
    ctrl_block_in = 2.6 if alone_controls else 0.0
    fig_h = len(hap2_vars) * cell_in + margin_h + ctrl_block_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax_ctrl = None  # placed after main heatmap layout is fixed (see below)

    cmap = mpl.colormaps[cmap_name].copy()
    cmap.set_bad("#dddddd")
    im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap=cmap, aspect="equal")

    ax.set_xticks(np.arange(len(dmp_vars)))
    ax.set_yticks(np.arange(len(hap2_vars)))
    ax.set_xticklabels(
        [format_dmp_tick(v) for v in dmp_vars],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=9,
    )
    ax.set_yticklabels(hap2_vars, fontsize=9)
    ax.set_xlabel("DMP variant (CL = CRISPR-Local on-target score)", fontsize=10.5)
    ax.set_ylabel("HAP2 variant", fontsize=10.5)
    ax.set_title(title, fontsize=11)

    annotate_thresholds = strong is not None and borderline is not None
    for i in range(len(hap2_vars)):
        for j in range(len(dmp_vars)):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", color="black", fontsize=8)
                continue
            marker = ""
            if annotate_thresholds:
                if v >= strong:
                    marker = "**"
                elif v >= borderline:
                    marker = "*"
            ax.text(
                j,
                i,
                f"{v:.2f}{marker}",
                ha="center",
                va="center",
                color="black",
                fontsize=9,
            )

    # Constrain the colorbar to match the heatmap height. shrink keeps it
    # short relative to the axes; fraction/pad keep it slim and close.
    cbar = fig.colorbar(
        im, ax=ax, label=metric_label,
        fraction=0.025, pad=0.02, shrink=0.65, aspect=18,
    )
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(metric_label, fontsize=9)
    if annotate_thresholds:
        cbar.ax.axhline(strong, color="black", lw=0.6)
        cbar.ax.axhline(borderline, color="black", lw=0.6, linestyle="--")

    # Main-heatmap layout first; placement of the controls strip below
    # depends on the main heatmap's final bbox.
    fig.tight_layout()

    if alone_controls:
        ctrl_names = [c[0] for c in alone_controls]
        ctrl_vals = []
        for _, metrics in alone_controls:
            v = metrics.get(metric_col)
            ctrl_vals.append(np.nan if v is None else float(v))

        # Cell dimensions in figure-fraction coords (inches / fig dim).
        cell_w_fig = cell_in / fig_w
        cell_h_fig = cell_in / fig_h
        ctrl_width_fig = cell_w_fig * len(ctrl_names)
        ctrl_height_fig = cell_h_fig

        # Main heatmap data area in figure coords (use imshow's image bbox so
        # we align under the actual cells, not the axes box).
        fig.canvas.draw()  # ensure transforms are finalised
        renderer = fig.canvas.get_renderer()
        img_bbox = im.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        main_left = img_bbox.x0
        main_width = img_bbox.width
        main_bottom = img_bbox.y0

        # Centre the controls strip horizontally under the heatmap data area.
        ctrl_left = main_left + (main_width - ctrl_width_fig) / 2
        # Vertical gap below main heatmap: room for x-tick labels + the
        # "DMP variant" xlabel + a small breathing space.
        gap_fig = 1.4 / fig_h
        ctrl_bottom = main_bottom - gap_fig - ctrl_height_fig

        ax_ctrl = fig.add_axes([ctrl_left, ctrl_bottom, ctrl_width_fig, ctrl_height_fig])
        ctrl_mat = np.array(ctrl_vals).reshape(1, -1)
        ax_ctrl.imshow(
            ctrl_mat, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto",
            extent=(-0.5, len(ctrl_names) - 0.5, -0.5, 0.5),
        )
        ax_ctrl.set_xticks(np.arange(len(ctrl_names)))
        ax_ctrl.set_xticklabels(
            ctrl_names, rotation=30, ha="right", rotation_mode="anchor", fontsize=9,
        )
        ax_ctrl.set_yticks([0])
        ax_ctrl.set_yticklabels(["Alone (no partner)"], fontsize=9)
        ax_ctrl.set_title(
            f"Reference controls -- {metric_label} for each protein folded alone "
            "(no interacting partner)",
            fontsize=10.5, pad=4,
        )
        for j, v in enumerate(ctrl_vals):
            txt = "n/a" if np.isnan(v) else f"{v:.2f}"
            ax_ctrl.text(j, 0, txt, ha="center", va="center", color="black", fontsize=9)

    # bbox_inches="tight" auto-trims top whitespace and captures the
    # controls strip's x-tick labels even if they overflow the figure.
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--complex-dir", required=True, help="02_Complexes/ root for this experiment")
    ap.add_argument("--stoich", required=True, help="stoichiometry label (e.g. monomeric)")
    ap.add_argument("--hap2-variants", required=True, help="space-separated HAP2 variant names")
    ap.add_argument("--dmp-variants", required=True, help="space-separated DMP variant names")
    ap.add_argument("--pair-labels", required=True, help="space-separated pair labels")
    ap.add_argument("--output-dir", required=True, help="07_Summary/ for this experiment")
    ap.add_argument("--format", default="pdf", choices=["pdf", "png", "svg", "jpeg"])
    ap.add_argument("--cmap", default="RdYlGn")
    ap.add_argument("--vmin", type=float, default=0.0)
    ap.add_argument("--vmax", type=float, default=1.0)
    ap.add_argument("--iptm-strong", type=float, default=0.80)
    ap.add_argument("--iptm-borderline", type=float, default=0.60)
    ap.add_argument("--model-index", type=int, default=0, help="AF3 rank index to read")
    ap.add_argument(
        "--alone-control",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Add a 'protein alone (no partner)' reference cell to the controls strip. "
             "PATH points to an AF3 summary_confidences_*.json file. Repeat per control.",
    )
    args = ap.parse_args()

    complex_root = Path(args.complex_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hap2_vars = args.hap2_variants.split()
    dmp_vars = args.dmp_variants.split()
    pair_labels = args.pair_labels.split()

    if not pair_labels:
        print("[ERROR] empty --pair-labels", file=sys.stderr)
        return 2

    # Parse alone-control JSONs into [(name, {iptm, ptm, ranking_score}), ...]
    alone_controls: list[tuple[str, dict]] = []
    for spec in args.alone_control:
        if "=" not in spec:
            print(f"[WARN] --alone-control '{spec}' missing '=NAME'; skipping", file=sys.stderr)
            continue
        name, json_path = spec.split("=", 1)
        jp = Path(json_path)
        if not jp.exists():
            print(f"[WARN] alone-control JSON not found: {jp}", file=sys.stderr)
            alone_controls.append((name, {}))
            continue
        try:
            data = json.loads(jp.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"[WARN] could not parse {jp}: {exc}", file=sys.stderr)
            alone_controls.append((name, {}))
            continue
        iptm = data.get("iptm")
        ptm = data.get("ptm")
        rs = data.get("ranking_score")
        if rs is None and iptm is not None and ptm is not None:
            rs = 0.8 * float(iptm) + 0.2 * float(ptm)
        alone_controls.append((name, {"iptm": iptm, "ptm": ptm, "ranking_score": rs}))
        print(f"[OK] alone-control '{name}': iptm={iptm} ptm={ptm} ranking_score={rs}")

    df = collect_metrics(complex_root, args.stoich, pair_labels, args.model_index)
    tsv_path = out_dir / f"iptm_table_{args.stoich}.tsv"
    df.to_csv(tsv_path, sep="\t", index=False, float_format="%.4f")
    print(
        f"[OK] wrote {tsv_path}  "
        f"(ipTM {df['iptm'].notna().sum()}/{len(df)}, "
        f"pTM {df['ptm'].notna().sum()}/{len(df)}, "
        f"ranking_score {df['ranking_score'].notna().sum()}/{len(df)})"
    )

    panels = [
        {
            "metric_col": "iptm",
            "metric_label": "ipTM",
            "title": f"AF3 ipTM heatmap ({args.stoich})",
            "fname": f"iptm_heatmap_{args.stoich}.{args.format}",
            "strong": args.iptm_strong,
            "borderline": args.iptm_borderline,
        },
        {
            "metric_col": "ptm",
            "metric_label": "pTM",
            "title": f"AF3 pTM heatmap ({args.stoich})",
            "fname": f"ptm_heatmap_{args.stoich}.{args.format}",
            "strong": None,
            "borderline": None,
        },
    ]

    for p in panels:
        fig_path = out_dir / p["fname"]
        render_heatmap(
            alone_controls=alone_controls or None,
            df=df,
            metric_col=p["metric_col"],
            metric_label=p["metric_label"],
            title=p["title"],
            hap2_vars=hap2_vars,
            dmp_vars=dmp_vars,
            out_path=fig_path,
            cmap_name=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            strong=p["strong"],
            borderline=p["borderline"],
        )
        print(f"[OK] wrote {fig_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
