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
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# CRISPR-P v2.0 on-target (Target) scores keyed by guide NUMBER -- shared
# across delGuide<N> + fsGuide<N>, since both variants come from the same
# gRNA design (and CRISPR-P v2.0 scores the gRNA, not the repair outcome).
# Source: 14_Interaction_Domain_MappingCONFIG.toml [dmp_variants].rows
# guide annotations. Embedded here so the heatmap is self-contained.
# Tool: Liu et al. (2017), CRISPR-P 2.0 (http://crispr.hzau.edu.cn/CRISPR2/).
CRISPR_P_V2_SCORES: dict[int, float] = {
    4:  0.1172,
    16: 0.4235,
    17: 0.8725,
    20: 0.4861,
    37: 0.5374,
    46: 0.4683,
    50: 0.4767,
}

# Matches a guide-anchored variant name and captures (kind, number).
# Both deletion (delGuide<N>) and frameshift (fsGuide<N>) variants of the
# same guide RNA share the same on-target score.
_GUIDE_RE = re.compile(r"^(del|fs)Guide(\d+)$")


def classify_hap2_variant(name: str) -> str:
    """Tag each HAP2 variant for the y-axis group annotation.

    - "WT"                              -> "control"
    - contains "And" (e.g. delPreTMDAndTMD, delEctoAndC)
                                        -> "combo"   (multi-region fusion)
    - everything else                   -> "single"  (one contiguous deletion)

    The heuristic matches the [hap2_variants].rows naming convention in
    14_Interaction_Domain_MappingCONFIG.toml (combo names use the "AndX"
    suffix; singles do not). Override by editing this function if a new
    naming scheme is introduced.
    """
    if name == "WT":
        return "control"
    if "And" in name:
        return "combo"
    return "single"


def hap2_group_spans(hap2_vars: list[str]) -> list[tuple[str, int, int]]:
    """Collapse the per-variant group tags into contiguous (label, start, end)
    spans (end inclusive) in input order. Used to draw the y-axis group
    brackets / labels. Groups must be contiguous in hap2_vars; if a future
    config interleaves them, this returns multiple spans with the same label
    (each rendered as its own bracket).
    """
    spans: list[tuple[str, int, int]] = []
    if not hap2_vars:
        return spans
    cur = classify_hap2_variant(hap2_vars[0])
    start = 0
    for i in range(1, len(hap2_vars)):
        tag = classify_hap2_variant(hap2_vars[i])
        if tag != cur:
            spans.append((cur, start, i - 1))
            cur = tag
            start = i
    spans.append((cur, start, len(hap2_vars) - 1))
    return spans


# Pretty labels for the y-axis group annotation.
_GROUP_DISPLAY = {
    "control": "WT",
    "single":  "Single domain deletions",
    "combo":   "Combined domain deletions",
}


def classify_dmp_variant(name: str) -> str:
    """Tag each DMP variant for the x-axis group annotation. Symmetric to
    classify_hap2_variant.
    - "WT"           -> "control"
    - delGuide<N> / fsGuide<N> (matches _GUIDE_RE) -> "guide"
    - everything else (delN, delC, delTMDcore, ...) -> "domain"
    """
    if name == "WT":
        return "control"
    if _GUIDE_RE.match(name):
        return "guide"
    return "domain"


def dmp_group_spans(dmp_vars: list[str]) -> list[tuple[str, int, int]]:
    """Collapse the per-variant DMP group tags into contiguous (label, start, end)
    spans (end inclusive) in input order. After pair_guide_variants() the
    expected layout is [WT, domain..., guide...]; the spans match that
    ordering and the divider falls between domain and guide segments.
    """
    spans: list[tuple[str, int, int]] = []
    if not dmp_vars:
        return spans
    cur = classify_dmp_variant(dmp_vars[0])
    start = 0
    for i in range(1, len(dmp_vars)):
        tag = classify_dmp_variant(dmp_vars[i])
        if tag != cur:
            spans.append((cur, start, i - 1))
            cur = tag
            start = i
    spans.append((cur, start, len(dmp_vars) - 1))
    return spans


_DMP_GROUP_DISPLAY = {
    "control": "WT",
    "domain":  "Domain deletions variant",
    "guide":   "Guide variants",
}


def format_hap2_tick(name: str) -> str:
    """Wrap multi-domain HAP2 variant names so each deletion segment lives
    on its own line. The convention in [hap2_variants].rows is that combo
    variants concatenate domain tags with "And" (e.g. delPreTMDAndTMD,
    delPreTMDAndTMDAndJuxtaMem, delEctoAndC). Splitting at every "And"
    yields a stack like
        delPreTMD
        AndTMD
        AndJuxtaMem
    which keeps each line short enough to fit in the y-axis margin without
    overlapping the rotated "HAP2 variant" axis label or bleeding into
    adjacent cells. Single-domain variants (no "And") are returned as-is.
    """
    if "And" in name:
        parts = name.split("And")
        return parts[0] + "".join("\nAnd" + p for p in parts[1:])
    return name


def format_dmp_tick(name: str) -> str:
    """Return a multi-line x-tick label. Guide variants get '(CP <score>)'
    appended on a second line (CP = CRISPR-P v2.0 Target Score); both
    delGuide<N> and fsGuide<N> resolve to the same on-target score
    (guide number is the lookup key).
    """
    m = _GUIDE_RE.match(name)
    if m:
        score = CRISPR_P_V2_SCORES.get(int(m.group(2)))
        if score is not None:
            return f"{name}\n(CP {score:.2f})"
    return name


def pair_guide_variants(dmp_vars: list[str]) -> list[str]:
    """Reorder DMP variant names so delGuide<N> and fsGuide<N> for the SAME
    guide-RNA number appear side-by-side on the x-axis. Non-guide variants
    (WT, delC, delN, delTMDcore, ...) keep their input order and are placed
    first; guide variants are then appended in ascending guide-number order,
    deletion before frameshift within each pair.
    """
    non_guides: list[str] = []
    guides: list[tuple[int, int, str]] = []
    for v in dmp_vars:
        m = _GUIDE_RE.match(v)
        if m:
            kind, num = m.group(1), int(m.group(2))
            type_rank = 0 if kind == "del" else 1   # del before fs
            guides.append((num, type_rank, v))
        else:
            non_guides.append(v)
    guides.sort(key=lambda t: (t[0], t[1]))
    return non_guides + [v for _, _, v in guides]


def find_summary_json(slot_dir: Path, model_idx: int) -> Path | None:
    """Return the summary_confidences_<model_idx>.json to read for this slot.

    Normally each slot holds exactly one AF3 bundle (the dedupe step in
    prepare_complexes enforces this). As a defensive fallback for slots
    that still have more than one bundle, prefer the SWISS-template-biased
    prediction: read each candidate's sibling job_request.json and pick
    the one whose 'name' field contains 'swiss' (case-insensitive). This
    codifies the project rule "if a SWISS-template-based model exists,
    prioritise using it" so downstream metrics consistently reflect the
    template-biased prediction even if cleanup ever lags behind.
    """
    cands = sorted(slot_dir.glob(f"*_summary_confidences_{model_idx}.json"))
    if not cands:
        return None
    if len(cands) == 1:
        return cands[0]
    for c in cands:
        prefix = c.name[: -len(f"_summary_confidences_{model_idx}.json")]
        jr = slot_dir / f"{prefix}_job_request.json"
        if not jr.exists():
            continue
        try:
            data = json.loads(jr.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        jobs = data if isinstance(data, list) else [data]
        name = (jobs[0].get("name") or "") if jobs else ""
        if "swiss" in name.lower():
            return c
    return cands[0]


def parse_pair_label(label: str) -> tuple[str, str]:
    if "__" not in label:
        return label, "WT"
    h, d = label.split("__", 1)
    return h, d


def _dmp_self_iptm(data: dict) -> float:
    """Extract the DMP chain's self-ipTM from AF3's chain_pair_iptm matrix.

    AF3 lists chains in input order; in Stage 14 pair-complex jobs the DMP
    is always submitted LAST (monomeric = [HAP2, DMP]; postfusion_like =
    [HAP2_A, HAP2_B, HAP2_C, DMP]). The DMP self-confidence is therefore
    chain_pair_iptm[-1][-1]. Values close to 0 (e.g. < 0.10) indicate
    AF3 found no internal structure for the chain - typical for very
    short frameshift peptides (e.g. fsGuide4 = 10 aa, B-self = 0.01),
    which inflates the overall ipTM with a meaningless point-dock score.
    """
    mat = data.get("chain_pair_iptm")
    if not isinstance(mat, list) or not mat:
        return float("nan")
    last = mat[-1]
    if not isinstance(last, list) or not last:
        return float("nan")
    try:
        return float(last[-1])
    except (TypeError, ValueError):
        return float("nan")


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
            rows.append((pl, h, d, np.nan, np.nan, np.nan, np.nan, ""))
            continue
        try:
            data = json.loads(js.read_text())
        except Exception as exc:  # noqa: BLE001
            print(f"[WARN] could not parse {js}: {exc}", file=sys.stderr)
            rows.append((pl, h, d, np.nan, np.nan, np.nan, np.nan, str(js.name)))
            continue
        iptm = float(data.get("iptm", np.nan))
        ptm = float(data.get("ptm", np.nan))
        rs_raw = data.get("ranking_score")
        if rs_raw is None and not (np.isnan(iptm) or np.isnan(ptm)):
            rs = 0.8 * iptm + 0.2 * ptm
        else:
            rs = float(rs_raw) if rs_raw is not None else np.nan
        dmp_self = _dmp_self_iptm(data)
        rows.append((pl, h, d, iptm, ptm, rs, dmp_self, js.name))
    return pd.DataFrame(
        rows,
        columns=["pair", "hap2", "dmp", "iptm", "ptm", "ranking_score",
                 "dmp_self_iptm", "source_json"],
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
    short_chain_threshold: float | None = None,
) -> None:
    """Render a single HAP2 x DMP heatmap for one metric column.

    strong / borderline (optional): if both provided, cells >= strong get **,
    cells >= borderline get *, and threshold lines are drawn on the colorbar.
    Pass None to omit threshold annotations (recommended for pTM / ranking_score
    where the ipTM 0.60/0.80 thresholds do not transfer cleanly).

    short_chain_threshold (optional, ipTM panel only): cells whose DMP-self
    ipTM (chain_pair_iptm[-1][-1] in the AF3 JSON) falls below this value
    are flagged as "null†" with the value parenthesised, because AF3's
    interface confidence on a chain it cannot fold internally is a
    point-dock artifact (e.g. fsGuide4 = 10 aa, B-self = 0.01, overall
    ipTM = 0.58). Footnote is rendered below the heatmap.
    """
    mat = np.full((len(hap2_vars), len(dmp_vars)), np.nan)
    short_mat = np.zeros((len(hap2_vars), len(dmp_vars)), dtype=bool)
    for _, row in df.iterrows():
        if row.hap2 in hap2_vars and row.dmp in dmp_vars:
            i = hap2_vars.index(row.hap2)
            j = dmp_vars.index(row.dmp)
            mat[i, j] = row[metric_col]
            if (short_chain_threshold is not None
                    and "dmp_self_iptm" in df.columns
                    and not np.isnan(row.dmp_self_iptm)
                    and row.dmp_self_iptm < short_chain_threshold):
                short_mat[i, j] = True

    # Fixed cell size so main-heatmap squares and control-strip squares
    # render identical on disk; height/width derived from cell counts.
    cell_in = 0.55
    margin_w = 3.2   # y-labels + outer group bracket + colorbar + tight side padding
    margin_h = 2.4   # title + 2-line rotated x-tick labels + x-group labels
                     # ("Domain deletions variant" / "Guide variants") +
                     # "DMP variant (CP = ...)" parent xlabel.
    fig_w = len(dmp_vars) * cell_in + margin_w
    # Vertical room reserved for the manually-placed controls strip below:
    # 2.4 in = small gap + strip title + 1 cell + strip x-labels.
    ctrl_block_in = 2.4 if alone_controls else 0.0
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
    ax.set_yticklabels([format_hap2_tick(v) for v in hap2_vars], fontsize=9)
    # "DMP variant (CP = CRISPR-P v2.0 Target Score)" xlabel acts as the
    # parent x-axis label; "Domain deletions variant" / "Guide variants"
    # group labels (placed via the x_spans block below) sit between the
    # rotated tick labels and this xlabel. Large labelpad opens the
    # space so the three text rows are visually separated.
    ax.set_xlabel("DMP variant (CP = CRISPR-P v2.0 Target Score)",
                  fontsize=11, fontweight="bold", labelpad=38)
    # The default ylabel is suppressed; "HAP2 variant" is placed manually
    # further left than the group brackets and centered between the singles
    # and combos group labels (see block after invert_yaxis below).
    ax.set_ylabel("")
    ax.set_title(title, fontsize=11)

    # WT__WT belongs in the bottom-left corner. Input order is
    # [WT, singles..., combos...]; matplotlib's default imshow origin is "upper"
    # (row 0 at top), so we invert the y-axis to flip WT (index 0) down.
    # Group dividers and labels are placed AFTER inversion so they sit on the
    # correct side of each span (top edge of each higher-index group).
    ax.invert_yaxis()

    # Horizontal dividers + outer-left group brackets (single / combo / control)
    spans = hap2_group_spans(hap2_vars)
    # Draw separators between consecutive groups. With y-axis inverted, axhline
    # still takes data y-coords, so boundary at index i.5 sits between rows
    # i and i+1 regardless of orientation.
    for k in range(1, len(spans)):
        boundary = spans[k][1] - 0.5  # = (end of previous span) + 0.5
        ax.axhline(boundary, color="black", lw=0.8, linestyle="-", alpha=0.55)
    # Vertical label for each span placed outside the y-tick labels. Uses
    # axes-fraction x (left of plot) and data-y (variant index) via
    # get_yaxis_transform(). Skip the singleton "control" row (label "WT"
    # is redundant with the y-tick label).
    for tag, start, end in spans:
        if tag == "control":
            continue
        mid = (start + end) / 2.0
        ax.text(
            -0.22, mid, _GROUP_DISPLAY.get(tag, tag),
            rotation=90,
            ha="center", va="center",
            transform=ax.get_yaxis_transform(),
            fontsize=10, fontweight="bold",
            color="#222222",
        )

    # "HAP2 variant" parent label, placed further left than the per-span
    # brackets and centered between the singles and combos group midpoints
    # so it visually sits one level above them.
    non_ctrl_centers = [(s[1] + s[2]) / 2.0 for s in spans if s[0] != "control"]
    if non_ctrl_centers:
        y_center = sum(non_ctrl_centers) / len(non_ctrl_centers)
        ax.text(
            -0.32, y_center, "HAP2 variant",
            rotation=90, ha="center", va="center",
            transform=ax.get_yaxis_transform(),
            fontsize=12, fontweight="bold", color="black",
        )

    # X-axis grouping: domain-deletion DMP variants vs guide variants.
    # Vertical dividers (axvline) mirror the axhline group separators on
    # the y-axis. Group labels are placed JUST BELOW the rotated x-tick
    # labels (and above the "DMP variant" xlabel, which set_xlabel pushed
    # further down via labelpad=22).
    x_spans = dmp_group_spans(dmp_vars)
    for k in range(1, len(x_spans)):
        boundary = x_spans[k][1] - 0.5
        ax.axvline(boundary, color="black", lw=0.8, linestyle="-", alpha=0.55)
    for tag, start, end in x_spans:
        if tag == "control":
            continue
        mid = (start + end) / 2.0
        # y=-0.14 places these just below the rotated x-tick labels (which
        # end around y=-0.12); the bold "DMP variant (CP = ...)" parent
        # xlabel sits ~0.10 below this thanks to labelpad, giving a clear
        # vertical gap between the two text rows.
        ax.text(
            mid, -0.14, _DMP_GROUP_DISPLAY.get(tag, tag),
            rotation=0,
            ha="center", va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=10, fontweight="bold",
            color="#222222",
        )

    annotate_thresholds = strong is not None and borderline is not None
    any_short = False
    for i in range(len(hap2_vars)):
        for j in range(len(dmp_vars)):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, "n/a", ha="center", va="center", color="black", fontsize=8)
                continue
            if short_mat[i, j]:
                # Short-chain artifact: AF3 cannot fold the DMP internally,
                # so the interface score is a point-dock artifact. Show the
                # raw value parenthesised plus a dagger so the reader sees
                # both the inflated number and the warning.
                any_short = True
                ax.text(
                    j, i, f"null†\n({v:.2f})",
                    ha="center", va="center",
                    color="black", fontsize=8, fontstyle="italic",
                )
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

    # Footnote is drawn after the controls strip (if any) so it can be
    # positioned below the strip's rotated tick labels instead of overlapping
    # with them at the figure bottom.
    footnote_text = (
        (f"† DMP self-ipTM < {short_chain_threshold:.2f} (chain too "
         "short for AF3 to fold internally); overall ipTM is a "
         "point-dock artifact and should not be read as interface affinity.")
        if (any_short and short_chain_threshold is not None) else None
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

        # Main heatmap axes bbox in figure coords. get_position() returns the
        # axes rectangle regardless of axis inversion, which im.get_window_extent
        # does not: after ax.invert_yaxis() the image bbox can land near the
        # top of the figure and the controls strip overlaps the heatmap.
        fig.canvas.draw()  # ensure transforms (and tight_layout) are finalised
        main_bbox = ax.get_position()
        main_left = main_bbox.x0
        main_width = main_bbox.width
        main_bottom = main_bbox.y0

        # Centre the controls strip horizontally under the heatmap data area.
        ctrl_left = main_left + (main_width - ctrl_width_fig) / 2
        # Vertical gap below main heatmap: room for x-tick labels + the
        # "Domain deletions variant" / "Guide variants" group labels + the
        # "DMP variant" parent xlabel + a small breathing gap before the
        # controls-strip title.
        gap_fig = 1.95 / fig_h
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

    if footnote_text is not None:
        # Place the footnote below everything else so it does not collide
        # with the controls-strip rotated tick labels. With a controls strip
        # present we anchor to the lowest rendered tick label; otherwise we
        # fall back to the figure-bottom slot.
        if ax_ctrl is not None:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            tick_lbls = [t for t in ax_ctrl.get_xticklabels() if t.get_text()]
            if tick_lbls:
                lowest_px = min(t.get_window_extent(renderer).y0 for t in tick_lbls)
                lowest_fig = lowest_px / fig.bbox.height
                footnote_y = lowest_fig - (0.18 / fig_h)
            else:
                footnote_y = ax_ctrl.get_position().y0 - (0.55 / fig_h)
        else:
            footnote_y = 0.01
        fig.text(
            0.01, footnote_y, footnote_text,
            fontsize=8, color="#555555", ha="left", va="top",
        )

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
    ap.add_argument("--short-chain-self-iptm", type=float, default=0.10,
                    help="DMP self-ipTM (chain_pair_iptm[-1][-1]) below this is flagged "
                         "as a short-chain artifact in the ipTM panel; set to 0 to disable. "
                         "Default 0.10 catches very short frameshift peptides like fsGuide4 "
                         "(B-self ~0.01) without touching real partial-fold variants like "
                         "fsGuide17/50 (B-self ~0.2-0.3).")
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
    dmp_vars = pair_guide_variants(args.dmp_variants.split())
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

    short_chain_thr = args.short_chain_self_iptm if args.short_chain_self_iptm > 0 else None
    panels = [
        {
            "metric_col": "iptm",
            "metric_label": "ipTM",
            "title": f"AF3 ipTM heatmap ({args.stoich})",
            "fname": f"iptm_heatmap_{args.stoich}.{args.format}",
            "strong": args.iptm_strong,
            "borderline": args.iptm_borderline,
            "short_chain_threshold": short_chain_thr,
        },
        {
            "metric_col": "ptm",
            "metric_label": "pTM",
            "title": f"AF3 pTM heatmap ({args.stoich})",
            "fname": f"ptm_heatmap_{args.stoich}.{args.format}",
            "strong": None,
            "borderline": None,
            "short_chain_threshold": None,  # pTM is whole-complex, not interface-only
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
            short_chain_threshold=p["short_chain_threshold"],
        )
        print(f"[OK] wrote {fig_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
