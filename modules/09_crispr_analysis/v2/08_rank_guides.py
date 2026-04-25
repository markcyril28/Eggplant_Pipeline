#!/usr/bin/env python3
"""
Module: 08_rank_guides.py
Stage [8] — Composite KO score computation and guide ranking.

Computes a weighted composite KO score from five component signals:
  w_ontarget   : normalised mean on-target efficiency
  w_frameshift : predicted frameshift fraction (inDelphi / Lindel)
  w_nmd        : NMD susceptibility score
  w_offtarget  : inverted off-target penalty (1 - normalised CFD sum)
  w_domain     : functional domain hit bonus

Outputs:
  - ranked_guides.tsv   : full table sorted by composite_ko_score (desc)
  - top_guides.tsv      : top-N guides only
  - ko_score_plot.png   : bar chart of top-N guides coloured by component

Usage:
    python3 08_rank_guides.py \\
        --input     <nmd.tsv>           \\
        --outdir    <output_dir>        \\
        --weights   '{"w_ontarget":0.30,"w_frameshift":0.25,"w_nmd":0.20,"w_offtarget":0.15,"w_domain":0.10}' \\
        --top-n     10                  \\
        --output-format tsv             \\
        --dpi       600                 \\
        --format    png
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


DEFAULT_WEIGHTS = {
    "w_ontarget":   0.30,
    "w_frameshift": 0.25,
    "w_nmd":        0.20,
    "w_offtarget":  0.15,
    "w_domain":     0.10,
}

# Absolute cap for CFD off-target sum normalisation.
# Rationale: per-batch max-normalisation made a guide's off-target penalty
# depend on the maximum CFD in its batch — two identical guides in different
# genes would receive different composite scores. Using an absolute cap
# (5.0 = "essentially unusable" per the Haeussler/CRISPOR off-target
# specificity guidance) yields reproducible, gene-independent scores.
# Override per-call with --cfd-sum-cap.
CFD_SUM_CAP_DEFAULT = 5.0

# Known on-target scorer columns grouped by output scale. Stage [8] normalises
# each scorer to [0, 1] independently, then averages across all available
# scores — the previous "break-on-first-hit" behaviour silently discarded
# CRISPRon/DeepSpCas9 rescores whenever a Doench/Moreno-Mateos score was
# present, wasting stage [2] compute and biasing the composite toward a
# single human-trained model.
# Order inside each tuple is informational only (no longer affects scoring).
ONTARGET_SCORE_COLS_0_100 = (
    "doench2016OnTarget", "Doench2016OnTarget",
    "Doench '16-Score", "DoenchScore",
    "Moreno-Mateos-Score",
    "crisprP_score",            # CRISPR-P v2.0 (modified Xu 2015, 0-100)
)
ONTARGET_SCORE_COLS_0_1 = (
    "crisprOn_score", "DeepSpCas9_score",
    # Plant-trained rescorers (stage [2]): DeepCRISPR emits 0-1 activity;
    # CRISPR-Local emits 0-1 efficacy in its Python helper output.
    "DeepCRISPR_score", "CRISPR_Local_score",
)
# Generic fallback column — scale detected heuristically at read time.
ONTARGET_SCORE_COL_GENERIC = "on_target_score"


def parse_args():
    p = argparse.ArgumentParser(description="Composite KO score & guide ranking")
    # --input is required only for per-gene mode (single input TSV).
    # --aggregate-dir switches to aggregate mode: read all per-gene ranked TSVs
    # in the directory, concatenate + re-sort, and emit aggregate tables/plot.
    p.add_argument("--input",          default=None)
    p.add_argument("--aggregate-dir",  default=None,
                   help="Directory containing per-gene *.ranked.tsv files. "
                        "When set, outputs aggregate ranked_guides/top_guides/plot.")
    p.add_argument("--outdir",         required=True)
    p.add_argument("--weights",        default=json.dumps(DEFAULT_WEIGHTS),
                   help="JSON dict of component weights (must sum to 1.0)")
    p.add_argument("--top-n",          type=int,   default=10)
    p.add_argument("--output-format",  default="tsv",  choices=["tsv", "csv", "xlsx"])
    p.add_argument("--dpi",            type=int,   default=600)
    p.add_argument("--format",         default="png",  choices=["png", "svg", "pdf"])
    p.add_argument("--cfd-sum-cap",    type=float, default=CFD_SUM_CAP_DEFAULT,
                   help=("Absolute CFD-sum cap for off-target normalisation (default 5.0). "
                         "Guides with cfd_sum >= cap receive c_offtarget = 0. Replaces the "
                         "previous per-batch max-normalisation, which made the penalty "
                         "depend on which siblings were in the run."))
    p.add_argument("--overwrite",      action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


# ─── Safe float parse ─────────────────────────────────────────────────────────

def safe_float(val, default=0.0) -> float:
    if val in (None, "", "NA", "N/A"):
        return default
    try:
        f = float(val)
        return f if math.isfinite(f) else default
    except (ValueError, TypeError):
        return default


# ─── Score computation ────────────────────────────────────────────────────────

def _detect_generic_denom(rows: list[dict]) -> float | None:
    """Column-wide scale detection for the `on_target_score` fallback column.

    The per-row heuristic (v > 1.5 ⇒ 0-100) is ambiguous: a legitimate 0-100
    score of 1.0 would be misread as 100 %, and a 0-1 score of 0.5 as 50 %.
    The scale is a property of the emitting scorer, so it is detected once
    per table by scanning the column's max value:

      - max > 1.5  ⇒ clearly 0-100 (some value exceeds the 0-1 ceiling)
      - max ≤ 1.0  ⇒ clearly 0-1
      - 1.0 < max ≤ 1.5 ⇒ ambiguous band — default to 0-1 (conservative:
        assumes high activity rather than inflating a 1 % score to 1.0)
        and log a WARN so the ambiguity is visible.

    Returns the denominator to divide by (100.0, 1.0) or None when no
    numeric value is present in the column.
    """
    vals: list[float] = []
    for r in rows:
        v = safe_float(r.get(ONTARGET_SCORE_COL_GENERIC))
        if v > 0:
            vals.append(v)
    if not vals:
        return None
    col_max = max(vals)
    if col_max > 1.5:
        return 100.0
    if col_max > 1.0:
        _log(f"[08_rank] '{ONTARGET_SCORE_COL_GENERIC}' column max={col_max:.4f} "
             "lies in the 1.0-1.5 ambiguous band; treating as 0-1 scale. "
             "Rename the column to an explicit scorer name (see "
             "ONTARGET_SCORE_COLS_0_100 / _0_1) to disambiguate.",
             level="WARN")
    return 1.0


def _normalised_ontarget_scores(row: dict, generic_denom: float | None = None) -> list[float]:
    """Return ALL available on-target scores from `row`, each scaled to [0, 1].

    Unlike the prior break-on-first-hit implementation, this collects every
    recognised scorer so CRISPRon / DeepSpCas9 rescores contribute even when
    a Doench/Moreno-Mateos column is present. Each scorer is normalised by
    its own native scale (0-100 vs 0-1), not by a batch-wide max — so the
    value of a guide's on-target component no longer depends on its
    siblings in the run.

    `generic_denom` is the column-wide scale for the ONTARGET_SCORE_COL_GENERIC
    fallback, produced by `_detect_generic_denom(rows)`. Pass None when
    calling this function in isolation (single-row contexts): a per-row
    heuristic kicks in with a WARN-level ambiguity note.
    """
    values: list[float] = []
    for col in ONTARGET_SCORE_COLS_0_100:
        v = safe_float(row.get(col))
        if v > 0:
            values.append(min(v / 100.0, 1.0))
    for col in ONTARGET_SCORE_COLS_0_1:
        v = safe_float(row.get(col))
        if v > 0:
            values.append(min(v, 1.0))
    # Generic fallback: accepted only when no other scorer fired.
    if not values:
        v = safe_float(row.get(ONTARGET_SCORE_COL_GENERIC))
        if v > 0:
            if generic_denom is not None:
                denom = generic_denom
            else:
                # Isolated per-row call — fall back to the old heuristic but
                # surface the ambiguity. Prefer the compute_scores() path.
                denom = 100.0 if v > 1.5 else 1.0
            values.append(min(v / denom, 1.0))
    return values


def compute_scores(rows: list[dict], cfd_sum_cap: float = CFD_SUM_CAP_DEFAULT) -> list[dict]:
    """Add component score columns to each row (normalised, 0-1 range).

    CFD-sum normalisation uses an absolute cap (cfd_sum_cap, default 5.0)
    rather than per-batch max, so composite scores are reproducible across
    runs and gene groups.
    """
    if cfd_sum_cap <= 0:
        _log(f"[08_rank] cfd_sum_cap must be > 0 (got {cfd_sum_cap}); "
             f"falling back to default {CFD_SUM_CAP_DEFAULT}.", level="WARN")
        cfd_sum_cap = CFD_SUM_CAP_DEFAULT

    # Determine the generic-fallback column's scale once from the whole table
    # rather than per row, so a 0-100 score of 1.0 is not misread as 100 %.
    generic_denom = _detect_generic_denom(rows)

    for row in rows:
        # ── on-target component: mean across ALL available, independently
        #    normalised scorers (replaces prior break-on-first-hit). ────────
        ot_values = _normalised_ontarget_scores(row, generic_denom=generic_denom)
        if ot_values:
            c_ontarget = sum(ot_values) / len(ot_values)
            # Stash scorer count so downstream auditing can see whether the
            # composite leaned on one model or several.
            row["c_ontarget_n_sources"] = str(len(ot_values))
        else:
            c_ontarget = 0.0
            row["c_ontarget_n_sources"] = "0"

        # ── frameshift component ───────────────────────────────────────────
        c_frameshift = safe_float(row.get("best_fs_frac"))
        c_frameshift = min(c_frameshift, 1.0)

        # ── NMD component ──────────────────────────────────────────────────
        nmd = row.get("nmd_summary", "")
        if nmd == "nmd_predicted":
            c_nmd = 1.0
        elif nmd == "nmd_escape":
            c_nmd = 0.5
        else:
            c_nmd = 0.0

        # ── off-target component (inverted penalty, absolute cap) ─────────
        cfd = safe_float(row.get("cfd_sum"))
        c_offtarget = max(0.0, 1.0 - min(cfd, cfd_sum_cap) / cfd_sum_cap)

        # ── domain component ──────────────────────────────────────────────
        c_domain = 1.0 if row.get("any_domain_hit", "").lower() == "true" else 0.0

        row["c_ontarget"]   = f"{c_ontarget:.4f}"
        row["c_frameshift"] = f"{c_frameshift:.4f}"
        row["c_nmd"]        = f"{c_nmd:.4f}"
        row["c_offtarget"]  = f"{c_offtarget:.4f}"
        row["c_domain"]     = f"{c_domain:.4f}"

    return rows


def _validate_weights(weights: dict, tol: float = 1e-6) -> dict:
    """Assert that weight values sum to 1.0 (within tolerance) and that every
    component has a non-negative weight. Returns the validated dict, falling
    back to DEFAULT_WEIGHTS for any missing keys.
    """
    w = {k: weights.get(k, DEFAULT_WEIGHTS[k]) for k in DEFAULT_WEIGHTS}
    # Coerce to float up front so a string in the JSON doesn't break later math.
    for k in list(w.keys()):
        try:
            w[k] = float(w[k])
        except (TypeError, ValueError):
            _log(f"[08_rank] Weight {k!r}={weights.get(k)!r} is not numeric; "
                 f"falling back to default {DEFAULT_WEIGHTS[k]}.", level="WARN")
            w[k] = DEFAULT_WEIGHTS[k]
    for k, v in w.items():
        if v < 0:
            _log(f"[08_rank] Weight {k}={v} is negative; clamped to 0.", level="WARN")
            w[k] = 0.0
    total = sum(w.values())
    if abs(total - 1.0) > tol:
        # Auto-renormalise and log loudly — refusing to run would abort a long
        # pipeline run over a config typo; renormalising lets the run finish
        # while making the drift obvious in the log and in a new column.
        _log(f"[08_rank] Composite weights sum to {total:.6f}, expected 1.0 "
             f"(±{tol:.0e}). Renormalising in place. Original weights: {w}",
             level="ERROR" if abs(total - 1.0) > 0.05 else "WARN")
        if total <= 0:
            _log("[08_rank] Weight sum ≤ 0; cannot renormalise. Reverting to "
                 "DEFAULT_WEIGHTS.", level="ERROR")
            return dict(DEFAULT_WEIGHTS)
        w = {k: v / total for k, v in w.items()}
    return w


def apply_weights(rows: list[dict], weights: dict) -> list[dict]:
    # Validated + renormalised weights (see _validate_weights).
    w = _validate_weights(weights)
    for row in rows:
        score = (
            w["w_ontarget"]   * safe_float(row.get("c_ontarget"))
            + w["w_frameshift"] * safe_float(row.get("c_frameshift"))
            + w["w_nmd"]        * safe_float(row.get("c_nmd"))
            + w["w_offtarget"]  * safe_float(row.get("c_offtarget"))
            + w["w_domain"]     * safe_float(row.get("c_domain"))
        )
        row["composite_ko_score"] = f"{score:.4f}"
    return rows


# ─── Visualisation ────────────────────────────────────────────────────────────

def plot_top_guides(rows: list[dict], outdir: Path, top_n: int,
                    dpi: int, fmt: str, stem: str = "",
                    filename: str | None = None) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        _log("[08_rank] matplotlib not installed — skipping plot.", level="WARN")
        return

    top = rows[:top_n]
    if not top:
        return

    guide_labels = []
    for r in top:
        # sgRNA_id is the canonical CRISPR-P v2.0 identifier; keep the CRISPOR
        # alternatives ahead of generic ones. Without sgRNA_id in this list the
        # log line "Best guide: ?" appears for every CRISPR-P-only run.
        for col in ("guideId", "guide_id", "sgRNA_id", "name", "ID", "targetSeq"):
            if r.get(col):
                guide_labels.append(str(r[col])[:25])
                break
        else:
            guide_labels.append("guide")

    components  = ["c_ontarget", "c_frameshift", "c_nmd", "c_offtarget", "c_domain"]
    comp_labels = ["On-target", "Frameshift", "NMD", "Off-target (inv)", "Domain"]
    colors      = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336"]

    data = [[safe_float(r.get(c)) for r in top] for c in components]

    x       = range(len(top))
    bar_w   = 0.6
    bottoms = [0.0] * len(top)

    fig, ax = plt.subplots(figsize=(max(8, len(top) * 0.9), 5))
    for comp_data, label, color in zip(data, comp_labels, colors):
        ax.bar(x, comp_data, bar_w, bottom=bottoms, label=label, color=color, alpha=0.85)
        bottoms = [b + d for b, d in zip(bottoms, comp_data)]

    ax.set_xticks(list(x))
    ax.set_xticklabels(guide_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Composite KO score (stacked components)")
    ax.set_title(f"Top {top_n} CRISPR guides by KO score")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()

    plot_path = outdir / (filename or f"{stem}.top{top_n}_ko_score.{fmt}")
    plt.savefig(plot_path, dpi=dpi)
    plt.close()
    _log(f"[08_rank] Plot saved: {plot_path}", level="INFO")


# ─── Main ─────────────────────────────────────────────────────────────────────

def _write_table(path: Path, data: list[dict], fields: list[str], fmt: str) -> str:
    """Write rows to path in the requested format; return the format actually used."""
    if fmt == "xlsx":
        try:
            import openpyxl
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.append(fields)
            for r in data:
                ws.append([r.get(f, "") for f in fields])
            wb.save(path)
            return "xlsx"
        except ImportError:
            _log("[08_rank] openpyxl not installed — falling back to TSV.", level="WARN")
            path = path.with_suffix(".tsv")
            fmt = "tsv"
    delim = "," if fmt == "csv" else "\t"
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter=delim,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)
    return fmt


def _run_per_gene(args, weights):
    """Per-gene ranking: writes {GENE}.ranked.tsv and {GENE}.top{N}.tsv to --outdir."""
    inpath = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Gene stem: prefix before the first dot (e.g. SMEL5_01g008730.nmd → SMEL5_01g008730).
    stem = inpath.stem.split(".")[0]
    ranked_path = outdir / f"{stem}.ranked.{args.output_format}"
    top_path    = outdir / f"{stem}.top{args.top_n}.{args.output_format}"

    if not args.overwrite and ranked_path.exists():
        _log(f"[08_rank] Skipping (overwrite=false): {ranked_path}", level="WARN")
        return

    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    rows = compute_scores(rows, cfd_sum_cap=args.cfd_sum_cap)
    rows = apply_weights(rows, weights)
    rows.sort(key=lambda r: safe_float(r.get("composite_ko_score")), reverse=True)

    new_cols = ["c_ontarget", "c_ontarget_n_sources",
                "c_frameshift", "c_nmd", "c_offtarget", "c_domain",
                "composite_ko_score"]
    new_fields = fields + [c for c in new_cols if c not in fields]

    _write_table(ranked_path, rows,                new_fields, args.output_format)
    _write_table(top_path,    rows[:args.top_n],   new_fields, args.output_format)

    plot_top_guides(rows, outdir, args.top_n, args.dpi, args.format,
                    stem=stem,
                    filename=f"{stem}.top{args.top_n}_ko_score.{args.format}")

    _log(f"[08_rank] {len(rows)} guides ranked; top-{args.top_n} -> {top_path}", level="INFO")
    if rows:
        best = rows[0]
        best_id = next((best.get(c) for c in ("guideId","guide_id","sgRNA_id","name") if best.get(c)), "?")
        _log(f"[08_rank] Best guide: {best_id}  KO score={best.get('composite_ko_score','?')}", level="INFO")


def _run_aggregate(args):
    """Aggregate mode: concatenate all per-gene *.ranked.<fmt> files under
    --aggregate-dir, re-sort, and emit ranked_guides.<fmt> / top_guides.<fmt> /
    top{N}_ko_score.png at --outdir (stage root)."""
    agg_dir = Path(args.aggregate_dir)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Sentinel check: skip when overwrite=false and the aggregate output already
    # exists. Mirrors the per-gene guard in _run_per_gene.
    ranked_path = outdir / f"ranked_guides.{args.output_format}"
    if not args.overwrite and ranked_path.exists():
        _log(f"[08_rank] Aggregate: skipping (overwrite=false): {ranked_path}", level="WARN")
        return

    # Glob matches the configured output format. If xlsx was requested but
    # _write_table fell back to .tsv (openpyxl absent), try .tsv as a fallback.
    fmt = args.output_format
    per_gene_files = sorted(agg_dir.glob(f"*.ranked.{fmt}"))
    if not per_gene_files and fmt == "xlsx":
        per_gene_files = sorted(agg_dir.glob("*.ranked.tsv"))
    if not per_gene_files:
        _log(f"[08_rank] No *.ranked.{fmt} under {agg_dir} — aggregate skipped.", level="WARN")
        return

    all_rows: list[dict] = []
    fields:   list[str]  = []
    for f in per_gene_files:
        with open(f, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            # Access fieldnames before the row loop so header-only files still
            # populate `fields` (DictReader only reads the header on first access).
            file_fields = reader.fieldnames
            for row in reader:
                all_rows.append(row)
            if file_fields and not fields:
                fields = list(file_fields)

    if not all_rows:
        _log(f"[08_rank] Per-gene files under {agg_dir} are all empty — aggregate skipped.",
             level="WARN")
        return

    all_rows.sort(key=lambda r: safe_float(r.get("composite_ko_score")), reverse=True)

    top_path = outdir / f"top{args.top_n}_guides.{args.output_format}"
    _write_table(ranked_path, all_rows,              fields, args.output_format)
    _write_table(top_path,    all_rows[:args.top_n], fields, args.output_format)

    plot_top_guides(all_rows, outdir, args.top_n, args.dpi, args.format,
                    filename=f"top{args.top_n}_ko_score.{args.format}")

    _log(f"[08_rank] Aggregate: {len(all_rows)} guides from {len(per_gene_files)} genes "
         f"-> {ranked_path}", level="INFO")


def main():
    args = parse_args()

    try:
        weights = json.loads(args.weights)
    except json.JSONDecodeError:
        _log("[08_rank] Invalid JSON for --weights; using defaults.", level="WARN")
        weights = DEFAULT_WEIGHTS

    if args.aggregate_dir:
        _run_aggregate(args)
    elif args.input:
        _run_per_gene(args, weights)
    else:
        _log("[08_rank] Either --input or --aggregate-dir is required.", level="ERROR")
        sys.exit(2)


if __name__ == "__main__":
    main()
