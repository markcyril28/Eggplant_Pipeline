#!/usr/bin/env python3
"""
combined_plots.py — cross-step aggregation for the GROMACS PPI pipeline.

Walks every completed simulation/analysis step under <run_dir>, joins
per-structure metrics by canonical structure stem, normalizes the values to
[0, 1] (with lower-is-better axes inverted so "outer = better" everywhere),
and emits ONE combined radar + ranking covering all structures across all
steps.

Outputs (under --output-dir):
  combined_radar_all_steps.<fmt>
  combined_ranking_all_steps.<fmt>
  combined_metrics.csv              (raw per-structure x per-step values)
  combined_metrics_normalized.csv   (0-1 values + composite score)

Step layouts understood:
  quick_stability         <step_dir>/<workdir>/structure_N/metrics.json
                          (PDB stem resolved via <workdir>/MANIFEST.txt)
  compare_chain_stability <step_dir>/<workdir>/structure_N/statistics/md_statistics.json
  production_md           same as compare_chain_stability
  batch_comparison        <step_dir>/<dataset>/<stem>/metrics.json
                          (stem is the PDB basename without extension)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Combined radar axis spec.
#
# Each entry: (axis_label, [step_priority], json_key, better_direction)
# The first step in `priority` that has the metric (with >=2 distinct values)
# wins. lower-is-better axes are inverted on normalization so the outer ring
# always corresponds to "better".
# ---------------------------------------------------------------------------

COMBINED_AXES: List[Tuple[str, List[str], str, str]] = [
    ("EM Energy",   ["batch_comparison", "quick_stability"], "potential",       "lower"),
    ("EM H-bonds",  ["batch_comparison", "quick_stability"], "hbonds",          "higher"),
    ("EM Contacts", ["batch_comparison", "quick_stability"], "contacts",        "higher"),
    ("EM SASA",     ["batch_comparison", "quick_stability"], "sasa",            "lower"),
    ("MD RMSD",     ["production_md", "compare_chain_stability"], "rmsd_mean",     "lower"),
    ("MD H-bonds",  ["production_md", "compare_chain_stability"], "hbonds_mean",   "higher"),
    ("MD Total IE", ["production_md", "compare_chain_stability"], "total_ie_mean", "lower"),
    ("MD Min Dist", ["production_md", "compare_chain_stability"], "mindist_mean",  "lower"),
]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _read_json(path: Path) -> Dict:
    try:
        with open(path) as f:
            txt = f.read().strip()
        return json.loads(txt) if txt else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _flatten(data: Dict) -> Dict[str, float]:
    """Flatten nested dicts into scalar key/value pairs (one level deep)."""
    flat: Dict[str, float] = {}
    for k, v in data.items():
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            flat[k] = float(v)
        elif isinstance(v, dict):
            for sk, sv in v.items():
                if isinstance(sv, bool):
                    continue
                if isinstance(sv, (int, float)):
                    flat[f"{k}_{sk}"] = float(sv)
    return flat


def _read_manifest(workdir: Path) -> Dict[str, str]:
    """Parse MANIFEST.txt (lines: 'structure_N  pdb_basename')."""
    out: Dict[str, str] = {}
    manifest = workdir / "MANIFEST.txt"
    if not manifest.exists():
        return out
    try:
        for line in manifest.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    except OSError:
        pass
    return out


def _stem(name: str) -> str:
    """Drop .pdb / .cif / .gz suffix to get a canonical structure name."""
    n = name
    for ext in (".gz",):
        if n.endswith(ext):
            n = n[: -len(ext)]
    for ext in (".pdb", ".cif"):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n


# ---------------------------------------------------------------------------
# Per-step collectors: each returns {structure_stem: {metric_key: value}}
# ---------------------------------------------------------------------------

def collect_em_workdir(step_dir: Path) -> Dict[str, Dict[str, float]]:
    """Step 1 layout: <step_dir>/<workdir>/structure_N/metrics.json."""
    out: Dict[str, Dict[str, float]] = {}
    if not step_dir.is_dir():
        return out
    for workdir in sorted(step_dir.iterdir()):
        if not workdir.is_dir():
            continue
        manifest = _read_manifest(workdir)
        if not manifest:
            continue
        for sn, pdb_basename in manifest.items():
            mfile = workdir / sn / "metrics.json"
            data = _read_json(mfile)
            if data:
                out[_stem(pdb_basename)] = _flatten(data)
    return out


def collect_md_workdir(step_dir: Path) -> Dict[str, Dict[str, float]]:
    """Steps 2/5 layout: <step_dir>/<workdir>/structure_N/statistics/md_statistics.json."""
    out: Dict[str, Dict[str, float]] = {}
    if not step_dir.is_dir():
        return out
    for workdir in sorted(step_dir.iterdir()):
        if not workdir.is_dir():
            continue
        manifest = _read_manifest(workdir)
        if not manifest:
            continue
        for sn, pdb_basename in manifest.items():
            stats = workdir / sn / "statistics" / "md_statistics.json"
            data = _read_json(stats)
            if data:
                out[_stem(pdb_basename)] = _flatten(data)
    return out


def collect_batch_dataset(step_dir: Path) -> Dict[str, Dict[str, float]]:
    """Step 4 layout: <step_dir>/<dataset>/<stem>/metrics.json."""
    out: Dict[str, Dict[str, float]] = {}
    if not step_dir.is_dir():
        return out
    for ds in sorted(step_dir.iterdir()):
        if not ds.is_dir():
            continue
        for struct_dir in sorted(ds.iterdir()):
            if not struct_dir.is_dir():
                continue
            mfile = struct_dir / "metrics.json"
            data = _read_json(mfile)
            if data:
                out[_stem(struct_dir.name)] = _flatten(data)
    return out


COLLECTORS: Dict[str, Callable[[Path], Dict[str, Dict[str, float]]]] = {
    "quick_stability": collect_em_workdir,
    "compare_chain_stability": collect_md_workdir,
    "batch_comparison": collect_batch_dataset,
    "production_md": collect_md_workdir,
}


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(
    run_dir: Path,
    step_prefixes: Dict[str, str],
    source_steps: Optional[List[str]],
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Returns {step_name: {structure_stem: {metric_key: value}}}."""
    auto = not source_steps
    per_step: Dict[str, Dict[str, Dict[str, float]]] = {}
    for step, fn in COLLECTORS.items():
        if not auto and step not in source_steps:
            continue
        prefix = step_prefixes.get(step, step)
        sd = run_dir / prefix
        if not sd.is_dir():
            continue
        result = fn(sd)
        if result:
            per_step[step] = result
    return per_step


def build_radar_data(
    per_step: Dict[str, Dict[str, Dict[str, float]]]
) -> Tuple[List[str], List[Tuple[str, str, str, str]],
          List[str], List[List[Optional[float]]],
          List[List[float]], List[float]]:
    """
    Returns (axis_labels, axis_specs, structures, raw, normalized, scores).
    Each `raw[i][j]` may be None for missing data; `normalized[i][j]` is 0
    when raw is None. Composite score is mean of normalized row.
    """
    available: List[Tuple[str, str, str, str]] = []  # (label, step, key, better)
    for label, priority, key, better in COMBINED_AXES:
        for step in priority:
            sdata = per_step.get(step, {})
            values = [m.get(key) for m in sdata.values() if m.get(key) is not None]
            if len({float(v) for v in values}) >= 2:
                available.append((label, step, key, better))
                break

    structs_set = set()
    for step_data in per_step.values():
        structs_set.update(step_data.keys())
    structs = sorted(structs_set)

    if not available or len(structs) < 2:
        return [], [], structs, [], [], []

    # Raw values (Optional[float])
    raw: List[List[Optional[float]]] = []
    for s in structs:
        row: List[Optional[float]] = []
        for _, step, key, _ in available:
            v = per_step.get(step, {}).get(s, {}).get(key)
            row.append(float(v) if v is not None else None)
        raw.append(row)

    # Normalize each column independently
    n_axes = len(available)
    norm: List[List[float]] = [[0.0] * n_axes for _ in structs]
    for j in range(n_axes):
        col = [raw[i][j] for i in range(len(structs)) if raw[i][j] is not None]
        if not col:
            continue
        mn = min(col)
        mx = max(col)
        rng = mx - mn if mx != mn else 1.0
        better = available[j][3]
        for i in range(len(structs)):
            v = raw[i][j]
            if v is None:
                norm[i][j] = 0.0  # missing = worst (clearly visible as a dip)
            else:
                n = (v - mn) / rng
                if better == "lower":
                    n = 1.0 - n
                norm[i][j] = n

    scores = [sum(row) / len(row) for row in norm]
    labels = [a[0] for a in available]
    return labels, available, structs, raw, norm, scores


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_radar(
    labels: List[str],
    structs: List[str],
    norm: List[List[float]],
    scores: List[float],
    output_path: Path,
    dpi: int,
    max_n: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    order = sorted(range(len(structs)), key=lambda i: scores[i], reverse=True)
    keep = order[: max_n if max_n > 0 else len(order)]
    structs_p = [structs[i] for i in keep]
    norm_p = [norm[i] for i in keep]

    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection="polar"))
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(structs_p), 1)))
    for data, name, color in zip(norm_p, structs_p, colors):
        d = list(data) + [data[0]]
        ax.plot(angles, d, "o-", linewidth=2, label=name, color=color)
        ax.fill(angles, d, alpha=0.15, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, size=9)
    ax.set_ylim(0, 1)
    ax.set_title(
        f"Combined Stability Comparison Across All Steps\n"
        f"({len(structs_p)} structure(s) shown, outer = better)",
        pad=20, fontsize=14, fontweight="bold",
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.0), fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"  Generated: {output_path}")


def render_ranking(
    labels: List[str],
    structs: List[str],
    norm: List[List[float]],
    scores: List[float],
    output_path: Path,
    dpi: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    order = sorted(range(len(structs)), key=lambda i: scores[i], reverse=True)
    sorted_names = [structs[i] for i in order]
    sorted_norm = [norm[i] for i in order]

    fig, ax = plt.subplots(figsize=(14, max(6, 0.4 * len(sorted_names))))
    metric_colors = plt.cm.Set2(np.linspace(0, 1, len(labels)))
    bottom = np.zeros(len(sorted_names))
    y = np.arange(len(sorted_names))
    for j, (label, mc) in enumerate(zip(labels, metric_colors)):
        vals = [sorted_norm[i][j] for i in range(len(sorted_names))]
        ax.barh(y, vals, 0.6, left=bottom, label=label, color=mc, edgecolor="white")
        bottom += vals

    ax.set_yticks(y)
    ax.set_yticklabels(sorted_names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Composite Score (sum of normalized axes; outer = better)")
    ax.set_title("Combined Ranking Across All Steps", fontsize=14, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8, ncol=min(4, len(labels)))
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"  Generated: {output_path}")


# ---------------------------------------------------------------------------
# CSV outputs
# ---------------------------------------------------------------------------

def write_csv_outputs(
    labels: List[str],
    axis_specs: List[Tuple[str, str, str, str]],
    structs: List[str],
    raw: List[List[Optional[float]]],
    norm: List[List[float]],
    scores: List[float],
    per_step: Dict[str, Dict[str, Dict[str, float]]],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    norm_csv = output_dir / "combined_metrics_normalized.csv"
    with open(norm_csv, "w", newline="") as f:
        w = csv.writer(f)
        # Header annotates which step + key fed each axis
        header_axes = [f"{lbl} [{step}.{key}]" for (lbl, step, key, _) in axis_specs]
        w.writerow(["structure", *header_axes, "composite_score"])
        for i, s in enumerate(structs):
            w.writerow([s, *(f"{v:.4f}" for v in norm[i]), f"{scores[i]:.4f}"])

    raw_csv = output_dir / "combined_metrics.csv"
    with open(raw_csv, "w", newline="") as f:
        w = csv.writer(f)
        cols: List[Tuple[str, str]] = []
        for step, sdata in per_step.items():
            keys = set()
            for m in sdata.values():
                keys.update(m.keys())
            for k in sorted(keys):
                cols.append((step, k))
        w.writerow(["structure", *(f"{step}.{k}" for step, k in cols)])
        for s in structs:
            row: List[str] = [s]
            for step, k in cols:
                v = per_step.get(step, {}).get(s, {}).get(k)
                row.append("" if v is None else f"{v:.6f}")
            w.writerow(row)

    print(f"  Generated: {norm_csv}")
    print(f"  Generated: {raw_csv}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_step_prefix_args(items: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            continue
        k, v = item.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def main() -> int:
    p = argparse.ArgumentParser(
        description="Aggregate per-structure metrics across GROMACS pipeline steps."
    )
    p.add_argument("--run-dir", required=True, type=Path,
                   help="Pipeline run directory (contains 1_quick_stability/, 2_compare_chain_stability/, ...)")
    p.add_argument("--output-dir", required=True, type=Path,
                   help="Where to write combined plots / CSVs")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--plot-format", default="png", choices=["png", "svg", "pdf"])
    p.add_argument("--max-structures", type=int, default=12,
                   help="Max structures to display in radar plot (top-N by composite score)")
    p.add_argument("--source-steps", nargs="*", default=None,
                   help="Limit aggregation to these step names (default: auto-detect all)")
    p.add_argument("--step-prefix", action="append", default=[],
                   help='Step-prefix override "step=prefix" (repeatable)')
    args = p.parse_args()

    run_dir: Path = args.run_dir
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    step_prefixes = parse_step_prefix_args(args.step_prefix)
    source_steps = args.source_steps if args.source_steps else None

    print(f"Aggregating combined metrics from: {run_dir}")
    if source_steps:
        print(f"Restricted to source steps: {source_steps}")

    per_step = aggregate(run_dir, step_prefixes, source_steps)
    if not per_step:
        print("No step outputs found - nothing to aggregate.")
        return 0

    print(f"Found data from steps: {sorted(per_step.keys())}")
    for step, sdata in per_step.items():
        print(f"  {step}: {len(sdata)} structure(s)")

    labels, axis_specs, structs, raw, norm, scores = build_radar_data(per_step)
    if not labels:
        print("Insufficient data for combined plots (need >=1 axis with >=2 distinct values).")
        return 0

    print(f"Combined axes ({len(labels)}): {labels}")
    print(f"Structures (union): {len(structs)}")

    fmt = args.plot_format
    radar_path = output_dir / f"combined_radar_all_steps.{fmt}"
    rank_path = output_dir / f"combined_ranking_all_steps.{fmt}"

    render_radar(labels, structs, norm, scores, radar_path,
                 dpi=args.dpi, max_n=args.max_structures)
    render_ranking(labels, structs, norm, scores, rank_path, dpi=args.dpi)
    write_csv_outputs(labels, axis_specs, structs, raw, norm, scores, per_step, output_dir)

    print(f"\nCombined plots: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
