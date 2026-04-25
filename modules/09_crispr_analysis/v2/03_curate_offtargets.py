#!/usr/bin/env python3
"""
Module: 03_curate_offtargets.py
Stage [3] — Paralog-aware off-target curation.

Reads the rescored guide TSV (stage [2]) plus the CRISPOR off-target report,
labels guides that hit DMP-family paralogs, and adds a CFD-sum penalty column.

Usage:
    python3 03_curate_offtargets.py \\
        --input         <rescored.tsv>          \\
        --offtarget-tsv <crispor_offtargets.tsv> \\
        --outdir        <output_dir>             \\
        --paralog-patterns DMP HAP2 GCS1        \\
        --paralog-gene-ids SmelDMP01 SmelDMP02  \\
        --paralog-hit-threshold 1               \\
        --cfd-sum-threshold 0.2
"""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


# ── Paralog pattern compilation ────────────────────────────────────────────
# Characters whose presence implies the user supplied a deliberate regex
# rather than a bare family tag. Bare tags get anchored to prevent substring
# false-positives like "HDMP" or "XDMP" matching the "DMP" family; deliberate
# regex patterns are passed through as-is.
_REGEX_METACHARS = re.compile(r'[\\\[\]\(\)\|\^\$\?\*\+\{\}]')


def _compile_paralog_pattern(pattern: str) -> re.Pattern:
    """Compile a paralog pattern.

    - If `pattern` contains regex metacharacters, treat it as a user-supplied
      regex and compile case-insensitively.
    - Otherwise, treat it as a family tag (e.g. "DMP") and anchor it so that:
        * it matches at start-of-string, after a separator (_ - .), or after
          a lowercase letter (species-prefix convention: "SmelDMP", "AtDMP");
        * it requires a digit, separator, or end-of-string after the tag.
      This accepts "DMP", "DMP1", "SmelDMP01", "AtDMP9" while rejecting
      "HDMP", "XDMP", "TDMP_unrelated" — a strict boundary \\b cannot do
      this because "SmelDMP" has no word-boundary between 'l' and 'D'.
    """
    if _REGEX_METACHARS.search(pattern):
        return re.compile(pattern, re.IGNORECASE)
    # Inline (?i:...) keeps the escaped tag case-insensitive while the
    # lookbehind character class [a-z_\-\.] stays case-sensitive — essential
    # for distinguishing "SmelDMP" (lowercase 'l' before 'D' → match) from
    # "HDMP" (uppercase 'H' before 'D' → reject).
    return re.compile(
        r'(?:^|(?<=[a-z_\-\.]))(?i:' + re.escape(pattern) + r')(?:\d+|[_\-.]|$)'
    )


def parse_args():
    p = argparse.ArgumentParser(description="Paralog-aware off-target curation")
    p.add_argument("--input",                required=True)
    p.add_argument("--offtarget-tsv",        required=True,
                   help="CRISPOR off-target report TSV (one row per off-target hit)")
    p.add_argument("--outdir",               required=True)
    p.add_argument("--paralog-patterns",     nargs="+", default=["DMP", "HAP2", "GCS1"],
                   help="Regex patterns (case-insensitive) to identify paralog genes")
    p.add_argument("--paralog-gene-ids",     nargs="*",  default=[],
                   help="Explicit gene IDs to flag as paralogs")
    p.add_argument("--paralog-hit-threshold",type=int,   default=1,
                   help="Min paralog hits to label guide 'paralog_risk'")
    p.add_argument("--cfd-sum-threshold",    type=float, default=0.2,
                   help="CFD off-target sum above which guide is flagged 'high_offtarget'")
    p.add_argument("--overwrite",            action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


def load_offtargets(path: str) -> dict[str, list[dict]]:
    """Return {guide_id: [hit_rows]} from a CRISPOR off-target TSV."""
    hits: dict[str, list[dict]] = {}
    path_obj = Path(path)
    if not path_obj.exists():
        # INFO not WARN: an absent off-target file is expected in the v3
        # plant-only pipeline when Cas-OFFinder has not yet been run.
        # Guides proceed with empty off-target annotation (paralog_hits=0).
        _log(f"[03_curate] Off-target TSV not found: {path} — skipping curation.", level="INFO")
        return hits

    with open(path_obj, newline="") as fh:
        # CRISPOR off-target headers: guideId, offtargetSeq, mismatchCount, cfdScore, gene, ...
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = row.get("guideId") or row.get("guide_id") or row.get("name") or ""
            hits.setdefault(gid, []).append(row)
    return hits


def get_guide_id(row: dict) -> str:
    for col in ("guideId", "guide_id", "name", "ID", "id"):
        if col in row and row[col]:
            return row[col].strip()
    return ""


def cfd_sum(hits: list[dict]) -> float:
    total = 0.0
    for h in hits:
        for col in ("cfdScore", "CFDScore", "cfd_score"):
            if col in h:
                try:
                    total += float(h[col])
                    break
                except ValueError:
                    pass
    return total


def main():
    args = parse_args()

    inpath  = Path(args.input)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.curated.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[03_curate] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    # Compile paralog regex patterns. Bare family tags (e.g. "DMP") are
    # auto-anchored to avoid "HDMP"/"XDMP" substring false-positives while
    # still matching "SmelDMP01", "AtDMP9", "DMP1" etc. See
    # _compile_paralog_pattern docstring for boundary semantics.
    paralog_re = [_compile_paralog_pattern(p) for p in args.paralog_patterns]
    paralog_ids = set(args.paralog_gene_ids)

    # Load guide table
    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    # Load off-target data
    offtargets = load_offtargets(args.offtarget_tsv)

    new_fields = fields + ["paralog_hits", "paralog_risk", "cfd_sum", "offtarget_flag"]

    for row in rows:
        gid  = get_guide_id(row)
        hits = offtargets.get(gid, [])

        # Count paralog hits
        ph_count = 0
        for h in hits:
            gene_name = h.get("gene") or h.get("geneName") or h.get("gene_name") or ""
            is_paralog = (
                gene_name in paralog_ids
                or any(rx.search(gene_name) for rx in paralog_re)
            )
            if is_paralog:
                ph_count += 1

        cfd = cfd_sum(hits)
        row["paralog_hits"]   = str(ph_count)
        row["paralog_risk"]   = "paralog_risk" if ph_count >= args.paralog_hit_threshold else ""
        row["cfd_sum"]        = f"{cfd:.4f}"
        row["offtarget_flag"] = "high_offtarget" if cfd > args.cfd_sum_threshold else ""

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    paralog_n  = sum(1 for r in rows if r.get("paralog_risk"))
    offtarget_n= sum(1 for r in rows if r.get("offtarget_flag"))
    _log(f"[03_curate] {len(rows)} guides curated; "
          f"{paralog_n} paralog_risk, {offtarget_n} high_offtarget -> {outpath}", level="INFO")


if __name__ == "__main__":
    main()
