#!/usr/bin/env python3
"""
Module: 01b_plant_sgrna_filter.py
Stage [1b] — Plant-specific sgRNA pre-filter.

Runs AFTER stage [1] CRISPOR design + filtering and BEFORE stage [2] on-target
rescoring. Enforces the three biology-mandated rules that mammalian-trained
scorers do not encode:

    (P1) Poly-T Pol III termination signal (TTTT+)
         Plant U6/U3 promoters use RNA Pol III, which terminates at any run
         of ≥ termination_run_length T residues. A guide with an internal
         TTTT yields a truncated sgRNA and no editing.
         Default action: REJECT.

    (P2) 5' nucleotide for Pol III promoter
         U6 promoter initiates only at G (some vectors auto-prepend a G);
         U3 prefers A. Misaligned 5' nucleotides drop expression severely.
         Default action: FLAG (annotate; user picks vector later).

    (P3) GC content window 30–70 %
         Extreme GC reduces Cas9 loading and strand invasion. Flagged,
         not rejected — some regulatory regions require extreme-GC guides
         and the biologist may want to keep them for follow-up.

Input : stage-[1] filtered TSV (*.filtered.tsv)
Output: plant-filtered TSV  (*.plant_filtered.tsv) at --outdir

Added columns:
    protospacer                : 20-nt guide extracted from the CRISPOR row
    plant_pass                 : True/False (False = rejected by stage 1b)
    plant_reject_reason        : "" when pass, else pipe-joined reasons
    poly_t_run                 : longest T run in the protospacer
    poly_t_flag                : "" | "polyT_terminator"
    gc_fraction                : 0–1
    gc_flag                    : "" | "low_gc" | "high_gc"
    five_prime_nt              : first nt of the protospacer (A/C/G/T)
    promoter_compat            : "" | "U6_ready" | "U3_ready" | "prepend_G" …
    prepended_protospacer      : 21-nt sequence when action="prepend"

Usage:
    python3 01b_plant_sgrna_filter.py \\
        --input   <gene>.filtered.tsv \\
        --outdir  <01b_Plant_Filter/> \\
        --termination-run-length 4 \\
        --termination-action     reject \\
        --promoter-type          U6 \\
        --promoter-action        flag \\
        --gc-min 0.30 --gc-max 0.70

References (plant CRISPR):
    Ma X et al. 2015. Mol Plant 8(8):1274-84.  doi:10.1016/j.molp.2015.04.007
    Gao Y, Zhao Y. 2014. J Integr Plant Biol 56(4):343-9.
        doi:10.1111/jipb.12152
    Nekrasov V et al. 2013. Nat Biotechnol 31:691-3. doi:10.1038/nbt.2655
    Shan Q et al. 2014. Nat Protoc 9(10):2395-410.   doi:10.1038/nprot.2014.157
    Liu X et al. 2016. Nucleic Acids Res 44(10):gkw223.
        doi:10.1093/nar/gkw223
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


# ─── Sequence helpers ────────────────────────────────────────────────────────

def _extract_protospacer(row: dict) -> str:
    """Return the 20-nt protospacer from a CRISPOR row."""
    for col in ("targetSeq", "guideSeq", "protospacer", "sequence", "Sequence"):
        if col in row and row[col]:
            seq = row[col].strip().upper().replace("U", "T")
            # CRISPOR often emits 23 nt (protospacer + PAM); strip to 20 from 5' end
            # UNLESS the column explicitly names the PAM-trimmed field.
            if col in ("protospacer", "guideSeq"):
                return seq[:20]
            # targetSeq/sequence variants typically put the protospacer first
            return seq[:20]
    return ""


def _longest_T_run(seq: str) -> int:
    """Length of the longest consecutive T stretch in seq."""
    if not seq:
        return 0
    best = cur = 0
    for ch in seq:
        if ch == "T":
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in "GC")
    return gc / len(seq)


# ─── Per-row plant-filter evaluation ─────────────────────────────────────────

def evaluate(row: dict,
             term_run: int, term_action: str,
             promoter_type: str, promoter_action: str,
             gc_min: float, gc_max: float) -> dict:
    """Return a dict of the new columns; does not mutate `row`."""
    protospacer = _extract_protospacer(row)
    out = {
        "protospacer":           protospacer,
        "plant_pass":            "True",
        "plant_reject_reason":   "",
        "poly_t_run":            "0",
        "poly_t_flag":           "",
        "gc_fraction":           "0.0000",
        "gc_flag":               "",
        "five_prime_nt":         "",
        "promoter_compat":       "",
        "prepended_protospacer": "",
    }
    if not protospacer or len(protospacer) < 20:
        out["plant_pass"] = "False"
        out["plant_reject_reason"] = "no_protospacer"
        return out

    reject_reasons = []

    # ── (P1) Poly-T terminator ────────────────────────────────────────────
    longest_t = _longest_T_run(protospacer)
    out["poly_t_run"] = str(longest_t)
    if longest_t >= term_run:
        out["poly_t_flag"] = "polyT_terminator"
        if term_action == "reject":
            reject_reasons.append(f"polyT_T{longest_t}")

    # ── (P3) GC content window ────────────────────────────────────────────
    gc = _gc_fraction(protospacer)
    out["gc_fraction"] = f"{gc:.4f}"
    if gc < gc_min:
        out["gc_flag"] = "low_gc"
    elif gc > gc_max:
        out["gc_flag"] = "high_gc"
    # GC is never a hard reject — only a flag (per config convention).

    # ── (P2) 5' nucleotide for Pol III promoter ──────────────────────────
    first_nt = protospacer[0]
    out["five_prime_nt"] = first_nt

    def _expected_5p(promoter: str) -> str:
        return {"U6": "G", "U3": "A"}.get(promoter.upper(), "")

    expected = _expected_5p(promoter_type)
    if expected:
        if first_nt == expected:
            out["promoter_compat"] = f"{promoter_type.upper()}_ready"
        else:
            if promoter_action == "require":
                reject_reasons.append(f"needs_5p_{expected}_for_{promoter_type.upper()}")
                out["promoter_compat"] = f"incompatible_{promoter_type.upper()}"
            elif promoter_action == "prepend":
                out["promoter_compat"] = f"prepend_{expected}_for_{promoter_type.upper()}"
                out["prepended_protospacer"] = expected + protospacer
            else:  # "flag"
                out["promoter_compat"] = f"needs_5p_{expected}"
    # promoter_type == "any" — no compatibility check

    if reject_reasons:
        out["plant_pass"] = "False"
        out["plant_reject_reason"] = "|".join(reject_reasons)

    return out


# ─── CLI ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Plant sgRNA pre-filter (U6/U3, TTTT, GC)")
    p.add_argument("--input",                      required=True,
                   help="Stage-[1] filtered TSV (*.filtered.tsv)")
    p.add_argument("--outdir",                     required=True)
    p.add_argument("--termination-run-length",     type=int,   default=4,
                   help="Minimum consecutive T count that triggers Pol III termination (default 4).")
    p.add_argument("--termination-action",         default="reject",
                   choices=("reject", "flag"),
                   help="'reject' drops guides with polyT runs; 'flag' annotates only.")
    p.add_argument("--promoter-type",              default="U6",
                   choices=("U6", "U3", "any"),
                   help="Pol III promoter the guide will be expressed from.")
    p.add_argument("--promoter-action",            default="flag",
                   choices=("require", "flag", "prepend"),
                   help="'require' drops incompatible 5' nt; 'flag' annotates; "
                        "'prepend' emits a candidate prepended sequence.")
    p.add_argument("--gc-min",                     type=float, default=0.30)
    p.add_argument("--gc-max",                     type=float, default=0.70)
    p.add_argument("--enabled",                    default="true",
                   help="'false'/'0'/'no' makes this stage a no-op pass-through.")
    p.add_argument("--overwrite",                  action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


def main():
    args = parse_args()

    inpath = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.plant_filtered.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[01b_plant] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    # CRISPR-P v2.0 exports comma-delimited .csv files; pipeline-internal
    # files use tab-delimited .tsv. Detect from extension so both work.
    delimiter = "," if inpath.suffix.lower() == ".csv" else "\t"
    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    new_cols = [
        "protospacer", "plant_pass", "plant_reject_reason",
        "poly_t_run", "poly_t_flag",
        "gc_fraction", "gc_flag",
        "five_prime_nt", "promoter_compat", "prepended_protospacer",
    ]
    new_fields = fields + [c for c in new_cols if c not in fields]

    enabled = str(args.enabled).lower() not in ("false", "0", "no", "off")

    if not enabled:
        # Pass-through: copy rows verbatim, still add the new columns as "" so
        # downstream stages see a consistent schema whether the filter ran or not.
        _log(f"[01b_plant] DISABLED — pass-through copy of {len(rows)} rows.", level="WARN")
        for row in rows:
            for c in new_cols:
                row.setdefault(c, "")
            row["plant_pass"] = "True"
        with open(outpath, "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t",
                           extrasaction="ignore").writeheader()
        with open(outpath, "a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t",
                           extrasaction="ignore").writerows(rows)
        return

    # Normal path: score each row and optionally drop rejected ones.
    kept: list[dict] = []
    rejected_poly_t = 0
    rejected_promoter = 0
    rejected_no_protospacer = 0
    flagged_gc = 0
    flagged_promoter = 0
    for row in rows:
        ann = evaluate(
            row,
            term_run=args.termination_run_length,
            term_action=args.termination_action,
            promoter_type=args.promoter_type,
            promoter_action=args.promoter_action,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
        )
        row.update(ann)
        if ann["plant_pass"] == "False":
            if "no_protospacer" in ann["plant_reject_reason"]:
                rejected_no_protospacer += 1
            elif "polyT" in ann["plant_reject_reason"]:
                rejected_poly_t += 1
            elif "needs_5p" in ann["plant_reject_reason"]:
                rejected_promoter += 1
            # Rejected guides are NOT written — they cannot be expressed in a
            # plant vector, so carrying them forward wastes stage [2]-[8] compute.
            continue
        if ann["gc_flag"]:
            flagged_gc += 1
        if ann["promoter_compat"] and "ready" not in ann["promoter_compat"]:
            flagged_promoter += 1
        kept.append(row)

    if rejected_no_protospacer:
        _log(
            f"[01b_plant] {rejected_no_protospacer}/{len(rows)} rows had no recognisable "
            f"protospacer column — check that the input file uses a known column name "
            f"(targetSeq, guideSeq, protospacer, sequence, Sequence) and the correct "
            f"delimiter (tab for .tsv, comma for .csv).",
            level="ERROR",
        )

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(kept)

    _log(
        f"[01b_plant] {len(rows)} in; {len(kept)} kept; "
        f"rejected: polyT={rejected_poly_t}, 5'-{args.promoter_type}={rejected_promoter}, "
        f"no_protospacer={rejected_no_protospacer}; "
        f"flagged: GC={flagged_gc}, promoter={flagged_promoter} -> {outpath}",
        level="INFO",
    )


if __name__ == "__main__":
    main()
