#!/usr/bin/env python3
"""
generate_variants.py -- build a single truncation-variant FASTA from a
source sequence by applying one or more deletion ranges.

Called by 14_interaction_Domain_Mapping.sh for each row in [variants] and
[dmp_variants].  All deletion coordinates refer to the ORIGINAL sequence
(1-indexed, inclusive).  When multiple ranges are given they are applied
simultaneously (positions to remove are unioned before slicing), so the
ranges in the TOML can be specified in any order and must not overlap.
"""

import argparse
import sys


def parse_deletions(deletion_str):
    """Return list of (start, end) tuples (1-indexed, inclusive)."""
    if not deletion_str.strip():
        return []
    ranges = []
    for part in deletion_str.split(","):
        part = part.strip()
        if not part:
            continue
        start_s, end_s = part.split("-", 1)
        ranges.append((int(start_s), int(end_s)))
    return ranges


def apply_deletions(seq, ranges):
    """Remove all positions covered by any range; coordinates are original."""
    if not ranges:
        return seq
    n = len(seq)
    to_remove = set()
    for start, end in ranges:
        if start < 1 or end > n or start > end:
            raise ValueError(
                f"Deletion range {start}-{end} out of bounds for "
                f"sequence of length {n}."
            )
        to_remove.update(range(start - 1, end))
    return "".join(c for i, c in enumerate(seq) if i not in to_remove)


def read_first_fasta(path):
    """Return (header_without_gt, sequence) for the first FASTA record."""
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            else:
                parts.append(line)
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    return header, "".join(parts)


def write_fasta(path, header, seq, line_width=80):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            fh.write(seq[i : i + line_width] + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Build a truncation-variant FASTA from a source sequence."
    )
    ap.add_argument("--input", required=True,
                    help="Source FASTA (first record used)")
    ap.add_argument("--name", required=True,
                    help="Variant identifier written into the FASTA header")
    ap.add_argument("--description", default="",
                    help="Human-readable variant description for the header")
    ap.add_argument("--deletions", default="",
                    help='Comma-separated "start-end" ranges to delete '
                         "(1-indexed, inclusive); empty = WT")
    ap.add_argument("--output", required=True, help="Output FASTA path")
    args = ap.parse_args()

    orig_header, orig_seq = read_first_fasta(args.input)
    ranges = parse_deletions(args.deletions)
    variant_seq = apply_deletions(orig_seq, ranges)

    header_parts = [args.name]
    if args.description:
        header_parts.append(args.description)
    header_parts.append(f"len={len(variant_seq)}")
    if ranges:
        header_parts.append(
            "deleted=" + ",".join(f"{s}-{e}" for s, e in ranges)
        )
    write_fasta(args.output, " | ".join(header_parts), variant_seq)

    suffix = (
        f" (deleted {args.deletions})" if ranges else " (WT, no deletion)"
    )
    print(
        f"  {args.name}: {len(orig_seq)} aa -> {len(variant_seq)} aa{suffix}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
