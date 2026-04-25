#!/usr/bin/env python3
"""
Extract DMP query FASTA records from a local genome's CDS/protein file by
matching regex patterns against the FASTA headers.

Usage:
    extract_dmp_from_local.py --fasta <path> --patterns <p1>[,<p2>...] \
                              --out <path> [--name <gene_name>]

Returns: number of records written.
"""

import argparse
import re
import sys
from pathlib import Path


def extract(fasta_path: Path, patterns: list[str], out_path: Path, name: str | None) -> int:
    compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    written = 0
    with fasta_path.open() as fin, out_path.open("a") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip()
                keep = any(p.search(header) for p in compiled)
                if keep:
                    if name:
                        fout.write(f">{name} | {header}\n")
                    else:
                        fout.write(line)
                    written += 1
            elif keep:
                fout.write(line)
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="Local CDS or protein FASTA")
    ap.add_argument("--patterns", required=True,
                    help="Pipe-separated regex patterns (e.g. ^AT1G09157|^AT5G39650)")
    ap.add_argument("--out", required=True, help="Output FASTA path (appended)")
    ap.add_argument("--name", default=None,
                    help="Optional gene name to prepend to extracted headers")
    args = ap.parse_args()

    fasta = Path(args.fasta)
    if not fasta.exists():
        print(f"[extract] source missing: {fasta}", file=sys.stderr)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    patterns = [p for p in args.patterns.split("|") if p]
    n = extract(fasta, patterns, out, args.name)
    print(f"[extract] {fasta.name}: {n} record(s) -> {out.name}")
    return 0 if n > 0 else 2


if __name__ == "__main__":
    sys.exit(main())
