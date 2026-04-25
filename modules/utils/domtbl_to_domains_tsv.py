#!/usr/bin/env python3
"""
Convert a HMMER --domtblout file into the flat per-gene domain TSV expected
by modules/09_crispr_analysis/v2/06_protein_consequence.py.

Output columns (tab-separated, header row):
    gene_id    domain_name    aa_start    aa_end

The domtbl target names usually include a transcript suffix (e.g.
SMEL5_01g008730.1). The CRISPR stage keys domain lookups by the gene_id
emitted in stage [05] output (SMEL5_01g008730, no suffix), so this script
strips the trailing ".<digits>" isoform suffix by default.

Usage:
    python3 modules/utils/domtbl_to_domains_tsv.py \\
        --domtbl III_RESULT/DMP/01_Identification/.../PF05078_DMP_hits_filtered.domtbl \\
        --output II_INPUTS/DMP/Pfam/DMP_domains.tsv
"""

import argparse
import re
from pathlib import Path


_GENE_SUFFIX_RE = re.compile(r"\.\d+$")


def parse_args():
    p = argparse.ArgumentParser(description="HMMER domtbl -> DMP_domains.tsv")
    p.add_argument("--domtbl", required=True, type=Path,
                   help="HMMER --domtblout file from hmmsearch/hmmscan.")
    p.add_argument("--output", required=True, type=Path,
                   help="Destination TSV (columns: gene_id domain_name aa_start aa_end).")
    p.add_argument("--keep-isoform-suffix", action="store_true",
                   help="Preserve the .N transcript suffix on target names "
                        "(default: strip so IDs match stage [05] gene_stem).")
    p.add_argument("--coord-column", choices=("ali", "env"), default="ali",
                   help="Which coord pair to emit (default: ali; HMMER domtbl "
                        "cols 18-19). Use 'env' for the more permissive "
                        "envelope coords (cols 20-21).")
    return p.parse_args()


def iter_domtbl_hits(domtbl_path: Path):
    """Yield (target_id, domain_name, ali_from, ali_to, env_from, env_to)
    tuples for every data row in a HMMER --domtblout file."""
    with open(domtbl_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            # domtbl is whitespace-delimited with 22+ columns; split on any run
            # of spaces since target IDs never contain whitespace.
            cols = line.split()
            if len(cols) < 21:
                continue
            target_id   = cols[0]
            domain_name = cols[3]
            try:
                ali_from = int(cols[17])
                ali_to   = int(cols[18])
                env_from = int(cols[19])
                env_to   = int(cols[20])
            except ValueError:
                continue
            yield target_id, domain_name, ali_from, ali_to, env_from, env_to


def main():
    args = parse_args()

    if not args.domtbl.exists():
        raise SystemExit(f"domtbl not found: {args.domtbl}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    coord_idx = (2, 3) if args.coord_column == "ali" else (4, 5)
    n_written = 0
    seen: set[tuple[str, str, int, int]] = set()

    with open(args.output, "w", newline="") as out:
        out.write("gene_id\tdomain_name\taa_start\taa_end\n")
        for hit in iter_domtbl_hits(args.domtbl):
            target_id = hit[0]
            gene_id = target_id if args.keep_isoform_suffix \
                else _GENE_SUFFIX_RE.sub("", target_id)
            coord_from = hit[coord_idx[0]]
            coord_to   = hit[coord_idx[1]]
            key = (gene_id, hit[1], coord_from, coord_to)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"{gene_id}\t{hit[1]}\t{coord_from}\t{coord_to}\n")
            n_written += 1

    print(f"Wrote {n_written} domain row(s) to {args.output}")


if __name__ == "__main__":
    main()
