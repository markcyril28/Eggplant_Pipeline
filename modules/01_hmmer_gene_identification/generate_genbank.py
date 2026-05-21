#!/usr/bin/env python3
"""
generate_genbank.py

Walk an e_GENE_Structures/ directory and produce structure.gb for every
mRNA folder that has structure.gff3 + structure.fa. Safe to re-run
idempotently; skips folders missing inputs and reports a summary.

Supports the with/without flank layout written by EXTRACT_GENE_STRUCTURES:
    {GENE}/without_up_and_downstream/{mRNA}/   -> flank forced to 0
    {GENE}/with_up_and_downstream/{mRNA}/      -> flank from --flank-bp
A legacy flat layout {GENE}/{mRNA}/ is also accepted for backward
compatibility (treated as the with-flanks side).

Reuses build_genbank() from extract_gene_structures.py so output is bit-for-bit
identical to the inline GenBank pass performed by EXTRACT_GENE_STRUCTURES.

Flanking sequence:
  With --flank-bp > 0 AND --genome-fasta pointing at a FASTA that has a
  sibling .fai, the script re-fetches the gene span plus flanks directly
  from the genome (recommended — matches EXTRACT_GENE_STRUCTURES output).
  Without --genome-fasta, flanks are forced to 0 and the structure.fa
  >*_gene record is used as the GenBank backbone.
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_gene_structures import (  # noqa: E402
    build_genbank,
    fetch_flanked_gene,
    load_fai,
)


def read_gene_fasta(path, mrna_id):
    """Return the '>{mrna_id}_gene' sequence from a structure.fa, or ''."""
    target = f"{mrna_id}_gene"
    current = None
    buf = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            if s.startswith(">"):
                if current == target:
                    return "".join(buf)
                current = s[1:].split()[0]
                buf = []
            elif s:
                buf.append(s)
    return "".join(buf) if current == target else ""


def partition_gff(path):
    """Return (gene_line, {feat_type: [raw_lines]}) from a structure.gff3."""
    by_type = defaultdict(list)
    gene_line = ""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2].lower()
            if feat == "gene":
                gene_line = line.rstrip("\n")
            elif feat in ("exon", "cds", "five_prime_utr",
                          "three_prime_utr", "intron"):
                by_type[feat].append(line.rstrip("\n"))
    return gene_line, by_type


def process_folder(mrna_dir, organism, overwrite, fasta_fh, fai, flank_bp):
    gff = mrna_dir / "structure.gff3"
    fa = mrna_dir / "structure.fa"
    gb = mrna_dir / "structure.gb"
    if not gff.is_file() or not fa.is_file():
        return "skip-missing-inputs"
    if gb.exists() and not overwrite:
        return "skip-exists"
    mrna_id = mrna_dir.name
    gene_line, by_type = partition_gff(gff)
    if not gene_line:
        return "skip-no-gene-feature"

    up = dn = 0
    if fasta_fh is not None and fai and flank_bp >= 0:
        backbone, up, dn = fetch_flanked_gene(fasta_fh, fai, gene_line, flank_bp)
        if not backbone:
            return "skip-fetch-empty"
    else:
        backbone = read_gene_fasta(fa, mrna_id)
        if not backbone:
            return "skip-no-gene-seq"

    content = build_genbank(
        mrna_id, gene_line,
        by_type["exon"], by_type["cds"],
        by_type["five_prime_utr"], by_type["three_prime_utr"],
        by_type["intron"], backbone,
        up_flank=up, dn_flank=dn,
        organism=organism,
    )
    if not content:
        return "skip-empty-build"
    gb.write_text(content, encoding="ascii")
    # Sibling copy at the gene-folder level for bulk GenBank import
    (mrna_dir.parent / f"{mrna_dir.name}.gb").write_text(content, encoding="ascii")
    return "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True,
                    help="e_GENE_Structures/ directory containing {GENE}/{mRNA}/ subfolders")
    ap.add_argument("--overwrite", default="true")
    ap.add_argument("--organism", default="Solanum melongena",
                    help="Organism name written into SOURCE/ORGANISM/source.organism")
    ap.add_argument("--genome-fasta", default="",
                    help="Genomic FASTA (+ sibling .fai). Required for flank extraction.")
    ap.add_argument("--flank-bp", type=int, default=1000,
                    help="5'/3' flank bp to include in structure.gb (default 1000). "
                         "Requires --genome-fasta; otherwise treated as 0.")
    args = ap.parse_args()

    overwrite = args.overwrite.lower() not in ("false", "0", "no")
    root = Path(args.root)
    if not root.is_dir():
        print(f"ERROR: --root not found: {root}", file=sys.stderr)
        sys.exit(1)

    fasta_fh = None
    fai = None
    effective_flank = 0
    if args.genome_fasta:
        if not os.path.isfile(args.genome_fasta):
            print(f"WARNING: --genome-fasta not found: {args.genome_fasta}", file=sys.stderr)
        else:
            fai_path = args.genome_fasta + ".fai"
            if not os.path.isfile(fai_path):
                print(f"WARNING: .fai index missing for {args.genome_fasta}", file=sys.stderr)
            else:
                fai = load_fai(fai_path)
                fasta_fh = open(args.genome_fasta, "rb")
                effective_flank = max(0, args.flank_bp)

    try:
        counts = defaultdict(int)
        for gene_dir in sorted(root.iterdir()):
            if not gene_dir.is_dir():
                continue

            # Preferred layout: gene_dir / {with,without}_up_and_downstream / mRNA
            no_flank_root = gene_dir / "without_up_and_downstream"
            with_flank_root = gene_dir / "with_up_and_downstream"
            split_layout = no_flank_root.is_dir() or with_flank_root.is_dir()

            if split_layout:
                for sub_root, sub_flank in (
                    (no_flank_root, 0),
                    (with_flank_root, effective_flank),
                ):
                    if not sub_root.is_dir():
                        continue
                    for mrna_dir in sorted(sub_root.iterdir()):
                        if not mrna_dir.is_dir():
                            continue
                        outcome = process_folder(mrna_dir, args.organism, overwrite,
                                                 fasta_fh, fai, sub_flank)
                        counts[outcome] += 1
                        if outcome.startswith("skip-") and outcome != "skip-exists":
                            print(f"  {outcome}: {mrna_dir}", file=sys.stderr)
            else:
                # Legacy flat layout: gene_dir / mRNA
                for mrna_dir in sorted(gene_dir.iterdir()):
                    if not mrna_dir.is_dir():
                        continue
                    outcome = process_folder(mrna_dir, args.organism, overwrite,
                                             fasta_fh, fai, effective_flank)
                    counts[outcome] += 1
                    if outcome.startswith("skip-") and outcome != "skip-exists":
                        print(f"  {outcome}: {mrna_dir}", file=sys.stderr)

        ok = counts.get("ok", 0)
        skipped = sum(v for k, v in counts.items() if k.startswith("skip-"))
        flank_note = f" flank={effective_flank}bp" if effective_flank else ""
        print(f"GenBank: wrote={ok} skipped={skipped}{flank_note} @ {root}",
              file=sys.stderr)
    finally:
        if fasta_fh:
            fasta_fh.close()


if __name__ == "__main__":
    main()
