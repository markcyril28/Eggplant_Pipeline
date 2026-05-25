#!/usr/bin/env python3
"""
Translate a CDS FASTA into a protein FASTA (standard genetic code, table 1).

Stop codons are translated to '*' and stripped from the trailing position.
Internal stops are kept as '*' so HMMER's hmmsearch can still process them
without crashing (it ignores '*' inside sequences). Sequences whose length
is not a multiple of 3 are trimmed at the 3' end with a warning.

Usage:
    python3 translate_cds_to_protein.py <input.cds.fa> <output.pep.fa>
    python3 translate_cds_to_protein.py <input.cds.fa> <output.pep.fa> --force

Designed to be idempotent: skips work if the output is newer than the input
unless --force is passed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    # Drop everything that is not a DNA base before chunking.
    seq = "".join(c for c in seq if c in "ACGTN")
    aa = []
    n = len(seq) - (len(seq) % 3)
    for i in range(0, n, 3):
        codon = seq[i:i + 3]
        if "N" in codon:
            aa.append("X")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa).rstrip("*")


def iter_fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(out_path: Path, records, line_width: int = 80) -> int:
    n = 0
    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        for header, aa in records:
            if not aa:
                continue
            out.write(f">{header}\n")
            for i in range(0, len(aa), line_width):
                out.write(aa[i:i + line_width] + "\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path, help="Input CDS FASTA")
    ap.add_argument("output", type=Path, help="Output protein FASTA")
    ap.add_argument("--force", action="store_true",
                    help="Re-translate even if output is newer than input")
    args = ap.parse_args()

    if not args.input.is_file():
        print(f"ERROR: input not found: {args.input}", file=sys.stderr)
        return 2

    if args.output.exists() and not args.force:
        if args.output.stat().st_mtime >= args.input.stat().st_mtime:
            print(f"translate_cds_to_protein: up-to-date, skipping ({args.output.name})",
                  file=sys.stderr)
            return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    records = ((h, translate(s)) for h, s in iter_fasta(args.input))
    n = write_fasta(args.output, records)
    print(f"translate_cds_to_protein: wrote {n} records -> {args.output}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
