#!/usr/bin/env python3
"""
extract_gene_structures.py

For each HMMER hit ID, extract the corresponding gene structure from a GFF/GFF3
annotation file and write per-feature subfiles into e_GENE_Structures/.

Output layout (relative to --output):
    {gene_name}/{mrna_id}/
        structure.gff3        all features for this gene locus (+ derived introns)
        exons.gff3
        cds.gff3
        five_prime_utr.gff3
        three_prime_utr.gff3
        introns.gff3          derived from gaps between consecutive exons
        structure.fa          per-feature DNA sequences (only with --genome-fasta)
        structure.gb          Benchling-ready GenBank with local coords
                              (only with --genome-fasta)

structure.fa layout (one record per feature, blank line between records):
        >{mrna_id}_gene                full genomic span
        >{mrna_id}_mRNA                spliced (exons concatenated in txpt order)
        >{mrna_id}_five_prime_UTR      concatenated 5' UTR pieces
        >{mrna_id}_three_prime_UTR     concatenated 3' UTR pieces
        >{mrna_id}_exon_N              individual exons (transcript order)
        >{mrna_id}_CDS_N               individual CDS segments (transcript order)
        >{mrna_id}_intron_N            introns derived from exon gaps
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def _parse_attrs(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def load_gff(gff_path):
    """
    Parse GFF2/GFF3. Returns:
      gene_records, mrna_records, child_records, mrna_to_gene, gene_to_mrnas
    """
    gene_records = {}
    mrna_records = {}
    child_records = defaultdict(list)
    mrna_to_gene = {}
    gene_to_mrnas = defaultdict(list)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2].lower()
            attrs = _parse_attrs(parts[8])
            feat_id = attrs.get("ID", "")
            parent = attrs.get("Parent", "")

            if feat_type == "gene":
                if feat_id:
                    gene_records[feat_id] = line.rstrip("\n")
            elif feat_type in ("mrna", "transcript"):
                if feat_id:
                    mrna_records[feat_id] = line.rstrip("\n")
                    if parent:
                        mrna_to_gene[feat_id] = parent
                        gene_to_mrnas[parent].append(feat_id)
            elif feat_type in ("exon", "cds", "five_prime_utr", "three_prime_utr"):
                if parent:
                    child_records[parent].append(line.rstrip("\n"))

    return gene_records, mrna_records, child_records, mrna_to_gene, gene_to_mrnas


def _derive_introns(exon_lines, anchor_line):
    """Return intron GFF lines derived from gaps between sorted exons."""
    if len(exon_lines) < 2 or not anchor_line:
        return []
    p = anchor_line.split("\t")
    seqid, source, strand = p[0], p[1], p[6]
    intervals = []
    for ln in exon_lines:
        ep = ln.split("\t")
        try:
            intervals.append((int(ep[3]), int(ep[4])))
        except (ValueError, IndexError):
            pass
    intervals.sort()
    introns = []
    for i in range(len(intervals) - 1):
        istart = intervals[i][1] + 1
        iend = intervals[i + 1][0] - 1
        if iend >= istart:
            introns.append(
                f"{seqid}\t{source}\tintron\t{istart}\t{iend}\t.\t{strand}\t.\t"
                f"ID=derived_intron_{i + 1};derived=true"
            )
    return introns


def _resolve_hits(hit_ids_file, mrna_records, gene_records, gene_to_mrnas):
    """Map HMMER hit IDs to mRNA IDs in the GFF. Returns list of (hit, [mrna_ids])."""
    with open(hit_ids_file) as fh:
        raw_ids = [ln.strip() for ln in fh if ln.strip()]

    results = []
    for hit in raw_ids:
        matched = []
        if hit in mrna_records:
            matched = [hit]
        elif hit in gene_records:
            matched = gene_to_mrnas.get(hit, [])
        else:
            stripped = re.sub(r"\.p\d*$|\.pep$", "", hit)
            if stripped != hit:
                if stripped in mrna_records:
                    matched = [stripped]
                elif stripped in gene_records:
                    matched = gene_to_mrnas.get(stripped, [])
            if not matched:
                gene_cand = re.sub(r"\.\d+$", "", hit)
                if gene_cand != hit:
                    if gene_cand in gene_records:
                        matched = gene_to_mrnas.get(gene_cand, [])
                    elif gene_cand in mrna_records:
                        matched = [gene_cand]
        if matched:
            results.append((hit, matched))
    return results


# ---------------------------------------------------------------------------
# FAI-based random-access FASTA reader
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                            "TGCANtgcanYRSWMKVHDBvhdbryswmk")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def load_fai(fai_path):
    """Parse .fai index. Returns {seqid: (length, offset, linebases, linewidth)}."""
    fai = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                fai[parts[0]] = (int(parts[1]), int(parts[2]),
                                 int(parts[3]), int(parts[4]))
    return fai


def fetch_seq(fh, fai, seqid, start, end):
    """1-based inclusive fetch. Returns uppercase ASCII string ('' on miss)."""
    if seqid not in fai:
        return ""
    length, offset, linebases, linewidth = fai[seqid]
    if start < 1:
        start = 1
    if end > length:
        end = length
    if end < start:
        return ""
    start0 = start - 1
    end0 = end  # exclusive in 0-based
    start_line = start0 // linebases
    start_col = start0 % linebases
    end_line = (end0 - 1) // linebases
    end_col = (end0 - 1) % linebases
    file_start = offset + start_line * linewidth + start_col
    file_end = offset + end_line * linewidth + end_col + 1
    fh.seek(file_start)
    raw = fh.read(file_end - file_start).decode("ascii", errors="replace")
    return raw.replace("\n", "").replace("\r", "").upper()


# ---------------------------------------------------------------------------
# Sequence assembly
# ---------------------------------------------------------------------------

def _line_coords(line):
    p = line.split("\t")
    return p[0], int(p[3]), int(p[4]), p[6]


def _concat_seq(lines, fh, fai, strand):
    """Fetch + concatenate in transcript order (reverse for - strand)."""
    if not lines:
        return ""
    sorted_lines = sorted(lines, key=lambda l: int(l.split("\t")[3]))
    pieces = []
    for ln in sorted_lines:
        seqid, s, e, _ = _line_coords(ln)
        pieces.append(fetch_seq(fh, fai, seqid, s, e))
    full = "".join(pieces)
    return revcomp(full) if strand == "-" else full


def _individual_seqs(lines, fh, fai, strand, label):
    """Return list of (header_suffix, seq) in transcript order."""
    if not lines:
        return []
    sorted_lines = sorted(lines, key=lambda l: int(l.split("\t")[3]))
    if strand == "-":
        sorted_lines = list(reversed(sorted_lines))
    out = []
    for i, ln in enumerate(sorted_lines, 1):
        seqid, s, e, _ = _line_coords(ln)
        seq = fetch_seq(fh, fai, seqid, s, e)
        if strand == "-":
            seq = revcomp(seq)
        out.append((f"{label}_{i}", seq))
    return out


def build_structure_fasta(mrna_id, gene_line, mrna_line,
                          exon_lines, cds_lines, utr5_lines, utr3_lines,
                          intron_lines, fh, fai):
    """Return list of (header, seq) for structure.fa."""
    strand = "+"
    anchor = gene_line or mrna_line
    if anchor:
        strand = anchor.split("\t")[6]

    entries = []

    if gene_line:
        seqid, s, e, _ = _line_coords(gene_line)
        seq = fetch_seq(fh, fai, seqid, s, e)
        if strand == "-":
            seq = revcomp(seq)
        entries.append((f"{mrna_id}_gene", seq))

    mrna_seq = _concat_seq(exon_lines, fh, fai, strand)
    if mrna_seq:
        entries.append((f"{mrna_id}_mRNA", mrna_seq))

    utr5_seq = _concat_seq(utr5_lines, fh, fai, strand)
    if utr5_seq:
        entries.append((f"{mrna_id}_five_prime_UTR", utr5_seq))

    utr3_seq = _concat_seq(utr3_lines, fh, fai, strand)
    if utr3_seq:
        entries.append((f"{mrna_id}_three_prime_UTR", utr3_seq))

    entries.extend(
        (f"{mrna_id}_{suf}", seq)
        for suf, seq in _individual_seqs(exon_lines, fh, fai, strand, "exon")
    )
    entries.extend(
        (f"{mrna_id}_{suf}", seq)
        for suf, seq in _individual_seqs(cds_lines, fh, fai, strand, "CDS")
    )
    entries.extend(
        (f"{mrna_id}_{suf}", seq)
        for suf, seq in _individual_seqs(intron_lines, fh, fai, strand, "intron")
    )
    return entries


def write_fasta(path, entries, line_width=80):
    with open(path, "w") as fh:
        for i, (header, seq) in enumerate(entries):
            if i > 0:
                fh.write("\n")
            fh.write(f">{header}\n")
            if not seq:
                continue
            for j in range(0, len(seq), line_width):
                fh.write(seq[j:j + line_width] + "\n")


# ---------------------------------------------------------------------------
# GenBank writer (Benchling-ready; coords are local 1-based inclusive)
# ---------------------------------------------------------------------------

import datetime


def _local_intervals(lines, gene_start, gene_end, strand):
    """Shift genomic intervals into local (gene-5'-origin) 1-based coords."""
    out = []
    for ln in lines:
        p = ln.split("\t")
        try:
            s, e = int(p[3]), int(p[4])
        except (ValueError, IndexError):
            continue
        if strand == "+":
            out.append((s - gene_start + 1, e - gene_start + 1))
        else:
            out.append((gene_end - e + 1, gene_end - s + 1))
    out.sort()
    return out


def _gb_location(intervals):
    """GenBank location string: single interval or join(a..b,c..d,...)."""
    if not intervals:
        return ""
    if len(intervals) == 1:
        return f"{intervals[0][0]}..{intervals[0][1]}"
    return "join(" + ",".join(f"{s}..{e}" for s, e in intervals) + ")"


def _gb_feature(lines_out, feat_type, location, qualifiers):
    """Emit a GenBank FEATURES entry with 21-char key column and /qual="val" lines."""
    if not location:
        return
    # Feature key starts at col 6, location at col 22 (21-char total left pad)
    lines_out.append(f"     {feat_type:<16}{location}")
    for k, v in qualifiers:
        lines_out.append(f"""                     /{k}="{v}" """.rstrip())


def fetch_flanked_gene(fh, fai, gene_line, flank_bp):
    """Fetch the gene sequence plus requested flanks in the gene's 5'→3' frame.
    Returns (seq, up_flank_actual, dn_flank_actual). Clamps at chromosome boundaries,
    so the effective flanks may be smaller than requested near the chromosome ends."""
    p = gene_line.split("\t")
    seqid, gene_start, gene_end, strand = p[0], int(p[3]), int(p[4]), p[6]
    chrom_len = fai.get(seqid, (0,))[0] if fai else 0
    fetch_start = max(1, gene_start - flank_bp)
    fetch_end = gene_end + flank_bp
    if chrom_len:
        fetch_end = min(chrom_len, fetch_end)
    seq = fetch_seq(fh, fai, seqid, fetch_start, fetch_end)
    up = gene_start - fetch_start     # actual upstream (genomic frame)
    dn = fetch_end - gene_end         # actual downstream (genomic frame)
    if strand == "-":
        seq = revcomp(seq)
        up, dn = dn, up               # in gene frame, swap upstream/downstream
    return seq, up, dn


def build_genbank(mrna_id, gene_line, exon_lines, cds_lines,
                  utr5_lines, utr3_lines, intron_lines, backbone_seq,
                  up_flank=0, dn_flank=0, organism="Solanum melongena"):
    """Return a GenBank-format string for this mRNA. backbone_seq must already be
    5'→3' on the gene strand (reverse-complemented for minus-strand genes) and
    include up_flank upstream + gene + dn_flank downstream bases, in that order."""
    if not gene_line or not backbone_seq:
        return ""
    p = gene_line.split("\t")
    seqid, gene_start, gene_end, strand = p[0], int(p[3]), int(p[4]), p[6]
    gene_len = gene_end - gene_start + 1
    total_len = len(backbone_seq)

    def shift(intervals):
        return [(s + up_flank, e + up_flank) for s, e in intervals]

    exons_loc = shift(_local_intervals(exon_lines, gene_start, gene_end, strand))
    cds_loc = shift(_local_intervals(cds_lines, gene_start, gene_end, strand))
    utr5_loc = shift(_local_intervals(utr5_lines, gene_start, gene_end, strand))
    utr3_loc = shift(_local_intervals(utr3_lines, gene_start, gene_end, strand))
    intr_loc = shift(_local_intervals(intron_lines, gene_start, gene_end, strand))

    date = datetime.date.today().strftime("%d-%b-%Y").upper()
    name_field = mrna_id if len(mrna_id) <= 16 else mrna_id
    pad = " " * max(1, 17 - len(name_field))
    locus = f"LOCUS       {name_field}{pad}{total_len:>9} bp    DNA     linear   PLN {date}"

    flank_desc = (f" +{up_flank}bp 5' / +{dn_flank}bp 3' flanks"
                  if (up_flank or dn_flank) else "")
    out = [
        locus,
        f"DEFINITION  Gene structure for {mrna_id} "
        f"({seqid}:{gene_start}..{gene_end} strand {strand}){flank_desc}.",
        f"ACCESSION   {mrna_id}",
        f"VERSION     {mrna_id}",
        "KEYWORDS    .",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
        "FEATURES             Location/Qualifiers",
    ]

    fetch_start = gene_start - up_flank
    fetch_end = gene_end + dn_flank
    _gb_feature(out, "source", f"1..{total_len}", [
        ("organism", organism),
        ("mol_type", "genomic DNA"),
        ("note", f"{seqid}:{fetch_start}..{fetch_end} strand={strand} "
                 f"gene={gene_start}..{gene_end} flanks=+{up_flank}/-{dn_flank}"),
    ])
    gene_loc_s = up_flank + 1
    gene_loc_e = up_flank + gene_len
    _gb_feature(out, "gene", f"{gene_loc_s}..{gene_loc_e}",
                [("locus_tag", mrna_id), ("label", "gene")])

    if up_flank > 0:
        _gb_feature(out, "misc_feature", f"1..{up_flank}",
                    [("label", f"upstream_{up_flank}bp"),
                     ("note", "5' flanking sequence")])
    if dn_flank > 0:
        _gb_feature(out, "misc_feature", f"{gene_loc_e + 1}..{total_len}",
                    [("label", f"downstream_{dn_flank}bp"),
                     ("note", "3' flanking sequence")])

    if exons_loc:
        _gb_feature(out, "mRNA", _gb_location(exons_loc),
                    [("product", "spliced mRNA"), ("label", "mRNA")])
    if cds_loc:
        _gb_feature(out, "CDS", _gb_location(cds_loc),
                    [("product", "protein"), ("label", "CDS")])
    for i, (s, e) in enumerate(exons_loc, 1):
        _gb_feature(out, "exon", f"{s}..{e}",
                    [("number", str(i)), ("label", f"exon {i}")])
    for s, e in utr5_loc:
        _gb_feature(out, "5'UTR", f"{s}..{e}", [("label", "5'UTR")])
    for s, e in utr3_loc:
        _gb_feature(out, "3'UTR", f"{s}..{e}", [("label", "3'UTR")])
    for i, (s, e) in enumerate(intr_loc, 1):
        _gb_feature(out, "intron", f"{s}..{e}",
                    [("number", str(i)), ("label", f"intron {i}")])

    out.append("ORIGIN")
    seq = backbone_seq.lower()
    for pos in range(0, len(seq), 60):
        chunk = seq[pos:pos + 60]
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        out.append(f"{pos + 1:>9} {' '.join(groups)}")
    out.append("//")
    return "\n".join(out) + "\n"


def write_genbank(path, *args, **kwargs):
    content = build_genbank(*args, **kwargs)
    if content:
        Path(path).write_text(content, encoding="ascii")


# ---------------------------------------------------------------------------
# Per-mRNA output writer
# ---------------------------------------------------------------------------

def _write_gff_file(path, header, lines):
    with open(path, "w") as fh:
        fh.write(header)
        for ln in lines:
            fh.write(ln + "\n")


def write_structure(mrna_id, gene_records, mrna_records, child_records,
                    mrna_to_gene, out_dir, overwrite,
                    fasta_fh=None, fai=None, flank_bp=1000,
                    organism="Solanum melongena"):
    safe_id = re.sub(r"[^\w.\-]", "_", mrna_id)
    mrna_dir = out_dir / safe_id
    if mrna_dir.exists() and not overwrite:
        return

    mrna_dir.mkdir(parents=True, exist_ok=True)

    gene_id = mrna_to_gene.get(mrna_id, "")
    gene_line = gene_records.get(gene_id, "")
    mrna_line = mrna_records.get(mrna_id, "")
    children = child_records.get(mrna_id, [])

    by_type = defaultdict(list)
    for ln in children:
        p = ln.split("\t")
        if len(p) >= 3:
            by_type[p[2].lower()].append(ln)

    exon_lines = by_type["exon"]
    cds_lines = by_type["cds"]
    utr5_lines = by_type["five_prime_utr"]
    utr3_lines = by_type["three_prime_utr"]
    intron_lines = _derive_introns(exon_lines, gene_line or mrna_line)

    header = f"##gff-version 3\n# Gene structure: {mrna_id}\n"

    all_lines = []
    if gene_line:
        all_lines.append(gene_line)
    if mrna_line:
        all_lines.append(mrna_line)
    all_lines.extend(children)
    all_lines.extend(intron_lines)

    _write_gff_file(mrna_dir / "structure.gff3", header, all_lines)
    _write_gff_file(mrna_dir / "exons.gff3", header, exon_lines)
    _write_gff_file(mrna_dir / "cds.gff3", header, cds_lines)
    _write_gff_file(mrna_dir / "five_prime_utr.gff3", header, utr5_lines)
    _write_gff_file(mrna_dir / "three_prime_utr.gff3", header, utr3_lines)
    _write_gff_file(mrna_dir / "introns.gff3", header, intron_lines)

    if fasta_fh is not None and fai:
        entries = build_structure_fasta(
            mrna_id, gene_line, mrna_line,
            exon_lines, cds_lines, utr5_lines, utr3_lines,
            intron_lines, fasta_fh, fai,
        )
        if entries:
            write_fasta(mrna_dir / "structure.fa", entries)
            if gene_line:
                backbone, up, dn = fetch_flanked_gene(fasta_fh, fai, gene_line, flank_bp)
                if backbone:
                    gb_content = build_genbank(
                        mrna_id, gene_line,
                        exon_lines, cds_lines, utr5_lines, utr3_lines,
                        intron_lines, backbone,
                        up_flank=up, dn_flank=dn,
                        organism=organism,
                    )
                    if gb_content:
                        (mrna_dir / "structure.gb").write_text(gb_content, encoding="ascii")
                        # Sibling copy next to the mRNA folder, named for the mRNA —
                        # convenient for bulk GenBank import (one file = one sequence).
                        (out_dir / f"{safe_id}.gb").write_text(gb_content, encoding="ascii")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True, help="GFF/GFF3 annotation file")
    ap.add_argument("--output", required=True, help="Output root (e_GENE_Structures/)")
    ap.add_argument("--overwrite", default="true", help="Overwrite existing outputs (true/false)")
    ap.add_argument("--genome-fasta", default="",
                    help="Genomic FASTA (must have sibling .fai). Enables structure.fa output.")
    ap.add_argument("--flank-bp", type=int, default=1000,
                    help="bp of 5'/3' flanking sequence to include in structure.gb (default 1000; "
                         "clamped at chromosome boundaries).")
    ap.add_argument("--organism", default="Solanum melongena",
                    help="Organism name written into GenBank SOURCE/ORGANISM/source.organism")
    ap.add_argument("--gene-spec", action="append", default=[],
                    help="GENE_NAME:HIT_IDS_PATH (repeatable). Parses GFF once for all specs.")
    # Backward-compat single-gene flags
    ap.add_argument("--gene-name", default="", help="(single-gene mode) gene family name")
    ap.add_argument("--hit-ids", default="", help="(single-gene mode) hit IDs file")
    args = ap.parse_args()

    overwrite = args.overwrite.lower() not in ("false", "0", "no")

    specs = []
    for spec in args.gene_spec:
        if ":" not in spec:
            print(f"ERROR: --gene-spec must be GENE:PATH, got: {spec}", file=sys.stderr)
            sys.exit(1)
        gene_name, hit_ids_path = spec.split(":", 1)
        specs.append((gene_name.strip(), hit_ids_path.strip()))
    if args.gene_name and args.hit_ids:
        specs.append((args.gene_name, args.hit_ids))
    if not specs:
        print("ERROR: provide --gene-spec GENE:PATH (repeatable) or --gene-name + --hit-ids",
              file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(args.gff3):
        print(f"ERROR: GFF not found: {args.gff3}", file=sys.stderr)
        sys.exit(1)
    for _gn, hp in specs:
        if not os.path.isfile(hp):
            print(f"ERROR: hit-ids not found: {hp}", file=sys.stderr)
            sys.exit(1)

    fasta_fh = None
    fai = None
    if args.genome_fasta:
        if not os.path.isfile(args.genome_fasta):
            print(f"WARNING: --genome-fasta not found: {args.genome_fasta} (skipping FASTA output)",
                  file=sys.stderr)
        else:
            fai_path = args.genome_fasta + ".fai"
            if not os.path.isfile(fai_path):
                print(f"WARNING: .fai index missing for {args.genome_fasta}; "
                      f"run 'samtools faidx' to enable structure.fa output",
                      file=sys.stderr)
            else:
                fai = load_fai(fai_path)
                fasta_fh = open(args.genome_fasta, "rb")

    try:
        print(f"Loading GFF: {args.gff3}", file=sys.stderr)
        gene_records, mrna_records, child_records, mrna_to_gene, gene_to_mrnas = load_gff(args.gff3)

        total_written = 0
        for gene_name, hit_ids_path in specs:
            resolved = _resolve_hits(hit_ids_path, mrna_records, gene_records, gene_to_mrnas)
            if not resolved:
                print(f"WARNING: no GFF matches for {gene_name} ({hit_ids_path})", file=sys.stderr)
                continue

            out_dir = Path(args.output) / gene_name
            out_dir.mkdir(parents=True, exist_ok=True)

            written = 0
            for _hit, mrna_ids in resolved:
                for mrna_id in mrna_ids:
                    write_structure(mrna_id, gene_records, mrna_records, child_records,
                                    mrna_to_gene, out_dir, overwrite,
                                    fasta_fh=fasta_fh, fai=fai,
                                    flank_bp=args.flank_bp,
                                    organism=args.organism)
                    written += 1
            total_written += written
            fa_note = " (+ structure.fa)" if fasta_fh else ""
            print(f"  {gene_name}: {written} mRNA(s){fa_note} → {out_dir}/", file=sys.stderr)

        print(f"Done: {total_written} total gene structure folder(s)", file=sys.stderr)
    finally:
        if fasta_fh:
            fasta_fh.close()


if __name__ == "__main__":
    main()
