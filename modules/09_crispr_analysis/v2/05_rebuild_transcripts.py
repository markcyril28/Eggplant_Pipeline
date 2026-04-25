#!/usr/bin/env python3
"""
Module: 05_rebuild_transcripts.py
Stage [5] — Mutant transcript reconstruction (Biopython + GTF).

For each guide's top indel outcomes, applies the predicted indel to the
reference CDS sequence and reconstructs the resulting mRNA transcript.

Outputs per guide:
  - mutant_cds_{indel_rank}.fa  : mutant CDS FASTA
  - mutant_tx_{indel_rank}.fa   : mutant mRNA (UTR + CDS) if UTR coords available
  - transcript_summary.tsv      : columns added to guide table for each indel

Usage:
    python3 05_rebuild_transcripts.py \\
        --input        <indels.tsv>          \\
        --outdir       <output_dir>          \\
        --genome-fasta <genome.fa>           \\
        --gtf          <annotation.gtf>      \\
        --gene-group   DMP                   \\
        --top-indels   3
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Mutant transcript reconstruction")
    p.add_argument("--input",         required=True)
    p.add_argument("--outdir",        required=True)
    p.add_argument("--genome-fasta",  required=True)
    p.add_argument("--gtf",           default="",
                   help="GTF annotation file for CDS boundary extraction (optional; skips rebuild if absent)")
    p.add_argument("--gene-group",    default="DMP",
                   help="Family tag used as a fallback regex when --gene-id does not match anything")
    p.add_argument("--gene-id",       default="",
                   help="Specific target gene id (e.g. SMEL5_01g008730). Matched against GTF gene_id "
                        "with prefix semantics; falls back to --gene-group regex if no hit.")
    p.add_argument("--top-indels",    type=int, default=3)
    p.add_argument("--overwrite",     action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


# ─── GTF utilities ────────────────────────────────────────────────────────────

def load_cds_coords(gtf_path: str, gene_group: str, gene_id: str = "") -> dict:
    """
    Parse GTF and return {transcript_id: [(chrom, start, end, strand), ...]}
    for CDS features belonging to the target gene.

    Matching priority:
      1. Exact match on gene_id (e.g. "SMEL5_01g008730").
      2. Prefix match on gene_id (to catch variants like "SMEL5_01g008730.v2").
      3. Fallback regex on gene_group when gene_id is empty OR no prefix hit.

    The previous pipeline relied only on (3), which failed silently for
    eggplant annotations where gene_id does not contain the family tag
    "DMP" (e.g. the SMEL5/GPE001970 nomenclature), producing an empty
    coords dict and downstream NA-filled outputs.
    """
    import re
    coords: dict = {}
    fallback_pattern = re.compile(gene_group, re.IGNORECASE) if gene_group else None
    gene_id_clean = (gene_id or "").strip()
    matched_by_id = 0
    matched_by_fallback = 0

    def _parse_attrs(attrs: str) -> tuple[str, str]:
        """Extract (gene_id, transcript_id) from GTF or GFF3 attribute strings.

        GTF  : gene_id "X"; transcript_id "Y";
        GFF3 : ID=X;Parent=Y  (on a CDS row, Parent is the transcript;
               gene_id is derived by stripping the trailing ".N" isoform suffix
               typical of Helixer/SMEL5/GPE001970 naming).
        """
        m = re.search(r'gene_id "([^"]+)"', attrs)
        gid = m.group(1) if m else ""
        m = re.search(r'transcript_id "([^"]+)"', attrs)
        tid = m.group(1) if m else ""
        if not tid:
            m = re.search(r'transcript_id=([^;]+)', attrs)
            if m:
                tid = m.group(1).strip()
        if not tid:
            m = re.search(r'Parent=([^;,]+)', attrs)
            if m:
                tid = m.group(1).strip()
        if not gid:
            m = re.search(r'gene_id=([^;]+)', attrs)
            if m:
                gid = m.group(1).strip()
        if not gid and tid:
            m = re.match(r'^(.+)\.\d+$', tid)
            gid = m.group(1) if m else tid
        return gid, tid

    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                attrs = parts[8]
                gid, tid = _parse_attrs(attrs)

                hit = False
                if gene_id_clean:
                    if gid == gene_id_clean or gid.startswith(gene_id_clean + ".") \
                       or gid.startswith(gene_id_clean + "_"):
                        hit = True
                        matched_by_id += 1
                if not hit and fallback_pattern is not None and not gene_id_clean:
                    # Only fall back to the family regex when no specific gene_id
                    # was supplied — otherwise an unrelated gene with the family
                    # tag in its id would pollute this gene's CDS map.
                    if fallback_pattern.search(gid):
                        hit = True
                        matched_by_fallback += 1
                if not hit:
                    continue

                coords.setdefault(tid or gid, []).append((
                    parts[0], int(parts[3]) - 1, int(parts[4]),
                    parts[6]  # strand
                ))
    except FileNotFoundError:
        _log(f"[05_rebuild] GTF not found: {gtf_path}", level="WARN")
        return coords

    # Second-chance fallback: the specific id was given but nothing matched
    # (e.g. gene_id column uses a different nomenclature than the target
    # FASTA name). Re-scan with the family regex so downstream stages still
    # produce non-empty output instead of silently collapsing to NA.
    if gene_id_clean and not coords and fallback_pattern is not None:
        _log(f"[05_rebuild] No GTF CDS matched gene_id={gene_id_clean!r}; "
             f"retrying with family regex {gene_group!r}.", level="WARN")
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                gid, tid = _parse_attrs(parts[8])
                if not fallback_pattern.search(gid):
                    continue
                coords.setdefault(tid or gid, []).append((
                    parts[0], int(parts[3]) - 1, int(parts[4]),
                    parts[6]
                ))
                matched_by_fallback += 1

    if gene_id_clean:
        _log(f"[05_rebuild] GTF match for gene_id={gene_id_clean!r}: "
             f"{matched_by_id} CDS rows by id, {matched_by_fallback} by fallback; "
             f"{len(coords)} transcript(s) loaded.", level="INFO")
    return coords


def extract_cds_sequence(genome_fasta: str, exons: list, strand: str) -> str:
    """Extract and concatenate CDS exon sequences from genome FASTA.

    Uses SeqIO.index() for random access without loading the entire genome
    into memory (important when multiple parallel instances run concurrently).
    """
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq

        genome = SeqIO.index(genome_fasta, "fasta")
        cds_seq = ""
        for chrom, start, end, _strand in sorted(exons, key=lambda x: x[1]):
            if chrom not in genome:
                continue
            cds_seq += str(genome[chrom].seq[start:end])
        genome.close()
        if strand == "-":
            cds_seq = str(Seq(cds_seq).reverse_complement())
        return cds_seq.upper()
    except ImportError:
        _log("[05_rebuild] Biopython not installed — cannot extract sequences.", level="WARN")
        return ""
    except Exception as exc:
        _log(f"[05_rebuild] Sequence extraction error: {exc}", level="ERROR")
        return ""


# ─── Indel application ────────────────────────────────────────────────────────

def apply_indel(cds: str, cut_position: int, indel_str: str) -> str:
    """
    Apply a simple indel to the CDS string.
    indel_str format (inDelphi/Lindel conventions):
      "+N:ACGT"  — insertion of ACGT after cut_position
      "-N"       — deletion of N bases starting at cut_position
    """
    if not indel_str or not cds:
        return cds
    cut = min(cut_position, len(cds))
    indel = indel_str.strip()
    if indel.startswith("+"):
        # Insertion
        parts = indel.split(":", 1)
        inserted = parts[1] if len(parts) > 1 else ""
        return cds[:cut] + inserted + cds[cut:]
    elif indel.startswith("-"):
        n_del = 0
        try:
            n_del = int(indel[1:].split(":")[0])
        except ValueError:
            pass
        return cds[:cut] + cds[cut + n_del:]
    return cds


def is_frameshift(indel_str: str) -> bool:
    """Return True if the indel length is not divisible by 3."""
    indel = indel_str.strip()
    n = 0
    if indel.startswith("+"):
        parts = indel.split(":", 1)
        n = len(parts[1]) if len(parts) > 1 else 0
    elif indel.startswith("-"):
        try:
            n = int(indel[1:].split(":")[0])
        except ValueError:
            pass
    return n % 3 != 0


# ─── Cut-site position estimate ───────────────────────────────────────────────

def estimate_cut_position(cds: str, guide_seq: str) -> int:
    """Locate guide in CDS and return SpCas9 cut position (3 nt from PAM).

    Returns -1 when the guide cannot be unambiguously located on either
    strand of the CDS. Callers MUST treat -1 as a hard skip and emit NA
    for the affected indel row — NEVER use a midpoint or arbitrary
    fallback. A synthetic cut position silently corrupts all downstream
    stages (mutant CDS, PTC position, NMD class, composite KO score)
    with numerically plausible but biologically meaningless predictions.
    """
    guide = guide_seq.upper()[:20]
    if not guide or not cds:
        return -1
    pos = cds.find(guide)
    if pos != -1:
        return pos + 17  # SpCas9 cuts between nt 17 and 18 of the protospacer
    # Try reverse complement (guide on opposite strand of CDS-sense mRNA)
    try:
        from Bio.Seq import Seq
        rc = str(Seq(guide).reverse_complement())
        pos = cds.find(rc)
        if pos != -1:
            return pos + 3  # 3 nt from PAM on the minus strand
    except ImportError:
        pass
    # Guide not findable on either strand — refuse to guess.
    return -1


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    inpath   = Path(args.input)
    outdir   = Path(args.outdir)
    fasta_dir= outdir / "mutant_fastas"
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(exist_ok=True)

    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.transcripts.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[05_rebuild] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    # Load guides
    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    # Derive the specific target gene id from the input filename stem when
    # --gene-id was not explicitly supplied. Each stage produces one TSV per
    # gene (e.g. SMEL5_01g008730.indels.tsv), so the stem is the gene id.
    gene_id_effective = args.gene_id or gene_stem

    # Load CDS coordinates from GTF
    if not args.gtf or not Path(args.gtf).exists():
        _log(f"[05_rebuild] GTF not provided or not found ({args.gtf!r}) — "
              "writing pass-through table with NA columns.", level="WARN")
        cds_coords = {}
    else:
        cds_coords = load_cds_coords(args.gtf, args.gene_group, gene_id_effective)

    # For each guide use the first matching transcript's CDS
    # (In production, map guide to specific transcript via locus coordinates)
    first_tx    = next(iter(cds_coords), None)
    ref_cds_seq = ""
    if first_tx and cds_coords[first_tx]:
        strand = cds_coords[first_tx][0][3]
        ref_cds_seq = extract_cds_sequence(
            args.genome_fasta, cds_coords[first_tx], strand)

    # New columns
    new_cols  = []
    for i in range(1, args.top_indels + 1):
        new_cols += [f"indel{i}_seq", f"indel{i}_is_frameshift",
                     f"indel{i}_ptc_pos", f"indel{i}_ptc_exon"]
    new_fields = fields + new_cols

    # Pre-filter: drop guides that fall outside the CDS (UTR / intron /
    # upstream-downstream flanking regions). These were previously kept in
    # the output with NA indel columns, which inflated downstream ranking
    # files with biologically meaningless rows. A reference CDS is required
    # for the check — if one is unavailable (no GTF), all rows pass through
    # so downstream stages can still run on best-effort context.
    filtered_rows: list[dict] = []
    outside_cds = 0
    for row in rows:
        if not ref_cds_seq:
            filtered_rows.append(row)
            continue
        guide_check = ""
        for col in ("targetSeq", "guideSeq", "sequence", "Sequence"):
            if col in row and row[col]:
                guide_check = row[col].strip().upper()[:20]
                break
        if estimate_cut_position(ref_cds_seq, guide_check) >= 0:
            filtered_rows.append(row)
        else:
            outside_cds += 1
    if outside_cds:
        _log(f"[05_rebuild] Filtered {outside_cds} guide(s) outside CDS "
             f"(stem={gene_stem}); {len(filtered_rows)} retained for "
             "transcript reconstruction.", level="INFO")
    rows = filtered_rows

    for row in rows:
        guide_seq = ""
        for col in ("targetSeq", "guideSeq", "sequence", "Sequence"):
            if col in row and row[col]:
                guide_seq = row[col].strip().upper()[:20]
                break

        # Parse top indels (try inDelphi first, then Lindel)
        top_indels_json = row.get("inDelphi_top_indels") or row.get("Lindel_top_indels") or "[]"
        try:
            top_indels = json.loads(top_indels_json)[:args.top_indels]
        except (json.JSONDecodeError, TypeError):
            top_indels = []

        cut_pos = estimate_cut_position(ref_cds_seq, guide_seq) if ref_cds_seq else -1
        # All retained rows are locatable on the CDS after the pre-filter
        # above; this branch is kept only for the no-GTF pass-through case.
        guide_locatable = cut_pos >= 0

        guide_id = ""
        for col in ("guideId", "guide_id", "name", "ID"):
            if col in row and row[col]:
                guide_id = row[col].strip()
                break

        for idx in range(1, args.top_indels + 1):
            col_seq    = f"indel{idx}_seq"
            col_fs     = f"indel{idx}_is_frameshift"
            col_ptc    = f"indel{idx}_ptc_pos"
            col_exon   = f"indel{idx}_ptc_exon"

            if idx - 1 < len(top_indels) and ref_cds_seq and guide_locatable:
                indel_str = top_indels[idx - 1].get("indel", "")
                mutant    = apply_indel(ref_cds_seq, cut_pos, indel_str)

                # Write mutant CDS FASTA
                fa_name = f"{guide_id}_indel{idx}.fa" if guide_id else f"guide{idx}_indel{idx}.fa"
                fa_path = fasta_dir / fa_name
                with open(fa_path, "w") as fa:
                    fa.write(f">{guide_id or 'guide'}_indel{idx} | {indel_str}\n")
                    for i in range(0, len(mutant), 80):
                        fa.write(mutant[i:i+80] + "\n")
                    fa.write("\n")

                # Find PTC: first stop codon in frame
                ptc_pos  = "NA"
                ptc_exon = "NA"
                for frame_start in range(0, len(mutant) - 2, 3):
                    codon = mutant[frame_start:frame_start + 3]
                    if codon in ("TAA", "TAG", "TGA"):
                        ptc_pos = str(frame_start)
                        # Rough exon number: map position back to sorted exons
                        # (simplified; uses cumulative exon lengths)
                        ptc_exon = "exon_unknown"
                        if first_tx and cds_coords[first_tx]:
                            cum = 0
                            for en, (_, s, e, _st) in enumerate(
                                    sorted(cds_coords[first_tx], key=lambda x: x[1]), 1):
                                cum += e - s
                                if frame_start < cum:
                                    ptc_exon = f"exon{en}"
                                    break
                        break

                # Store path relative to outdir so the TSV is portable across
                # machines and PIPELINE_DIR moves (stage 6 resolves it back).
                row[col_seq]  = str(fa_path.relative_to(outdir))
                row[col_fs]   = str(is_frameshift(indel_str))
                row[col_ptc]  = ptc_pos
                row[col_exon] = ptc_exon
            else:
                row[col_seq]  = "NA"
                row[col_fs]   = "NA"
                row[col_ptc]  = "NA"
                row[col_exon] = "NA"

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    _log(f"[05_rebuild] {len(rows)} guides processed -> {outpath}", level="INFO")
    _log(f"[05_rebuild] Mutant FASTA directory: {fasta_dir}", level="INFO")


if __name__ == "__main__":
    main()
