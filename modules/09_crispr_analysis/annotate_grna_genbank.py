#!/usr/bin/env python3
"""
annotate_grna_genbank.py

Annotate per-gene GenBank gene-structure files with two complementary tracks:

  (A) Secondary-structure features  -- N-terminal domain, 4 transmembrane
      helices (TM1..TM4), and the extracellular beta-pleated sheet -- whose
      amino-acid boundaries are taken from the motif-analysis reference
      (config/colors_config/meme_motif_colors.toml :: dmp_interaction.
      domain_ranges_json; boundaries derived by PyMOL DSS on the AlphaFold3
      model of the reference SmelDMP01.030 / SMEL5_01g026030.1, 178 aa).
      For non-reference paralogs the reference AA boundaries are projected
      onto the paralog's own CDS by proportional scaling
      (paralog_aa = round(ref_aa * paralog_len / reference_len)) and the
      feature /note records the projection.

  (B) gRNA features derived from CRISPR-P V2 filtered score tables, augmented
      with genome-wide off-target mismatch counts from Cas-OFFinder.

For each gene in --gb-dir, two GenBank outputs are written under --output-dir:
  01_Score_<HIGH>_above/<gene>_score_ge<HIGH>_mm<=<MAX_MM>.gb
      strict shortlist (score >= HIGH only)
  02_Score_<MOD>_above_incl_<HIGH>/<gene>_score_ge<MOD>_incl_ge<HIGH>_mm<=<MAX_MM>.gb
      combined list (score >= MOD; guides also passing the HIGH threshold are
      tagged "_HIGH>=<HIGH>" in the feature label so the strict subset stays
      visible in the combined file)

Each gRNA becomes a misc_feature with:
  /label="<guide>_S<score>_mm0=<n0>_mm1=<n1>_mm2=<n2>_mm3=<n3>[_HIGH>=<HIGH>]"
  /sgRNA_score, /sequence, /strand, /off_target_summary, /note

Off-target hits with strictly more than --max-mm mismatches are dropped before
counting (default max_mm = 3).

Inputs are intentionally robust to gene-name vs. Cas-OFFinder-filename mismatch:
all Cas-OFFinder output rows are pooled into a single sequence -> mm-count
dictionary, so on-target/off-target accounting is by gRNA *sequence*, not by
filename.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# ---------------------------------------------------------------------------
# DMP secondary-structure reference
# ---------------------------------------------------------------------------
# Reference: SmelDMP01.030 / SMEL5_01g026030.1, 178 aa.
# Boundaries derived from PyMOL DSS on the AlphaFold3 model of the reference,
# transcribed verbatim from config/colors_config/meme_motif_colors.toml ::
# [dmp_interaction].domain_ranges_json (0-based, end-inclusive on the reference
# protein). Loops (gray) are intentionally omitted; only the structurally
# named features requested by the user are emitted.
DMP_REFERENCE_GENE  = "SMEL5_01g026030.1"
DMP_REFERENCE_AA_LEN = 178

# mRNA-view companion is intentionally emitted ONLY for the SmelDMPv5_10.610
# transcript (SMEL5_10g017610.1) per author decision (2026-05-25). Other
# paralogs are not flanked with a spliced .gb. To add another gene, append its
# transcript id here.
MRNA_VIEW_TRANSCRIPTS = frozenset({"SMEL5_10g017610.1"})
DMP_SS_REFERENCE: List[Tuple[str, int, int, str]] = [
    # (label,                              aa_start_0based, aa_end_0based_incl, color_hex)
    ("N-terminal Domain",                    0,    9, "#D9A621"),
    ("TM1 (alpha4)",                        10,   50, "#D43005"),
    ("TM2 (alpha5 N-half)",                 54,   68, "#D43005"),
    ("Extracellular beta-pleated sheet",    69,   89, "#F5DEB3"),
    ("TM3 (alpha5 C-half)",                 98,  138, "#D43005"),
    ("TM4 (alpha6)",                       143,  171, "#D43005"),
]


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# GenBank parsing (minimal; we only need to split header / features / ORIGIN)
# ---------------------------------------------------------------------------

def split_genbank(text: str) -> Tuple[str, str, str, str]:
    """Return (locus_line, header_to_features_inclusive, features_block, origin_block).

    header_to_features_inclusive ends with the 'FEATURES ...' line (with newline).
    features_block contains the existing feature entries (between FEATURES line
    and ORIGIN line, exclusive).
    origin_block starts at the 'ORIGIN' line through the closing '//'.
    """
    lines = text.splitlines(keepends=True)
    features_start: Optional[int] = None
    origin_start: Optional[int] = None
    locus_line = ""
    for i, line in enumerate(lines):
        if line.startswith("LOCUS"):
            locus_line = line.rstrip("\n")
        if features_start is None and line.startswith("FEATURES"):
            features_start = i
        if origin_start is None and line.startswith("ORIGIN"):
            origin_start = i
            break
    if features_start is None or origin_start is None:
        raise ValueError("GenBank file missing FEATURES or ORIGIN block")
    header = "".join(lines[: features_start + 1])
    features = "".join(lines[features_start + 1 : origin_start])
    origin = "".join(lines[origin_start:])
    return locus_line, header, features, origin


def extract_origin_sequence(origin_block: str) -> str:
    """Concatenate the lowercase nucleotide letters from an ORIGIN block."""
    seq_chars = []
    for line in origin_block.splitlines():
        if line.startswith("ORIGIN") or line.startswith("//"):
            continue
        seq_chars.extend(ch for ch in line if ch.isalpha())
    return "".join(seq_chars).upper()


def extract_locus_length(locus_line: str) -> int:
    m = re.search(r"\s(\d+)\s+bp", locus_line)
    return int(m.group(1)) if m else 0


def extract_cds_intervals(features_block: str) -> List[Tuple[int, int]]:
    """Return the CDS as an ordered list of (start, end) 1-based inclusive
    gene-coordinate intervals.

    Handles both flat CDS (`CDS  179..865`) and spliced join CDS
    (`CDS  join(752..1194,2030..2255)`). Returns [] if no usable CDS is found.

    A leading `complement(...)` wrapper is tolerated but ignored: the gene-
    structure GenBank files already serve the sense strand, so the listed
    coordinates are directly usable.
    """
    # Match the CDS feature line plus any continuation lines that belong to
    # the same location (qualifier lines start with "/", location-continuation
    # lines don't). Stop at the next feature line (5 spaces + non-space) or at
    # a top-level keyword (ORIGIN / //). Use literal spaces in the lookahead so
    # the wider `\s` class doesn't gobble qualifier continuations.
    pat = re.compile(
        r"^ {5}CDS {2,}([^\n]*(?:\n {21}(?!/)[^\n]*)*)",
        re.MULTILINE,
    )
    m = pat.search(features_block)
    if not m:
        return []
    raw = re.sub(r"\s+", "", m.group(1))
    if raw.startswith("complement(") and raw.endswith(")"):
        raw = raw[len("complement("):-1]
    if raw.startswith("join(") and raw.endswith(")"):
        raw = raw[len("join("):-1]
    intervals: List[Tuple[int, int]] = []
    for piece in raw.split(","):
        m2 = re.match(r"(\d+)\.\.(\d+)", piece)
        if not m2:
            continue
        intervals.append((int(m2.group(1)), int(m2.group(2))))
    return intervals


def cds_total_length(intervals: List[Tuple[int, int]]) -> int:
    return sum(e - s + 1 for s, e in intervals)


def cds_offset_to_gene(intervals: List[Tuple[int, int]], offset: int
                       ) -> Optional[int]:
    """Map a 0-based offset within the spliced CDS to a 1-based gene-coord
    position. Returns None if the offset falls past the last exon.
    """
    remaining = offset
    for s, e in intervals:
        length = e - s + 1
        if remaining < length:
            return s + remaining
        remaining -= length
    return None


def aa_range_to_gene_intervals(aa_start_0: int, aa_end_0_incl: int,
                               cds_intervals: List[Tuple[int, int]]
                               ) -> List[Tuple[int, int]]:
    """Project an AA range onto gene coordinates, returning one or more
    contiguous gene-coordinate intervals (multiple when the AA range crosses
    one or more introns).
    """
    nt_start_off = aa_start_0 * 3            # 0-based inclusive into CDS
    nt_end_off   = (aa_end_0_incl + 1) * 3 - 1  # 0-based inclusive into CDS
    out: List[Tuple[int, int]] = []
    cursor = 0
    for s, e in cds_intervals:
        ex_len = e - s + 1
        ex_first = cursor
        ex_last  = cursor + ex_len - 1
        # Intersect [nt_start_off, nt_end_off] with [ex_first, ex_last]
        lo = max(nt_start_off, ex_first)
        hi = min(nt_end_off,   ex_last)
        if lo <= hi:
            gene_lo = s + (lo - ex_first)
            gene_hi = s + (hi - ex_first)
            out.append((gene_lo, gene_hi))
        cursor += ex_len
        if cursor > nt_end_off:
            break
    return out


def format_location(intervals: List[Tuple[int, int]]) -> str:
    """Render intervals as a GenBank location string."""
    if not intervals:
        return ""
    if len(intervals) == 1:
        s, e = intervals[0]
        return f"{s}..{e}"
    return "join(" + ",".join(f"{s}..{e}" for s, e in intervals) + ")"


def project_aa_range(aa_start_0: int, aa_end_0_incl: int,
                     paralog_aa_len: int, reference_aa_len: int
                     ) -> Tuple[int, int]:
    """Scale a reference AA range onto a paralog of different length.
    Endpoints are rounded; the result is clamped to the paralog's AA span.
    """
    if reference_aa_len <= 0:
        return aa_start_0, aa_end_0_incl
    ratio = paralog_aa_len / reference_aa_len
    s = int(round(aa_start_0 * ratio))
    e = int(round((aa_end_0_incl + 1) * ratio)) - 1
    s = max(0, min(s, paralog_aa_len - 1))
    e = max(s, min(e, paralog_aa_len - 1))
    return s, e


# ---------------------------------------------------------------------------
# Cas-OFFinder output parsing
# ---------------------------------------------------------------------------

def load_casoffinder_mm_counts(
    casoff_dir: Path, max_mm: int
) -> Dict[str, Dict[int, int]]:
    """Aggregate per-sequence mismatch-count histogram from every *_output.txt
    file found under casoff_dir.

    Returns:
        {sgRNA_sequence (uppercase, 23 nt incl PAM) -> {mm: hit_count}}
    Hits with mm > max_mm are excluded.
    """
    counts: Dict[str, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
    if not casoff_dir.is_dir():
        return counts

    for txt in sorted(casoff_dir.rglob("*_output.txt")):
        try:
            with txt.open() as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    query_seq = parts[0].strip().upper()
                    try:
                        mm = int(parts[5].strip())
                    except ValueError:
                        continue
                    if mm > max_mm:
                        continue
                    counts[query_seq][mm] += 1
        except OSError as exc:
            print(f"[WARN] could not read {txt}: {exc}", file=sys.stderr)
    return counts


# ---------------------------------------------------------------------------
# Filtered CSV parsing
# ---------------------------------------------------------------------------

GuideRow = Dict[str, str]


def load_filtered_csv(csv_path: Path) -> List[GuideRow]:
    rows: List[GuideRow] = []
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if not row.get("Sequence"):
                continue
            rows.append({
                "sgRNA_id": row.get("sgRNA_id", "").strip(),
                "Score":    row.get("Score", "0").strip(),
                "Sequence": row["Sequence"].strip().upper(),
                "strand":   row.get("strand", "+").strip(),
                "pos":      row.get("pos", "0").strip(),
            })
    return rows


def gene_id_from_csv_name(name: str) -> str:
    """'SMEL5_01g008730.1_filtered_0.7.csv' -> 'SMEL5_01g008730.1'."""
    base = Path(name).stem
    for suf in ("_filtered_0.7", "_filtered_0.5", "_filtered_0.4",
                "_filtered_0.3999", "_filtered_0.4999"):
        if base.endswith(suf):
            return base[: -len(suf)]
    base = re.sub(r"_filtered_[0-9.]+$", "", base)
    return base


# ---------------------------------------------------------------------------
# Locate a gRNA on the gene sequence
# ---------------------------------------------------------------------------

def locate_guide(gene_seq: str, guide_seq: str, strand: str,
                 pos_hint: int) -> Optional[Tuple[int, int, str]]:
    """Return (start, end, strand) in 1-based inclusive GenBank coords.

    Searches the gene sequence first for the guide on the indicated strand;
    falls back to the opposite strand if not found; finally falls back to the
    CSV pos_hint with our best-guess offset table:
      + strand: start = pos_hint - 19  (pos = last base of 20-nt protospacer)
      - strand: start = pos_hint + 1   (pos = +strand start of 23-mer - 1)
    """
    seq_len = len(gene_seq)
    glen = len(guide_seq)

    forward_hits = [m.start() for m in re.finditer(re.escape(guide_seq), gene_seq)]
    revcomp_seq = revcomp(guide_seq)
    reverse_hits = [m.start() for m in re.finditer(re.escape(revcomp_seq), gene_seq)]

    if strand == "+":
        primary, secondary = forward_hits, reverse_hits
        primary_strand, secondary_strand = "+", "-"
    else:
        primary, secondary = reverse_hits, forward_hits
        primary_strand, secondary_strand = "-", "+"

    def pick_closest(hits: List[int], hint_start_0based: int) -> int:
        if len(hits) == 1:
            return hits[0]
        return min(hits, key=lambda h: abs(h - hint_start_0based))

    if primary:
        if strand == "+":
            hint_start_0 = max(0, pos_hint - 20)
        else:
            hint_start_0 = max(0, pos_hint)
        start0 = pick_closest(primary, hint_start_0)
        return (start0 + 1, start0 + glen, primary_strand)

    if secondary:
        if secondary_strand == "+":
            hint_start_0 = max(0, pos_hint - 20)
        else:
            hint_start_0 = max(0, pos_hint)
        start0 = pick_closest(secondary, hint_start_0)
        return (start0 + 1, start0 + glen, secondary_strand)

    # Final fallback: trust CSV pos
    if strand == "+":
        start = pos_hint - 19
    else:
        start = pos_hint + 1
    end = start + glen - 1
    if start < 1 or end > seq_len:
        return None
    return (start, end, strand)


# ---------------------------------------------------------------------------
# Feature emission
# ---------------------------------------------------------------------------

def format_feature(
    guide_id: str,
    score: float,
    sequence: str,
    strand_actual: str,
    start: int,
    end: int,
    mm_counts: Dict[int, int],
    max_mm: int,
    high_tag: Optional[str],
) -> str:
    """Return a fixed-column GenBank misc_feature block (with trailing newline)."""
    loc = f"{start}..{end}" if strand_actual == "+" else f"complement({start}..{end})"

    mm_pairs = [f"mm{i}={mm_counts.get(i, 0)}" for i in range(max_mm + 1)]
    mm_short = "_".join(mm_pairs)
    mm_long = ", ".join(f"{i}mm={mm_counts.get(i, 0)}" for i in range(max_mm + 1))

    label_extra = f"_HIGH>={high_tag}" if high_tag else ""
    label = f"{guide_id}_S{score:.2f}_{mm_short}{label_extra}"

    note = (
        f"sgRNA {guide_id}; score={score:.4f}; strand={strand_actual}; "
        f"sequence={sequence}; off_target_counts(<= {max_mm}mm): {mm_long}"
    )
    if high_tag:
        note += f"; tier=HIGH (score>={high_tag})"

    qual_prefix = " " * 21 + "/"
    lines = [
        f"     {'misc_feature':<16}{loc}",
        f'{qual_prefix}label="{label}"',
        f'{qual_prefix}sgRNA_score="{score:.4f}"',
        f'{qual_prefix}sequence="{sequence}"',
        f'{qual_prefix}strand="{strand_actual}"',
        f'{qual_prefix}off_target_summary="{mm_short}"',
        f'{qual_prefix}note="{note}"',
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Secondary-structure feature emission
# ---------------------------------------------------------------------------

def build_ss_features(features_block: str, gene_id: str,
                      seq_len: int) -> Tuple[str, int]:
    """Return (ss_features_string, n_emitted) for one gene.

    Skips silently (returns ("", 0)) if no CDS is parseable from the
    existing FEATURES block. Spliced CDSs (join(...)) are supported: an SS
    range that crosses an intron is emitted as a GenBank join() location so
    the highlight ends at the donor splice site and resumes at the acceptor.
    """
    cds_intervals = extract_cds_intervals(features_block)
    if not cds_intervals:
        return "", 0

    cds_nt_len = cds_total_length(cds_intervals)
    paralog_aa_len = cds_nt_len // 3
    if cds_nt_len % 3 == 0 and paralog_aa_len > 0:
        paralog_aa_len -= 1  # drop the trailing stop codon

    is_reference = (gene_id == DMP_REFERENCE_GENE)

    qual_prefix = " " * 21 + "/"
    lines: List[str] = []
    emitted = 0
    for label, ref_s0, ref_e0, color_hex in DMP_SS_REFERENCE:
        if is_reference:
            aa_s0, aa_e0 = ref_s0, ref_e0
            projection_note = (
                f"reference coordinates from "
                f"meme_motif_colors.toml::dmp_interaction.domain_ranges_json "
                f"(PyMOL DSS on AlphaFold3 model of {DMP_REFERENCE_GENE}, "
                f"{DMP_REFERENCE_AA_LEN} aa)"
            )
        else:
            aa_s0, aa_e0 = project_aa_range(
                ref_s0, ref_e0, paralog_aa_len, DMP_REFERENCE_AA_LEN
            )
            projection_note = (
                f"projected from {DMP_REFERENCE_GENE} ({DMP_REFERENCE_AA_LEN} aa) "
                f"onto {gene_id} ({paralog_aa_len} aa) by proportional scaling; "
                f"reference range AA {ref_s0+1}-{ref_e0+1} (1-based)"
            )

        if aa_s0 >= paralog_aa_len:
            continue

        nt_pieces = aa_range_to_gene_intervals(aa_s0, aa_e0, cds_intervals)
        if not nt_pieces:
            continue
        # Clamp to gene span (paranoia; CDS already inside the gene span).
        nt_pieces = [(max(1, s), min(seq_len, e) if seq_len else e)
                     for s, e in nt_pieces if s <= e]
        if not nt_pieces:
            continue
        loc = format_location(nt_pieces)

        feature_label = (
            f"{label.replace(' ', '_')}_AA{aa_s0+1}-{aa_e0+1}"
        )
        spliced_note = (" (spliced across "
                        f"{len(nt_pieces)} exons)" if len(nt_pieces) > 1 else "")
        note = (
            f"{label}; aa {aa_s0+1}-{aa_e0+1} (1-based on {gene_id}){spliced_note}; "
            f"color={color_hex}; {projection_note}"
        )

        lines.append(f"     {'misc_feature':<16}{loc}")
        lines.append(f'{qual_prefix}label="{feature_label}"')
        lines.append(f'{qual_prefix}structural_element="{label}"')
        lines.append(f'{qual_prefix}color="{color_hex}"')
        lines.append(f'{qual_prefix}aa_range="{aa_s0+1}..{aa_e0+1}"')
        lines.append(f'{qual_prefix}note="{note}"')
        emitted += 1

    if not lines:
        return "", 0
    return "\n".join(lines) + "\n", emitted


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def build_annotated_gb(
    gb_text: str,
    guide_rows: List[GuideRow],
    mm_table: Dict[str, Dict[int, int]],
    max_mm: int,
    high_threshold: Optional[float],
    gene_id_hint: str = "",
) -> Tuple[str, int, int, int]:
    """Return (annotated_gb_text, n_added, n_skipped).

    gene_id_hint identifies the gene to the secondary-structure projector so
    the reference paralog (SMEL5_01g026030.1) keeps verbatim coordinates and
    other paralogs get a proportional projection.
    """
    locus_line, header, features_block, origin_block = split_genbank(gb_text)
    gene_seq = extract_origin_sequence(origin_block)
    seq_len = len(gene_seq) or extract_locus_length(locus_line)

    # Secondary-structure features come first so they render below the
    # built-in gene/CDS/exon stack and above the per-guide misc_features.
    ss_block, n_ss = build_ss_features(features_block, gene_id_hint, seq_len)

    new_features = []
    if ss_block:
        new_features.append(ss_block)
    n_added = 0
    n_skipped = 0
    for row in guide_rows:
        try:
            score = float(row["Score"])
        except ValueError:
            n_skipped += 1
            continue
        try:
            pos_hint = int(row["pos"])
        except ValueError:
            n_skipped += 1
            continue

        guide_seq = row["Sequence"]
        loc = locate_guide(gene_seq, guide_seq, row["strand"], pos_hint)
        if loc is None:
            print(f"[WARN] could not place guide {row['sgRNA_id']} ({guide_seq}) "
                  f"in gene (len={seq_len}); skipped", file=sys.stderr)
            n_skipped += 1
            continue
        start, end, strand_actual = loc

        mm_counts = mm_table.get(guide_seq, {})
        high_tag = None
        if high_threshold is not None and score >= high_threshold:
            high_tag = f"{high_threshold:g}"

        new_features.append(format_feature(
            guide_id=row["sgRNA_id"],
            score=score,
            sequence=guide_seq,
            strand_actual=strand_actual,
            start=start,
            end=end,
            mm_counts=mm_counts,
            max_mm=max_mm,
            high_tag=high_tag,
        ))
        n_added += 1

    annotated = header + features_block + "".join(new_features) + origin_block
    return annotated, n_added, n_skipped, n_ss


# ---------------------------------------------------------------------------
# mRNA-view companion writer
# ---------------------------------------------------------------------------
# For every annotated genomic GenBank we emit, also emit a sibling spliced
# mRNA view: introns removed, sequence = concatenated exons, every feature
# remapped from gene coords onto the mature transcript. Features that fall
# fully inside an intron are dropped; spliced CDS / SS join() locations
# collapse to a single contiguous interval in mRNA coords.

_FEATURE_KEY_COL = " " * 5
_QUAL_COL        = " " * 21


def iter_features(features_block: str):
    """Yield (key, location_str, qualifiers_text) for each feature.

    qualifiers_text preserves the original 21-space-indented lines verbatim
    (with trailing newline if present in source). Location-continuation lines
    are merged into location_str (whitespace collapsed).
    """
    lines = features_block.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if (line.startswith(_FEATURE_KEY_COL)
                and len(line) > 5 and line[5] != " "):
            m = re.match(r"^ {5}(\S+) +(.*)$", line)
            if not m:
                i += 1
                continue
            key = m.group(1)
            loc_parts = [m.group(2)]
            i += 1
            # Location-continuation lines: 21 spaces, NOT starting with '/'
            while (i < n and lines[i].startswith(_QUAL_COL)
                   and not lines[i][21:22].startswith("/")):
                loc_parts.append(lines[i].strip())
                i += 1
            location = "".join(loc_parts)
            qual_lines: List[str] = []
            # Qualifier lines: 21 spaces, starting with '/' (with continuations)
            while (i < n and lines[i].startswith(_QUAL_COL)
                   and lines[i][21:22].startswith("/")):
                qual_lines.append(lines[i])
                i += 1
                while (i < n and lines[i].startswith(_QUAL_COL)
                       and not lines[i][21:22].startswith("/")):
                    qual_lines.append(lines[i])
                    i += 1
            quals_text = "\n".join(qual_lines)
            if qual_lines:
                quals_text += "\n"
            yield key, location, quals_text
        else:
            i += 1


def parse_location(loc_str: str) -> Tuple[List[Tuple[int, int]], bool]:
    """Parse a GenBank location. Returns (intervals_1based_inclusive, complement).

    Tolerated: 'a..b', 'complement(...)', 'join(...)', nesting of the two,
    single positions ('123'), and fuzzy modifiers ('<', '>') which are stripped.
    """
    s = re.sub(r"\s+", "", loc_str)
    is_complement = False
    if s.startswith("complement(") and s.endswith(")"):
        is_complement = True
        s = s[len("complement("):-1]
    if s.startswith("join(") and s.endswith(")"):
        s = s[len("join("):-1]
    if s.startswith("order(") and s.endswith(")"):
        s = s[len("order("):-1]
    intervals: List[Tuple[int, int]] = []
    for piece in s.split(","):
        piece = piece.replace("<", "").replace(">", "")
        m = re.match(r"^(\d+)\.\.(\d+)$", piece)
        if m:
            intervals.append((int(m.group(1)), int(m.group(2))))
            continue
        m2 = re.match(r"^(\d+)$", piece)
        if m2:
            p = int(m2.group(1))
            intervals.append((p, p))
    return intervals, is_complement


def format_location_with_complement(intervals: List[Tuple[int, int]],
                                    complement: bool) -> str:
    if not intervals:
        return ""
    if len(intervals) == 1:
        s, e = intervals[0]
        loc = f"{s}..{e}"
    else:
        loc = "join(" + ",".join(f"{s}..{e}" for s, e in intervals) + ")"
    return f"complement({loc})" if complement else loc


def gene_to_mrna_intervals(gene_intervals: List[Tuple[int, int]],
                           exons: List[Tuple[int, int]]
                           ) -> List[Tuple[int, int]]:
    """Map gene-coord intervals onto mRNA coords (exons concatenated in order).

    Pieces of an interval that fall in introns are silently dropped. The result
    is post-merged so adjacent mRNA intervals (i.e. a feature that spans an
    exon-exon junction in genomic coords) collapse to one contiguous range.
    """
    out: List[Tuple[int, int]] = []
    for gs, ge in gene_intervals:
        cursor = 0
        for ex_s, ex_e in exons:
            ex_len = ex_e - ex_s + 1
            lo = max(gs, ex_s)
            hi = min(ge, ex_e)
            if lo <= hi:
                mrna_lo = cursor + (lo - ex_s) + 1
                mrna_hi = cursor + (hi - ex_s) + 1
                out.append((mrna_lo, mrna_hi))
            cursor += ex_len
    if not out:
        return []
    out.sort()
    merged = [out[0]]
    for s, e in out[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def format_origin_block(seq: str) -> str:
    """Return a standard 60-bp / 10-bp-group ORIGIN block ending with '//'."""
    seq = seq.lower()
    out_lines = ["ORIGIN\n"]
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60]
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        out_lines.append(f"{i + 1:>9} " + " ".join(groups) + "\n")
    out_lines.append("//\n")
    return "".join(out_lines)


def rewrite_locus_line(locus_line: str, new_length: int) -> str:
    """Update the LOCUS line: replace length and DNA->mRNA, preserving layout."""
    m = re.match(
        r"^(LOCUS\s+\S+\s+)(\d+)(\s+bp\s+)(\S+)(\s+.*)$",
        locus_line,
    )
    if m:
        return m.group(1) + str(new_length) + m.group(3) + "mRNA" + m.group(5)
    # Fallback: best-effort regex swap
    out = re.sub(r"\s\d+\s+bp", f"  {new_length} bp", locus_line, count=1)
    out = re.sub(r"\bDNA\b", "mRNA", out, count=1)
    return out


def _patch_source_qualifiers(quals_text: str, mrna_len: int) -> str:
    """Patch mol_type to mRNA and rewrite the /note coord span if present."""
    quals_text = re.sub(r'/mol_type="[^"]*"', '/mol_type="mRNA"', quals_text)
    # Keep the source /note as-is (it documents the gene's genomic span,
    # which is still useful provenance on the mRNA view).
    return quals_text


def build_mrna_view(annotated_gb_text: str) -> Optional[str]:
    """Return a spliced mRNA-view GenBank text, or None if exons can't be
    determined from the input."""
    locus_line, header, features_block, origin_block = split_genbank(annotated_gb_text)
    gene_seq = extract_origin_sequence(origin_block)

    # Source of exon intervals: prefer the parent mRNA feature; fall back to
    # the concatenation of exon features.
    exons: List[Tuple[int, int]] = []
    for key, loc, _ in iter_features(features_block):
        if key == "mRNA":
            ivs, _ = parse_location(loc)
            exons = ivs
            break
    if not exons:
        exon_list: List[Tuple[int, int]] = []
        for key, loc, _ in iter_features(features_block):
            if key == "exon":
                ivs, _ = parse_location(loc)
                exon_list.extend(ivs)
        if exon_list:
            exons = sorted(exon_list)
    if not exons:
        return None

    # Drop zero-length / out-of-range exons defensively.
    exons = [(s, e) for s, e in exons if 1 <= s <= e <= max(1, len(gene_seq))]
    if not exons:
        return None

    mrna_seq = "".join(gene_seq[s - 1:e] for s, e in exons)
    mrna_len = len(mrna_seq)
    if mrna_len == 0:
        return None

    new_feature_blocks: List[str] = []
    for key, loc, quals in iter_features(features_block):
        if key == "intron":
            continue
        intervals, complement = parse_location(loc)

        if key == "source":
            new_loc = f"1..{mrna_len}"
            new_quals = _patch_source_qualifiers(quals, mrna_len)
            new_feature_blocks.append(
                f"{_FEATURE_KEY_COL}{key:<16}{new_loc}\n{new_quals}"
            )
            continue
        if key in ("gene", "mRNA"):
            new_loc = f"1..{mrna_len}"
            new_feature_blocks.append(
                f"{_FEATURE_KEY_COL}{key:<16}{new_loc}\n{quals}"
            )
            continue

        mrna_intervals = gene_to_mrna_intervals(intervals, exons)
        if not mrna_intervals:
            # Feature lives entirely in an intron (or off the mRNA).
            continue
        new_loc = format_location_with_complement(mrna_intervals, complement)
        new_feature_blocks.append(
            f"{_FEATURE_KEY_COL}{key:<16}{new_loc}\n{quals}"
        )

    # Header LOCUS line update.
    header_lines = header.splitlines(keepends=True)
    for i, ln in enumerate(header_lines):
        if ln.startswith("LOCUS"):
            header_lines[i] = rewrite_locus_line(ln.rstrip("\n"), mrna_len) + "\n"
            break
    new_header = "".join(header_lines)

    return new_header + "".join(new_feature_blocks) + format_origin_block(mrna_seq)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gb-dir", required=True, type=Path,
                    help="Directory of source gene-structure .gb files")
    ap.add_argument("--filtered-high-dir", required=True, type=Path,
                    help="Directory of high-tier filtered score CSVs (e.g. 0.7)")
    ap.add_argument("--filtered-mod-dir", required=True, type=Path,
                    help="Directory of moderate-tier filtered score CSVs (e.g. 0.5)")
    ap.add_argument("--casoff-output-dir", required=True, type=Path,
                    help="Cas-OFFinder output directory (recursively scanned for *_output.txt)")
    ap.add_argument("--output-dir", required=True, type=Path,
                    help="Root output directory (07_GenBank_with_gRNA_Annotations)")
    ap.add_argument("--high-threshold", type=float, default=0.7)
    ap.add_argument("--mod-threshold", type=float, default=0.5)
    ap.add_argument("--max-mm", type=int, default=3,
                    help="Cap off-target mismatches at this many (default 3)")
    args = ap.parse_args()

    mm_table = load_casoffinder_mm_counts(args.casoff_output_dir, args.max_mm)
    print(f"[INFO] loaded mm counts for {len(mm_table)} unique sgRNA sequences "
          f"(<= {args.max_mm} mm)")

    high_label = f"{args.high_threshold:g}"
    mod_label  = f"{args.mod_threshold:g}"
    high_out = args.output_dir / f"01_Score_{high_label}_above"
    mod_out  = args.output_dir / f"02_Score_{mod_label}_above_incl_{high_label}"
    high_out.mkdir(parents=True, exist_ok=True)
    mod_out.mkdir(parents=True, exist_ok=True)

    high_csvs = {gene_id_from_csv_name(p.name): p
                 for p in sorted(args.filtered_high_dir.glob("*_filtered_*.csv"))
                 if not p.name.startswith("CRISPR_results_combined")}
    mod_csvs  = {gene_id_from_csv_name(p.name): p
                 for p in sorted(args.filtered_mod_dir.glob("*_filtered_*.csv"))
                 if not p.name.startswith("CRISPR_results_combined")}

    gb_files = sorted(args.gb_dir.glob("*.gb"))
    n_genes_processed = 0
    for gb_path in gb_files:
        gene_id = gb_path.stem
        try:
            gb_text = gb_path.read_text()
        except OSError as exc:
            print(f"[WARN] cannot read {gb_path}: {exc}", file=sys.stderr)
            continue

        # ---- High-only variant ---------------------------------------------
        if gene_id in high_csvs:
            rows = load_filtered_csv(high_csvs[gene_id])
            annotated, added, skipped, n_ss = build_annotated_gb(
                gb_text, rows, mm_table, args.max_mm,
                high_threshold=None, gene_id_hint=gene_id,
            )
            out = high_out / f"{gene_id}_score_ge{high_label}_mm{args.max_mm}.gb"
            out.write_text(annotated)
            print(f"[OK] {gene_id} >= {high_label}: {added} gRNA features + "
                  f"{n_ss} SS features (skipped {skipped}) -> {out}")
            if gene_id in MRNA_VIEW_TRANSCRIPTS:
                mrna_text = build_mrna_view(annotated)
                if mrna_text is not None:
                    mrna_out = out.with_name(out.stem + "_mRNA.gb")
                    mrna_out.write_text(mrna_text)
                    print(f"[OK]   mRNA view -> {mrna_out}")
                else:
                    print(f"[WARN] could not build mRNA view for {gene_id} "
                          f"(no exon/mRNA feature); skipped mRNA companion",
                          file=sys.stderr)
        else:
            print(f"[INFO] {gene_id}: no high-tier CSV (>= {high_label}); skipped")

        # ---- Combined variant (mod + high tag) ----------------------------
        if gene_id in mod_csvs:
            rows = load_filtered_csv(mod_csvs[gene_id])
            annotated, added, skipped, n_ss = build_annotated_gb(
                gb_text, rows, mm_table, args.max_mm,
                high_threshold=args.high_threshold,
                gene_id_hint=gene_id,
            )
            out = mod_out / f"{gene_id}_score_ge{mod_label}_incl_ge{high_label}_mm{args.max_mm}.gb"
            out.write_text(annotated)
            print(f"[OK] {gene_id} >= {mod_label} (HIGH tagged >= {high_label}): "
                  f"{added} gRNA features + {n_ss} SS features "
                  f"(skipped {skipped}) -> {out}")
            if gene_id in MRNA_VIEW_TRANSCRIPTS:
                mrna_text = build_mrna_view(annotated)
                if mrna_text is not None:
                    mrna_out = out.with_name(out.stem + "_mRNA.gb")
                    mrna_out.write_text(mrna_text)
                    print(f"[OK]   mRNA view -> {mrna_out}")
                else:
                    print(f"[WARN] could not build mRNA view for {gene_id} "
                          f"(no exon/mRNA feature); skipped mRNA companion",
                          file=sys.stderr)
        else:
            print(f"[INFO] {gene_id}: no moderate-tier CSV (>= {mod_label}); skipped")

        n_genes_processed += 1

    print(f"[DONE] processed {n_genes_processed} gene GenBank files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
