#!/usr/bin/env python3
"""
generate_variants.py -- build a single Stage-14 variant FASTA from a source
sequence. Two modes:

  deletion (default, used by HAP2 + paper-anchored DMP truncations):
      Apply one or more deletion ranges to a PROTEIN source sequence and
      write the truncated protein. Coordinates are 1-indexed inclusive,
      referring to the original sequence.

  frameshift (used by guide-anchored DMP NHEJ models):
      Locate a Cas9 protospacer + PAM in a TRANSCRIPT (DNA), apply a +1 bp
      insertion at the predicted DSB (3 bp 5' of PAM on the protospacer),
      then re-translate the CDS in the shifted frame, stopping at the first
      new in-frame stop codon. The output is a Met-initiated protein that
      matches WT up to the cut and ends with the NHEJ-shifted tail.

Called by 14_Interaction_Domain_Mapping.sh for each row in [hap2_variants]
and [dmp_variants]; the orchestrator picks the mode based on which TOML
section the variant is declared under (coords vs. guides).
"""

import argparse
import re
import sys


# --------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------

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


def read_named_fasta(path, record_id):
    """Return (header_without_gt, sequence) for the FASTA record whose header
    contains record_id (substring match). Falls back to the first record if
    record_id is empty / None.
    """
    if not record_id:
        return read_first_fasta(path)
    header = None
    parts = []
    capture = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if capture:
                    break
                if record_id in line:
                    capture = True
                    header = line[1:]
                continue
            if capture:
                parts.append(line)
    if header is None:
        raise ValueError(f"Record '{record_id}' not found in {path}")
    return header, "".join(parts)


def write_fasta(path, header, seq, line_width=80):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            fh.write(seq[i : i + line_width] + "\n")


# --------------------------------------------------------------------------
# Mode: deletion
# --------------------------------------------------------------------------

def parse_deletions(deletion_str):
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


def run_deletion_mode(args):
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


# --------------------------------------------------------------------------
# Mode: frameshift
# --------------------------------------------------------------------------

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

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def locate_cut(transcript, guide23, strand):
    """Locate the Cas9 cut site for a 23-mer guide (= 20 nt protospacer + 3 nt PAM).

    Returns (cut_0idx, binding_0idx). cut_0idx is the 0-indexed position
    BEFORE which a +1 insertion goes (= the 1-indexed position of the last
    preserved bp on the + strand). Cas9 cleaves 3 bp 5' of the PAM on the
    protospacer, so the blunt DSB sits between 23-mer indices 16 and 17 on
    the bound strand.
    """
    t = transcript.upper()
    g = guide23.upper()
    if strand == "+":
        idx = t.find(g)
        if idx < 0:
            raise ValueError(f"Guide '{guide23}' not found on + strand of transcript")
        return idx + 17, idx
    elif strand == "-":
        rc = revcomp(g)
        idx = t.find(rc)
        if idx < 0:
            raise ValueError(
                f"Guide '{guide23}' (RC '{rc}') not found on - strand of transcript"
            )
        # On + strand: PAM-RC (CCN) at idx..idx+2, protospacer at idx+3..idx+22.
        # Cut between t[idx+5] and t[idx+6] (mirror of + strand case via blunt DSB).
        return idx + 6, idx
    else:
        raise ValueError(f"strand must be '+' or '-': got '{strand}'")


def translate(dna):
    """Translate DNA in frame 0; stop at first stop codon. Returns (aa, hit_stop)."""
    aa = []
    n = len(dna)
    i = 0
    while i + 2 < n:
        codon = dna[i:i + 3].upper()
        i += 3
        res = CODON_TABLE.get(codon, "X")
        if res == "*":
            return "".join(aa), True
        aa.append(res)
    return "".join(aa), False


def parse_cds_from_header(header):
    """Return (cds_start, cds_end) parsed from "CDS=N-M" in a FASTA header, or (None, None)."""
    m = re.search(r"CDS=(\d+)-(\d+)", header)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def run_frameshift_mode(args):
    if not args.dna_input:
        raise SystemExit("--dna-input is required in frameshift mode")
    if not args.guide_23mer:
        raise SystemExit("--guide-23mer is required in frameshift mode")
    if args.guide_strand not in ("+", "-"):
        raise SystemExit("--guide-strand must be '+' or '-'")

    header, transcript = read_named_fasta(args.dna_input, args.dna_record)
    cds_start = args.cds_start
    cds_end = args.cds_end
    if cds_start is None or cds_end is None:
        h_start, h_end = parse_cds_from_header(header)
        if cds_start is None:
            cds_start = h_start
        if cds_end is None:
            cds_end = h_end
    if cds_start is None or cds_end is None:
        raise SystemExit(
            "CDS coords not provided and not parseable from header "
            "(expected 'CDS=N-M'); pass --cds-start / --cds-end."
        )

    cut_0idx, binding_0idx = locate_cut(transcript, args.guide_23mer, args.guide_strand)
    cut_t_1idx = cut_0idx  # last preserved 1-indexed bp on + strand
    cds_pos = cut_t_1idx - cds_start + 1
    insert_nt = args.insert_nt or "A"

    # Apply +1 insertion at the cut site on the transcript, then extract CDS.
    # If the cut precedes the CDS start the CDS itself is unchanged (frame
    # preserved); otherwise the indel sits inside the CDS and shifts the
    # downstream reading frame.
    mod_transcript = transcript[:cut_0idx] + insert_nt + transcript[cut_0idx:]
    if cut_0idx < cds_start - 1:
        # Cut in 5'UTR; CDS slides by 1 in transcript coords but its
        # internal frame is preserved -> identical to WT translation.
        new_cds_start = cds_start + 1
        scope_note = "cut in 5'UTR (CDS frame preserved); translation == WT"
    elif cut_0idx >= cds_end:
        new_cds_start = cds_start
        scope_note = "cut in 3'UTR (CDS frame preserved); translation == WT"
    else:
        new_cds_start = cds_start
        scope_note = f"cut inside CDS at cds_pos {cds_pos} (codon {(cds_pos - 1) // 3 + 1})"
    mod_cds = mod_transcript[new_cds_start - 1:]
    fs_protein, hit_stop = translate(mod_cds)
    if not fs_protein:
        raise SystemExit(
            "Frameshift translation yielded empty protein; "
            f"start codon may be disrupted. Note: {scope_note}"
        )
    if fs_protein[0] != "M":
        # Should not happen if CDS is intact upstream; flag for visibility.
        print(
            f"  WARNING: first codon of frameshift product is "
            f"'{fs_protein[0]}' (expected 'M'); start may be disrupted.",
            file=sys.stderr,
        )

    # Compose header. Carry guide + cut metadata so AF3 / downstream
    # consumers can trace which CRISPR outcome this peptide represents.
    header_parts = [args.name]
    if args.description:
        header_parts.append(args.description)
    header_parts.append(f"len={len(fs_protein)}")
    header_parts.append(
        f"frameshift=guide:{args.guide_23mer},strand:{args.guide_strand},"
        f"insert:+1{insert_nt},cds_pos:{cds_pos},hit_stop:{hit_stop}"
    )
    write_fasta(args.output, " | ".join(header_parts), fs_protein)

    print(
        f"  {args.name}: NHEJ +1{insert_nt} at transcript {cut_t_1idx}, "
        f"{scope_note}, protein {len(fs_protein)} aa "
        f"({'hit' if hit_stop else 'no'} stop codon)",
        file=sys.stderr,
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Build a Stage-14 variant FASTA (deletion or frameshift)."
    )
    ap.add_argument("--mode", choices=("deletion", "frameshift"), default="deletion",
                    help="Variant construction mode (default: deletion)")
    ap.add_argument("--name", required=True,
                    help="Variant identifier written into the FASTA header")
    ap.add_argument("--description", default="",
                    help="Human-readable variant description for the header")
    ap.add_argument("--output", required=True, help="Output FASTA path")

    # deletion-mode inputs
    ap.add_argument("--input",
                    help="Source PROTEIN FASTA (deletion mode; first record used)")
    ap.add_argument("--deletions", default="",
                    help='Comma-separated "start-end" ranges to delete '
                         "(deletion mode; 1-indexed, inclusive; empty = WT)")

    # frameshift-mode inputs
    ap.add_argument("--dna-input",
                    help="Source TRANSCRIPT FASTA (frameshift mode; DNA)")
    ap.add_argument("--dna-record", default="",
                    help="Record ID inside --dna-input (substring match; "
                         "default = first record)")
    ap.add_argument("--cds-start", type=int, default=None,
                    help="1-indexed CDS start in the transcript (default: parse "
                         "'CDS=N-M' from the record header)")
    ap.add_argument("--cds-end", type=int, default=None,
                    help="1-indexed CDS end in the transcript (default: parse "
                         "'CDS=N-M' from the record header)")
    ap.add_argument("--guide-23mer", default="",
                    help="Cas9 guide sequence with PAM (23 nt; 20 nt protospacer "
                         "+ 3 nt NGG PAM)")
    ap.add_argument("--guide-strand", default="",
                    help="Guide strand: '+' or '-'")
    ap.add_argument("--insert-nt", default="A",
                    help="Nucleotide inserted at the DSB to model +1 NHEJ (default: A)")

    args = ap.parse_args()

    if args.mode == "deletion":
        if not args.input:
            ap.error("--input is required in deletion mode")
        run_deletion_mode(args)
    else:
        run_frameshift_mode(args)


if __name__ == "__main__":
    main()
