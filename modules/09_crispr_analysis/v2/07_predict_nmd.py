#!/usr/bin/env python3
"""
Module: 07_predict_nmd.py
Stage [7] — NMD prediction (rule-based heuristics, plant-aware).

Applies two complementary heuristics — the first canonical (mammalian), the
second plant-specific — and labels a guide as NMD-susceptible when EITHER
rule fires:

  (a) 50-nucleotide rule (Maquat 2004) — mammalian canon:
      A premature termination codon (PTC) is susceptible to NMD if it lies
      ≥ `--ptc-distance-threshold` upstream of the last exon-exon junction.

  (b) Long-3'UTR rule (Kerényi et al. 2008) — plant-specific, EJC-independent:
      A plant transcript whose (induced + original) 3' UTR is
      ≥ `--long-3utr-threshold` nt is NMD-susceptible even without a
      downstream exon-exon junction. The induced 3' UTR length is
      approximated as the CDS distance from the PTC to the natural stop
      codon plus the original 3' UTR length parsed from the GTF (zero if
      no three_prime_utr / UTR features are annotated — conservative
      lower bound, so plant NMD hits are undercounted, never overcounted,
      in the common plant-GTF case where only CDS rows are present).

For each guide's top indel outcomes, reads the PTC position (stage [5]) and
classifies each indel as:
  - nmd_predicted : either rule (a) or rule (b) triggers
  - nmd_escape    : PTC in last exon AND induced 3' UTR < long_3utr_threshold
  - no_ptc        : no in-frame stop codon introduced

Adds columns: indel{N}_nmd_class, indel{N}_induced_3utr_len, nmd_summary
(worst-case across indels).

Usage:
    python3 07_predict_nmd.py \\
        --input                <protein.tsv>   \\
        --outdir               <output_dir>    \\
        --gtf                  <annotation.gtf> \\
        --gene-group           DMP             \\
        --ptc-distance-threshold 50 \\
        --long-3utr-threshold    350

References:
    Maquat LE. 2004. Nat Rev Mol Cell Biol 5:89-99. doi:10.1038/nrm1310
    Kerényi Z et al. 2008. EMBO J 27(11):1585-1595.
        doi:10.1038/emboj.2008.189
    Shaul O. 2015. Trends Plant Sci 20(12):767-779.
        doi:10.1016/j.tplants.2015.08.011
"""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="NMD prediction via 50-nt heuristic")
    p.add_argument("--input",                    required=True)
    p.add_argument("--outdir",                   required=True)
    p.add_argument("--gtf",                      default="",
                   help="GTF annotation for EJC position calculation (optional; defaults to no_ptc if absent)")
    p.add_argument("--gene-group",               default="DMP",
                   help="Family tag used as a fallback regex when --gene-id does not match anything")
    p.add_argument("--gene-id",                  default="",
                   help="Specific target gene id (e.g. SMEL5_01g008730). Matched against GTF gene_id "
                        "with prefix semantics; falls back to --gene-group regex if no hit.")
    p.add_argument("--ptc-distance-threshold",   type=int, default=50,
                   help="Mammalian 50-nt EJC rule: PTC ≥ this many nt upstream of "
                        "the last EJC triggers NMD (Maquat 2004).")
    p.add_argument("--long-3utr-threshold",      type=int, default=350,
                   help="Plant long-3'UTR rule: induced 3' UTR ≥ this many nt "
                        "triggers NMD even without a downstream EJC "
                        "(Kerényi et al. 2008). Set 0 to disable.")
    p.add_argument("--overwrite",                action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


# ─── GTF: last EJC position per transcript ───────────────────────────────────

def last_ejc_position(gtf_path: str, gene_group: str, gene_id: str = "") -> dict:
    """
    Return {transcript_id: (last_ejc_nt_pos_in_CDS, total_exon_count,
                            total_cds_len, original_3utr_len)}.

        last_ejc        : position (in spliced-CDS coordinates) of the last
                          exon-exon junction, approximated as total_cds_len
                          minus last_exon_len.
        total_exon_count: number of CDS-containing exons for that transcript.
        total_cds_len   : cumulative CDS length (used by the plant long-3'UTR
                          rule to compute induced 3' UTR length).
        original_3utr_len: sum of three_prime_utr / UTR features in the GTF
                          for this transcript. 0 when the GTF has no UTR rows
                          (common in SMEL5/GPE001970 annotations).

    Matching priority (same contract as load_cds_coords in stage 5):
      1. Exact / prefix match on gene_id (e.g. "SMEL5_01g008730")
      2. Fallback regex on gene_group when gene_id is empty OR no hit
    Returns empty dict if GTF not found.
    """
    ejc: dict = {}
    fallback_pattern = re.compile(gene_group, re.IGNORECASE) if gene_group else None
    gene_id_clean = (gene_id or "").strip()
    exons_by_tx: dict = {}
    utr3_by_tx: dict = {}   # {transcript_id: summed_3utr_length}
    matched_by_id = 0
    matched_by_fallback = 0

    def _matches(gid: str) -> bool:
        if gene_id_clean:
            return (gid == gene_id_clean
                    or gid.startswith(gene_id_clean + ".")
                    or gid.startswith(gene_id_clean + "_"))
        if fallback_pattern is not None:
            return bool(fallback_pattern.search(gid))
        return False

    def _parse_attrs(attrs: str) -> tuple[str, str]:
        """Extract (gene_id, transcript_id) from GTF or GFF3 attributes.

        Mirrors the parser in 05_rebuild_transcripts.py so both stages agree
        on the same transcript id for any given CDS row.
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

    def _scan(use_fallback_only: bool):
        nonlocal matched_by_id, matched_by_fallback
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                feat = parts[2]
                # Capture CDS for exon counting / cumulative length, plus any
                # three_prime_utr variant the annotator may have used.
                if feat not in ("CDS", "three_prime_utr", "three_prime_UTR",
                                "3UTR", "UTR3"):
                    continue
                gid, tid = _parse_attrs(parts[8])

                if use_fallback_only:
                    if fallback_pattern is None or not fallback_pattern.search(gid):
                        continue
                    matched_by_fallback += 1
                else:
                    if not _matches(gid):
                        continue
                    if gene_id_clean:
                        matched_by_id += 1
                    else:
                        matched_by_fallback += 1

                tid = tid or gid
                length = int(parts[4]) - int(parts[3]) + 1
                if feat == "CDS":
                    exons_by_tx.setdefault(tid, []).append((int(parts[3]), length))
                else:
                    utr3_by_tx[tid] = utr3_by_tx.get(tid, 0) + length

    try:
        _scan(use_fallback_only=False)
    except FileNotFoundError:
        _log(f"[07_nmd] GTF not found: {gtf_path}", level="WARN")
        return ejc

    # Second-chance: gene_id was given but missed the GTF — retry with family regex.
    if gene_id_clean and not exons_by_tx and fallback_pattern is not None:
        _log(f"[07_nmd] No GTF CDS matched gene_id={gene_id_clean!r}; "
             f"retrying with family regex {gene_group!r}.", level="WARN")
        _scan(use_fallback_only=True)

    for tid, exon_list in exons_by_tx.items():
        sorted_exons  = sorted(exon_list, key=lambda x: x[0])
        total_cds     = sum(e[1] for e in sorted_exons)
        last_exon_len = sorted_exons[-1][1]
        n_exons       = len(sorted_exons)
        # EJC is deposited ~20-24 nt upstream of each junction; use junction position
        # approximation: last_ejc = total_cds - last_exon_len
        ejc[tid] = (
            max(0, total_cds - last_exon_len),
            n_exons,
            total_cds,
            utr3_by_tx.get(tid, 0),   # 0 when GTF has no UTR rows (common in plants)
        )

    return ejc


# ─── NMD classification ───────────────────────────────────────────────────────

def induced_3utr_length(ptc_pos: int, total_cds: int, original_3utr: int) -> int:
    """Length of the 3' UTR that results after an indel introduces a PTC.

    The induced 3' UTR spans from the PTC stop codon (exclusive) to the end
    of the CDS plus any original 3' UTR that remains annealed downstream.
    A PTC at or past the natural stop collapses the CDS contribution to 0;
    the original 3' UTR length is always added (or 0 when the GTF lacks
    three_prime_utr rows, as is typical for plant annotations — in which
    case this is a conservative lower bound).
    """
    remaining_cds = max(0, total_cds - ptc_pos - 3)  # -3 to exclude the stop itself
    return remaining_cds + max(0, original_3utr)


def classify_nmd(ptc_pos_str: str, ptc_exon_str: str,
                 last_ejc: int, total_exons: int,
                 total_cds: int, original_3utr: int,
                 ejc_threshold: int, long_3utr_threshold: int) -> tuple:
    """
    Returns (class, induced_3utr_len).

    class ∈ {'nmd_predicted', 'nmd_escape', 'no_ptc'}.

    A guide is 'nmd_predicted' if EITHER
       (a) PTC ≥ ejc_threshold nt upstream of last EJC (Maquat 2004), OR
       (b) induced 3' UTR ≥ long_3utr_threshold nt (Kerényi 2008 plant rule;
           skipped when long_3utr_threshold <= 0).
    'nmd_escape' requires PTC in the last exon AND a short induced 3' UTR.
    """
    if ptc_pos_str in ("NA", "", None):
        return ("no_ptc", 0)

    try:
        ptc_pos = int(ptc_pos_str)
    except ValueError:
        return ("no_ptc", 0)

    induced = induced_3utr_length(ptc_pos, total_cds, original_3utr) if total_cds else 0

    # Rule (b) — plant long-3'UTR trigger. Checked first because it fires
    # even when the PTC is in the last exon (EJC-independent).
    if long_3utr_threshold > 0 and induced >= long_3utr_threshold:
        return ("nmd_predicted", induced)

    # Rule (a) — mammalian 50-nt EJC rule.
    # PTC in the last exon escapes (no downstream EJC to license NMD)
    # UNLESS the plant long-3'UTR rule already fired above.
    if ptc_exon_str and "exon" in ptc_exon_str.lower() and total_exons > 0:
        m = re.search(r"exon(\d+)", ptc_exon_str, re.IGNORECASE)
        if m and int(m.group(1)) >= total_exons:
            return ("nmd_escape", induced)

    distance_from_last_ejc = last_ejc - ptc_pos
    if distance_from_last_ejc >= ejc_threshold:
        return ("nmd_predicted", induced)
    return ("nmd_escape", induced)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    inpath  = Path(args.input)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.nmd.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[07_nmd] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    # Derive the specific target gene id from the input filename stem when
    # --gene-id was not explicitly supplied (e.g. SMEL5_01g008730.nmd.tsv).
    gene_id_effective = args.gene_id or gene_stem

    if not args.gtf or not Path(args.gtf).exists():
        _log(f"[07_nmd] GTF not provided or not found ({args.gtf!r}) — "
              "all guides will be classified as no_ptc.", level="WARN")
        ejc_map = {}
    else:
        ejc_map = last_ejc_position(args.gtf, args.gene_group, gene_id_effective)
    # Single representative EJC value from the first transcript in the GTF.
    # This matches the single-transcript assumption in stage [5] (05_rebuild_transcripts.py
    # also uses next(iter(cds_coords))). Both stages must use the same transcript
    # for PTC placement and NMD classification to be consistent. If the GTF
    # contains multiple isoforms, the first encountered depends on GTF line order.
    # Tuple: (last_ejc, total_exons, total_cds, original_3utr_len)
    _default_tuple = next(iter(ejc_map.values()), (0, 0, 0, 0))
    default_ejc    = _default_tuple[0]
    default_exons  = _default_tuple[1]
    default_cds    = _default_tuple[2] if len(_default_tuple) >= 3 else 0
    default_utr3   = _default_tuple[3] if len(_default_tuple) >= 4 else 0

    if default_cds == 0 and args.long_3utr_threshold > 0:
        _log("[07_nmd] Plant long-3'UTR rule ENABLED but CDS length unknown "
             "(no GTF CDS rows matched) — rule (b) will be skipped and only "
             "the 50-nt EJC rule will fire.", level="WARN")
    elif default_utr3 == 0 and args.long_3utr_threshold > 0:
        _log(f"[07_nmd] Plant long-3'UTR rule: using conservative induced-3'UTR "
             f"proxy (CDS {default_cds} nt; no three_prime_utr rows in GTF — "
             f"original 3' UTR treated as 0 nt). Kerényi 2008 threshold = "
             f"{args.long_3utr_threshold} nt.", level="INFO")

    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    n_indels = sum(1 for f in fields if f.startswith("indel") and f.endswith("_ptc_pos"))
    new_cols = []
    for i in range(1, n_indels + 1):
        new_cols += [f"indel{i}_nmd_class", f"indel{i}_induced_3utr_len"]
    new_cols += ["nmd_summary", "nmd_trigger"]
    new_fields = fields + new_cols

    for row in rows:
        nmd_classes  = []
        nmd_triggers = []
        for i in range(1, n_indels + 1):
            ptc_pos  = row.get(f"indel{i}_ptc_pos",  "NA")
            ptc_exon = row.get(f"indel{i}_ptc_exon", "NA")
            nmd_cls, induced = classify_nmd(
                ptc_pos, ptc_exon,
                default_ejc, default_exons,
                default_cds, default_utr3,
                args.ptc_distance_threshold,
                args.long_3utr_threshold,
            )
            row[f"indel{i}_nmd_class"]        = nmd_cls
            row[f"indel{i}_induced_3utr_len"] = str(induced)
            nmd_classes.append(nmd_cls)
            # Record WHICH rule triggered so the thesis can report both signals.
            if nmd_cls == "nmd_predicted":
                if args.long_3utr_threshold > 0 and induced >= args.long_3utr_threshold:
                    nmd_triggers.append("long_3utr")
                else:
                    nmd_triggers.append("ejc_50nt")

        # Summary: worst-case across all indel outcomes
        if "nmd_predicted" in nmd_classes:
            row["nmd_summary"] = "nmd_predicted"
        elif "nmd_escape" in nmd_classes:
            row["nmd_summary"] = "nmd_escape"
        else:
            row["nmd_summary"] = "no_ptc"
        # Trigger summary: pipe-joined unique rule names (e.g. "long_3utr|ejc_50nt")
        row["nmd_trigger"] = "|".join(sorted(set(nmd_triggers))) if nmd_triggers else ""

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    nmd_n    = sum(1 for r in rows if r.get("nmd_summary") == "nmd_predicted")
    escape_n = sum(1 for r in rows if r.get("nmd_summary") == "nmd_escape")
    long_3utr_n = sum(1 for r in rows if "long_3utr" in (r.get("nmd_trigger") or ""))
    ejc_n       = sum(1 for r in rows if "ejc_50nt"  in (r.get("nmd_trigger") or ""))
    _log(f"[07_nmd] {len(rows)} guides; "
          f"{nmd_n} nmd_predicted ({ejc_n} via EJC-50nt, {long_3utr_n} via long-3'UTR), "
          f"{escape_n} nmd_escape -> {outpath}", level="INFO")


if __name__ == "__main__":
    main()
