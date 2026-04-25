#!/usr/bin/env python3
"""
Module: 06_protein_consequence.py
Stage [6] — Protein-level consequence assessment.

For each guide's mutant CDS FASTA (stage [5]):
  1. Translates the mutant CDS to protein using Biopython.
  2. Checks whether the predicted truncation falls within a known Pfam domain
     (using a pre-computed domain annotation TSV).
  3. Flags whether an ESMFold/AF3 structure prediction is recommended.

Outputs a TSV with columns added to the guide table:
  mutant_protein_len, domain_hit, domain_name, structure_flag

Usage:
    python3 06_protein_consequence.py \\
        --input         <transcripts.tsv>   \\
        --outdir        <output_dir>        \\
        --domain-tsv    <pfam_domains.tsv>  \\
        --flag-domain-hits true
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Protein-level consequence assessment")
    p.add_argument("--input",             required=True)
    p.add_argument("--outdir",            required=True)
    p.add_argument("--domain-tsv",        default="",
                   help="TSV: gene_id<tab>domain_name<tab>aa_start<tab>aa_end")
    p.add_argument("--flag-domain-hits",  default="true")
    p.add_argument("--overwrite",         action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


# ─── Domain loading ───────────────────────────────────────────────────────────

def load_domains(tsv_path: str) -> dict:
    """Return {gene_id: [(domain_name, aa_start, aa_end), ...]}."""
    domains: dict = {}
    if not tsv_path or not Path(tsv_path).exists():
        return domains
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = row.get("gene_id") or row.get("gene") or ""
            dm  = row.get("domain_name") or row.get("domain") or ""
            try:
                s = int(row.get("aa_start", 0))
                e = int(row.get("aa_end", 0))
            except ValueError:
                s, e = 0, 0
            domains.setdefault(gid, []).append((dm, s, e))
    return domains


# ─── Translation ──────────────────────────────────────────────────────────────

def translate_cds(fasta_path: str) -> tuple[str, int]:
    """
    Translate the first sequence in fasta_path.
    Returns (protein_sequence, original_nt_length).
    Stops at first stop codon.
    """
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq

        rec = next(SeqIO.parse(fasta_path, "fasta"), None)
        if rec is None:
            return "", 0
        nt = str(rec.seq).upper()
        # Trim to nearest codon
        trim = len(nt) - (len(nt) % 3)
        protein = str(Seq(nt[:trim]).translate(to_stop=True))
        return protein, len(nt)
    except ImportError:
        _log("[06_protein] Biopython not installed — translation skipped.", level="WARN")
        return "", 0
    except Exception as exc:
        _log(f"[06_protein] Translation error for {fasta_path}: {exc}", level="ERROR")
        return "", 0


# ─── Domain intersection ──────────────────────────────────────────────────────

def check_domain_hit(protein_len: int, domains: list) -> tuple[bool, str]:
    """
    Return (hit_bool, domain_name_str) if the truncation point (protein_len)
    falls within any annotated domain.
    """
    hits = []
    for dm_name, dm_start, dm_end in domains:
        if dm_start <= protein_len <= dm_end:
            hits.append(dm_name)
    if hits:
        return True, "|".join(hits)
    return False, ""


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    inpath   = Path(args.input)
    outdir   = Path(args.outdir)
    # Stage [5] stores FASTA paths relative to its own outdir (inpath.parent).
    fasta_base = inpath.parent
    outdir.mkdir(parents=True, exist_ok=True)
    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.protein.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[06_protein] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    flag_domains = args.flag_domain_hits.lower() not in ("false", "0", "no")
    domains_db   = load_domains(args.domain_tsv)

    # Surface the most common cause of silent "0 flagged for structure review"
    # runs: the domain TSV was left unset, does not exist, or contains no
    # parseable rows. Without this log the pipeline appeared healthy (stage
    # [6b] reported ok=0 skipped=N error=0) while actually running with no
    # domain annotation whatsoever.
    if flag_domains:
        if not args.domain_tsv:
            _log("[06_protein] No --domain-tsv provided; falling back to "
                 "frameshift+early-truncation heuristic for structure_flag. "
                 "Provide a Pfam TSV (gene_id\\tdomain_name\\taa_start\\taa_end) "
                 "for domain-aware flagging.", level="WARN")
        elif not Path(args.domain_tsv).exists():
            _log(f"[06_protein] Domain TSV not found: {args.domain_tsv} — "
                 "falling back to frameshift+early-truncation heuristic for "
                 "structure_flag. Generate with modules/utils/"
                 "domtbl_to_domains_tsv.py from the stage-[01] HMMER domtbl.",
                 level="WARN")
        elif not domains_db:
            _log(f"[06_protein] Domain TSV {args.domain_tsv} loaded 0 rows "
                 "— falling back to frameshift+early-truncation heuristic.",
                 level="WARN")

    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    # Determine how many indels were annotated in stage [5]
    n_indels = sum(1 for f in fields if f.startswith("indel") and f.endswith("_seq"))

    new_cols: list[str] = []
    for i in range(1, n_indels + 1):
        new_cols += [f"indel{i}_protein_len", f"indel{i}_domain_hit", f"indel{i}_domain_name"]
    new_cols += ["any_domain_hit", "structure_flag"]
    new_fields = fields + new_cols

    # Fallback heuristic: when no Pfam domain data is available for a gene
    # the only usable signal from stage [5] is (a) whether the top indel
    # produces a frameshift and (b) how early the premature termination
    # codon lands. Early frameshifts almost always disrupt function, so
    # flag them for structural review when we can't consult Pfam.
    _early_truncation_frac = 0.9  # PTC before 90% of the full-length protein

    for row in rows:
        # Determine gene-group key for domain lookup. Prefer an explicit
        # gene_id column on the row; fall back to the file-derived
        # gene_stem (e.g. SMEL5_01g008730) which matches the id schema
        # emitted by modules/utils/domtbl_to_domains_tsv.py.
        gid_key = ""
        for col in ("gene_id", "gene", "guideId"):
            if col in row and row[col]:
                gid_key = row[col].strip()
                break
        if not gid_key or gid_key not in domains_db:
            gid_key = gene_stem
        gene_domains = domains_db.get(gid_key, [])

        # Capture reference protein length from the longest indel translation
        # on this row so the fallback heuristic can compare PTC position to
        # something meaningful. Computed inside the per-indel loop below.
        ref_protein_len = 0

        any_hit = False
        for i in range(1, n_indels + 1):
            fa_rel = row.get(f"indel{i}_seq", "NA")
            fa_path = str(fasta_base / fa_rel) if (fa_rel and fa_rel != "NA") else "NA"
            if fa_path and fa_path != "NA" and Path(fa_path).exists():
                protein, nt_len = translate_cds(fa_path)
                p_len = len(protein)
                ref_protein_len = max(ref_protein_len, nt_len // 3)
                if flag_domains and gene_domains:
                    hit, dm_name = check_domain_hit(p_len, gene_domains)
                else:
                    hit, dm_name = False, ""
                if hit:
                    any_hit = True
                row[f"indel{i}_protein_len"]  = str(p_len)
                row[f"indel{i}_domain_hit"]   = str(hit)
                row[f"indel{i}_domain_name"]  = dm_name
            else:
                row[f"indel{i}_protein_len"]  = "NA"
                row[f"indel{i}_domain_hit"]   = "NA"
                row[f"indel{i}_domain_name"]  = "NA"

        row["any_domain_hit"]  = str(any_hit)
        # Prefer the domain-based flag when Pfam data exists for this gene.
        # Otherwise consult the frameshift+early-truncation fallback using
        # stage [5] columns that are always present.
        if any_hit:
            row["structure_flag"] = "recommend_structure"
        elif flag_domains and not gene_domains:
            frameshifted = any(
                str(row.get(f"indel{i}_is_frameshift", "")).lower() == "true"
                for i in range(1, n_indels + 1)
            )
            early_ptc = False
            if ref_protein_len > 0:
                for i in range(1, n_indels + 1):
                    ptc_str = row.get(f"indel{i}_ptc_pos", "NA")
                    if ptc_str in (None, "", "NA", "N/A"):
                        continue
                    try:
                        ptc_nt = int(ptc_str)
                    except ValueError:
                        continue
                    # ptc_pos is a 0-based nucleotide offset on the mutant
                    # CDS; divide by 3 to get amino-acid position, then
                    # compare to ref_protein_len.
                    if (ptc_nt / 3.0) < _early_truncation_frac * ref_protein_len:
                        early_ptc = True
                        break
            row["structure_flag"] = ("recommend_structure"
                                     if frameshifted and early_ptc else "")
        else:
            row["structure_flag"] = ""

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    struct_n = sum(1 for r in rows if r.get("structure_flag"))
    _log(f"[06_protein] {len(rows)} guides assessed; "
          f"{struct_n} flagged for structure review -> {outpath}", level="INFO")


if __name__ == "__main__":
    main()
