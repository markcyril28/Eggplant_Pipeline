#!/usr/bin/env python3
"""
Build per-paralog MSA input FASTAs from BLAST lollipop groups.

For each anchor SmelDMP paralog, the curated BLAST CSVs are filtered to keep
only hits with Bit Score >= --bitscore. The matched Query IDs are then
resolved against the union of [ortholog_blast].query_fastas and
query_protein_fastas, and written alongside the anchor sequence to:

    <out_dir>/<paralog_id>_NUCLEOTIDE_Sequence.fasta   (from blastn CSV)
    <out_dir>/<paralog_id>_AMINO_ACID_Sequence.fasta   (from blastp CSV)

Usage:
    python3 build_v4_blast_groups.py \\
        --pipeline-dir <root> \\
        --gene-group DMP \\
        --bitscore 200 \\
        --out-dir <III_RESULT/DMP/04_MSA/merged_input/Selected/v4_BLAST_Groups_bitscore200>

When invoked without --pipeline-dir, the script auto-detects it as the parent
of `modules/`.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    import tomllib  # py311+
except ModuleNotFoundError:
    import tomli as tomllib  # type: ignore

# Make modules/utils/ importable so we can prettify FASTA headers using the
# same naming logic the BLAST y-axis uses.
_UTILS = Path(__file__).resolve().parent.parent / "utils"
if str(_UTILS) not in sys.path:
    sys.path.insert(0, str(_UTILS))
from dmp_query_labels import short_label  # noqa: E402


_DEFAULT_PARALOGS = (
    "SMEL5_01g008730.1",
    "SMEL5_01g026030.1",
    "SMEL5_02g013320.1",
    "SMEL5_04g005390.1",
    "SMEL5_10g003660.1",
    "SMEL5_10g017610.1",
    "SMEL5_12g005350.1",
)


def parse_fasta(path: Path):
    header, chunks = None, []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def first_token(header: str) -> str:
    return header.split()[0] if header else header


def load_fasta_dict(paths: list[Path]) -> dict[str, tuple[str, str]]:
    """Read multiple FASTAs and return {first_token: (full_header, sequence)}."""
    out: dict[str, tuple[str, str]] = {}
    for p in paths:
        if not p.exists():
            print(f"  [warn] missing FASTA: {p}", file=sys.stderr)
            continue
        for header, seq in parse_fasta(p):
            key = first_token(header)
            if key not in out:
                out[key] = (header, seq)
    return out


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def prettify_header(header: str) -> str:
    """Replace the first whitespace token of a FASTA header with its
    short_label() form, dropping any trailing annotation. Matches the BLAST
    y-axis label convention so MSA outputs and BLAST figures share names."""
    if not header:
        return header
    first_token = header.split(None, 1)[0]
    return short_label(first_token)


def write_fasta(records: list[tuple[str, str]], out_path: Path) -> None:
    """Write FASTA with prettified headers. When the same short label maps to
    multiple sequences in this file (e.g. three Musa paralogs all collapse to
    'MaDMP8-like'), append _2, _3, ... so every record stays distinct."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    with out_path.open("w", encoding="utf-8", newline="\n") as fh:
        first = True
        for header, seq in records:
            new_id = prettify_header(header)
            counts[new_id] = counts.get(new_id, 0) + 1
            display = new_id if counts[new_id] == 1 else f"{new_id}_{counts[new_id]}"
            if not first:
                fh.write("\n")
            first = False
            fh.write(f">{display}\n{wrap(seq)}\n")


def filter_csv(csv_path: Path, paralog_ids: set[str], bitscore_min: float) -> dict[str, list[str]]:
    """Return {paralog_id: [Query ID, ...]} after filtering."""
    if not csv_path.exists():
        print(f"  [warn] CSV missing: {csv_path}", file=sys.stderr)
        return {pid: [] for pid in paralog_ids}
    out: dict[str, list[str]] = {pid: [] for pid in paralog_ids}
    seen_per_pid: dict[str, set[str]] = {pid: set() for pid in paralog_ids}
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid = row.get("Subject ID", "").strip()
            if sid not in paralog_ids:
                continue
            try:
                bs = float(row.get("Bit Score", "") or 0.0)
            except ValueError:
                continue
            if bs < bitscore_min:
                continue
            qid = row.get("Query ID", "").strip()
            if not qid or qid in seen_per_pid[sid]:
                continue
            seen_per_pid[sid].add(qid)
            out[sid].append(qid)
    return out


def autodiscover_csv(blast_dir: Path, kind: str) -> Path | None:
    """Find <blast_dir>/curated_results/merged_blast{n,p}_*_plant_only.csv,
    preferring the most recent date in the filename."""
    pat = re.compile(rf"^merged_blast{kind}_\d{{4}}-\d{{2}}-\d{{2}}_.*_plant_only\.csv$")
    candidates = sorted(
        (p for p in (blast_dir / "curated_results").glob(f"merged_blast{kind}_*_plant_only.csv") if pat.match(p.name)),
        reverse=True,
    )
    return candidates[0] if candidates else None


def safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


def build_for_type(
    csv_path: Path,
    anchor_fasta: Path,
    query_fastas: list[Path],
    paralog_ids: list[str],
    bitscore_min: float,
    out_dir: Path,
    type_suffix: str,
) -> None:
    print(f"\n=== {type_suffix} ===")
    print(f"  CSV         : {csv_path}")
    print(f"  Anchor FASTA: {anchor_fasta}")
    print(f"  Query FASTAs: {len(query_fastas)} files")

    anchor_dict = load_fasta_dict([anchor_fasta])
    query_dict = load_fasta_dict(query_fastas)
    hits_by_paralog = filter_csv(csv_path, set(paralog_ids), bitscore_min)

    n_total_written = 0
    for pid in paralog_ids:
        anchor = anchor_dict.get(pid)
        if anchor is None:
            print(f"  [warn] anchor sequence not found in {anchor_fasta.name}: {pid} — skipping {type_suffix} for this paralog", file=sys.stderr)
            continue
        records = [anchor]
        missed = []
        for qid in hits_by_paralog.get(pid, []):
            entry = query_dict.get(qid)
            if entry is None:
                missed.append(qid)
                continue
            records.append(entry)
        out_path = out_dir / f"{safe_filename(pid)}_{type_suffix}_Sequence.fasta"
        write_fasta(records, out_path)
        n_hits = len(hits_by_paralog.get(pid, []))
        n_kept = len(records) - 1
        n_missed = len(missed)
        print(f"  {pid}: {n_kept}/{n_hits} hits resolved (missed={n_missed})  ->  {out_path.name}")
        if missed:
            for m in missed[:5]:
                print(f"      [no FASTA] {m}", file=sys.stderr)
            if len(missed) > 5:
                print(f"      ... and {len(missed)-5} more", file=sys.stderr)
        n_total_written += 1
    print(f"  wrote {n_total_written} {type_suffix} FASTA(s) to {out_dir}/")


def load_query_fastas(pipeline_dir: Path, group_cfg: dict) -> tuple[list[Path], list[Path]]:
    """Resolve [ortholog_blast].query_fastas and query_protein_fastas as Paths."""
    ob = group_cfg.get("ortholog_blast", {})
    nt = [pipeline_dir / p for p in ob.get("query_fastas", [])]
    aa = [pipeline_dir / p for p in ob.get("query_protein_fastas", [])]
    return nt, aa


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pipeline-dir", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--gene-group", default="DMP")
    parser.add_argument("--genome", default="GPE001970_SMEL5", help="Genome subdirectory under <gene-group>/02_BLAST_Alignment/")
    parser.add_argument("--bitscore", type=float, default=None, help="Minimum Bit Score (default: read from [blast_groups_msa].bitscore_threshold or 200)")
    parser.add_argument("--paralogs", default=None, help="Comma-separated paralog IDs (default: read from config or hardcoded 7 SmelDMPs)")
    parser.add_argument("--blast-csv-nt", type=Path, default=None, help="Override blastn CSV path")
    parser.add_argument("--blast-csv-aa", type=Path, default=None, help="Override blastp CSV path")
    parser.add_argument("--anchor-nt", type=Path, default=None, help="Override anchor NT FASTA path")
    parser.add_argument("--anchor-aa", type=Path, default=None, help="Override anchor AA FASTA path")
    parser.add_argument("--out-dir", type=Path, default=None, help="Output directory (default: read from config)")
    parser.add_argument("--types", default="nucleotide,amino_acid", help="Comma-separated subset of {nucleotide, amino_acid} to build")
    parser.add_argument("--overwrite", action="store_true", default=True)
    args = parser.parse_args()

    pipeline_dir: Path = args.pipeline_dir.resolve()
    group = args.gene_group
    group_cfg_dir = pipeline_dir / "config" / group
    if not group_cfg_dir.is_dir():
        print(f"Error: config dir not found: {group_cfg_dir}", file=sys.stderr)
        return 1

    # Merge group TOMLs (read query_fastas from blast config; bitscore from msa config)
    blast_toml = group_cfg_dir / "02_blast_ortholog_alignment.toml"
    msa_toml = group_cfg_dir / "04_multiple_sequence_alignment.toml"
    blast_cfg = tomllib.loads(blast_toml.read_text(encoding="utf-8")) if blast_toml.exists() else {}
    msa_cfg = tomllib.loads(msa_toml.read_text(encoding="utf-8")) if msa_toml.exists() else {}

    bg_cfg = msa_cfg.get("blast_groups_msa", {})
    bitscore = float(args.bitscore) if args.bitscore is not None else float(bg_cfg.get("bitscore_threshold", 200))
    paralog_ids = (
        [p.strip() for p in args.paralogs.split(",") if p.strip()]
        if args.paralogs is not None
        else list(bg_cfg.get("paralog_ids", _DEFAULT_PARALOGS))
    )

    base_dir = pipeline_dir / f"III_RESULT/{group}"
    blast_dir = base_dir / "02_BLAST_Alignment" / args.genome

    nt_csv = args.blast_csv_nt or autodiscover_csv(blast_dir, "n")
    aa_csv = args.blast_csv_aa or autodiscover_csv(blast_dir, "p")

    anchor_nt = args.anchor_nt or (base_dir / "02_BLAST_Alignment" / "hmmer_identified_transcripts.fa")
    anchor_aa = args.anchor_aa or (base_dir / "02_BLAST_Alignment" / "hmmer_identified_proteins.fa")

    out_rel = bg_cfg.get("output_dir", f"04_MSA/merged_input/Selected/v4_BLAST_Groups_bitscore{int(bitscore)}")
    out_dir = args.out_dir or (base_dir / out_rel)

    nt_queries, aa_queries = load_query_fastas(pipeline_dir, blast_cfg)

    print("=== build_v4_blast_groups ===")
    print(f"  pipeline_dir : {pipeline_dir}")
    print(f"  gene_group   : {group}")
    print(f"  genome       : {args.genome}")
    print(f"  bitscore_min : {bitscore}")
    print(f"  paralogs     : {paralog_ids}")
    print(f"  out_dir      : {out_dir}")

    types = {t.strip() for t in args.types.split(",") if t.strip()}
    if "nucleotide" in types:
        if nt_csv is None:
            print("  [skip] no blastn CSV located; pass --blast-csv-nt to override", file=sys.stderr)
        else:
            build_for_type(nt_csv, anchor_nt, nt_queries, paralog_ids, bitscore, out_dir, "NUCLEOTIDE")
    if "amino_acid" in types:
        if aa_csv is None:
            print("  [skip] no blastp CSV located; pass --blast-csv-aa to override", file=sys.stderr)
        else:
            build_for_type(aa_csv, anchor_aa, aa_queries, paralog_ids, bitscore, out_dir, "AMINO_ACID")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
