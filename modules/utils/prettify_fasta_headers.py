#!/usr/bin/env python3
"""
Prettify FASTA headers using the DMP BLAST short-label rules.

Rewrites every ``>...`` header so the first whitespace token (the sequence ID)
is replaced by the publication-friendly label produced by
``dmp_query_labels.short_label()``. Sequences are left untouched.

Header forms handled (see ``dmp_query_labels`` for the full list):
  ``>lcl|XM_003621193.2_cds_XP_003621241.1_1 [gene=...]``  ->  ``>MtDMP9``
  ``>CDX74441.``                                            ->  ``>CDX74441``
  ``>mRNA.BjuA04g10430S.``                                  ->  ``>BjuA04g10430S``
  ``>SMEL5_10g017610.1``                                    ->  ``>SmelDMP10.610``

Usage:
  python3 prettify_fasta_headers.py <fasta>...
      Rewrite each FASTA in place (after creating a .bak backup).

  python3 prettify_fasta_headers.py --inplace=false --out OUT.fa IN.fa
      Write a single prettified FASTA to a new path.

  python3 prettify_fasta_headers.py --report <fasta>
      Print the (old, new) header pairs without modifying the file.

  python3 prettify_fasta_headers.py --keep-original-as-comment <fasta>
      Append the original ID after the new label as a free-form comment:
      ``>MtDMP9  lcl|XM_003621193.2_cds_XP_003621241.1_1``
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Make the shared label module importable regardless of how we're run.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from dmp_query_labels import short_label  # noqa: E402


def parse_fasta(path: Path):
    """Yield (header_no_gt, sequence_lines_list) tuples."""
    header, chunks = None, []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, chunks
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line)
        if header is not None:
            yield header, chunks


def rewrite_header(header: str, *, keep_original: bool = False) -> tuple[str, bool]:
    """Return (new_header, changed) for the given header (without the leading '>').

    The first whitespace token is mapped through ``short_label()``. Annotation
    suffixes (like ``[gene=...]``) are dropped, since they bloat downstream
    tree-tip labels. If ``keep_original`` is True, the original first-token is
    appended as ``"<new>  <orig_token>"`` so the audit trail is preserved.
    """
    if not header:
        return header, False
    parts = header.split(None, 1)
    orig_id = parts[0]
    new_id = short_label(orig_id)
    if new_id == orig_id and len(parts) == 1:
        return header, False
    if keep_original and new_id != orig_id:
        new_header = f"{new_id}  {orig_id}"
    else:
        new_header = new_id
    return new_header, new_header != header


def write_fasta(records, path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for i, (header, chunks) in enumerate(records):
            if i > 0:
                fh.write("\n")
            fh.write(f">{header}\n")
            seq = "".join(chunks)
            if seq:
                for j in range(0, len(seq), width):
                    fh.write(seq[j : j + width] + "\n")


def process_file(
    path: Path,
    *,
    inplace: bool = True,
    backup_suffix: str = ".bak_pretty",
    keep_original: bool = False,
    report_only: bool = False,
    out_path: Path | None = None,
) -> tuple[int, int]:
    """Returns (n_records, n_renamed). When two records map to the same short
    label, append _2, _3, ... to preserve uniqueness for downstream MSA/tree
    tools that require distinct sequence names."""
    records_in = list(parse_fasta(path))
    records_out = []
    n_renamed = 0
    label_counts: dict[str, int] = {}
    for header, chunks in records_in:
        new_header, changed = rewrite_header(header, keep_original=keep_original)
        # Disambiguate collisions on the new ID (first whitespace token).
        new_id = new_header.split()[0] if new_header else ""
        label_counts[new_id] = label_counts.get(new_id, 0) + 1
        if label_counts[new_id] > 1:
            suffix = f"_{label_counts[new_id]}"
            if " " in new_header:
                new_header = new_header.replace(new_id, new_id + suffix, 1)
            else:
                new_header = new_id + suffix
            changed = True
        if changed:
            n_renamed += 1
        if report_only:
            old_id = header.split()[0] if header else ""
            shown_id = new_header.split()[0] if new_header else ""
            marker = "*" if changed else " "
            print(f"  {marker} {old_id}  ->  {shown_id}")
        records_out.append((new_header, chunks))

    if report_only:
        return len(records_in), n_renamed

    if out_path is None:
        out_path = path

    if inplace and n_renamed > 0:
        backup = path.with_suffix(path.suffix + backup_suffix)
        if not backup.exists():
            backup.write_bytes(path.read_bytes())

    write_fasta(records_out, out_path)
    return len(records_in), n_renamed


def is_clustal(path: Path) -> bool:
    """Sniff the first non-empty line for the CLUSTAL signature."""
    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if line:
                    return line.startswith("CLUSTAL")
    except OSError:
        pass
    return False


def process_clustal(
    path: Path,
    *,
    inplace: bool = True,
    backup_suffix: str = ".bak_pretty",
    report_only: bool = False,
    out_path: Path | None = None,
) -> tuple[int, int]:
    """Rewrite the first whitespace token of every taxon-sequence row in a
    CLUSTAL .aln file using short_label(), preserving the column at which the
    sequence starts. Same-name collisions get _2, _3 disambiguators applied
    consistently across all blocks (so every block keeps the same row order)."""
    raw = path.read_text(encoding="utf-8").splitlines(keepends=False)
    # Pass 1: build the rename map by scanning all unique first-tokens that
    # appear in taxon rows (lines starting with a non-whitespace character that
    # contain at least one whitespace gap to a sequence column).
    rename: dict[str, str] = {}
    label_counts: dict[str, int] = {}
    for line in raw:
        if not line or line.startswith("CLUSTAL") or line[0].isspace():
            continue
        m = re.match(r"^(\S+)(\s+)(\S.*)$", line)
        if not m:
            continue
        old_id = m.group(1)
        if old_id in rename:
            continue
        new_id = short_label(old_id)
        # Disambiguate collisions across distinct old_ids that map to the same
        # short label.
        label_counts[new_id] = label_counts.get(new_id, 0) + 1
        if label_counts[new_id] > 1:
            new_id = f"{new_id}_{label_counts[new_id]}"
        rename[old_id] = new_id

    n_renamed = sum(1 for k, v in rename.items() if k != v)

    # Pass 2: rewrite each row, preserving the sequence column.
    out_lines = []
    for line in raw:
        if not line or line.startswith("CLUSTAL") or line[0].isspace():
            out_lines.append(line)
            continue
        m = re.match(r"^(\S+)(\s+)(\S.*)$", line)
        if not m:
            out_lines.append(line)
            continue
        old_id, gap, rest = m.groups()
        new_id = rename.get(old_id, old_id)
        target_col = len(old_id) + len(gap)
        new_gap_len = max(1, target_col - len(new_id))
        out_lines.append(f"{new_id}{' ' * new_gap_len}{rest}")

    if report_only:
        for k, v in rename.items():
            marker = "*" if k != v else " "
            print(f"  {marker} {k}  ->  {v}")
        return len(rename), n_renamed

    if out_path is None:
        out_path = path
    if inplace and n_renamed > 0:
        backup = path.with_suffix(path.suffix + backup_suffix)
        if not backup.exists():
            backup.write_bytes(path.read_bytes())
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8", newline="\n")
    return len(rename), n_renamed


def main() -> int:
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    p.add_argument("fastas", nargs="+", type=Path)
    p.add_argument("--out", type=Path, default=None, help="Write to this path instead of in-place. Only valid with one input.")
    p.add_argument("--report", action="store_true", help="Print rename pairs without modifying files")
    p.add_argument("--keep-original-as-comment", action="store_true", help="Append original ID after the new label, separated by two spaces")
    p.add_argument("--no-inplace", dest="inplace", action="store_false", default=True, help="Skip the .bak_pretty backup and don't modify the source")
    args = p.parse_args()

    if args.out is not None and len(args.fastas) != 1:
        print("Error: --out requires exactly one input FASTA", file=sys.stderr)
        return 1

    total_renamed = 0
    for fasta in args.fastas:
        if not fasta.exists():
            print(f"  [skip missing] {fasta}", file=sys.stderr)
            continue
        if is_clustal(fasta):
            n_in, n_ren = process_clustal(
                fasta,
                inplace=args.inplace,
                report_only=args.report,
                out_path=args.out,
            )
            fmt = "clustal"
        else:
            n_in, n_ren = process_file(
                fasta,
                inplace=args.inplace,
                keep_original=args.keep_original_as_comment,
                report_only=args.report,
                out_path=args.out,
            )
            fmt = "fasta"
        total_renamed += n_ren
        action = "would-rename" if args.report else "renamed"
        print(f"  [{fmt}] {fasta}: {n_ren}/{n_in} {action}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
