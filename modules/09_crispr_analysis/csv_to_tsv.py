#!/usr/bin/env python3
"""Convert every *.csv under <dir> to a sibling *.tsv (idempotent).

CRISPR-P v2.0 exports comma-delimited .csv. The downstream pipeline globs
for both .tsv and .csv, but tab is the canonical pipeline delimiter, so
materializing TSV siblings up-front keeps every stage on a single format
and removes the comma-vs-tab branch from hot paths.

A .tsv is regenerated only when missing or older than its .csv source,
so reruns are cheap. Originals are preserved.
"""
import csv
import sys
from pathlib import Path


def convert(directory: Path) -> tuple[int, int]:
    converted = skipped = 0
    for csv_path in sorted(directory.glob("*.csv")):
        tsv_path = csv_path.with_suffix(".tsv")
        if tsv_path.exists() and tsv_path.stat().st_mtime >= csv_path.stat().st_mtime:
            skipped += 1
            continue
        with open(csv_path, newline="") as fin, open(tsv_path, "w", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
            writer.writerows(csv.reader(fin))
        converted += 1
    return converted, skipped


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: csv_to_tsv.py <directory>", file=sys.stderr)
        return 2
    target = Path(sys.argv[1])
    if not target.is_dir():
        print(f"[csv_to_tsv] not a directory: {target}", file=sys.stderr)
        return 1
    converted, skipped = convert(target)
    print(f"[csv_to_tsv] {target}: converted={converted}, up-to-date={skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
