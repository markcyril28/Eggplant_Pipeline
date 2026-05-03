#!/usr/bin/env python3
"""Merge MFEprimer + isPcr per-genome amplicon TSVs into one summary table.

Called by 12_In_Silico_PCR.sh after both engines complete. Combines:
  PCR_DIR/{genome}/02_MFEprimer/{set}__{genome}.tsv
  PCR_DIR/{genome}/03_isPcr/{set}__{genome}.tsv
into a long-format TSV with one row per (genome, engine, primer_id, hit).
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pcr-dir", required=True, type=Path)
    ap.add_argument("--set-name", required=True)
    ap.add_argument("--genomes", required=True, nargs="+")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    fields = [
        "engine", "genome", "primer_id", "chrom", "start", "end",
        "strand", "size", "tm_f", "tm_r", "dg_f", "dg_r", "product_gc",
        "forward", "reverse",
    ]

    rows: list[dict] = []
    for genome in args.genomes:
        mfe = read_tsv(args.pcr_dir / genome / "02_MFEprimer" /
                       f"{args.set_name}__{genome}.tsv")
        for r in mfe:
            r["engine"] = "mfeprimer"
            rows.append(r)
        isp = read_tsv(args.pcr_dir / genome / "03_isPcr" /
                       f"{args.set_name}__{genome}.tsv")
        for r in isp:
            r["engine"] = "ispcr"
            rows.append(r)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    n_mfe = sum(1 for r in rows if r["engine"] == "mfeprimer")
    n_isp = sum(1 for r in rows if r["engine"] == "ispcr")
    print(f"[summarize] {args.set_name}: MFEprimer={n_mfe}, isPcr={n_isp} hits "
          f"-> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
