#!/usr/bin/env python3
"""Merge MFEprimer + isPcr per-genome amplicon TSVs into one summary table.

Called by 12_In_Silico_PCR.sh after both engines complete. Reads:
  PCR_DIR/{genome}/02_MFEprimer/{set}__{genome}.tsv
  PCR_DIR/{genome}/02_MFEprimer/{set}__{genome}.bands_per_pair.tsv
  PCR_DIR/{genome}/03_isPcr/{set}__{genome}.tsv

Writes (under PCR_DIR/04_Summary/):
  {set}_amplicons.tsv                    long-form, one row per amplicon hit
  {set}_amplicons_bands_per_pair.tsv     long-form, one row per (engine,
                                          genome, primer_pair) with band
                                          counts; includes 0-band pairs
  {set}_bands_per_pair_wide.tsv          wide pivot, one row per primer pair
                                          with total_bands per genome and
                                          a row total for at-a-glance triage

Off-target conventions:
  - clean_bands    : amplicons where F and R primer IDs share the same pair
  - chimeric_bands : F from one pair paired with R from another (cross-talk).
                     A chimeric amplicon credits BOTH constituent pairs.
  - total_bands    : clean + chimeric per pair per genome
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict
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

    # ── 1. Long-form amplicons table (one row per hit) ────────────────────
    amp_fields = [
        "engine", "genome", "primer_id", "f_primer_id", "r_primer_id",
        "is_chimeric", "chrom", "start", "end", "strand", "size",
        "tm_f", "tm_r", "dg_f", "dg_r", "product_gc", "ppc",
        "forward", "reverse",
    ]
    rows: list[dict] = []
    for genome in args.genomes:
        for r in read_tsv(args.pcr_dir / genome / "02_MFEprimer" /
                          f"{args.set_name}__{genome}.tsv"):
            r["engine"] = "mfeprimer"
            rows.append(r)
        for r in read_tsv(args.pcr_dir / genome / "03_isPcr" /
                          f"{args.set_name}__{genome}.tsv"):
            r["engine"] = "ispcr"
            r.setdefault("is_chimeric", "no")
            rows.append(r)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=amp_fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in amp_fields})

    # ── 2. Long-form bands-per-pair table ─────────────────────────────────
    # Prefer the per-genome MFEprimer .bands_per_pair.tsv (already lists every
    # input pair including 0-band ones, and credits chimeric bands to both
    # constituent pairs). Fall back to recomputing from the amplicon rows
    # when the rollup is absent (e.g., isPcr, which never writes one).
    band_long: list[dict] = []
    pair_universe: set[str] = set()
    per_genome_isp: dict[tuple, dict] = defaultdict(
        lambda: {"clean": 0, "chimeric": 0, "chroms": set()})

    for genome in args.genomes:
        mfe_rollup = (args.pcr_dir / genome / "02_MFEprimer" /
                      f"{args.set_name}__{genome}.bands_per_pair.tsv")
        for r in read_tsv(mfe_rollup):
            pid = (r.get("pair_name") or r.get("primer_id") or "").strip()
            if not pid:
                continue
            pair_universe.add(pid)
            band_long.append({
                "engine": "mfeprimer",
                "genome": r.get("genome", genome),
                "primer_id": pid,
                "clean_bands": int(r.get("clean_bands") or 0),
                "chimeric_bands": int(r.get("chimeric_bands") or 0),
                "total_bands": int(r.get("total_bands") or 0),
                "chrom_count": int(r.get("chrom_count") or 0),
                "chroms": r.get("chroms", ""),
            })

    # isPcr: derive band counts from the merged amplicon rows (no rollup).
    for r in rows:
        if r.get("engine") != "ispcr":
            continue
        pid = (r.get("primer_id") or "").strip()
        if not pid:
            continue
        pair_universe.add(pid)
        key = (r.get("genome", ""), pid)
        is_chim = (r.get("is_chimeric", "no").lower() == "yes")
        per_genome_isp[key]["chimeric" if is_chim else "clean"] += 1
        chrom = (r.get("chrom") or "").strip()
        if chrom:
            per_genome_isp[key]["chroms"].add(chrom)

    for (genome, pid), v in per_genome_isp.items():
        band_long.append({
            "engine": "ispcr", "genome": genome, "primer_id": pid,
            "clean_bands": v["clean"],
            "chimeric_bands": v["chimeric"],
            "total_bands": v["clean"] + v["chimeric"],
            "chrom_count": len(v["chroms"]),
            "chroms": ";".join(sorted(v["chroms"])),
        })

    # Backfill: ensure every (engine, genome, pair) appears even with zeros so
    # downstream consumers don't have to detect missing rows.
    seen = {(r["engine"], r["genome"], r["primer_id"]) for r in band_long}
    for engine in ("mfeprimer", "ispcr"):
        for genome in args.genomes:
            for pid in pair_universe:
                if (engine, genome, pid) not in seen:
                    band_long.append({
                        "engine": engine, "genome": genome, "primer_id": pid,
                        "clean_bands": 0, "chimeric_bands": 0,
                        "total_bands": 0, "chrom_count": 0, "chroms": "",
                    })

    bands_path = args.output.with_name(
        args.output.stem + "_bands_per_pair.tsv")
    band_fields = ["engine", "genome", "primer_id",
                   "clean_bands", "chimeric_bands", "total_bands",
                   "chrom_count", "chroms"]
    band_long.sort(key=lambda r: (r["engine"], r["primer_id"], r["genome"]))
    with bands_path.open("w", newline="") as fh:
        bw = csv.DictWriter(fh, fieldnames=band_fields, delimiter="\t")
        bw.writeheader()
        bw.writerows(band_long)

    # # ── 3. Wide pivot: pair × genome (one row per pair, totals per genome) ─
    # # Combines mfeprimer + ispcr counts so a single column shows the full
    # # off-target burden for that pair on that genome. The "all_genomes_*"
    # # columns are sums across genomes for quick triage.
    # wide_path = args.output.with_name(
    #     args.output.stem + "_bands_per_pair_wide.tsv")
    # sorted_pairs = sorted(pair_universe)
    # sorted_genomes = list(args.genomes)
    # pair_genome_total: dict[tuple, int] = defaultdict(int)
    # pair_genome_clean: dict[tuple, int] = defaultdict(int)
    # pair_genome_chim: dict[tuple, int] = defaultdict(int)
    # for r in band_long:
    #     k = (r["primer_id"], r["genome"])
    #     pair_genome_total[k] += r["total_bands"]
    #     pair_genome_clean[k] += r["clean_bands"]
    #     pair_genome_chim[k] += r["chimeric_bands"]

    # wide_header = ["primer_id"]
    # for g in sorted_genomes:
    #     wide_header += [f"{g}__total", f"{g}__clean", f"{g}__chimeric"]
    # wide_header += ["all_genomes_total", "all_genomes_clean",
    #                 "all_genomes_chimeric"]
    # with wide_path.open("w", newline="") as fh:
    #     ww = csv.writer(fh, delimiter="\t")
    #     ww.writerow(wide_header)
    #     for pid in sorted_pairs:
    #         row = [pid]
    #         tot = clean = chim = 0
    #         for g in sorted_genomes:
    #             t = pair_genome_total[(pid, g)]
    #             c = pair_genome_clean[(pid, g)]
    #             x = pair_genome_chim[(pid, g)]
    #             row += [t, c, x]
    #             tot += t; clean += c; chim += x
    #         row += [tot, clean, chim]
    #         ww.writerow(row)

    n_mfe = sum(1 for r in rows if r["engine"] == "mfeprimer")
    n_isp = sum(1 for r in rows if r["engine"] == "ispcr")
    print(f"[summarize] {args.set_name}: MFEprimer={n_mfe}, isPcr={n_isp} hits "
          f"-> {args.output}")
    print(f"[summarize] {args.set_name}: bands-per-pair (long) -> {bands_path}")
    # print(f"[summarize] {args.set_name}: bands-per-pair (wide) -> {wide_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
