#!/usr/bin/env python3
"""remap_hap2_coords.py -- remap AtHAP2 (or any reference) deletion coordinates
onto an orthologous HAP2 sequence via a pairwise MAFFT L-INS-i alignment.

Used by Stage 14 (Interaction Domain Mapping) to derive [hap2_variants.coords.sm]
(SmelHAP2) from [hap2_variants.coords.at] (AtHAP2 / Wang 2022 Fig. 3A numbering).
Endpoint policy is conservative: target-side deletion endpoints are pulled
inward to the nearest aligned target residue, so target-specific insertions
relative to the reference are NEVER folded into a deletion implicitly.

Usage (from repo root, inside the `egg` conda env):

    conda run -n egg python3 modules/14_special_pipeline/remap_hap2_coords.py \\
        --ref-fasta    II_INPUTS/HAP2/AtHAP2_At4g11720.fasta \\
        --target-fasta II_INPUTS/HAP2/SmelHAP2.fasta \\
        --ref-coords   14_Interaction_Domain_MappingCONFIG.toml \\
        --ref-section  hap2_variants.coords.at \\
        --target-section hap2_variants.coords.sm \\
        --workdir      /tmp/hap2_remap

The script writes <workdir>/pair.fasta and <workdir>/aligned.fasta and prints
the remapped TOML block to stdout. The block is intended to be pasted back
into the per-stage config (it does not edit the TOML in place).
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


def read_fasta(path: Path) -> tuple[str, str]:
    name, parts = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
        elif name is not None:
            parts.append(line.strip())
    if name is None:
        raise ValueError(f"no records in {path}")
    return name, "".join(parts)


def parse_toml_section(toml_path: Path, section: str) -> dict[str, str]:
    """Read [section] as name -> coord-string. Tolerates inline comments."""
    body, in_section = [], False
    header = f"[{section}]"
    for line in toml_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_section = (stripped == header)
            continue
        if in_section and "=" in stripped and not stripped.startswith("#"):
            body.append(stripped)
    coords: dict[str, str] = {}
    for line in body:
        key, _, rhs = line.partition("=")
        rhs = rhs.split("#", 1)[0].strip()
        m = re.match(r'^"([^"]*)"', rhs)
        if not m:
            continue
        val = m.group(1).strip()
        if val and val.upper() != "TBD":
            coords[key.strip()] = val
    return coords


def parse_ranges(spec: str) -> list[tuple[int, int]]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        a, _, b = part.partition("-")
        out.append((int(a), int(b)))
    return out


def run_mafft(pair_fa: Path, aligned_fa: Path, threads: int = 4) -> None:
    mafft = shutil.which("mafft")
    if mafft is None:
        raise RuntimeError("mafft not on PATH (activate the `egg` conda env)")
    cmd = [mafft, "--localpair", "--maxiterate", "1000", "--thread", str(threads), str(pair_fa)]
    with open(aligned_fa, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)


def build_col_maps(aln_ref: str, aln_tgt: str) -> tuple[list[int], list[int], list[int | None], list[int | None]]:
    assert len(aln_ref) == len(aln_tgt), "alignment length mismatch"
    ncols = len(aln_ref)
    col2ref = [0] * ncols
    col2tgt = [0] * ncols
    ri = ti = 0
    for c in range(ncols):
        if aln_ref[c] != "-":
            ri += 1
            col2ref[c] = ri
        if aln_tgt[c] != "-":
            ti += 1
            col2tgt[c] = ti
    res2col_ref: list[int | None] = [None] * (ri + 2)
    res2col_tgt: list[int | None] = [None] * (ti + 2)
    for c in range(ncols):
        if col2ref[c]:
            res2col_ref[col2ref[c]] = c
        if col2tgt[c]:
            res2col_tgt[col2tgt[c]] = c
    return col2ref, col2tgt, res2col_ref, res2col_tgt


def map_range(a: int, b: int, res2col_ref, col2tgt) -> tuple[int, int] | None:
    """Conservative: pull both endpoints inward to the nearest aligned target residue."""
    col_a, col_b = res2col_ref[a], res2col_ref[b]
    tgt_start = next((col2tgt[c] for c in range(col_a, col_b + 1) if col2tgt[c]), None)
    if tgt_start is None:
        return None
    tgt_end = next((col2tgt[c] for c in range(col_b, col_a - 1, -1) if col2tgt[c]), None)
    return tgt_start, tgt_end


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref-fasta", type=Path, required=True)
    p.add_argument("--target-fasta", type=Path, required=True)
    p.add_argument("--ref-coords", type=Path, required=True, help="TOML config with the reference coord block.")
    p.add_argument("--ref-section", default="hap2_variants.coords.at")
    p.add_argument("--target-section", default="hap2_variants.coords.sm")
    p.add_argument("--workdir", type=Path, default=Path("/tmp/hap2_remap"))
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args(argv)

    args.workdir.mkdir(parents=True, exist_ok=True)
    pair_fa = args.workdir / "pair.fasta"
    aligned_fa = args.workdir / "aligned.fasta"
    pair_fa.write_text(args.ref_fasta.read_text().rstrip() + "\n" + args.target_fasta.read_text().rstrip() + "\n")

    print(f"[remap] aligning {args.ref_fasta.name} + {args.target_fasta.name} with MAFFT L-INS-i ...", file=sys.stderr)
    run_mafft(pair_fa, aligned_fa, threads=args.threads)

    records: list[tuple[str, str]] = []
    name, parts = None, []
    for line in aligned_fa.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].split()[0], []
        elif name is not None:
            parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    if len(records) != 2:
        print(f"[remap] expected 2 aligned records, got {len(records)}", file=sys.stderr)
        return 1
    (ref_name, ref_aln), (tgt_name, tgt_aln) = records
    col2ref, col2tgt, res2col_ref, res2col_tgt = build_col_maps(ref_aln, tgt_aln)
    ref_len = sum(1 for v in col2ref if v)
    tgt_len = sum(1 for v in col2tgt if v)
    print(f"[remap] ungapped: ref={ref_len} aa, target={tgt_len} aa", file=sys.stderr)

    ref_coords = parse_toml_section(args.ref_coords, args.ref_section)
    if not ref_coords:
        print(f"[remap] no usable rows in [{args.ref_section}] of {args.ref_coords}", file=sys.stderr)
        return 1

    print(f"# Remapped from [{args.ref_section}] via MAFFT L-INS-i pairwise alignment.")
    print(f"# Source FASTAs: {args.ref_fasta.name} ({ref_len} aa) + {args.target_fasta.name} ({tgt_len} aa).")
    print(f"[{args.target_section}]")
    for name, spec in ref_coords.items():
        ranges = parse_ranges(spec)
        out_parts, lost = [], False
        for (a, b) in ranges:
            mapped = map_range(a, b, res2col_ref, col2tgt)
            if mapped is None:
                lost = True
                break
            out_parts.append(f"{mapped[0]}-{mapped[1]}")
        if lost:
            print(f'{name:<11} = "TBD"  # <- {spec} (no target coverage)')
        else:
            print(f'{name:<11} = "{",".join(out_parts)}"  # <- {spec}')
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
