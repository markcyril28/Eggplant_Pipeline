#!/usr/bin/env python3
"""Detect double-extracted AF3 bundles in 02_Complexes/<stoich>/<pair>/
slots and dedupe the safe ones.

A pair folder normally holds exactly one AF3 bundle: a set of files
sharing a common prefix (e.g. `fold_WT__delN_monomeric_*`). When the
same prediction is unzipped twice under two different job names (or the
orchestrator extracts a zip whose top-level files use the AF3 server's
default lowercase name while a prior camelCase extraction still sits in
place), the folder ends up with two parallel prefixes. This module:

  1. Identifies bundles by their `*_summary_confidences_0.json` filename
     (always present in an AF3 download, unlike `*_model_0.cif` which a
     user might remove).
  2. For each folder with >1 distinct prefix, compares the rank-CIFs the
     bundles share (intersection of available ranks).
  3. All overlapping ranks match (sha1) -> safe dedupe; keep the most-
     complete bundle (largest rank set), delete the others.
  4. Any rank disagrees -> divergent; flag, do NOT touch.

`*_canonical/` and `*_SWISS/` sibling folders (intentional archive splits)
are NOT auto-merged; each is treated as an independent pair slot. Pass
`--ignore-pair` to additionally skip specific pair names.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path


SUMMARY_RE = re.compile(r"^(?P<prefix>.+?)_summary_confidences_0\.json$")
RANK_RE = re.compile(r"_model_(?P<rank>\d+)\.cif$")


def sha1_of(p: Path) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def bundle_files(pair_dir: Path, prefix: str) -> list[Path]:
    """Every file (top level + msas/ + templates/) belonging to one prefix."""
    out: list[Path] = []
    if not pair_dir.is_dir():
        return out
    for f in pair_dir.rglob("*"):
        if f.is_file() and f.name.startswith(prefix + "_"):
            out.append(f)
    return out


def bundle_rank_cifs(pair_dir: Path, prefix: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for p in pair_dir.iterdir():
        if not p.is_file() or not p.name.startswith(prefix + "_"):
            continue
        m = RANK_RE.search(p.name)
        if m:
            out[int(m.group("rank"))] = p
    return out


def scan(complex_root: Path, ignore_pairs: set[str]) -> tuple[list, list]:
    """Return (safe_actions, divergent).

    safe_actions: list of (slot, keep_prefix, drop_prefixes, drop_files, ranks_map)
    divergent:    list of (slot, ref_prefix, ranks_map, [(other_prefix, reason)...])
    """
    safe_actions: list = []
    divergent: list = []
    if not complex_root.is_dir():
        return safe_actions, divergent
    for slab_dir in sorted(complex_root.iterdir()):
        if not slab_dir.is_dir():
            continue
        for pair_dir in sorted(slab_dir.iterdir()):
            if not pair_dir.is_dir():
                continue
            if pair_dir.name in ignore_pairs:
                continue
            prefixes = sorted({
                SUMMARY_RE.match(f.name).group("prefix")
                for f in pair_dir.iterdir()
                if f.is_file() and SUMMARY_RE.match(f.name)
            })
            if len(prefixes) <= 1:
                continue
            ranks = {p: bundle_rank_cifs(pair_dir, p) for p in prefixes}
            ref = max(prefixes, key=lambda p: len(ranks[p]))
            slot = f"{slab_dir.name}/{pair_dir.name}"
            divergent_with: list[tuple[str, str]] = []
            for p in prefixes:
                if p == ref:
                    continue
                shared = set(ranks[ref]) & set(ranks[p])
                if not shared:
                    divergent_with.append((p, "no shared ranks to compare"))
                    continue
                for r in sorted(shared):
                    if sha1_of(ranks[ref][r]) != sha1_of(ranks[p][r]):
                        divergent_with.append((p, f"rank {r} differs"))
                        break
            if divergent_with:
                divergent.append((slot, ref, ranks, divergent_with))
            else:
                drop_prefixes = [p for p in prefixes if p != ref]
                drop_files = [f for p in drop_prefixes for f in bundle_files(pair_dir, p)]
                safe_actions.append((slot, ref, drop_prefixes, drop_files, ranks))
    return safe_actions, divergent


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--complex-dir", required=True,
                    help="02_Complexes root to scan (per-stoich subfolders are walked)")
    ap.add_argument("--apply", action="store_true",
                    help="actually delete redundant files (default: dry-run)")
    ap.add_argument("--ignore-pair", default="",
                    help="comma-separated pair names to skip entirely")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress per-folder output; only print summary line")
    args = ap.parse_args()

    complex_root = Path(args.complex_dir)
    ignore_pairs = {p.strip() for p in args.ignore_pair.split(",") if p.strip()}
    safe_actions, divergent = scan(complex_root, ignore_pairs)

    if divergent and not args.quiet:
        print(f"[dedupe] DIVERGENT bundles (NOT touched; review manually):")
        for slot, ref, ranks, dw in divergent:
            print(f"  {slot}")
            print(f"    ref  : {ref}  (ranks {sorted(ranks[ref])})")
            for p, why in dw:
                print(f"    vs   : {p}  (ranks {sorted(ranks[p])})  -> {why}")

    if not safe_actions:
        msg = f"[dedupe] {len(divergent)} divergent; 0 safe-dedupe candidates"
        print(msg)
        return 0

    if not args.quiet:
        print(f"[dedupe] {len(safe_actions)} safe-dedupe candidate(s):")
        for slot, keep, drop, drop_files, ranks in safe_actions:
            print(f"  {slot}")
            print(f"    KEEP : {keep}  (ranks {sorted(ranks[keep])})")
            for p in drop:
                cnt = sum(1 for f in drop_files if f.name.startswith(p + "_"))
                print(f"    DROP : {p}  (ranks {sorted(ranks[p])}, {cnt} files)")

    if not args.apply:
        n_files = sum(len(a[3]) for a in safe_actions)
        print(f"[dedupe] DRY-RUN: would delete {n_files} files across "
              f"{len(safe_actions)} folder(s). Re-run with --apply to execute.")
        return 0

    deleted = 0
    for slot, keep, drop, drop_files, ranks in safe_actions:
        for fp in drop_files:
            try:
                fp.unlink()
                deleted += 1
            except OSError as exc:
                print(f"  [WARN] failed to delete {fp}: {exc}", file=sys.stderr)
    print(f"[dedupe] APPLIED: deleted {deleted} file(s) across "
          f"{len(safe_actions)} folder(s); {len(divergent)} divergent left as-is.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
