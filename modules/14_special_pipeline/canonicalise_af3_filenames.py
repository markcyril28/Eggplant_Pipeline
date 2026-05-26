#!/usr/bin/env python3
"""
Canonicalise AF3 zip-extracted filenames inside 02_Complexes/{stoich}/{pair}/.

Background:
  `distribute_af3_zips.py` matches a zip's filename and routes it to the
  correct slot folder (e.g. WT__fsGuide4) but does NOT rename the interior
  CIFs / JSONs. AF3 jobs submitted with the default job-name template land
  as `fold_<stoich_token>_<hap2_lower>_<dmp_lower>_eggplant[_<tag>]_<rest>`,
  whereas manually-submitted jobs (delN, delC, delTMDcore) land as
  `fold_<HAP2>__<DMP>_<canonical_stoich>_<rest>`. Downstream code
  (`refresh_all_v2.py`, `iptm_heatmap.py`, etc.) expects the canonical form.

This script walks 02_Complexes/{stoich}/{pair}/ and renames the AF3-raw
files to the canonical convention. It is idempotent:
  * Files already in canonical form are skipped.
  * Files whose AF3 raw form cannot be confidently parsed are reported and
    left alone (no destructive moves).

The slot folder name (e.g. "WT__fsGuide4") is the source of truth for the
target pair label; the parent folder name (e.g. "monomeric",
"postfusion_like") is the canonical stoich. The script does NOT trust the
AF3 job-name tokens beyond using them to locate the rest-of-filename split.

Usage:
  python canonicalise_af3_filenames.py --apply        # rename
  python canonicalise_af3_filenames.py                # dry-run (default)
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

COMPLEX_ROOT = Path(r"c:\PIPELINE\Eggplant_Pipeline"
                    r"\III_RESULT\DMP_x_SmelHAP2\14_Domain_Mapping"
                    r"\deletion_ladder\02_Complexes")

# Markers that separate the AF3 job-name prefix from the file-specific tail.
# All AF3 outputs in a job folder end in one of these:
TAIL_MARKERS = (
    "_model_",
    "_summary_confidences_",
    "_full_data_",
    "_job_request",
    "_confidences_",
    "_data_",
)


def canonical_rest(name: str) -> str | None:
    """Strip the AF3 job-name prefix and return the file-specific tail.

    Example:
      "fold_monomeric_wt_fsguide4_eggplant_redo_model_0.cif"
        -> "model_0.cif"
      "fold_monomeric_wt_fsguide4_eggplant_redo_summary_confidences_2.json"
        -> "summary_confidences_2.json"
    """
    for marker in TAIL_MARKERS:
        pos = name.find(marker)
        if pos > 0:
            return name[pos + 1:]  # +1 to drop the leading underscore
    return None


def is_canonical(name: str, pair_label: str, canonical_stoich: str) -> bool:
    return name.startswith(f"fold_{pair_label}_{canonical_stoich}_") \
        or name == f"fold_{pair_label}_{canonical_stoich}"


def plan_renames(root: Path) -> list[tuple[Path, Path, str]]:
    """Return list of (src, dst, reason). reason = "rename" | "skip-canonical" | "skip-unknown"."""
    plan: list[tuple[Path, Path, str]] = []
    if not root.exists():
        print(f"[err] complex root not found: {root}")
        return plan
    for stoich_dir in sorted(root.iterdir()):
        if not stoich_dir.is_dir():
            continue
        canonical_stoich = stoich_dir.name
        for pair_dir in sorted(stoich_dir.iterdir()):
            if not pair_dir.is_dir() or "__" not in pair_dir.name:
                continue
            pair_label = pair_dir.name
            for f in sorted(pair_dir.iterdir()):
                if not f.is_file():
                    continue
                name = f.name
                if not name.startswith("fold_"):
                    continue
                if is_canonical(name, pair_label, canonical_stoich):
                    continue  # already correct
                tail = canonical_rest(name)
                if tail is None:
                    plan.append((f, f, "skip-unknown"))
                    continue
                new_name = f"fold_{pair_label}_{canonical_stoich}_{tail}"
                if new_name == name:
                    continue
                plan.append((f, f.with_name(new_name), "rename"))
    return plan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="actually rename (default = dry-run)")
    ap.add_argument("--root", default=str(COMPLEX_ROOT),
                    help="02_Complexes root")
    ap.add_argument("--verbose", action="store_true",
                    help="print every rename, not just one sample per slot")
    ap.add_argument("--exclude", action="append", default=[],
                    help="pair-label substring to skip (repeatable). Match is on folder name only.")
    args = ap.parse_args()

    root = Path(args.root)
    print(f"[scan] {root}")
    plan = plan_renames(root)
    if not plan:
        print("[done] nothing to do (all files already canonical or no fold_* files found).")
        return 0

    renames_raw = [p for p in plan if p[2] == "rename"]
    unknown = [p for p in plan if p[2] == "skip-unknown"]

    # Apply --exclude filter on folder name.
    def excluded(src: Path) -> bool:
        return any(ex in src.parent.name for ex in args.exclude)

    renames = [p for p in renames_raw if not excluded(p[0])]
    skipped_excluded = [p for p in renames_raw if excluded(p[0])]

    print(f"\n[plan] {len(renames)} rename(s), "
          f"{len(skipped_excluded)} excluded by filter, "
          f"{len(unknown)} unparseable")
    by_pair: dict[str, list[tuple[Path, Path]]] = {}
    for src, dst, _ in renames:
        key = f"{src.parent.parent.name}/{src.parent.name}"
        by_pair.setdefault(key, []).append((src, dst))
    for key in sorted(by_pair):
        files = by_pair[key]
        print(f"  {key}/  ({len(files)} files)")
        if args.verbose:
            for src, dst in files:
                print(f"    {src.name}")
                print(f"      -> {dst.name}")
        else:
            sample = files[0]
            print(f"    e.g. {sample[0].name}")
            print(f"      -> {sample[1].name}")
    if skipped_excluded:
        skipped_pairs = sorted({p[0].parent.name for p in skipped_excluded})
        print(f"\n[excluded] skipped {len(skipped_excluded)} files in: {', '.join(skipped_pairs)}")

    if unknown:
        print("\n[warn] unparseable files (left alone):")
        for src, _, _ in unknown[:10]:
            print(f"  {src}")
        if len(unknown) > 10:
            print(f"  ... and {len(unknown) - 10} more")

    if not args.apply:
        print("\n[dry-run] re-run with --apply to perform renames.")
        return 0

    n_ok = n_collide = 0
    for src, dst, _ in renames:
        if dst.exists():
            print(f"  [SKIP-COLLIDE] {src.name}  (target {dst.name} already exists)")
            n_collide += 1
            continue
        src.rename(dst)
        n_ok += 1
    print(f"\n[apply] renamed={n_ok} collisions={n_collide}")
    return 0 if n_collide == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
