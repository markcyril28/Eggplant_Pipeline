#!/usr/bin/env python3
"""Extract AF3 backup zips sitting at _AF3_Backup/ (or _AF3_Backup/<stoich>/)
into 02_Complexes/<stoich>/<pair>/ slots.

Complements the orchestrator's inline SCATTER pass, which handles zips
dropped at $EXP_OUT_DIR/, $EXP_OUT_DIR/02_Complexes/, and
$EXP_OUT_DIR/02_Complexes/<stoich>/. SCATTER deletes the source zip on
success; this script does NOT (the backup mirror is immutable).

Filename grammar (lowercase tokens, separated by underscores):
    fold_<stoich>_<variant>_wt_eggplant[_<tag>].zip   -> <variant>__WT pair
    fold_<stoich>_wt_<variant>_eggplant[_<tag>].zip   -> WT__<variant> pair
    fold_<stoich>_hap2_and_<gene>_<tag>.zip           -> AF3 default; matched
                                                         against --pair-labels
                                                         (if provided) when
                                                         exactly one empty slot
                                                         is plausible
<stoich> token mapping is configurable via --stoich-tokens (default:
"monomeric:monomeric,trimeric:postfusion_like,dimeric:dimeric") so the
filename "trimeric" routes to the 02_Complexes/postfusion_like/ subfolder
matching the project's stoichiometry labels.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


# Lowercase tokens used in [hap2_variants].rows / [dmp_variants].rows ->
# canonical camelCase variant name. Mirrors the rows defined in
# 14_Interaction_Domain_MappingCONFIG.toml. Keep in sync when new rows are added.
_CAMEL = {
    "wt": "WT",
    "delc": "delC", "deln": "delN", "delecto": "delEcto",
    "delectod2": "delEctoD2", "delfl": "delFL",
    "delpretmd": "delPreTMD", "deltmd": "delTMD", "deljuxtamem": "delJuxtaMem",
    "delpretmdandtmd": "delPreTMDAndTMD",
    "delpretmdandtmdandjuxtamem": "delPreTMDAndTMDAndJuxtaMem",
    "delectoandc": "delEctoAndC", "deltmdcore": "delTMDcore",
}
_GUIDE_RE = re.compile(r"^(del|fs)guide(\d+)$")


def camel_variant(token: str) -> str | None:
    """delguide4 -> delGuide4 ; fsguide17 -> fsGuide17 ; fall back to _CAMEL map."""
    m = _GUIDE_RE.match(token)
    if m:
        kind = "delGuide" if m.group(1) == "del" else "fsGuide"
        return f"{kind}{m.group(2)}"
    return _CAMEL.get(token)


def parse_stoich_tokens(spec: str) -> dict[str, str]:
    """\"monomeric:monomeric,trimeric:postfusion_like\" -> {\"monomeric\":\"monomeric\", ...}."""
    out: dict[str, str] = {}
    for pair in spec.split(","):
        pair = pair.strip()
        if not pair:
            continue
        k, _, v = pair.partition(":")
        out[k.strip().lower()] = (v.strip() or k.strip()).lower()
    return out


def parse_zip(name: str, stoich_map: dict[str, str]) -> tuple[str, str] | None:
    """Return (target_stoich, pair_label) or None if filename does not match.

    Recognises the two domain-mapping naming patterns:
        fold_<stoich>_<variant>_wt_eggplant[_<tag>].zip
        fold_<stoich>_wt_<variant>_eggplant[_<tag>].zip
    where <variant> may itself contain underscores (e.g., "delpretmdandtmd")
    and <tag> is an optional trailing suffix ("redo", "swiss_template", ...).
    """
    if not name.endswith(".zip") or not name.startswith("fold_"):
        return None
    stem = name[len("fold_"):-len(".zip")]
    parts = stem.split("_")
    try:
        egg_idx = parts.index("eggplant")
    except ValueError:
        return None
    if not parts:
        return None
    stoich_tok = parts[0]
    if stoich_tok not in stoich_map:
        return None
    middle = parts[1:egg_idx]
    if len(middle) < 2 or "wt" not in middle:
        return None
    if middle[-1] == "wt":            # <hap2_variant>_wt
        hap2 = camel_variant("".join(middle[:-1]))
        dmp = "WT"
    elif middle[0] == "wt":           # wt_<dmp_variant>
        hap2 = "WT"
        dmp = camel_variant("".join(middle[1:]))
    else:
        return None
    if hap2 is None or dmp is None:
        return None
    return stoich_map[stoich_tok], f"{hap2}__{dmp}"


def discover_zips(backup_dir: Path) -> list[Path]:
    """All zips anywhere under backup_dir (recursive). Supports flexible
    organisational layouts: parent-level, per-stoich subdirs, batch_NN
    subdirs (one per submission session), or arbitrarily nested combos.
    The zip's filename is the only thing that controls the target slot, so
    physical location is purely organisational."""
    if not backup_dir.exists():
        return []
    return sorted(backup_dir.rglob("*.zip"))


def extract_zip(zip_path: Path, target: Path, clean: bool, dry_run: bool) -> tuple[bool, int]:
    """Extract `zip_path` flat into `target`, optionally cleaning target first.
    Returns (ok, n_files_after)."""
    if dry_run:
        return True, -1
    if clean and target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)
    rc = subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", str(target)],
                        check=False).returncode
    if rc != 0:
        return False, 0
    n = sum(1 for _ in target.rglob("*") if _.is_file())
    return True, n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backup-dir", required=True,
                    help="path to _AF3_Backup root (parent-level zips + per-stoich subfolders are both swept)")
    ap.add_argument("--complex-dir", required=True,
                    help="path to 02_Complexes root (target tree)")
    ap.add_argument("--stoich-tokens",
                    default="monomeric:monomeric,trimeric:postfusion_like,dimeric:dimeric",
                    help="comma-separated FILENAME_TOKEN:OUTPUT_STOICH mappings (default maps trimeric->postfusion_like)")
    ap.add_argument("--no-clean", action="store_true",
                    help="do NOT wipe target slot before extracting (default: clean+extract)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report planned actions without touching disk")
    args = ap.parse_args()

    backup = Path(args.backup_dir)
    complex_root = Path(args.complex_dir)
    if not backup.exists():
        print(f"[distribute] no _AF3_Backup dir at {backup}; nothing to do.")
        return 0

    stoich_map = parse_stoich_tokens(args.stoich_tokens)
    zips = discover_zips(backup)
    if not zips:
        print(f"[distribute] no zips found under {backup}")
        return 0

    print(f"[distribute] scanning {len(zips)} zip(s) under {backup}")
    n_ok = n_skip = n_fail = 0
    for zp in zips:
        parsed = parse_zip(zp.name, stoich_map)
        if parsed is None:
            print(f"  [SKIP] {zp.name}  (unrecognised filename pattern)")
            n_skip += 1
            continue
        stoich, pair = parsed
        target = complex_root / stoich / pair
        rel_zip = zp.relative_to(backup) if zp.is_relative_to(backup) else zp
        if args.dry_run:
            print(f"  [DRY ] {rel_zip}  ->  {stoich}/{pair}/")
            continue
        ok, n_files = extract_zip(zp, target, clean=(not args.no_clean), dry_run=False)
        if ok:
            print(f"  [OK  ] {rel_zip}  ->  {stoich}/{pair}/  ({n_files} files)")
            n_ok += 1
        else:
            print(f"  [FAIL] {rel_zip}  ->  {stoich}/{pair}/  (unzip exited non-zero)")
            n_fail += 1

    print(f"\n[distribute] extracted={n_ok} skipped={n_skip} failed={n_fail}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
