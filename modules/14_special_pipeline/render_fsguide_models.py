#!/usr/bin/env python3
"""
Render the missing fsGuide model PSE + structure PNGs into
SmelDMP_variants/frameshift_model_{monomeric,postfusion}/.

For each fsGuide variant this writes (per stoichiometry):
  hap2_dmp_<v>_model_<stoich>.pse                 (PyMOL session)
  hap2_dmp_<v>_model_<stoich>_structure.png       (rendered image)
  hap2_dmp_<v>_model_<stoich>_structure_cropped.png

Uses the WT Figure-0 view (fig0_monomeric.json / fig0_trimer.json) as the
starting orientation; the user is expected to re-open each .pse in PyMOL,
reorient as needed, save, then re-run render_from_pse.py to refresh the PNGs.

Inputs (canonical naming after canonicalise_af3_filenames.py):
  02_Complexes/monomeric/WT__<v>/fold_WT__<v>_monomeric_model_0.cif
  02_Complexes/postfusion_like/WT__<v>/fold_WT__<v>_postfusion_like_model_0.cif

Idempotent: skips a target whose .pse already exists (so user re-orientations
are not clobbered) UNLESS --force is supplied.
"""
from pathlib import Path
import argparse
import subprocess
import sys

PYMOL_ENV = r"C:\ProgramData\pymol"
PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
DMP_DIR = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/07_Summary/deletion_variants/SmelDMP_variants"
FS_MOD_MO = DMP_DIR / "frameshift_model_monomeric"
FS_MOD_PF = DMP_DIR / "frameshift_model_postfusion"
COMPLEXES = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/02_Complexes"
HAP2_SCRIPT = PROJECT / "modules/14_special_pipeline/render_hap2_domain_map.py"
CROP = PROJECT / "modules/14_special_pipeline/crop_to_content.py"
VIEW_TMP = PROJECT / "modules/14_special_pipeline/.views"
MONO_VIEW = VIEW_TMP / "fig0_monomeric.json"
TRI_VIEW = VIEW_TMP / "fig0_trimer.json"

# (variant, last_kept_aa, tail_len) - mirrors refresh_smeldmp_panels.py
DMP_FRAMESHIFT = [
    ("fsGuide4",  7,   3),
    ("fsGuide17", 118, 23),
    ("fsGuide16", 117, 24),
    ("fsGuide20", 128, 13),
    ("fsGuide37", 215, 25),
    ("fsGuide46", 27,  15),
    ("fsGuide50", 80,  3),
]


def run_pymol(*args: str) -> int:
    cmd = ["conda", "run", "-p", PYMOL_ENV, "pymol", "-cq", *args]
    r = subprocess.run(cmd, cwd=PROJECT, capture_output=True, text=True)
    for line in (r.stdout + r.stderr).splitlines():
        if any(tag in line for tag in ("[OK]", "[INFO]", "[WARN]", "[ERROR]")):
            print(" ", line)
    return r.returncode


def run_py(*args: str) -> int:
    r = subprocess.run([sys.executable, *args], cwd=PROJECT,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
    for line in (r.stdout + r.stderr).splitlines():
        if any(tag in line for tag in ("[OK]", "[WARN]", "[ERROR]", "cropped")):
            print(" ", line)
    return r.returncode


def cif_path(v: str, stoich: str) -> Path:
    return COMPLEXES / stoich / f"WT__{v}" / f"fold_WT__{v}_{stoich}_model_0.cif"


def render_one(v: str, stoich: str, out_dir: Path, view_json: Path,
               dmp_chains: str, hap2_chains: str, last_kept: int,
               tail_len: int, force: bool) -> None:
    # File naming convention (matches the delN / delC / delTMDcore deletion
    # variants in model_{monomeric,postfusion}/):
    #   <base>.pse              <- PyMOL session
    #   <base>_structure.png    <- rendered image (renderer auto-appends _structure)
    #   <base>_structure_cropped.png
    # `base` carries NO _structure suffix; render_hap2_domain_map.py adds it.
    suffix = "monomeric" if stoich == "monomeric" else "postfusion"
    base = out_dir / f"hap2_dmp_{v}_model_{suffix}"
    pse = base.with_suffix(".pse")
    out_arg = base.with_suffix(".png")  # what we pass to --out
    structure = base.with_name(base.name + "_structure.png")
    cropped = base.with_name(base.name + "_structure_cropped.png")

    cif = cif_path(v, stoich)
    if not cif.exists():
        print(f"[SKIP] {v} {stoich}: missing CIF {cif.name}")
        return

    if pse.exists() and not force:
        print(f"[SKIP] {v} {stoich}: .pse already exists (pass --force to redo)")
        if not structure.exists():
            # Render PNG from existing pse so the folder is not partial.
            run_pymol(str(PROJECT / "modules/14_special_pipeline/render_from_pse.py"),
                      "--", "--pse", str(pse), "--out", str(structure))
        if structure.exists() and not cropped.exists():
            run_py(str(CROP), "--image", str(structure),
                   "--out", str(cropped), "--pad", "12")
        return

    # Fresh render from the CIF using the WT Figure-0 view as the starting
    # orientation. render_hap2_domain_map.py saves the .pse next to the
    # output PNG so the user can reopen, reorient, and re-render later.
    resi_lost = f"{last_kept + 1}-222"
    args = [
        str(HAP2_SCRIPT), "--",
        "--cif", str(cif),
        "--out", str(out_arg),
        "--hap2-chains", hap2_chains,
        "--dmp-chains", dmp_chains,
        "--view-json", str(view_json),
        "--dmp-deleted-residues", resi_lost,
    ]
    if tail_len > 0:
        args.extend(["--dmp-novel-tail", str(tail_len)])
    print(f"[RUN]  {v} {stoich}: rendering {structure.name}")
    rc = run_pymol(*args)
    if rc != 0:
        print(f"[FAIL] {v} {stoich}: PyMOL exited {rc}")
        return
    if structure.exists():
        run_py(str(CROP), "--image", str(structure),
               "--out", str(cropped), "--pad", "12")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing .pse / .png (default: skip if .pse exists)")
    ap.add_argument("--only", default="", help="comma-separated subset of variants (e.g. fsGuide4,fsGuide17)")
    args = ap.parse_args()

    FS_MOD_MO.mkdir(parents=True, exist_ok=True)
    FS_MOD_PF.mkdir(parents=True, exist_ok=True)

    keep = set(s.strip() for s in args.only.split(",") if s.strip())
    todo = [t for t in DMP_FRAMESHIFT if not keep or t[0] in keep]

    print(f"[plan] {len(todo)} variants x 2 stoich = {2*len(todo)} renders")
    for v, last_kept, tail_len in todo:
        render_one(v, "monomeric", FS_MOD_MO, MONO_VIEW,
                   dmp_chains="B", hap2_chains="A",
                   last_kept=last_kept, tail_len=tail_len, force=args.force)
        render_one(v, "postfusion_like", FS_MOD_PF, TRI_VIEW,
                   dmp_chains="D", hap2_chains="A,B,C",
                   last_kept=last_kept, tail_len=tail_len, force=args.force)

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
