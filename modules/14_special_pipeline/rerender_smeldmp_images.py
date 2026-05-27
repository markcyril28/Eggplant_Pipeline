#!/usr/bin/env python3
"""
Re-render every PSE session file under SmelDMP_variants/ to refresh its
PNG output, re-crop, re-label red panels, and rebuild the four grouped
figures (deletion mono/postfusion + frameshift mono/postfusion).

Use case: the user re-oriented one or more .pse files in PyMOL and wants
the corresponding PNGs (and the grouped panels that aggregate them)
refreshed without re-running the full HAP2_variants subtree or
re-rendering domain_map.

Scope (this script touches ONLY SmelDMP_variants/<sub>/...):
  deletion_red_monomeric/        delN / delC / delTMDcore
  deletion_red_postfusion/
  model_monomeric/
  model_postfusion/
  frameshift_red_monomeric/      fsGuide{4,16,17,20,37,46,50}
  frameshift_red_postfusion/
  frameshift_model_monomeric/
  frameshift_model_postfusion/

For each .pse we (re-)produce, per folder convention:
  red folders:    <base>.pse -> <base>.png -> <base>_cropped.png
                                -> <base>_labeled.png  (overlay_red_legend)
  model folders:  <base>.pse -> <base>_structure.png -> <base>_structure_cropped.png

Idempotent w.r.t. session files: each PSE renders the view stored inside.
PNG/labeled targets are overwritten.

Finally calls build_smeldmp_grouped_panels.py to refresh the four
grouped figures from the new labeled/cropped PNGs, plus the side-by-side
deletion composite (smeldmp_deletion_variants_composite.png) that pairs
Monomeric (left) with Trimeric / post-fusion (right) row-by-row.
"""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys

PYMOL_ENV = r"C:\ProgramData\pymol"
PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
DMP_DIR = PROJECT / ("III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/"
                    "deletion_ladder/07_Summary/deletion_variants/SmelDMP_variants")
RENDER_PSE = PROJECT / "modules/14_special_pipeline/render_from_pse.py"
CROP = PROJECT / "modules/14_special_pipeline/crop_to_content.py"
OVERLAY = PROJECT / "modules/14_special_pipeline/overlay_red_legend.py"
PANEL_BUILDER = PROJECT / "modules/14_special_pipeline/build_smeldmp_grouped_panels.py"

# Mirror of the relevant tables from refresh_smeldmp_panels.py.
# (variant, mono_iptm, pf_iptm, SmelDMP-deleted-range, AtDMP-deleted-range,
#  description, frameshift_tail_length [0 for deletion variants])
DMP_DELETION = [
    ("delN",       0.18, 0.52, "1-83",    "1-92",    "N-terminal cytoplasmic / specificity region", 0),
    ("delC",       0.19, 0.49, "221-222", "242-243", "C-terminal cytoplasmic tail",                 0),
    ("delTMDcore", 0.14, 0.51, "84-220",  "93-241",  "four-TMD bundle (MEME-1 to MEME-4)",          0),
]
# (variant, mono_iptm, pf_iptm, last_kept_aa, tail_len, description)
DMP_FRAMESHIFT = [
    ("fsGuide4",  0.47, 0.56, 7,   3,  "+1 NHEJ frameshift @ aa 7 (immediate post-ATG; null)"),
    ("fsGuide17", 0.23, 0.48, 118, 23, "+1 NHEJ frameshift @ aa 118 (TMD1-TMD2 linker)"),
    ("fsGuide16", 0.23, 0.48, 117, 24, "+1 NHEJ frameshift @ aa 118 (1 bp upstream of Guide17)"),
    ("fsGuide20", 0.26, 0.52, 128, 13, "+1 NHEJ frameshift @ aa 129 (TMD2 onset)"),
    ("fsGuide37", 0.18, 0.48, 215, 25, "+1 NHEJ frameshift @ aa 216 (near WT stop; near-WT length)"),
    ("fsGuide46", 0.22, 0.56, 27,  15, "+1 NHEJ frameshift @ aa 27 (early N-term; severe)"),
    ("fsGuide50", 0.23, 0.53, 80,  3,  "+1 NHEJ frameshift @ aa 80 (pre-TMD1)"),
]

# AF3 baselines for the WT__WT reference (used by classify()).
WT_MONO = 0.22
WT_POSTFUSION = 0.52
SHORT_CHAIN_SELF_IPTM = 0.10


def classify(iptm: float, stoich: str, dmp_self_iptm: float | None = None) -> str:
    if dmp_self_iptm is not None and dmp_self_iptm < SHORT_CHAIN_SELF_IPTM:
        return "null (short chain)"
    if iptm <= 0:
        return "n/a"
    baseline = WT_MONO if stoich == "monomeric" else WT_POSTFUSION
    frac = iptm / baseline
    if frac >= 0.80: return "tolerated"
    if frac >= 0.50: return "reduced"
    if frac >= 0.25: return "strongly reduced"
    return "catastrophic"


def read_dmp_self_iptm(pair_dir: Path, stoich: str) -> float | None:
    import json
    js = pair_dir / f"fold_{pair_dir.name}_{stoich}_summary_confidences_0.json"
    if not js.exists():
        return None
    try:
        data = json.loads(js.read_text())
    except (OSError, ValueError):
        return None
    mat = data.get("chain_pair_iptm")
    if not isinstance(mat, list) or not mat or not isinstance(mat[-1], list) or not mat[-1]:
        return None
    try:
        return float(mat[-1][-1])
    except (TypeError, ValueError):
        return None


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
    return r.returncode


def render_from_pse(pse: Path, png: Path) -> bool:
    if png.exists():
        png.unlink()
    rc = run_pymol(str(RENDER_PSE), "--", "--pse", str(pse), "--out", str(png))
    return rc == 0 and png.exists()


def crop(src: Path, dst: Path) -> bool:
    if dst.exists():
        dst.unlink()
    rc = run_py(str(CROP), "--image", str(src), "--out", str(dst), "--pad", "12")
    return rc == 0 and dst.exists()


def label_red(cropped: Path, labeled: Path, label: str, iptm: float, cat: str) -> bool:
    if labeled.exists():
        labeled.unlink()
    rc = run_py(str(OVERLAY),
                "--image", str(cropped), "--out", str(labeled),
                "--label", label,
                "--iptm", f"{iptm:.2f}", "--iptm-category", cat)
    return rc == 0 and labeled.exists()


def rerender_red(folder: Path, variant: str, red_pat: str, label: str,
                 iptm: float, category: str) -> None:
    """For red panels: render <base>.png from .pse, crop, label."""
    pse_name = red_pat.replace(".png", ".pse")
    base = folder / red_pat[:-len(".png")]  # strip .png
    pse = folder / pse_name
    png = folder / red_pat
    cropped = base.with_name(base.name + "_cropped.png")
    labeled = base.with_name(base.name + "_labeled.png")
    if not pse.exists():
        print(f"  [SKIP] {variant} ({folder.name}): no .pse")
        return
    if not render_from_pse(pse, png):
        print(f"  [FAIL] {variant} ({folder.name}): render failed")
        return
    if not crop(png, cropped):
        print(f"  [FAIL] {variant} ({folder.name}): crop failed")
        return
    if not label_red(cropped, labeled, label, iptm, category):
        print(f"  [FAIL] {variant} ({folder.name}): label failed")
        return
    print(f"  [OK]   {variant} ({folder.name}): {labeled.name}")


def rerender_model(folder: Path, variant: str, base_name: str) -> None:
    """For model panels: render <base>_structure.png from .pse, crop."""
    pse = folder / f"{base_name}.pse"
    structure = folder / f"{base_name}_structure.png"
    cropped = folder / f"{base_name}_structure_cropped.png"
    if not pse.exists():
        print(f"  [SKIP] {variant} ({folder.name}): no .pse")
        return
    if not render_from_pse(pse, structure):
        print(f"  [FAIL] {variant} ({folder.name}): render failed")
        return
    if not crop(structure, cropped):
        print(f"  [FAIL] {variant} ({folder.name}): crop failed")
        return
    print(f"  [OK]   {variant} ({folder.name}): {cropped.name}")


def main() -> int:
    DEL_RED_MO = DMP_DIR / "deletion_red_monomeric"
    DEL_RED_PF = DMP_DIR / "deletion_red_postfusion"
    DEL_MOD_MO = DMP_DIR / "model_monomeric"
    DEL_MOD_PF = DMP_DIR / "model_postfusion"
    FS_RED_MO  = DMP_DIR / "frameshift_red_monomeric"
    FS_RED_PF  = DMP_DIR / "frameshift_red_postfusion"
    FS_MOD_MO  = DMP_DIR / "frameshift_model_monomeric"
    FS_MOD_PF  = DMP_DIR / "frameshift_model_postfusion"

    COMPLEXES = PROJECT / ("III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/"
                           "deletion_ladder/02_Complexes")

    # ---- Deletion variants (red panels: mono + postfusion) ----
    print("\n## Deletion red panels ##")
    for v, iptm_mo, iptm_pf, sm, at, descr, _tail in DMP_DELETION:
        label = f"Red: {v} (SmelDMP {sm} / AtDMP8 {at}) -- {descr}"
        rerender_red(DEL_RED_MO, v,
                     red_pat=f"hap2_dmp_{v}_red_monomeric.png",
                     label=label, iptm=iptm_mo,
                     category=classify(iptm_mo, "monomeric"))
        rerender_red(DEL_RED_PF, v,
                     red_pat=f"hap2_dmp_{v}_red.png",
                     label=label, iptm=iptm_pf,
                     category=classify(iptm_pf, "postfusion_like"))

    # ---- Deletion variants (model panels: mono + postfusion) ----
    print("\n## Deletion model panels ##")
    for v, *_ in DMP_DELETION:
        rerender_model(DEL_MOD_MO, v, base_name=f"hap2_dmp_{v}_model_monomeric")
        rerender_model(DEL_MOD_PF, v, base_name=f"hap2_dmp_{v}_model_postfusion")

    # ---- Frameshift variants (red panels: mono + postfusion) ----
    print("\n## Frameshift red panels ##")
    for v, iptm_mo, iptm_pf, last_kept, tail_len, descr in DMP_FRAMESHIFT:
        label = f"Red: {v} (lost SmelDMP {last_kept+1}-222; +{tail_len} aa novel tail) -- {descr}"
        self_mo = read_dmp_self_iptm(COMPLEXES / f"monomeric/WT__{v}", "monomeric")
        self_pf = read_dmp_self_iptm(COMPLEXES / f"postfusion_like/WT__{v}", "postfusion_like")
        rerender_red(FS_RED_MO, v,
                     red_pat=f"hap2_dmp_{v}_red_monomeric.png",
                     label=label, iptm=iptm_mo,
                     category=classify(iptm_mo, "monomeric", self_mo))
        rerender_red(FS_RED_PF, v,
                     red_pat=f"hap2_dmp_{v}_red.png",
                     label=label, iptm=iptm_pf,
                     category=classify(iptm_pf, "postfusion_like", self_pf))

    # ---- Frameshift variants (model panels: mono + postfusion) ----
    print("\n## Frameshift model panels ##")
    for v, *_ in DMP_FRAMESHIFT:
        rerender_model(FS_MOD_MO, v, base_name=f"hap2_dmp_{v}_model_monomeric")
        rerender_model(FS_MOD_PF, v, base_name=f"hap2_dmp_{v}_model_postfusion")

    # ---- Rebuild the four grouped panels ----
    print("\n## Rebuilding grouped panels ##")
    run_py(str(PANEL_BUILDER))

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
