#!/usr/bin/env python3
"""
End-to-end refresh:
  1. HAP2_variants/ subtree -- existing HAP2 deletion variants
  2. SmelDMP_variants/ subtree -- DMP deletion + frameshift variants (NEW)
  3. domain_map/ -- Figure 0 (mono + trimer)

For every .pse already present, re-render its PNG from the saved view.
For variants that don't have a hand-tuned .pse yet, render fresh with the
WT Figure-0 view as fallback orientation.

Applies the new baseline-relative ipTM classification:
    >= 80% of WT baseline  -> tolerated
    50-80% of WT          -> reduced
    25-50% of WT          -> strongly reduced
    < 25% of WT           -> catastrophic
(WT baselines: monomeric 0.22, postfusion 0.52.)

Output structure:
  07_Summary/
    hap2_dmp_domain_map.png                              (Fig 0B, trimer)
    hap2_dmp_domain_map_monomeric.png                    (Fig 0A, mono)
    deletion_variants/
      domain_map/                                        (Fig 0 raw + .pse)
      HAP2_variants/
        deletion_red_{monomeric,postfusion}/             (.pse + variants)
        model_{monomeric,postfusion}/                    (.pse + variants)
        hap2_dmp_ectodomain_deletions_grouped_*.png      (combo top, single bottom)
      SmelDMP_variants/
        deletion_red_{monomeric,postfusion}/             (DMP red)
        model_{monomeric,postfusion}/                    (DMP variant models)
        smeldmp_deletion_variants_{monomeric,postfusion}.png  (delN/delC/delTMDcore)
        smeldmp_frameshift_variants_{monomeric,postfusion}.png  (fsGuide*)
"""
from pathlib import Path
import subprocess
import sys

PYMOL_ENV = r"C:\ProgramData\pymol"
PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
OUT = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/07_Summary"
DV = OUT / "deletion_variants"
DOMAIN_DIR = DV / "domain_map"
HAP2_DIR = DV / "HAP2_variants"
DMP_DIR  = DV / "SmelDMP_variants"
COMPLEXES = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/02_Complexes"

# Module scripts.
HAP2_SCRIPT = PROJECT / "modules/14_special_pipeline/render_hap2_domain_map.py"
RED_SCRIPT  = PROJECT / "modules/14_special_pipeline/render_deletion_red.py"
COMPOSE_LEG = PROJECT / "modules/14_special_pipeline/compose_legend.py"
CROP        = PROJECT / "modules/14_special_pipeline/crop_to_content.py"
OVERLAY     = PROJECT / "modules/14_special_pipeline/overlay_red_legend.py"
COMBINE     = PROJECT / "modules/14_special_pipeline/combine_panels.py"
RENDER_PSE  = PROJECT / "modules/14_special_pipeline/render_from_pse.py"
EXTRACT_VIEW = PROJECT / "modules/14_special_pipeline/extract_view.py"

# View JSONs (PyMOL camera orientation cache) live next to the script tree
# under modules/14_special_pipeline/.views/ so they travel with the repo.
VIEW_TMP = PROJECT / "modules/14_special_pipeline/.views"; VIEW_TMP.mkdir(exist_ok=True)
MONO_VIEW = VIEW_TMP / "fig0_monomeric.json"
TRI_VIEW  = VIEW_TMP / "fig0_trimer.json"

# ipTM classification (baseline-relative).
WT_MONO = 0.22
WT_POSTFUSION = 0.52
# Below this DMP self-ipTM (chain_pair_iptm[-1][-1]) the DMP chain has no
# internal structure confidence (typical for very short frameshift peptides,
# e.g. fsGuide4 = 10 aa with B-self = 0.01). The overall ipTM in that case
# is a point-dock artifact, not an interface-affinity signal, so we mark
# the variant as "null (short chain)" instead of classifying it on the
# baseline-relative scale (where 0.58 would otherwise read as "tolerated").
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
    """Read chain_pair_iptm[-1][-1] from the rank-0 summary JSON for this pair.
    Returns None if the JSON is missing or the field cannot be parsed.
    """
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


# HAP2 variants and metadata.
HAP2_VARIANTS = ["delEctoAndC", "delEcto", "delEctoD2", "delFL", "delPreTMD",
                 "delPreTMDAndTMD", "delPreTMDAndTMDAndJuxtaMem"]
HAP2_IPTM = {
    "delEcto":                    (0.12, 0.13),
    "delEctoAndC":                (0.17, 0.14),
    "delPreTMDAndTMDAndJuxtaMem": (0.17, 0.49),
    "delPreTMD":                  (0.19, 0.48),
    "delFL":                      (0.20, 0.43),
    "delPreTMDAndTMD":            (0.20, 0.48),
    "delEctoD2":                  (0.22, 0.38),
}
HAP2_LABEL_META = {
    "delEcto":                    ("22-589",            "66-425",            "HAP2 ectodomain", "22-589"),
    "delEctoAndC":                ("22-589 + 655-804",  "66-425 + 596-705",  "ectodomain + C-tail", "22-589,655-804"),
    "delEctoD2":                  ("186-325",           "66-304",            "ectodomain D2 subdomain", "186-325"),
    "delFL":                      ("231-247",           "280-295",           "fusion loop / cd-loop", "231-247"),
    "delPreTMD":                  ("590-619",           "531-560",           "pre-TMD linker / stem base", "590-619"),
    "delPreTMDAndTMD":            ("590-641",           "531-582",           "pre-TMD linker + TMD", "590-641"),
    "delPreTMDAndTMDAndJuxtaMem": ("590-654",           "531-595",           "pre-TMD linker + TMD + juxtamembrane", "590-654"),
}
HAP2_SINGLE = ["delEcto", "delEctoD2", "delFL", "delPreTMD"]
HAP2_COMBO  = ["delEctoAndC", "delPreTMDAndTMD", "delPreTMDAndTMDAndJuxtaMem"]

# SmelDMP variants and metadata.
# (variant, mono_iptm, pf_iptm, SmelDMP-delete-range, AtDMP-delete-range,
#  descr, frameshift_tail_length [0 for deletion variants])
DMP_DELETION = [
    ("delN",       0.18, 0.52, "1-83",    "1-92",    "N-terminal cytoplasmic / specificity region", 0),
    ("delC",       0.19, 0.49, "221-222", "242-243", "C-terminal cytoplasmic tail",                 0),
    ("delTMDcore", 0.14, 0.51, "84-220",  "93-241",  "four-TMD bundle (MEME-1 to MEME-4)",          0),
]
# (variant, mono_iptm, pf_iptm, last_kept_aa, tail_len, descr)
DMP_FRAMESHIFT = [
    ("fsGuide4",  0.47, 0.56, 7,   3,  "+1 NHEJ frameshift @ aa 7 (immediate post-ATG; null)"),
    ("fsGuide17", 0.23, 0.48, 118, 23, "+1 NHEJ frameshift @ aa 118 (TMD1-TMD2 linker)"),
    ("fsGuide16", 0.23, 0.48, 117, 24, "+1 NHEJ frameshift @ aa 118 (1 bp upstream of Guide17)"),
    ("fsGuide20", 0.26, 0.52, 128, 13, "+1 NHEJ frameshift @ aa 129 (TMD2 onset)"),
    ("fsGuide37", 0.18, 0.48, 215, 25, "+1 NHEJ frameshift @ aa 216 (near WT stop; near-WT length)"),
    ("fsGuide46", 0.22, 0.56, 27,  15, "+1 NHEJ frameshift @ aa 27 (early N-term; severe)"),
    ("fsGuide50", 0.23, 0.53, 80,  3,  "+1 NHEJ frameshift @ aa 80 (pre-TMD1)"),
]


def run_pymol(*args):
    cmd = ["conda", "run", "-p", PYMOL_ENV, "pymol", "-cq", *args]
    res = subprocess.run(cmd, cwd=PROJECT, capture_output=True, text=True)
    for line in (res.stdout + res.stderr).splitlines():
        if "[OK]" in line or "[INFO]" in line or "[WARN]" in line or "[ERROR]" in line:
            print("  ", line)
    return res

def run_py(*args):
    res = subprocess.run([sys.executable, *args], cwd=PROJECT, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout); print(res.stderr)
    for line in (res.stdout + res.stderr).splitlines():
        if "[OK]" in line or "[WARN]" in line or "[ERROR]" in line or "cropped" in line:
            print(" ", line)


def render_from_pse_if_exists(pse: Path, png: Path):
    if pse.exists():
        run_pymol(str(RENDER_PSE), "--", "--pse", str(pse), "--out", str(png))


def crop(src: Path, dst: Path):
    if src.exists():
        run_py(str(CROP), "--image", str(src), "--out", str(dst), "--pad", "12")


def label_red(cropped: Path, labeled: Path, label: str, iptm: float, cat: str):
    run_py(str(OVERLAY),
           "--image", str(cropped), "--out", str(labeled),
           "--label", label,
           "--iptm", f"{iptm:.2f}", "--iptm-category", cat)


def combine_grid(out_path: Path, pairs: list[Path], panel_height: float = 2.6, dpi: int = 150):
    cmd = [sys.executable, str(COMBINE), "--out", str(out_path),
           "--cols", "2", "--panel-height", str(panel_height), "--dpi", str(dpi)]
    for p in pairs:
        cmd.extend(["--image", str(p)])
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout); print(res.stderr)
    print(f"  {out_path.name}: rc={res.returncode}")


# ---------- Step 1: extract Figure 0 views and refresh domain map ----------
print("\n## Figure 0 domain maps ##")
run_pymol(str(EXTRACT_VIEW), "--", "--pse", str(DOMAIN_DIR / "hap2_dmp_domain_map_monomeric.pse"),
          "--out", str(MONO_VIEW))
run_pymol(str(EXTRACT_VIEW), "--", "--pse", str(DOMAIN_DIR / "hap2_dmp_domain_map.pse"),
          "--out", str(TRI_VIEW))
render_from_pse_if_exists(DOMAIN_DIR / "hap2_dmp_domain_map_monomeric.pse",
                          DOMAIN_DIR / "hap2_dmp_domain_map_monomeric_structure.png")
render_from_pse_if_exists(DOMAIN_DIR / "hap2_dmp_domain_map.pse",
                          DOMAIN_DIR / "hap2_dmp_domain_map_structure.png")
crop(DOMAIN_DIR / "hap2_dmp_domain_map_monomeric_structure.png",
     DOMAIN_DIR / "hap2_dmp_domain_map_monomeric_structure_cropped.png")
crop(DOMAIN_DIR / "hap2_dmp_domain_map_structure.png",
     DOMAIN_DIR / "hap2_dmp_domain_map_structure_cropped.png")
run_py(str(COMPOSE_LEG),
       "--structure", str(DOMAIN_DIR / "hap2_dmp_domain_map_monomeric_structure_cropped.png"),
       "--out",       str(OUT / "hap2_dmp_domain_map_monomeric.png"),
       "--hap2-title", "SmelHAP2 (chain A) -- deletion-variant bands",
       "--dmp-title",  "SmelDMPv5_10.610 (chain B) -- DMP topology palette (08 config)")
run_py(str(COMPOSE_LEG),
       "--structure", str(DOMAIN_DIR / "hap2_dmp_domain_map_structure_cropped.png"),
       "--out",       str(OUT / "hap2_dmp_domain_map.png"))


# ---------- Step 2: HAP2_variants subtree -- refresh + rebuild grouped panels ----------
RED_MO = HAP2_DIR / "deletion_red_monomeric"
RED_PF = HAP2_DIR / "deletion_red_postfusion"
MOD_MO = HAP2_DIR / "model_monomeric"
MOD_PF = HAP2_DIR / "model_postfusion"

print("\n## HAP2_variants: refreshing renders from .pse ##")
for v in HAP2_VARIANTS:
    render_from_pse_if_exists(RED_MO / f"hap2_dmp_{v}_red_monomeric.pse",
                              RED_MO / f"hap2_dmp_{v}_red_monomeric.png")
    render_from_pse_if_exists(RED_PF / f"hap2_dmp_{v}_red.pse",
                              RED_PF / f"hap2_dmp_{v}_red.png")
    render_from_pse_if_exists(MOD_MO / f"hap2_dmp_{v}_model_monomeric.pse",
                              MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure.png")
    render_from_pse_if_exists(MOD_PF / f"hap2_dmp_{v}_model_postfusion.pse",
                              MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure.png")

print("\n## HAP2_variants: cropping ##")
for v in HAP2_VARIANTS:
    crop(RED_MO / f"hap2_dmp_{v}_red_monomeric.png",
         RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png")
    crop(RED_PF / f"hap2_dmp_{v}_red.png",
         RED_PF / f"hap2_dmp_{v}_red_cropped.png")
    crop(MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure.png",
         MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
    crop(MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure.png",
         MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")

print("\n## HAP2_variants: labelling red panels (new ipTM classifications) ##")
for v in HAP2_VARIANTS:
    sm, at, descr, _ = HAP2_LABEL_META[v]
    iptm_mo, iptm_pf = HAP2_IPTM[v]
    label = f"Red: {v} (SmelHAP2 {sm} / AtHAP2 {at}) -- {descr}"
    label_red(RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png",
              RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png",
              label, iptm_mo, classify(iptm_mo, "monomeric"))
    label_red(RED_PF / f"hap2_dmp_{v}_red_cropped.png",
              RED_PF / f"hap2_dmp_{v}_red_labeled.png",
              label, iptm_pf, classify(iptm_pf, "postfusion_like"))

print("\n## HAP2_variants: rebuilding grouped panels ##")
# The manuscript references only the split combined+single grouped panels
# (built below via split_hap2_grouped.py). The "all" version stays around
# as reference but is written into _legacy/ to keep the working folder clean.
LEGACY_DIR = OUT / "_legacy"
LEGACY_DIR.mkdir(exist_ok=True)
mo_grouped = sorted(HAP2_COMBO, key=lambda v: HAP2_IPTM[v][0]) + sorted(HAP2_SINGLE, key=lambda v: HAP2_IPTM[v][0])
pf_grouped = sorted(HAP2_COMBO, key=lambda v: HAP2_IPTM[v][1]) + sorted(HAP2_SINGLE, key=lambda v: HAP2_IPTM[v][1])
mo_pairs = []
for v in mo_grouped:
    mo_pairs.append(RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png")
    mo_pairs.append(MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
combine_grid(LEGACY_DIR / "hap2_dmp_ectodomain_deletions_grouped_monomeric.png", mo_pairs)
pf_pairs = []
for v in pf_grouped:
    pf_pairs.append(RED_PF / f"hap2_dmp_{v}_red_labeled.png")
    pf_pairs.append(MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
combine_grid(LEGACY_DIR / "hap2_dmp_ectodomain_deletions_grouped_postfusion.png", pf_pairs)

# Sorted-by-ipTM legacy panels. Written into 07_Summary/_legacy/ alongside
# the "all" grouped versions above; the manuscript references only the
# split combined/single grouped panels under HAP2_variants/.
mo_sort = sorted(HAP2_VARIANTS, key=lambda v: HAP2_IPTM[v][0])
pf_sort = sorted(HAP2_VARIANTS, key=lambda v: HAP2_IPTM[v][1])
mo_pairs = []
for v in mo_sort:
    mo_pairs.append(RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png")
    mo_pairs.append(MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
combine_grid(LEGACY_DIR / "hap2_dmp_ectodomain_deletions_panel_monomeric.png", mo_pairs)
pf_pairs = []
for v in pf_sort:
    pf_pairs.append(RED_PF / f"hap2_dmp_{v}_red_labeled.png")
    pf_pairs.append(MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
combine_grid(LEGACY_DIR / "hap2_dmp_ectodomain_deletions_panel_postfusion.png", pf_pairs)


# ---------- Step 3: SmelDMP_variants subtree (NEW) ----------
# Deletion-variant assets (delN / delC / delTMDcore) live in the original
# `deletion_*` and `model_*` folders. Frameshift-variant assets (fsGuide*)
# now live in their own parallel sibling folders so the two variant classes
# are easy to find separately. The PSE session files and rendered PNGs for
# each fsGuide variant live alongside each other in:
#   frameshift_red_monomeric/    <- WT__WT red-highlight (mono) of lost residues
#   frameshift_red_postfusion/   <- WT__WT red-highlight (trimer) of lost residues
#   frameshift_model_monomeric/  <- AF3 WT__fsGuide* model (mono)
#   frameshift_model_postfusion/ <- AF3 WT__fsGuide* model (trimer)
DMP_RED_MO = DMP_DIR / "deletion_red_monomeric"
DMP_RED_PF = DMP_DIR / "deletion_red_postfusion"
DMP_MOD_MO = DMP_DIR / "model_monomeric"
DMP_MOD_PF = DMP_DIR / "model_postfusion"
FS_RED_MO  = DMP_DIR / "frameshift_red_monomeric"
FS_RED_PF  = DMP_DIR / "frameshift_red_postfusion"
FS_MOD_MO  = DMP_DIR / "frameshift_model_monomeric"
FS_MOD_PF  = DMP_DIR / "frameshift_model_postfusion"
for d in [DMP_RED_MO, DMP_RED_PF, DMP_MOD_MO, DMP_MOD_PF,
          FS_RED_MO, FS_RED_PF, FS_MOD_MO, FS_MOD_PF]:
    d.mkdir(parents=True, exist_ok=True)

# Render deletion + frameshift DMP variants. For each variant:
#   * WT red highlight on the WT__WT complex (target=dmp); the WT-DMP CIF
#     residues highlight the deletion-range or "post-cut to end" frame.
#   * AF3 variant model with WT-aware DMP topology + frameshift-tail marker.
print("\n## SmelDMP_variants: rendering DMP variants (this populates fresh .pse files) ##")
WT_MONO_CIF = COMPLEXES / "monomeric/WT__WT/fold_WT__WT_monomeric_model_0.cif"
WT_PF_CIF   = COMPLEXES / "postfusion_like/WT__WT/fold_WT__WT_postfusion_like_model_0.cif"

def dmp_variant_cif(v: str, stoich: str) -> Path:
    return COMPLEXES / f"{stoich}/WT__{v}/fold_WT__{v}_{stoich}_model_0.cif"

# Convert deletion specs:
#   Deletion variants: explicit ranges given.
#   Frameshift variants: deletion = (last_kept+1)..222 (residues LOST from WT).
def deletion_range(v_entry):
    if len(v_entry) == 7:  # deletion tuple: (v, mono, pf, sm, at, descr, tail)
        return v_entry[3]
    # frameshift: (v, mono, pf, last_kept, tail, descr)
    last_kept = v_entry[3]
    return f"{last_kept+1}-222"

# Render fn -- if a .pse already exists at the target path, prefer rendering
# from it (preserves the user's hand-tuned orientation); otherwise render
# from the CIF using the WT Figure-0 view.
def render_red_dmp(variant: str, stoich: str, wt_cif: Path, view_json: Path, out_png: Path, resi: str):
    pse = out_png.with_suffix(".pse")
    if pse.exists():
        render_from_pse_if_exists(pse, out_png)
    else:
        run_pymol(str(RED_SCRIPT), "--",
                  "--cif", str(wt_cif),
                  "--residues", resi,
                  "--target", "dmp",
                  "--out", str(out_png),
                  "--view-json", str(view_json))

def render_model_dmp(variant_cif: Path, stoich: str, view_json: Path,
                     out_png: Path, dmp_chains: str, hap2_chains: str,
                     dmp_deleted: str, tail_len: int):
    pse = out_png.with_suffix(".pse")
    if pse.exists():
        struct = out_png.with_name(out_png.stem + "_structure.png")
        render_from_pse_if_exists(pse, struct)
    else:
        args = [str(HAP2_SCRIPT), "--",
                "--cif", str(variant_cif),
                "--out", str(out_png),
                "--hap2-chains", hap2_chains,
                "--dmp-chains", dmp_chains,
                "--view-json", str(view_json)]
        if dmp_deleted:
            args.extend(["--dmp-deleted-residues", dmp_deleted])
        if tail_len > 0:
            args.extend(["--dmp-novel-tail", str(tail_len)])
        run_pymol(*args)

# Deletion variants
print("\n## DMP deletion variants ##")
for entry in DMP_DELETION:
    v, iptm_mo, iptm_pf, sm, at, descr, _tail = entry
    # WT red highlights (mono, pf)
    render_red_dmp(v, "monomeric", WT_MONO_CIF, MONO_VIEW,
                   DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric.png", sm)
    render_red_dmp(v, "postfusion_like", WT_PF_CIF, TRI_VIEW,
                   DMP_RED_PF / f"hap2_dmp_{v}_red.png", sm)
    # Variant models (mono = chain A HAP2 + chain B DMP; postfusion = A,B,C + D)
    render_model_dmp(dmp_variant_cif(v, "monomeric"), "monomeric", MONO_VIEW,
                     DMP_MOD_MO / f"hap2_dmp_{v}_model_monomeric.png",
                     dmp_chains="B", hap2_chains="A",
                     dmp_deleted=sm, tail_len=0)
    render_model_dmp(dmp_variant_cif(v, "postfusion_like"), "postfusion_like", TRI_VIEW,
                     DMP_MOD_PF / f"hap2_dmp_{v}_model_postfusion.png",
                     dmp_chains="D", hap2_chains="A,B,C",
                     dmp_deleted=sm, tail_len=0)

# Frameshift variants - write to the dedicated frameshift_* siblings so the
# fsGuide assets stay separate from the delN/delC/delTMDcore deletions.
print("\n## DMP frameshift variants ##")
for entry in DMP_FRAMESHIFT:
    v, iptm_mo, iptm_pf, last_kept, tail_len, descr = entry
    resi_lost = f"{last_kept+1}-222"
    render_red_dmp(v, "monomeric", WT_MONO_CIF, MONO_VIEW,
                   FS_RED_MO / f"hap2_dmp_{v}_red_monomeric.png", resi_lost)
    render_red_dmp(v, "postfusion_like", WT_PF_CIF, TRI_VIEW,
                   FS_RED_PF / f"hap2_dmp_{v}_red.png", resi_lost)
    # For frameshift models: deleted-residues = resi_lost (post-cut to 222),
    # but the variant chain has the cut-1 retained WT + appended tail.
    render_model_dmp(dmp_variant_cif(v, "monomeric"), "monomeric", MONO_VIEW,
                     FS_MOD_MO / f"hap2_dmp_{v}_model_monomeric.png",
                     dmp_chains="B", hap2_chains="A",
                     dmp_deleted=resi_lost, tail_len=tail_len)
    render_model_dmp(dmp_variant_cif(v, "postfusion_like"), "postfusion_like", TRI_VIEW,
                     FS_MOD_PF / f"hap2_dmp_{v}_model_postfusion.png",
                     dmp_chains="D", hap2_chains="A,B,C",
                     dmp_deleted=resi_lost, tail_len=tail_len)

print("\n## SmelDMP_variants: cropping ##")
# Deletion variants crop in-place under deletion_*/model_*; frameshift
# variants crop in-place under frameshift_*/.
for entry in DMP_DELETION:
    v = entry[0]
    crop(DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric.png",
         DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png")
    crop(DMP_RED_PF / f"hap2_dmp_{v}_red.png",
         DMP_RED_PF / f"hap2_dmp_{v}_red_cropped.png")
    crop(DMP_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure.png",
         DMP_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
    crop(DMP_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure.png",
         DMP_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
for entry in DMP_FRAMESHIFT:
    v = entry[0]
    crop(FS_RED_MO / f"hap2_dmp_{v}_red_monomeric.png",
         FS_RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png")
    crop(FS_RED_PF / f"hap2_dmp_{v}_red.png",
         FS_RED_PF / f"hap2_dmp_{v}_red_cropped.png")
    crop(FS_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure.png",
         FS_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
    crop(FS_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure.png",
         FS_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")

print("\n## SmelDMP_variants: labelling red panels ##")
for entry in DMP_DELETION:
    v, iptm_mo, iptm_pf, sm, at, descr, _ = entry
    label = f"Red: {v} (SmelDMP {sm} / AtDMP8 {at}) -- {descr}"
    label_red(DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png",
              DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png",
              label, iptm_mo, classify(iptm_mo, "monomeric"))
    label_red(DMP_RED_PF / f"hap2_dmp_{v}_red_cropped.png",
              DMP_RED_PF / f"hap2_dmp_{v}_red_labeled.png",
              label, iptm_pf, classify(iptm_pf, "postfusion_like"))

for entry in DMP_FRAMESHIFT:
    v, iptm_mo, iptm_pf, last_kept, tail_len, descr = entry
    label = f"Red: {v} (lost SmelDMP {last_kept+1}-222; +{tail_len} aa novel tail) -- {descr}"
    # Frameshift peptides can be short enough that AF3 cannot fold the chain
    # internally; in that case the overall ipTM is a point-dock artifact and
    # classify() should return "null (short chain)" instead of a baseline-
    # relative bucket. The DMP self-ipTM comes from chain_pair_iptm[-1][-1]
    # in the rank-0 AF3 summary JSON for the WT__<variant> pair.
    pair_mo = COMPLEXES / f"monomeric/WT__{v}"
    pair_pf = COMPLEXES / f"postfusion_like/WT__{v}"
    self_mo = read_dmp_self_iptm(pair_mo, "monomeric")
    self_pf = read_dmp_self_iptm(pair_pf, "postfusion_like")
    label_red(FS_RED_MO / f"hap2_dmp_{v}_red_monomeric_cropped.png",
              FS_RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png",
              label, iptm_mo, classify(iptm_mo, "monomeric", self_mo))
    label_red(FS_RED_PF / f"hap2_dmp_{v}_red_cropped.png",
              FS_RED_PF / f"hap2_dmp_{v}_red_labeled.png",
              label, iptm_pf, classify(iptm_pf, "postfusion_like", self_pf))

print("\n## SmelDMP_variants: building grouped panels ##")
DMP_DEL_NAMES = [e[0] for e in DMP_DELETION]
DMP_FS_NAMES  = [e[0] for e in DMP_FRAMESHIFT]
DMP_IPTM = {e[0]: (e[1], e[2]) for e in DMP_DELETION + [(e[0], e[1], e[2], None, None, None) for e in DMP_FRAMESHIFT]}

def dmp_sort_by_iptm(names, idx):
    return sorted(names, key=lambda v: DMP_IPTM[v][idx])

# Deletion subset (mono + postfusion).
del_mo = dmp_sort_by_iptm(DMP_DEL_NAMES, 0)
pairs = []
for v in del_mo:
    pairs.append(DMP_RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png")
    pairs.append(DMP_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
combine_grid(DMP_DIR / "smeldmp_deletion_variants_monomeric.png", pairs)

del_pf = dmp_sort_by_iptm(DMP_DEL_NAMES, 1)
pairs = []
for v in del_pf:
    pairs.append(DMP_RED_PF / f"hap2_dmp_{v}_red_labeled.png")
    pairs.append(DMP_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
combine_grid(DMP_DIR / "smeldmp_deletion_variants_postfusion.png", pairs)

# Frameshift subset - read from the frameshift_* siblings.
fs_mo = dmp_sort_by_iptm(DMP_FS_NAMES, 0)
pairs = []
for v in fs_mo:
    pairs.append(FS_RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png")
    pairs.append(FS_MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
combine_grid(DMP_DIR / "smeldmp_frameshift_variants_monomeric.png", pairs)

fs_pf = dmp_sort_by_iptm(DMP_FS_NAMES, 1)
pairs = []
for v in fs_pf:
    pairs.append(FS_RED_PF / f"hap2_dmp_{v}_red_labeled.png")
    pairs.append(FS_MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
combine_grid(DMP_DIR / "smeldmp_frameshift_variants_postfusion.png", pairs)

print("\n[DONE]")
