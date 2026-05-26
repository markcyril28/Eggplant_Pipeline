#!/usr/bin/env python3
"""
Rebuild the SmelDMP_variants grouped panels with per-panel captions, so each
subpanel letter sits next to its variant name and panel type (red/model),
matching the HAP2 grouped panel style:

    (A) fsGuide4 (red)
    (B) fsGuide4 (model)
    ...

Output files (overwrites in place):
  smeldmp_deletion_variants_monomeric.png      # delN / delC / delTMDcore
  smeldmp_deletion_variants_postfusion.png
  smeldmp_frameshift_variants_monomeric.png    # fsGuide{4,16,17,20,37,46,50}
  smeldmp_frameshift_variants_postfusion.png

Sort within each subset: ipTM ascending (worst-first), same as refresh_all_v2.py.
"""
from pathlib import Path
import subprocess
import sys

PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
DMP_DIR = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/07_Summary/deletion_variants/SmelDMP_variants"
# Deletion-variant assets (delN / delC / delTMDcore)
DEL_RED_MO = DMP_DIR / "deletion_red_monomeric"
DEL_RED_PF = DMP_DIR / "deletion_red_postfusion"
DEL_MOD_MO = DMP_DIR / "model_monomeric"
DEL_MOD_PF = DMP_DIR / "model_postfusion"
# Frameshift-variant assets (fsGuide*)
FS_RED_MO  = DMP_DIR / "frameshift_red_monomeric"
FS_RED_PF  = DMP_DIR / "frameshift_red_postfusion"
FS_MOD_MO  = DMP_DIR / "frameshift_model_monomeric"
FS_MOD_PF  = DMP_DIR / "frameshift_model_postfusion"
COMBINE = PROJECT / "modules/14_special_pipeline/combine_panels.py"

# (variant, mono_iptm, postfusion_iptm) - mirrors refresh_all_v2.py.
DMP_DELETION_IPTM = {
    "delN":       (0.18, 0.52),
    "delC":       (0.19, 0.49),
    "delTMDcore": (0.14, 0.51),
}
DMP_FRAMESHIFT_IPTM = {
    "fsGuide4":  (0.47, 0.56),
    "fsGuide17": (0.23, 0.48),
    "fsGuide16": (0.23, 0.48),
    "fsGuide20": (0.26, 0.52),
    "fsGuide37": (0.18, 0.48),
    "fsGuide46": (0.22, 0.56),
    "fsGuide50": (0.23, 0.53),
}

# Variants present in DMP_FRAMESHIFT_IPTM but intentionally OMITTED from the
# smeldmp_frameshift_variants_*.png grouped panels. They are still scored in
# the ipTM heatmap (the orchestrator builds PAIR_LABEL from the TOML rows);
# only the figure panel is suppressed. fsGuide16 produces a frameshift that
# is structurally too similar to fsGuide17 and fsGuide20 (same TMD1-TMD2
# linker cut window, +1 NHEJ, 23-24 aa tail), so showing all three rows
# adds noise without adding information. Documented in
# 14_Interaction_Domain_Mapping_In_Silico_Experiments.md (Limitations).
_PANEL_EXCLUDE = {"fsGuide16"}

# CRISPR-P v2.0 on-target (Target) scores keyed by guide NUMBER -- shared with
# iptm_heatmap.CRISPR_P_V2_SCORES (kept in lockstep; update both when guides
# are re-scored). Source: Liu et al. (2017), CRISPR-P 2.0
# (http://crispr.hzau.edu.cn/CRISPR2/). Used to annotate the fsGuide panel
# captions: "(A) fsGuide4 (CP 0.12) (red)".
CRISPR_P_V2_SCORES: dict[int, float] = {
    4:  0.1172,
    16: 0.4235,
    17: 0.8725,
    20: 0.4861,
    37: 0.5374,
    46: 0.4683,
    50: 0.4767,
}


def fsguide_cp_suffix(variant: str) -> str:
    """Return ' (CP <score>)' for a fsGuide<N> variant, '' for anything else."""
    if not variant.startswith("fsGuide"):
        return ""
    try:
        num = int(variant[len("fsGuide"):])
    except ValueError:
        return ""
    score = CRISPR_P_V2_SCORES.get(num)
    return f" (CP {score:.2f})" if score is not None else ""


def combine(out: Path, pairs, captions, panel_h=3.2, dpi=150, cols=2,
            caption_fontsize: float | None = None):
    cmd = [sys.executable, str(COMBINE), "--out", str(out),
           "--cols", str(cols), "--panel-height", str(panel_h), "--dpi", str(dpi)]
    if caption_fontsize is not None:
        cmd.extend(["--caption-fontsize", str(caption_fontsize)])
    for p in pairs:
        cmd.extend(["--image", str(p)])
    for c in captions:
        cmd.extend(["--caption", c])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
    print(f"  {out.name}: rc={r.returncode}")


def interleave_4col(names, red_dir, mod_dir, red_pat, mod_pat):
    """Reorder N=6 variants for a 4-col (red,model,red,model) x 3-row grid where
    the left super-column holds names[0:3] (top-to-bottom) and the right
    super-column holds names[3:6] (top-to-bottom). combine_panels.py lays out
    images in row-major order, so the input list must be:
        [n0_red, n0_model, n3_red, n3_model,   # row 0
         n1_red, n1_model, n4_red, n4_model,   # row 1
         n2_red, n2_model, n5_red, n5_model]   # row 2
    """
    if len(names) != 6:
        raise ValueError(f"interleave_4col expects exactly 6 names, got {len(names)}")
    pairs, caps = [], []
    for row in range(3):
        for col_idx in (row, row + 3):
            v = names[col_idx]
            cp = fsguide_cp_suffix(v)   # " (CP 0.12)" for fsGuide* else ""
            pairs.append(red_dir / red_pat.format(v=v))
            caps.append(f"{v}{cp} (red)")
            pairs.append(mod_dir / mod_pat.format(v=v))
            caps.append(f"{v}{cp} (model)")
    return pairs, caps


def pairs_and_caps_mo(names, red_dir, mod_dir):
    pairs, caps = [], []
    for v in names:
        pairs.append(red_dir / f"hap2_dmp_{v}_red_monomeric_labeled.png")
        caps.append(f"{v} (red)")
        pairs.append(mod_dir / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
        caps.append(f"{v} (model)")
    return pairs, caps


def pairs_and_caps_pf(names, red_dir, mod_dir):
    pairs, caps = [], []
    for v in names:
        pairs.append(red_dir / f"hap2_dmp_{v}_red_labeled.png")
        caps.append(f"{v} (red)")
        pairs.append(mod_dir / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
        caps.append(f"{v} (model)")
    return pairs, caps


# Deletion subset (delN / delC / delTMDcore) - lives in deletion_*/model_*.
del_names = list(DMP_DELETION_IPTM.keys())
del_mo = sorted(del_names, key=lambda v: DMP_DELETION_IPTM[v][0])
del_pf = sorted(del_names, key=lambda v: DMP_DELETION_IPTM[v][1])
p, c = pairs_and_caps_mo(del_mo, DEL_RED_MO, DEL_MOD_MO)
combine(DMP_DIR / "smeldmp_deletion_variants_monomeric.png", p, c)
p, c = pairs_and_caps_pf(del_pf, DEL_RED_PF, DEL_MOD_PF)
combine(DMP_DIR / "smeldmp_deletion_variants_postfusion.png", p, c)

# Frameshift subset (fsGuide{4,16,17,20,37,46,50}) - lives in frameshift_*.
# _PANEL_EXCLUDE drops fsGuide16 from the figure (kept in the ipTM heatmap).
# Layout: 4-col x 3-row (two side-by-side super-columns, each holding 3
# variants as (red, model) pair-rows). interleave_4col() rearranges the
# 6 sorted names so combine_panels.py renders them in the right cells.
fs_names = [v for v in DMP_FRAMESHIFT_IPTM if v not in _PANEL_EXCLUDE]
fs_mo = sorted(fs_names, key=lambda v: DMP_FRAMESHIFT_IPTM[v][0])
fs_pf = sorted(fs_names, key=lambda v: DMP_FRAMESHIFT_IPTM[v][1])
# Smaller caption font here than the deletion panels above: the fsGuide
# captions now carry a CRISPR-P v2.0 score ("(A) fsGuide4 (CP 0.12) (red)"),
# which would otherwise overrun the narrow tile width at the default size.
FS_CAPTION_FONTSIZE = 5.5
p, c = interleave_4col(fs_mo, FS_RED_MO, FS_MOD_MO,
                       red_pat="hap2_dmp_{v}_red_monomeric_labeled.png",
                       mod_pat="hap2_dmp_{v}_model_monomeric_structure_cropped.png")
combine(DMP_DIR / "smeldmp_frameshift_variants_monomeric.png", p, c, cols=4,
        caption_fontsize=FS_CAPTION_FONTSIZE)
p, c = interleave_4col(fs_pf, FS_RED_PF, FS_MOD_PF,
                       red_pat="hap2_dmp_{v}_red_labeled.png",
                       mod_pat="hap2_dmp_{v}_model_postfusion_structure_cropped.png")
combine(DMP_DIR / "smeldmp_frameshift_variants_postfusion.png", p, c, cols=4,
        caption_fontsize=FS_CAPTION_FONTSIZE)

print("\n[DONE]")
