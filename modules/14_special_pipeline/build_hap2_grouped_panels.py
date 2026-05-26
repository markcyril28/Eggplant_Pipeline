#!/usr/bin/env python3
"""
Build "Single domain deletions only" and "Combined domain deletions only"
versions of the HAP2_variants grouped panel, with per-panel captions
("{variant} (red)" / "{variant} (model)") so the variant name appears
next to each subpanel letter.

Output files:
  hap2_dmp_ectodomain_deletions_grouped_monomeric_prefusion_{single,combined}.png
  hap2_dmp_ectodomain_deletions_grouped_trimeric_postfusion_{single,combined}.png

Sort within each subset is by ipTM ascending (worst-first).
"""
from pathlib import Path
import subprocess, sys

PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
HAP2_DIR = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder/07_Summary/deletion_variants/HAP2_variants"
RED_MO = HAP2_DIR / "deletion_red_monomeric"
RED_PF = HAP2_DIR / "deletion_red_postfusion"
MOD_MO = HAP2_DIR / "model_monomeric"
MOD_PF = HAP2_DIR / "model_postfusion"
COMBINE = PROJECT / "modules/14_special_pipeline/combine_panels.py"

HAP2_IPTM = {
    "delEcto":                    (0.12, 0.13),
    "delEctoAndC":                (0.17, 0.14),
    "delPreTMDAndTMDAndJuxtaMem": (0.17, 0.49),
    "delPreTMD":                  (0.19, 0.48),
    "delFL":                      (0.20, 0.43),
    "delPreTMDAndTMD":            (0.20, 0.48),
    "delEctoD2":                  (0.22, 0.38),
}
HAP2_SINGLE = ["delEcto", "delEctoD2", "delFL", "delPreTMD"]
HAP2_COMBO  = ["delEctoAndC", "delPreTMDAndTMD", "delPreTMDAndTMDAndJuxtaMem"]


def combine(out: Path, pairs, captions, panel_h=3.2, dpi=150):
    cmd = [sys.executable, str(COMBINE), "--out", str(out),
           "--cols", "2", "--panel-height", str(panel_h), "--dpi", str(dpi)]
    for p in pairs:
        cmd.extend(["--image", str(p)])
    for c in captions:
        cmd.extend(["--caption", c])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
    print(f"  {out.name}: rc={r.returncode}")


def pairs_and_caps_mo(names):
    pairs, caps = [], []
    for v in names:
        pairs.append(RED_MO / f"hap2_dmp_{v}_red_monomeric_labeled.png")
        caps.append(f"{v} (red)")
        pairs.append(MOD_MO / f"hap2_dmp_{v}_model_monomeric_structure_cropped.png")
        caps.append(f"{v} (model)")
    return pairs, caps

def pairs_and_caps_pf(names):
    pairs, caps = [], []
    for v in names:
        pairs.append(RED_PF / f"hap2_dmp_{v}_red_labeled.png")
        caps.append(f"{v} (red)")
        pairs.append(MOD_PF / f"hap2_dmp_{v}_model_postfusion_structure_cropped.png")
        caps.append(f"{v} (model)")
    return pairs, caps


# Single-domain only.
single_mo = sorted(HAP2_SINGLE, key=lambda v: HAP2_IPTM[v][0])
single_pf = sorted(HAP2_SINGLE, key=lambda v: HAP2_IPTM[v][1])
p, c = pairs_and_caps_mo(single_mo)
combine(HAP2_DIR / "hap2_dmp_ectodomain_deletions_grouped_monomeric_prefusion_single.png", p, c)
p, c = pairs_and_caps_pf(single_pf)
combine(HAP2_DIR / "hap2_dmp_ectodomain_deletions_grouped_trimeric_postfusion_single.png", p, c)

# Combined-domain only.
combo_mo = sorted(HAP2_COMBO, key=lambda v: HAP2_IPTM[v][0])
combo_pf = sorted(HAP2_COMBO, key=lambda v: HAP2_IPTM[v][1])
p, c = pairs_and_caps_mo(combo_mo)
combine(HAP2_DIR / "hap2_dmp_ectodomain_deletions_grouped_monomeric_prefusion_combined.png", p, c)
p, c = pairs_and_caps_pf(combo_pf)
combine(HAP2_DIR / "hap2_dmp_ectodomain_deletions_grouped_trimeric_postfusion_combined.png", p, c)

print("\n[DONE]")
