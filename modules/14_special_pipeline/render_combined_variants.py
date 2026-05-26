#!/usr/bin/env python3
"""
Render the "special" combined-variant pairs (HAP2 deletion x DMP deletion
double-knockouts) into a dedicated Combined_Variants/ subtree under
07_Summary/deletion_variants/.

Special pairs come from [pairing.special_pairs].pairs in the orchestrator
TOML; they are pairs that sit OUTSIDE the orthogonal ladder (e.g. double
deletions that probe the cumulative HAP2-side + DMP-side knockout effect).
Default set hard-coded here mirrors the current TOML; pass --pairs to override.

Output layout (idempotent w.r.t. existing .pse files):
  Combined_Variants/
    model_monomeric/
      hap2_dmp_<pair>_model_monomeric.pse              (PyMOL session)
      hap2_dmp_<pair>_model_monomeric_structure.png    (rendered image)
      hap2_dmp_<pair>_model_monomeric_structure_cropped.png
    model_postfusion/
      hap2_dmp_<pair>_model_postfusion.{pse,_structure.png,_structure_cropped.png}
    combined_variants_monomeric.png      (grouped panel, captioned per variant)
    combined_variants_postfusion.png

Highlights BOTH deletions on the rendered model (HAP2 chain region in one
colour band, DMP chain region in another) via render_hap2_domain_map.py's
--deleted-residues + --dmp-deleted-residues flags. Uses the WT Figure-0
view (.views/fig0_{monomeric,trimer}.json) as the starting orientation;
re-open each .pse in PyMOL to reorient and re-save, then re-run this
script (or call render_from_pse.py directly) to refresh the PNG.

Run after canonicalise_af3_filenames.py so the CIFs are in the canonical
fold_{pair}_{stoich}_model_0.cif form.
"""
from __future__ import annotations
from pathlib import Path
import argparse
import subprocess
import sys

PYMOL_ENV = r"C:\ProgramData\pymol"
PROJECT = Path(r"c:\PIPELINE\Eggplant_Pipeline")
EXP_ROOT = PROJECT / "III_RESULT/DMP_x_SmelHAP2/14_Domain_Mapping/deletion_ladder"
COMPLEXES = EXP_ROOT / "02_Complexes"
COMBINED_DIR = EXP_ROOT / "07_Summary/deletion_variants/Combined_Variants"
RENDER_HAP2 = PROJECT / "modules/14_special_pipeline/render_hap2_domain_map.py"
RENDER_RED = PROJECT / "modules/14_special_pipeline/render_deletion_red.py"
RENDER_PSE = PROJECT / "modules/14_special_pipeline/render_from_pse.py"
EXTRACT_VIEW = PROJECT / "modules/14_special_pipeline/extract_view.py"
CROP = PROJECT / "modules/14_special_pipeline/crop_to_content.py"
COMBINE = PROJECT / "modules/14_special_pipeline/combine_panels.py"
OVERLAY = PROJECT / "modules/14_special_pipeline/overlay_red_legend.py"
VIEW_TMP = PROJECT / "modules/14_special_pipeline/.views"
MONO_VIEW = VIEW_TMP / "fig0_monomeric.json"
TRI_VIEW = VIEW_TMP / "fig0_trimer.json"

# ipTM classification baselines (WT__WT iptm; matches refresh_smeldmp_panels.py).
WT_MONO_IPTM = 0.22
WT_POSTFUSION_IPTM = 0.52


def classify_iptm(iptm: float, stoich: str) -> str:
    """Baseline-relative ipTM classification, mirroring
    refresh_smeldmp_panels.classify (without the short-chain branch -- DMP
    chains in combined pairs are >= 83 aa so DMP-self is never artifact-low).
    """
    if iptm <= 0:
        return "n/a"
    baseline = WT_MONO_IPTM if stoich == "monomeric" else WT_POSTFUSION_IPTM
    frac = iptm / baseline
    if frac >= 0.80: return "tolerated"
    if frac >= 0.50: return "reduced"
    if frac >= 0.25: return "strongly reduced"
    return "catastrophic"


def read_pair_iptm(pair: str, stoich: str) -> float | None:
    """Pull the overall iptm (rank-0) for one pair from the AF3 summary JSON."""
    import json
    js = (COMPLEXES / stoich / pair
          / f"fold_{pair}_{stoich}_summary_confidences_0.json")
    if not js.exists():
        return None
    try:
        return float(json.loads(js.read_text()).get("iptm"))
    except Exception:  # noqa: BLE001
        return None

# Default special pairs - mirror [pairing.special_pairs].pairs in the TOML.
# Each entry: (pair_label, hap2_deleted_residues, dmp_deleted_residues,
#              monomeric_iptm, postfusion_iptm, description)
# Coordinates are SmelHAP2 / SmelDMP. HAP2 delPreTMDAndTMD = 590-641 (from
# HAP2_LABEL_META in refresh_smeldmp_panels.py). DMP delN = 1-83, delTMDcore
# = 84-220 (from DMP_DELETION in the same file).
SPECIAL_PAIRS = [
    ("delPreTMDAndTMD__delN",
     "590-641", "1-83",
     "HAP2 pre-TMD+TMD removed AND DMP N-terminal cytoplasmic removed"),
    ("delPreTMDAndTMD__delTMDcore",
     "590-641", "84-220",
     "HAP2 pre-TMD+TMD removed AND DMP four-TMD bundle removed"),
]

# Human-readable region labels for the red-legend caption. Kept in lockstep
# with the residue ranges in SPECIAL_PAIRS above; AtHAP2 / AtDMP equivalents
# (in parens) come from the [hap2_variants.coords.at] / [dmp_variants.coords.at]
# tables in 14_Interaction_Domain_MappingCONFIG.toml.
HAP2_REGION_LABEL = {
    "delPreTMDAndTMD": "HAP2 pre-TMD + TMD (SmelHAP2 590-641 / AtHAP2 531-651)",
}
DMP_REGION_LABEL = {
    "delN":        "DMP N-terminal cytoplasmic (SmelDMP 1-83 / AtDMP 1-92)",
    "delTMDcore":  "DMP four-TMD bundle (SmelDMP 84-220 / AtDMP 93-241)",
}


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


def pair_cif(pair: str, stoich: str) -> Path:
    return COMPLEXES / stoich / pair / f"fold_{pair}_{stoich}_model_0.cif"


def render_pair(pair: str, stoich: str, hap2_del: str, dmp_del: str,
                out_dir: Path, force: bool) -> None:
    """Render the model PSE + structure PNG for one special pair."""
    suffix = "monomeric" if stoich == "monomeric" else "postfusion"
    base = out_dir / f"hap2_dmp_{pair}_model_{suffix}"
    pse = base.with_suffix(".pse")
    out_arg = base.with_suffix(".png")           # what we pass to --out
    structure = base.with_name(base.name + "_structure.png")
    cropped = base.with_name(base.name + "_structure_cropped.png")
    cif = pair_cif(pair, stoich)
    if not cif.exists():
        print(f"[SKIP] {pair} {stoich}: missing CIF {cif.name}")
        return

    if pse.exists() and not force:
        print(f"[KEEP] {pair} {stoich}: .pse exists; re-rendering PNG from it")
        if structure.exists():
            structure.unlink()
        run_pymol(str(RENDER_PSE), "--", "--pse", str(pse), "--out", str(structure))
    else:
        # Fresh CIF render. Stoich-specific chain layout:
        #   monomeric        = A (HAP2) + B (DMP)
        #   postfusion_like  = A,B,C (HAP2) + D (DMP)
        if stoich == "monomeric":
            hap2_chains, dmp_chains, view = "A", "B", MONO_VIEW
        else:
            hap2_chains, dmp_chains, view = "A,B,C", "D", TRI_VIEW
        args = [
            str(RENDER_HAP2), "--",
            "--cif", str(cif),
            "--out", str(out_arg),
            "--hap2-chains", hap2_chains,
            "--dmp-chains", dmp_chains,
            "--view-json", str(view),
            "--deleted-residues", hap2_del,
            "--dmp-deleted-residues", dmp_del,
        ]
        print(f"[RUN]  {pair} {stoich}: fresh render with HAP2 del={hap2_del} DMP del={dmp_del}")
        rc = run_pymol(*args)
        if rc != 0:
            print(f"[FAIL] {pair} {stoich}: PyMOL exit {rc}")
            return

    if structure.exists():
        if cropped.exists():
            cropped.unlink()
        run_py(str(CROP), "--image", str(structure),
               "--out", str(cropped), "--pad", "12")
        print(f"[OK]   {pair} {stoich}: {cropped.name}")


def render_pair_red(pair: str, stoich: str, hap2_del: str, dmp_del: str,
                    out_dir: Path, force: bool) -> Path | None:
    """Render the deletion-red variant (HAP2 + DMP zones both painted red) and
    its labeled overlay. Returns the labeled PNG path on success, None on skip.

    The red overlay is painted on the WT__WT complex CIF -- not the variant's
    own CIF -- because AF3 renumbers each variant chain starting from 1, so
    selecting `resi 590-641` on a delPreTMDAndTMD chain (which only has 752
    residues numbered 1-752) hits the wrong band. Painting on WT__WT keeps
    SmelHAP2 / SmelDMP residue numbering intact, so the red zones correspond
    to "what gets removed in this variant" -- matching the convention used
    by render_deletion_red.py for the HAP2_variants / SmelDMP_variants
    single-side red renders.

    Mirrors the HAP2_variants/SmelDMP_variants red-render layout so the
    Combined_Variants subtree can be grouped the same way:
      hap2_dmp_<pair>_red_<suffix>.pse
      hap2_dmp_<pair>_red_<suffix>.png             (raw PyMOL render)
      hap2_dmp_<pair>_red_<suffix>_cropped.png     (auto-cropped)
      hap2_dmp_<pair>_red_<suffix>_labeled.png     (legend + ipTM badge)
    """
    suffix = "monomeric" if stoich == "monomeric" else "postfusion"
    base = out_dir / f"hap2_dmp_{pair}_red_{suffix}"
    pse = base.with_suffix(".pse")
    raw = base.with_suffix(".png")
    cropped = base.with_name(base.name + "_cropped.png")
    labeled = base.with_name(base.name + "_labeled.png")

    if pse.exists() and not force:
        # User-edited .pse takes precedence: render straight from the session
        # so any colour / orientation / view tweaks survive. This is safe
        # because the .pse was built from the WT__WT CIF on the previous
        # full run (correct residue numbering, correct red ranges); we are
        # only refreshing the PNG. Pass --force to start over from CIF.
        print(f"[KEEP] {pair} {stoich} red: .pse exists; "
              f"re-rendering PNG from it (preserves color/orientation edits)")
        if raw.exists():
            raw.unlink()
        rc = run_pymol(str(RENDER_PSE), "--", "--pse", str(pse), "--out", str(raw))
        if rc != 0 or not raw.exists():
            print(f"[FAIL] {pair} {stoich} red: render_from_pse exit {rc}")
            return None
    else:
        cif = pair_cif("WT__WT", stoich)
        if not cif.exists():
            print(f"[SKIP] {pair} {stoich} red: missing WT__WT CIF {cif.name}")
            return None
        # Fresh render: use the WT__WT CIF (variant CIFs have AF3-renumbered
        # chains starting at 1, so red selections on WT residue ranges would
        # land on the wrong band). View baseline = canonical fig0_*.json.
        view_path = MONO_VIEW if stoich == "monomeric" else TRI_VIEW
        args = [
            str(RENDER_RED), "--",
            "--cif", str(cif),
            "--out", str(raw),
            "--target", "hap2",
            "--residues", hap2_del,
            "--also-dmp-residues", dmp_del,
            "--view-json", str(view_path),
        ]
        print(f"[RUN]  {pair} {stoich} red: HAP2 red={hap2_del}, DMP red={dmp_del} "
              f"(WT__WT CIF, view={view_path.name})")
        rc = run_pymol(*args)
        if rc != 0:
            print(f"[FAIL] {pair} {stoich} red: PyMOL exit {rc}")
            return None

    if not raw.exists():
        return None

    if cropped.exists():
        cropped.unlink()
    run_py(str(CROP), "--image", str(raw), "--out", str(cropped), "--pad", "12")
    if not cropped.exists():
        return None

    hap2_var, dmp_var = pair.split("__", 1)
    hap2_label = HAP2_REGION_LABEL.get(hap2_var, f"HAP2 {hap2_var}")
    dmp_label = DMP_REGION_LABEL.get(dmp_var, f"DMP {dmp_var}")
    legend = f"Red: {hap2_label} + {dmp_label}"

    iptm = read_pair_iptm(pair, stoich)
    overlay_args = [sys.executable, str(OVERLAY),
                    "--image", str(cropped), "--out", str(labeled),
                    "--label", legend]
    if iptm is not None:
        cat = classify_iptm(iptm, stoich)
        overlay_args += ["--iptm", f"{iptm:.4f}", "--iptm-category", cat]
    r = subprocess.run(overlay_args, cwd=PROJECT, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
        return None
    print(f"[OK]   {pair} {stoich} red: {labeled.name}"
          + (f"  (ipTM={iptm:.2f}, {cat})" if iptm is not None else ""))
    return labeled


def combine_panel(out_path: Path, pairs_pngs, captions, panel_h: float = 3.4) -> None:
    cmd = [sys.executable, str(COMBINE), "--out", str(out_path),
           "--cols", "2", "--panel-height", str(panel_h), "--dpi", "150"]
    for p in pairs_pngs:
        cmd.extend(["--image", str(p)])
    for c in captions:
        cmd.extend(["--caption", c])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
    print(f"  {out_path.name}: rc={r.returncode}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="re-render even if .pse exists (default: render PNG from existing .pse)")
    args = ap.parse_args()

    mod_mo = COMBINED_DIR / "model_monomeric"
    mod_pf = COMBINED_DIR / "model_postfusion"
    red_mo = COMBINED_DIR / "deletion_red_monomeric"
    red_pf = COMBINED_DIR / "deletion_red_postfusion"
    for d in (mod_mo, mod_pf, red_mo, red_pf):
        d.mkdir(parents=True, exist_ok=True)

    print(f"[plan] {len(SPECIAL_PAIRS)} special pairs x 2 stoich x 2 views "
          f"= {4*len(SPECIAL_PAIRS)} renders (model + red)")
    for pair, hap2_del, dmp_del, _descr in SPECIAL_PAIRS:
        render_pair(pair, "monomeric",       hap2_del, dmp_del, mod_mo, args.force)
        render_pair(pair, "postfusion_like", hap2_del, dmp_del, mod_pf, args.force)
        render_pair_red(pair, "monomeric",       hap2_del, dmp_del, red_mo, args.force)
        render_pair_red(pair, "postfusion_like", hap2_del, dmp_del, red_pf, args.force)

    # Grouped panels: one row per stoich containing all special pairs side by
    # side. The "combined" panels pair red + model views per variant; the
    # original "model-only" panels are retained for backward-compatibility.
    print("\n## Combined_Variants grouped panels ##")
    for stoich, mod_dir, red_dir, suffix in [
        ("monomeric",       mod_mo, red_mo, "monomeric"),
        ("postfusion_like", mod_pf, red_pf, "postfusion"),
    ]:
        # Model-only panel (legacy layout).
        pngs, caps = [], []
        for pair, _hd, _dd, _descr in SPECIAL_PAIRS:
            p = mod_dir / f"hap2_dmp_{pair}_model_{suffix}_structure_cropped.png"
            if p.exists():
                pngs.append(p)
                caps.append(f"{pair} (combined model)")
        if pngs:
            out_panel = COMBINED_DIR / f"combined_variants_{suffix}.png"
            combine_panel(out_panel, pngs, caps)

        # Red + model paired panel.
        pngs2, caps2 = [], []
        for pair, _hd, _dd, _descr in SPECIAL_PAIRS:
            red = red_dir / f"hap2_dmp_{pair}_red_{suffix}_labeled.png"
            model = mod_dir / f"hap2_dmp_{pair}_model_{suffix}_structure_cropped.png"
            if red.exists():
                pngs2.append(red); caps2.append(f"{pair} (red)")
            if model.exists():
                pngs2.append(model); caps2.append(f"{pair} (model)")
        if pngs2:
            out_panel2 = COMBINED_DIR / f"combined_variants_red_and_model_{suffix}.png"
            combine_panel(out_panel2, pngs2, caps2)

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
