#!/usr/bin/env python3
"""
Render the trimeric HAP2 + DMP from a Stage 14 WT__WT AF3 model, colour-
coded by the deletion-variant boundaries used in this thesis. The bands
are chosen to match the variants defined in 14_Interaction_Domain_Mapping
CONFIG.toml so the reader can map a heatmap row/column directly onto the
structural region removed by that variant.

All boundaries are in SmelHAP2 / SmelDMPv5_10.610 numbering (the AtHAP2/AtDMP8
ranges from Wang et al. 2022 already remapped onto the eggplant orthologs;
source: [hap2_variants.coords.sm] and [dmp_variants.coords.sm] in
14_Interaction_Domain_MappingCONFIG.toml).

HAP2 bands (SmelHAP2, 804 aa; class II fusion protein with C-anchored TM):
    1-21        Pre-ectodomain                              (not deleted)
    22-185      Ectodomain D1 proximal                      (delEcto removes 22-589)
    186-230     Ectodomain D2 N-flank                       (delEctoD2 removes 186-325)
    231-247     Fusion loop / cd-loop                       (delFL removes; nested in D2)
    248-325     Ectodomain D2 C-flank                       (delEctoD2 cont.)
    326-589     Ectodomain D1 distal / stem                 (delEcto cont.)
    590-619     Pre-TMD linker                              (not deleted)
    620-641     Transmembrane domain                        (delTMD removes; C-anchored)
    642-654     Post-TMD linker                             (not deleted)
    655-804     C-terminal cytoplasmic tail                 (delC removes)

DMP colouring (topology-based, NOT deletion-variant ranges): mirrors the
`[DMP]` palette in config/colors_config/protein_structure_colors.toml
referenced by 08_protein_structureCONFIG.toml.
    Amphipathic helices     golden amber (#D9A621)   helices in residues 1-83
    Transmembrane helices   deep orange-red (#D43005)  helices in residues 84-220
    Beta sheets             warm peach (#FFB366)     all sheets in the chain
    Loops                   gray70                    everything else

Usage:
    pymol -cq render_hap2_domain_map.py -- \
        --cif <path/to/fold_WT__WT_postfusion_like_model_0.cif> \
        --out <path/to/hap2_dmp_domain_map.png> \
        [--ray-width 1800] [--ray-height 1200]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# pymol is imported lazily inside main() so this module can be imported
# (for its band tables) from environments where pymol is not installed
# (e.g. the matplotlib-only compose_legend.py companion).


# HAP2 domain bands on SmelHAP2 (804 aa). (start, end, hex_colour, legend_label).
# The three short "linker / not deleted" segments (pre-ectodomain, pre-TMD,
# post-TMD) each get a DISTINCT muted colour so they can be told apart in
# the rendered structure -- previously all three shared one light grey and
# blended together.
HAP2_BANDS: list[tuple[int, int, str, str]] = [
    (1,   21,  "#F4A582", "Pre-ectodomain (not deleted; light salmon)"),
    (22,  185, "#228B22", "Ectodomain D1 proximal (delEcto: 22-589)"),
    (186, 230, "#4A7CC4", "Ectodomain D2 N-flank (delEctoD2: 186-325)"),
    (231, 247, "#FFD700", "Fusion loop (delFL: 231-247)"),
    (248, 325, "#4A7CC4", "Ectodomain D2 C-flank (delEctoD2 cont.)"),
    (326, 589, "#228B22", "Ectodomain D1 distal / stem (delEcto cont.)"),
    (590, 619, "#80CDC1", "Pre-TMD linker (not deleted; pale teal)"),
    (620, 641, "#7B3F99", "Transmembrane domain (delTMD: 620-641)"),
    (642, 654, "#E78AC3", "Post-TMD linker (not deleted; pale magenta)"),
    (655, 804, "#DC143C", "C-terminal cytoplasmic tail (delC: 655-804)"),
]

# DMP topology palette (mirrors [DMP] in protein_structure_colors.toml).
# (start, end, hex, label). start=None means "secondary-structure rule, no
# fixed residue range"; compose_legend.py renders those rows without the
# "start-end:" prefix.
DMP_BANDS: list[tuple[int | None, int | None, str, str]] = [
    (None, None, "#D9A621", "Amphipathic helices (helices in residues 1-83; delN zone)"),
    (None, None, "#D43005", "Transmembrane helices (helices in residues 84-220; delTMDcore zone)"),
    (None, None, "#FFB366", "Extracellular beta sheets"),
    (None, None, "#B3B3B3", "Loops (everything else)"),
]
# Selection ranges used at render time (must match the legend descriptions).
DMP_N_ZONE_END = 83          # residues 1..end inclusive are the N-terminal zone
DMP_TM_ZONE_START = 84       # residues start..end inclusive are the TM core zone
DMP_TM_ZONE_END = 220

HAP2_LENGTH_THRESHOLD = 400     # chains >= this length classed as HAP2
BACKGROUND = "black"


def main() -> int:
    from pymol import cmd  # local import: this module is also imported by
                           # compose_legend.py in a non-pymol env.

    def hex_to_pymol(name: str, hex_code: str) -> None:
        h = hex_code.lstrip("#")
        r, g, b = (int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
        cmd.set_color(name, [r, g, b])

    def chain_lengths(obj: str) -> dict[str, int]:
        lengths: dict[str, int] = {}
        for ch in cmd.get_chains(obj):
            n = cmd.count_atoms(f"{obj} and chain {ch} and name CA")
            lengths[ch] = n
        return lengths

    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ray-width", type=int, default=1800)
    ap.add_argument("--ray-height", type=int, default=1200)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--hap2-chains", default=None,
                    help='Comma-separated HAP2 chain IDs, e.g. "A,B,C". '
                         "Bypasses auto-classification by chain length (needed "
                         "for truncated variants where HAP2 < DMP in size).")
    ap.add_argument("--dmp-chains", default=None,
                    help='Comma-separated DMP chain IDs, e.g. "D".')
    ap.add_argument("--view-json", default=None,
                    help='Optional JSON file with {"rotation": [9 floats]} '
                         "from a saved PyMOL view (cmd.get_view()[:9]). "
                         "Overrides the auto SVD orientation so a hand-tuned "
                         "pose from a .pse can be reused across variants.")
    ap.add_argument("--deleted-residues", default=None,
                    help='Comma-separated SmelHAP2-numbering ranges removed '
                         "by the HAP2 variant, e.g. \"22-589\" or "
                         "\"22-589,655-804\". When set, HAP2_BANDS colours "
                         "are applied by mapping the variant\'s 1..N residues "
                         "back to their WT identities (variant residue i = "
                         "i-th surviving WT residue). Omit for WT HAP2.")
    ap.add_argument("--hap2-wt-length", type=int, default=804,
                    help="Length of full-length SmelHAP2 in residues "
                         "(default 804). Used with --deleted-residues for the "
                         "WT-aware residue mapping.")
    ap.add_argument("--dmp-deleted-residues", default=None,
                    help='Comma-separated SmelDMP-numbering ranges removed by '
                         "the DMP variant, e.g. \"1-83\" (delN) or "
                         "\"84-220\" (delTMDcore). When set, the DMP topology "
                         "palette is applied by mapping the variant DMP chain\'s "
                         "1..N residues back to their WT DMP identities. Omit "
                         "for WT DMP.")
    ap.add_argument("--dmp-wt-length", type=int, default=222,
                    help="Length of full-length SmelDMPv5_10.610 in residues "
                         "(default 222). Used with --dmp-deleted-residues for "
                         "the WT-aware residue mapping.")
    ap.add_argument("--dmp-novel-tail", type=int, default=0,
                    help="Number of NEW residues appended to the DMP variant "
                         "chain after the last surviving WT residue (frameshift "
                         "tail length). These trailing residues are coloured "
                         "with a distinct 'frameshift' colour so they are not "
                         "confused with WT-mapped DMP topology.")
    args = ap.parse_args()

    cif = Path(args.cif).resolve()
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    cmd.reinitialize()
    cmd.bg_color(BACKGROUND)
    cmd.load(str(cif), "complex")

    # Classify chains. Manual override wins; otherwise auto-detect by CA
    # count (works for full-length variants but fails for severe truncations
    # where HAP2 drops below DMP size — pass --hap2-chains / --dmp-chains).
    lens = chain_lengths("complex")
    if args.hap2_chains or args.dmp_chains:
        hap2_chains = sorted(c.strip() for c in (args.hap2_chains or "").split(",") if c.strip())
        dmp_chains  = sorted(c.strip() for c in (args.dmp_chains  or "").split(",") if c.strip())
    else:
        hap2_chains = sorted(c for c, n in lens.items() if n >= HAP2_LENGTH_THRESHOLD)
        dmp_chains  = sorted(c for c, n in lens.items() if n <  HAP2_LENGTH_THRESHOLD)
    if not hap2_chains:
        print(f"[ERROR] no HAP2 chains identified (auto-threshold {HAP2_LENGTH_THRESHOLD} aa). "
              "Pass --hap2-chains A,B,C for truncated variants.", file=sys.stderr)
        return 2
    print(f"[INFO] HAP2 chains: {hap2_chains}  (lengths {[lens.get(c, 0) for c in hap2_chains]})")
    print(f"[INFO] DMP  chains: {dmp_chains}  (lengths {[lens.get(c, 0) for c in dmp_chains]})")

    # Register HAP2 band colours (deletion-variant ranges).
    for i, (_, _, hex_code, _) in enumerate(HAP2_BANDS):
        hex_to_pymol(f"hap2_band{i}", hex_code)

    # Register DMP topology palette colours, indexed by their role in DMP_BANDS:
    #   0 = amphipathic, 1 = TM, 2 = sheet, 3 = loops
    for i, (_, _, hex_code, _) in enumerate(DMP_BANDS):
        hex_to_pymol(f"dmp_band{i}", hex_code)

    # Colour HAP2 bands on every HAP2 chain.
    #
    # AF3 renumbers each chain's residues 1..N regardless of which WT
    # residues are present, so a deletion variant's residue 22 is NOT WT
    # residue 22 -- it is the 22nd surviving residue. Without remapping,
    # HAP2_BANDS (which is in WT numbering) would shift all colours toward
    # the N-terminus. When --deleted-residues is set we therefore build a
    # variant->WT residue map and group variant residues by which WT band
    # they belong to before painting.
    hap2_sel = "complex and chain " + "+".join(hap2_chains)
    if args.deleted_residues:
        deleted = set()
        for chunk in args.deleted_residues.split(","):
            a, b = chunk.strip().split("-")
            deleted.update(range(int(a), int(b) + 1))
        # variant_to_wt[i] = WT residue number for variant residue i+1
        variant_to_wt = [wt for wt in range(1, args.hap2_wt_length + 1)
                         if wt not in deleted]
        # Group variant residues by their assigned band index. Use a
        # default band of -1 (skip) for any residue with no band.
        def band_of_wt(wt: int) -> int:
            for idx, (s, e, _, _) in enumerate(HAP2_BANDS):
                if s <= wt <= e:
                    return idx
            return -1
        per_band: dict[int, list[int]] = {}
        for vi, wt in enumerate(variant_to_wt, start=1):
            b = band_of_wt(wt)
            if b >= 0:
                per_band.setdefault(b, []).append(vi)
        # Apply each band's colour to its variant residues. Compress
        # consecutive runs into PyMOL "a-b" tokens to keep selections short.
        def compress(nums: list[int]) -> str:
            if not nums:
                return ""
            nums = sorted(nums)
            parts = []
            start = prev = nums[0]
            for n in nums[1:]:
                if n == prev + 1:
                    prev = n
                else:
                    parts.append(f"{start}-{prev}" if start != prev else f"{start}")
                    start = prev = n
            parts.append(f"{start}-{prev}" if start != prev else f"{start}")
            return "+".join(parts)
        for i in range(len(HAP2_BANDS)):
            sel = compress(per_band.get(i, []))
            if sel:
                cmd.color(f"hap2_band{i}", f"{hap2_sel} and resi {sel}")
        print(f"[INFO] WT-aware coloring applied; deleted ranges {args.deleted_residues} "
              f"-> {len(variant_to_wt)} surviving residues mapped to bands")
    else:
        for i, (start, end, _, _) in enumerate(HAP2_BANDS):
            cmd.color(f"hap2_band{i}", f"{hap2_sel} and resi {start}-{end}")

    # Colour DMP by secondary-structure topology (matches stage 08's DMP
    # palette). PyMOL auto-assigns ss on load: 'H' = helix, 'S' = sheet,
    # 'L' (or unassigned) = loop.
    if dmp_chains:
        dmp_sel = "complex and chain " + "+".join(dmp_chains)
        # Register the frameshift-tail colour (vivid pink) once per call so
        # it shows up distinct from any band/topology hue.
        hex_to_pymol("dmp_frameshift_tail", "#FF1493")

        # Default every DMP residue to the loops colour first, then overlay
        # helices and sheets so any residue PyMOL did not flag as H or S
        # falls through to grey.
        cmd.color("dmp_band3", dmp_sel)

        if args.dmp_deleted_residues or args.dmp_novel_tail:
            # WT-aware DMP colouring: build a variant->WT residue map for
            # the DMP chain so zone tests (1-83 amphipathic, 84-220 TM)
            # apply against the WT residue identity, not the variant\'s
            # renumbered position. Frameshifted tail residues (no WT
            # mapping) are coloured with a distinct frameshift hue.
            deleted_dmp: set[int] = set()
            if args.dmp_deleted_residues:
                for chunk in args.dmp_deleted_residues.split(","):
                    a, b = chunk.strip().split("-")
                    deleted_dmp.update(range(int(a), int(b) + 1))
            wt_kept = [wt for wt in range(1, args.dmp_wt_length + 1)
                       if wt not in deleted_dmp]
            tail = max(0, int(args.dmp_novel_tail))
            variant_total = len(wt_kept) + tail
            # variant residue i (1-indexed) -> WT residue, or None for tail
            def variant_to_wt(i: int) -> int | None:
                if i <= len(wt_kept):
                    return wt_kept[i - 1]
                return None  # frameshifted-tail residue
            amphi_var: list[int] = []
            tm_var:    list[int] = []
            for vi in range(1, variant_total + 1):
                wt = variant_to_wt(vi)
                if wt is None:
                    continue
                if wt <= DMP_N_ZONE_END:
                    amphi_var.append(vi)
                elif DMP_TM_ZONE_START <= wt <= DMP_TM_ZONE_END:
                    tm_var.append(vi)
            def compress(nums: list[int]) -> str:
                if not nums: return ""
                nums = sorted(nums); parts = []
                s = p = nums[0]
                for n in nums[1:]:
                    if n == p + 1: p = n
                    else: parts.append(f"{s}-{p}" if s != p else f"{s}"); s = p = n
                parts.append(f"{s}-{p}" if s != p else f"{s}")
                return "+".join(parts)
            if amphi_var:
                cmd.color("dmp_band0", f"{dmp_sel} and ss H and resi {compress(amphi_var)}")
            if tm_var:
                cmd.color("dmp_band1", f"{dmp_sel} and ss H and resi {compress(tm_var)}")
            cmd.color("dmp_band2", f"{dmp_sel} and ss S")
            if tail > 0:
                tail_start = len(wt_kept) + 1
                tail_end = variant_total
                cmd.color("dmp_frameshift_tail",
                          f"{dmp_sel} and resi {tail_start}-{tail_end}")
            print(f"[INFO] DMP WT-aware coloring: deleted={args.dmp_deleted_residues or 'none'} "
                  f"tail={tail} amphipathic={len(amphi_var)} TM={len(tm_var)}")
        else:
            cmd.color("dmp_band0", f"{dmp_sel} and ss H and resi 1-{DMP_N_ZONE_END}")
            cmd.color("dmp_band1", f"{dmp_sel} and ss H and resi {DMP_TM_ZONE_START}-{DMP_TM_ZONE_END}")
            cmd.color("dmp_band2", f"{dmp_sel} and ss S")

    # Upright orientation, unified across monomeric and trimer cases.
    # HAP2 is a class II fusion protein anchored at its C-terminus (the
    # TMD spans residues 620-641 in SmelHAP2). We always want the membrane
    # normal -- HAP2's long principal axis -- vertical on the screen with
    # the TMD at the bottom, and (where DMP is positioned off-axis) the
    # HAP2->DMP direction pointing toward the viewer so the contact patch
    # faces forward. Implemented via a principal-axis (SVD) computation
    # because cmd.orient lands the model in an arbitrary tilt and rotate()
    # about a screen axis cannot recover a deterministic upright pose.
    tmd_sel = f"{hap2_sel} and resi 620-641"
    import numpy as np
    hap2_ca = cmd.get_coords(f"{hap2_sel} and name CA")
    tmd_ca  = cmd.get_coords(f"{tmd_sel} and name CA")
    dmp_ca  = cmd.get_coords(f"{dmp_sel} and name CA") if dmp_chains else None
    if hap2_ca is not None and len(hap2_ca) >= 20:
        hap2_ca = np.asarray(hap2_ca)
        hap2_com = hap2_ca.mean(axis=0)
        # Principal long axis of HAP2 via SVD on the centred CA cloud.
        _, _, vh = np.linalg.svd(hap2_ca - hap2_com, full_matrices=False)
        long_axis = vh[0] / np.linalg.norm(vh[0])
        # Orient long_axis so it points FROM the TMD TOWARD the rest of
        # HAP2 (i.e. "up" once installed as the camera +Y axis). For severe
        # truncations where the TMD is the only ectodomain element left,
        # the dot-product flip still produces a deterministic pose because
        # the non-TMD residues anchor the comparison.
        if tmd_ca is not None and len(tmd_ca) > 0:
            tmd_com = np.asarray(tmd_ca).mean(axis=0)
            if np.dot(long_axis, tmd_com - hap2_com) > 0:
                long_axis = -long_axis
        up = long_axis
        # Forward = HAP2->DMP direction, projected orthogonal to up. For
        # the trimer + DMP postfusion-like case DMP sits near the membrane
        # along the threefold (almost collinear with up), so the
        # perpendicular component is small; we fall back to a deterministic
        # perpendicular direction in that case.
        if dmp_ca is not None and len(dmp_ca) > 0:
            dmp_com = np.asarray(dmp_ca).mean(axis=0)
            to_dmp = dmp_com - hap2_com
            forward = to_dmp - np.dot(to_dmp, up) * up
            n_fwd = np.linalg.norm(forward)
        else:
            forward = np.zeros(3)
            n_fwd = 0.0
        if n_fwd < 1e-3:
            tmp = np.array([1.0, 0.0, 0.0]) if abs(up[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            forward = tmp - np.dot(tmp, up) * up
            forward /= np.linalg.norm(forward)
        else:
            forward /= n_fwd
        # Right-handed camera frame: +X right, +Y up, +Z toward viewer.
        right = np.cross(up, forward)
        right /= np.linalg.norm(right)
        forward = np.cross(right, up)
        # PyMOL get_view stores the rotation row-major; rows are the
        # camera axes expressed in model coordinates.
        R = np.array([right, up, forward])
        new_view = list(R.flatten()) + list(cmd.get_view())[9:]
        cmd.set_view(new_view)
        stoich = "monomeric" if len(hap2_chains) == 1 else f"{len(hap2_chains)}-chain"
        print(f"[INFO] {stoich}: HAP2 long axis vertical (TMD bottom)"
              + (", HAP2->DMP facing viewer" if n_fwd >= 1e-3 else ", arbitrary azimuth (DMP collinear or absent)"))
    elif cmd.count_atoms(tmd_sel) >= 30:
        cmd.orient(tmd_sel)
        cmd.rotate("x", 90)
        print("[WARN] not enough HAP2 CA atoms for SVD orient; fell back to TMD-based trimer orient")
    else:
        cmd.orient(hap2_sel)
        print("[WARN] not enough HAP2 CA atoms for SVD orient; fell back to cmd.orient(hap2)")

    # Hand-tuned rotation override: replace the rotation portion of the
    # current view matrix with the one extracted from a saved .pse.
    if args.view_json:
        import json
        with open(args.view_json) as fh:
            saved = json.load(fh)
        rot = list(saved["rotation"])
        assert len(rot) == 9, f"--view-json rotation must have 9 floats, got {len(rot)}"
        new_view = rot + list(cmd.get_view())[9:]
        cmd.set_view(new_view)
        print(f"[INFO] applied saved rotation from {args.view_json}")

    cmd.zoom("complex", buffer=4, complete=1)
    cmd.set("ray_opaque_background", 1)

    # Save the reloadable PyMOL session (colours, view, selections) so the
    # exact figure state can be reopened in the GUI for inspection or
    # touch-ups: File > Open... > <session>.pse.
    session_pse = out.with_suffix(".pse")
    cmd.save(str(session_pse))
    print(f"[OK] wrote {session_pse}")

    # Render structure to a temporary file first; the legend composite is
    # done in a second step by the base-Python wrapper since matplotlib is
    # typically not installed in the PyMOL conda env.
    structure_png = out.with_name(out.stem + "_structure.png")
    cmd.ray(args.ray_width, args.ray_height)
    cmd.png(str(structure_png), dpi=args.dpi, ray=0)
    print(f"[OK] wrote {structure_png}")

    # If matplotlib is available in this env, composite directly. Otherwise
    # the caller is expected to run compose_legend.py against structure_png.
    try:
        compose_with_legend(structure_png, out, args.dpi)
        print(f"[OK] wrote {out}")
    except ImportError:
        print("[INFO] matplotlib not in this env; skipping legend composite. "
              "Run compose_legend.py separately to build the final figure.")
    return 0


def compose_with_legend(structure_png: Path, out: Path, dpi: int) -> None:
    """Paste the PyMOL structure render into a matplotlib figure and add a
    coloured-band legend panel listing HAP2 and DMP variant boundaries."""
    import matplotlib.image as mpimg
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    img = mpimg.imread(str(structure_png))
    h_in = 8
    w_in = h_in * (img.shape[1] / img.shape[0]) + 5.5  # extra width for legend

    fig = plt.figure(figsize=(w_in, h_in), dpi=dpi)
    gs = GridSpec(1, 2, width_ratios=[img.shape[1] / img.shape[0] * h_in, 5.5], wspace=0.02)

    ax_img = fig.add_subplot(gs[0, 0])
    ax_img.imshow(img)
    ax_img.set_axis_off()

    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.set_axis_off()

    def patches(bands, title, y_start):
        ax_leg.text(0.0, y_start, title, fontsize=11, fontweight="bold",
                    transform=ax_leg.transAxes, va="top")
        y = y_start - 0.04
        for start, end, hex_code, label in bands:
            ax_leg.add_patch(mpatches.Rectangle(
                (0.0, y - 0.025), 0.06, 0.025,
                transform=ax_leg.transAxes, facecolor=hex_code,
                edgecolor="black", linewidth=0.4,
            ))
            ax_leg.text(0.08, y - 0.012, f"{start}-{end}: {label}",
                        fontsize=8.5, transform=ax_leg.transAxes, va="center")
            y -= 0.035
        return y - 0.02

    y_after_hap2 = patches(HAP2_BANDS, "SmelHAP2 trimer (chains A/B/C)", 0.98)
    patches(DMP_BANDS, "SmelDMPv5_10.610 (chain D)", y_after_hap2)

    fig.savefig(str(out), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__" or __name__ == "pymol":
    sys.exit(main())
