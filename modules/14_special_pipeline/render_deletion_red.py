#!/usr/bin/env python3
"""
Render a HAP2 + DMP AF3 complex with one or more deletion-variant residue
ranges painted bright red on every HAP2 chain, with the rest of the
structure rendered in a dim grey. Used to visualise the regions removed by
the catastrophic-loss variants delEcto (22-589) and delEctoAndC
(22-589, 655-804) on the SmelHAP2 postfusion-like 3:1 complex.

Residue ranges are given in SmelHAP2 numbering (the remapped values from
14_Interaction_Domain_MappingCONFIG.toml [hap2_variants.coords.sm]). The
AtHAP2 ranges that the manuscript tables quote (66-425, 596-705 etc.) are
remapped equivalents - do NOT pass AtHAP2 numbering directly here.

Usage:
    pymol -cq render_deletion_red.py -- \
        --cif <path/to/fold_WT__WT_postfusion_like_model_0.cif> \
        --residues "22-589" \
        --out <path/to/output.png> \
        [--label "delEcto (SmelHAP2 22-589)"]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


HAP2_LENGTH_THRESHOLD = 400
BACKGROUND = "black"


def parse_ranges(spec: str) -> list[tuple[int, int]]:
    """Parse "22-589,655-804" into [(22, 589), (655, 804)]."""
    out: list[tuple[int, int]] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        a, b = chunk.split("-")
        out.append((int(a), int(b)))
    return out


def main() -> int:
    from pymol import cmd

    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True)
    ap.add_argument("--residues", required=True,
                    help='SmelHAP2 ranges, e.g. "22-589" or "22-589,655-804"')
    ap.add_argument("--out", required=True)
    ap.add_argument("--ray-width", type=int, default=1800)
    ap.add_argument("--ray-height", type=int, default=1200)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--label", default=None,
                    help="Optional label burned into PyMOL view (legend uses --label).")
    ap.add_argument("--view-json", default=None,
                    help='Optional JSON file with {"rotation": [9 floats]} '
                         "from a saved PyMOL view (cmd.get_view()[:9]). "
                         "Overrides the auto SVD orientation so a hand-tuned "
                         "pose from a .pse can be reused across variants.")
    args = ap.parse_args()

    cif = Path(args.cif).resolve()
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    ranges = parse_ranges(args.residues)
    if not ranges:
        print("[ERROR] no residue ranges parsed from --residues", file=sys.stderr)
        return 2

    cmd.reinitialize()
    cmd.bg_color(BACKGROUND)
    cmd.load(str(cif), "complex")

    # Classify chains.
    lens = {ch: cmd.count_atoms(f"complex and chain {ch} and name CA")
            for ch in cmd.get_chains("complex")}
    hap2_chains = sorted(c for c, n in lens.items() if n >= HAP2_LENGTH_THRESHOLD)
    dmp_chains  = sorted(c for c, n in lens.items() if n <  HAP2_LENGTH_THRESHOLD)
    print(f"[INFO] HAP2 chains: {hap2_chains}  (lengths {[lens[c] for c in hap2_chains]})")
    print(f"[INFO] DMP  chains: {dmp_chains}  (lengths {[lens[c] for c in dmp_chains]})")
    if not hap2_chains:
        print(f"[ERROR] no chains >= {HAP2_LENGTH_THRESHOLD} aa found", file=sys.stderr)
        return 2

    # Register colours.
    def hex_to_pymol(name: str, hex_code: str) -> None:
        h = hex_code.lstrip("#")
        r, g, b = (int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
        cmd.set_color(name, [r, g, b])
    hex_to_pymol("dim_grey",    "#5A5A5A")
    hex_to_pymol("dmp_grey",    "#3A3A3A")
    hex_to_pymol("deletion_red", "#E60000")

    # Paint baseline.
    cmd.color("dim_grey", "complex")
    if dmp_chains:
        dmp_sel = "complex and chain " + "+".join(dmp_chains)
        cmd.color("dmp_grey", dmp_sel)

    # Paint deletion zones red on every HAP2 chain.
    hap2_sel = "complex and chain " + "+".join(hap2_chains)
    range_sel_parts = [f"resi {a}-{b}" for a, b in ranges]
    red_sel = f"({hap2_sel}) and ({' or '.join(range_sel_parts)})"
    n_red = cmd.count_atoms(red_sel)
    print(f"[INFO] painting {n_red} atoms red across ranges {ranges} on chains {hap2_chains}")
    cmd.color("deletion_red", red_sel)

    # Upright orientation via SVD on HAP2 CA cloud (same logic as
    # render_hap2_domain_map.py): HAP2 long axis -> vertical (TMD at
    # bottom); HAP2 -> DMP direction projected perpendicular to the long
    # axis -> facing the viewer. Works for both monomeric (1 HAP2 chain)
    # and trimer (3 HAP2 chains).
    tmd_sel = f"{hap2_sel} and resi 620-641"
    import numpy as np
    hap2_ca = cmd.get_coords(f"{hap2_sel} and name CA")
    tmd_ca  = cmd.get_coords(f"{tmd_sel} and name CA")
    dmp_ca  = cmd.get_coords(f"{dmp_sel} and name CA") if dmp_chains else None
    if hap2_ca is not None and len(hap2_ca) >= 20:
        hap2_ca = np.asarray(hap2_ca)
        hap2_com = hap2_ca.mean(axis=0)
        _, _, vh = np.linalg.svd(hap2_ca - hap2_com, full_matrices=False)
        long_axis = vh[0] / np.linalg.norm(vh[0])
        if tmd_ca is not None and len(tmd_ca) > 0:
            tmd_com = np.asarray(tmd_ca).mean(axis=0)
            if np.dot(long_axis, tmd_com - hap2_com) > 0:
                long_axis = -long_axis
        up = long_axis
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
        right = np.cross(up, forward)
        right /= np.linalg.norm(right)
        forward = np.cross(right, up)
        R = np.array([right, up, forward])
        cmd.set_view(list(R.flatten()) + list(cmd.get_view())[9:])
        stoich = "monomeric" if len(hap2_chains) == 1 else f"{len(hap2_chains)}-chain"
        print(f"[INFO] {stoich}: HAP2 long axis vertical (TMD bottom)"
              + (", HAP2->DMP facing viewer" if n_fwd >= 1e-3 else ", arbitrary azimuth (DMP collinear/absent)"))
    else:
        cmd.orient(hap2_sel)
        print("[WARN] SVD orient skipped (insufficient CA); fell back to cmd.orient(hap2)")

    # Hand-tuned rotation override: replace the rotation portion of the
    # current view matrix with one extracted from a saved .pse.
    if args.view_json:
        import json
        with open(args.view_json) as fh:
            saved = json.load(fh)
        rot = list(saved["rotation"])
        assert len(rot) == 9, f"--view-json rotation must have 9 floats, got {len(rot)}"
        cmd.set_view(rot + list(cmd.get_view())[9:])
        print(f"[INFO] applied saved rotation from {args.view_json}")

    cmd.zoom("complex", buffer=4, complete=1)
    cmd.set("ray_opaque_background", 1)

    # Save reloadable session (so user can re-orient in GUI).
    session_pse = out.with_suffix(".pse")
    cmd.save(str(session_pse))
    print(f"[OK] wrote {session_pse}")

    # Ray-traced structure render.
    cmd.ray(args.ray_width, args.ray_height)
    cmd.png(str(out), dpi=args.dpi, ray=0)
    print(f"[OK] wrote {out}")
    return 0


if __name__ == "__main__" or __name__ == "pymol":
    sys.exit(main())
