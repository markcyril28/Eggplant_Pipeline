#!/usr/bin/env python3
"""Op 8: Comparative overlay render of AlphaFold3 vs SWISS-MODEL structures.

Superimposes each matched AF3/SWISS gene pair and renders an overlay image
with the two structures coloured distinctly.  All visual parameters are
controlled via the [comparison] section of protein_structure_colors.toml and
the shared [rendering] section.

Orientation matches Op 5 (render_structures_dmp.py):
  1. reference.pse in AlphaFold3_Results/ → extract rotation matrix.
  2. Superpose each AF3 structure onto reference.cif, apply that rotation.
  3. PCA fallback (helix axis → Y-up, sheet centre → Z-forward) when
     reference files are absent.
  4. Pass 1 computes a uniform camera Z across all AF3 structures so every
     overlay shares the same zoom level.  AF3 is oriented first; SWISS is
     then superposed onto AF3's final (post-orientation) coordinates.

Outputs:
  Comparison_Results/{gene}/overlay_{bg}.jpg

Usage:
    python3 comparative_render.py --run-dir /path/to/genome_dir \
        --color-config /path/to/protein_structure_colors.toml
"""

import argparse
import os
import re
from pathlib import Path

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

try:
    from PIL import Image
except ImportError:
    Image = None

# Matches raw AlphaFold3 timestamped extract folders: YYYY_MM_DD_HH_MM[_SS]
_AF3_TIMESTAMP_RE = re.compile(r"^\d{4}_\d{2}_\d{2}_\d{2}_\d{2}")


def af3_to_swiss_name(af3_name: str) -> str:
    idx = af3_name.rfind("_")
    if idx == -1:
        return af3_name
    base = af3_name[:idx]
    version = af3_name[idx + 1:]
    first_under = base.find("_")
    if first_under == -1:
        return base.upper() + "." + version
    prefix = base[:first_under].upper()
    return prefix + base[first_under:] + "." + version


# ── Rendering constants (defaults; overridden from [rendering] in colors config) ─
IMAGE_WIDTH  = 1200
IMAGE_HEIGHT = 1200
DPI          = 900
JPEG_QUALITY = 95
ANTIALIAS    = 2
RAY_SHADOWS  = 0
SPECULAR     = 0.3
DEPTH_CUE    = 1
_DEFAULT_BACKGROUNDS = ("black", "white")

# ── Overlay constants (defaults; overridden from [comparison] in colors config) ─
_AF3_COLOR    = "salmon"   # PyMOL colour for AF3 structure
_SWISS_COLOR  = "palecyan" # PyMOL colour for SWISS-MODEL structure
_SWISS_TRANSP = 0.20       # Cartoon transparency for SWISS chain (0=opaque, 1=invisible)
_ZOOM_BUFFER  = 10         # Buffer pixels around structures during Pass 1 zoom

# ── Orientation constants (defaults; overridden from [rendering] in colors config) ─
_ORIENTATION_METHOD = "reference"   # "reference" or "pca"
_REF_PSE_FILE       = "reference.pse"  # filename within AlphaFold3_Results/
_REF_CIF_FILE       = "reference.cif"  # filename within AlphaFold3_Results/
_HELIX_UP           = False         # PCA fallback: False = helix axis → −Y (down; ectodomain up)
_SLAB_DEPTH         = 200           # clipping half-depth from camera centre (Å)


def _init_pymol():
    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd
    return cmd


def _apply_ray_settings(cmd):
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias",  ANTIALIAS)
    cmd.set("ray_shadows", RAY_SHADOWS)
    cmd.set("depth_cue",  DEPTH_CUE)
    cmd.set("specular",   SPECULAR)
    cmd.set("max_threads", 0)
    cmd.set("hash_max",   250)


def _ray_and_save(cmd, out_path: Path, bg: str):
    """Ray-trace and save as JPEG (via PNG intermediate if Pillow is available)."""
    cmd.bg_color(bg)
    _apply_ray_settings(cmd)
    cmd.ray(IMAGE_WIDTH, IMAGE_HEIGHT)

    if Image is not None:
        png_path = out_path.with_suffix(".png")
        cmd.png(str(png_path), width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=DPI, quiet=1)
        img = Image.open(str(png_path)).convert("RGB")
        img.save(str(out_path), "JPEG", quality=JPEG_QUALITY, dpi=(DPI, DPI))
        png_path.unlink()
    else:
        cmd.png(str(out_path.with_suffix(".png")),
                width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=DPI, quiet=1)


# ── Orientation helpers (mirrored from render_structures_dmp.py) ─────────────

def _get_ca_coords(cmd, selection: str):
    """Return Nx3 NumPy array of Cα coordinates for *selection*, or None."""
    from pymol import stored
    stored.xyz = []
    cmd.iterate_state(1, f"({selection}) and name CA",
                      "stored.xyz.append((x, y, z))")
    if not stored.xyz:
        return None
    return np.array(stored.xyz)


def _principal_axis(coords: np.ndarray) -> np.ndarray:
    centered = coords - coords.mean(axis=0)
    _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    return Vt[0]


def _zoom_fit(cmd, gene: str, cam_z: float | None, slab: tuple | None) -> None:
    """Centre view on Cα centre of mass; apply uniform camera Z when provided."""
    cmd.zoom(gene, buffer=0)
    all_ca = _get_ca_coords(cmd, gene)
    if all_ca is not None and len(all_ca) >= 1:
        center = all_ca.mean(axis=0)
        view = list(cmd.get_view())
        view[12] = float(center[0])
        view[13] = float(center[1])
        view[14] = float(center[2])
        cmd.set_view(view)
    if cam_z is not None and slab is not None:
        view = list(cmd.get_view())
        view[11] = cam_z
        view[15] = slab[0]
        view[16] = slab[1]
        cmd.set_view(view)


def _orient_pca(cmd, gene: str, cam_z: float | None, slab: tuple | None) -> None:
    """PCA-based orientation: helix axis → Y-up, sheet centre → Z-forward."""
    helix_ca = _get_ca_coords(cmd, f"{gene} and ss h")
    if helix_ca is None or len(helix_ca) < 4:
        cmd.orient(gene)
        _zoom_fit(cmd, gene, cam_z, slab)
        return

    all_ca = _get_ca_coords(cmd, gene)
    center = all_ca.mean(axis=0) if all_ca is not None else np.zeros(3)

    Y = _principal_axis(helix_ca)
    Y = Y / np.linalg.norm(Y)

    # Flip Y so it points toward sheets (ectodomain up on screen)
    sheet_ca = _get_ca_coords(cmd, f"{gene} and ss s")
    if sheet_ca is not None and len(sheet_ca) >= 3:
        sheet_dir = sheet_ca.mean(axis=0) - center
        if np.dot(Y, sheet_dir) < 0:
            Y = -Y
    else:
        if (Y[1] < 0) == _HELIX_UP:
            Y = -Y
        centered_all = all_ca - center
        _, _, Vt = np.linalg.svd(centered_all, full_matrices=False)
        sheet_dir = Vt[1]

    sheet_dir = sheet_dir - np.dot(sheet_dir, Y) * Y
    norm = np.linalg.norm(sheet_dir)
    if norm < 1e-6:
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(perp, Y)) > 0.9:
            perp = np.array([0.0, 0.0, 1.0])
        sheet_dir = perp - np.dot(perp, Y) * Y
        norm = np.linalg.norm(sheet_dir)
    Z = sheet_dir / norm
    X = np.cross(Y, Z)
    X = X / np.linalg.norm(X)
    Z = np.cross(X, Y)

    R = np.column_stack([X, Y, Z])
    cmd.orient(gene)
    view = list(cmd.get_view())
    view[0:9] = R.flatten().tolist()
    cmd.set_view(view)
    _zoom_fit(cmd, gene, cam_z, slab)


def _apply_orientation(cmd, gene: str, ref_cif: Path | None,
                       ref_rotation: list | None,
                       cam_z: float | None, slab: tuple | None) -> None:
    """Superpose *gene* onto reference.cif and apply ref_rotation; fall back to PCA."""
    cmd.dss(gene)  # ensure secondary structure is assigned for PCA fallback
    if ref_rotation is not None and ref_cif is not None and ref_cif.exists():
        try:
            cmd.load(str(ref_cif), "_ref_struct")
            cmd.super(gene, "_ref_struct")
            cmd.delete("_ref_struct")
            cmd.orient(gene)
            view = list(cmd.get_view())
            view[0:9] = list(ref_rotation)
            cmd.set_view(view)
            _zoom_fit(cmd, gene, cam_z, slab)
            return
        except Exception as e:
            print(f"    Warning: superposition onto reference failed ({e}), using PCA")
            try:
                cmd.delete("_ref_struct")
            except Exception:
                pass
    _orient_pca(cmd, gene, cam_z, slab)


# ── Render ────────────────────────────────────────────────────────────────────

def render_overlay(cmd, af3_pdb: Path, swiss_pdb: Path, out_dir: Path,
                   backgrounds: tuple[str, ...], orient_fn) -> None:
    """Orient AF3 first, then superpose SWISS onto the oriented AF3, render per background.

    Order matters: reference orientation calls cmd.super("af3", ref_cif) which moves
    AF3's atom coordinates in world space.  SWISS must be aligned *after* that move,
    otherwise it ends up in AF3's pre-orientation coordinate frame and appears misaligned.
    """
    for bg in (backgrounds or _DEFAULT_BACKGROUNDS):
        cmd.reinitialize()
        cmd.load(str(af3_pdb), "af3")
        cmd.hide("everything")
        cmd.show("cartoon", "af3")

        # Orient AF3 first — this may transform AF3's world coordinates
        # (reference method superposes AF3 onto reference.cif).
        orient_fn("af3")

        # Save view before loading SWISS: PyMOL's auto_zoom fires on cmd.load()
        # and overwrites the camera position set by orient_fn.
        _saved_view = cmd.get_view()

        # Superpose SWISS onto AF3's final (post-orientation) coordinates.
        cmd.load(str(swiss_pdb), "swiss")
        cmd.super("swiss", "af3")

        # Restore the oriented view (SWISS is colocated with AF3 after super).
        cmd.set_view(_saved_view)

        cmd.show("cartoon", "swiss")
        cmd.color(_AF3_COLOR,   "af3")
        cmd.color(_SWISS_COLOR, "swiss")
        cmd.set("cartoon_transparency", _SWISS_TRANSP, "swiss")

        _ray_and_save(cmd, out_dir / f"overlay_{bg}.jpg", bg)


def render_confidence(cmd, af3_pdb: Path, out_dir: Path,
                      conf_colors: dict, backgrounds: tuple[str, ...],
                      orient_fn) -> None:
    """Render AF3 coloured by pLDDT confidence (B-factor column) per background."""
    thresholds = conf_colors.get("thresholds", [50, 70, 90])
    for bg in (backgrounds or _DEFAULT_BACKGROUNDS):
        cmd.reinitialize()
        cmd.load(str(af3_pdb), "af3")
        cmd.dss("af3")
        cmd.hide("everything")
        cmd.show("cartoon", "af3")

        cmd.set_color("plddt_vhigh", conf_colors.get("very_high", [0.00, 0.33, 0.84]))
        cmd.set_color("plddt_high",  conf_colors.get("high",      [0.40, 0.80, 0.95]))
        cmd.set_color("plddt_low",   conf_colors.get("low",       [1.00, 0.86, 0.07]))
        cmd.set_color("plddt_vlow",  conf_colors.get("very_low",  [1.00, 0.49, 0.27]))
        cmd.color("plddt_vlow",  "af3")
        cmd.color("plddt_low",   f"(af3) and b > {thresholds[0]}")
        cmd.color("plddt_high",  f"(af3) and b > {thresholds[1]}")
        cmd.color("plddt_vhigh", f"(af3) and b > {thresholds[2]}")

        orient_fn("af3")
        _ray_and_save(cmd, out_dir / f"confidence_{bg}.jpg", bg)


def render_deviation(cmd, af3_pdb: Path, swiss_pdb: Path, out_dir: Path,
                     backgrounds: tuple[str, ...], orient_fn) -> None:
    """Render AF3 coloured by per-residue Cα distance vs SWISS-MODEL.

    After superposition, the Cα–Cα distance for each matched residue is
    written into the B-factor column of af3, then coloured with a
    blue→white→red spectrum (low deviation = blue, high = red).
    """
    import math
    from pymol import stored

    for bg in (backgrounds or _DEFAULT_BACKGROUNDS):
        cmd.reinitialize()
        cmd.load(str(af3_pdb),   "af3")
        cmd.load(str(swiss_pdb), "swiss")
        cmd.super("swiss", "af3")
        cmd.dss("af3")

        # Gather Cα coordinates keyed by (chain, resi)
        stored.af3_ca = {}
        cmd.iterate_state(1, "af3 and name CA",
                          "stored.af3_ca[(chain, resi)] = (x, y, z)")
        stored.swiss_ca = {}
        cmd.iterate_state(1, "swiss and name CA",
                          "stored.swiss_ca[(chain, resi)] = (x, y, z)")

        # Per-residue Cα–Cα distance → stored.dev_map for use in cmd.alter
        stored.dev_map = {}
        for key, (ax, ay, az) in stored.af3_ca.items():
            if key in stored.swiss_ca:
                sx, sy, sz = stored.swiss_ca[key]
                stored.dev_map[key] = math.sqrt(
                    (ax - sx) ** 2 + (ay - sy) ** 2 + (az - sz) ** 2
                )

        if not stored.dev_map:
            print(f"    Warning: no matching residues for deviation render of "
                  f"{af3_pdb.parent.name} — skipping {bg}")
            continue

        # Write deviations into B-factor column of af3 (all atoms in each residue)
        cmd.alter("af3", "b = stored.dev_map.get((chain, resi), 0.0)")

        cmd.hide("everything")
        cmd.show("cartoon", "af3")
        cmd.spectrum("b", "blue_white_red", "af3")

        orient_fn("af3")
        _ray_and_save(cmd, out_dir / f"deviation_{bg}.jpg", bg)


# ── Main ──────────────────────────────────────────────────────────────────────

def render_all(run_dir: Path, overwrite: bool = False,
               backgrounds: tuple[str, ...] = (),
               render_types: tuple[str, ...] = ("overlay",),
               conf_colors: dict | None = None) -> int:
    af3_dir   = run_dir / "AlphaFold3_Results"
    swiss_dir = run_dir / "SWISS_Results"
    comp_dir  = run_dir / "Comparison_Results"

    if not af3_dir.is_dir():
        print(f"  Missing AlphaFold3_Results/ in {run_dir.name}")
        return 0
    needs_swiss = "overlay" in render_types or "deviation" in render_types
    if needs_swiss and not swiss_dir.is_dir():
        print(f"  Missing SWISS_Results/ in {run_dir.name} "
              "(required for overlay/deviation)")
        return 0

    cmd = _init_pymol()

    gene_dirs = sorted(
        d for d in af3_dir.iterdir()
        if d.is_dir()
        and d.name not in ("data", "__pycache__")
        and not _AF3_TIMESTAMP_RE.match(d.name)
    )

    # ── Load reference orientation (mirrors Op 5) ─────────────────────────────
    ref_cif = af3_dir / _REF_CIF_FILE
    ref_pse = af3_dir / _REF_PSE_FILE
    ref_rotation = None
    if _ORIENTATION_METHOD == "reference" and ref_pse.exists():
        try:
            cmd.load(str(ref_pse))
            ref_view = list(cmd.get_view())
            ref_rotation = ref_view[0:9]
            print(f"  Loaded reference rotation from {ref_pse.name}")
            cmd.reinitialize()
        except Exception as e:
            print(f"  Warning: could not load reference.pse: {e}")

    ref_cif_arg = ref_cif if ref_cif.exists() else None

    # ── Pass 1: compute uniform camera Z across all AF3 structures ────────────
    print("  Pass 1: computing uniform zoom ...")
    cam_z_vals = []
    for gene_dir in gene_dirs:
        gene_name = gene_dir.name
        af3_pdb   = gene_dir / f"{gene_name}.pdb"
        if not af3_pdb.exists():
            continue
        cmd.reinitialize()
        cmd.load(str(af3_pdb), gene_name)
        cmd.hide("everything")
        cmd.show("cartoon", gene_name)
        _apply_orientation(cmd, gene_name, ref_cif_arg, ref_rotation, None, None)
        cmd.zoom(gene_name, buffer=_ZOOM_BUFFER)
        v = cmd.get_view()
        cam_z_vals.append(v[11])
        cmd.reinitialize()

    if cam_z_vals:
        uniform_cam_z  = min(cam_z_vals)
        cam_dist       = abs(uniform_cam_z)
        uniform_slab   = (cam_dist - _SLAB_DEPTH, cam_dist + _SLAB_DEPTH)
        print(f"  Uniform camera Z = {uniform_cam_z:.1f}, "
              f"slab = {uniform_slab[0]:.1f} .. {uniform_slab[1]:.1f}")
    else:
        uniform_cam_z = None
        uniform_slab  = None

    def orient_fn(gene: str) -> None:
        _apply_orientation(cmd, gene, ref_cif_arg, ref_rotation,
                           uniform_cam_z, uniform_slab)

    # ── Render loop ───────────────────────────────────────────────────────────
    processed = 0
    for gene_dir in gene_dirs:
        gene_name  = gene_dir.name
        swiss_name = af3_to_swiss_name(gene_name)

        af3_pdb        = gene_dir / f"{gene_name}.pdb"
        swiss_gene_dir = swiss_dir / swiss_name if swiss_dir.is_dir() else None
        swiss_pdbs     = sorted(swiss_gene_dir.glob("*_model_*.pdb")) \
                         if swiss_gene_dir and swiss_gene_dir.is_dir() else []

        if not af3_pdb.exists():
            print(f"  {gene_name}: AF3 PDB not found — skipping")
            continue
        if not swiss_pdbs and needs_swiss:
            print(f"  {gene_name}: no SWISS match for '{swiss_name}' — "
                  "skipping overlay/deviation")

        gene_out = comp_dir / gene_name
        first_bg = (backgrounds or _DEFAULT_BACKGROUNDS)[0]
        did_render = False

        if "overlay" in render_types and swiss_pdbs:
            marker = gene_out / f"overlay_{first_bg}.jpg"
            if not marker.exists() or overwrite:
                os.makedirs(gene_out, exist_ok=True)
                print(f"  {gene_name}: overlay")
                render_overlay(cmd, af3_pdb, swiss_pdbs[0], gene_out,
                               backgrounds, orient_fn)
            did_render = True

        if "confidence" in render_types:
            marker = gene_out / f"confidence_{first_bg}.jpg"
            if not marker.exists() or overwrite:
                os.makedirs(gene_out, exist_ok=True)
                print(f"  {gene_name}: confidence")
                render_confidence(cmd, af3_pdb, gene_out,
                                  conf_colors or {}, backgrounds, orient_fn)
            did_render = True

        if "deviation" in render_types and swiss_pdbs:
            marker = gene_out / f"deviation_{first_bg}.jpg"
            if not marker.exists() or overwrite:
                os.makedirs(gene_out, exist_ok=True)
                print(f"  {gene_name}: deviation")
                render_deviation(cmd, af3_pdb, swiss_pdbs[0], gene_out,
                                 backgrounds, orient_fn)
            did_render = True

        if did_render:
            processed += 1

    return processed


def main():
    parser = argparse.ArgumentParser(description="Comparative overlay renders AF3 vs SWISS")
    parser.add_argument("--run-dir", required=True,
                        help="Genome run directory (contains AlphaFold3_Results/ and SWISS_Results/)")
    parser.add_argument("--overwrite", default="false",
                        help="Overwrite existing outputs (true/false)")
    parser.add_argument("--backgrounds", default="",
                        help="Space-separated background colours (e.g. 'black' or 'black white'). "
                             "Defaults to [rendering].backgrounds from the colors config.")
    parser.add_argument("--color-config", default=None,
                        help="Path to protein_structure_colors.toml. "
                             "Rendering and overlay settings are read from [rendering] and "
                             "[comparison] sections respectively.")
    parser.add_argument("--render-types", default="",
                        help="Space-separated Op 8 sub-types to run: "
                             "overlay confidence deviation. "
                             "Defaults to 'overlay' when omitted.")
    args = parser.parse_args()

    # ── Load settings from colors config ─────────────────────────────────────
    global IMAGE_WIDTH, IMAGE_HEIGHT, DPI, JPEG_QUALITY, ANTIALIAS, \
           RAY_SHADOWS, SPECULAR, DEPTH_CUE, _DEFAULT_BACKGROUNDS, \
           _AF3_COLOR, _SWISS_COLOR, _SWISS_TRANSP, _ZOOM_BUFFER, \
           _ORIENTATION_METHOD, _REF_PSE_FILE, _REF_CIF_FILE, _HELIX_UP, _SLAB_DEPTH
    conf_colors: dict = {}
    if args.color_config:
        _cfg_path = Path(args.color_config)
        if _cfg_path.exists():
            with open(_cfg_path, "rb") as _f:
                _ccfg = tomllib.load(_f)

            _r = _ccfg.get("rendering", {})
            IMAGE_WIDTH         = _r.get("width",               IMAGE_WIDTH)
            IMAGE_HEIGHT        = _r.get("height",              IMAGE_HEIGHT)
            DPI                 = _r.get("dpi",                 DPI)
            JPEG_QUALITY        = _r.get("jpeg_quality",        JPEG_QUALITY)
            ANTIALIAS           = _r.get("antialias",           ANTIALIAS)
            RAY_SHADOWS         = 1 if _r.get("ray_shadows", bool(RAY_SHADOWS)) else 0
            SPECULAR            = _r.get("specular",            SPECULAR)
            DEPTH_CUE           = 1 if _r.get("depth_cue", bool(DEPTH_CUE)) else 0
            _ORIENTATION_METHOD = _r.get("orientation_method",  _ORIENTATION_METHOD)
            _REF_PSE_FILE       = _r.get("reference_pse_file",  _REF_PSE_FILE)
            _REF_CIF_FILE       = _r.get("reference_cif_file",  _REF_CIF_FILE)
            _HELIX_UP           = _r.get("helix_up",            _HELIX_UP)
            _SLAB_DEPTH         = _r.get("slab_depth",          _SLAB_DEPTH)
            _cfg_bgs = _r.get("backgrounds", None)
            if _cfg_bgs:
                _DEFAULT_BACKGROUNDS = tuple(_cfg_bgs)

            _ov = _ccfg.get("comparison", {})
            _AF3_COLOR    = _ov.get("af3_color",         _AF3_COLOR)
            _SWISS_COLOR  = _ov.get("swiss_color",        _SWISS_COLOR)
            _SWISS_TRANSP = _ov.get("swiss_transparency", _SWISS_TRANSP)
            _ZOOM_BUFFER  = _ov.get("zoom_buffer",        _ZOOM_BUFFER)
            conf_colors   = _ccfg.get("confidence", {})
        else:
            print(f"Warning: color config not found at {args.color_config}, using defaults")

    run_dir   = Path(args.run_dir)
    overwrite = args.overwrite.lower() in ("true", "1", "yes")
    # --backgrounds CLI arg overrides config default; empty string → use config default.
    if args.backgrounds.strip():
        backgrounds = tuple(b.strip() for b in args.backgrounds.split() if b.strip())
    else:
        backgrounds = _DEFAULT_BACKGROUNDS
    # --render-types CLI arg; default to "overlay" when absent.
    if args.render_types.strip():
        render_types = tuple(t.strip() for t in args.render_types.split() if t.strip())
    else:
        render_types = ("overlay",)

    n = render_all(run_dir, overwrite, backgrounds, render_types, conf_colors)
    print(f"  Rendered {n} gene(s) — types: {', '.join(render_types)}")


if __name__ == "__main__":
    main()
