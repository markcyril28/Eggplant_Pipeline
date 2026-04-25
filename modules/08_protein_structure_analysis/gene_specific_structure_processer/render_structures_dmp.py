#!/usr/bin/env python3
"""
DMP-specific: Render all PDB structures as JPEG images using PyMOL.
Two color versions per structure (original + contrast), each on black and
white backgrounds.

Coloring: topology-aware partitioning of helices relative to the central
sheet region (TM helices, amphipathic helices, extracellular domain,
loops).  Colors read from protein_structure_colors.toml [DMP] section.

Orientation: structural superposition onto reference.cif using camera from
reference.pse — gives consistent orientation across all DMP proteins.
Fallback: PCA-based helix-vertical / sheet-forward.

Usage (standalone):
    python3 render_structures_dmp.py --input-dir /path/to/Protein_Structures

Usage (orchestrated):
    python3 render_structures_dmp.py \
        --input-dir /path/to/Protein_Structures \
        --color-config /path/to/protein_structure_colors.toml
"""

import argparse
import os
import sys
from pathlib import Path

os.environ["PYMOL_QUIET"] = "1"

import numpy as np
import pymol
from pymol import cmd, stored
from PIL import Image

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Render PDB structures as publication-quality JPEG images via PyMOL."
    )
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Directory with gene sub-folders containing .pdb files, "
                        "plus optional reference.cif / reference.pse.")
    p.add_argument("--color-config", type=Path, default=None,
                   help="Path to protein_structure_colors.toml. "
                        "If omitted, built-in defaults are used.")
    p.add_argument("--workers", type=int, default=1,
                   help="Parallel render subprocesses (default: 1).")
    p.add_argument("--pdb-subset", type=str, default=None,
                   help="Comma-separated PDB paths (worker mode: skip Pass 1).")
    p.add_argument("--uniform-cam-z", type=float, default=None,
                   help="Pre-computed camera Z (worker mode).")
    p.add_argument("--uniform-cam-slab", type=str, default=None,
                   help="Pre-computed slab 'near,far' (worker mode).")
    p.add_argument("--color-versions", type=str, default=None,
                   help="Space-separated color versions to render (e.g. 'original contrast'). "
                        "If omitted, all versions from the color config are rendered.")
    p.add_argument("--backgrounds", type=str, default=None,
                   help="Space-separated background colours (e.g. 'black' or 'black white'). "
                        "Overrides the color config [rendering].backgrounds when provided.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(color_toml: Path | None, gene_group: str):
    """Return (render_cfg, palette, topo, fallback, confidence) dicts.

    *palette* is the full gene-group section (includes sub-keys like
    'contrast' and 'color_versions').
    *confidence* is the [confidence] section (pLDDT threshold colours).
    """
    if color_toml and color_toml.exists():
        with open(color_toml, "rb") as f:
            cfg = tomllib.load(f)
        render_cfg = cfg.get("rendering", {})
        palette = cfg.get(gene_group, {})
        topo = palette.get("topology", {})
        fallback = topo.get("fallback", {})
        confidence = cfg.get("confidence", {})
        print(f"Loaded color config from {color_toml}")
    else:
        if color_toml:
            print(f"Warning: color config not found at {color_toml}, using defaults")
        render_cfg, palette, topo, fallback, confidence = {}, {}, {}, {}, {}
    return render_cfg, palette, topo, fallback, confidence


def _build_palette(base_palette: dict, version_tag: str, conf: dict | None = None) -> dict:
    """Build a flat colour dict for *version_tag* ('original' → base).

    For 'confidence', return the pLDDT threshold colour dict instead.
    """
    if version_tag == "confidence":
        return conf or {}
    src = base_palette if version_tag == "original" else \
          base_palette.get(version_tag, base_palette)
    fb_src = src.get("topology", {}).get("fallback", {})
    return {
        "tm":       src.get("transmembrane_helices", src.get("tm_helices",           "red")),
        "amphi":    src.get("N-terminal_domain",      src.get("amphipathic_helices",   "slate")),
        "extra":    src.get("extracellular_beta_pleated_sheet", src.get("extracellular_domain",  "yellow")),
        "loop":     src.get("loops",                 "gray70"),
        "fb_helix": fb_src.get("all_helices",        "red"),
        "fb_sheet": fb_src.get("all_sheets",         "yellow"),
    }


def _register_palette(pal: dict, tag: str = "") -> None:
    """Register RGB-list colours in PyMOL with a tag prefix."""
    pfx = f"{tag}_" if tag else ""
    for key, val in pal.items():
        if isinstance(val, list) and len(val) == 3:
            cmd.set_color(f"{pfx}cfg_{key}", val)


def _resolved(pal: dict, key: str, tag: str = ""):
    """Return a PyMOL-usable colour name for *key* (registered name or string)."""
    val = pal[key]
    if isinstance(val, list) and len(val) == 3:
        pfx = f"{tag}_" if tag else ""
        return f"{pfx}cfg_{key}"
    return val


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    BASE = args.input_dir.resolve()

    if not BASE.is_dir():
        print(f"ERROR: directory does not exist: {BASE}", file=sys.stderr)
        sys.exit(1)

    if args.pdb_subset:
        pdb_files = [Path(p) for p in args.pdb_subset.split(",")]
    else:
        # One canonical PDB per gene subdirectory: prefer {dir}/{dir}.pdb,
        # fall back to the first *model_N.pdb found.
        pdb_files = []
        for d in sorted(BASE.iterdir()):
            if not d.is_dir():
                continue
            gene_pdb = d / f"{d.name}.pdb"
            if gene_pdb.exists():
                pdb_files.append(gene_pdb)
            else:
                candidates = sorted(d.glob("*model_[0-9].pdb"))
                if candidates:
                    pdb_files.append(candidates[0])
    if not pdb_files:
        print("No PDB files found — nothing to render.")
        return

    # ── Load config ────────────────────────────────────────────────────────
    _render, _palette, _topo, _fallback, _conf = load_config(args.color_config, "DMP")

    # Rendering settings (TOML → fallback)
    IMAGE_WIDTH  = _render.get("width", 1200)
    IMAGE_HEIGHT = _render.get("height", 1200)
    DPI          = _render.get("dpi", 600)
    JPEG_QUALITY = _render.get("jpeg_quality", 95)
    RAY_SHADOWS  = 1 if _render.get("ray_shadows", False) else 0
    SPECULAR     = _render.get("specular", 0.3)
    ANTIALIAS    = _render.get("antialias", 2)
    DEPTH_CUE    = 1 if _render.get("depth_cue", True) else 0
    ZOOM_BUFFER  = _render.get("zoom_buffer", 20)
    BACKGROUNDS  = _render.get("backgrounds", ["black", "white"])
    if args.backgrounds:
        BACKGROUNDS = [b.strip() for b in args.backgrounds.split() if b.strip()]

    # Color versions ──────────────────────────────────────────────────────
    version_names = _palette.get("color_versions", {}).get("names", ["original"])
    if args.color_versions:
        requested = [v.strip() for v in args.color_versions.split() if v.strip()]
        version_names = [v for v in version_names if v in requested]
    palettes = {v: _build_palette(_palette, v, conf=_conf) for v in version_names}

    # Topology parameters (TOML → fallback)
    HELIX_MERGE_GAP  = _topo.get("helix_merge_gap", 2)
    TM_BEFORE_SHEETS = _topo.get("tm_before_sheets", 2)
    TM_AFTER_SHEETS  = _topo.get("tm_after_sheets", 2)

    # Orientation parameters (TOML → fallback)
    ORIENTATION_METHOD = _render.get("orientation_method", "reference")
    REF_PSE_FILE       = _render.get("reference_pse_file", "reference.pse")
    REF_CIF_FILE       = _render.get("reference_cif_file", "reference.cif")
    HELIX_UP           = _render.get("helix_up", False)
    SLAB_DEPTH         = _render.get("slab_depth", 200)

    REFERENCE_CIF = BASE / REF_CIF_FILE
    REFERENCE_PSE = BASE / REF_PSE_FILE

    # ── Launch PyMOL headless and extract reference view ───────────────────
    pymol.finish_launching(["pymol", "-cq"])

    ref_rotation = None
    if ORIENTATION_METHOD == "reference" and REFERENCE_PSE.exists():
        try:
            cmd.load(str(REFERENCE_PSE))
            ref_view = list(cmd.get_view())
            ref_rotation = ref_view[0:9]
            print(f"Reference rotation extracted from {REFERENCE_PSE.name}")
            cmd.reinitialize()
        except Exception as e:
            print(f"Warning: could not load reference.pse: {e}")

    # -------------------------------------------------------------------
    # Coloring helpers
    # -------------------------------------------------------------------

    def color_by_confidence(gene: str, conf: dict, tag: str = "") -> None:
        """Colour by AlphaFold pLDDT confidence (B-factor column)."""
        pfx = f"{tag}_" if tag else ""
        cmd.set_color(f"{pfx}plddt_vhigh", conf.get("very_high", [0.00, 0.33, 0.84]))
        cmd.set_color(f"{pfx}plddt_high",  conf.get("high",      [0.40, 0.80, 0.95]))
        cmd.set_color(f"{pfx}plddt_low",   conf.get("low",       [1.00, 0.86, 0.07]))
        cmd.set_color(f"{pfx}plddt_vlow",  conf.get("very_low",  [1.00, 0.49, 0.27]))
        thresholds = conf.get("thresholds", [50, 70, 90])
        cmd.color(f"{pfx}plddt_vlow",  gene)
        cmd.color(f"{pfx}plddt_low",   f"({gene}) and b > {thresholds[0]}")
        cmd.color(f"{pfx}plddt_high",  f"({gene}) and b > {thresholds[1]}")
        cmd.color(f"{pfx}plddt_vhigh", f"({gene}) and b > {thresholds[2]}")

    def color_by_structure(gene: str, pal: dict, tag: str = "") -> None:
        COL_TM    = _resolved(pal, "tm",    tag)
        COL_AMPHI = _resolved(pal, "amphi",  tag)
        COL_EXTRA = _resolved(pal, "extra",  tag)
        COL_LOOP  = _resolved(pal, "loop",   tag)
        FB_HELIX  = _resolved(pal, "fb_helix", tag)
        FB_SHEET  = _resolved(pal, "fb_sheet", tag)

        cmd.dss(gene)
        cmd.color(COL_LOOP, gene)

        stored.resi_ss = []
        cmd.iterate(f"{gene} and name CA",
                    "stored.resi_ss.append((int(resi), ss))")
        if not stored.resi_ss:
            return

        raw_segments = []
        current = []
        for resi, ss in stored.resi_ss:
            if ss == 'H':
                if current and resi - current[-1] > 1:
                    raw_segments.append(current[:])
                    current = [resi]
                else:
                    current.append(resi)
            else:
                if current:
                    raw_segments.append(current[:])
                    current = []
        if current:
            raw_segments.append(current[:])

        if not raw_segments:
            cmd.color(FB_SHEET, f"{gene} and ss s")
            return

        merged = [[raw_segments[0]]]
        for i in range(1, len(raw_segments)):
            if raw_segments[i][0] - raw_segments[i - 1][-1] <= HELIX_MERGE_GAP:
                merged[-1].append(raw_segments[i])
            else:
                merged.append([raw_segments[i]])

        groups = []
        for m in merged:
            total = sum(len(s) for s in m)
            min_r = min(s[0] for s in m)
            max_r = max(s[-1] for s in m)
            groups.append({'segs': m, 'len': total, 'min': min_r, 'max': max_r})

        groups.sort(key=lambda g: g['min'])

        sheet_resi = [r for r, s in stored.resi_ss if s == 'S']
        if not sheet_resi:
            cmd.color(FB_HELIX, f"{gene} and ss h")
            return

        sheet_min = min(sheet_resi)
        sheet_max = max(sheet_resi)

        before = [g for g in groups if g['max'] < sheet_min]
        after  = [g for g in groups if g['min'] > sheet_max]

        # Recover helix groups that straddle the sheet boundary — assign by centre
        sheet_center = (sheet_min + sheet_max) / 2.0
        classified_ids = {id(g) for g in before + after}
        for g in groups:
            if id(g) not in classified_ids:
                if (g['min'] + g['max']) / 2.0 <= sheet_center:
                    before.append(g)
                else:
                    after.append(g)
        before.sort(key=lambda g: g['min'])
        after.sort(key=lambda g: g['min'])

        if len(before) < TM_BEFORE_SHEETS or len(after) < TM_AFTER_SHEETS:
            cmd.color(FB_HELIX, f"{gene} and ss h")
            cmd.color(FB_SHEET, f"{gene} and ss s")
            return

        tm_before = before[-TM_BEFORE_SHEETS:]
        non_tm_before = before[:-TM_BEFORE_SHEETS]

        tm_after = after[:TM_AFTER_SHEETS]
        tm_groups = tm_before + tm_after

        for g in tm_groups:
            for seg in g['segs']:
                cmd.color(COL_TM, f"{gene} and resi {seg[0]}-{seg[-1]} and ss h")

        for g in non_tm_before:
            for seg in g['segs']:
                cmd.color(COL_AMPHI, f"{gene} and resi {seg[0]}-{seg[-1]} and ss h")

        yellow_start = tm_before[-1]['max'] + 1
        yellow_end = tm_after[0]['min'] - 1

        if yellow_start <= yellow_end:
            cmd.color(COL_EXTRA, f"{gene} and resi {yellow_start}-{yellow_end}")

        # Sweep: color remaining helix groups (C-terminal / split TM) as TM
        # so no helix stays at the default loop-gray colour
        classified_ids = {id(g) for g in tm_groups + non_tm_before}
        for g in groups:
            if id(g) not in classified_ids:
                for seg in g['segs']:
                    cmd.color(COL_TM,
                              f"{gene} and resi {seg[0]}-{seg[-1]} and ss h")

    # -------------------------------------------------------------------
    # Orientation helpers
    # -------------------------------------------------------------------

    def _get_ca_coords(selection: str):
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

    UNIFORM_CAM_Z = None
    UNIFORM_SLAB = None

    def _zoom_fit_all(gene: str) -> None:
        cmd.zoom(gene, buffer=0)

        core_sel = f"{gene} and (ss h or ss s) and name CA"
        core_ca = _get_ca_coords(core_sel)
        if core_ca is not None and len(core_ca) >= 4:
            core_center = core_ca.mean(axis=0)
            view = list(cmd.get_view())
            view[12] = core_center[0]
            view[13] = core_center[1]
            view[14] = core_center[2]
            cmd.set_view(view)

        if UNIFORM_CAM_Z is not None:
            view = list(cmd.get_view())
            view[11] = UNIFORM_CAM_Z
            view[15] = UNIFORM_SLAB[0]
            view[16] = UNIFORM_SLAB[1]
            cmd.set_view(view)

    def orient_pca_fallback(gene: str) -> None:
        helix_ca = _get_ca_coords(f"{gene} and ss h")
        if helix_ca is None or len(helix_ca) < 4:
            cmd.orient(gene)
            return

        all_ca = _get_ca_coords(gene)
        center = all_ca.mean(axis=0) if all_ca is not None else np.zeros(3)

        Y = _principal_axis(helix_ca)
        Y = Y / np.linalg.norm(Y)

        # Flip Y so it points toward sheets (ectodomain up on screen)
        sheet_ca = _get_ca_coords(f"{gene} and ss s")
        if sheet_ca is not None and len(sheet_ca) >= 3:
            sheet_dir = sheet_ca.mean(axis=0) - center
            if np.dot(Y, sheet_dir) < 0:
                Y = -Y
        else:
            if (Y[1] < 0) == HELIX_UP:
                Y = -Y
            centered = all_ca - center
            _, _, Vt = np.linalg.svd(centered, full_matrices=False)
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
        _zoom_fit_all(gene)

    def apply_orientation(gene: str) -> None:
        if ref_rotation is not None and REFERENCE_CIF.exists():
            try:
                cmd.load(str(REFERENCE_CIF), "ref_struct")
                cmd.super(gene, "ref_struct")
                cmd.delete("ref_struct")

                cmd.orient(gene)
                view = list(cmd.get_view())
                view[0:9] = list(ref_rotation)
                cmd.set_view(view)
                _zoom_fit_all(gene)
                return
            except Exception as e:
                print(f"  Warning: superposition failed ({e}), using PCA fallback")
                try:
                    cmd.delete("ref_struct")
                except Exception:
                    pass

        orient_pca_fallback(gene)

    # -------------------------------------------------------------------
    # Rendering
    # -------------------------------------------------------------------

    def render(pdb_path: Path, bg_color: str, pal: dict, version_tag: str) -> None:
        gene = pdb_path.parent.name
        suffix = "black" if bg_color == "black" else "white"
        ver_sfx = f"_{version_tag}" if version_tag != "original" else ""
        out_png = pdb_path.parent / f"{gene}_{suffix}{ver_sfx}.png"
        out_jpg = pdb_path.parent / f"{gene}_{suffix}{ver_sfx}.jpg"

        cmd.reinitialize()
        cmd.load(str(pdb_path), gene)

        cmd.hide("everything")
        cmd.show("cartoon", gene)
        cmd.dss(gene)  # required for ss-based selections in orientation + coloring

        if version_tag == "confidence":
            color_by_confidence(gene, pal, version_tag)
        else:
            _register_palette(pal, version_tag)
            color_by_structure(gene, pal, version_tag)
        apply_orientation(gene)

        cmd.bg_color(bg_color)

        cmd.set("ray_opaque_background", 1)
        cmd.set("antialias", ANTIALIAS)
        cmd.set("ray_shadows", RAY_SHADOWS)
        cmd.set("depth_cue", DEPTH_CUE)
        cmd.set("specular", SPECULAR)
        cmd.set("max_threads", 0)
        cmd.set("hash_max", 250)

        cmd.ray(IMAGE_WIDTH, IMAGE_HEIGHT)
        cmd.png(str(out_png), width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=DPI, quiet=1)

        img = Image.open(str(out_png)).convert("RGB")
        img.save(str(out_jpg), "JPEG", quality=JPEG_QUALITY, dpi=(DPI, DPI))
        out_png.unlink()

        print(f"  Saved: {out_jpg.relative_to(BASE)}")

    # -------------------------------------------------------------------
    # Pass 1: compute uniform camera distance (skip in worker mode)
    # -------------------------------------------------------------------
    if args.pdb_subset and args.uniform_cam_z is not None:
        UNIFORM_CAM_Z = args.uniform_cam_z
        near, far = args.uniform_cam_slab.split(",")
        UNIFORM_SLAB = (float(near), float(far))
        print(f"  Worker: cam_z={UNIFORM_CAM_Z:.1f}, slab={UNIFORM_SLAB}")
    else:
        print("Pass 1: computing uniform zoom (based on structured core) ...")
        cam_z_vals = []
        slab_vals = []

        for pdb in pdb_files:
            gene = pdb.parent.name
            cmd.reinitialize()
            cmd.load(str(pdb), gene)
            cmd.hide("everything")
            cmd.show("cartoon", gene)
            cmd.dss(gene)
            apply_orientation(gene)
            cmd.zoom(gene, buffer=ZOOM_BUFFER)
            v = cmd.get_view()
            cam_z_vals.append(v[11])
            slab_vals.append((v[15], v[16]))
            cmd.delete(gene)

        UNIFORM_CAM_Z = min(cam_z_vals)
        cam_dist = abs(UNIFORM_CAM_Z)
        UNIFORM_SLAB = (cam_dist - SLAB_DEPTH, cam_dist + SLAB_DEPTH)
        print(f"  Uniform camera Z = {UNIFORM_CAM_Z:.1f}")
        print(f"  Clip slab = {UNIFORM_SLAB[0]:.1f} .. {UNIFORM_SLAB[1]:.1f}")

    # -------------------------------------------------------------------
    # Fork workers (main process only, when --workers > 1)
    # -------------------------------------------------------------------
    if args.workers > 1 and not args.pdb_subset:
        import subprocess, math
        n = min(args.workers, len(pdb_files))
        chunk_size = math.ceil(len(pdb_files) / n)
        chunks = [pdb_files[i:i + chunk_size]
                  for i in range(0, len(pdb_files), chunk_size)]

        slab_str = f"{UNIFORM_SLAB[0]},{UNIFORM_SLAB[1]}"
        procs = []
        for idx, chunk in enumerate(chunks):
            subset = ",".join(str(p) for p in chunk)
            worker_cmd = [
                sys.executable, str(Path(__file__).resolve()),
                "--input-dir", str(BASE),
                "--workers", "1",
                "--uniform-cam-z", str(UNIFORM_CAM_Z),
                "--uniform-cam-slab", slab_str,
                "--pdb-subset", subset,
            ]
            if args.color_config:
                worker_cmd += ["--color-config", str(args.color_config)]
            if args.color_versions:
                worker_cmd += ["--color-versions", args.color_versions]
            if args.backgrounds:
                worker_cmd += ["--backgrounds", args.backgrounds]
            print(f"  Spawning worker {idx+1}/{len(chunks)} "
                  f"({len(chunk)} PDB files)")
            procs.append(subprocess.Popen(worker_cmd))

        failures = sum(1 for p in procs if p.wait() != 0)
        cmd.quit()
        if failures:
            print(f"Warning: {failures}/{len(procs)} render worker(s) failed")
            sys.exit(1)
        print(f"\nDone. All images saved ({n} workers).")
        return

    # -------------------------------------------------------------------
    # Pass 2: render with uniform scale (all versions × backgrounds)
    # -------------------------------------------------------------------
    for vtag, pal in palettes.items():
        print(f"\n── Color version: {vtag} ──")
        for pdb in pdb_files:
            gene = pdb.parent.name
            print(f"Rendering {gene} ...")
            for bg in BACKGROUNDS:
                render(pdb, bg, pal, vtag)

    cmd.quit()
    print("\nDone. All images saved.")


if __name__ == "__main__":
    main()
