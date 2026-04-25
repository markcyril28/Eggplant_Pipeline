#!/usr/bin/env python3
"""
Render HAP2 PDB structures as publication-quality JPEG images using PyMOL.
Two versions per structure: black background and white background.

Orientation: PCA-based helix-vertical / sheet-forward.

Coloring strategy (per-chain):
    TM helices   = 4 longest helix segments per chain
    other helices = remaining shorter helices
    beta sheets  = beta sheets / strands
    loops        = loops / coils

Colors are read from config/colors_config/protein_structure_colors.toml
under the [HAP2] section.  Falls back to built-in defaults when no config
is given.

Usage (standalone):
    python3 render_structures_hap2.py --input-dir /path/to/Protein_Structure_HAP2

Usage (orchestrated):
    python3 render_structures_hap2.py \
        --input-dir /path/to/Protein_Structure_HAP2 \
        --color-config /path/to/protein_structure_colors.toml \
        --gene-group HAP2
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
        description="Render HAP2 PDB structures as JPEG images via PyMOL."
    )
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Directory with gene sub-folders containing .pdb files.")
    p.add_argument("--color-config", type=Path, default=None,
                   help="Path to protein_structure_colors.toml.")
    p.add_argument("--gene-group", type=str, default="HAP2",
                   help="Gene-group key inside the color TOML (default: HAP2).")
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
    p.add_argument("--tbtools-heatmap", action="store_true",
                   help="Extract structural features (helix/sheet/loop counts, residues, chains) "
                        "per gene during Pass 1 and write a TBTools-style gene×feature matrix TSV "
                        "plus a publication JPEG heatmap to Structural_Feature_Matrix/.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(color_toml, gene_group):
    """Return (render_cfg, palette, confidence) dicts."""
    if color_toml and color_toml.exists():
        with open(color_toml, "rb") as f:
            cfg = tomllib.load(f)
        render_cfg = cfg.get("rendering", {})
        palette = cfg.get(gene_group, {})
        confidence = cfg.get("confidence", {})
        print(f"Loaded color config from {color_toml}  [{gene_group}]")
    else:
        if color_toml:
            print(f"Warning: color config not found at {color_toml}, using defaults")
        render_cfg, palette, confidence = {}, {}, {}
    return render_cfg, palette, confidence


def resolve_color(value, color_name):
    """Register an RGB list as a custom PyMOL color; return the color name."""
    if isinstance(value, list) and len(value) == 3:
        cmd.set_color(color_name, value)
        return color_name
    return value  # already a PyMOL named color string


def _build_palette(base_palette, version_tag, conf=None):
    """Build a flat palette dict for *version_tag* ('original' or 'contrast').

    For 'confidence', return the pLDDT threshold colour dict instead.
    """
    if version_tag == "confidence":
        return conf or {}
    if version_tag == "original":
        src = base_palette
    else:
        src = base_palette.get(version_tag, base_palette)
    return {
        "tm":    src.get("tm_helices",    "red"),
        "helix": src.get("other_helices", "pink"),
        "sheet": src.get("beta_sheets",   "yellow"),
        "loop":  src.get("loops",         "gray70"),
    }


def _register_palette(pal, tag=""):
    """Register any RGB-list colors in *pal* as PyMOL custom colors."""
    pfx = f"{tag}_" if tag else ""
    resolve_color(pal["tm"],    f"{pfx}cfg_tm")
    resolve_color(pal["helix"], f"{pfx}cfg_helix")
    resolve_color(pal["sheet"], f"{pfx}cfg_sheet")
    resolve_color(pal["loop"],  f"{pfx}cfg_loop")


def _resolved(pal, key, tag=""):
    """Return the PyMOL color name for *key* (custom name for RGB lists)."""
    val = pal[key]
    pfx = f"{tag}_" if tag else ""
    if isinstance(val, list) and len(val) == 3:
        return f"{pfx}cfg_{key}"
    return val


# ---------------------------------------------------------------------------
# Structural feature extraction  (--tbtools-heatmap)
# ---------------------------------------------------------------------------

def _extract_ss_metrics(selection: str) -> dict:
    """Extract secondary-structure counts from *selection* (loaded in PyMOL).

    Returns a dict of integer-valued features suitable for a TBTools-style
    gene × feature matrix (analogous to a PlantCARE motif-count matrix).
    Helix classification mirrors color_by_structure: the 4 longest helix
    segments **per chain** are counted as TM helices; the rest as other helices.
    """
    chains = cmd.get_chains(selection) or [""]
    n_chains = len(chains)

    empty = {
        "Chains": n_chains, "Residues": 0, "TM_Helices": 0,
        "Other_Helices": 0, "Beta_Segments": 0, "Helix_Residues": 0,
        "Sheet_Residues": 0, "Loop_Residues": 0, "Longest_Helix_aa": 0,
    }

    def _chain_segments(chain_data, ss_char):
        """Split *chain_data* (sorted (resi,ss) pairs) into contiguous runs of *ss_char*.
        Uses gap > 1 to match color_by_structure segmentation."""
        segs, cur = [], []
        for resi, ss in chain_data:
            if ss == ss_char:
                if cur and resi - cur[-1] > 1:
                    segs.append(cur)
                    cur = [resi]
                else:
                    cur.append(resi)
            else:
                if cur:
                    segs.append(cur)
                    cur = []
        if cur:
            segs.append(cur)
        return segs

    total_residues  = 0
    total_helix_res = 0
    total_sheet_res = 0
    total_tm        = 0
    total_other_hx  = 0
    total_beta_segs = 0
    longest_helix   = 0

    for chain in chains:
        chain_sel = f"{selection} and chain {chain}" if chain else selection
        stored.resi_ss_ext = []
        cmd.iterate(
            f"({chain_sel}) and name CA",
            "stored.resi_ss_ext.append((int(resi), ss))",
        )
        if not stored.resi_ss_ext:
            continue

        chain_data = sorted(stored.resi_ss_ext, key=lambda x: x[0])
        total_residues  += len(chain_data)
        total_helix_res += sum(1 for _, ss in chain_data if ss == "H")
        total_sheet_res += sum(1 for _, ss in chain_data if ss == "S")

        helix_segs = sorted(_chain_segments(chain_data, "H"), key=len, reverse=True)
        beta_segs  = _chain_segments(chain_data, "S")

        n_tm_chain = min(4, len(helix_segs))  # top 4 per chain, matching color_by_structure
        total_tm        += n_tm_chain
        total_other_hx  += max(0, len(helix_segs) - n_tm_chain)
        total_beta_segs += len(beta_segs)
        if helix_segs:
            longest_helix = max(longest_helix, len(helix_segs[0]))

    if total_residues == 0:
        return empty

    return {
        "Chains":           n_chains,
        "Residues":         total_residues,
        "TM_Helices":       total_tm,
        "Other_Helices":    total_other_hx,
        "Beta_Segments":    total_beta_segs,
        "Helix_Residues":   total_helix_res,
        "Sheet_Residues":   total_sheet_res,
        "Loop_Residues":    total_residues - total_helix_res - total_sheet_res,
        "Longest_Helix_aa": longest_helix,
    }


def _write_structural_matrix(base_dir: Path, gene_metrics: dict, gene_group: str) -> None:
    """Write gene × structural-feature matrix TSV and a publication JPEG heatmap.

    Output layout (mirrors PlantCARE TBTools heatmap convention):
      Structural_Feature_Matrix/{gene_group}_structural_matrix.tsv
      Structural_Feature_Matrix/{gene_group}_structural_heatmap.jpg
    Rows = genes (sorted), columns = structural features, cells = integer counts.
    Column-wise normalisation is applied only for coloring; raw counts are shown
    as text inside cells (matching the style of the attached reference image).
    """
    out_dir = base_dir / "Structural_Feature_Matrix"
    out_dir.mkdir(exist_ok=True)

    genes = sorted(gene_metrics.keys())
    if not genes:
        print("  No metrics collected — skipping structural matrix.")
        return
    features = list(next(iter(gene_metrics.values())).keys())

    # ── TSV ──────────────────────────────────────────────────────────────
    tsv_path = out_dir / f"{gene_group.lower()}_structural_matrix.tsv"
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(["Gene"] + features) + "\n")
        for gene in genes:
            row = [gene] + [str(gene_metrics[gene].get(f, 0)) for f in features]
            fh.write("\t".join(row) + "\n")
    print(f"  Structural matrix  → {tsv_path.relative_to(base_dir)}")

    # ── Heatmap ───────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available — skipping heatmap.")
        return

    values   = np.array(
        [[gene_metrics[g].get(f, 0) for f in features] for g in genes],
        dtype=float,
    )
    col_max  = values.max(axis=0)
    col_max[col_max == 0] = 1          # avoid divide-by-zero for constant-zero columns
    norm_vals = values / col_max

    cell_w = 1.25
    cell_h = 0.80
    fig_w  = max(8,  len(features) * cell_w + 3.0)
    fig_h  = max(4,  len(genes)    * cell_h + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_facecolor("#f0f0f0")

    im = ax.imshow(norm_vals, aspect="auto", cmap="YlGn", vmin=0, vmax=1,
                   interpolation="nearest")

    # Axis ticks
    ax.set_xticks(range(len(features)))
    ax.set_xticklabels(
        [f.replace("_", "\n") for f in features],
        rotation=45, ha="right", fontsize=9, fontweight="bold",
    )
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([g.upper() for g in genes], fontsize=9)
    ax.set_xlabel("Structural Features", fontsize=11, labelpad=8)
    ax.set_ylabel("Genes of Interest",   fontsize=11, labelpad=8)
    ax.set_title(f"{gene_group} — Structural Feature Profile", fontsize=13, pad=12)

    # Integer count labels inside cells
    for i, gene in enumerate(genes):
        for j, feat in enumerate(features):
            raw        = int(gene_metrics[gene].get(feat, 0))
            text_color = "white" if norm_vals[i, j] > 0.55 else "black"
            ax.text(j, i, str(raw), ha="center", va="center",
                    fontsize=8.5, fontweight="bold", color=text_color)

    # Cell grid lines (minor ticks)
    ax.set_xticks(np.arange(-0.5, len(features), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes),    1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Relative Frequency", rotation=270, labelpad=14, fontsize=9)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0", "0.25", "0.50", "0.75", "1.00"])

    plt.tight_layout()
    png_path = out_dir / f"{gene_group.lower()}_structural_heatmap.png"
    jpg_path = out_dir / f"{gene_group.lower()}_structural_heatmap.jpg"
    fig.savefig(str(png_path), dpi=300, bbox_inches="tight")
    img = Image.open(str(png_path)).convert("RGB")
    img.save(str(jpg_path), "JPEG", quality=95, dpi=(300, 300))
    png_path.unlink(missing_ok=True)
    plt.close(fig)
    print(f"  Structural heatmap → {jpg_path.relative_to(base_dir)}")


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
        pdb_files = sorted(p for p in BASE.rglob("*.pdb") if p.name != "reference.pdb")
    if not pdb_files:
        print("No PDB files found — nothing to render.")
        return

    # ── Load config ────────────────────────────────────────────────────────
    _render, _palette, _conf = load_config(args.color_config, args.gene_group)

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

    # Orientation parameters (TOML → fallback)
    ORIENTATION_METHOD = _render.get("orientation_method", "reference")
    REF_PSE_FILE       = _render.get("reference_pse_file", "reference.pse")
    REF_CIF_FILE       = _render.get("reference_cif_file", "reference.cif")
    HELIX_UP           = _render.get("helix_up", False)
    SLAB_DEPTH         = _render.get("slab_depth", 200)

    REFERENCE_CIF = BASE / REF_CIF_FILE
    REFERENCE_PSE = BASE / REF_PSE_FILE

    # ── Launch PyMOL headless and extract reference view ──────────────────
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

    # ── Build color versions ───────────────────────────────────────────────
    version_names = _palette.get("color_versions", {}).get("names", ["original"])
    if args.color_versions:
        requested = [v.strip() for v in args.color_versions.split() if v.strip()]
        version_names = [v for v in version_names if v in requested]
    palettes = {v: _build_palette(_palette, v, conf=_conf) for v in version_names}

    for vname, pal in palettes.items():
        if vname == "confidence":
            print(f"  [{vname}] pLDDT confidence colouring")
        else:
            print(f"  [{vname}] TM={pal['tm']}  helix={pal['helix']}  "
                  f"sheet={pal['sheet']}  loop={pal['loop']}")

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

    def color_by_structure(gene, pal, tag):
        """
        Per-chain coloring for consistency:
          loops color  → loops (base)
          sheet color  → beta sheets
          TM color     → 4 longest helix segments per chain
          helix color  → remaining shorter helices
        """
        cmd.dss(gene)

        COL_TM    = _resolved(pal, "tm",    tag)
        COL_HELIX = _resolved(pal, "helix", tag)
        COL_SHEET = _resolved(pal, "sheet", tag)
        COL_LOOP  = _resolved(pal, "loop",  tag)

        # Base color: loops for all (covers loops/coils)
        cmd.color(COL_LOOP, gene)
        # Sheets
        cmd.color(COL_SHEET, f"{gene} and ss s")

        # Process each chain independently so monomeric and trimeric get
        # the same per-chain TM helix assignment (top 4 per chain).
        chains = cmd.get_chains(gene)
        if not chains:
            chains = [""]

        for chain in chains:
            if chain:
                chain_sel = f"{gene} and chain {chain}"
            else:
                chain_sel = gene

            stored.resi_ss = []
            cmd.iterate(f"{chain_sel} and name CA",
                        "stored.resi_ss.append((int(resi), ss))")
            if not stored.resi_ss:
                cmd.color(COL_TM, f"{chain_sel} and ss h")
                continue

            helix_segments = []
            current = []
            for resi, ss in stored.resi_ss:
                if ss == 'H':
                    if current and resi - current[-1] > 1:
                        helix_segments.append(current)
                        current = [resi]
                    else:
                        current.append(resi)
                else:
                    if current:
                        helix_segments.append(current)
                        current = []
            if current:
                helix_segments.append(current)

            if not helix_segments:
                continue

            # Top 4 longest → TM helices, rest → other helices
            helix_segments.sort(key=len, reverse=True)
            for i, seg in enumerate(helix_segments):
                color = COL_TM if i < 4 else COL_HELIX
                cmd.color(color, f"{chain_sel} and resi {seg[0]}-{seg[-1]} and ss h")

    # -------------------------------------------------------------------
    # Orientation helpers
    # -------------------------------------------------------------------

    def _get_ca_coords(selection):
        stored.xyz = []
        cmd.iterate_state(1, f"({selection}) and name CA",
                           "stored.xyz.append((x, y, z))")
        if not stored.xyz:
            return None
        return np.array(stored.xyz)

    def _principal_axis(coords):
        centered = coords - coords.mean(axis=0)
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        return Vt[0]

    def orient_pca_fallback(gene):
        """PCA-based orientation: helix axis → Y (up/down per HELIX_UP), sheet → Z."""
        helix_ca = _get_ca_coords(f"{gene} and ss h")
        if helix_ca is None or len(helix_ca) < 4:
            cmd.orient(gene)
            _zoom_fit(gene)
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
        _zoom_fit(gene)

    def apply_orientation(gene):
        if ref_rotation is not None and REFERENCE_CIF.exists():
            try:
                cmd.load(str(REFERENCE_CIF), "ref_struct")
                cmd.super(gene, "ref_struct")
                cmd.delete("ref_struct")
                cmd.orient(gene)
                view = list(cmd.get_view())
                view[0:9] = list(ref_rotation)
                cmd.set_view(view)
                _zoom_fit(gene)
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

    UNIFORM_CAM_Z = None
    UNIFORM_SLAB  = None
    gene_metrics: dict = {}  # populated in Pass 1 when --tbtools-heatmap is set

    def _zoom_fit(gene):
        cmd.zoom(gene, buffer=0)
        # Center on structured core
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

    def render(pdb_path, bg_color, pal, version_tag):
        gene = pdb_path.parent.name
        bg_suffix = "black" if bg_color == "black" else "white"
        ver_suffix = f"_{version_tag}" if version_tag != "original" else ""
        out_png = pdb_path.parent / f"{gene}_{bg_suffix}{ver_suffix}.png"
        out_jpg = pdb_path.parent / f"{gene}_{bg_suffix}{ver_suffix}.jpg"

        cmd.reinitialize()

        cmd.load(str(pdb_path), gene)

        cmd.hide("everything")
        cmd.show("cartoon", gene)
        cmd.dss(gene)  # required for ss-based selections in orientation + coloring

        if version_tag == "confidence":
            color_by_confidence(gene, pal, version_tag)
        else:
            # Re-register custom colors after reinitialize
            _register_palette(pal, version_tag)
            color_by_structure(gene, pal, version_tag)
        apply_orientation(gene)

        cmd.bg_color(bg_color)

        cmd.set("ray_opaque_background", 1)
        cmd.set("antialias", ANTIALIAS)
        cmd.set("ray_shadows", RAY_SHADOWS)
        cmd.set("depth_cue", DEPTH_CUE)
        cmd.set("specular", SPECULAR)

        cmd.ray(IMAGE_WIDTH, IMAGE_HEIGHT)
        cmd.png(str(out_png), width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=DPI, quiet=1)

        img = Image.open(str(out_png)).convert("RGB")
        img.save(str(out_jpg), "JPEG", quality=JPEG_QUALITY, dpi=(DPI, DPI))
        out_png.unlink()

        print(f"  Saved: {out_jpg.relative_to(BASE)}")

    # -------------------------------------------------------------------
    # Pass 1: compute uniform zoom (skip in worker mode)
    # -------------------------------------------------------------------
    if args.pdb_subset and args.uniform_cam_z is not None:
        UNIFORM_CAM_Z = args.uniform_cam_z
        near, far = args.uniform_cam_slab.split(",")
        UNIFORM_SLAB = (float(near), float(far))
        print(f"  Worker: cam_z={UNIFORM_CAM_Z:.1f}, slab={UNIFORM_SLAB}")
    else:
        print("\nPass 1: computing uniform zoom ...")
        cam_z_vals = []

        for pdb in pdb_files:
            gene = pdb.parent.name
            cmd.reinitialize()
            cmd.load(str(pdb), gene)
            cmd.hide("everything")
            cmd.show("cartoon", gene)
            cmd.dss(gene)
            if args.tbtools_heatmap:
                gene_metrics[gene] = _extract_ss_metrics(gene)
            apply_orientation(gene)
            cmd.zoom(gene, buffer=ZOOM_BUFFER)
            v = cmd.get_view()
            cam_z_vals.append(v[11])
            cmd.delete(gene)

        UNIFORM_CAM_Z = min(cam_z_vals)
        cam_dist = abs(UNIFORM_CAM_Z)
        UNIFORM_SLAB = (cam_dist - SLAB_DEPTH, cam_dist + SLAB_DEPTH)
        print(f"  Uniform camera Z = {UNIFORM_CAM_Z:.1f}")
        print(f"  Clip slab = {UNIFORM_SLAB[0]:.1f} .. {UNIFORM_SLAB[1]:.1f}")

    # ── Structural feature matrix (data collected during Pass 1) ──────────
    if args.tbtools_heatmap and gene_metrics:
        print("\nGenerating structural feature matrix and heatmap ...")
        _write_structural_matrix(BASE, gene_metrics, args.gene_group)

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
                "--gene-group", args.gene_group,
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
    # Pass 2: render with uniform scale — all versions × all backgrounds
    # -------------------------------------------------------------------
    for vtag, pal in palettes.items():
        print(f"\n-- Rendering color version: {vtag} --")
        for pdb in pdb_files:
            gene = pdb.parent.name
            print(f"  {gene} ...")
            for bg in BACKGROUNDS:
                render(pdb, bg, pal, vtag)

    cmd.quit()
    print("\nDone. All images saved.")


if __name__ == "__main__":
    main()
