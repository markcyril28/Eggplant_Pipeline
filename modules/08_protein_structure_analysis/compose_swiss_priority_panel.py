#!/usr/bin/env python3
"""Compose the 2x2 priority SmelDMP SWISS-MODEL panel.

Panels (v5 naming, priority order):
  (i)   SmelDMPv5_10.610  Tier 1  top-left
  (ii)  SmelDMPv5_01.030  Tier 1  top-right
  (iii) SmelDMPv5_02       Tier 2  bottom-left
  (iv)  SmelDMPv5_01.730  Tier 2  bottom-right

Each source image is cropped to its non-white bounding box, expanded to
a square with uniform padding so the structure occupies ~75% of each cell,
then all four cells are placed at identical sizes.

Output: SWISS_Results/priority_smeldmp_swiss_model_composite.png

Usage:
    python3 compose_swiss_priority_panel.py [--run-dir /path/to/GPE001970_SMEL5]
    python3 compose_swiss_priority_panel.py --cell-size 1800 --pad-frac 0.15
"""

from __future__ import annotations
import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

# ── Panel definitions ─────────────────────────────────────────────────────────
# Each entry: gene directory under SWISS_Results/, a glob pattern to find the
# white-background PNG, v5 display name, gene accession, tier label, tier color.
PANELS = [
    {
        "panel":      "i",
        "gene_dir":   "SMEL5_10g017610.1",
        "src_glob":   "*white*BG*.png",
        "name":       "SmelDMPv5_10.610",
        "accession":  "SMEL5_10g017610.1",
        "tier":       "Tier 1",
        "tier_color": (200, 150, 0),    # dark gold (readable on white)
    },
    {
        "panel":      "ii",
        "gene_dir":   "SMEL5_01g026030.1",
        "src_glob":   "*white*BG.png",   # excludes " - Copy.png"
        "name":       "SmelDMPv5_01.030",
        "accession":  "SMEL5_01g026030.1",
        "tier":       "Tier 1",
        "tier_color": (200, 150, 0),
    },
    {
        "panel":      "iii",
        "gene_dir":   "SMEL5_02g013320.1",
        "src_glob":   "*white*BG*.png",
        "name":       "SmelDMPv5_02",
        "accession":  "SMEL5_02g013320.1",
        "tier":       "Tier 2",
        "tier_color": (50, 130, 200),   # mid-blue (readable on white)
    },
    {
        "panel":      "iv",
        "gene_dir":   "SMEL5_01g008730.1",
        "src_glob":   "*white*BG*.png",
        "name":       "SmelDMPv5_01.730",
        "accession":  "SMEL5_01g008730.1",
        "tier":       "Tier 2",
        "tier_color": (50, 130, 200),
    },
]

# ── Layout constants ──────────────────────────────────────────────────────────
CELL_SIZE    = 1800    # pixels per panel (square)
LABEL_H      = 160     # height of label bar below each cell
CELL_GAP     = 12      # gap between cells
BG_COLOR     = (255, 255, 255)
TEXT_COLOR   = (30,  30,  30)
SUBTEXT_COLOR= (100, 100, 100)
BORDER_COLOR = (220, 220, 220)
BORDER_W     = 2
FONT_NAME_SZ = 52
FONT_ACC_SZ  = 38
FONT_TIER_SZ = 40
FONT_PANEL_SZ= 42
PAD_FRAC     = 0.12    # fraction of bbox side added as padding on each side
WHITE_THRESH = 240     # luminance >= this is "white background"
BLACK_THRESH = 15      # luminance <= this is "black background"


def get_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold
        else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size)
        except (OSError, IOError):
            continue
    return ImageFont.load_default()


def _bg_is_black(img: Image.Image) -> bool:
    """Return True if the image has a black background (corner pixels near 0)."""
    import numpy as np
    arr = np.array(img.convert("RGB"))
    corners = [arr[0, 0], arr[0, -1], arr[-1, 0], arr[-1, -1]]
    return float(sum(c.mean() for c in corners) / 4) < 30


def crop_to_content(img: Image.Image,
                    pad_frac: float = PAD_FRAC) -> Image.Image:
    """Crop to structure bounding box, expand to square, add padding.

    Auto-detects black vs white background.  Returns the cropped square
    placed on a white background so all panels share the same canvas color.
    """
    import numpy as np
    rgb = np.array(img.convert("RGB"))

    black_bg = _bg_is_black(img)
    if black_bg:
        # Content = any channel significantly above black
        mask_arr = (rgb > BLACK_THRESH).any(axis=2)
    else:
        # Content = any channel significantly below white
        mask_arr = (rgb < WHITE_THRESH).any(axis=2)

    rows = np.any(mask_arr, axis=1)
    cols = np.any(mask_arr, axis=0)
    if not rows.any():
        return img
    top    = int(np.argmax(rows))
    bottom = int(len(rows) - np.argmax(rows[::-1]))
    left   = int(np.argmax(cols))
    right  = int(len(cols) - np.argmax(cols[::-1]))

    content_h = bottom - top
    content_w = right  - left
    side = max(content_h, content_w)
    pad  = int(side * pad_frac)
    side += 2 * pad

    cx = (left + right)  // 2
    cy = (top  + bottom) // 2
    half = side // 2

    sq_l = cx - half
    sq_t = cy - half
    sq_r = sq_l + side
    sq_b = sq_t + side

    W, H = img.size
    if sq_l < 0:   sq_r -= sq_l;  sq_l = 0
    if sq_t < 0:   sq_b -= sq_t;  sq_t = 0
    if sq_r > W:   sq_l -= sq_r - W; sq_r = W
    if sq_b > H:   sq_t -= sq_b - H; sq_b = H
    sq_l = max(0, sq_l)
    sq_t = max(0, sq_t)

    cropped = img.crop((sq_l, sq_t, sq_r, sq_b))

    # For black-BG renders: replace background (near-black) pixels with white
    # so all panels share a consistent white canvas.
    if black_bg:
        cr_arr = np.array(cropped.convert("RGB"))
        bg_mask = ~(cr_arr > BLACK_THRESH).any(axis=2)  # True = background pixel
        cr_arr[bg_mask] = [255, 255, 255]
        return Image.fromarray(cr_arr, "RGB")

    return cropped


def find_source(swiss_dir: Path, gene_dir: str, glob_pat: str) -> Path | None:
    d = swiss_dir / gene_dir
    if not d.is_dir():
        return None
    matches = sorted(d.glob(glob_pat))
    if not matches:
        # fallback: any PNG
        matches = sorted(d.glob("*.png"))
    return matches[0] if matches else None


def compose(swiss_dir: Path, out_path: Path,
            cell_size: int = CELL_SIZE, pad_frac: float = PAD_FRAC) -> None:
    import numpy as np

    cell_w  = cell_size
    cell_h  = cell_size
    row_h   = cell_h + LABEL_H
    total_w = cell_w * 2 + CELL_GAP
    total_h = row_h  * 2 + CELL_GAP

    canvas = Image.new("RGB", (total_w, total_h), BG_COLOR)
    draw   = ImageDraw.Draw(canvas)

    fn  = get_font(FONT_NAME_SZ,  bold=True)
    fa  = get_font(FONT_ACC_SZ,   bold=False)
    ft  = get_font(FONT_TIER_SZ,  bold=True)
    fp  = get_font(FONT_PANEL_SZ, bold=False)

    positions = [
        (0,                       0),
        (cell_w + CELL_GAP,       0),
        (0,                       row_h + CELL_GAP),
        (cell_w + CELL_GAP,       row_h + CELL_GAP),
    ]

    for p, (px, py) in zip(PANELS, positions):
        src_path = find_source(swiss_dir, p["gene_dir"], p["src_glob"])

        # ── Place structure image ──────────────────────────────────────────
        if src_path is None or not src_path.exists():
            print(f"  WARNING: source not found for {p['name']} "
                  f"(dir={p['gene_dir']}, glob={p['src_glob']})")
            placeholder = Image.new("RGB", (cell_w, cell_h), (240, 240, 240))
            canvas.paste(placeholder, (px, py))
        else:
            print(f"  {p['name']}: loading {src_path.name}")
            raw = Image.open(src_path).convert("RGB")
            cropped = crop_to_content(raw, pad_frac=pad_frac)
            cell_img = cropped.resize((cell_w, cell_h), Image.LANCZOS)
            canvas.paste(cell_img, (px, py))

        # ── Optional thin border ───────────────────────────────────────────
        if BORDER_W > 0:
            draw.rectangle(
                [px, py, px + cell_w - 1, py + cell_h - 1],
                outline=BORDER_COLOR, width=BORDER_W,
            )

        # ── Label bar ─────────────────────────────────────────────────────
        label_y = py + cell_h
        draw.rectangle([px, label_y, px + cell_w, label_y + LABEL_H], fill=BG_COLOR)

        # Panel letter "(i)" at far left
        panel_str = f"({p['panel']})"
        pb = draw.textbbox((0, 0), panel_str, font=fp)
        panel_w = pb[2] - pb[0]

        # Gene name
        nb = draw.textbbox((0, 0), p["name"], font=fn)
        name_h = nb[3] - nb[1]

        # Accession
        ab = draw.textbbox((0, 0), p["accession"], font=fa)
        acc_h = ab[3] - ab[1]

        # Tier
        tb = draw.textbbox((0, 0), p["tier"], font=ft)
        tier_w = tb[2] - tb[0]
        tier_h = tb[3] - tb[1]

        total_text_h = name_h + 10 + acc_h
        text_top = label_y + (LABEL_H - total_text_h) // 2

        x_name = px + 16 + panel_w + 12
        draw.text((px + 16, text_top + 4), panel_str, fill=SUBTEXT_COLOR, font=fp)
        draw.text((x_name, text_top), p["name"], fill=TEXT_COLOR, font=fn)
        draw.text((x_name, text_top + name_h + 10), p["accession"],
                  fill=SUBTEXT_COLOR, font=fa)

        # Tier at right
        tier_x = px + cell_w - tier_w - 16
        tier_y = label_y + (LABEL_H - tier_h) // 2
        draw.text((tier_x, tier_y), p["tier"], fill=p["tier_color"], font=ft)

    canvas.save(str(out_path), "PNG", optimize=False)
    print(f"  Saved: {out_path}")


def main() -> None:
    global CELL_SIZE, PAD_FRAC
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        default="III_RESULT/DMP/08_Protein_Structure/GPE001970_SMEL5",
        help="Path to genome run directory (contains SWISS_Results/)",
    )
    parser.add_argument(
        "--cell-size", type=int, default=CELL_SIZE,
        help=f"Pixel size of each panel cell (square, default {CELL_SIZE})",
    )
    parser.add_argument(
        "--pad-frac", type=float, default=PAD_FRAC,
        help=f"Padding fraction around cropped structure (default {PAD_FRAC})",
    )
    args = parser.parse_args()
    CELL_SIZE = args.cell_size
    PAD_FRAC  = args.pad_frac

    run_dir    = Path(args.run_dir).resolve()
    swiss_dir  = run_dir / "SWISS_Results"
    out_path   = swiss_dir / "priority_smeldmp_swiss_model_composite.png"

    if not swiss_dir.is_dir():
        print(f"ERROR: SWISS_Results not found at {swiss_dir}")
        return

    print(f"Composing SWISS-MODEL priority panel from: {swiss_dir}")
    compose(swiss_dir, out_path, cell_size=CELL_SIZE, pad_frac=PAD_FRAC)
    print("Done.")


if __name__ == "__main__":
    main()
