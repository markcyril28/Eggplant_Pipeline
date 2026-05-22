#!/usr/bin/env python3
"""Compose the 2×2 reference.png panel from four priority SmelDMP protein renders.

Panels (v5 naming):
  (i)   SmelDMPv5_10.610  Tier 1  (primary candidate)
  (ii)  SmelDMPv5_01.030  Tier 1
  (iii) SmelDMPv5_02       Tier 2
  (iv)  SmelDMPv5_01.730  Tier 2

Output: AlphaFold3_Results/reference.png

Usage:
    python3 compose_reference_panel.py [--run-dir /path/to/GPE001970_SMEL5]
"""

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

# ── Panel definitions (v5 gene names, authoritative) ────────────────────────
PANELS = [
    {
        "panel": "i",
        "subdir": "smel5_10g017610_1",
        "stem":   "smel5_10g017610_1_black",
        "name":   "SmelDMPv5_10.610",
        "tier":   "Tier 1",
        "tier_color": (255, 215, 0),   # gold
    },
    {
        "panel": "ii",
        "subdir": "smel5_01g026030_1",
        "stem":   "smel5_01g026030_1_black",
        "name":   "SmelDMPv5_01.030",
        "tier":   "Tier 1",
        "tier_color": (255, 215, 0),
    },
    {
        "panel": "iii",
        "subdir": "smel5_02g013320_1",
        "stem":   "smel5_02g013320_1_black",
        "name":   "SmelDMPv5_02",
        "tier":   "Tier 2",
        "tier_color": (180, 220, 255),  # light blue
    },
    {
        "panel": "iv",
        "subdir": "smel5_01g008730_1",
        "stem":   "smel5_01g008730_1_black",
        "name":   "SmelDMPv5_01.730",
        "tier":   "Tier 2",
        "tier_color": (180, 220, 255),
    },
]

# ── Layout constants ─────────────────────────────────────────────────────────
CELL_SIZE    = 600      # pixels; source renders (1200×1200) scaled to this
LABEL_H      = 80       # height of the label bar below each cell
CELL_GAP     = 4        # gap between cells (horizontal and vertical)
BG_COLOR     = (0, 0, 0)
TEXT_COLOR   = (255, 255, 255)
FONT_SIZE    = 30
PANEL_FONT_SIZE = 26
CROP_PAD_FRAC = 0.06    # margin kept around the protein when cropping (closeup)


def get_font(size: int) -> ImageFont.FreeTypeFont:
    for name in [
        "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "arial.ttf",
    ]:
        try:
            return ImageFont.truetype(name, size)
        except (OSError, IOError):
            continue
    return ImageFont.load_default()


def crop_to_content(img: Image.Image, threshold: int = 16,
                    pad_frac: float = CROP_PAD_FRAC) -> Image.Image:
    """Crop *img* to the bounding box of non-background (non-black) pixels,
    expand that box to a square (centred on the box) so the structure is not
    distorted when resized, then add a small padding margin.

    *threshold* is the per-pixel luminance below which a pixel counts as
    background. *pad_frac* is the margin added on each side as a fraction of
    the square side length. Returns the cropped square sub-image.
    """
    gray = img.convert("L")
    mask = gray.point(lambda v: 255 if v > threshold else 0)
    bbox = mask.getbbox()
    if bbox is None:                       # fully black — nothing to crop
        return img

    left, upper, right, lower = bbox
    side = max(right - left, lower - upper)
    side += 2 * int(side * pad_frac)
    cx, cy = (left + right) // 2, (upper + lower) // 2
    half = side // 2

    sq_left, sq_upper = cx - half, cy - half
    sq_right, sq_lower = sq_left + side, sq_upper + side

    # Shift the square back inside the image if it overhangs an edge.
    W, H = img.size
    if sq_left < 0:
        sq_right -= sq_left; sq_left = 0
    if sq_upper < 0:
        sq_lower -= sq_upper; sq_upper = 0
    if sq_right > W:
        sq_left -= (sq_right - W); sq_right = W
    if sq_lower > H:
        sq_upper -= (sq_lower - H); sq_lower = H
    sq_left, sq_upper = max(0, sq_left), max(0, sq_upper)

    return img.crop((sq_left, sq_upper, sq_right, sq_lower))


def compose(af3_dir: Path, out_path: Path) -> None:
    cell_w    = CELL_SIZE
    cell_h    = CELL_SIZE
    row_h     = cell_h + LABEL_H
    total_w   = cell_w * 2 + CELL_GAP
    total_h   = row_h * 2 + CELL_GAP

    canvas = Image.new("RGB", (total_w, total_h), BG_COLOR)
    draw   = ImageDraw.Draw(canvas)
    font   = get_font(FONT_SIZE)
    pfont  = get_font(PANEL_FONT_SIZE)

    positions = [
        (0, 0),            # panel (i)  top-left
        (cell_w + CELL_GAP, 0),        # panel (ii) top-right
        (0, row_h + CELL_GAP),         # panel (iii) bottom-left
        (cell_w + CELL_GAP, row_h + CELL_GAP),   # panel (iv) bottom-right
    ]

    for p, (px, py) in zip(PANELS, positions):
        src = af3_dir / p["subdir"] / f"{p['stem']}.jpg"
        if not src.exists():
            print(f"  WARNING: source image not found — {src}")
            placeholder = Image.new("RGB", (cell_w, cell_h), (20, 20, 20))
            canvas.paste(placeholder, (px, py))
        else:
            img = Image.open(src).convert("RGB")
            img = crop_to_content(img)          # closeup: drop black margin
            img = img.resize((cell_w, cell_h), Image.LANCZOS)
            canvas.paste(img, (px, py))

        # Label bar (below the image)
        label_y = py + cell_h
        label_rect = [px, label_y, px + cell_w, label_y + LABEL_H]
        draw.rectangle(label_rect, fill=BG_COLOR)

        # Build label: "(i)  SmelDMPv5_10.610    Tier 1"
        panel_str = f"({p['panel']})"
        name_str  = f"  {p['name']}"

        # Measure widths to position tier label at right side
        pb = draw.textbbox((0, 0), panel_str, font=pfont)
        nb = draw.textbbox((0, 0), name_str,  font=font)
        tb = draw.textbbox((0, 0), p["tier"], font=pfont)

        text_y = label_y + (LABEL_H - (nb[3] - nb[1])) // 2

        x_cursor = px + 12
        draw.text((x_cursor, text_y + 2), panel_str, fill=(180, 180, 180), font=pfont)
        x_cursor += pb[2] - pb[0]
        draw.text((x_cursor, text_y), name_str, fill=TEXT_COLOR, font=font)

        # Tier label — right-aligned within the cell
        tier_x = px + cell_w - (tb[2] - tb[0]) - 12
        draw.text((tier_x, text_y + 2), p["tier"], fill=p["tier_color"], font=pfont)

    canvas.save(str(out_path), "PNG")
    print(f"  Saved: {out_path}")


def main() -> None:
    global CROP_PAD_FRAC
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        default="III_RESULT/DMP/08_Protein_Structure/GPE001970_SMEL5",
        help="Path to genome run directory (contains AlphaFold3_Results/)",
    )
    parser.add_argument(
        "--pad-frac", type=float, default=CROP_PAD_FRAC,
        help="Margin kept around each protein when cropping. Smaller = tighter "
             f"closeup (default: {CROP_PAD_FRAC}).",
    )
    args = parser.parse_args()
    CROP_PAD_FRAC = args.pad_frac

    run_dir  = Path(args.run_dir).resolve()
    af3_dir  = run_dir / "AlphaFold3_Results"
    out_path = af3_dir / "reference.png"

    if not af3_dir.is_dir():
        print(f"ERROR: AlphaFold3_Results not found at {af3_dir}")
        return

    print(f"Composing reference panel from: {af3_dir}")
    compose(af3_dir, out_path)
    print("Done.")


if __name__ == "__main__":
    main()
