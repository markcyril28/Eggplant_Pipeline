#!/usr/bin/env python3
"""
Compose side-by-side figures: pLDDT confidence (left) vs interaction domain (right)
for each DMP-HAP2 AlphaFold 3 complex structure.

Auto-discovers all subfolders under AlphaFold3_Results/ that contain
both *_black_confidence.jpg and *_black_interaction.jpg images.

Chain assignment:  A/B/C = HAP2 trimer,  D = DMP monomer

Outputs:
  - Individual per-complex composites  (Combined_Confidence_Interaction/)
  - Combined grid composite            (grid.jpg / grid.png)

CLI args (all optional):
  --af3-dir     Path to AlphaFold3_Results directory
  --output-dir  Output subdirectory name (default: Combined_Confidence_Interaction)
  --background  black | white (default: black)
"""

import argparse
import math
import re
import sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ── Default AlphaFold3 results directory (relative to pipeline root) ──────────
DEFAULT_AF3_DIR = Path(
    "3_RESULT/DMP-HAP2/08_Protein_Structure/"
    "GPE001970_SMEL5/AlphaFold3_Results"
)
DEFAULT_OUTPUT_SUBDIR = "Combined_Confidence_Interaction"

# ── Layout constants ──────────────────────────────────────────────────────────
BG_COLOR           = (0, 0, 0)        # black canvas
TEXT_COLOR          = (255, 255, 255)  # white title text
SUBTITLE_COLOR      = (200, 200, 200)  # light-gray subtitle
IMAGE_SCALE        = 1.3              # scale factor applied to source protein renders
TITLE_H            = 135              # pixels for title + subtitles
LEGEND_H           = 165              # pixels for legend below panels (two rows)
GAP                = 30               # horizontal gap between panels
GRID_GAP           = 40               # gap between rows/cols in grid
FONT_SIZE          = 44
SUBTITLE_FONT_SIZE = 30
LEGEND_FONT_SIZE   = 28
LEGEND_BOX_SIZE    = 24
LEGEND_SPACING     = 13
LEGEND_ROW_SPACING = 10               # vertical gap between legend rows

# ── pLDDT Confidence legend (AlphaFold standard) ─────────────────────────────
PLDDT_LEGEND = [
    ("Very high (>90)",   (0, 84, 214)),       # deep blue
    ("Confident (70\u201390)", (102, 204, 242)),   # cyan
    ("Low (50\u201370)",       (255, 219, 18)),    # yellow
    ("Very low (<50)",    (255, 125, 69)),      # orange
]

# ── DMP-HAP2 Interaction Domain legend ────────────────────────────────────────
# HAP2 chains (A/B/C) — cool palette
# DMP chain (D) — warm palette
INTERACTION_LEGEND_HAP2 = [
    ("HAP2 TM helices",     (0, 140, 140)),      # deep teal   [0.00, 0.55, 0.55]
    ("HAP2 other helices",  (128, 89, 194)),      # periwinkle  [0.50, 0.35, 0.76]
    ("HAP2 \u03b2-sheets",       (92, 135, 191)),     # steel blue  [0.36, 0.53, 0.75]
    ("HAP2 loops",          (179, 179, 179)),     # gray70
]

INTERACTION_LEGEND_DMP = [
    ("DMP TM helices",             (212, 48, 5)),    # deep orange-red
    ("DMP N-terminal domain",      (217, 166, 33)),   # golden amber
    ("DMP extracellular \u03B2-sheet", (245, 222, 179)),  # beige
    ("DMP loops",                  (179, 179, 179)),  # gray70
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--af3-dir", default=None,
        help="Path to AlphaFold3_Results directory (auto-detected from pipeline root if omitted)"
    )
    p.add_argument(
        "--output-dir", default=DEFAULT_OUTPUT_SUBDIR,
        dest="output_dir",
        help=f"Output sub-directory under AF3 dir (default: {DEFAULT_OUTPUT_SUBDIR})"
    )
    p.add_argument(
        "--background", default="black",
        help="Background colour of source renders (default: black)"
    )
    return p.parse_args()


def get_font(size):
    """Try to load a TrueType font, fall back to PIL default."""
    for name in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "C:/Windows/Fonts/arialbd.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "arial.ttf", "Arial.ttf", "DejaVuSans.ttf",
    ]:
        try:
            return ImageFont.truetype(name, size)
        except (OSError, IOError):
            continue
    return ImageFont.load_default()


def measure_legend_row(draw, items, font):
    """Measure the total width of one row of legend items."""
    total_w = 0
    for label, _colour in items:
        bbox = draw.textbbox((0, 0), label, font=font)
        tw = bbox[2] - bbox[0]
        total_w += LEGEND_BOX_SIZE + LEGEND_SPACING + tw
    total_w += LEGEND_SPACING * max(0, len(items) - 1)
    return total_w


def draw_legend_row(draw, items, x_start, y_center, max_width, font):
    """Draw a horizontal colour legend centred within *max_width* starting at *x_start*."""
    total_w = measure_legend_row(draw, items, font)
    x = x_start + max(0, (max_width - total_w)) // 2
    box_y = y_center - LEGEND_BOX_SIZE // 2

    for label, colour in items:
        draw.rectangle(
            [x, box_y, x + LEGEND_BOX_SIZE, box_y + LEGEND_BOX_SIZE],
            fill=colour, outline=(180, 180, 180),
        )
        draw.text(
            (x + LEGEND_BOX_SIZE + LEGEND_SPACING, box_y - 1),
            label, fill=SUBTITLE_COLOR, font=font,
        )
        bbox = draw.textbbox((0, 0), label, font=font)
        tw = bbox[2] - bbox[0]
        x += LEGEND_BOX_SIZE + LEGEND_SPACING + tw + LEGEND_SPACING


def extract_gene_name(folder_name):
    """Extract a readable gene name from folder like 'hap2_and_smel5_01g008730_1'."""
    m = re.search(r'smel5_(\d+g\d+)_(\d+)', folder_name)
    if m:
        return f"SMEL5_{m.group(1)}.{m.group(2)}"
    return folder_name


def extract_model_info(folder_name):
    """Extract model number if present (e.g. '_model_4')."""
    m = re.search(r'_model_(\d+)$', folder_name)
    if m:
        return int(m.group(1))
    return 0


def discover_pairs(af3_dir, bg):
    """Find all subfolders with confidence+interaction image pairs."""
    pairs = []
    for subdir in sorted(af3_dir.iterdir()):
        if not subdir.is_dir():
            continue
        if subdir.name.startswith("Combined") or subdir.name.startswith("reference"):
            continue

        conf_files = list(subdir.glob(f"*_{bg}_confidence.jpg"))
        inter_files = list(subdir.glob(f"*_{bg}_interaction.jpg"))

        if conf_files and inter_files:
            gene_name = extract_gene_name(subdir.name)
            model_num = extract_model_info(subdir.name)
            model_str = f" (model {model_num})" if model_num > 0 else ""
            pairs.append({
                "folder": subdir.name,
                "confidence": conf_files[0],
                "interaction": inter_files[0],
                "gene": gene_name,
                "model_str": model_str,
                "model_num": model_num,
            })
    return pairs


def make_pair(pair, out_dir, panel_idx):
    """Create a single side-by-side composite for one complex."""
    conf_img  = Image.open(pair["confidence"]).convert("RGB")
    inter_img = Image.open(pair["interaction"]).convert("RGB")

    # Scale up protein images for readability
    if IMAGE_SCALE != 1.0:
        sw = int(conf_img.width * IMAGE_SCALE)
        sh = int(conf_img.height * IMAGE_SCALE)
        conf_img = conf_img.resize((sw, sh), Image.LANCZOS)

    w, h = conf_img.size
    if inter_img.size != (w, h):
        inter_img = inter_img.resize((w, h), Image.LANCZOS)

    canvas_w = w * 2 + GAP
    canvas_h = h + TITLE_H + LEGEND_H
    canvas   = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)

    canvas.paste(conf_img,  (0,        TITLE_H))
    canvas.paste(inter_img, (w + GAP,  TITLE_H))

    draw    = ImageDraw.Draw(canvas)
    font    = get_font(FONT_SIZE)
    subfont = get_font(SUBTITLE_FONT_SIZE)
    legfont = get_font(LEGEND_FONT_SIZE)

    # ── Main title ────────────────────────────────────────────────────────
    panel_label = chr(ord('a') + panel_idx)
    title = (
        f"({panel_label})  HAP2\u2013{pair['gene']} Complex"
        f"{pair['model_str']}  \u2014  AlphaFold 3 Predicted Structure"
    )
    bbox = draw.textbbox((0, 0), title, font=font)
    draw.text(
        ((canvas_w - (bbox[2] - bbox[0])) / 2, 6),
        title, fill=TEXT_COLOR, font=font
    )

    # ── Panel subtitles ───────────────────────────────────────────────────
    conf_label  = "pLDDT Confidence"
    inter_label = "Interaction Domain"
    cb = draw.textbbox((0, 0), conf_label,  font=subfont)
    ib = draw.textbbox((0, 0), inter_label, font=subfont)
    subtitle_y = TITLE_H - SUBTITLE_FONT_SIZE - 10
    draw.text(
        ((w - (cb[2] - cb[0])) / 2, subtitle_y),
        conf_label, fill=SUBTITLE_COLOR, font=subfont
    )
    draw.text(
        (w + GAP + (w - (ib[2] - ib[0])) / 2, subtitle_y),
        inter_label, fill=SUBTITLE_COLOR, font=subfont
    )

    # ── Chain info under interaction subtitle ─────────────────────────────
    chain_info = "Chains A/B/C = HAP2 trimer  |  Chain D = DMP"
    chain_font = get_font(LEGEND_FONT_SIZE)
    ci_bbox = draw.textbbox((0, 0), chain_info, font=chain_font)
    draw.text(
        (w + GAP + (w - (ci_bbox[2] - ci_bbox[0])) / 2, subtitle_y + SUBTITLE_FONT_SIZE + 2),
        chain_info, fill=(160, 160, 160), font=chain_font
    )

    # ── Legends (below panels) ────────────────────────────────────────────
    # Left: pLDDT legend (single row)
    legend_y1 = TITLE_H + h + LEGEND_H // 3
    draw_legend_row(draw, PLDDT_LEGEND, 0, legend_y1, w, legfont)

    # Right: Interaction legend (two rows — HAP2 top, DMP bottom)
    legend_y_hap2 = TITLE_H + h + LEGEND_H // 3 - LEGEND_ROW_SPACING
    legend_y_dmp  = TITLE_H + h + 2 * LEGEND_H // 3 + LEGEND_ROW_SPACING
    draw_legend_row(draw, INTERACTION_LEGEND_HAP2, w + GAP, legend_y_hap2, w, legfont)
    draw_legend_row(draw, INTERACTION_LEGEND_DMP,  w + GAP, legend_y_dmp,  w, legfont)

    # Save
    safe_name = pair["folder"].replace("/", "_")
    out_name = f"{safe_name}_confidence_vs_interaction.jpg"
    out_path = out_dir / out_name
    canvas.save(out_path, quality=95)
    # Also save PNG
    out_png = out_dir / out_name.replace(".jpg", ".png")
    canvas.save(out_png)
    print(f"  Saved: {out_path.name}")
    return canvas


def make_grid(composites, pairs, out_dir):
    """Arrange composites in a grid (auto-sized)."""
    if not composites:
        return

    n = len(composites)
    cols = min(2, n)
    rows = math.ceil(n / cols)

    cw, ch = composites[0].size
    grid_w = cols * cw + (cols - 1) * GRID_GAP
    top_h  = 110
    grid_h = rows * ch + (rows - 1) * GRID_GAP + top_h
    grid   = Image.new("RGB", (grid_w, grid_h), BG_COLOR)

    draw  = ImageDraw.Draw(grid)
    font  = get_font(FONT_SIZE + 4)
    title = (
        "DMP\u2013HAP2 Protein Complex Structures: "
        "pLDDT Confidence vs Interaction Domain Coloring"
    )
    bbox = draw.textbbox((0, 0), title, font=font)
    draw.text(
        ((grid_w - (bbox[2] - bbox[0])) / 2, 14),
        title, fill=TEXT_COLOR, font=font
    )

    # Subtitle
    subfont = get_font(SUBTITLE_FONT_SIZE)
    subtitle = "AlphaFold 3 predicted structures  |  Chains A/B/C = HAP2 trimer, Chain D = DMP monomer"
    sb = draw.textbbox((0, 0), subtitle, font=subfont)
    draw.text(
        ((grid_w - (sb[2] - sb[0])) / 2, 50),
        subtitle, fill=SUBTITLE_COLOR, font=subfont
    )

    for idx, comp in enumerate(composites):
        row, col = divmod(idx, cols)
        x = col * (cw + GRID_GAP)
        y = top_h + row * (ch + GRID_GAP)
        grid.paste(comp, (x, y))

    out_jpg = out_dir / "DMP-HAP2_confidence_vs_interaction_grid.jpg"
    out_png = out_dir / "DMP-HAP2_confidence_vs_interaction_grid.png"
    grid.save(out_jpg, quality=95)
    grid.save(out_png)
    print(f"  Saved grid: {out_jpg.name}")
    print(f"  Saved grid: {out_png.name}")


def main():
    args = parse_args()

    # Resolve AF3 directory
    if args.af3_dir:
        af3_dir = Path(args.af3_dir).resolve()
    else:
        script_dir = Path(__file__).resolve().parent
        pipeline_dir = script_dir.parent.parent
        af3_dir = pipeline_dir / DEFAULT_AF3_DIR

    if not af3_dir.is_dir():
        print(f"ERROR: AlphaFold3 results directory not found: {af3_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir = af3_dir / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"AF3 results dir: {af3_dir}")
    print(f"Output dir:      {out_dir}")
    print(f"Background:      {args.background}")
    print()

    # Discover all pairs
    pairs = discover_pairs(af3_dir, args.background)
    if not pairs:
        print("No confidence/interaction image pairs found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(pairs)} complexes:\n")
    for i, p in enumerate(pairs):
        print(f"  {i+1}. {p['gene']}{p['model_str']}  ({p['folder']})")
    print()

    # Generate individual composites
    composites = []
    for i, pair in enumerate(pairs):
        print(f"Processing ({chr(ord('a') + i)}) {pair['gene']}{pair['model_str']}...")
        comp = make_pair(pair, out_dir, i)
        if comp is not None:
            composites.append(comp)

    # Generate combined grid
    print()
    if composites:
        print(f"Building {math.ceil(len(composites)/2)}\u00d72 combined grid...")
        make_grid(composites, pairs, out_dir)
    else:
        print("No composites generated — skipping grid.")

    print(f"\nDone. {len(composites)} composites generated.")


if __name__ == "__main__":
    main()
