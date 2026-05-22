#!/usr/bin/env python3
"""
Compose side-by-side figures: pLDDT confidence (left) vs interaction domain (right)
for each SmelDMP AlphaFold 3 structure.

Left panel  → confidence render from the canonical (timestamped) run directory.
Right panel → interaction render from the latest run directory (may differ).

Outputs:
  - Individual per-protein composites  (Combined_Confidence_Interaction/)
  - Combined 3×2 grid composite        (Figure_5cd_confidence_vs_interaction_grid.jpg/.png)

CLI args (all optional — defaults match h_protein_structureCONFIG.toml):
  --background          black | white        Background colour of source renders
  --interaction-version interaction | ...    Color version used for right panel
  --output-subdir       Sub-directory name under AlphaFold3_Results/ for output
"""

import argparse
import sys
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ── Result root (relative to pipeline root) ───────────────────────────────────
AF3_DIR = Path(
    "III_RESULT/DMP/08_Protein_Structure/GPE001970_SMEL5/AlphaFold3_Results"
)
DEFAULT_OUTPUT_SUBDIR = "Combined_Confidence_Interaction"

# ── Protein table ─────────────────────────────────────────────────────────────
# confidence_dir / confidence_stem → source for LEFT  (pLDDT confidence render)
# interaction_dir / interaction_stem → source for RIGHT (interaction color render)
# Ordered by descending mean pLDDT (matches manuscript Table 3).
PROTEINS = [
    {
        "confidence_dir":  "2026_04_07_13_18_smel5_04g005390_1",
        "confidence_stem": "2026_04_07_13_18_smel5_04g005390_1",
        "interaction_dir":  "2026_04_07_13_18_smel5_04g005390_1",
        "interaction_stem": "2026_04_07_13_18_smel5_04g005390_1",
        "name": "SmelDMPv5_04",
        "gene": "SMEL5_04g005390.1",
        "panel": "i",
        "tier": "",
    },
    {
        "confidence_dir":  "2026_04_07_13_17_smel5_02g013320_1",
        "confidence_stem": "2026_04_07_13_17_smel5_02g013320_1",
        "interaction_dir":  "2026_04_07_13_17_smel5_02g013320_1",
        "interaction_stem": "2026_04_07_13_17_smel5_02g013320_1",
        "name": "SmelDMPv5_02",
        "gene": "SMEL5_02g013320.1",
        "panel": "ii",
        "tier": "",
    },
    {
        "confidence_dir":  "2026_04_07_13_19_smel5_10g003660_1",
        "confidence_stem": "2026_04_07_13_19_smel5_10g003660_1",
        "interaction_dir":  "smel5_10g003660_1",   # latest re-render
        "interaction_stem": "smel5_10g003660_1",
        "name": "SmelDMPv5_10.660",
        "gene": "SMEL5_10g003660.1",
        "panel": "iii",
        "tier": "",
    },
    {
        "confidence_dir":  "2026_04_07_13_17_smel5_01g026030_1",
        "confidence_stem": "2026_04_07_13_17_smel5_01g026030_1",
        "interaction_dir":  "2026_04_07_13_17_smel5_01g026030_1",
        "interaction_stem": "2026_04_07_13_17_smel5_01g026030_1",
        "name": "SmelDMPv5_01.030",
        "gene": "SMEL5_01g026030.1",
        "panel": "iv",
        "tier": "Tier 1",
    },
    {
        "confidence_dir":  "2026_04_07_13_16_smel5_01g008730_1",
        "confidence_stem": "2026_04_07_13_16_smel5_01g008730_1",
        "interaction_dir":  "smel5_01g008730_1",   # latest re-render
        "interaction_stem": "smel5_01g008730_1",
        "name": "SmelDMPv5_01.730",
        "gene": "SMEL5_01g008730.1",
        "panel": "v",
        "tier": "Tier 2",
    },
    {
        "confidence_dir":  "2026_04_07_13_29_smel5_12g005350_1",
        "confidence_stem": "2026_04_07_13_29_smel5_12g005350_1",
        "interaction_dir":  "smel5_12g005350_1",   # latest re-render
        "interaction_stem": "smel5_12g005350_1",
        "name": "SmelDMPv5_12",
        "gene": "SMEL5_12g005350.1",
        "panel": "vi",
        "tier": "",
    },
    # Primary HI candidate. Placed in row 4 of the grid (centered, spanning both columns)
    # via the in_grid_solo flag; still rendered as a per-gene composite as well.
    {
        "confidence_dir":  "fold_smel5_10g017610_1",
        "confidence_stem": "fold_smel5_10g017610_1",
        "interaction_dir":  "smel5_10g017610_1",   # latest re-render
        "interaction_stem": "smel5_10g017610_1",
        "name": "SmelDMPv5_10.610",
        "gene": "SMEL5_10g017610.1",
        "panel": "vii",
        "tier": "Tier 1, primary",
        "in_grid_solo": True,
    },
]

# ── Layout constants ───────────────────────────────────────────────────────────
BG_COLOR         = (0, 0, 0)        # black canvas
TEXT_COLOR       = (255, 255, 255)  # white title text
SUBTITLE_COLOR   = (200, 200, 200)  # light-gray subtitle text
IMAGE_SCALE      = 1.3              # scale factor applied to source protein renders
TITLE_H          = 120              # pixels reserved for title above each pair
LEGEND_H         = 120              # pixels reserved for legend below each pair
GAP              = 30               # horizontal gap between left/right panels
GRID_GAP         = 40               # gap between rows/cols in combined grid
FONT_SIZE        = 44
SUBTITLE_FONT_SIZE = 30
LEGEND_FONT_SIZE = 28
LEGEND_BOX_SIZE  = 24               # colour swatch side length
LEGEND_SPACING   = 13               # gap between swatch and label / between items

# ── Colour legend definitions (RGB tuples) ─────────────────────────────────────
# pLDDT Confidence (from config/colors_config/protein_structure_colors.toml [confidence])
PLDDT_LEGEND = [
    ("Very high (>90)",  (0, 84, 214)),      # deep blue   [0.00, 0.33, 0.84]
    ("Confident (70–90)", (102, 204, 242)),   # cyan        [0.40, 0.80, 0.95]
    ("Low (50–70)",      (255, 219, 18)),     # yellow      [1.00, 0.86, 0.07]
    ("Very low (<50)",   (255, 125, 69)),     # orange      [1.00, 0.49, 0.27]
]

# Interaction Domain (from [DMP.interaction] in same TOML)
INTERACTION_LEGEND = [
    ("Transmembrane helices",          (212, 48, 5)),    # deep orange-red [0.83, 0.19, 0.02]
    ("N-terminal domain",             (217, 166, 33)),   # golden amber    [0.85, 0.65, 0.13]
    ("Extracellular \u03B2-pleated sheet", (245, 222, 179)),  # beige      [0.96, 0.87, 0.70]
    ("Loops",                         (179, 179, 179)),  # gray70
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--background", default="black",
        help="Background colour of source renders (default: black)"
    )
    p.add_argument(
        "--interaction-version", default="interaction",
        dest="interaction_version",
        help="Color version used for right panel filename suffix (default: interaction)"
    )
    p.add_argument(
        "--output-subdir", default=DEFAULT_OUTPUT_SUBDIR,
        dest="output_subdir",
        help="Sub-directory under AlphaFold3_Results/ for composite output"
    )
    return p.parse_args()


def get_font(size):
    """Try to load a TrueType font, fall back to PIL default."""
    for name in [
        "arial.ttf", "Arial.ttf", "DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "C:/Windows/Fonts/arial.ttf",
    ]:
        try:
            return ImageFont.truetype(name, size)
        except (OSError, IOError):
            continue
    return ImageFont.load_default()


def draw_legend(draw, items, x_start, y_center, max_width, font):
    """Draw a horizontal colour legend centred within *max_width* starting at *x_start*.

    ``items`` is a list of (label, rgb_tuple) pairs.
    Swatches are drawn at *y_center* - LEGEND_BOX_SIZE/2.
    """
    # First pass: measure total width
    total_w = 0
    segments = []
    for label, colour in items:
        bbox = draw.textbbox((0, 0), label, font=font)
        tw = bbox[2] - bbox[0]
        seg_w = LEGEND_BOX_SIZE + LEGEND_SPACING + tw
        segments.append((label, colour, seg_w, tw))
        total_w += seg_w
    total_w += LEGEND_SPACING * (len(items) - 1)  # gaps between segments

    # Centre within the panel
    x = x_start + max(0, (max_width - total_w)) // 2
    box_y = y_center - LEGEND_BOX_SIZE // 2

    for label, colour, seg_w, tw in segments:
        # Swatch
        draw.rectangle(
            [x, box_y, x + LEGEND_BOX_SIZE, box_y + LEGEND_BOX_SIZE],
            fill=colour, outline=(180, 180, 180),
        )
        # Label
        draw.text(
            (x + LEGEND_BOX_SIZE + LEGEND_SPACING, box_y - 1),
            label, fill=SUBTITLE_COLOR, font=font,
        )
        x += seg_w + LEGEND_SPACING


def make_pair(pipeline_dir, protein, out_dir, bg, interaction_version):
    """Create a single side-by-side composite for one protein.

    Left  panel: <confidence_dir>/<confidence_stem>_<bg>_confidence.jpg
    Right panel: <interaction_dir>/<interaction_stem>_<bg>_<interaction_version>.jpg
    """
    conf_path = (
        pipeline_dir / AF3_DIR
        / protein["confidence_dir"]
        / f"{protein['confidence_stem']}_{bg}_confidence.jpg"
    )
    inter_path = (
        pipeline_dir / AF3_DIR
        / protein["interaction_dir"]
        / f"{protein['interaction_stem']}_{bg}_{interaction_version}.jpg"
    )

    if not conf_path.exists():
        print(f"  SKIP (confidence missing): {conf_path}")
        return None
    if not inter_path.exists():
        print(f"  SKIP (interaction missing): {inter_path}")
        return None

    conf_img  = Image.open(conf_path).convert("RGB")
    inter_img = Image.open(inter_path).convert("RGB")

    # Scale up protein images for readability
    if IMAGE_SCALE != 1.0:
        sw = int(conf_img.width * IMAGE_SCALE)
        sh = int(conf_img.height * IMAGE_SCALE)
        conf_img = conf_img.resize((sw, sh), Image.LANCZOS)

    # Resize interaction panel to match confidence dimensions if needed
    w, h = conf_img.size
    if inter_img.size != (w, h):
        inter_img = inter_img.resize((w, h), Image.LANCZOS)

    # Canvas: two panels side-by-side + gap + title bar + legend bar
    canvas_w = w * 2 + GAP
    canvas_h = h + TITLE_H + LEGEND_H
    canvas   = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)

    canvas.paste(conf_img,  (0,        TITLE_H))
    canvas.paste(inter_img, (w + GAP,  TITLE_H))

    draw    = ImageDraw.Draw(canvas)
    font    = get_font(FONT_SIZE)
    subfont = get_font(SUBTITLE_FONT_SIZE)
    legfont = get_font(LEGEND_FONT_SIZE)

    # Main title (centered)
    tier_str = f"  [{protein['tier']}]" if protein["tier"] else ""
    title    = f"({protein['panel']})  {protein['name']} ({protein['gene']}){tier_str}"
    bbox     = draw.textbbox((0, 0), title, font=font)
    draw.text(
        ((canvas_w - (bbox[2] - bbox[0])) / 2, 8),
        title, fill=TEXT_COLOR, font=font
    )

    # Panel subtitles
    conf_label  = "pLDDT Confidence"
    inter_label = "Interaction Domain"
    cb = draw.textbbox((0, 0), conf_label,  font=subfont)
    ib = draw.textbbox((0, 0), inter_label, font=subfont)
    draw.text(
        ((w - (cb[2] - cb[0])) / 2, TITLE_H - SUBTITLE_FONT_SIZE - 8),
        conf_label, fill=SUBTITLE_COLOR, font=subfont
    )
    draw.text(
        (w + GAP + (w - (ib[2] - ib[0])) / 2, TITLE_H - SUBTITLE_FONT_SIZE - 8),
        inter_label, fill=SUBTITLE_COLOR, font=subfont
    )

    # Legends (below panels)
    legend_y = TITLE_H + h + LEGEND_H // 2
    draw_legend(draw, PLDDT_LEGEND,       0,        legend_y, w,  legfont)
    draw_legend(draw, INTERACTION_LEGEND,  w + GAP,  legend_y, w,  legfont)

    out_name = f"{protein['name']}_{protein['gene']}_confidence_vs_interaction.jpg"
    out_path = out_dir / out_name
    canvas.save(out_path, quality=95)
    print(f"  Saved: {out_path.relative_to(pipeline_dir)}")
    return canvas


def make_grid(composites, pipeline_dir, out_dir, solo_composite=None):
    """Arrange paired composites in a 2-column grid; optionally append a single
    solo_composite centered in its own row at the bottom."""
    if not composites:
        return

    cw, ch     = composites[0].size
    cols       = 2
    pair_rows  = (len(composites) + cols - 1) // cols
    extra_rows = 1 if solo_composite is not None else 0
    total_rows = pair_rows + extra_rows
    grid_w     = cols * cw + (cols - 1) * GRID_GAP
    grid_h     = total_rows * ch + (total_rows - 1) * GRID_GAP
    top_h      = 110
    grid       = Image.new("RGB", (grid_w, grid_h + top_h), BG_COLOR)

    draw  = ImageDraw.Draw(grid)
    font  = get_font(FONT_SIZE + 4)
    title = "Figure 5c\u2013d.  AlphaFold 3 SmelDMP Structures: pLDDT Confidence vs Interaction Domain Coloring"
    bbox  = draw.textbbox((0, 0), title, font=font)
    draw.text(((grid_w - (bbox[2] - bbox[0])) / 2, 12), title, fill=TEXT_COLOR, font=font)

    for idx, comp in enumerate(composites):
        row, col = divmod(idx, cols)
        x = col * (cw + GRID_GAP)
        y = top_h + row * (ch + GRID_GAP)
        grid.paste(comp, (x, y))

    if solo_composite is not None:
        x = (grid_w - cw) // 2
        y = top_h + pair_rows * (ch + GRID_GAP)
        grid.paste(solo_composite, (x, y))

    out_jpg = out_dir / "Figure_5cd_confidence_vs_interaction_grid.jpg"
    out_png = out_dir / "Figure_5cd_confidence_vs_interaction_grid.png"
    grid.save(out_jpg, quality=95)
    grid.save(out_png)
    print(f"  Saved grid: {out_jpg.relative_to(pipeline_dir)}")
    print(f"  Saved grid: {out_png.relative_to(pipeline_dir)}")


def main():
    args = parse_args()

    script_dir    = Path(__file__).resolve().parent
    pipeline_dir  = script_dir.parent.parent  # modules/08_.../script → pipeline root

    out_dir = pipeline_dir / AF3_DIR / args.output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Pipeline root:        {pipeline_dir}")
    print(f"Output dir:           {out_dir}")
    print(f"Background:           {args.background}")
    print(f"Interaction version:  {args.interaction_version}")
    print()

    grid_composites = []
    solo_composite  = None
    for p in PROTEINS:
        print(f"Processing {p['name']} ({p['gene']})...")
        comp = make_pair(pipeline_dir, p, out_dir, args.background, args.interaction_version)
        if comp is None:
            continue
        if p.get("in_grid_solo"):
            solo_composite = comp
        elif p.get("in_grid", True):
            grid_composites.append(comp)

    print()
    if len(grid_composites) == 6:
        layout = "3\xd72 + 1 centered solo row" if solo_composite is not None else "3\xd72"
        print(f"Building {layout} combined grid...")
        make_grid(grid_composites, pipeline_dir, out_dir, solo_composite=solo_composite)
    else:
        print(f"Only {len(grid_composites)}/6 paired grid composites generated, skipping grid.")

    print("\nDone.")


if __name__ == "__main__":
    main()
