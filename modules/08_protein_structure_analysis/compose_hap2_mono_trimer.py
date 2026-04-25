#!/usr/bin/env python3
"""
Compose a side-by-side figure of monomeric (left) vs trimeric (right)
HAP2/GCS1 AlphaFold 3 predicted structures, using interaction-domain
coloring only (no pLDDT panel).

Output:
  Figure_6a_hap2_monomer_vs_trimer.jpg / .png
  under AlphaFold3_Results/Combined_HAP2_Mono_Trimer/

CLI args (all optional — defaults match the pipeline layout):
  --background      black | white  (default: black)
  --output-subdir   sub-directory under AlphaFold3_Results/
"""

import argparse
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# ── Source paths (relative to pipeline root) ──────────────────────────────────
HAP2_AF3_DIR = Path(
    "3_RESULT/HAP2/08_Protein_Structure/"
    "GPE001970_SMEL5/AlphaFold3_Results"
)
MONO_SUBDIR  = "fold_hap2_monomeric"
TRIMER_SUBDIR = "fold_hap2_trimeric"

DEFAULT_OUTPUT_SUBDIR = "Combined_HAP2_Mono_Trimer"

# ── Layout constants (matching compose_confidence_interaction.py) ──────────────
BG_COLOR           = (0, 0, 0)
TEXT_COLOR         = (255, 255, 255)
SUBTITLE_COLOR     = (200, 200, 200)
IMAGE_SCALE        = 1.3
TITLE_H            = 120
LEGEND_H           = 120
GAP                = 30
FONT_SIZE          = 44
SUBTITLE_FONT_SIZE = 30
LEGEND_FONT_SIZE   = 28
LEGEND_BOX_SIZE    = 24
LEGEND_SPACING     = 13

# ── Interaction-domain legend (HAP2 palette) ──────────────────────────────────
HAP2_LEGEND = [
    ("Transmembrane helices", (0, 140, 140)),     # deep teal
    ("Other helices",         (128, 89, 194)),    # periwinkle
    ("\u03b2-sheets",         (92, 135, 191)),    # steel blue
    ("Loops",                 (179, 179, 179)),   # gray70
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--background", default="black",
                   help="Background colour of source renders (default: black)")
    p.add_argument("--output-subdir", default=DEFAULT_OUTPUT_SUBDIR,
                   dest="output_subdir",
                   help="Sub-directory under AlphaFold3_Results/ for output")
    return p.parse_args()


def get_font(size):
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
    total_w = 0
    segments = []
    for label, colour in items:
        bbox = draw.textbbox((0, 0), label, font=font)
        tw = bbox[2] - bbox[0]
        seg_w = LEGEND_BOX_SIZE + LEGEND_SPACING + tw
        segments.append((label, colour, seg_w))
        total_w += seg_w
    total_w += LEGEND_SPACING * (len(items) - 1)

    x = x_start + max(0, (max_width - total_w)) // 2
    box_y = y_center - LEGEND_BOX_SIZE // 2

    for label, colour, seg_w in segments:
        draw.rectangle(
            [x, box_y, x + LEGEND_BOX_SIZE, box_y + LEGEND_BOX_SIZE],
            fill=colour, outline=(180, 180, 180),
        )
        draw.text(
            (x + LEGEND_BOX_SIZE + LEGEND_SPACING, box_y - 1),
            label, fill=SUBTITLE_COLOR, font=font,
        )
        x += seg_w + LEGEND_SPACING


def main():
    args = parse_args()

    script_dir   = Path(__file__).resolve().parent
    pipeline_dir = script_dir.parent.parent

    af3_dir  = pipeline_dir / HAP2_AF3_DIR
    out_dir  = af3_dir / args.output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    bg       = args.background

    mono_path   = af3_dir / MONO_SUBDIR  / f"fold_hap2_monomeric_{bg}_interaction.jpg"
    trimer_path = af3_dir / TRIMER_SUBDIR / f"fold_hap2_trimeric_{bg}_interaction.jpg"

    for p in (mono_path, trimer_path):
        if not p.exists():
            print(f"ERROR: file not found: {p}")
            return

    mono_img   = Image.open(mono_path).convert("RGB")
    trimer_img = Image.open(trimer_path).convert("RGB")

    # Scale up for readability
    if IMAGE_SCALE != 1.0:
        sw = int(mono_img.width * IMAGE_SCALE)
        sh = int(mono_img.height * IMAGE_SCALE)
        mono_img = mono_img.resize((sw, sh), Image.LANCZOS)

    w, h = mono_img.size
    if trimer_img.size != (w, h):
        trimer_img = trimer_img.resize((w, h), Image.LANCZOS)

    # Canvas
    canvas_w = w * 2 + GAP
    canvas_h = h + TITLE_H + LEGEND_H
    canvas   = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)
    canvas.paste(mono_img,   (0,       TITLE_H))
    canvas.paste(trimer_img, (w + GAP, TITLE_H))

    draw    = ImageDraw.Draw(canvas)
    font    = get_font(FONT_SIZE)
    subfont = get_font(SUBTITLE_FONT_SIZE)
    legfont = get_font(LEGEND_FONT_SIZE)

    # Main title
    title = "HAP2/GCS1 GPE001970 \u2014 Interaction-Domain Coloring"
    bbox  = draw.textbbox((0, 0), title, font=font)
    draw.text(
        ((canvas_w - (bbox[2] - bbox[0])) / 2, 8),
        title, fill=TEXT_COLOR, font=font,
    )

    # Panel subtitles
    for label, x_off in [
        ("(i)  Monomeric (pre-fusion)", 0),
        ("(ii)  Trimeric post-fusion",  w + GAP),
    ]:
        bb = draw.textbbox((0, 0), label, font=subfont)
        draw.text(
            (x_off + (w - (bb[2] - bb[0])) / 2,
             TITLE_H - SUBTITLE_FONT_SIZE - 10),
            label, fill=SUBTITLE_COLOR, font=subfont,
        )

    # Shared legend centred across full canvas
    legend_y = TITLE_H + h + LEGEND_H // 2
    draw_legend(draw, HAP2_LEGEND, 0, legend_y, canvas_w, legfont)

    out_jpg = out_dir / "Figure_6a_hap2_monomer_vs_trimer.jpg"
    out_png = out_dir / "Figure_6a_hap2_monomer_vs_trimer.png"
    canvas.save(out_jpg, quality=95)
    canvas.save(out_png)
    print(f"Saved: {out_jpg.relative_to(pipeline_dir)}")
    print(f"Saved: {out_png.relative_to(pipeline_dir)}")
    print("Done.")


if __name__ == "__main__":
    main()
