#!/usr/bin/env python3
"""Three-panel composite of monomeric (left), dimeric (middle), and trimeric
(right) HAP2/GCS1 AlphaFold 3 structures with HAP2 [interaction] coloring.

Extends compose_hap2_mono_trimer.py to a third panel covering the AF3 dimer
prediction stored under HAP2_Monomeric_to_Quatermeric/. Outputs go under
III_RESULT/DMP-HAP2/.../HAP2_Monomeric_to_Quatermeric/Combined_HAP2_Mono_Dimer_Trimer/
so the composite lives alongside the source folder containing the dimer +
quaternary AF3 predictions.

CLI args (all optional):
  --background      black | white  (default: black)
  --output-subdir   sub-directory name for the composite output (default:
                    Combined_HAP2_Mono_Dimer_Trimer)
"""

import argparse
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# Source render layout (relative to pipeline root) - same monomer/trimer JPGs
# used by compose_hap2_mono_trimer.py; the dimer JPG is rendered into a new
# sibling folder under the same AlphaFold3_Results/ tree.
HAP2_AF3_DIR = Path(
    "III_RESULT/HAP2/08_Protein_Structure/"
    "GPE001970_SMEL5/AlphaFold3_Results"
)
MONO_SUBDIR   = "fold_hap2_monomeric"
DIMER_SUBDIR  = "fold_hap2_dimeric"
TRIMER_SUBDIR = "fold_hap2_trimeric"

# Composite output goes under the DMP-HAP2 tree (source of the dimer/quaternary AF3 runs)
COMPOSITE_PARENT = Path(
    "III_RESULT/DMP-HAP2/08_Protein_Structure/"
    "GPE001970_SMEL5/HAP2_Monomeric_to_Quatermeric"
)
DEFAULT_OUTPUT_SUBDIR = "Combined_HAP2_Mono_Dimer_Trimer"
OUTPUT_STEM = "hap2_monomer_dimer_trimer"

# Layout constants (matched to compose_hap2_mono_trimer.py)
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

HAP2_LEGEND = [
    ("Transmembrane helices", (0, 140, 140)),
    ("Other helices",         (128, 89, 194)),
    ("β-sheets",         (92, 135, 191)),
    ("Loops",                 (179, 179, 179)),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--background", default="black",
                   help="Background colour of source renders (default: black)")
    p.add_argument("--output-subdir", default=DEFAULT_OUTPUT_SUBDIR,
                   dest="output_subdir",
                   help="Sub-directory name under HAP2_Monomeric_to_Quatermeric/")
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

    af3_dir   = pipeline_dir / HAP2_AF3_DIR
    out_dir   = pipeline_dir / COMPOSITE_PARENT / args.output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    bg        = args.background

    panel_paths = [
        ("(i)  Monomeric (pre-fusion)",
         af3_dir / MONO_SUBDIR  / f"fold_hap2_monomeric_{bg}_interaction.jpg"),
        ("(ii)  Dimeric assembly",
         af3_dir / DIMER_SUBDIR / f"fold_hap2_dimeric_{bg}_interaction.jpg"),
        ("(iii)  Trimeric post-fusion",
         af3_dir / TRIMER_SUBDIR / f"fold_hap2_trimeric_{bg}_interaction.jpg"),
    ]

    for _, p in panel_paths:
        if not p.exists():
            print(f"ERROR: panel image not found: {p}")
            return

    imgs = [Image.open(p).convert("RGB") for _, p in panel_paths]

    if IMAGE_SCALE != 1.0:
        sw = int(imgs[0].width  * IMAGE_SCALE)
        sh = int(imgs[0].height * IMAGE_SCALE)
        imgs[0] = imgs[0].resize((sw, sh), Image.LANCZOS)

    w, h = imgs[0].size
    for i in range(1, len(imgs)):
        if imgs[i].size != (w, h):
            imgs[i] = imgs[i].resize((w, h), Image.LANCZOS)

    canvas_w = w * 3 + GAP * 2
    canvas_h = h + TITLE_H + LEGEND_H
    canvas   = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)
    for i, img in enumerate(imgs):
        canvas.paste(img, (i * (w + GAP), TITLE_H))

    draw    = ImageDraw.Draw(canvas)
    font    = get_font(FONT_SIZE)
    subfont = get_font(SUBTITLE_FONT_SIZE)
    legfont = get_font(LEGEND_FONT_SIZE)

    title = "HAP2/GCS1 GPE001970 — Interaction-Domain Coloring"
    bbox  = draw.textbbox((0, 0), title, font=font)
    draw.text(
        ((canvas_w - (bbox[2] - bbox[0])) / 2, 8),
        title, fill=TEXT_COLOR, font=font,
    )

    for i, (label, _) in enumerate(panel_paths):
        x_off = i * (w + GAP)
        bb = draw.textbbox((0, 0), label, font=subfont)
        draw.text(
            (x_off + (w - (bb[2] - bb[0])) / 2,
             TITLE_H - SUBTITLE_FONT_SIZE - 10),
            label, fill=SUBTITLE_COLOR, font=subfont,
        )

    legend_y = TITLE_H + h + LEGEND_H // 2
    draw_legend(draw, HAP2_LEGEND, 0, legend_y, canvas_w, legfont)

    out_jpg = out_dir / f"{OUTPUT_STEM}.jpg"
    out_png = out_dir / f"{OUTPUT_STEM}.png"
    canvas.save(out_jpg, quality=95)
    canvas.save(out_png)
    print(f"Saved: {out_jpg.relative_to(pipeline_dir)}")
    print(f"Saved: {out_png.relative_to(pipeline_dir)}")
    print("Done.")


if __name__ == "__main__":
    main()
