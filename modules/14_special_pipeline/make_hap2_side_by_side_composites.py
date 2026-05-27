#!/usr/bin/env python3
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

BASE = Path(__file__).resolve().parents[2] / "III_RESULT" / "DMP_x_SmelHAP2" / "14_Domain_Mapping" / "deletion_ladder" / "07_Summary" / "deletion_variants" / "HAP2_variants"

PAIRS = [
    (
        BASE / "hap2_dmp_ectodomain_deletions_grouped_monomeric_prefusion_single.png",
        BASE / "hap2_dmp_ectodomain_deletions_grouped_trimeric_postfusion_single.png",
        BASE / "hap2_single_domain_deletions_side_by_side.png",
    ),
    (
        BASE / "hap2_dmp_ectodomain_deletions_grouped_monomeric_prefusion_combined.png",
        BASE / "hap2_dmp_ectodomain_deletions_grouped_trimeric_postfusion_combined.png",
        BASE / "hap2_combined_domain_deletions_side_by_side.png",
    ),
]

GAP = 0
BG = (0, 0, 0, 255)
FG = (255, 255, 255, 255)
LEFT_LABEL = "Monomeric (1:1 HAP2-DMP)"
RIGHT_LABEL = "Trimeric (3:1 HAP2-DMP, postfusion-like)"

def load_font(size):
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()

for left_path, right_path, out_path in PAIRS:
    left = Image.open(left_path).convert("RGBA")
    right = Image.open(right_path).convert("RGBA")
    h = max(left.height, right.height)

    def pad(img):
        if img.height == h:
            return img
        canvas = Image.new("RGBA", (img.width, h), BG)
        canvas.paste(img, (0, (h - img.height) // 2), img)
        return canvas

    left = pad(left)
    right = pad(right)

    w = left.width + GAP + right.width
    font_size = max(28, w // 60)
    font = load_font(font_size)
    label_band = int(font_size * 2.2)

    out = Image.new("RGBA", (w, h + label_band), BG)
    out.paste(left, (0, label_band), left)
    out.paste(right, (left.width + GAP, label_band), right)

    draw = ImageDraw.Draw(out)
    def center_text(text, x0, x1, y_band):
        bbox = draw.textbbox((0, 0), text, font=font)
        tw = bbox[2] - bbox[0]
        th = bbox[3] - bbox[1]
        x = x0 + ((x1 - x0) - tw) // 2 - bbox[0]
        y = (y_band - th) // 2 - bbox[1]
        draw.text((x, y), text, font=font, fill=FG)

    center_text(LEFT_LABEL, 0, left.width, label_band)
    center_text(RIGHT_LABEL, left.width + GAP, w, label_band)

    out.save(out_path, format="PNG", optimize=True)
    print(f"wrote {out_path.name}: {out.size}")
