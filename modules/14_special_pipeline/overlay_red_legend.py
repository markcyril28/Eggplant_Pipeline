#!/usr/bin/env python3
"""
Overlay a single-row "red = <what>" legend onto a deletion-highlight render
produced by render_deletion_red.py. The legend is drawn as a small swatch
plus a caption in the bottom-left corner of the figure on a black canvas,
so the structure render is not cropped or rescaled.

Usage:
    python overlay_red_legend.py \
        --image hap2_dmp_delEcto_red.png \
        --out   hap2_dmp_delEcto_red_labeled.png \
        --label "Red: delEcto deletion zone (SmelHAP2 22-589 / AtHAP2 66-425) -- HAP2 ectodomain (D1, D2, fusion loop, stem)"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def overlay(image_path: Path, out_path: Path, label: str,
            dpi: int = 300, background: str = "black",
            text_colour: str = "white", swatch: str = "#E60000",
            iptm: float | None = None,
            iptm_category: str | None = None) -> None:
    img = mpimg.imread(str(image_path))
    h, w = img.shape[:2]
    h_in = 8
    w_in = h_in * (w / h)

    fig = plt.figure(figsize=(w_in, h_in), dpi=dpi, facecolor=background)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=(0, 1, 0, 1), aspect="auto")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    # Red swatch + caption in the bottom-left corner, anchored to axes coords.
    # The swatch height is set in axes coords; the width is derived so the
    # swatch is square in display space (axes coords span w_in horizontally
    # and h_in vertically, so sw_w_axes = sw_h_axes * (h_in / w_in)).
    pad_x, pad_y = 0.018, 0.025
    sw_h = 0.030
    sw_w = sw_h * (h_in / w_in)
    ax.add_patch(mpatches.Rectangle(
        (pad_x, pad_y), sw_w, sw_h,
        facecolor=swatch, edgecolor="white", linewidth=0.6,
        transform=ax.transAxes, zorder=10,
    ))
    ax.text(pad_x + sw_w + 0.010, pad_y + sw_h / 2, label,
            fontsize=9, color=text_colour, va="center", ha="left",
            transform=ax.transAxes, zorder=10)

    # Prominent ipTM badge in the top-right corner (the value of ipTM when
    # the red region is deleted from the AF3 model). Colour-coded by
    # category for at-a-glance reading.
    if iptm is not None:
        cat = (iptm_category or "").lower()
        # Baseline-relative classification (see classify_iptm.py):
        #   tolerated         green   (>= 80% of WT baseline)
        #   reduced           amber   (50-80% of WT)
        #   strongly reduced  orange  (25-50% of WT)
        #   catastrophic      dark red (< 25% of WT)
        if "catastrophic" in cat:
            badge_face = "#7a0000"
        elif "strongly reduced" in cat:
            badge_face = "#a85a00"
        elif "reduced" in cat:
            badge_face = "#8a6500"
        elif "tolerated" in cat:
            badge_face = "#1f5f1f"
        else:
            badge_face = "#3a3a3a"
        badge_text = f"ipTM = {iptm:.2f}"
        if iptm_category:
            badge_text += f"\n({iptm_category})"
        ax.text(0.98, 0.97, badge_text,
                fontsize=13, color="white", va="top", ha="right",
                fontweight="bold",
                transform=ax.transAxes, zorder=11,
                bbox=dict(facecolor=badge_face, edgecolor="white",
                          linewidth=0.8, boxstyle="round,pad=0.45"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, facecolor=background,
                bbox_inches=None, pad_inches=0)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--iptm", type=float, default=None,
                    help="optional ipTM badge value (top-right of panel)")
    ap.add_argument("--iptm-category", default=None,
                    help='optional category, e.g. "catastrophic", "tolerated", "reduced but tolerated"')
    args = ap.parse_args()
    overlay(Path(args.image), Path(args.out), args.label,
            dpi=args.dpi, iptm=args.iptm, iptm_category=args.iptm_category)
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
