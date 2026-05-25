#!/usr/bin/env python3
"""
Composite a PyMOL structure render with a matplotlib legend panel describing
the HAP2/DMP variant band colours. Runs in any conda env with matplotlib;
the band definitions are imported from render_hap2_domain_map.py to keep
them in sync.

Usage:
    python compose_legend.py --structure <structure.png> --out <final.png>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent))
from render_hap2_domain_map import HAP2_BANDS, DMP_BANDS


def compose(structure_png: Path, out_png: Path, dpi: int = 300,
            background: str = "black",
            hap2_title: str = "SmelHAP2 trimer (chains A/B/C) -- deletion-variant bands",
            dmp_title: str = "SmelDMPv5_10.610 (chain D) -- DMP topology palette (08 config)",
            ) -> None:
    img = mpimg.imread(str(structure_png))
    h_in = 8
    img_w_in = h_in * (img.shape[1] / img.shape[0])
    legend_w_in = 6.0
    fig = plt.figure(figsize=(img_w_in + legend_w_in, h_in), dpi=dpi,
                     facecolor=background)
    fig.patch.set_facecolor(background)
    gs = GridSpec(1, 2, width_ratios=[img_w_in, legend_w_in], wspace=0.02)

    text_colour  = "white" if background == "black" else "black"
    patch_edge   = "white" if background == "black" else "black"

    ax_img = fig.add_subplot(gs[0, 0])
    ax_img.imshow(img)
    ax_img.set_axis_off()
    # Force background colour with an explicit full-axes rectangle: set_axis_off
    # hides the axes patch, so figure facecolor would otherwise show through.
    ax_img.add_patch(mpatches.Rectangle(
        (0, 0), 1, 1, transform=ax_img.transAxes,
        facecolor=background, edgecolor="none", zorder=-10,
    ))

    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.set_axis_off()
    ax_leg.add_patch(mpatches.Rectangle(
        (0, 0), 1, 1, transform=ax_leg.transAxes,
        facecolor=background, edgecolor="none", zorder=-10,
    ))

    def patches(bands, title, y_start):
        ax_leg.text(0.0, y_start, title, fontsize=11, fontweight="bold",
                    color=text_colour,
                    transform=ax_leg.transAxes, va="top")
        y = y_start - 0.045
        for start, end, hex_code, label in bands:
            ax_leg.add_patch(mpatches.Rectangle(
                (0.0, y - 0.025), 0.06, 0.025,
                transform=ax_leg.transAxes, facecolor=hex_code,
                edgecolor=patch_edge, linewidth=0.4,
            ))
            # Topology-based rows (e.g. DMP secondary-structure colouring)
            # pass start=None to skip the residue-range prefix.
            text = f"{start}-{end}: {label}" if start is not None else label
            ax_leg.text(0.08, y - 0.012, text,
                        fontsize=8.5, color=text_colour,
                        transform=ax_leg.transAxes, va="center")
            y -= 0.038
        return y - 0.02

    y_after_hap2 = patches(HAP2_BANDS, hap2_title, 0.98)
    patches(DMP_BANDS, dmp_title, y_after_hap2)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_png), dpi=dpi, bbox_inches="tight", facecolor=background)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--background", default="black", choices=["black", "white"])
    ap.add_argument("--hap2-title", default=None, help="override the HAP2 legend section title")
    ap.add_argument("--dmp-title", default=None, help="override the DMP legend section title")
    args = ap.parse_args()
    kwargs = {}
    if args.hap2_title is not None:
        kwargs["hap2_title"] = args.hap2_title
    if args.dmp_title is not None:
        kwargs["dmp_title"] = args.dmp_title
    compose(Path(args.structure), Path(args.out), args.dpi, args.background, **kwargs)
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
