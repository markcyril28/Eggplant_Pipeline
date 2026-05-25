#!/usr/bin/env python3
"""
Combine several deletion-highlight renders into one figure with subpanel
labels (A, B, C, ...) along the top of each tile. Designed to stitch the
delEcto / delEctoD2 / delEctoAndC red-highlight PNGs into a single
Figure 3 image. Panels are laid out left-to-right on a black canvas.

Usage:
    python combine_panels.py \
        --image hap2_dmp_delEcto_red_labeled.png \
        --image hap2_dmp_delEctoD2_red_labeled.png \
        --image hap2_dmp_delEctoAndC_red_labeled.png \
        --out   hap2_dmp_catastrophic_panel.png
"""

from __future__ import annotations

import argparse
import string
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def combine(images: list[Path], out: Path, dpi: int = 300,
            background: str = "black", text_colour: str = "white",
            captions: list[str] | None = None,
            panel_height_in: float = 6.0,
            cols: int = 0) -> None:
    """Combine images in a grid of `cols` columns. cols=0 means 1xN single row."""
    if not images:
        raise ValueError("no images supplied")
    arrs = [mpimg.imread(str(p)) for p in images]
    if cols <= 0:
        cols = len(arrs)
    rows = (len(arrs) + cols - 1) // cols
    # Same height per panel; widths follow aspect. For a clean grid, use the
    # max width across the row's panels for every column-cell.
    aspects = [a.shape[1] / a.shape[0] for a in arrs]
    panel_widths = [panel_height_in * r for r in aspects]
    # Per-column width = max of the panels falling in that column.
    col_widths = [0.0] * cols
    for i, w in enumerate(panel_widths):
        c = i % cols
        col_widths[c] = max(col_widths[c], w)
    fig_w = sum(col_widths)
    fig_h = panel_height_in * rows

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi, facecolor=background)
    gs = GridSpec(rows, cols, figure=fig, width_ratios=col_widths,
                  wspace=0.01, hspace=0.02,
                  left=0, right=1, top=1 - 0.04 / rows, bottom=0)

    for i, (img, p) in enumerate(zip(arrs, images)):
        r, c = divmod(i, cols)
        ax = fig.add_subplot(gs[r, c])
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_facecolor(background)
        letter = string.ascii_uppercase[i]
        cap = f"({letter})"
        if captions and i < len(captions) and captions[i]:
            cap = f"({letter}) {captions[i]}"
        # Subpanel label, top-left of each tile.
        ax.text(0.02, 0.98, cap, transform=ax.transAxes,
                fontsize=14, fontweight="bold", color=text_colour,
                va="top", ha="left",
                bbox=dict(facecolor=background, edgecolor="none",
                          alpha=0.55, pad=3.0))

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=dpi, facecolor=background,
                bbox_inches=None, pad_inches=0.05)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", action="append", required=True,
                    help="Panel image (repeat for each panel; order = A, B, C, ...).")
    ap.add_argument("--caption", action="append", default=[],
                    help="Optional caption per panel, in the same order as --image.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--panel-height", type=float, default=6.0)
    ap.add_argument("--cols", type=int, default=0,
                    help="grid columns (default 0 = 1xN single row).")
    args = ap.parse_args()
    combine([Path(p) for p in args.image], Path(args.out),
            dpi=args.dpi, captions=args.caption,
            panel_height_in=args.panel_height, cols=args.cols)
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
