#!/usr/bin/env python3
"""
Crop an image to its tight non-black bounding box. Useful for trimming the
generous black margins PyMOL leaves around ray-traced renders before
stitching the renders into a multi-panel figure.

Usage:
    python crop_to_content.py --image in.png --out out.png [--pad 12] [--threshold 8]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def crop(image: Path, out: Path, pad: int = 12, threshold: int = 8) -> None:
    img = Image.open(str(image))
    arr = np.asarray(img.convert("RGB"))
    # Non-black = any channel above threshold (0-255).
    mask = arr.max(axis=2) > threshold
    if not mask.any():
        # Nothing to crop; just copy through.
        img.save(str(out))
        print(f"[WARN] {image.name} is entirely below threshold; copied as-is")
        return
    ys, xs = np.where(mask)
    y0, y1 = ys.min(), ys.max() + 1
    x0, x1 = xs.min(), xs.max() + 1
    # Apply pad inside the original frame.
    h, w = arr.shape[:2]
    y0 = max(0, y0 - pad)
    y1 = min(h, y1 + pad)
    x0 = max(0, x0 - pad)
    x1 = min(w, x1 + pad)
    cropped = img.crop((x0, y0, x1, y1))
    out.parent.mkdir(parents=True, exist_ok=True)
    cropped.save(str(out))
    print(f"[OK] cropped {image.name}: ({w}x{h}) -> ({x1-x0}x{y1-y0}) -> {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=12,
                    help="pixels of black padding kept around the content bbox")
    ap.add_argument("--threshold", type=int, default=8,
                    help="max RGB value (0-255) treated as 'black' background")
    args = ap.parse_args()
    crop(Path(args.image), Path(args.out), pad=args.pad, threshold=args.threshold)
    return 0


if __name__ == "__main__":
    sys.exit(main())
