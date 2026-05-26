#!/usr/bin/env python3
"""
Load a saved PyMOL session (.pse) and re-render its current view to PNG.
Preserves the orientation/colours that were saved in the .pse, so the user's
hand-tuned pose is reproduced exactly. Used to refresh the
deletion_variants/model_*/  PNGs after the user re-oriented in PyMOL GUI.

Usage:
    pymol -cq render_from_pse.py -- --pse path/to/session.pse --out path/to/out.png
"""
import argparse
import sys
from pathlib import Path


def main() -> int:
    from pymol import cmd
    ap = argparse.ArgumentParser()
    ap.add_argument("--pse", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ray-width", type=int, default=1800)
    ap.add_argument("--ray-height", type=int, default=1200)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd.reinitialize()
    cmd.load(args.pse)
    cmd.set("ray_opaque_background", 1)
    cmd.ray(args.ray_width, args.ray_height)
    cmd.png(str(out), dpi=args.dpi, ray=0)
    print(f"[OK] wrote {out}")
    return 0


if __name__ == "__main__" or __name__ == "pymol":
    sys.exit(main())
