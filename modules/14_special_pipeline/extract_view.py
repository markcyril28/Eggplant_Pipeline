#!/usr/bin/env python3
"""
Extract the rotation matrix from a saved PyMOL session (.pse) and write it
to a JSON file so other renders can re-apply the same orientation.

Usage:
    pymol -cq extract_view.py -- --pse path/to/session.pse --out path/to/view.json
"""
import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    from pymol import cmd
    ap = argparse.ArgumentParser()
    ap.add_argument("--pse", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    cmd.load(args.pse)
    view = list(cmd.get_view())
    rotation = view[:9]  # first 9 floats = 3x3 rotation, row-major
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"rotation": rotation}, indent=2), encoding="utf-8")
    print(f"[OK] wrote view rotation to {out}")
    print("rotation =", rotation)
    return 0


if __name__ == "__main__" or __name__ == "pymol":
    sys.exit(main())
