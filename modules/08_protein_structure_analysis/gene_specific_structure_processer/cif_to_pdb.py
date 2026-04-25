#!/usr/bin/env python3
"""Convert all .cif files in each gene subfolder to .pdb format using gemmi.

Usage (standalone):
    python3 cif_to_pdb.py --input-dir /path/to/Protein_Structures

Usage (orchestrated):
    Called by process_structures.sh with --input-dir flag.
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gemmi


def _convert_one(cif_str: str, base_str: str) -> str:
    """Convert a single CIF to PDB.  Returns status message."""
    cif_path = Path(cif_str)
    pdb_path = cif_path.with_suffix(".pdb")
    try:
        structure = gemmi.read_structure(str(cif_path))
        structure.write_pdb(str(pdb_path))
        return f"Converted: {cif_path.relative_to(base_str)}  →  {pdb_path.name}"
    except Exception as e:
        return f"ERROR: {cif_path.relative_to(base_str)}: {e}"


def convert_cif_to_pdb(input_dir: Path, workers: int = 1) -> int:
    """Convert all .cif files under *input_dir* to .pdb.  Returns count."""
    cif_files = sorted(
        p for p in input_dir.rglob("*.cif") if p.name != "reference.cif"
    )
    if not cif_files:
        return 0

    converted = 0
    base_str = str(input_dir)

    if workers > 1 and len(cif_files) > 1:
        actual_workers = min(workers, len(cif_files))
        print(f"  Converting {len(cif_files)} CIF files with {actual_workers} parallel workers")
        with ProcessPoolExecutor(max_workers=actual_workers) as pool:
            futures = {
                pool.submit(_convert_one, str(p), base_str): p
                for p in cif_files
            }
            for future in as_completed(futures):
                msg = future.result()
                print(f"  {msg}")
                if not msg.startswith("ERROR"):
                    converted += 1
    else:
        for cif_path in cif_files:
            msg = _convert_one(str(cif_path), base_str)
            print(f"  {msg}")
            if not msg.startswith("ERROR"):
                converted += 1
    return converted


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert CIF structure files to PDB format via gemmi."
    )
    parser.add_argument(
        "--input-dir", required=True, type=Path,
        help="Root directory containing gene sub-folders with .cif files."
    )
    parser.add_argument(
        "--workers", type=int,
        default=max(1, os.cpu_count() or 1),
        help="Number of parallel worker processes (default: CPU count)."
    )
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        print(f"ERROR: directory does not exist: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    n = convert_cif_to_pdb(args.input_dir, workers=args.workers)
    if n == 0:
        print("No CIF files found — nothing to convert.")
    else:
        print(f"\nConverted {n} CIF file(s).")


if __name__ == "__main__":
    main()
