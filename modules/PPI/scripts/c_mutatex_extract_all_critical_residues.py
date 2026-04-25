#!/usr/bin/env python3
"""
Automatically extract critical residues from all MutateX results
Processes all results in run_results/ directory
"""

import argparse
import os
import sys
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# ==============================================================================
# CONFIGURATION - EDIT THESE PARAMETERS
# ==============================================================================

# Directory containing MutateX result folders
DEFAULT_RUN_RESULTS_DIR = "3_RESULT/DMP-HAP2/11_PPI_MutateX"

# Directory containing input PDB files
DEFAULT_INPUTS_DIR = "3_RESULT/DMP-HAP2/08_Protein_Structure/GPE001970_SMEL5/Protein_Structures"

# Output subdirectory name within each result folder
DEFAULT_OUTPUT_SUBDIR = "critical_residues"

# DDG threshold for destabilizing mutations (kcal/mol)
# Residues with DDG > this value are considered destabilizing
DEFAULT_DESTABILIZING_THRESHOLD = 2.0

# DDG threshold for stabilizing mutations (kcal/mol)
# Residues with DDG < this value are considered stabilizing
DEFAULT_STABILIZING_THRESHOLD = -2.0

# Number of top residues to extract per category
DEFAULT_TOP_N = 30

# ==============================================================================
# FUNCTIONS - DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU'RE DOING
# ==============================================================================

def find_pdb_file(pdb_name, inputs_dir):
    """Find the corresponding PDB file in inputs directory.
    
    Searches both flat layout (inputs_dir/name.pdb) and
    subdirectory layout (inputs_dir/name/name.pdb) used by
    AlphaFold3 structure directories.
    """
    inputs = Path(inputs_dir)
    
    # Flat layout: inputs_dir/name.pdb
    pdb_path = inputs / f"{pdb_name}.pdb"
    if pdb_path.exists():
        return pdb_path
    
    # Subdirectory layout: inputs_dir/name/name.pdb
    pdb_path_sub = inputs / pdb_name / f"{pdb_name}.pdb"
    if pdb_path_sub.exists():
        return pdb_path_sub
    
    # Try with _model0_checked suffix (FoldX repair output)
    pdb_path_checked = inputs / f"{pdb_name}_model0_checked.pdb"
    if pdb_path_checked.exists():
        return pdb_path_checked
    
    # Search recursively as last resort
    matches = list(inputs.rglob(f"{pdb_name}*.pdb"))
    if matches:
        return matches[0]
    
    return None


def _process_one_pdb(result_dir, inputs_path, output_subdir, thresholds):
    """Process a single PDB result directory. Returns (pdb_name, success, message)."""
    pdb_name = result_dir.name

    # Locate the actual data directory: may be result_dir itself
    # or a nruns_N subdirectory created by the orchestrator
    data_dir = result_dir
    results_subdir = data_dir / "results"
    if not results_subdir.exists():
        nruns_dirs = sorted(result_dir.glob("nruns_*/"))
        if nruns_dirs:
            data_dir = nruns_dirs[-1]
            results_subdir = data_dir / "results"

    if not results_subdir.exists():
        return (pdb_name, False, f"No results directory found in {result_dir}")

    # Find corresponding PDB file
    pdb_file = find_pdb_file(pdb_name, inputs_path)
    if not pdb_file:
        return (pdb_name, False, f"PDB file not found for {pdb_name}")

    # Create output directory within the top-level result directory
    output_dir = result_dir / output_subdir

    # Run extraction script (pass data_dir where mutatex wrote results)
    cmd = [
        "python3",
        str(Path(__file__).parent / "extract_critical_residues.py"),
        str(data_dir),
        str(pdb_file),
        "-o", str(output_dir),
        "--destabilizing", str(thresholds['destabilizing']),
        "--stabilizing", str(thresholds['stabilizing']),
        "--top-n", str(thresholds['top_n'])
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return (pdb_name, True, result.stdout)
    except subprocess.CalledProcessError as e:
        return (pdb_name, False, f"Failed: {e.stderr}")
    except Exception as e:
        return (pdb_name, False, f"Unexpected error: {e}")


def process_all_results(run_results_dir, inputs_dir, output_subdir, thresholds, max_workers=1):
    """Process all MutateX results in run_results directory."""
    run_results_path = Path(run_results_dir)
    inputs_path = Path(inputs_dir)

    if not run_results_path.exists():
        print(f"[ERROR] Results directory not found: {run_results_dir}")
        return

    # Find all result directories
    result_dirs = [d for d in run_results_path.iterdir() if d.is_dir()]

    if not result_dirs:
        print(f"[WARN] No result directories found in {run_results_dir}")
        return

    print(f"[INFO] Found {len(result_dirs)} result directories")
    print(f"[INFO] Thresholds: destabilizing >{thresholds['destabilizing']}, "
          f"stabilizing <{thresholds['stabilizing']}, top N={thresholds['top_n']}")
    print("="*80)

    success_count = 0
    failed_count = 0

    effective_workers = min(max_workers, len(result_dirs))
    if effective_workers > 1:
        print(f"[INFO] Processing {len(result_dirs)} PDBs with {effective_workers} parallel workers")

    if effective_workers > 1:
        with ThreadPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                executor.submit(_process_one_pdb, rd, inputs_path, output_subdir, thresholds): rd
                for rd in result_dirs
            }
            for future in as_completed(futures):
                pdb_name, success, msg = future.result()
                if success:
                    print(f"\n[PROCESSING] {pdb_name}")
                    print(msg)
                    success_count += 1
                else:
                    print(f"\n[PROCESSING] {pdb_name}")
                    print(f"  [SKIP/ERROR] {msg}")
                    failed_count += 1
    else:
        for result_dir in result_dirs:
            pdb_name = result_dir.name
            print(f"\n[PROCESSING] {pdb_name}")
            pdb_name, success, msg = _process_one_pdb(
                result_dir, inputs_path, output_subdir, thresholds
            )
            if success:
                print(msg)
                success_count += 1
            else:
                print(f"  [SKIP/ERROR] {msg}")
                failed_count += 1

    print("\n" + "="*80)
    print(f"[SUMMARY] Processed {len(result_dirs)} directories")
    print(f"  Success: {success_count}")
    print(f"  Failed:  {failed_count}")
    print(f"[OUTPUT] Critical residue files saved in each run_results subdirectory under '{output_subdir}/'")


def main():
    parser = argparse.ArgumentParser(
        description="Automatically extract critical residues from all MutateX results"
    )
    parser.add_argument(
        "--run-results",
        default=DEFAULT_RUN_RESULTS_DIR,
        help=f"Directory containing MutateX result folders (default: {DEFAULT_RUN_RESULTS_DIR})"
    )
    parser.add_argument(
        "--inputs",
        default=DEFAULT_INPUTS_DIR,
        help=f"Directory containing PDB files (default: {DEFAULT_INPUTS_DIR})"
    )
    parser.add_argument(
        "-o", "--output",
        default=DEFAULT_OUTPUT_SUBDIR,
        help=f"Output subdirectory name within each result folder (default: {DEFAULT_OUTPUT_SUBDIR})"
    )
    parser.add_argument(
        "--destabilizing",
        type=float,
        default=DEFAULT_DESTABILIZING_THRESHOLD,
        help=f"DDG threshold for destabilizing mutations (default: {DEFAULT_DESTABILIZING_THRESHOLD} kcal/mol)"
    )
    parser.add_argument(
        "--stabilizing",
        type=float,
        default=DEFAULT_STABILIZING_THRESHOLD,
        help=f"DDG threshold for stabilizing mutations (default: {DEFAULT_STABILIZING_THRESHOLD} kcal/mol)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_TOP_N,
        help=f"Number of top residues to extract per category (default: {DEFAULT_TOP_N})"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers for processing PDB results (default: 1)"
    )

    args = parser.parse_args()

    thresholds = {
        'destabilizing': args.destabilizing,
        'stabilizing': args.stabilizing,
        'top_n': args.top_n
    }

    process_all_results(args.run_results, args.inputs, args.output, thresholds,
                        max_workers=args.workers)


if __name__ == "__main__":
    main()
