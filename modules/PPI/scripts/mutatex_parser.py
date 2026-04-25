#!/usr/bin/env python3
"""Shared parser for MutateX selfmutation .dat files.

Used by both b_mutatex_generate_mutatex_visualizations.py and
extract_critical_residues.py to avoid code duplication.
"""

import os
import pandas as pd


def parse_selfmutation_dat_file(dat_file_path, data_type='interface'):
    """Parse a selfmutation .dat file (format: residue_id avg std min max).

    Args:
        dat_file_path: Path to the .dat file
        data_type: 'interface' or 'folding'

    Returns:
        DataFrame with columns: original_aa, position, position_label, chain,
        ddg, ddg_std, ddg_min, ddg_max, data_type, and optionally chain_pair.
    """
    data = []

    try:
        with open(dat_file_path, 'r') as f:
            lines = f.readlines()

        # Extract chain pair from filename if interface data
        chain_pair = None
        if 'selfmutation_energies_' in os.path.basename(dat_file_path):
            chain_pair = os.path.basename(dat_file_path).replace('selfmutation_energies_', '').replace('.dat', '')

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            residue_id = parts[0]
            avg_ddg = float(parts[1])
            std_ddg = float(parts[2])
            min_ddg = float(parts[3])
            max_ddg = float(parts[4])

            # Parse residue_id: e.g., "MA1" -> M (aa), A (chain), 1 (position)
            if len(residue_id) < 3:
                continue

            original_aa = residue_id[0]
            chain_id = residue_id[1]
            position_str = residue_id[2:]

            try:
                position = int(position_str)
            except ValueError:
                continue

            record = {
                "residue_id": residue_id,
                "original_aa": original_aa,
                "position": position,
                "position_label": f"{original_aa}{chain_id}{position}",
                "chain": chain_id,
                "ddg": avg_ddg,
                "ddg_std": std_ddg,
                "ddg_min": min_ddg,
                "ddg_max": max_ddg,
                "data_type": data_type,
            }

            if chain_pair:
                record["chain_pair"] = chain_pair

            data.append(record)

    except (IOError, ValueError) as e:
        print(f"[WARNING] Error parsing {dat_file_path}: {e}")

    return pd.DataFrame(data)
