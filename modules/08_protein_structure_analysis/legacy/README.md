# Module 08 — Protein Structure Analysis

Orchestrator: `h_protein_structure.sh`

## Production Scripts

| Script | Purpose |
|--------|---------|
| `translate.sh` | Translates nucleotide FASTA to protein using EMBOSS `transeq` with parallel chunking |

## `legacy/` — PyMOL Visualization Ecosystem

One-off and semi-automated scripts for rendering AlphaFold3 predictions in PyMOL.
These are not called by the orchestrator; they are run manually.

### Directory Layout

```
legacy/
├── alphafold/              # AlphaFold output handling
│   └── unzip_all.sh        # Batch-unzip AlphaFold result archives
├── create_folders_dmp.sh   # One-off: create per-gene output folders (PLA gene IDs)
├── create_folders_smel_plas.sh  # Same as above (identical gene list)
├── pymol_modules/          # Python config + utilities (imported by visualization scripts)
│   ├── pymol_config_*.py   # Per-interaction-type color/style/sequence configs
│   ├── pymol_utils*.py     # Chain ID, rendering, bubble surface utilities
│   └── visualize_*.py      # Main visualization entry points (run under `pymol -c`)
├── pymol_fold_scripts/     # Fold-specific rendering scripts
│   ├── pymol_script_fold.py        # Canonical generic fold renderer
│   └── visualize_fold_*.py         # Per-fold visualization configs (unique)
└── pymol_runners/          # Bash wrappers that invoke PyMOL with the right scripts
    ├── a_run_PyMOL_for_SmelGRF-GIF_residue.sh
    └── b_run_PyMOL_for_SmelGRF-GIF_and_SmelGIF-SWI2.sh
```

### PyMOL Module Pattern

Each interaction type follows a **config → utils → visualize** pattern:

| Interaction | Config | Utils | Visualize |
|-------------|--------|-------|-----------|
| GRF-GIF | `pymol_config_gif_grf.py` | `pymol_utils.py` | `visualize_gif_grf.py` |
| GRF-GIF residue | `pymol_config_gif_grf_residue.py` | `pymol_utils_residue.py` | `visualize_gif_grf_residue.py` |
| GIF-SWI2 | `pymol_config_gif_swi2.py` | `pymol_utils_gif_swi2.py` | `visualize_gif_swi2.py` |
| Generic protein | `pymol_config_protein.py` | `pymol_utils_protein.py` | `pymol_script_fold.py` |
| Generic residue | `pymol_config_residue.py` | `pymol_utils_residue.py` | — |

### Known Duplication

- `pymol_utils.py` and `pymol_utils_gif_swi2.py` share most functions (hex_to_rgb, chain ID, highlight).
- `visualize_gif_grf.py` and `visualize_gif_swi2.py` differ by ~15 lines (string labels, config import).
- Config files duplicate `CARTOON_CONFIG`, `AMINO_ACID_CODES`, and `SEQUENCES` dicts.
- Future consolidation could extract shared code into a base module.
