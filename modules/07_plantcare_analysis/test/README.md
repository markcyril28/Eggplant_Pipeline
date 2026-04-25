# PlantCARE Test Suite

Self-contained test scripts for early-stage PlantCARE matrix generation.

**These scripts use their own local copies** of `plantCARE_to_matrix.py` and
`visualize_plantCARE_matrix.R` (ggplot2/pheatmap version). They are **not** the
production versions used by the main pipeline (`run_pipeline.sh`).

## Production scripts (parent directory)

| Script | Purpose |
|--------|---------|
| `run_pipeline.sh` | Main pipeline entry point (orchestrated) |
| `plantCARE_to_matrix.py` | Production matrix generator (v1 motif + v2 function) |
| `visualize_plantCARE_matrix.R` | ComplexHeatmap-based heatmap |
| `post_process_plantcare.py` | Raw PlantCARE post-processor (.tab / .tar.gz) |

## Files in this directory

| File | Notes |
|------|-------|
| `plantCARE_to_matrix.py` | Older standalone 320-line version with more output types |
| `visualize_plantCARE_matrix.R` | ggplot2/pheatmap version (not ComplexHeatmap) |
| `test_pipeline.py` | Quick test importing from local `plantCARE_to_matrix.py` |
| `example_workflow.sh` | End-to-end demo workflow |
| `batch_process_plantCARE.sh` | Batch processing wrapper |
| `run_plantCARE_to_matrix.sh` | Single-file wrapper |
