#!/bin/bash
set -euo pipefail

################################################################################
# PyMOL Visualization Pipeline for SmelGRF-GIF Residue Analysis
# Target: fold_1_x_1 (SmelGRF-GIF complex with detailed residue labels)
#
# This script:
#   - Processes fold_1_x_1 with residue-level detail
#   - Uses leader_lines for residue labels
#   - Highlights hydrophobic interface residues (GLN, LEU, VAL, PRO, TYR, ILE in GRF; LEU, ILE in GIF)
#   - Creates both full and helix-only PSE sessions
#   - Generates multi-format outputs in fold_1_x_1/OUTPUT_RESIDUE/
# 
# Usage: 
#   bash a_run_PyMOL_for_SmelGRF-GIF_residue.sh
#   bash a_run_PyMOL_for_SmelGRF-GIF_residue.sh clear-logs
#   bash a_run_PyMOL_for_SmelGRF-GIF_residue.sh session-only
################################################################################

SKIP_EXISTING="false"
CLEAR_LOGS="true"
SESSION_ONLY="false"

if [ $# -gt 0 ]; then
    if [[ "$1" == "clear-logs" ]] || [[ "$1" == "--clear-logs" ]] || [[ "$1" == "-cl" ]]; then
        CLEAR_LOGS="true"
    elif [[ "$1" == "session-only" ]] || [[ "$1" == "--session-only" ]] || [[ "$1" == "-so" ]]; then
        SESSION_ONLY="true"
    fi
fi

echo "============================================================"
echo "PyMOL Visualization - SmelGRF-GIF Residue Analysis"
echo "============================================================"
echo "Target: fold_1_x_1"
echo "Configuration: CLEAR_LOGS=${CLEAR_LOGS}"
echo "Configuration: SESSION_ONLY=${SESSION_ONLY}"
echo ""

if [ "$CLEAR_LOGS" = "true" ]; then
    echo "Clearing log files..."
    if [ -d "logs" ]; then
        rm -f logs/*.log
        echo "[OK] Log files cleared"
    else
        echo "[OK] No logs directory found"
    fi
    echo ""
fi

if [ "$SESSION_ONLY" = "true" ]; then
    echo "Session-only mode: Disabling image rendering..."
    cp modules/pymol_config_gif_grf_residue.py modules/pymol_config_gif_grf_residue.py.backup
    sed -i 's/^RENDER_IMAGES = True/RENDER_IMAGES = False/' modules/pymol_config_gif_grf_residue.py
    echo "[OK] Image rendering disabled"
    echo ""
fi

if ! command -v pymol &> /dev/null; then
    echo "[ERROR] PyMOL not found!"
    exit 1
fi

echo "[OK] PyMOL found"
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VISUALIZE_SCRIPT="${SCRIPT_DIR}/modules/visualize_gif_grf_residue.py"

echo "Running residue-focused visualization..."
cd "${SCRIPT_DIR}" && pymol -c "${VISUALIZE_SCRIPT}"

EXIT_STATUS=$?

if [ "$SESSION_ONLY" = "true" ] && [ -f "modules/pymol_config_gif_grf_residue.py.backup" ]; then
    mv modules/pymol_config_gif_grf_residue.py.backup modules/pymol_config_gif_grf_residue.py
    echo "[OK] Configuration restored"
fi

if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "[OK] SmelGRF-GIF residue visualization completed!"
    echo "============================================================"
else
    echo ""
    echo "============================================================"
    echo "[ERROR] Visualization failed. Check logs for details."
    echo "============================================================"
    exit 1
fi
