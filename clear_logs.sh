#!/bin/bash
# clear_logs.sh — Remove all contents of logs/ except z_logs_prompt/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGS_DIR="$SCRIPT_DIR/logs"
KEEP_DIR="z_logs_prompt"

if [[ ! -d "$LOGS_DIR" ]]; then
    echo "ERROR: logs directory not found: $LOGS_DIR"
    exit 1
fi

echo "The following will be deleted from: $LOGS_DIR"
echo "--------------------------------------------"
find "$LOGS_DIR" -maxdepth 1 ! -name "$(basename "$LOGS_DIR")" ! -name "$KEEP_DIR" -print
echo "--------------------------------------------"
echo "Directory preserved: $KEEP_DIR"
echo ""
read -r -p "Proceed? [y/N] " confirm
if [[ "${confirm,,}" != "y" ]]; then
    echo "Aborted."
    exit 0
fi

find "$LOGS_DIR" -maxdepth 1 ! -name "$(basename "$LOGS_DIR")" ! -name "$KEEP_DIR" -exec rm -rf {} +

echo "Done. logs/ cleared (z_logs_prompt preserved)."
