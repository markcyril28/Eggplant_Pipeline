#!/bin/bash
set -euo pipefail
# Explore PlantCARE data structure

RAW_DIR="PlantCARE_Results"

echo "=== PlantCARE Data Analysis ==="
echo ""

# Count files
FILE_COUNT=$(ls -1 "$RAW_DIR"/*.tab 2>/dev/null | wc -l)
echo "Tab files found: $FILE_COUNT"

# Get all unique major groups
echo ""
echo "Major Groups:"
awk -F'\t' '{print $9}' "$RAW_DIR"/*.tab | grep -v '^$' | sort -u | nl

# Count functions per group
echo ""
echo "Functions by Major Group:"
while IFS= read -r group; do
    if [ -n "$group" ]; then
        echo ""
        echo "$group:"
        awk -F'\t' -v grp="$group" '$9 == grp {print $8}' "$RAW_DIR"/*.tab | grep -v '^$' | sort -u | sed 's/^/  - /'
    fi
done < <(awk -F'\t' '{print $9}' "$RAW_DIR"/*.tab | grep -v '^$' | sort -u)

# Count unique motifs
echo ""
echo "Unique Motifs: $(awk -F'\t' '{print $2}' "$RAW_DIR"/*.tab | sort -u | wc -l)"
echo "Unique Functions: $(awk -F'\t' '{print $8}' "$RAW_DIR"/*.tab | grep -v '^$' | sort -u | wc -l)"
