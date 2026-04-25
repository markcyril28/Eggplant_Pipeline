#!/bin/bash
set -euo pipefail
#
# Batch processing script for multiple PlantCARE files
# Processes all .tab files in a directory and generates matrices
#

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Help function
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Batch process multiple PlantCARE tab files.

OPTIONS:
    -i, --input-dir DIR      Input directory containing .tab files (default: current)
    -o, --output-dir DIR     Output directory (default: ./plantCARE_matrices)
    -p, --pattern PATTERN    File pattern to match (default: *.tab)
    -v, --visualize          Generate R visualizations (requires R)
    -m, --min-freq NUM       Minimum maximum frequency for CAREs in heatmap (default: 8)
    -h, --help               Show this help message

EXAMPLES:
    # Process all .tab files in current directory
    $(basename "$0")

    # Process files from specific directory
    $(basename "$0") -i raw_data/ -o results/

    # Process with custom pattern and visualization
    $(basename "$0") -p "plantCARE_*.tab" -v

    # Process with lower frequency threshold
    $(basename "$0") -v -m 5

EOF
}

# Default values
INPUT_DIR="."
OUTPUT_DIR="./plantCARE_matrices"
PATTERN="*.tab"
VISUALIZE=false
MIN_FREQ=8

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -v|--visualize)
            VISUALIZE=true
            shift
            ;;
        -m|--min-freq)
            MIN_FREQ="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Error: Unknown option: $1${NC}" >&2
            show_help
            exit 1
            ;;
    esac
done

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo -e "${RED}Error: Input directory '$INPUT_DIR' not found${NC}" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Find all matching files
FILES=("$INPUT_DIR"/$PATTERN)

if [ ! -e "${FILES[0]}" ]; then
    echo -e "${RED}Error: No files matching pattern '$PATTERN' found in '$INPUT_DIR'${NC}" >&2
    exit 1
fi

# Count files
FILE_COUNT=${#FILES[@]}

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Batch PlantCARE Processing${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Input directory: $INPUT_DIR"
echo -e "Output directory: $OUTPUT_DIR"
echo -e "Pattern: $PATTERN"
echo -e "Files found: $FILE_COUNT"
echo -e "Visualization: $VISUALIZE"
echo -e "Min frequency threshold: $MIN_FREQ"
echo ""

# Process each file
SUCCESS_COUNT=0
FAIL_COUNT=0

for file in "${FILES[@]}"; do
    # Get filename without path
    filename=$(basename "$file")
    # Get prefix (remove .tab extension)
    prefix="${filename%.tab}"
    
    echo -e "${YELLOW}Processing: $filename${NC}"
    
    # Run the Python script
    if python3 "$SCRIPT_DIR/plantCARE_to_matrix.py" \
        "$file" \
        -o "$OUTPUT_DIR" \
        -p "$prefix"; then
        
        echo -e "${GREEN}✓ Successfully processed: $filename${NC}"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        
        # Generate visualizations if requested
        if [ "$VISUALIZE" = true ]; then
            if command -v Rscript &> /dev/null; then
                echo -e "${YELLOW}  Generating visualizations...${NC}"
                cd "$OUTPUT_DIR"
                # Note: This uses the Test directory visualization script (ggplot2/pheatmap version)
                # For ComplexHeatmap version with min_freq filtering, use the modules version instead
                if Rscript "$SCRIPT_DIR/visualize_plantCARE_matrix.R" "$prefix"; then
                    echo -e "${GREEN}  ✓ Visualizations created${NC}"
                else
                    echo -e "${YELLOW}  ⚠ Visualization failed (continuing)${NC}"
                fi
                cd - > /dev/null
            else
                echo -e "${YELLOW}  ⚠ Rscript not found, skipping visualization${NC}"
            fi
        fi
        
    else
        echo -e "${RED}✗ Failed to process: $filename${NC}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    echo ""
done

# Final summary
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Batch Processing Complete${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Total files: $FILE_COUNT"
echo -e "${GREEN}Successful: $SUCCESS_COUNT${NC}"
if [ $FAIL_COUNT -gt 0 ]; then
    echo -e "${RED}Failed: $FAIL_COUNT${NC}"
fi
echo -e "Output directory: $OUTPUT_DIR"
echo ""

# List output files
echo -e "${GREEN}Generated files:${NC}"
ls -lh "$OUTPUT_DIR" | tail -n +2

exit 0
