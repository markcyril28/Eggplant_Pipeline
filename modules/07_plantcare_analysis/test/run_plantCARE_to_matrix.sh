#!/bin/bash
set -euo pipefail
#
# PlantCARE to Matrix Pipeline Wrapper
# Converts PlantCARE tab files to matrix format
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Help function
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] <input_file>

Convert PlantCARE tab-delimited output to matrix format.

OPTIONS:
    -o, --output-dir DIR    Output directory (default: current directory)
    -p, --prefix PREFIX     Output file prefix (default: plantCARE_matrix)
    -h, --help              Show this help message

EXAMPLES:
    # Process single file
    $(basename "$0") plantCARE_output_PlantCARE_5913.tab

    # Process with custom output directory
    $(basename "$0") -o results/ plantCARE_output_PlantCARE_5913.tab

    # Process with custom prefix
    $(basename "$0") -p gene_analysis plantCARE_output_PlantCARE_5913.tab

    # Process all .tab files in current directory
    for file in *.tab; do
        $(basename "$0") -o matrix_results/ -p "\${file%.tab}" "\$file"
    done

OUTPUT FILES:
    *_count_matrix.tsv           - Count of each motif type per sequence
    *_position_matrix.tsv        - Positions of each motif (comma-separated)
    *_strand_matrix.tsv          - Strand distribution (+/-)
    *_detailed.tsv               - Detailed information per motif
    *_functional_categories.tsv  - Counts by functional category
    *_summary.tsv                - Summary statistics

EOF
}

# Default values
OUTPUT_DIR="."
PREFIX="plantCARE_matrix"
INPUT_FILE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -*)
            echo -e "${RED}Error: Unknown option: $1${NC}" >&2
            echo "Use -h or --help for usage information"
            exit 1
            ;;
        *)
            INPUT_FILE="$1"
            shift
            ;;
    esac
done

# Check if input file provided
if [ -z "$INPUT_FILE" ]; then
    echo -e "${RED}Error: No input file specified${NC}" >&2
    echo "Use -h or --help for usage information"
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file '$INPUT_FILE' not found${NC}" >&2
    exit 1
fi

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if Python script exists
PYTHON_SCRIPT="$SCRIPT_DIR/plantCARE_to_matrix.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo -e "${RED}Error: Python script not found at $PYTHON_SCRIPT${NC}" >&2
    exit 1
fi

# Check if pandas is installed
if ! python3 -c "import pandas" 2>/dev/null; then
    echo -e "${YELLOW}Warning: pandas not installed. Installing...${NC}"
    pip3 install pandas numpy
fi

# Run the Python script
echo -e "${GREEN}Processing PlantCARE file: $INPUT_FILE${NC}"
echo -e "${GREEN}Output directory: $OUTPUT_DIR${NC}"
echo -e "${GREEN}Output prefix: $PREFIX${NC}"
echo ""

python3 "$PYTHON_SCRIPT" "$INPUT_FILE" -o "$OUTPUT_DIR" -p "$PREFIX"

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Pipeline completed successfully!${NC}"
    echo -e "${GREEN}Output files are in: $OUTPUT_DIR/${NC}"
else
    echo -e "${RED}✗ Pipeline failed${NC}" >&2
    exit 1
fi
