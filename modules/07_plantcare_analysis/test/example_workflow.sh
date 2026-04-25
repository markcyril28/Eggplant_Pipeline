#!/bin/bash
set -euo pipefail
#
# Complete example workflow for PlantCARE matrix pipeline
# This script demonstrates all features of the pipeline
#

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  PlantCARE to Matrix Pipeline - Complete Example          ║${NC}"
echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo -e "${GREEN}Step 1: Verify input file${NC}"
echo "--------------------------------------------------------------"
INPUT_FILE="plantCARE_output_PlantCARE_5913.tab"

if [ -f "$INPUT_FILE" ]; then
    echo "✓ Input file found: $INPUT_FILE"
    LINE_COUNT=$(wc -l < "$INPUT_FILE")
    echo "  Lines: $LINE_COUNT"
else
    echo "✗ Input file not found: $INPUT_FILE"
    echo "  Please ensure the PlantCARE output file exists"
    exit 1
fi
echo ""

echo -e "${GREEN}Step 2: Check dependencies${NC}"
echo "--------------------------------------------------------------"

# Check Python
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    echo "✓ Python: $PYTHON_VERSION"
else
    echo "✗ Python 3 not found"
    exit 1
fi

# Check pandas
if python3 -c "import pandas" 2>/dev/null; then
    PANDAS_VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
    echo "✓ pandas: $PANDAS_VERSION"
else
    echo "⚠ pandas not installed"
    echo "  Installing pandas..."
    pip3 install pandas numpy
fi

# Check R (optional)
if command -v Rscript &> /dev/null; then
    R_VERSION=$(Rscript --version 2>&1 | head -1)
    echo "✓ R: $R_VERSION"
    HAS_R=true
else
    echo "⚠ R not found (visualization will be skipped)"
    HAS_R=false
fi
echo ""

echo -e "${GREEN}Step 3: Test the pipeline${NC}"
echo "--------------------------------------------------------------"
echo "Running quick test..."
python3 test_pipeline.py
echo ""

echo -e "${GREEN}Step 4: Process single file${NC}"
echo "--------------------------------------------------------------"
OUTPUT_DIR="example_output"
mkdir -p "$OUTPUT_DIR"

echo "Processing: $INPUT_FILE"
python3 plantCARE_to_matrix.py "$INPUT_FILE" \
    -o "$OUTPUT_DIR" \
    -p "example_gene"

if [ $? -eq 0 ]; then
    echo "✓ Processing completed successfully"
else
    echo "✗ Processing failed"
    exit 1
fi
echo ""

echo -e "${GREEN}Step 5: Verify output files${NC}"
echo "--------------------------------------------------------------"
EXPECTED_FILES=(
    "example_gene_count_matrix.tsv"
    "example_gene_position_matrix.tsv"
    "example_gene_strand_matrix.tsv"
    "example_gene_detailed.tsv"
    "example_gene_functional_categories.tsv"
    "example_gene_summary.tsv"
)

ALL_FOUND=true
for file in "${EXPECTED_FILES[@]}"; do
    if [ -f "$OUTPUT_DIR/$file" ]; then
        SIZE=$(du -h "$OUTPUT_DIR/$file" | cut -f1)
        echo "✓ $file ($SIZE)"
    else
        echo "✗ $file (NOT FOUND)"
        ALL_FOUND=false
    fi
done
echo ""

if [ "$ALL_FOUND" = false ]; then
    echo "⚠ Some output files are missing"
    exit 1
fi

echo -e "${GREEN}Step 6: Display sample results${NC}"
echo "--------------------------------------------------------------"
echo -e "${YELLOW}Summary Statistics:${NC}"
cat "$OUTPUT_DIR/example_gene_summary.tsv"
echo ""

echo -e "${YELLOW}Top 5 rows from Count Matrix:${NC}"
head -6 "$OUTPUT_DIR/example_gene_count_matrix.tsv" | column -t -s $'\t'
echo ""

echo -e "${YELLOW}Functional Categories:${NC}"
cat "$OUTPUT_DIR/example_gene_functional_categories.tsv" | column -t -s $'\t'
echo ""

if [ "$HAS_R" = true ]; then
    echo -e "${GREEN}Step 7: Generate visualizations${NC}"
    echo "--------------------------------------------------------------"
    cd "$OUTPUT_DIR"
    
    # Note: This uses the Test directory visualization script (ggplot2/pheatmap version)
    # For ComplexHeatmap version with min_freq filtering, use:
    # Rscript ../modules/visualize_plantCARE_matrix.R -i example_gene_count_matrix.tsv -o heatmap.png --min_freq 8
    if Rscript ../visualize_plantCARE_matrix.R "example_gene" 2>/dev/null; then
        echo "✓ Visualizations created successfully"
        
        # List visualization files
        echo ""
        echo "Generated plots:"
        ls -lh example_gene_plots*.pdf 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
    else
        echo "⚠ Visualization generation encountered issues"
    fi
    cd ..
    echo ""
fi

echo -e "${GREEN}Step 8: Create analysis report${NC}"
echo "--------------------------------------------------------------"

REPORT_FILE="$OUTPUT_DIR/analysis_report.txt"

cat > "$REPORT_FILE" << EOF
================================================================================
PlantCARE Cis-Regulatory Element Analysis Report
================================================================================

Input File: $INPUT_FILE
Analysis Date: $(date)
Output Directory: $OUTPUT_DIR

================================================================================
SUMMARY STATISTICS
================================================================================

$(cat "$OUTPUT_DIR/example_gene_summary.tsv")

================================================================================
FUNCTIONAL CATEGORY DISTRIBUTION
================================================================================

$(cat "$OUTPUT_DIR/example_gene_functional_categories.tsv")

================================================================================
TOP 10 MOST FREQUENT MOTIFS
================================================================================

$(awk -F'\t' 'NR>1 {motif[$2]++} END {for (m in motif) print motif[m], m}' \
   "$OUTPUT_DIR/example_gene_detailed.tsv" | sort -rn | head -10 | \
   awk '{printf "%-40s %5d\n", $2, $1}')

================================================================================
OUTPUT FILES
================================================================================

Count Matrix:              example_gene_count_matrix.tsv
Position Matrix:           example_gene_position_matrix.tsv
Strand Matrix:             example_gene_strand_matrix.tsv
Detailed Matrix:           example_gene_detailed.tsv
Functional Categories:     example_gene_functional_categories.tsv
Summary Statistics:        example_gene_summary.tsv

$(if [ "$HAS_R" = true ]; then
    echo "Visualizations:"
    ls example_gene_plots*.pdf 2>/dev/null | sed 's/^/    /'
fi)

================================================================================
NOTES
================================================================================

- All matrices are tab-delimited and can be opened in Excel, R, or Python
- Position matrix contains comma-separated genomic positions
- Strand matrix shows distribution as +:count/-:count
- Functional categories are automatically assigned based on motif descriptions

For detailed usage information, see:
- README_PlantCARE_Pipeline.md
- QUICK_START.md

================================================================================
EOF

echo "✓ Analysis report created: $REPORT_FILE"
echo ""

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  Pipeline Execution Complete!                              ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}Summary:${NC}"
echo "  ✓ Input file processed successfully"
echo "  ✓ 6 matrix files generated"
if [ "$HAS_R" = true ]; then
    echo "  ✓ Visualizations created"
fi
echo "  ✓ Analysis report generated"
echo ""
echo -e "${YELLOW}Output Location:${NC}"
echo "  $OUTPUT_DIR/"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Review the analysis report: cat $REPORT_FILE"
echo "  2. Open matrices in your preferred tool (Excel, R, Python)"
if [ "$HAS_R" = true ]; then
    echo "  3. View the PDF plots for visual insights"
fi
echo "  4. Use the matrices for downstream analysis"
echo ""
echo -e "${GREEN}For more examples, see QUICK_START.md${NC}"
echo ""
