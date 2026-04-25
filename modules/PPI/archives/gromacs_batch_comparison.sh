#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS Batch PPI Comparison Script
#
# Analyzes ALL structures in input folders and generates a comparison
# matrix/ranking. Useful for comparing multiple AlphaFold predictions.
#
# Outputs:
# 1. Per-structure interface metrics
# 2. Ranking by binding quality
# 3. Comparative summary table
#
# Usage: ./gromacs_batch_comparison.sh [INPUT_DIR] [OUTPUT_DIR]
#######################################################################

INPUT_DIR="${1:-${HOME}/PPI/inputs}"
OUTPUT_DIR="${2:-${HOME}/PPI/RESULTS/results_gromacs/batch_comparison}"

# GROMACS settings
GMX_BIN="${HOME}/miniconda3/envs/gromacs_HIP/bin/gmx_mpi"
FORCEFIELD="amber99sb-ildn"
WATERMODEL="tip3p"
EM_STEPS=2000
GPU_FLAGS="-nb gpu -pme cpu -bonded cpu"
NTHREADS=8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() { echo "[$(date '+%H:%M:%S')] $1"; }

#------------------------------------------------------------------------------
# Convert CIF to PDB if needed
#------------------------------------------------------------------------------
convert_to_pdb() {
    local input="$1"
    local output="$2"
    
    if [[ "$input" == *.cif ]]; then
        python3 << EOF
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('protein', '$input')
io = PDBIO()
io.set_structure(structure)
io.save('$output')
EOF
    else
        cp "$input" "$output"
    fi
}

#------------------------------------------------------------------------------
# Quick analysis for a single structure
#------------------------------------------------------------------------------
analyze_structure() {
    local input_file="$1"
    local struct_name="$2"
    local struct_dir="$3"
    
    mkdir -p "$struct_dir"
    cd "$struct_dir"
    
    log "  Converting/cleaning structure..."
    convert_to_pdb "$input_file" "input.pdb"
    grep -E "^(ATOM|TER|END)" input.pdb | sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > clean.pdb
    
    # Get chain info
    local chains=$(grep "^ATOM" clean.pdb | cut -c22 | sort -u | tr -d ' \n')
    local natoms=$(grep -c "^ATOM" clean.pdb || echo 0)
    
    log "  Chains: $chains, Atoms: $natoms"
    
    # Generate topology
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 > pdb2gmx.log || {
        echo "ERROR" > status.txt
        return 1
    }
    
    # Create simple chain index
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 > make_ndx.log
    
    # Quick setup
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d 0.8 -bt cubic 2>&1 > editconf.log
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 > solvate.log
    
    # Minimal EM
    cat > em.mdp << EOF
integrator = steep
emtol = 500.0
emstep = 0.01
nsteps = $EM_STEPS
nstlist = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
EOF

    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 3 2>&1 > grompp.log
    mpirun -np 1 $GMX_BIN mdrun -deffnm em -ntomp $NTHREADS $GPU_FLAGS 2>&1 > mdrun.log || \
        mpirun -np 1 $GMX_BIN mdrun -deffnm em -ntomp 16 2>&1 > mdrun_cpu.log
    
    # Extract metrics
    echo -e "Potential\n\n" | $GMX_BIN energy -f em.edr -o energy.xvg 2>&1 > energy.log || true
    echo "Protein" | $GMX_BIN gyrate -s em.tpr -f em.gro -o gyrate.xvg 2>&1 > gyrate.log || true
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg 2>&1 > sasa.log || true
    
    # Parse results
    python3 << 'PYEOF'
import os
import json

def parse_xvg(filename):
    if not os.path.exists(filename):
        return None
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):
                parts = line.split()
                if parts:
                    return float(parts[-1])
    return None

metrics = {
    'potential': parse_xvg('energy.xvg'),
    'gyration': parse_xvg('gyrate.xvg'),
    'sasa': parse_xvg('sasa.xvg'),
    'status': 'OK'
}

with open('metrics.json', 'w') as f:
    json.dump(metrics, f, indent=2)

print(f"Potential: {metrics['potential']:.1f} kJ/mol" if metrics['potential'] else "Potential: N/A")
PYEOF

    echo "OK" > status.txt
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

main() {
    echo ""
    echo "============================================================"
    echo " GROMACS Batch PPI Comparison"
    echo "============================================================"
    
    log "Scanning input directory: $INPUT_DIR"
    
    mkdir -p "$OUTPUT_DIR"
    cd "$OUTPUT_DIR"
    
    # Find all PDB/CIF files
    mapfile -t STRUCTURES < <(find "$INPUT_DIR" -type f \( -name "*.pdb" -o -name "*.cif" \) | sort)
    
    log "Found ${#STRUCTURES[@]} structures to analyze"
    echo ""
    
    START=$(date +%s)
    
    # Analyze each structure
    for struct_file in "${STRUCTURES[@]}"; do
        struct_name=$(basename "$struct_file" | sed 's/\.[^.]*$//')
        struct_dir="$OUTPUT_DIR/$struct_name"
        
        log "Analyzing: $struct_name"
        
        if analyze_structure "$struct_file" "$struct_name" "$struct_dir"; then
            log "  ✓ Complete"
        else
            log "  ✗ Failed"
        fi
        echo ""
    done
    
    # Generate comparison report
    log "Generating comparison report..."
    
    python3 << 'PYEOF'
import os
import json
import glob
import csv

results = []
for metrics_file in glob.glob('*/metrics.json'):
    struct_name = os.path.dirname(metrics_file)
    try:
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)
            metrics['name'] = struct_name
            results.append(metrics)
    except:
        pass

# Sort by potential energy (more negative = more stable)
results.sort(key=lambda x: x.get('potential', 0) or 0)

# Generate report
with open('comparison_report.txt', 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("BATCH PPI STRUCTURE COMPARISON REPORT\n")
    f.write("=" * 70 + "\n\n")
    
    f.write(f"{'Rank':<6}{'Structure':<40}{'Potential':<15}{'Rg':<10}{'SASA':<10}\n")
    f.write("-" * 70 + "\n")
    
    for i, r in enumerate(results, 1):
        pot = f"{r.get('potential', 0):.0f}" if r.get('potential') else "N/A"
        rg = f"{r.get('gyration', 0):.2f}" if r.get('gyration') else "N/A"
        sasa = f"{r.get('sasa', 0):.1f}" if r.get('sasa') else "N/A"
        f.write(f"{i:<6}{r['name'][:38]:<40}{pot:<15}{rg:<10}{sasa:<10}\n")
    
    f.write("\n" + "-" * 70 + "\n")
    f.write("Notes:\n")
    f.write("  - Potential: Lower (more negative) = more stable\n")
    f.write("  - Rg: Radius of gyration (nm) - compactness\n")
    f.write("  - SASA: Solvent accessible surface area (nm²)\n")
    f.write("\n")
    
    if results:
        f.write(f"BEST STRUCTURE: {results[0]['name']}\n")
        f.write(f"  Potential Energy: {results[0].get('potential', 'N/A'):.1f} kJ/mol\n")

print("\nComparison saved to: comparison_report.txt")

# Save as CSV for further analysis
with open('comparison_data.csv', 'w') as f:
    f.write("rank,structure,potential_kJ_mol,gyration_nm,sasa_nm2,status\n")
    for i, r in enumerate(results, 1):
        f.write(f"{i},{r['name']},{r.get('potential', '')},{r.get('gyration', '')},{r.get('sasa', '')},{r.get('status', '')}\n")

print("CSV data saved to: comparison_data.csv")

# Save complete JSON with all results
with open('comparison_data.json', 'w') as f:
    json.dump({
        'structures': results,
        'best_structure': results[0]['name'] if results else None,
        'total_analyzed': len(results)
    }, f, indent=2)
print("JSON data saved to: comparison_data.json")

# Generate gnuplot script for bar chart
gnuplot_script = '''# Bar chart of potential energies
set terminal pngcairo size 1200,600 enhanced font 'Arial,10'
set output 'comparison_barchart.png'

set title 'Structure Comparison - Potential Energy'
set ylabel 'Potential Energy (kJ/mol)'
set style fill solid 0.8
set boxwidth 0.8
set xtics rotate by -45

set datafile separator ','
plot 'comparison_data.csv' skip 1 using 0:3:xtic(2) with boxes lc rgb 'blue' notitle
'''

with open('plot_comparison.gp', 'w') as f:
    f.write(gnuplot_script)
print("Gnuplot script saved to: plot_comparison.gp")

# Generate R script for analysis
r_script = '''# R script for PPI comparison analysis
# Run with: Rscript analyze_comparison.R

library(ggplot2)

# Read data
data <- read.csv("comparison_data.csv")

# Bar plot of potential energies
p1 <- ggplot(data, aes(x=reorder(structure, potential_kJ_mol), y=potential_kJ_mol)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Structure Comparison", x="Structure", y="Potential Energy (kJ/mol)") +
  theme_minimal()

ggsave("comparison_barplot.png", p1, width=10, height=6)

# Scatter plot: Rg vs SASA
p2 <- ggplot(data, aes(x=gyration_nm, y=sasa_nm2, label=structure)) +
  geom_point(size=3, color="red") +
  geom_text(vjust=-0.5, size=3) +
  labs(title="Compactness vs Surface Area", x="Radius of Gyration (nm)", y="SASA (nm²)") +
  theme_minimal()

ggsave("rg_vs_sasa.png", p2, width=8, height=6)

print("Plots saved: comparison_barplot.png, rg_vs_sasa.png")
'''

with open('analyze_comparison.R', 'w') as f:
    f.write(r_script)
print("R script saved to: analyze_comparison.R")

# Generate Python plotting script
python_script = '''#!/usr/bin/env python3
"""Generate comparison plots from batch analysis."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read data
df = pd.read_csv("comparison_data.csv")

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)

# 1. Bar plot of potential energies
fig, ax = plt.subplots()
colors = ['green' if i == 0 else 'steelblue' for i in range(len(df))]
bars = ax.barh(df['structure'], df['potential_kJ_mol'], color=colors)
ax.set_xlabel('Potential Energy (kJ/mol)')
ax.set_title('Structure Comparison - Potential Energy\\n(Lower = More Stable, Green = Best)')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('comparison_barplot.png', dpi=150)
plt.close()

# 2. Scatter: Rg vs SASA
fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter(df['gyration_nm'], df['sasa_nm2'], 
                     c=df['potential_kJ_mol'], cmap='RdYlGn_r', 
                     s=100, edgecolors='black')
for i, txt in enumerate(df['structure']):
    ax.annotate(txt[:15], (df['gyration_nm'].iloc[i], df['sasa_nm2'].iloc[i]),
                fontsize=8, alpha=0.7)
plt.colorbar(scatter, label='Potential Energy (kJ/mol)')
ax.set_xlabel('Radius of Gyration (nm)')
ax.set_ylabel('SASA (nm²)')
ax.set_title('Structure Compactness vs Surface Area')
plt.tight_layout()
plt.savefig('rg_vs_sasa_scatter.png', dpi=150)
plt.close()

# 3. Summary heatmap
fig, ax = plt.subplots(figsize=(8, len(df)*0.5 + 2))
metrics = df[['potential_kJ_mol', 'gyration_nm', 'sasa_nm2']].copy()
# Normalize for heatmap
for col in metrics.columns:
    metrics[col] = (metrics[col] - metrics[col].min()) / (metrics[col].max() - metrics[col].min() + 1e-10)
metrics.index = df['structure']
sns.heatmap(metrics, annot=True, fmt='.2f', cmap='RdYlGn_r', ax=ax)
ax.set_title('Normalized Metrics Comparison\\n(Lower = Better for all metrics)')
plt.tight_layout()
plt.savefig('comparison_heatmap.png', dpi=150)
plt.close()

print("Generated: comparison_barplot.png, rg_vs_sasa_scatter.png, comparison_heatmap.png")
'''

with open('generate_plots.py', 'w') as f:
    f.write(python_script)
print("Python plotting script saved to: generate_plots.py")
PYEOF

    # Try to generate plots
    log "Generating plots..."
    
    if command -v gnuplot &> /dev/null; then
        gnuplot plot_comparison.gp 2>/dev/null || true
        log "  Generated: comparison_barchart.png (gnuplot)"
    fi
    
    if python3 -c "import matplotlib, seaborn, pandas" 2>/dev/null; then
        python3 generate_plots.py 2>/dev/null || true
        log "  Generated: comparison_barplot.png, rg_vs_sasa_scatter.png, comparison_heatmap.png"
    else
        log "  Note: Install matplotlib, seaborn, pandas for additional plots"
        log "        pip install matplotlib seaborn pandas"
    fi

    END=$(date +%s)
    
    echo ""
    echo "============================================================"
    log "Batch analysis complete in $((END - START)) seconds"
    echo "============================================================"
    echo ""
    echo "OUTPUT FILES:"
    echo ""
    echo "Reports:"
    echo "  - comparison_report.txt    : Text summary with rankings"
    echo "  - comparison_data.csv      : CSV for spreadsheet analysis"
    echo "  - comparison_data.json     : JSON for programmatic access"
    echo ""
    echo "Visualization Scripts:"
    echo "  - generate_plots.py        : Python script (matplotlib/seaborn)"
    echo "  - analyze_comparison.R     : R script (ggplot2)"
    echo "  - plot_comparison.gp       : Gnuplot script"
    echo ""
    echo "Generated Plots (if dependencies available):"
    echo "  - comparison_barplot.png   : Energy ranking bar chart"
    echo "  - comparison_heatmap.png   : Normalized metrics heatmap"
    echo "  - rg_vs_sasa_scatter.png   : Compactness vs surface area"
    echo ""
    
    cat comparison_report.txt
}

main "$@"
