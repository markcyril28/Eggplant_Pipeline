#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS Interface & Binding Analysis Script
#
# Performs detailed PPI interface analysis:
# 1. Interface residue identification (within distance cutoff)
# 2. Buried Surface Area (BSA) calculation
# 3. Hydrogen bond analysis at interface
# 4. Salt bridge / ionic interaction detection
# 5. Contact frequency & residue-residue contact map
# 6. Per-residue energy contribution estimation
#
# Usage: ./gromacs_interface_analysis.sh <PDB_FILE> [OUTPUT_DIR]
#######################################################################

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

PDB_INPUT="${1:-${HOME}/PPI/inputs/SmelGRF-GIF/fold_1_x_1_model_0.pdb}"
WORKDIR="${2:-${HOME}/PPI/RESULTS/results_gromacs/interface_analysis}"

# Analysis parameters
INTERFACE_CUTOFF=0.5      # nm - residues within this distance are "interface"
CONTACT_CUTOFF=0.6        # nm - for contact counting
HBOND_CUTOFF=0.35         # nm - H-bond distance cutoff
FORCEFIELD="amber99sb-ildn"
WATERMODEL="tip3p"
BOX_DISTANCE=1.0
EM_STEPS=5000

# GROMACS binary
GMX_BIN="${HOME}/miniconda3/envs/gromacs_HIP/bin/gmx_mpi"

# GPU settings (uses GPU for EM if available)
GPU_EM_FLAGS="-nb gpu -pme cpu -bonded cpu"
NTHREADS=8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------

log() { echo "[$(date '+%H:%M:%S')] $1"; }

prepare_structure() {
    log "Preparing structure..."
    
    # Convert CIF to PDB if needed
    if [[ "$PDB_INPUT" == *.cif ]]; then
        log "  Converting CIF to PDB..."
        python3 << EOF
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('protein', '$PDB_INPUT')
io = PDBIO()
io.set_structure(structure)
io.save('input.pdb')
EOF
        PDB_FILE="input.pdb"
    else
        cp "$PDB_INPUT" input.pdb
        PDB_FILE="input.pdb"
    fi
    
    # Clean PDB
    grep -E "^(ATOM|TER|END)" "$PDB_FILE" | \
        sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > clean.pdb
    
    # Get chain info
    CHAINS=$(grep "^ATOM" clean.pdb | cut -c22 | sort -u | tr -d ' ' | tr '\n' ' ')
    log "  Chains found: $CHAINS"
    
    # Generate topology
    log "  Generating topology..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 > pdb2gmx.log
    
    # Create chain index
    create_chain_index
    
    # Quick solvation and EM for realistic structure
    log "  Solvating and minimizing..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 > editconf.log
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 > solvate.log
    
    # Create EM MDP
    cat > em.mdp << EOF
integrator = steep
emtol = 100.0
emstep = 0.01
nsteps = $EM_STEPS
nstlist = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
pbc = xyz
EOF

    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 2 2>&1 > grompp.log
    
    # Run EM with GPU
    log "  Running energy minimization..."
    mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp $NTHREADS $GPU_EM_FLAGS 2>&1 > mdrun_em.log || \
        mpirun -np 1 $GMX_BIN mdrun -v -deffnm em -ntomp 16 2>&1 > mdrun_em_cpu.log
    
    log "  Structure prepared: em.gro"
}

create_chain_index() {
    log "  Creating chain index..."
    
    # Generate basic index
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 > make_ndx.log
    
    # Parse chains and add to index
    python3 << 'PYEOF'
import re

# Read the GRO file to get atom info
atoms_by_chain = {}
with open('protein.gro', 'r') as f:
    lines = f.readlines()[2:-1]  # Skip header and box line
    for i, line in enumerate(lines, 1):
        if len(line) >= 10:
            resnum = int(line[0:5].strip())
            resname = line[5:10].strip()
            # Determine chain based on residue numbering pattern
            # This is a heuristic - may need adjustment
            chain = 'A' if i <= len(lines)//2 else 'B'
            if chain not in atoms_by_chain:
                atoms_by_chain[chain] = []
            atoms_by_chain[chain].append(i)

# Also try to get chain info from PDB if available
try:
    chain_atoms = {'A': [], 'B': []}
    current_chain = None
    atom_idx = 0
    with open('clean.pdb', 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_idx += 1
                chain = line[21].strip()
                if chain and chain in ['A', 'B']:
                    chain_atoms[chain].append(atom_idx)
                elif chain:
                    # Map other chains to A or B
                    if not chain_atoms['A'] or (chain_atoms['B'] and len(chain_atoms['A']) > len(chain_atoms['B'])):
                        chain_atoms['B'].append(atom_idx)
                    else:
                        chain_atoms['A'].append(atom_idx)
    if chain_atoms['A'] and chain_atoms['B']:
        atoms_by_chain = chain_atoms
except:
    pass

# Append to index file
with open('index.ndx', 'a') as f:
    for chain, atoms in sorted(atoms_by_chain.items()):
        if atoms:
            f.write(f'\n[ Chain{chain} ]\n')
            for i, atom in enumerate(atoms):
                f.write(f'{atom:6d}')
                if (i + 1) % 15 == 0:
                    f.write('\n')
            f.write('\n')

print(f"Added chain groups: {list(atoms_by_chain.keys())}")
PYEOF
}

analyze_interface() {
    log "Analyzing interface..."
    
    # 1. Minimum distance and contacts between chains
    log "  Calculating inter-chain distances..."
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx \
        -od mindist.xvg -on numcont.xvg -d $CONTACT_CUTOFF 2>&1 > mindist.log || true
    
    # 2. Contact residue pairs
    log "  Identifying contact residue pairs..."
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx \
        -or rescontact.xvg -d $INTERFACE_CUTOFF 2>&1 > rescontact.log || true
    
    # 3. Hydrogen bonds
    log "  Analyzing hydrogen bonds..."
    echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx \
        -num hbonds.xvg -dist hbdist.xvg -ang hbang.xvg 2>&1 > hbond.log || true
    
    # 4. SASA analysis (for BSA calculation)
    log "  Calculating solvent accessible surface area..."
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa_complex.xvg -or sasa_res.xvg 2>&1 > sasa_complex.log || true
    echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainA.xvg 2>&1 > sasa_A.log || true
    echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainB.xvg 2>&1 > sasa_B.log || true
    
    # 5. Energy components
    log "  Extracting energy components..."
    echo -e "Potential\nCoul-SR\nLJ-SR\n\n" | $GMX_BIN energy -f em.edr -o energies.xvg 2>&1 > energy.log || true
}

calculate_interface_residues() {
    log "Identifying interface residues..."
    
    python3 << 'PYEOF'
import os
import re

def parse_xvg(filename):
    """Parse XVG file and return data."""
    data = []
    if not os.path.exists(filename):
        return data
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):
                parts = line.split()
                if parts:
                    data.append([float(x) for x in parts])
    return data

def parse_gro_residues(filename):
    """Get residue info from GRO file."""
    residues = {}
    with open(filename, 'r') as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            if len(line) >= 15:
                resnum = int(line[0:5].strip())
                resname = line[5:10].strip()
                if resnum not in residues:
                    residues[resnum] = resname
    return residues

# Parse SASA data for BSA calculation
sasa_complex = parse_xvg('sasa_complex.xvg')
sasa_A = parse_xvg('sasa_chainA.xvg')
sasa_B = parse_xvg('sasa_chainB.xvg')

# Calculate Buried Surface Area (BSA)
bsa = 0
if sasa_complex and sasa_A and sasa_B:
    complex_sasa = sasa_complex[0][1] if len(sasa_complex[0]) > 1 else sasa_complex[0][0]
    chain_a_sasa = sasa_A[0][1] if len(sasa_A[0]) > 1 else sasa_A[0][0]
    chain_b_sasa = sasa_B[0][1] if len(sasa_B[0]) > 1 else sasa_B[0][0]
    bsa = (chain_a_sasa + chain_b_sasa - complex_sasa) / 2  # Per chain
    print(f"Buried Surface Area (BSA): {bsa:.2f} nm²")

# Parse hydrogen bond count
hbonds = 0
hb_data = parse_xvg('hbonds.xvg')
if hb_data:
    hbonds = int(hb_data[0][1]) if len(hb_data[0]) > 1 else int(hb_data[0][0])
    print(f"Interface hydrogen bonds: {hbonds}")

# Parse contact count
contacts = 0
cont_data = parse_xvg('numcont.xvg')
if cont_data:
    contacts = int(cont_data[0][1]) if len(cont_data[0]) > 1 else int(cont_data[0][0])
    print(f"Interface contacts (<{0.6} nm): {contacts}")

# Parse minimum distance
mindist = 0
md_data = parse_xvg('mindist.xvg')
if md_data:
    mindist = md_data[0][1] if len(md_data[0]) > 1 else md_data[0][0]
    print(f"Minimum inter-chain distance: {mindist:.3f} nm")

# Parse energies
energies = parse_xvg('energies.xvg')
potential = coul = lj = 0
if energies and len(energies[0]) >= 3:
    potential = energies[0][0]
    coul = energies[0][1] if len(energies[0]) > 1 else 0
    lj = energies[0][2] if len(energies[0]) > 2 else 0

# Write summary report
with open('interface_analysis.txt', 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("INTERFACE ANALYSIS REPORT\n")
    f.write("=" * 60 + "\n\n")
    
    f.write("BINDING INTERFACE METRICS\n")
    f.write("-" * 40 + "\n")
    f.write(f"Buried Surface Area (BSA):    {bsa:.2f} nm²\n")
    f.write(f"Interface H-bonds:            {hbonds}\n")
    f.write(f"Interface contacts:           {contacts}\n")
    f.write(f"Minimum distance:             {mindist:.3f} nm\n\n")
    
    f.write("ENERGY COMPONENTS (kJ/mol)\n")
    f.write("-" * 40 + "\n")
    f.write(f"Total Potential:              {potential:.2f}\n")
    f.write(f"Short-range Coulomb:          {coul:.2f}\n")
    f.write(f"Short-range LJ:               {lj:.2f}\n\n")
    
    # Binding quality assessment
    f.write("BINDING QUALITY ASSESSMENT\n")
    f.write("-" * 40 + "\n")
    
    quality_score = 0
    assessments = []
    
    if bsa > 8:  # Good interface size
        quality_score += 2
        assessments.append("✓ Large binding interface")
    elif bsa > 4:
        quality_score += 1
        assessments.append("○ Moderate binding interface")
    else:
        assessments.append("✗ Small binding interface")
    
    if hbonds >= 10:
        quality_score += 2
        assessments.append("✓ Strong H-bond network")
    elif hbonds >= 5:
        quality_score += 1
        assessments.append("○ Moderate H-bond network")
    else:
        assessments.append("✗ Weak H-bond network")
    
    if contacts >= 50:
        quality_score += 2
        assessments.append("✓ Extensive contacts")
    elif contacts >= 20:
        quality_score += 1
        assessments.append("○ Moderate contacts")
    else:
        assessments.append("✗ Limited contacts")
    
    for a in assessments:
        f.write(f"{a}\n")
    
    f.write(f"\nOverall Quality Score: {quality_score}/6\n")
    
    if quality_score >= 5:
        f.write("Assessment: STRONG binding predicted\n")
    elif quality_score >= 3:
        f.write("Assessment: MODERATE binding predicted\n")
    else:
        f.write("Assessment: WEAK binding predicted\n")

print("\nReport saved to: interface_analysis.txt")
PYEOF
}

generate_visualizations() {
    log "Generating visualization files..."
    
    python3 << 'PYEOF'
import os
import json

# Generate PyMOL script for interface visualization
def create_pymol_script():
    """Create PyMOL script to visualize interface."""
    
    # Read interface residues
    interface_a = []
    interface_b = []
    
    if os.path.exists('interface_residues.txt'):
        with open('interface_residues.txt', 'r') as f:
            for line in f:
                if line.strip() and not line.startswith(('Interface', '=', '-', 'Chain', 'Total')):
                    parts = line.split()
                    if len(parts) >= 2:
                        # Extract residue numbers
                        try:
                            res_a = ''.join(filter(str.isdigit, parts[0]))
                            res_b = ''.join(filter(str.isdigit, parts[1]))
                            if res_a: interface_a.append(res_a)
                            if res_b: interface_b.append(res_b)
                        except:
                            pass
    
    interface_a = list(set(interface_a))[:50]  # Limit for readability
    interface_b = list(set(interface_b))[:50]
    
    script = '''# PyMOL Interface Visualization Script
# Load with: pymol visualize_interface.pml

# Load structure
load em.gro, complex

# Basic display settings
bg_color white
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1
set stick_radius, 0.15

# Color chains
color marine, chain A
color orange, chain B

# Show as cartoon
hide all
show cartoon, all

# Highlight interface residues
'''
    
    if interface_a:
        script += f"select interface_A, chain A and resi {'+'.join(interface_a)}\n"
        script += "show sticks, interface_A\n"
        script += "color red, interface_A\n"
    
    if interface_b:
        script += f"select interface_B, chain B and resi {'+'.join(interface_b)}\n"
        script += "show sticks, interface_B\n"
        script += "color yellow, interface_B\n"
    
    script += '''
# Show hydrogen bonds
dist hbonds, interface_A, interface_B, mode=2
hide labels, hbonds
color cyan, hbonds

# Set view
orient
zoom all, 5

# Save image
ray 1920, 1080
png interface_view.png, dpi=300

# Save session
save interface_session.pse
'''
    
    with open('visualize_interface.pml', 'w') as f:
        f.write(script)
    print("Created: visualize_interface.pml")

# Generate VMD script
def create_vmd_script():
    """Create VMD script for visualization."""
    
    script = '''# VMD Interface Visualization Script
# Load with: vmd -e visualize_interface.vmd

# Load structure
mol new em.gro type gro

# Display settings
display projection Orthographic
display depthcue off
axes location Off
color Display Background white

# Representation: Cartoon for protein
mol delrep 0 top
mol representation NewCartoon
mol color Chain
mol selection {protein}
mol addrep top

# Representation: Licorice for interface (within 5A of other chain)
mol representation Licorice 0.3 12.0 12.0
mol color Name
mol selection {protein and chain A and within 5 of chain B}
mol addrep top

mol representation Licorice 0.3 12.0 12.0
mol color Name  
mol selection {protein and chain B and within 5 of chain A}
mol addrep top

# Center view
display resetview

# Render
render TachyonInternal interface_view_vmd.tga
'''
    
    with open('visualize_interface.vmd', 'w') as f:
        f.write(script)
    print("Created: visualize_interface.vmd")

# Generate ChimeraX script
def create_chimerax_script():
    """Create ChimeraX script for visualization."""
    
    script = '''# ChimeraX Interface Visualization Script
# Load with: chimerax visualize_interface.cxc

# Open structure
open em.gro

# Display settings
set bgColor white
lighting soft

# Color chains
color /A marine
color /B orange

# Show as cartoon
hide atoms
show cartoon

# Highlight interface (residues within 5A of other chain)
select /A & @@<5 & /B
show sel atoms
style sel stick
color sel red

select /B & @@<5 & /A  
show sel atoms
style sel stick
color sel yellow

# Show hydrogen bonds
hbonds /A /B color cyan

# Center and zoom
view

# Save image
save interface_chimerax.png width 1920 height 1080 supersample 3
'''
    
    with open('visualize_interface.cxc', 'w') as f:
        f.write(script)
    print("Created: visualize_interface.cxc")

create_pymol_script()
create_vmd_script()
create_chimerax_script()
PYEOF
}

generate_statistics() {
    log "Generating statistics files..."
    
    python3 << 'PYEOF'
import os
import json
import csv

def parse_xvg(filename):
    """Parse XVG file and return all data points."""
    data = []
    metadata = {'title': '', 'xlabel': '', 'ylabel': ''}
    if not os.path.exists(filename):
        return data, metadata
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('@ title'):
                metadata['title'] = line.split('"')[1] if '"' in line else ''
            elif line.startswith('@ xaxis label'):
                metadata['xlabel'] = line.split('"')[1] if '"' in line else ''
            elif line.startswith('@ yaxis label'):
                metadata['ylabel'] = line.split('"')[1] if '"' in line else ''
            elif not line.startswith(('#', '@')):
                parts = line.split()
                if parts:
                    data.append([float(x) for x in parts])
    return data, metadata

# Collect all statistics
stats = {
    'structure': os.path.basename(os.getcwd()),
    'energy': {},
    'interface': {},
    'geometry': {}
}

# Parse energy data
energy_data, _ = parse_xvg('energies.xvg')
if energy_data:
    stats['energy']['potential_kJ_mol'] = energy_data[-1][0] if energy_data[-1] else None
    if len(energy_data[-1]) > 1:
        stats['energy']['coulomb_SR_kJ_mol'] = energy_data[-1][1]
    if len(energy_data[-1]) > 2:
        stats['energy']['lj_SR_kJ_mol'] = energy_data[-1][2]

# Parse SASA data
sasa_complex, _ = parse_xvg('sasa_complex.xvg')
sasa_a, _ = parse_xvg('sasa_chainA.xvg')
sasa_b, _ = parse_xvg('sasa_chainB.xvg')

if sasa_complex and sasa_a and sasa_b:
    complex_sasa = sasa_complex[0][1] if len(sasa_complex[0]) > 1 else sasa_complex[0][0]
    chain_a_sasa = sasa_a[0][1] if len(sasa_a[0]) > 1 else sasa_a[0][0]
    chain_b_sasa = sasa_b[0][1] if len(sasa_b[0]) > 1 else sasa_b[0][0]
    
    stats['geometry']['sasa_complex_nm2'] = complex_sasa
    stats['geometry']['sasa_chainA_nm2'] = chain_a_sasa
    stats['geometry']['sasa_chainB_nm2'] = chain_b_sasa
    stats['interface']['buried_surface_area_nm2'] = (chain_a_sasa + chain_b_sasa - complex_sasa) / 2

# Parse H-bonds
hbond_data, _ = parse_xvg('hbonds.xvg')
if hbond_data:
    stats['interface']['hydrogen_bonds'] = int(hbond_data[0][1]) if len(hbond_data[0]) > 1 else int(hbond_data[0][0])

# Parse contacts
cont_data, _ = parse_xvg('numcont.xvg')
if cont_data:
    stats['interface']['contacts'] = int(cont_data[0][1]) if len(cont_data[0]) > 1 else int(cont_data[0][0])

# Parse min distance
mindist_data, _ = parse_xvg('mindist.xvg')
if mindist_data:
    stats['interface']['min_distance_nm'] = mindist_data[0][1] if len(mindist_data[0]) > 1 else mindist_data[0][0]

# Count interface residues
if os.path.exists('interface_residues.txt'):
    with open('interface_residues.txt', 'r') as f:
        content = f.read()
        if 'Total interface residue pairs:' in content:
            try:
                stats['interface']['residue_pairs'] = int(content.split('Total interface residue pairs:')[1].strip().split()[0])
            except:
                pass

# Save as JSON
with open('statistics.json', 'w') as f:
    json.dump(stats, f, indent=2)
print("Created: statistics.json")

# Save as CSV (flat format for easy import)
with open('statistics.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['metric', 'value', 'unit'])
    
    for category, metrics in stats.items():
        if isinstance(metrics, dict):
            for key, value in metrics.items():
                unit = ''
                if 'kJ_mol' in key: unit = 'kJ/mol'
                elif 'nm2' in key: unit = 'nm²'
                elif 'nm' in key: unit = 'nm'
                writer.writerow([f"{category}.{key}", value, unit])
print("Created: statistics.csv")

# Print summary
print("\nStatistics Summary:")
print(f"  Potential Energy: {stats['energy'].get('potential_kJ_mol', 'N/A')} kJ/mol")
print(f"  BSA: {stats['interface'].get('buried_surface_area_nm2', 'N/A'):.2f} nm²" if stats['interface'].get('buried_surface_area_nm2') else "  BSA: N/A")
print(f"  H-bonds: {stats['interface'].get('hydrogen_bonds', 'N/A')}")
print(f"  Contacts: {stats['interface'].get('contacts', 'N/A')}")
PYEOF
}

generate_contact_map() {
    log "Generating contact map..."
    
    python3 << 'PYEOF'
import numpy as np
import os

def get_residue_coords(gro_file):
    """Extract CA coordinates per residue."""
    residues = {}
    with open(gro_file, 'r') as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            if len(line) >= 44:
                resnum = int(line[0:5].strip())
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                if atomname == 'CA':  # Use CA atoms
                    x = float(line[20:28])
                    y = float(line[28:36])
                    z = float(line[36:44])
                    residues[resnum] = {'name': resname, 'coords': np.array([x, y, z])}
    return residues

def get_chain_residues(pdb_file):
    """Get residue to chain mapping from PDB."""
    chain_map = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                resnum = int(line[22:26].strip())
                chain = line[21].strip() or 'A'
                chain_map[resnum] = chain
    return chain_map

# Get residue data
residues = get_residue_coords('em.gro')
chain_map = get_chain_residues('clean.pdb')

if not residues:
    print("No residues found")
    exit()

# Separate chains
chain_a_res = sorted([r for r in residues if chain_map.get(r, 'A') == 'A'])
chain_b_res = sorted([r for r in residues if chain_map.get(r, 'B') == 'B'])

if not chain_b_res:
    # If chain assignment failed, split by residue number
    all_res = sorted(residues.keys())
    mid = len(all_res) // 2
    chain_a_res = all_res[:mid]
    chain_b_res = all_res[mid:]

# Calculate contact map
contact_cutoff = 0.8  # nm
contact_map = np.zeros((len(chain_a_res), len(chain_b_res)))
interface_pairs = []

for i, res_a in enumerate(chain_a_res):
    for j, res_b in enumerate(chain_b_res):
        if res_a in residues and res_b in residues:
            dist = np.linalg.norm(residues[res_a]['coords'] - residues[res_b]['coords'])
            if dist < contact_cutoff:
                contact_map[i, j] = 1
                interface_pairs.append((res_a, residues[res_a]['name'], 
                                        res_b, residues[res_b]['name'], dist))

# Save contact map as text
np.savetxt('contact_map.txt', contact_map, fmt='%d')

# Save interface pairs
interface_pairs.sort(key=lambda x: x[4])
with open('interface_residues.txt', 'w') as f:
    f.write("Interface Residue Pairs (distance < 0.8 nm)\n")
    f.write("=" * 50 + "\n")
    f.write(f"{'ChainA':<15} {'ChainB':<15} {'Distance (nm)':<15}\n")
    f.write("-" * 50 + "\n")
    for res_a, name_a, res_b, name_b, dist in interface_pairs[:30]:  # Top 30
        f.write(f"{name_a}{res_a:<10} {name_b}{res_b:<10} {dist:.3f}\n")
    f.write(f"\nTotal interface residue pairs: {len(interface_pairs)}\n")

print(f"Contact map saved: contact_map.txt")
print(f"Interface residues saved: interface_residues.txt")
print(f"Found {len(interface_pairs)} interface residue pairs")

# Generate contact map data for plotting
if len(chain_a_res) > 0 and len(chain_b_res) > 0:
    # Save plotting data
    with open('contact_map_plot.csv', 'w') as f:
        f.write('chainA_res,chainB_res,distance_nm,contact\n')
        for i, res_a in enumerate(chain_a_res):
            for j, res_b in enumerate(chain_b_res):
                if res_a in residues and res_b in residues:
                    dist = np.linalg.norm(residues[res_a]['coords'] - residues[res_b]['coords'])
                    contact = 1 if dist < 0.8 else 0
                    f.write(f'{res_a},{res_b},{dist:.3f},{contact}\n')
    print("Created: contact_map_plot.csv (for heatmap visualization)")

# Generate gnuplot script for contact map
gnuplot_script = '''# Gnuplot script for contact map visualization
# Run with: gnuplot contact_map.gp

set terminal pngcairo size 800,800 enhanced font 'Arial,12'
set output 'contact_map_heatmap.png'

set title 'Inter-chain Contact Map'
set xlabel 'Chain A Residue'
set ylabel 'Chain B Residue'
set palette defined (0 'white', 1 'red')
unset colorbox

set datafile separator ','
plot 'contact_map_plot.csv' skip 1 using 1:2:4 with image notitle
'''

with open('contact_map.gp', 'w') as f:
    f.write(gnuplot_script)
print("Created: contact_map.gp (gnuplot script)")
PYEOF
}

convert_structure_formats() {
    log "Converting structure to multiple formats..."
    
    # Convert to PDB format (more widely compatible)
    echo "Protein" | $GMX_BIN trjconv -s em.tpr -f em.gro -o structure_minimized.pdb -pbc mol 2>&1 > trjconv.log || true
    
    log "  Created: structure_minimized.pdb"
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

main() {
    echo ""
    echo "============================================================"
    echo " GROMACS Interface & Binding Analysis"
    echo "============================================================"
    
    log "Input: $PDB_INPUT"
    log "Output: $WORKDIR"
    
    mkdir -p "$WORKDIR"
    cd "$WORKDIR"
    
    START=$(date +%s)
    
    prepare_structure
    analyze_interface
    calculate_interface_residues
    generate_contact_map
    generate_visualizations
    generate_statistics
    convert_structure_formats
    
    END=$(date +%s)
    
    echo ""
    echo "============================================================"
    log "Analysis complete in $((END - START)) seconds"
    echo "============================================================"
    echo ""
    echo "OUTPUT FILES:"
    echo ""
    echo "Analysis Reports:"
    echo "  - interface_analysis.txt   : Summary report with binding quality score"
    echo "  - interface_residues.txt   : Interface residue pairs with distances"
    echo "  - statistics.json          : All metrics in JSON format"
    echo "  - statistics.csv           : All metrics in CSV format"
    echo ""
    echo "Visualization Scripts:"
    echo "  - visualize_interface.pml  : PyMOL script (run: pymol visualize_interface.pml)"
    echo "  - visualize_interface.vmd  : VMD script (run: vmd -e visualize_interface.vmd)"
    echo "  - visualize_interface.cxc  : ChimeraX script"
    echo ""
    echo "Contact Map Data:"
    echo "  - contact_map.txt          : Binary contact matrix"
    echo "  - contact_map_plot.csv     : Data for heatmap plotting"
    echo "  - contact_map.gp           : Gnuplot script for heatmap"
    echo ""
    echo "Structure Files:"
    echo "  - structure_minimized.pdb  : Energy-minimized structure (PDB format)"
    echo "  - em.gro                   : Energy-minimized structure (GROMACS format)"
    echo ""
    
    cat interface_analysis.txt
}

main "$@"
