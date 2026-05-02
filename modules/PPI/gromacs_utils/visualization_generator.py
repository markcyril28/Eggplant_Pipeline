#!/usr/bin/env python3
"""
Visualization Generator Module for GROMACS PPI Analysis

Generates visualization scripts for PyMOL, VMD, ChimeraX, and plotting tools.
"""

import os
from pathlib import Path
from typing import List, Dict, Optional, Any


class VisualizationGenerator:
    """Generate visualization scripts for molecular structures."""
    
    def __init__(self, workdir: Path, structure_file: str = "em.gro"):
        """
        Initialize visualization generator.
        
        Args:
            workdir: Working directory
            structure_file: Name of the structure file (GRO or PDB)
        """
        self.workdir = Path(workdir)
        self.structure_file = structure_file
        self.interface_residues_a = []
        self.interface_residues_b = []
        
        # Try to load interface residues
        self._load_interface_residues()
    
    def _load_interface_residues(self):
        """Load interface residues from file if available."""
        interface_file = self.workdir / 'interface_residues.txt'
        if interface_file.exists():
            with open(interface_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith(('Interface', '=', '-', 'Chain', 'Total')):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                res_a = ''.join(filter(str.isdigit, parts[0]))
                                res_b = ''.join(filter(str.isdigit, parts[1]))
                                if res_a:
                                    self.interface_residues_a.append(res_a)
                                if res_b:
                                    self.interface_residues_b.append(res_b)
                            except (ValueError, IndexError):
                                pass
            
            # Remove duplicates and limit
            self.interface_residues_a = list(set(self.interface_residues_a))[:50]
            self.interface_residues_b = list(set(self.interface_residues_b))[:50]
    
    def generate_pymol_interface(self, output_file: Optional[Path] = None) -> str:
        """
        Generate PyMOL script for interface visualization.
        
        Args:
            output_file: Optional path to save script
            
        Returns:
            Script content as string
        """
        script = f'''# PyMOL Interface Visualization Script
# Load with: pymol visualize_interface.pml

# Load structure
load {self.structure_file}, complex

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
        if self.interface_residues_a:
            script += f"select interface_A, chain A and resi {'+'.join(self.interface_residues_a)}\n"
            script += "show sticks, interface_A\n"
            script += "color red, interface_A\n"
        
        if self.interface_residues_b:
            script += f"select interface_B, chain B and resi {'+'.join(self.interface_residues_b)}\n"
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
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        
        return script
    
    def generate_pymol_trajectory(self, trajectory_file: str = "md_protein.xtc",
                                  pdb_file: str = "final_structure.pdb",
                                  output_file: Optional[Path] = None) -> str:
        """
        Generate PyMOL script for trajectory visualization.
        
        Args:
            trajectory_file: Name of trajectory file
            pdb_file: Name of PDB structure file
            output_file: Optional path to save script
            
        Returns:
            Script content as string
        """
        script = f'''# PyMOL MD Trajectory Visualization
# Load with: pymol visualize_trajectory.pml

# Load trajectory
load {pdb_file}, protein
load_traj {trajectory_file}, protein, state=1

# Display settings
bg_color white
set cartoon_fancy_helices, 1
set movie_fps, 15

# Color by chain
color marine, chain A
color orange, chain B

# Show as cartoon
hide all
show cartoon

# Center and zoom
orient
zoom all, 5

# Movie controls
mset 1 x100
mview store, 1, state=1
mview store, 100, state=0
mplay

# To save movie (uncomment):
# set ray_trace_frames, 1
# mpng frame_, width=1280, height=720
'''
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        
        return script
    
    def generate_vmd_interface(self, output_file: Optional[Path] = None) -> str:
        """
        Generate VMD script for interface visualization.
        
        Args:
            output_file: Optional path to save script
            
        Returns:
            Script content as string
        """
        script = f'''# VMD Interface Visualization Script
# Load with: vmd -e visualize_interface.vmd

# Load structure
mol new {self.structure_file} type gro

# Display settings
display projection Orthographic
display depthcue off
axes location Off
color Display Background white

# Representation: Cartoon for protein
mol delrep 0 top
mol representation NewCartoon
mol color Chain
mol selection {{protein}}
mol addrep top

# Representation: Licorice for interface (within 5A of other chain)
mol representation Licorice 0.3 12.0 12.0
mol color Name
mol selection {{protein and chain A and within 5 of chain B}}
mol addrep top

mol representation Licorice 0.3 12.0 12.0
mol color Name  
mol selection {{protein and chain B and within 5 of chain A}}
mol addrep top

# Center view
display resetview

# Render
render TachyonInternal interface_view_vmd.tga
'''
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        
        return script
    
    def generate_vmd_movie_script(self,
                                  trajectory_file: str = "trajectories/md_center.xtc",
                                  pdb_file: str = "structures/final_structure.pdb",
                                  frames_dir: str = "video/frames",
                                  output_file: Optional[Path] = None,
                                  width: int = 1280,
                                  height: int = 720,
                                  stride: int = 1,
                                  rotate_per_frame: float = 0.0,
                                  render_mode: str = "TachyonInternal") -> str:
        """
        Generate a self-contained VMD script that renders every frame to a TGA
        file and quits. Designed for headless use:

            vmd -dispdev text -e visualize_movie.vmd

        Args:
            trajectory_file:  .xtc/.trr to load
            pdb_file:         .pdb topology reference
            frames_dir:       directory (relative to workdir) for frame TGAs
            output_file:      optional path to write the script
            width, height:    render canvas size
            stride:           render every Nth frame (>=1)
            rotate_per_frame: degrees of y-axis rotation between frames
            render_mode:      VMD render method (TachyonInternal | snapshot)

        Returns:
            Script content as string.
        """
        stride = max(int(stride), 1)
        script = f'''# VMD Movie Render Script (headless)
# Run with: vmd -dispdev text -e visualize_movie.vmd

# Load structure and trajectory
mol new {pdb_file} type pdb waitfor all
mol addfile {trajectory_file} type xtc waitfor all

# Display / canvas
display projection Orthographic
display depthcue off
display resize {width} {height}
color Display Background white
axes location Off

# Cartoon rep coloured by chain
mol delrep 0 top
mol representation NewCartoon
mol color Chain
mol selection {{protein}}
mol addrep top

# Wrap and centre on protein
package require pbctools
pbc wrap -compound res -centersel "protein" -center com -all

# Frame output dir
file mkdir {frames_dir}

set num_frames [molinfo top get numframes]
puts "Rendering $num_frames frames at stride {stride} -> {frames_dir}/"

set out 0
for {{set i 0}} {{$i < $num_frames}} {{incr i {stride}}} {{
    animate goto $i
    if {{ {rotate_per_frame:.4f} != 0.0 }} {{
        rotate y by {rotate_per_frame:.4f}
    }}
    set fname [format "{frames_dir}/frame_%05d.tga" $out]
    render {render_mode} $fname
    incr out
}}

puts "Done. Wrote $out frames."
quit
'''

        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)

        return script

    def generate_vmd_trajectory(self, trajectory_file: str = "md_protein.xtc",
                                pdb_file: str = "final_structure.pdb",
                                output_file: Optional[Path] = None) -> str:
        """
        Generate VMD script for trajectory visualization.
        
        Args:
            trajectory_file: Name of trajectory file
            pdb_file: Name of PDB structure file
            output_file: Optional path to save script
            
        Returns:
            Script content as string
        """
        script = f'''# VMD Trajectory Visualization
# Load with: vmd -e visualize_trajectory.vmd

# Load structure and trajectory
mol new {pdb_file} type pdb
mol addfile {trajectory_file} type xtc waitfor all

# Display settings
display projection Orthographic
display depthcue off
color Display Background white
axes location Off

# Cartoon representation
mol delrep 0 top
mol representation NewCartoon
mol color Chain
mol selection {{protein}}
mol addrep top

# Align trajectory
package require pbctools
pbc wrap -compound res -centersel "protein" -center com

# Animation settings
animate style Loop
animate speed 0.5
animate goto 0

# To render movie (uncomment):
# for {{set i 0}} {{$i < [molinfo top get numframes]}} {{incr i}} {{
#     animate goto $i
#     render TachyonInternal frame_[format %04d $i].tga
# }}
'''
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        
        return script
    
    def generate_chimerax_interface(self, output_file: Optional[Path] = None) -> str:
        """
        Generate ChimeraX script for interface visualization.
        
        Args:
            output_file: Optional path to save script
            
        Returns:
            Script content as string
        """
        script = f'''# ChimeraX Interface Visualization Script
# Load with: chimerax visualize_interface.cxc

# Open structure
open {self.structure_file}

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
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        
        return script
    
    def generate_all_interface_scripts(self) -> Dict[str, Path]:
        """
        Generate all interface visualization scripts.
        
        Returns:
            Dictionary mapping script type to file path
        """
        viz_dir = self.workdir / 'visualization'
        viz_dir.mkdir(exist_ok=True)
        
        files = {}
        
        pymol_file = viz_dir / 'visualize_interface.pml'
        self.generate_pymol_interface(pymol_file)
        files['pymol'] = pymol_file
        
        vmd_file = viz_dir / 'visualize_interface.vmd'
        self.generate_vmd_interface(vmd_file)
        files['vmd'] = vmd_file
        
        chimerax_file = viz_dir / 'visualize_interface.cxc'
        self.generate_chimerax_interface(chimerax_file)
        files['chimerax'] = chimerax_file
        
        return files
    
    def generate_all_trajectory_scripts(self, 
                                        trajectory_file: str = "md_protein.xtc",
                                        pdb_file: str = "final_structure.pdb") -> Dict[str, Path]:
        """
        Generate all trajectory visualization scripts.
        
        Returns:
            Dictionary mapping script type to file path
        """
        viz_dir = self.workdir / 'visualization'
        viz_dir.mkdir(exist_ok=True)
        
        files = {}
        
        pymol_file = viz_dir / 'visualize_trajectory.pml'
        self.generate_pymol_trajectory(trajectory_file, pdb_file, pymol_file)
        files['pymol'] = pymol_file

        vmd_file = viz_dir / 'visualize_trajectory.vmd'
        self.generate_vmd_trajectory(trajectory_file, pdb_file, vmd_file)
        files['vmd'] = vmd_file

        # Headless render-and-quit script for video pipelines
        movie_file = viz_dir / 'visualize_movie.vmd'
        self.generate_vmd_movie_script(
            trajectory_file=trajectory_file,
            pdb_file=pdb_file,
            output_file=movie_file,
        )
        files['vmd_movie'] = movie_file

        return files


class PlotGenerator:
    """Generate plotting scripts for analysis data."""
    
    def __init__(self, workdir: Path):
        """
        Initialize plot generator.
        
        Args:
            workdir: Working directory
        """
        self.workdir = Path(workdir)
    
    def generate_gnuplot_rmsd(self, output_file: Optional[Path] = None) -> str:
        """Generate gnuplot script for RMSD plot."""
        script = '''# RMSD plot
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'rmsd_plot.png'

set title 'Backbone RMSD over Time'
set xlabel 'Time (ps)'
set ylabel 'RMSD (nm)'
set grid

plot 'rmsd.xvg' using 1:2 with lines lw 2 lc rgb 'blue' title 'RMSD'
'''
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        return script
    
    def generate_gnuplot_rmsf(self, output_file: Optional[Path] = None) -> str:
        """Generate gnuplot script for RMSF plot."""
        script = '''# RMSF plot
set terminal pngcairo size 1000,500 enhanced font 'Arial,12'
set output 'rmsf_plot.png'

set title 'Per-Residue RMSF (Flexibility)'
set xlabel 'Residue Number'
set ylabel 'RMSF (nm)'
set grid

plot 'rmsf.xvg' using 1:2 with lines lw 2 lc rgb 'red' title 'RMSF'
'''
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        return script
    
    def generate_gnuplot_gyration(self, output_file: Optional[Path] = None) -> str:
        """Generate gnuplot script for radius of gyration plot."""
        script = '''# Radius of Gyration plot
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'gyration_plot.png'

set title 'Radius of Gyration over Time'
set xlabel 'Time (ps)'
set ylabel 'Rg (nm)'
set grid

plot 'gyrate.xvg' using 1:2 with lines lw 2 lc rgb 'green' title 'Rg'
'''
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        return script
    
    def generate_gnuplot_contact_map(self, output_file: Optional[Path] = None) -> str:
        """Generate gnuplot script for contact map heatmap."""
        script = '''# Contact map heatmap
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
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        return script
    
    def generate_gnuplot_comparison(self, output_file: Optional[Path] = None) -> str:
        """Generate gnuplot script for batch comparison bar chart."""
        script = '''# Bar chart of potential energies
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
        if output_file:
            with open(output_file, 'w') as f:
                f.write(script)
        return script
    
    def generate_all_gnuplot_scripts(self) -> Dict[str, Path]:
        """Generate all gnuplot scripts."""
        plots_dir = self.workdir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        files = {}
        
        # Only generate scripts for files that exist
        if (self.workdir / 'rmsd.xvg').exists():
            files['rmsd'] = plots_dir / 'plot_rmsd.gp'
            self.generate_gnuplot_rmsd(files['rmsd'])
        
        if (self.workdir / 'rmsf.xvg').exists():
            files['rmsf'] = plots_dir / 'plot_rmsf.gp'
            self.generate_gnuplot_rmsf(files['rmsf'])
        
        if (self.workdir / 'gyrate.xvg').exists():
            files['gyration'] = plots_dir / 'plot_gyration.gp'
            self.generate_gnuplot_gyration(files['gyration'])
        
        if (self.workdir / 'contact_map_plot.csv').exists():
            files['contact_map'] = plots_dir / 'plot_contact_map.gp'
            self.generate_gnuplot_contact_map(files['contact_map'])
        
        return files


def generate_python_plotting_script(output_file: Path) -> None:
    """
    Generate Python script for comprehensive plotting.
    
    Args:
        output_file: Path to save the script
    """
    script = '''#!/usr/bin/env python3
"""Generate comparison plots from batch analysis."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

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
if 'gyration_nm' in df.columns and 'sasa_nm2' in df.columns:
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

print("Plots generated successfully!")
'''
    
    with open(output_file, 'w') as f:
        f.write(script)
    os.chmod(output_file, 0o755)


def generate_r_analysis_script(output_file: Path) -> None:
    """
    Generate R script for analysis.
    
    Args:
        output_file: Path to save the script
    """
    script = '''# R script for PPI comparison analysis
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

ggsave("comparison_barplot_r.png", p1, width=10, height=6)

# Scatter plot: Rg vs SASA
if("gyration_nm" %in% names(data) && "sasa_nm2" %in% names(data)) {
  p2 <- ggplot(data, aes(x=gyration_nm, y=sasa_nm2, label=structure)) +
    geom_point(size=3, color="red") +
    geom_text(vjust=-0.5, size=3) +
    labs(title="Compactness vs Surface Area", x="Radius of Gyration (nm)", y="SASA (nm²)") +
    theme_minimal()
  
  ggsave("rg_vs_sasa_r.png", p2, width=8, height=6)
}

print("Plots saved!")
'''
    
    with open(output_file, 'w') as f:
        f.write(script)


def main():
    """CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate visualization scripts")
    parser.add_argument("--workdir", required=True, help="Working directory")
    parser.add_argument("--type", choices=['interface', 'trajectory', 'plots', 'movie', 'all'],
                        default='all', help="Type of visualization to generate")
    parser.add_argument("--trajectory-file", default="trajectories/md_center.xtc",
                        help="Trajectory file for movie/trajectory scripts (relative to workdir)")
    parser.add_argument("--pdb-file", default="structures/final_structure.pdb",
                        help="PDB topology for movie/trajectory scripts (relative to workdir)")
    parser.add_argument("--frames-dir", default="video/frames",
                        help="Where the movie script writes frame TGAs")
    parser.add_argument("--width", type=int, default=1280)
    parser.add_argument("--height", type=int, default=720)
    parser.add_argument("--stride", type=int, default=1,
                        help="Render every Nth trajectory frame")
    parser.add_argument("--rotate-per-frame", type=float, default=0.0,
                        help="Y-axis rotation between frames (degrees)")
    parser.add_argument("--render-mode", default="TachyonInternal",
                        help="VMD render method (TachyonInternal, snapshot, etc.)")
    args = parser.parse_args()
    
    workdir = Path(args.workdir)
    
    if args.type in ['interface', 'all']:
        viz_gen = VisualizationGenerator(workdir)
        files = viz_gen.generate_all_interface_scripts()
        print("Generated interface scripts:")
        for name, path in files.items():
            print(f"  {name}: {path}")
    
    if args.type in ['trajectory', 'all']:
        viz_gen = VisualizationGenerator(workdir)
        files = viz_gen.generate_all_trajectory_scripts()
        print("Generated trajectory scripts:")
        for name, path in files.items():
            print(f"  {name}: {path}")
    
    if args.type in ['plots', 'all']:
        plot_gen = PlotGenerator(workdir)
        files = plot_gen.generate_all_gnuplot_scripts()
        print("Generated gnuplot scripts:")
        for name, path in files.items():
            print(f"  {name}: {path}")

    if args.type in ['movie', 'all']:
        viz_gen = VisualizationGenerator(workdir)
        viz_dir = workdir / 'visualization'
        viz_dir.mkdir(exist_ok=True)
        movie_file = viz_dir / 'visualize_movie.vmd'
        viz_gen.generate_vmd_movie_script(
            trajectory_file=args.trajectory_file,
            pdb_file=args.pdb_file,
            frames_dir=args.frames_dir,
            output_file=movie_file,
            width=args.width,
            height=args.height,
            stride=args.stride,
            rotate_per_frame=args.rotate_per_frame,
            render_mode=args.render_mode,
        )
        print(f"Generated movie script: {movie_file}")


if __name__ == "__main__":
    main()
