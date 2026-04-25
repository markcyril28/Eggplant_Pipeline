#!/usr/bin/env python3
"""
Output Organizer Module for GROMACS PPI Analysis

Manages standardized output folder structure and file organization.
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class OutputStructure:
    """Defines the output directory structure."""
    
    # Main subdirectories
    structures: str = "structures"       # Prepared structure files
    analysis: str = "analysis"           # Analysis results (metrics, contacts)
    visualization: str = "visualization" # Visualization scripts
    statistics: str = "statistics"       # JSON/CSV statistics
    plots: str = "plots"                 # Generated plots
    logs: str = "logs"                   # Log files
    trajectories: str = "trajectories"   # MD trajectory files
    
    # File patterns for organization
    structure_patterns: tuple = ('.gro', '.pdb', '.tpr')
    log_patterns: tuple = ('.log',)
    trajectory_patterns: tuple = ('.xtc', '.trr', '.edr')
    analysis_patterns: tuple = ('.xvg', '.xpm')
    visualization_patterns: tuple = ('.pml', '.vmd', '.cxc', '.tcl')
    plot_patterns: tuple = ('.png', '.svg', '.pdf', '.gp')
    statistics_patterns: tuple = ('.json', '.csv')


class OutputOrganizer:
    """Organizes GROMACS output files into standardized structure."""
    
    def __init__(self, base_dir: Path, structure: Optional[OutputStructure] = None):
        """
        Initialize output organizer.
        
        Args:
            base_dir: Base output directory
            structure: Output structure definition (uses default if None)
        """
        self.base_dir = Path(base_dir)
        self.structure = structure or OutputStructure()
        self._dirs = {}
    
    def create_structure(self) -> Dict[str, Path]:
        """
        Create all output directories.
        
        Returns:
            Dictionary mapping directory names to paths
        """
        self._dirs = {}
        
        for attr in ['structures', 'analysis', 'visualization', 
                     'statistics', 'plots', 'logs', 'trajectories']:
            dir_name = getattr(self.structure, attr)
            path = self.base_dir / dir_name
            path.mkdir(parents=True, exist_ok=True)
            self._dirs[attr] = path
        
        return self._dirs
    
    def get_dir(self, name: str) -> Path:
        """Get path for a specific directory."""
        if not self._dirs:
            self.create_structure()
        return self._dirs.get(name, self.base_dir / name)
    
    def organize_files(self, source_dir: Optional[Path] = None) -> Dict[str, List[Path]]:
        """
        Organize files from source directory into proper structure.
        
        Args:
            source_dir: Source directory (uses base_dir if None)
            
        Returns:
            Dictionary mapping categories to lists of organized files
        """
        if not self._dirs:
            self.create_structure()
        
        source = Path(source_dir) if source_dir else self.base_dir
        organized = {cat: [] for cat in self._dirs.keys()}
        
        for file_path in source.iterdir():
            if file_path.is_file():
                dest = self._categorize_file(file_path)
                if dest:
                    category, dest_dir = dest
                    dest_path = dest_dir / file_path.name
                    if file_path != dest_path:
                        shutil.copy2(file_path, dest_path)
                        organized[category].append(dest_path)
        
        return organized
    
    def _categorize_file(self, file_path: Path) -> Optional[tuple]:
        """
        Determine which category a file belongs to.
        
        Args:
            file_path: Path to file
            
        Returns:
            Tuple of (category_name, destination_path) or None
        """
        suffix = file_path.suffix.lower()
        
        if suffix in self.structure.log_patterns:
            return ('logs', self._dirs['logs'])
        elif suffix in self.structure.trajectory_patterns:
            return ('trajectories', self._dirs['trajectories'])
        elif suffix in self.structure.visualization_patterns:
            return ('visualization', self._dirs['visualization'])
        elif suffix in self.structure.plot_patterns:
            return ('plots', self._dirs['plots'])
        elif suffix in self.structure.statistics_patterns:
            # Distinguish between analysis stats and comparison stats
            if 'comparison' in file_path.name or 'batch' in file_path.name:
                return ('statistics', self._dirs['statistics'])
            return ('statistics', self._dirs['statistics'])
        elif suffix in self.structure.analysis_patterns:
            return ('analysis', self._dirs['analysis'])
        elif suffix in self.structure.structure_patterns:
            return ('structures', self._dirs['structures'])
        
        return None
    
    def collect_logs(self, source_dir: Path, prefix: str = "") -> List[Path]:
        """
        Collect log files from source directory to logs folder.
        
        Args:
            source_dir: Source directory
            prefix: Optional prefix for log files
            
        Returns:
            List of collected log files
        """
        if not self._dirs:
            self.create_structure()
        
        logs_dir = self._dirs['logs']
        collected = []
        
        for log_file in source_dir.glob('*.log'):
            if prefix:
                dest_name = f"{prefix}_{log_file.name}"
            else:
                dest_name = log_file.name
            
            dest_path = logs_dir / dest_name
            shutil.copy2(log_file, dest_path)
            collected.append(dest_path)
        
        return collected
    
    def create_summary(self) -> Path:
        """
        Create a summary of the output structure.
        
        Returns:
            Path to summary file
        """
        if not self._dirs:
            self.create_structure()
        
        summary_file = self.base_dir / 'OUTPUT_STRUCTURE.txt'
        
        with open(summary_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("OUTPUT DIRECTORY STRUCTURE\n")
            f.write("=" * 60 + "\n\n")
            
            for name, path in sorted(self._dirs.items()):
                files = list(path.glob('*'))
                file_count = len([f for f in files if f.is_file()])
                f.write(f"{name}/ ({file_count} files)\n")
                
                # List a few files
                for file_path in sorted(files)[:5]:
                    if file_path.is_file():
                        f.write(f"  - {file_path.name}\n")
                if file_count > 5:
                    f.write(f"  ... and {file_count - 5} more files\n")
                f.write("\n")
        
        return summary_file


def setup_output_structure(base_dir: Path, 
                          structure_name: str = "") -> OutputOrganizer:
    """
    Set up output structure for an analysis.
    
    Args:
        base_dir: Base output directory
        structure_name: Optional structure name for subdirectory
        
    Returns:
        Configured OutputOrganizer
    """
    if structure_name:
        work_dir = base_dir / structure_name
    else:
        work_dir = base_dir
    
    organizer = OutputOrganizer(work_dir)
    organizer.create_structure()
    
    return organizer


def get_standard_paths(base_dir: Path) -> Dict[str, Path]:
    """
    Get dictionary of standard output paths.
    
    Args:
        base_dir: Base output directory
        
    Returns:
        Dictionary mapping category names to paths
    """
    structure = OutputStructure()
    return {
        'structures': base_dir / structure.structures,
        'analysis': base_dir / structure.analysis,
        'visualization': base_dir / structure.visualization,
        'statistics': base_dir / structure.statistics,
        'plots': base_dir / structure.plots,
        'logs': base_dir / structure.logs,
        'trajectories': base_dir / structure.trajectories
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Organize GROMACS output files")
    parser.add_argument("--base-dir", required=True, help="Base output directory")
    parser.add_argument("--source", help="Source directory to organize (if different)")
    parser.add_argument("--create", action="store_true", help="Create directory structure")
    parser.add_argument("--organize", action="store_true", help="Organize existing files")
    parser.add_argument("--summary", action="store_true", help="Create summary file")
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    organizer = OutputOrganizer(base_dir)
    
    if args.create:
        dirs = organizer.create_structure()
        print("Created directories:")
        for name, path in dirs.items():
            print(f"  {name}: {path}")
    
    if args.organize:
        source = Path(args.source) if args.source else None
        organized = organizer.organize_files(source)
        print("\nOrganized files:")
        for category, files in organized.items():
            if files:
                print(f"  {category}: {len(files)} files")
    
    if args.summary:
        summary = organizer.create_summary()
        print(f"\nSummary: {summary}")
