#!/usr/bin/env python3
"""
Configuration Module for GROMACS PPI Analysis

Centralizes all configuration settings, paths, and constants.
"""

import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Dict, List

# Base paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
MODULES_DIR = PROJECT_ROOT / "modules" / "gromacs_utils"

# Default GROMACS settings
DEFAULT_FORCEFIELD = "amber99sb-ildn"
DEFAULT_WATERMODEL = "tip3p"
DEFAULT_BOX_DISTANCE = 1.0  # nm
DEFAULT_EM_STEPS = 5000
DEFAULT_NTHREADS = 8

# GPU settings for NVIDIA CUDA
# EM: steep integrator cannot use -pme gpu or -update gpu; -bonded gpu OK on NVIDIA
GPU_EM_FLAGS = "-nb gpu -pme cpu -bonded gpu"
GPU_MD_FLAGS = "-nb gpu -pme gpu -bonded gpu -update gpu"

# Analysis parameters
INTERFACE_CUTOFF = 0.5   # nm - residues within this distance are "interface"
CONTACT_CUTOFF = 0.6     # nm - for contact counting  
HBOND_CUTOFF = 0.35      # nm - H-bond distance cutoff

# Visualization settings
VIZ_FORMATS = ['pymol', 'vmd', 'chimerax']
PLOT_DPI = 150
MOVIE_FPS = 15


@dataclass
class GromacsConfig:
    """GROMACS simulation configuration."""
    
    # Paths
    gmx_bin: str = field(default_factory=lambda: str(Path.home() / "miniconda3/envs/gromacs_CUDA/bin/gmx"))
    workdir: Path = field(default_factory=lambda: Path.cwd())
    
    # Force field
    forcefield: str = DEFAULT_FORCEFIELD
    watermodel: str = DEFAULT_WATERMODEL
    
    # Box and solvation
    box_distance: float = DEFAULT_BOX_DISTANCE
    box_type: str = "dodecahedron"
    ion_concentration: float = 0.15  # M NaCl
    
    # Energy minimization
    em_steps: int = DEFAULT_EM_STEPS
    em_tolerance: float = 100.0  # kJ/mol/nm
    
    # MD parameters
    nvt_steps: int = 50000   # 100 ps with 2fs timestep
    npt_steps: int = 50000   # 100 ps
    md_steps: int = 250000   # 500 ps
    timestep: float = 0.002  # ps (2 fs)
    temperature: float = 300  # K
    pressure: float = 1.0    # bar
    
    # Performance
    nthreads: int = DEFAULT_NTHREADS
    gpu_id: int = 0
    use_gpu: bool = True
    
    # GPU flags
    gpu_em_flags: str = GPU_EM_FLAGS
    gpu_md_flags: str = GPU_MD_FLAGS
    
    def get_duration_ps(self, phase: str) -> float:
        """Get simulation duration in picoseconds."""
        steps = getattr(self, f"{phase}_steps", 0)
        return steps * self.timestep


@dataclass 
class AnalysisConfig:
    """Analysis configuration."""
    
    interface_cutoff: float = INTERFACE_CUTOFF
    contact_cutoff: float = CONTACT_CUTOFF
    hbond_cutoff: float = HBOND_CUTOFF
    
    # Binding quality thresholds
    bsa_good: float = 8.0       # nm² - good interface
    bsa_moderate: float = 4.0   # nm² - moderate interface
    hbond_good: int = 10        # good H-bond network
    hbond_moderate: int = 5     # moderate H-bond network
    contact_good: int = 50      # extensive contacts
    contact_moderate: int = 20  # moderate contacts


@dataclass
class OutputConfig:
    """Output folder configuration."""
    
    base_dir: Path = field(default_factory=lambda: Path.cwd())
    
    # Subdirectory names
    structures_dir: str = "structures"
    analysis_dir: str = "analysis"
    visualization_dir: str = "visualization"
    statistics_dir: str = "statistics"
    plots_dir: str = "plots"
    logs_dir: str = "logs"
    trajectories_dir: str = "trajectories"
    
    def create_dirs(self) -> Dict[str, Path]:
        """Create all output directories and return paths."""
        dirs = {}
        for name in ['structures', 'analysis', 'visualization', 
                     'statistics', 'plots', 'logs', 'trajectories']:
            dir_name = getattr(self, f"{name}_dir")
            path = self.base_dir / dir_name
            path.mkdir(parents=True, exist_ok=True)
            dirs[name] = path
        return dirs
    
    def get_path(self, category: str) -> Path:
        """Get path for a specific category."""
        dir_name = getattr(self, f"{category}_dir", category)
        return self.base_dir / dir_name


def get_default_config() -> GromacsConfig:
    """Get default GROMACS configuration."""
    return GromacsConfig()


def get_analysis_config() -> AnalysisConfig:
    """Get default analysis configuration."""
    return AnalysisConfig()


def get_output_config(base_dir: Path) -> OutputConfig:
    """Get output configuration for a given base directory."""
    return OutputConfig(base_dir=base_dir)


# Environment setup for NVIDIA CUDA GPUs
NVIDIA_GPU_ENV = {
    'GMX_ENABLE_DIRECT_GPU_COMM': '1',
    'CUDA_VISIBLE_DEVICES': '0',
}


def setup_nvidia_gpu_env():
    """Set environment variables for NVIDIA GPU."""
    for key, value in NVIDIA_GPU_ENV.items():
        os.environ[key] = value


if __name__ == "__main__":
    # Print configuration info
    config = get_default_config()
    print(f"Project root: {PROJECT_ROOT}")
    print(f"GMX binary: {config.gmx_bin}")
    print(f"Force field: {config.forcefield}")
    print(f"GPU EM flags: {config.gpu_em_flags}")
    print(f"GPU MD flags: {config.gpu_md_flags}")
