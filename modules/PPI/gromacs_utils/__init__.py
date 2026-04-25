#!/usr/bin/env python3
# GROMACS Utility Modules for PPI Analysis
"""
This package contains utility functions for GROMACS-based
protein-protein interaction (PPI) analysis.

Modules:
- config: Configuration settings and paths
- xvg_parser: Parse GROMACS XVG output files
- chain_parser: Parse chain information from PDB/GRO files
- chain_indexer: Create GROMACS index groups for chains
- interface_analyzer: Analyze PPI interfaces (BSA, H-bonds, contacts)
- contact_map: Generate inter-chain contact maps
- md_statistics: Extract MD simulation statistics
- visualization_generator: Generate PyMOL/VMD/ChimeraX scripts
- batch_analyzer: Batch structure comparison
- output_organizer: Organize output files
- metrics_extractor: Extract metrics from GROMACS outputs
- stability_analyzer: Analyze structure stability
- results_comparator: Compare quick stability results
- md_results_comparator: Compare MD simulation results
"""

# Configuration
from .config import (
    GromacsConfig, 
    AnalysisConfig, 
    OutputConfig,
    get_default_config,
    get_analysis_config,
    get_output_config,
    setup_nvidia_gpu_env
)

# XVG parsing
from .xvg_parser import (
    parse_xvg,
    parse_xvg_last,
    parse_xvg_column,
    calculate_statistics,
    xvg_to_csv
)

# Chain handling
from .chain_parser import get_chain_info, parse_chains_for_index
from .chain_indexer import create_chain_index, parse_gro_for_chains

# Interface analysis
from .interface_analyzer import (
    InterfaceMetrics,
    extract_interface_metrics,
    assess_binding_quality,
    generate_interface_report,
    save_metrics_json
)

# Contact map
from .contact_map import (
    generate_contact_map,
    get_residue_coords,
    get_chain_residues,
    calculate_contact_map
)

# MD statistics
from .md_statistics import (
    MDStatistics,
    extract_md_statistics,
    generate_md_statistics,
    save_statistics_json as save_md_statistics_json,
    save_statistics_csv as save_md_statistics_csv
)

# Visualization
from .visualization_generator import (
    VisualizationGenerator,
    PlotGenerator,
    generate_python_plotting_script,
    generate_r_analysis_script
)

# Batch analysis
from .batch_analyzer import (
    StructureMetrics,
    collect_structure_metrics,
    analyze_all_structures,
    generate_comparison_report,
    generate_plotting_scripts,
    run_batch_analysis
)

# Output organization
from .output_organizer import (
    OutputOrganizer,
    OutputStructure,
    setup_output_structure,
    get_standard_paths
)

# Structure conversion
from .structure_converter import (
    convert_cif_to_pdb,
    clean_pdb_for_gromacs,
    prepare_structure,
    get_structure_info
)

# MDP generation
from .mdp_generator import (
    MDPConfig,
    generate_em_mdp,
    generate_nvt_mdp,
    generate_npt_mdp,
    generate_md_mdp,
    generate_ie_mdp,
    generate_all_mdp_files
)

# Plotting
from .plotting import (
    check_matplotlib,
    plot_comparison_barplot,
    plot_scatter_comparison,
    plot_timeseries,
    plot_rmsd,
    plot_rmsf,
    plot_gyration,
    plot_hbonds,
    plot_contact_map,
    plot_all_md_analysis,
    generate_batch_plots
)

# Legacy modules
from .metrics_extractor import extract_metrics
from .stability_analyzer import analyze_stability, read_xvg
from .results_comparator import compare_stability_results
from .md_results_comparator import compare_md_results

__all__ = [
    # Configuration
    'GromacsConfig',
    'AnalysisConfig', 
    'OutputConfig',
    'get_default_config',
    'get_analysis_config',
    'get_output_config',
    'setup_nvidia_gpu_env',
    
    # XVG parsing
    'parse_xvg',
    'parse_xvg_last',
    'parse_xvg_column',
    'calculate_statistics',
    'xvg_to_csv',
    
    # Chain handling
    'get_chain_info',
    'parse_chains_for_index',
    'create_chain_index',
    'parse_gro_for_chains',
    
    # Interface analysis
    'InterfaceMetrics',
    'extract_interface_metrics',
    'assess_binding_quality',
    'generate_interface_report',
    'save_metrics_json',
    
    # Contact map
    'generate_contact_map',
    'get_residue_coords',
    'get_chain_residues',
    'calculate_contact_map',
    
    # MD statistics
    'MDStatistics',
    'extract_md_statistics',
    'generate_md_statistics',
    'save_md_statistics_json',
    'save_md_statistics_csv',
    
    # Visualization
    'VisualizationGenerator',
    'PlotGenerator',
    'generate_python_plotting_script',
    'generate_r_analysis_script',
    
    # Batch analysis
    'StructureMetrics',
    'collect_structure_metrics',
    'analyze_all_structures',
    'generate_comparison_report',
    'generate_plotting_scripts',
    'run_batch_analysis',
    
    # Output organization
    'OutputOrganizer',
    'OutputStructure',
    'setup_output_structure',
    'get_standard_paths',
    
    # Structure conversion
    'convert_cif_to_pdb',
    'clean_pdb_for_gromacs',
    'prepare_structure',
    'get_structure_info',
    
    # MDP generation
    'MDPConfig',
    'generate_em_mdp',
    'generate_nvt_mdp',
    'generate_npt_mdp',
    'generate_md_mdp',
    'generate_ie_mdp',
    'generate_all_mdp_files',
    
    # Plotting
    'check_matplotlib',
    'plot_comparison_barplot',
    'plot_scatter_comparison',
    'plot_timeseries',
    'plot_rmsd',
    'plot_rmsf',
    'plot_gyration',
    'plot_hbonds',
    'plot_contact_map',
    'plot_all_md_analysis',
    'generate_batch_plots',
    
    # Legacy
    'extract_metrics',
    'analyze_stability',
    'read_xvg',
    'compare_stability_results',
    'compare_md_results'
]
