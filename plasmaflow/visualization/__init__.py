"""
Visualization module for PlasmaFlow

This module provides functionality for creating visualizations and matrices
from differential analysis results, including heatmaps, volcano plots,
and deepTools matrix generation.
"""

from .bed_utils import create_anchor_beds, parse_loop_coordinates
from .deeptools import DeepToolsInterface, validate_bigwig_files
from .heatmaps import (HeatmapGenerator, create_clustered_heatmap,
                       create_heatmap)
from .matrix_generator import (MatrixGenerator, create_bed_files,
                               generate_deeptools_matrix)
from .plots import create_comparison_plots, create_summary_plots
from .volcano import create_enhanced_volcano_plot, create_volcano_plot

__all__ = [
    "MatrixGenerator",
    "create_bed_files",
    "generate_deeptools_matrix",
    "HeatmapGenerator",
    "create_heatmap",
    "create_clustered_heatmap",
    "create_volcano_plot",
    "create_enhanced_volcano_plot",
    "create_comparison_plots",
    "create_summary_plots",
    "parse_loop_coordinates",
    "create_anchor_beds",
    "DeepToolsInterface",
    "validate_bigwig_files",
]
