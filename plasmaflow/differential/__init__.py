"""
Differential analysis module for PlasmaFlow

This module provides functionality for differential chromatin loop analysis
using multiple statistical methods including diffHic, DESeq2, and edgeR through R integration.
"""

from .analyzer import DifferentialAnalyzer, DifferentialResult
from .comparison import ComparisonManager, create_comparison_matrix
from .filtering import (FilterManager, apply_count_filters,
                        apply_distance_filters)
from .methods import DESeq2Analyzer, DiffHicAnalyzer, EdgeRAnalyzer
from .results import ResultsProcessor, export_results, merge_results
from .visualization import (create_comparison_plots, create_ma_plot,
                            create_volcano_plot)

__all__ = [
    "DifferentialAnalyzer",
    "DifferentialResult",
    "DiffHicAnalyzer",
    "DESeq2Analyzer",
    "EdgeRAnalyzer",
    "ComparisonManager",
    "create_comparison_matrix",
    "create_volcano_plot",
    "create_ma_plot",
    "create_comparison_plots",
    "FilterManager",
    "apply_distance_filters",
    "apply_count_filters",
    "ResultsProcessor",
    "merge_results",
    "export_results",
]
