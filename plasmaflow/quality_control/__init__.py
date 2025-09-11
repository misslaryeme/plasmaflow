"""
Quality control module for PlasmaFlow

This module provides functionality for quality control analysis of chromatin loops,
including Aggregate Peak Analysis (APA) and comparative APA analysis across samples.
"""

from .apa import APAAnalyzer, APAResult, calculate_apa_scores
from .caching import APACache, load_cached_apa, save_apa_cache
from .comparative import APAComparison, ComparativeAPAAnalyzer
from .metrics import QualityMetrics, calculate_quality_metrics
from .visualization import APAPlotter, create_apa_plots

__all__ = [
    "APAAnalyzer",
    "APAResult",
    "calculate_apa_scores",
    "ComparativeAPAAnalyzer",
    "APAComparison",
    "QualityMetrics",
    "calculate_quality_metrics",
    "APAPlotter",
    "create_apa_plots",
    "APACache",
    "load_cached_apa",
    "save_apa_cache",
]
