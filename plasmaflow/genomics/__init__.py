"""
Genomics analysis module for PlasmaFlow

This module provides functionality for:
- Gene proximity analysis for differential chromatin loops
- TSS-based distance calculations
- Regulatory region classification
- Cluster and category statistical comparisons
- Gene annotation and filtering
- Differential gene expression analysis
- Volcano plot visualization with key gene labeling
"""

from .annotations import GeneAnnotationManager
from .expression import GeneExpressionAnalyzer, analyze_expression_overlap
from .proximity import GeneProximityAnalyzer, ProximityResult
from .statistical import GeneDistanceStatistics

__all__ = [
    "GeneProximityAnalyzer",
    "ProximityResult",
    "GeneAnnotationManager",
    "GeneDistanceStatistics",
    "GeneExpressionAnalyzer",
    "analyze_expression_overlap",
]
