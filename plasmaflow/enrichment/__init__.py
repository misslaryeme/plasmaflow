"""
Pathway enrichment analysis module for PlasmaFlow

This module provides comprehensive pathway enrichment functionality including:
- Gene Ontology (GO) enrichment analysis with GSEA and ORA methods
- KEGG pathway analysis
- MSigDB C2 curated gene set analysis for B-cell specific pathways
- Publication-ready visualization with organized plot output
- R integration via rpy2 for clusterProfiler analysis
- Statistical testing with proper multiple testing correction
"""

from .clusterprofiler import ClusterProfilerAnalyzer, PathwayEnrichmentResult
from .msigdb import MSigDBAnalyzer, prepare_bcell_msigdb_gene_sets
from .r_interface import RPathwayInterface, validate_r_environment
from .visualization import PathwayPlotter, create_comprehensive_pathway_plots

__all__ = [
    "ClusterProfilerAnalyzer",
    "PathwayEnrichmentResult",
    "MSigDBAnalyzer",
    "prepare_bcell_msigdb_gene_sets",
    "PathwayPlotter",
    "create_comprehensive_pathway_plots",
    "RPathwayInterface",
    "validate_r_environment",
]
