"""
MSigDB integration for specialized gene set analysis

This module provides access to Molecular Signatures Database (MSigDB) gene sets,
with particular focus on B-cell and immunological pathways relevant to plasma cell differentiation.
"""

import gzip
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from urllib.request import urlretrieve

import numpy as np
import pandas as pd

from ..core.config import Config
from ..core.exceptions import PlasmaFlowError

logger = logging.getLogger(__name__)


class MSigDBAnalyzer:
    """Analysis interface for MSigDB gene sets"""

    def __init__(self, config: Config):
        self.config = config
        self.cache_dir = Path(config.analysis_config.output_dir) / "msigdb_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # MSigDB collections relevant to B-cell biology
        self.bcell_collections = {
            "C2_CP_REACTOME": "Reactome canonical pathways",
            "C2_CP_KEGG": "KEGG canonical pathways",
            "C2_CGP": "Chemical and genetic perturbations",
            "C7_IMMUNESIGDB": "Immunological signatures",
            "C8_CELL_TYPE": "Cell type signatures",
        }

        # B-cell specific keywords for filtering
        self.bcell_keywords = [
            "B_CELL",
            "BCELL",
            "PLASMA",
            "PLASMABLAST",
            "MEMORY_B",
            "GERMINAL_CENTER",
            "IMMUNOGLOBULIN",
            "ANTIBODY",
            "BCR",
            "CD19",
            "CD20",
            "CD27",
            "CD38",
            "CD138",
            "PAX5",
            "BLIMP1",
            "XBP1",
            "IRF4",
            "PRDM1",
        ]

    def prepare_bcell_msigdb_gene_sets(self) -> Dict[str, Set[str]]:
        """
        Prepare B-cell relevant MSigDB gene sets

        Returns:
        --------
        Dict[str, Set[str]]
            Dictionary mapping pathway names to gene sets
        """
        logger.info("Preparing B-cell specific MSigDB gene sets")

        bcell_gene_sets = {}

        # Load and filter gene sets from each collection
        for collection, description in self.bcell_collections.items():
            try:
                gene_sets = self._load_msigdb_collection(collection)
                filtered_sets = self._filter_bcell_relevant_sets(gene_sets)
                bcell_gene_sets.update(filtered_sets)

                logger.info(
                    f"Loaded {len(filtered_sets)} B-cell relevant sets from {collection}"
                )

            except Exception as e:
                logger.warning(f"Could not load {collection}: {str(e)}")
                continue

        # Add custom B-cell differentiation gene sets
        custom_sets = self._get_custom_bcell_gene_sets()
        bcell_gene_sets.update(custom_sets)

        logger.info(f"Total B-cell gene sets prepared: {len(bcell_gene_sets)}")

        # Save gene sets for R analysis
        self._save_gene_sets_for_r(bcell_gene_sets)

        return bcell_gene_sets

    def _load_msigdb_collection(self, collection: str) -> Dict[str, Set[str]]:
        """Load MSigDB gene sets from a specific collection"""
        cache_file = self.cache_dir / f"{collection}.json"

        # Check cache first
        if cache_file.exists():
            logger.debug(f"Loading {collection} from cache")
            with open(cache_file, "r") as f:
                cached_data = json.load(f)
                return {k: set(v) for k, v in cached_data.items()}

        # For demonstration, return hardcoded B-cell relevant pathways
        # In production, this would interface with MSigDB API or files
        gene_sets = self._get_demo_msigdb_sets(collection)

        # Cache the results
        cache_data = {k: list(v) for k, v in gene_sets.items()}
        with open(cache_file, "w") as f:
            json.dump(cache_data, f, indent=2)

        return gene_sets

    def _get_demo_msigdb_sets(self, collection: str) -> Dict[str, Set[str]]:
        """Get demonstration MSigDB gene sets relevant to B-cell biology"""

        if collection == "C2_CP_REACTOME":
            return {
                "REACTOME_B_CELL_RECEPTOR_SIGNALING_PATHWAY": {
                    "CD19",
                    "CD79A",
                    "CD79B",
                    "BTK",
                    "SYK",
                    "LYN",
                    "FYN",
                    "BLNK",
                    "PLCG2",
                    "PIK3CG",
                    "AKT1",
                    "NFKB1",
                    "NFKB2",
                },
                "REACTOME_IMMUNOGLOBULIN_PRODUCTION": {
                    "IGHM",
                    "IGHD",
                    "IGHG1",
                    "IGHG2",
                    "IGHG3",
                    "IGHG4",
                    "IGHA1",
                    "IGHA2",
                    "IGHE",
                    "IGLC1",
                    "IGLC2",
                    "IGKC",
                },
                "REACTOME_CLASS_SWITCHING": {
                    "AID",
                    "UNG",
                    "MSH2",
                    "MSH6",
                    "PCNA",
                    "RFC1",
                    "POLD1",
                    "POLD2",
                    "POLE",
                    "LIG1",
                    "XRCC1",
                },
                "REACTOME_PLASMA_CELL_DIFFERENTIATION": {
                    "PRDM1",
                    "XBP1",
                    "IRF4",
                    "PAX5",
                    "BCL6",
                    "BACH2",
                    "CD27",
                    "CD38",
                    "CD138",
                    "TNFRSF17",
                    "IL6R",
                },
            }

        elif collection == "C2_CP_KEGG":
            return {
                "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY": {
                    "CD19",
                    "CD22",
                    "CD79A",
                    "CD79B",
                    "BTK",
                    "SYK",
                    "LYN",
                    "VAV1",
                    "VAV3",
                    "PLCG2",
                    "PRKCB",
                    "MAP2K1",
                    "MAPK1",
                },
                "KEGG_PRIMARY_IMMUNODEFICIENCY": {
                    "BTK",
                    "IGHM",
                    "CD79A",
                    "CD79B",
                    "BLNK",
                    "LRRC8A",
                    "ADA",
                    "PNP",
                    "RAG1",
                    "RAG2",
                    "DCLRE1C",
                },
            }

        elif collection == "C7_IMMUNESIGDB":
            return {
                "PLASMA_CELLS_VS_B_CELLS_UP": {
                    "PRDM1",
                    "XBP1",
                    "IRF4",
                    "CD27",
                    "CD38",
                    "CD138",
                    "TNFRSF17",
                    "IGHG1",
                    "IGHG3",
                    "JCHAIN",
                    "MZB1",
                },
                "MEMORY_B_CELLS_VS_NAIVE_B_CELLS_UP": {
                    "CD27",
                    "IGHG1",
                    "IGHG3",
                    "IGHD",
                    "FCRL4",
                    "FCRL5",
                    "ITGAX",
                    "TBX21",
                    "FCGR2B",
                    "CD86",
                    "CD95",
                },
                "GERMINAL_CENTER_B_CELLS_UP": {
                    "BCL6",
                    "AICDA",
                    "CD95",
                    "CD83",
                    "CXCR4",
                    "CXCR5",
                    "CD86",
                    "ICOS",
                    "SLAMF1",
                    "MEF2B",
                    "EZH2",
                },
            }

        elif collection == "C8_CELL_TYPE":
            return {
                "NAIVE_B_CELL_SIGNATURE": {
                    "CD19",
                    "CD20",
                    "PAX5",
                    "MS4A1",
                    "CR2",
                    "FCER2",
                    "IL4R",
                    "CXCR5",
                    "CCR7",
                    "SELL",
                    "TCL1A",
                },
                "PLASMA_CELL_SIGNATURE": {
                    "CD138",
                    "PRDM1",
                    "XBP1",
                    "IRF4",
                    "TNFRSF17",
                    "JCHAIN",
                    "MZB1",
                    "DERL3",
                    "FKBP11",
                    "SSR4",
                },
            }

        else:
            return {}

    def _filter_bcell_relevant_sets(
        self, gene_sets: Dict[str, Set[str]]
    ) -> Dict[str, Set[str]]:
        """Filter gene sets for B-cell relevance based on keywords"""
        filtered_sets = {}

        for pathway_name, gene_set in gene_sets.items():
            # Check if pathway name contains B-cell keywords
            pathway_upper = pathway_name.upper()
            if any(keyword in pathway_upper for keyword in self.bcell_keywords):
                filtered_sets[pathway_name] = gene_set

            # Also check if gene set contains key B-cell genes
            elif self._contains_bcell_genes(gene_set):
                filtered_sets[pathway_name] = gene_set

        return filtered_sets

    def _contains_bcell_genes(self, gene_set: Set[str]) -> bool:
        """Check if gene set contains key B-cell marker genes"""
        bcell_markers = {
            "CD19",
            "CD20",
            "PAX5",
            "PRDM1",
            "XBP1",
            "IRF4",
            "CD27",
            "CD38",
            "CD138",
            "BCL6",
            "AICDA",
        }

        # Require at least 2 B-cell markers
        overlap = len(gene_set.intersection(bcell_markers))
        return overlap >= 2

    def _get_custom_bcell_gene_sets(self) -> Dict[str, Set[str]]:
        """Get custom B-cell differentiation gene sets"""
        return {
            "PLASMA_CELL_DIFFERENTIATION_CORE": {
                "PRDM1",
                "XBP1",
                "IRF4",
                "CD27",
                "CD38",
                "CD138",
                "TNFRSF17",
                "JCHAIN",
                "MZB1",
                "DERL3",
                "FKBP11",
                "SSR4",
                "CALR",
                "CANX",
                "HSPA5",
                "DDIT3",
            },
            "B_CELL_ACTIVATION_SIGNATURE": {
                "CD19",
                "CD79A",
                "CD79B",
                "BTK",
                "SYK",
                "LYN",
                "BLNK",
                "PLCG2",
                "NFKB1",
                "REL",
                "NFATC1",
                "CD86",
                "CD80",
                "ICOS",
                "CD40",
                "IL4R",
            },
            "IMMUNOGLOBULIN_HEAVY_CHAIN": {
                "IGHM",
                "IGHD",
                "IGHG1",
                "IGHG2",
                "IGHG3",
                "IGHG4",
                "IGHA1",
                "IGHA2",
                "IGHE",
            },
            "IMMUNOGLOBULIN_LIGHT_CHAIN": {
                "IGLC1",
                "IGLC2",
                "IGLC3",
                "IGLC6",
                "IGLC7",
                "IGKC",
                "IGKV1",
                "IGKV2",
                "IGKV3",
                "IGKV4",
            },
            "GERMINAL_CENTER_REACTION": {
                "BCL6",
                "AICDA",
                "CD95",
                "CXCR4",
                "CXCR5",
                "SLAMF1",
                "MEF2B",
                "EZH2",
                "PCNA",
                "MKI67",
                "TOP2A",
                "BIRC5",
                "CCNB1",
                "CDK1",
            },
        }

    def _save_gene_sets_for_r(self, gene_sets: Dict[str, Set[str]]) -> None:
        """Save gene sets in format suitable for R analysis"""

        # Save as GMT format (standard for GSEA)
        gmt_file = self.cache_dir / "bcell_msigdb_gene_sets.gmt"

        with open(gmt_file, "w") as f:
            for pathway_name, gene_set in gene_sets.items():
                # GMT format: pathway_name\tdescription\tgene1\tgene2\t...
                genes_str = "\t".join(sorted(gene_set))
                f.write(f"{pathway_name}\t{pathway_name}\t{genes_str}\n")

        logger.info(f"Gene sets saved to {gmt_file} for R analysis")

        # Also save as JSON for Python analysis
        json_file = self.cache_dir / "bcell_msigdb_gene_sets.json"
        json_data = {k: list(v) for k, v in gene_sets.items()}

        with open(json_file, "w") as f:
            json.dump(json_data, f, indent=2)

    def analyze_gene_set_overlap(
        self, query_genes: List[str], gene_sets: Optional[Dict[str, Set[str]]] = None
    ) -> pd.DataFrame:
        """
        Analyze overlap between query genes and MSigDB gene sets

        Parameters:
        -----------
        query_genes : List[str]
            List of query gene symbols
        gene_sets : Optional[Dict[str, Set[str]]]
            Gene sets to analyze (if None, uses B-cell specific sets)

        Returns:
        --------
        pd.DataFrame
            Overlap analysis results
        """
        if gene_sets is None:
            gene_sets = self.prepare_bcell_msigdb_gene_sets()

        query_set = set(query_genes)
        overlap_results = []

        for pathway_name, pathway_genes in gene_sets.items():
            overlap_genes = query_set.intersection(pathway_genes)
            overlap_count = len(overlap_genes)

            if overlap_count > 0:
                overlap_results.append(
                    {
                        "pathway": pathway_name,
                        "pathway_size": len(pathway_genes),
                        "query_size": len(query_set),
                        "overlap_count": overlap_count,
                        "overlap_genes": ",".join(sorted(overlap_genes)),
                        "overlap_fraction": overlap_count / len(pathway_genes),
                        "query_fraction": overlap_count / len(query_set),
                    }
                )

        results_df = pd.DataFrame(overlap_results)

        if not results_df.empty:
            # Sort by overlap significance
            results_df = results_df.sort_values("overlap_count", ascending=False)

        return results_df

    def get_pathway_genes(self, pathway_name: str) -> Set[str]:
        """Get genes for a specific pathway"""
        gene_sets = self.prepare_bcell_msigdb_gene_sets()
        return gene_sets.get(pathway_name, set())

    def search_pathways(self, keyword: str) -> List[str]:
        """Search for pathways containing a keyword"""
        gene_sets = self.prepare_bcell_msigdb_gene_sets()
        matching_pathways = []

        keyword_upper = keyword.upper()
        for pathway_name in gene_sets.keys():
            if keyword_upper in pathway_name.upper():
                matching_pathways.append(pathway_name)

        return matching_pathways


def prepare_bcell_msigdb_gene_sets(config: Config) -> Dict[str, Set[str]]:
    """
    Convenience function to prepare B-cell MSigDB gene sets

    Parameters:
    -----------
    config : Config
        PlasmaFlow configuration

    Returns:
    --------
    Dict[str, Set[str]]
        B-cell relevant gene sets
    """
    analyzer = MSigDBAnalyzer(config)
    return analyzer.prepare_bcell_msigdb_gene_sets()
