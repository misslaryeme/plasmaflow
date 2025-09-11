"""
ClusterProfiler integration for pathway enrichment analysis

This module provides Python wrappers for R's clusterProfiler package,
enabling comprehensive pathway analysis including GO, KEGG, and MSigDB enrichment.
"""

import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..core.config import Config
from ..core.exceptions import PlasmaFlowError, RIntegrationError
from .r_interface import RPathwayInterface

logger = logging.getLogger(__name__)


@dataclass
class PathwayEnrichmentResult:
    """Container for pathway enrichment analysis results"""

    method: str
    database: str
    gene_count: int
    significant_pathways: int
    results_df: pd.DataFrame
    plot_paths: Dict[str, Path]
    parameters: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary"""
        result_dict = asdict(self)
        result_dict["results_df"] = self.results_df.to_dict("records")
        result_dict["plot_paths"] = {k: str(v) for k, v in self.plot_paths.items()}
        return result_dict

    def save_results(self, output_dir: Path) -> None:
        """Save results to files"""
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save results table
        results_file = output_dir / f"{self.method}_{self.database}_results.csv"
        self.results_df.to_csv(results_file, index=False)

        # Save summary
        summary_file = output_dir / f"{self.method}_{self.database}_summary.txt"
        with open(summary_file, "w") as f:
            f.write(f"Pathway Enrichment Analysis Summary\n")
            f.write(f"Method: {self.method}\n")
            f.write(f"Database: {self.database}\n")
            f.write(f"Input genes: {self.gene_count}\n")
            f.write(
                f"Significant pathways (p.adjust < 0.05): {self.significant_pathways}\n"
            )


class ClusterProfilerAnalyzer:
    """Python interface for R clusterProfiler pathway analysis"""

    def __init__(self, config: Config):
        self.config = config
        self.r_interface = RPathwayInterface(config)
        self.output_dir = Path(config.analysis_config.output_dir)

    def prepare_gene_lists(
        self, gene_data: pd.DataFrame
    ) -> Tuple[List[str], List[str]]:
        """
        Prepare gene lists for enrichment analysis

        Parameters:
        -----------
        gene_data : pd.DataFrame
            DataFrame with gene expression data including log2FoldChange and padj columns

        Returns:
        --------
        Tuple[List[str], List[str]]
            (significant_genes, background_genes)
        """
        if "log2FoldChange" not in gene_data.columns or "padj" not in gene_data.columns:
            raise PlasmaFlowError(
                "Gene data must contain 'log2FoldChange' and 'padj' columns"
            )

        # Filter for significant genes
        significant_mask = (gene_data["padj"] < 0.05) & (
            np.abs(gene_data["log2FoldChange"]) > 1
        )

        significant_genes = gene_data[significant_mask].index.tolist()
        background_genes = gene_data.dropna(subset=["padj"]).index.tolist()

        logger.info(
            f"Prepared {len(significant_genes)} significant genes from {len(background_genes)} background"
        )

        return significant_genes, background_genes

    def run_ora_analysis(
        self,
        gene_list: List[str],
        background_genes: Optional[List[str]] = None,
        databases: List[str] = None,
        organism: str = "org.Hs.eg.db",
        output_prefix: str = "ora",
    ) -> Dict[str, PathwayEnrichmentResult]:
        """
        Run Over-Representation Analysis (ORA)

        Parameters:
        -----------
        gene_list : List[str]
            List of significant gene symbols
        background_genes : Optional[List[str]]
            Background gene universe
        databases : List[str]
            Databases to query ('GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'MSigDB_C2')
        organism : str
            Organism database
        output_prefix : str
            Prefix for output files

        Returns:
        --------
        Dict[str, PathwayEnrichmentResult]
            Results for each database
        """
        if databases is None:
            databases = ["GO_BP", "GO_MF", "GO_CC", "KEGG"]

        if not gene_list:
            raise PlasmaFlowError("Gene list cannot be empty")

        logger.info(
            f"Running ORA analysis for {len(gene_list)} genes across {len(databases)} databases"
        )

        results = {}
        output_dir = self.output_dir / "enrichment" / output_prefix
        output_dir.mkdir(parents=True, exist_ok=True)

        for database in databases:
            try:
                logger.info(f"Analyzing {database} database")

                # Run R analysis
                result_data = self.r_interface.run_ora_analysis(
                    gene_list=gene_list,
                    background_genes=background_genes,
                    database=database,
                    organism=organism,
                    output_dir=output_dir,
                    output_prefix=f"{output_prefix}_{database}",
                )

                # Load results
                results_file = output_dir / f"{output_prefix}_{database}_results.csv"
                if results_file.exists():
                    results_df = pd.read_csv(results_file)
                    significant_count = len(results_df[results_df["p.adjust"] < 0.05])

                    # Collect plot paths
                    plot_paths = {}
                    for plot_type in ["dotplot", "barplot", "cnetplot"]:
                        plot_file = (
                            output_dir / f"{output_prefix}_{database}_{plot_type}.pdf"
                        )
                        if plot_file.exists():
                            plot_paths[plot_type] = plot_file

                    results[database] = PathwayEnrichmentResult(
                        method="ORA",
                        database=database,
                        gene_count=len(gene_list),
                        significant_pathways=significant_count,
                        results_df=results_df,
                        plot_paths=plot_paths,
                        parameters={
                            "organism": organism,
                            "pvalue_cutoff": 0.05,
                            "qvalue_cutoff": 0.2,
                        },
                    )

                    logger.info(
                        f"{database}: Found {significant_count} significant pathways"
                    )
                else:
                    logger.warning(f"No results file found for {database}")

            except Exception as e:
                logger.error(f"Error analyzing {database}: {str(e)}")
                continue

        return results

    def run_gsea_analysis(
        self,
        ranked_genes: Dict[str, float],
        databases: List[str] = None,
        organism: str = "org.Hs.eg.db",
        output_prefix: str = "gsea",
    ) -> Dict[str, PathwayEnrichmentResult]:
        """
        Run Gene Set Enrichment Analysis (GSEA)

        Parameters:
        -----------
        ranked_genes : Dict[str, float]
            Dictionary of gene symbols to fold change values
        databases : List[str]
            Databases to query
        organism : str
            Organism database
        output_prefix : str
            Prefix for output files

        Returns:
        --------
        Dict[str, PathwayEnrichmentResult]
            Results for each database
        """
        if databases is None:
            databases = ["GO_BP", "GO_MF", "GO_CC", "KEGG"]

        if not ranked_genes:
            raise PlasmaFlowError("Ranked gene list cannot be empty")

        logger.info(
            f"Running GSEA analysis for {len(ranked_genes)} genes across {len(databases)} databases"
        )

        results = {}
        output_dir = self.output_dir / "enrichment" / output_prefix
        output_dir.mkdir(parents=True, exist_ok=True)

        for database in databases:
            try:
                logger.info(f"Analyzing {database} database")

                # Run R analysis
                result_data = self.r_interface.run_gsea_analysis(
                    ranked_genes=ranked_genes,
                    database=database,
                    organism=organism,
                    output_dir=output_dir,
                    output_prefix=f"{output_prefix}_{database}",
                )

                # Load results
                results_file = output_dir / f"{output_prefix}_{database}_results.csv"
                if results_file.exists():
                    results_df = pd.read_csv(results_file)
                    significant_count = len(results_df[results_df["p.adjust"] < 0.05])

                    # Collect plot paths
                    plot_paths = {}
                    for plot_type in ["dotplot", "ridgeplot", "gseaplot"]:
                        plot_file = (
                            output_dir / f"{output_prefix}_{database}_{plot_type}.pdf"
                        )
                        if plot_file.exists():
                            plot_paths[plot_type] = plot_file

                    results[database] = PathwayEnrichmentResult(
                        method="GSEA",
                        database=database,
                        gene_count=len(ranked_genes),
                        significant_pathways=significant_count,
                        results_df=results_df,
                        plot_paths=plot_paths,
                        parameters={
                            "organism": organism,
                            "pvalue_cutoff": 0.05,
                            "eps": 0,
                        },
                    )

                    logger.info(
                        f"{database}: Found {significant_count} significant pathways"
                    )
                else:
                    logger.warning(f"No results file found for {database}")

            except Exception as e:
                logger.error(f"Error analyzing {database}: {str(e)}")
                continue

        return results

    def run_comprehensive_analysis(
        self,
        gene_data: pd.DataFrame,
        output_prefix: str = "comprehensive",
        include_msigdb: bool = True,
    ) -> Dict[str, Dict[str, PathwayEnrichmentResult]]:
        """
        Run comprehensive pathway analysis including both ORA and GSEA

        Parameters:
        -----------
        gene_data : pd.DataFrame
            DataFrame with gene expression data
        output_prefix : str
            Prefix for output files
        include_msigdb : bool
            Whether to include MSigDB analysis

        Returns:
        --------
        Dict[str, Dict[str, PathwayEnrichmentResult]]
            Nested dictionary with method -> database -> results
        """
        logger.info("Running comprehensive pathway enrichment analysis")

        # Prepare gene lists
        significant_genes, background_genes = self.prepare_gene_lists(gene_data)

        # Prepare ranked gene list for GSEA
        ranked_genes = gene_data["log2FoldChange"].dropna().to_dict()

        # Define databases
        databases = ["GO_BP", "GO_MF", "GO_CC", "KEGG"]
        if include_msigdb:
            databases.append("MSigDB_C2")

        results = {}

        # Run ORA analysis
        if significant_genes:
            logger.info("Running ORA analysis")
            ora_results = self.run_ora_analysis(
                gene_list=significant_genes,
                background_genes=background_genes,
                databases=databases,
                output_prefix=f"{output_prefix}_ora",
            )
            results["ORA"] = ora_results
        else:
            logger.warning("No significant genes found for ORA analysis")

        # Run GSEA analysis
        if ranked_genes:
            logger.info("Running GSEA analysis")
            gsea_results = self.run_gsea_analysis(
                ranked_genes=ranked_genes,
                databases=databases,
                output_prefix=f"{output_prefix}_gsea",
            )
            results["GSEA"] = gsea_results
        else:
            logger.warning("No ranked genes found for GSEA analysis")

        # Save summary
        self._save_comprehensive_summary(results, output_prefix)

        return results

    def _save_comprehensive_summary(
        self, results: Dict[str, Dict[str, PathwayEnrichmentResult]], output_prefix: str
    ) -> None:
        """Save comprehensive analysis summary"""
        output_dir = self.output_dir / "enrichment"
        summary_file = output_dir / f"{output_prefix}_summary.txt"

        with open(summary_file, "w") as f:
            f.write("Comprehensive Pathway Enrichment Analysis Summary\n")
            f.write("=" * 50 + "\n\n")

            for method, method_results in results.items():
                f.write(f"{method} Analysis:\n")
                f.write("-" * 20 + "\n")

                total_significant = 0
                for database, result in method_results.items():
                    f.write(
                        f"  {database}: {result.significant_pathways} significant pathways\n"
                    )
                    total_significant += result.significant_pathways

                f.write(f"  Total: {total_significant} significant pathways\n\n")

        logger.info(f"Summary saved to {summary_file}")
