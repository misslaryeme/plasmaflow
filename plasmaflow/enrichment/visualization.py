"""
Visualization functions for pathway enrichment analysis

This module provides comprehensive plotting functions for pathway analysis results,
including both Python-based matplotlib plots and R-based ggplot2 visualizations.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from ..core.config import Config
from ..core.exceptions import PlasmaFlowError
from .clusterprofiler import PathwayEnrichmentResult
from .r_interface import RPathwayInterface

logger = logging.getLogger(__name__)


class PathwayPlotter:
    """Comprehensive plotting for pathway enrichment results"""

    def __init__(self, config: Config):
        self.config = config
        self.r_interface = RPathwayInterface(config)
        self.output_dir = (
            Path(config.analysis_config.output_dir) / "plots" / "enrichment"
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Set plotting style
        plt.style.use("seaborn-v0_8-whitegrid")
        sns.set_palette("husl")

        # Color schemes for different databases
        self.database_colors = {
            "GO_BP": "#1f77b4",
            "GO_MF": "#ff7f0e",
            "GO_CC": "#2ca02c",
            "KEGG": "#d62728",
            "MSigDB_C2": "#9467bd",
            "Reactome": "#8c564b",
        }

    def plot_enrichment_dotplot(
        self,
        results: Dict[str, PathwayEnrichmentResult],
        top_n: int = 20,
        output_file: Optional[str] = None,
    ) -> Path:
        """
        Create dotplot visualization of enrichment results

        Parameters:
        -----------
        results : Dict[str, PathwayEnrichmentResult]
            Enrichment results for multiple databases
        top_n : int
            Number of top pathways to show per database
        output_file : Optional[str]
            Output filename

        Returns:
        --------
        Path
            Path to saved plot
        """
        if output_file is None:
            output_file = "enrichment_dotplot.pdf"

        plot_path = self.output_dir / output_file

        # Prepare data for plotting
        plot_data = []

        for database, result in results.items():
            if result.results_df.empty:
                continue

            # Get top pathways
            top_pathways = result.results_df.nsmallest(top_n, "p.adjust")

            for _, row in top_pathways.iterrows():
                plot_data.append(
                    {
                        "database": database,
                        "pathway": row.get("Description", row.get("ID", "Unknown")),
                        "pvalue": -np.log10(row["p.adjust"]),
                        "count": row.get("Count", row.get("setSize", 0)),
                        "gene_ratio": row.get("GeneRatio", "0/0"),
                    }
                )

        if not plot_data:
            logger.warning("No data available for dotplot")
            return plot_path

        plot_df = pd.DataFrame(plot_data)

        # Create plot
        fig, ax = plt.subplots(figsize=(12, max(8, len(plot_df) * 0.3)))

        # Create dotplot
        for i, (database, group_data) in enumerate(plot_df.groupby("database")):
            y_positions = np.arange(len(group_data)) + i * (top_n + 1)

            scatter = ax.scatter(
                group_data["pvalue"],
                y_positions,
                s=group_data["count"] * 10,
                c=self.database_colors.get(database, "#666666"),
                alpha=0.7,
                label=database,
                edgecolors="black",
                linewidth=0.5,
            )

        # Customize plot
        ax.set_xlabel("-log10(adjusted p-value)", fontsize=12)
        ax.set_ylabel("Pathways", fontsize=12)
        ax.set_title("Pathway Enrichment Results", fontsize=14, fontweight="bold")

        # Set y-tick labels
        y_labels = []
        y_positions = []
        current_pos = 0

        for database, group_data in plot_df.groupby("database"):
            group_positions = np.arange(len(group_data)) + current_pos
            y_positions.extend(group_positions)
            y_labels.extend(
                [
                    (
                        f"{database}: {pathway[:50]}..."
                        if len(pathway) > 50
                        else f"{database}: {pathway}"
                    )
                    for pathway in group_data["pathway"]
                ]
            )
            current_pos += len(group_data) + 1

        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=8)

        # Add legend
        ax.legend(title="Database", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Dotplot saved to {plot_path}")
        return plot_path

    def plot_enrichment_barplot(
        self,
        results: Dict[str, PathwayEnrichmentResult],
        top_n: int = 10,
        output_file: Optional[str] = None,
    ) -> Path:
        """
        Create barplot visualization of enrichment results
        """
        if output_file is None:
            output_file = "enrichment_barplot.pdf"

        plot_path = self.output_dir / output_file

        # Prepare data
        plot_data = []

        for database, result in results.items():
            if result.results_df.empty:
                continue

            top_pathways = result.results_df.nsmallest(top_n, "p.adjust")

            for _, row in top_pathways.iterrows():
                plot_data.append(
                    {
                        "database": database,
                        "pathway": row.get("Description", row.get("ID", "Unknown"))[
                            :50
                        ],
                        "pvalue": -np.log10(row["p.adjust"]),
                        "count": row.get("Count", row.get("setSize", 0)),
                    }
                )

        if not plot_data:
            logger.warning("No data available for barplot")
            return plot_path

        plot_df = pd.DataFrame(plot_data)

        # Create subplots for each database
        databases = plot_df["database"].unique()
        n_databases = len(databases)

        fig, axes = plt.subplots(
            n_databases, 1, figsize=(12, 4 * n_databases), constrained_layout=True
        )

        if n_databases == 1:
            axes = [axes]

        for i, database in enumerate(databases):
            ax = axes[i]
            db_data = plot_df[plot_df["database"] == database].sort_values("pvalue")

            bars = ax.barh(
                db_data["pathway"],
                db_data["pvalue"],
                color=self.database_colors.get(database, "#666666"),
                alpha=0.7,
                edgecolor="black",
                linewidth=0.5,
            )

            ax.set_xlabel("-log10(adjusted p-value)")
            ax.set_title(f"{database} Enrichment", fontweight="bold")

            # Add value labels on bars
            for j, (bar, count) in enumerate(zip(bars, db_data["count"])):
                ax.text(
                    bar.get_width() + 0.1,
                    bar.get_y() + bar.get_height() / 2,
                    f"n={count}",
                    va="center",
                    fontsize=8,
                )

        plt.suptitle("Pathway Enrichment Results", fontsize=16, fontweight="bold")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Barplot saved to {plot_path}")
        return plot_path

    def plot_enrichment_heatmap(
        self,
        results: Dict[str, Dict[str, PathwayEnrichmentResult]],
        top_n: int = 15,
        output_file: Optional[str] = None,
    ) -> Path:
        """
        Create heatmap comparing enrichment across methods and databases

        Parameters:
        -----------
        results : Dict[str, Dict[str, PathwayEnrichmentResult]]
            Nested results: method -> database -> results
        """
        if output_file is None:
            output_file = "enrichment_heatmap.pdf"

        plot_path = self.output_dir / output_file

        # Collect all pathways and their p-values
        pathway_data = {}

        for method, method_results in results.items():
            for database, result in method_results.items():
                if result.results_df.empty:
                    continue

                top_pathways = result.results_df.nsmallest(top_n, "p.adjust")

                for _, row in top_pathways.iterrows():
                    pathway_name = row.get("Description", row.get("ID", "Unknown"))
                    key = f"{method}_{database}"

                    if pathway_name not in pathway_data:
                        pathway_data[pathway_name] = {}

                    pathway_data[pathway_name][key] = -np.log10(row["p.adjust"])

        if not pathway_data:
            logger.warning("No data available for heatmap")
            return plot_path

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(pathway_data).T.fillna(0)

        # Create heatmap
        plt.figure(figsize=(12, max(8, len(heatmap_df) * 0.3)))

        sns.heatmap(
            heatmap_df,
            cmap="YlOrRd",
            annot=True,
            fmt=".1f",
            cbar_kws={"label": "-log10(adjusted p-value)"},
            linewidths=0.5,
        )

        plt.title("Pathway Enrichment Heatmap", fontsize=14, fontweight="bold")
        plt.xlabel("Method_Database")
        plt.ylabel("Pathways")
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)

        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Heatmap saved to {plot_path}")
        return plot_path

    def plot_gene_concept_network(
        self,
        result: PathwayEnrichmentResult,
        top_n: int = 10,
        output_file: Optional[str] = None,
    ) -> Path:
        """
        Create gene-concept network plot using R
        """
        if output_file is None:
            output_file = f"gene_concept_network_{result.database}.pdf"

        plot_path = self.output_dir / output_file

        try:
            # Use R interface for advanced network plotting
            self.r_interface.create_gene_concept_network(
                results_file=self.output_dir.parent
                / "enrichment"
                / f"{result.method.lower()}_{result.database}_results.csv",
                output_path=plot_path,
                top_n=top_n,
            )

            logger.info(f"Gene-concept network saved to {plot_path}")

        except Exception as e:
            logger.error(f"Failed to create gene-concept network: {str(e)}")

        return plot_path

    def create_comprehensive_pathway_plots(
        self,
        results: Dict[str, Dict[str, PathwayEnrichmentResult]],
        output_prefix: str = "comprehensive",
    ) -> Dict[str, Path]:
        """
        Create comprehensive set of pathway visualization plots

        Parameters:
        -----------
        results : Dict[str, Dict[str, PathwayEnrichmentResult]]
            Nested results: method -> database -> results
        output_prefix : str
            Prefix for output files

        Returns:
        --------
        Dict[str, Path]
            Dictionary mapping plot type to file path
        """
        logger.info("Creating comprehensive pathway plots")

        plot_paths = {}

        # Create plots for each method separately
        for method, method_results in results.items():
            if not method_results:
                continue

            method_prefix = f"{output_prefix}_{method.lower()}"

            # Dotplot
            try:
                dotplot_path = self.plot_enrichment_dotplot(
                    method_results, output_file=f"{method_prefix}_dotplot.pdf"
                )
                plot_paths[f"{method}_dotplot"] = dotplot_path
            except Exception as e:
                logger.error(f"Failed to create {method} dotplot: {str(e)}")

            # Barplot
            try:
                barplot_path = self.plot_enrichment_barplot(
                    method_results, output_file=f"{method_prefix}_barplot.pdf"
                )
                plot_paths[f"{method}_barplot"] = barplot_path
            except Exception as e:
                logger.error(f"Failed to create {method} barplot: {str(e)}")

        # Combined heatmap
        try:
            heatmap_path = self.plot_enrichment_heatmap(
                results, output_file=f"{output_prefix}_heatmap.pdf"
            )
            plot_paths["combined_heatmap"] = heatmap_path
        except Exception as e:
            logger.error(f"Failed to create combined heatmap: {str(e)}")

        # Summary plot
        try:
            summary_path = self.plot_enrichment_summary(
                results, output_file=f"{output_prefix}_summary.pdf"
            )
            plot_paths["summary"] = summary_path
        except Exception as e:
            logger.error(f"Failed to create summary plot: {str(e)}")

        logger.info(f"Created {len(plot_paths)} pathway plots")
        return plot_paths

    def plot_enrichment_summary(
        self,
        results: Dict[str, Dict[str, PathwayEnrichmentResult]],
        output_file: str = "enrichment_summary.pdf",
    ) -> Path:
        """Create summary overview of enrichment results"""

        plot_path = self.output_dir / output_file

        # Collect summary statistics
        summary_data = []

        for method, method_results in results.items():
            for database, result in method_results.items():
                summary_data.append(
                    {
                        "Method": method,
                        "Database": database,
                        "Total_Pathways": (
                            len(result.results_df) if not result.results_df.empty else 0
                        ),
                        "Significant_Pathways": result.significant_pathways,
                        "Input_Genes": result.gene_count,
                    }
                )

        if not summary_data:
            logger.warning("No summary data available")
            return plot_path

        summary_df = pd.DataFrame(summary_data)

        # Create summary plots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle(
            "Pathway Enrichment Analysis Summary", fontsize=16, fontweight="bold"
        )

        # Plot 1: Significant pathways by method and database
        pivot_sig = summary_df.pivot(
            index="Database", columns="Method", values="Significant_Pathways"
        )
        pivot_sig.plot(kind="bar", ax=ax1, color=["skyblue", "lightcoral"])
        ax1.set_title("Significant Pathways by Method")
        ax1.set_ylabel("Number of Significant Pathways")
        ax1.tick_params(axis="x", rotation=45)

        # Plot 2: Total vs significant pathways
        for method in summary_df["Method"].unique():
            method_data = summary_df[summary_df["Method"] == method]
            ax2.scatter(
                method_data["Total_Pathways"],
                method_data["Significant_Pathways"],
                label=method,
                s=100,
                alpha=0.7,
            )

        ax2.set_xlabel("Total Pathways")
        ax2.set_ylabel("Significant Pathways")
        ax2.set_title("Total vs Significant Pathways")
        ax2.legend()

        # Plot 3: Pathway count by database
        db_counts = (
            summary_df.groupby("Database")["Significant_Pathways"].sum().sort_values()
        )
        db_counts.plot(kind="barh", ax=ax3, color="lightgreen")
        ax3.set_title("Significant Pathways by Database")
        ax3.set_xlabel("Number of Significant Pathways")

        # Plot 4: Input genes distribution
        summary_df.boxplot(column="Input_Genes", by="Method", ax=ax4)
        ax4.set_title("Input Gene Distribution by Method")
        ax4.set_ylabel("Number of Input Genes")
        plt.suptitle("")  # Remove automatic title

        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Summary plot saved to {plot_path}")
        return plot_path


def create_comprehensive_pathway_plots(
    results: Dict[str, Dict[str, PathwayEnrichmentResult]],
    config: Config,
    output_prefix: str = "comprehensive",
) -> Dict[str, Path]:
    """
    Convenience function to create comprehensive pathway plots

    Parameters:
    -----------
    results : Dict[str, Dict[str, PathwayEnrichmentResult]]
        Nested enrichment results
    config : Config
        PlasmaFlow configuration
    output_prefix : str
        Prefix for output files

    Returns:
    --------
    Dict[str, Path]
        Dictionary mapping plot type to file path
    """
    plotter = PathwayPlotter(config)
    return plotter.create_comprehensive_pathway_plots(results, output_prefix)
