"""
Statistical analysis for gene distance data

This module provides comprehensive statistical testing and visualization
for gene proximity analysis results.
"""

import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from ..config import Config
from ..utils import get_logger

logger = get_logger(__name__)


class GeneDistanceStatistics:
    """Statistical analysis for gene distance data"""

    def __init__(self, config: Config):
        """
        Initialize statistical analyzer

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config

        # Get visualization parameters
        self.viz_params = config.genomics.get("visualization", {})

        # Color schemes
        self.category_colors = self.viz_params.get(
            "category_colors",
            {
                "up": "#E69F00",  # Orange
                "down": "#56B4E9",  # Sky blue
                "common": "#009E73",  # Bluish green
            },
        )

        self.cluster_colors = self.viz_params.get(
            "cluster_colors",
            {
                "1": "#E69F00",
                "2": "#56B4E9",
                "3": "#009E73",
                "4": "#F0E442",
                "5": "#0072B2",
                "6": "#D55E00",
                "7": "#CC79A7",
            },
        )

        # Statistical parameters
        self.stat_params = self.viz_params.get(
            "statistical",
            {
                "confidence_level": 0.95,
                "n_bootstrap": 1000,
                "alpha": 0.05,
                "multiple_testing_correction": "fdr_bh",
            },
        )

        # Set publication quality plotting parameters
        self._set_plot_params()

    def _set_plot_params(self):
        """Set publication-quality plotting parameters"""
        plt.rcParams.update(
            {
                "font.size": 12,
                "axes.titlesize": 14,
                "axes.labelsize": 12,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10,
                "legend.fontsize": 10,
                "figure.titlesize": 16,
                "figure.dpi": 100,
            }
        )

    def create_comprehensive_analysis(
        self,
        results: List,  # List of ProximityResult objects
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Dict[str, Path]:
        """
        Create comprehensive statistical analysis and visualizations

        Args:
            results: List of ProximityResult objects
            comparison_name: Name of comparison
            method_name: Analysis method name
            output_dir: Output directory

        Returns:
            Dictionary of created visualization files
        """

        logger.info(
            f"Creating comprehensive statistical analysis for {comparison_name}"
        )

        # Organize data by category and cluster
        category_data, cluster_data = self._organize_data(results)

        visualization_files = {}

        # Category-level analysis
        if category_data:
            category_files = self._create_category_analysis(
                category_data, comparison_name, method_name, output_dir
            )
            visualization_files.update(category_files)

        # Cluster-level analysis
        for category, clusters in cluster_data.items():
            if len(clusters) >= 2:  # Need at least 2 clusters for comparison
                cluster_files = self._create_cluster_analysis(
                    clusters, category, comparison_name, method_name, output_dir
                )
                visualization_files.update(
                    {f"{category}_{k}": v for k, v in cluster_files.items()}
                )

        logger.info(f"Created {len(visualization_files)} visualization files")
        return visualization_files

    def _organize_data(self, results: List) -> Tuple[Dict, Dict]:
        """Organize results by category and cluster"""

        category_data = {}
        cluster_data = defaultdict(dict)

        for result in results:
            if len(result.distances_df) > 0:
                if result.file_type == "category":
                    category_data[result.category] = result.distances_df
                    logger.debug(
                        f"Category data: {result.category} - {len(result.distances_df)} genes"
                    )
                else:
                    cluster_data[result.category][result.cluster] = result.distances_df
                    logger.debug(
                        f"Cluster data: {result.category} cluster {result.cluster} - {len(result.distances_df)} genes"
                    )

        return category_data, cluster_data

    def _create_category_analysis(
        self,
        category_data: Dict[str, pd.DataFrame],
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Dict[str, Path]:
        """Create category-level statistical analysis"""

        logger.info("Creating category-level statistical analysis")

        # Calculate statistical tests
        test_results = self._calculate_statistical_tests(category_data, "category")

        # Save statistical results
        stats_file = (
            output_dir
            / f"{comparison_name}_{method_name}_category_statistical_tests.csv"
        )
        if test_results:
            stats_df = pd.DataFrame(test_results)
            stats_df.to_csv(stats_file, index=False)

        # Create comprehensive plot
        plot_file = self._create_category_plot(
            category_data, test_results, comparison_name, method_name, output_dir
        )

        return {"category_plot": plot_file, "category_stats": stats_file}

    def _create_cluster_analysis(
        self,
        cluster_data: Dict[str, pd.DataFrame],
        category: str,
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Dict[str, Path]:
        """Create cluster-level statistical analysis"""

        logger.info(f"Creating {category} cluster statistical analysis")

        # Calculate statistical tests
        test_results = self._calculate_statistical_tests(cluster_data, "cluster")

        # Save statistical results
        stats_file = (
            output_dir
            / f"{comparison_name}_{method_name}_{category}_cluster_statistical_tests.csv"
        )
        if test_results:
            stats_df = pd.DataFrame(test_results)
            stats_df.to_csv(stats_file, index=False)

        # Create comprehensive plot
        plot_file = self._create_cluster_plot(
            cluster_data,
            test_results,
            category,
            comparison_name,
            method_name,
            output_dir,
        )

        # Create overview plot
        overview_file = self._create_cluster_overview(
            cluster_data, category, comparison_name, method_name, output_dir
        )

        return {
            "cluster_plot": plot_file,
            "cluster_overview": overview_file,
            "cluster_stats": stats_file,
        }

    def _calculate_statistical_tests(
        self, data: Dict[str, pd.DataFrame], level: str
    ) -> List[Dict[str, Any]]:
        """Calculate pairwise statistical tests"""

        groups = list(data.keys())
        test_results = []
        p_values = []

        logger.debug(f"Performing statistical tests between {len(groups)} {level}s")

        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                group1, group2 = groups[i], groups[j]
                data1 = data[group1]["distance_to_loop"]
                data2 = data[group2]["distance_to_loop"]

                # Mann-Whitney U test (non-parametric)
                u_stat, p_val = mannwhitneyu(data1, data2, alternative="two-sided")

                test_results.append(
                    {
                        f"{level}1": group1,
                        f"{level}2": group2,
                        "u_statistic": u_stat,
                        "p_value": p_val,
                        "median1": np.median(data1),
                        "median2": np.median(data2),
                        "n1": len(data1),
                        "n2": len(data2),
                        "mean1": np.mean(data1),
                        "mean2": np.mean(data2),
                    }
                )
                p_values.append(p_val)

        # Apply multiple testing correction
        if p_values:
            correction_method = self.stat_params.get(
                "multiple_testing_correction", "fdr_bh"
            )
            corrected_p = multipletests(p_values, method=correction_method)[1]
            for i, result in enumerate(test_results):
                result["p_adjusted"] = corrected_p[i]
                result["significant"] = corrected_p[i] < self.stat_params.get(
                    "alpha", 0.05
                )

        # Log results
        for result in test_results:
            significance = ""
            if result["p_adjusted"] < 0.001:
                significance = " ***"
            elif result["p_adjusted"] < 0.01:
                significance = " **"
            elif result["p_adjusted"] < 0.05:
                significance = " *"

            logger.debug(
                f"{result[f'{level}1']} vs {result[f'{level}2']}: "
                f"p = {result['p_value']:.4f}, "
                f"p_adj = {result['p_adjusted']:.4f}{significance}"
            )

        return test_results

    def _create_category_plot(
        self,
        category_data: Dict[str, pd.DataFrame],
        test_results: List[Dict],
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Path:
        """Create comprehensive category comparison plot"""

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        categories = list(category_data.keys())
        x_pos = np.arange(len(categories))

        # Plot 1: Median distances with confidence intervals
        medians = []
        ci_lowers = []
        ci_uppers = []

        for category in categories:
            data = category_data[category]["distance_to_loop"]
            median_val = np.median(data)
            ci_lower, ci_upper = self._calculate_confidence_intervals(data)

            medians.append(median_val)
            ci_lowers.append(median_val - ci_lower)
            ci_uppers.append(ci_upper - median_val)

        bars = axes[0, 0].bar(
            x_pos,
            medians,
            yerr=[ci_lowers, ci_uppers],
            capsize=5,
            color=[
                self.category_colors.get(cat, f"C{i}")
                for i, cat in enumerate(categories)
            ],
            alpha=0.8,
            edgecolor="black",
            linewidth=1,
        )

        # Add significance bars
        self._add_significance_bars(
            axes[0, 0], category_data, test_results, x_pos, "category"
        )

        axes[0, 0].set_xlabel("Gene Categories")
        axes[0, 0].set_ylabel("Median TSS Distance to Loop Anchors (bp)")
        axes[0, 0].set_title(
            "TSS Distance by Category\\n(with 95% Confidence Intervals)"
        )
        axes[0, 0].set_xticks(x_pos)
        axes[0, 0].set_xticklabels([cat.upper() for cat in categories])
        axes[0, 0].grid(True, alpha=0.3)

        # Plot 2: Cumulative distributions
        for category, df in category_data.items():
            sorted_dist = np.sort(df["distance_to_loop"])
            cumulative = np.arange(1, len(sorted_dist) + 1) / len(sorted_dist)
            axes[0, 1].plot(
                sorted_dist,
                cumulative,
                label=f"{category.upper()} (n={len(df)})",
                color=self.category_colors.get(category, "gray"),
                linewidth=2,
            )

            # Mark median
            p50 = np.percentile(sorted_dist, 50)
            axes[0, 1].axvline(
                p50,
                color=self.category_colors.get(category, "gray"),
                linestyle="--",
                alpha=0.7,
            )

        axes[0, 1].set_xlabel("TSS Distance to Loop Anchors (bp)")
        axes[0, 1].set_ylabel("Cumulative Probability")
        axes[0, 1].set_title("Cumulative Distribution Functions")
        axes[0, 1].set_xscale("log")
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        # Plot 3: Regulatory zone enrichment
        regulatory_zones = [
            "TSS-proximal",
            "Promoter",
            "Local-enhancer",
            "Nearby-enhancer",
            "Distal-enhancer",
        ]

        zone_data = {}
        for category, df in category_data.items():
            if "regulatory_category" in df.columns:
                zone_counts = df["regulatory_category"].value_counts()
                zone_data[category] = zone_counts.reindex(
                    regulatory_zones, fill_value=0
                )

        if zone_data:
            x = np.arange(len(regulatory_zones))
            width = 0.25

            for i, (category, counts) in enumerate(zone_data.items()):
                axes[1, 0].bar(
                    x + i * width,
                    counts.values,
                    width,
                    label=category.upper(),
                    color=self.category_colors.get(category, f"C{i}"),
                    alpha=0.8,
                )

            axes[1, 0].set_xlabel("Regulatory Zones")
            axes[1, 0].set_ylabel("Number of Genes")
            axes[1, 0].set_title("Gene Distribution Across Regulatory Zones")
            axes[1, 0].set_xticks(x + width)
            axes[1, 0].set_xticklabels(regulatory_zones, rotation=45, ha="right")
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].text(
                0.5,
                0.5,
                "No regulatory category data",
                ha="center",
                va="center",
                transform=axes[1, 0].transAxes,
            )

        # Plot 4: Statistics table
        axes[1, 1].axis("off")
        stats_data = []

        for category, df in category_data.items():
            within_loop = (
                df["is_within_loop"].sum() if "is_within_loop" in df.columns else 0
            )
            stats_data.append(
                [
                    category.upper(),
                    len(df),
                    f"{np.median(df['distance_to_loop']):.0f}",
                    f"{np.mean(df['distance_to_loop']):.0f}",
                    f"{within_loop}",
                    f"{within_loop / len(df) * 100:.1f}%" if len(df) > 0 else "0%",
                ]
            )

        if stats_data:
            table = axes[1, 1].table(
                cellText=stats_data,
                colLabels=[
                    "Category",
                    "Total Genes",
                    "Median Dist (bp)",
                    "Mean Dist (bp)",
                    "Within Loop",
                    "% Within Loop",
                ],
                cellLoc="center",
                loc="center",
                bbox=[0, 0, 1, 1],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1, 2)

            # Color table cells
            for i in range(1, len(stats_data) + 1):
                category = list(category_data.keys())[i - 1]
                for j in range(len(stats_data[0])):
                    table[(i, j)].set_facecolor(
                        self.category_colors.get(category, "lightgray")
                    )
                    table[(i, j)].set_alpha(0.3)

        plt.suptitle(
            f"{comparison_name} - {method_name}: Category Statistical Analysis",
            fontsize=16,
        )
        plt.tight_layout()

        # Save plot
        plot_file = (
            output_dir
            / f"{comparison_name}_{method_name}_category_statistical_analysis.png"
        )
        plt.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Category statistical plot saved: {plot_file}")
        return plot_file

    def _create_cluster_plot(
        self,
        cluster_data: Dict[str, pd.DataFrame],
        test_results: List[Dict],
        category: str,
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Path:
        """Create comprehensive cluster comparison plot"""

        if len(cluster_data) < 2:
            logger.warning(
                f"Cannot create cluster plot for {category}: need at least 2 clusters"
            )
            return None

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        clusters = sorted(cluster_data.keys())
        x_pos = np.arange(len(clusters))

        # Plot 1: Median distances with confidence intervals
        medians = []
        ci_lowers = []
        ci_uppers = []

        for cluster in clusters:
            data = cluster_data[cluster]["distance_to_loop"]
            median_val = np.median(data)
            ci_lower, ci_upper = self._calculate_confidence_intervals(data)

            medians.append(median_val)
            ci_lowers.append(median_val - ci_lower)
            ci_uppers.append(ci_upper - median_val)

        bars = axes[0, 0].bar(
            x_pos,
            medians,
            yerr=[ci_lowers, ci_uppers],
            capsize=5,
            color=[
                self.cluster_colors.get(str(c), f"C{i}") for i, c in enumerate(clusters)
            ],
            alpha=0.8,
            edgecolor="black",
            linewidth=1,
        )

        # Add significance bars
        self._add_significance_bars(
            axes[0, 0], cluster_data, test_results, x_pos, "cluster"
        )

        axes[0, 0].set_xlabel("Clusters")
        axes[0, 0].set_ylabel("Median TSS Distance to Loop Anchors (bp)")
        axes[0, 0].set_title(
            f"{category.upper()} - TSS Distance by Cluster\\n(with 95% Confidence Intervals)"
        )
        axes[0, 0].set_xticks(x_pos)
        axes[0, 0].set_xticklabels([f"Cluster {c}" for c in clusters])
        axes[0, 0].grid(True, alpha=0.3)

        # Plot 2: Cumulative distributions
        for cluster, df in cluster_data.items():
            sorted_dist = np.sort(df["distance_to_loop"])
            cumulative = np.arange(1, len(sorted_dist) + 1) / len(sorted_dist)
            axes[0, 1].plot(
                sorted_dist,
                cumulative,
                label=f"Cluster {cluster} (n={len(df)})",
                color=self.cluster_colors.get(str(cluster), f"C{cluster}"),
                linewidth=2,
            )

            # Mark median
            p50 = np.percentile(sorted_dist, 50)
            axes[0, 1].axvline(
                p50,
                color=self.cluster_colors.get(str(cluster), f"C{cluster}"),
                linestyle="--",
                alpha=0.7,
            )

        axes[0, 1].set_xlabel("TSS Distance to Loop Anchors (bp)")
        axes[0, 1].set_ylabel("Cumulative Probability")
        axes[0, 1].set_title("Cumulative Distribution Functions")
        axes[0, 1].set_xscale("log")
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)

        # Plot 3: Regulatory zone enrichment
        regulatory_zones = [
            "TSS-proximal",
            "Promoter",
            "Local-enhancer",
            "Nearby-enhancer",
            "Distal-enhancer",
        ]

        zone_data = {}
        for cluster, df in cluster_data.items():
            if "regulatory_category" in df.columns:
                zone_counts = df["regulatory_category"].value_counts()
                zone_data[cluster] = zone_counts.reindex(regulatory_zones, fill_value=0)

        if zone_data:
            x = np.arange(len(regulatory_zones))
            width = 0.8 / len(clusters)

            for i, cluster in enumerate(clusters):
                counts = zone_data.get(cluster, pd.Series([0] * len(regulatory_zones)))
                axes[1, 0].bar(
                    x + i * width,
                    counts.values,
                    width,
                    label=f"Cluster {cluster}",
                    color=self.cluster_colors.get(str(cluster), f"C{i}"),
                    alpha=0.8,
                )

            axes[1, 0].set_xlabel("Regulatory Zones")
            axes[1, 0].set_ylabel("Number of Genes")
            axes[1, 0].set_title("Gene Distribution Across Regulatory Zones")
            axes[1, 0].set_xticks(x + width * (len(clusters) - 1) / 2)
            axes[1, 0].set_xticklabels(regulatory_zones, rotation=45, ha="right")
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].text(
                0.5,
                0.5,
                "No regulatory category data",
                ha="center",
                va="center",
                transform=axes[1, 0].transAxes,
            )

        # Plot 4: Statistics table
        axes[1, 1].axis("off")
        stats_data = []

        for cluster in clusters:
            df = cluster_data[cluster]
            within_loop = (
                df["is_within_loop"].sum() if "is_within_loop" in df.columns else 0
            )
            stats_data.append(
                [
                    f"Cluster {cluster}",
                    len(df),
                    f"{np.median(df['distance_to_loop']):.0f}",
                    f"{np.mean(df['distance_to_loop']):.0f}",
                    f"{within_loop}",
                    f"{within_loop / len(df) * 100:.1f}%" if len(df) > 0 else "0%",
                ]
            )

        if stats_data:
            table = axes[1, 1].table(
                cellText=stats_data,
                colLabels=[
                    "Cluster",
                    "Total Genes",
                    "Median Dist (bp)",
                    "Mean Dist (bp)",
                    "Within Loop",
                    "% Within Loop",
                ],
                cellLoc="center",
                loc="center",
                bbox=[0, 0, 1, 1],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1, 1.8)

            # Color table cells
            for i in range(1, len(stats_data) + 1):
                cluster = clusters[i - 1]
                for j in range(len(stats_data[0])):
                    table[(i, j)].set_facecolor(
                        self.cluster_colors.get(str(cluster), "lightgray")
                    )
                    table[(i, j)].set_alpha(0.3)

        plt.suptitle(
            f"{category.upper()} - Detailed Cluster Statistical Analysis", fontsize=16
        )
        plt.tight_layout()

        # Save plot
        plot_file = (
            output_dir
            / f"{comparison_name}_{method_name}_{category}_cluster_statistical_analysis.png"
        )
        plt.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"{category} cluster statistical plot saved: {plot_file}")
        return plot_file

    def _create_cluster_overview(
        self,
        cluster_data: Dict[str, pd.DataFrame],
        category: str,
        comparison_name: str,
        method_name: str,
        output_dir: Path,
    ) -> Path:
        """Create cluster overview plot"""

        plt.figure(figsize=(10, 6))

        for cluster, df in cluster_data.items():
            plt.hist(
                np.log10(df["distance_to_loop"] + 1),
                bins=30,
                alpha=0.6,
                label=f"Cluster {cluster} (n={len(df)})",
                color=self.cluster_colors.get(str(cluster), f"C{cluster}"),
                density=True,
            )

        plt.xlabel("Log10(TSS Distance + 1)")
        plt.ylabel("Density")
        plt.title(f"{category.upper()} - Gene Distance Distribution by Cluster")
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Save plot
        overview_file = (
            output_dir
            / f"{comparison_name}_{method_name}_{category}_cluster_overview.png"
        )
        plt.savefig(overview_file, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"{category} cluster overview plot saved: {overview_file}")
        return overview_file

    def _calculate_confidence_intervals(
        self, data: np.ndarray, confidence: float = None
    ) -> Tuple[float, float]:
        """Calculate bootstrap confidence intervals"""

        if confidence is None:
            confidence = self.stat_params.get("confidence_level", 0.95)

        n_bootstrap = self.stat_params.get("n_bootstrap", 1000)
        bootstrap_medians = []

        for _ in range(n_bootstrap):
            sample = np.random.choice(data, size=len(data), replace=True)
            bootstrap_medians.append(np.median(sample))

        alpha = 1 - confidence
        lower = np.percentile(bootstrap_medians, (alpha / 2) * 100)
        upper = np.percentile(bootstrap_medians, (1 - alpha / 2) * 100)

        return lower, upper

    def _add_significance_bars(
        self,
        ax,
        data: Dict[str, pd.DataFrame],
        test_results: List[Dict],
        x_positions: np.ndarray,
        level: str,
    ):
        """Add significance bars to plots"""

        # Calculate max height for positioning bars
        max_height = 0
        for group, df in data.items():
            values = df["distance_to_loop"]
            median_val = np.median(values)
            sem_val = np.std(values) / np.sqrt(len(values))
            max_height = max(max_height, median_val + sem_val)

        y_offset = max_height * 0.1
        current_y = max_height + y_offset

        groups = list(data.keys())

        for result in test_results:
            if result.get("significant", False):
                group1, group2 = result[f"{level}1"], result[f"{level}2"]

                try:
                    i = groups.index(group1)
                    j = groups.index(group2)
                except ValueError:
                    continue

                # Draw significance bar
                h = max_height * 0.03

                ax.plot(
                    [x_positions[i], x_positions[i], x_positions[j], x_positions[j]],
                    [current_y, current_y + h, current_y + h, current_y],
                    lw=1.5,
                    c="black",
                )

                # Add significance annotation
                if result["p_adjusted"] < 0.001:
                    sig_text = "***"
                elif result["p_adjusted"] < 0.01:
                    sig_text = "**"
                else:
                    sig_text = "*"

                ax.text(
                    (x_positions[i] + x_positions[j]) * 0.5,
                    current_y + h,
                    sig_text,
                    ha="center",
                    va="bottom",
                    fontweight="bold",
                    fontsize=10,
                )

                current_y += h + max_height * 0.05


def calculate_gene_distance_statistics(
    distances_df: pd.DataFrame, grouping_column: str = "category"
) -> Dict[str, Any]:
    """
    Calculate comprehensive statistics for gene distances

    Args:
        distances_df: DataFrame with gene distance data
        grouping_column: Column to group by for statistics

    Returns:
        Dictionary with statistical summaries
    """

    if len(distances_df) == 0:
        return {"error": "Empty distance DataFrame"}

    statistics = {"total_genes": len(distances_df), "groups": {}}

    # Overall statistics
    if "distance_to_loop" in distances_df.columns:
        statistics["overall"] = {
            "median_distance": distances_df["distance_to_loop"].median(),
            "mean_distance": distances_df["distance_to_loop"].mean(),
            "std_distance": distances_df["distance_to_loop"].std(),
            "min_distance": distances_df["distance_to_loop"].min(),
            "max_distance": distances_df["distance_to_loop"].max(),
            "quartiles": {
                "q25": distances_df["distance_to_loop"].quantile(0.25),
                "q75": distances_df["distance_to_loop"].quantile(0.75),
            },
        }

    # Group-wise statistics
    if grouping_column in distances_df.columns:
        for group in distances_df[grouping_column].unique():
            group_data = distances_df[distances_df[grouping_column] == group]

            if "distance_to_loop" in group_data.columns:
                within_loop = group_data.get(
                    "is_within_loop", pd.Series([False] * len(group_data))
                ).sum()

                statistics["groups"][group] = {
                    "n_genes": len(group_data),
                    "median_distance": group_data["distance_to_loop"].median(),
                    "mean_distance": group_data["distance_to_loop"].mean(),
                    "std_distance": group_data["distance_to_loop"].std(),
                    "genes_within_loop": within_loop,
                    "percent_within_loop": (
                        (within_loop / len(group_data)) * 100
                        if len(group_data) > 0
                        else 0
                    ),
                }

                # Regulatory category breakdown
                if "regulatory_category" in group_data.columns:
                    reg_counts = group_data["regulatory_category"].value_counts()
                    statistics["groups"][group][
                        "regulatory_categories"
                    ] = reg_counts.to_dict()

    return statistics
