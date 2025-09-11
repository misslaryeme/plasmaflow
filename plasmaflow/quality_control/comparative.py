"""
Comparative APA analysis functionality
"""

import logging
import os
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from coolpuppy import coolpup

from ..config import Config
from ..utils import get_logger
from .apa import APAAnalyzer, APAResult, calculate_apa_scores
from .caching import load_cached_apa, save_apa_cache

logger = get_logger(__name__)


class ComparativeAPAAnalyzer:
    """Analyzer for comparative APA analysis between cell types"""

    def __init__(self, config: Config):
        """
        Initialize comparative APA analyzer

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config
        self.apa_analyzer = APAAnalyzer(config)

        # Create output directory
        self.output_dir = (
            Path(config.output_dir) / "quality_control" / "comparative_apa"
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Cache directory
        self.cache_dir = self.output_dir / "cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def perform_pairwise_apa_comparison(
        self,
        cell_type_data: Dict[str, Dict[str, Union[str, Path]]],
        comparison_pairs: Optional[List[Tuple[str, str]]] = None,
        **apa_params,
    ) -> Dict[str, Any]:
        """
        Perform pairwise APA comparison between cell types

        Args:
            cell_type_data: Dict of cell_type -> {"cool_file": path, "loops_file": path}
            comparison_pairs: List of (cell_type1, cell_type2) tuples
            **apa_params: APA analysis parameters

        Returns:
            Dictionary with comparison results
        """
        logger.info("Starting pairwise APA comparison analysis")

        # Generate all pairs if not specified
        if comparison_pairs is None:
            cell_types = list(cell_type_data.keys())
            comparison_pairs = list(combinations(cell_types, 2))

        # Load loops for each cell type
        cell_type_loops = {}
        for cell_type, data in cell_type_data.items():
            loops_file = data["loops_file"]
            loops_df = self._load_loops_file(loops_file)
            cell_type_loops[cell_type] = self._loops_to_set(loops_df)
            logger.info(
                f"Loaded {len(cell_type_loops[cell_type])} loops for {cell_type}"
            )

        results = {}

        # Process each comparison pair
        for cell_type1, cell_type2 in comparison_pairs:
            pair_name = f"{cell_type1}_vs_{cell_type2}"
            logger.info(f"Processing comparison: {pair_name}")

            # Calculate loop sets
            loops1 = cell_type_loops[cell_type1]
            loops2 = cell_type_loops[cell_type2]

            common_loops = loops1 & loops2
            gained_loops = loops2 - loops1  # loops gained in cell_type2
            lost_loops = loops1 - loops2  # loops lost in cell_type2

            logger.info(f"{pair_name} statistics:")
            logger.info(f"  Common loops: {len(common_loops)}")
            logger.info(f"  Gained loops: {len(gained_loops)}")
            logger.info(f"  Lost loops: {len(lost_loops)}")

            # Perform APA analysis for each loop set
            pair_result = self._analyze_pairwise_loop_sets(
                pair_name,
                cell_type1,
                cell_type2,
                cell_type_data[cell_type1]["cool_file"],
                cell_type_data[cell_type2]["cool_file"],
                common_loops,
                gained_loops,
                lost_loops,
                **apa_params,
            )

            results[pair_name] = pair_result

        # Create comprehensive summary
        summary = self._create_pairwise_summary(results)

        # Save summary
        summary_file = self.output_dir / "pairwise_apa_summary.csv"
        summary.to_csv(summary_file, index=False)
        logger.info(f"Pairwise APA summary saved to {summary_file}")

        return {
            "pairwise_results": results,
            "summary": summary,
            "cell_type_loops": {k: len(v) for k, v in cell_type_loops.items()},
        }

    def _analyze_pairwise_loop_sets(
        self,
        pair_name: str,
        cell_type1: str,
        cell_type2: str,
        cool_file1: Union[str, Path],
        cool_file2: Union[str, Path],
        common_loops: Set,
        gained_loops: Set,
        lost_loops: Set,
        **apa_params,
    ) -> Dict[str, Any]:
        """Analyze loop sets for a single pairwise comparison"""

        loop_categories = {
            "common": common_loops,
            "gained": gained_loops,
            "lost": lost_loops,
        }

        # Process each loop category in both cell types
        apa_results = {}

        for category, loop_set in loop_categories.items():
            if len(loop_set) == 0:
                logger.warning(f"No loops in {category} category for {pair_name}")
                continue

            # Process in cell_type1
            result_key1 = f"{category}_in_{cell_type1}"
            apa_result1 = self._process_loop_set_apa(
                cool_file1, loop_set, result_key1, pair_name, **apa_params
            )
            apa_results[result_key1] = apa_result1

            # Process in cell_type2
            result_key2 = f"{category}_in_{cell_type2}"
            apa_result2 = self._process_loop_set_apa(
                cool_file2, loop_set, result_key2, pair_name, **apa_params
            )
            apa_results[result_key2] = apa_result2

        # Create visualization
        visualization_results = self._create_pairwise_visualizations(
            pair_name, cell_type1, cell_type2, apa_results, loop_categories
        )

        return {
            "pair_name": pair_name,
            "cell_type1": cell_type1,
            "cell_type2": cell_type2,
            "loop_counts": {
                "common": len(common_loops),
                "gained": len(gained_loops),
                "lost": len(lost_loops),
            },
            "apa_results": apa_results,
            "visualizations": visualization_results,
        }

    def _process_loop_set_apa(
        self,
        cool_file: Union[str, Path],
        loop_set: Set,
        result_key: str,
        pair_name: str,
        **apa_params,
    ) -> Optional[APAResult]:
        """Process APA for a single loop set"""

        if len(loop_set) == 0:
            return None

        # Create cache key
        cell_type = Path(cool_file).stem.split(".")[0]
        cache_key = f"{pair_name}_{result_key}"
        cache_file = self.cache_dir / f"{cache_key}.pkl"

        # Check cache
        try:
            cached_result = load_cached_apa(cache_file, apa_params)
            if cached_result:
                logger.debug(f"Using cached APA result for {result_key}")
                return cached_result
        except Exception as e:
            logger.debug(f"Could not load cached result: {e}")

        logger.info(f"Computing APA for {result_key}: {len(loop_set)} loops")

        try:
            # Load cooler
            clr = cooler.Cooler(str(cool_file))

            # Generate expected values
            cis_exp = cooltools.expected_cis(clr=clr, nproc=apa_params.get("nproc", 4))

            # Convert loop set to DataFrame
            loops_df = self._set_to_dataframe(loop_set)

            # Perform pileup
            pileup = coolpup.pileup(
                clr,
                loops_df,
                features_format="bedpe",
                expected_df=cis_exp,
                flank=apa_params.get("flank", 100000),
                min_diag=apa_params.get("min_diag", 3),
                maxdist=apa_params.get("maxdist", 1000000),
                nproc=apa_params.get("nproc", 4),
            )

            # Extract matrix
            apa_matrix = self._extract_apa_matrix(pileup)

            # Calculate scores
            center_score, corner_scores = calculate_apa_scores(apa_matrix)

            # Create result
            result = APAResult(
                sample_name=result_key,
                success=True,
                apa_matrix=apa_matrix,
                pileup_data=pileup,
                center_score=center_score,
                corner_scores=corner_scores,
                cache_file=cache_file,
                **apa_params,
            )

            # Cache result
            try:
                save_apa_cache(result, cache_file, apa_params)
            except Exception as e:
                logger.warning(f"Could not cache result: {e}")

            return result

        except Exception as e:
            logger.error(f"APA analysis failed for {result_key}: {e}")
            return None

    def _create_pairwise_visualizations(
        self,
        pair_name: str,
        cell_type1: str,
        cell_type2: str,
        apa_results: Dict[str, APAResult],
        loop_categories: Dict[str, Set],
    ) -> Dict[str, Path]:
        """Create visualizations for pairwise comparison"""

        visualization_files = {}

        # Create APA matrix comparison plot
        try:
            matrix_plot = self._create_apa_matrix_plot(
                pair_name, cell_type1, cell_type2, apa_results, loop_categories
            )
            visualization_files["matrix_plot"] = matrix_plot
        except Exception as e:
            logger.error(f"Failed to create matrix plot for {pair_name}: {e}")

        # Create violin/distribution plot
        try:
            violin_plot = self._create_apa_distribution_plot(
                pair_name, cell_type1, cell_type2, apa_results
            )
            visualization_files["violin_plot"] = violin_plot
        except Exception as e:
            logger.error(f"Failed to create violin plot for {pair_name}: {e}")

        return visualization_files

    def _create_apa_matrix_plot(
        self,
        pair_name: str,
        cell_type1: str,
        cell_type2: str,
        apa_results: Dict[str, APAResult],
        loop_categories: Dict[str, Set],
    ) -> Path:
        """Create APA matrix comparison plot"""

        # Create figure with 2 rows (cell types) x 3 columns (categories)
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))

        category_names = ["common", "gained", "lost"]

        for col, category in enumerate(category_names):
            # Skip if no loops in this category
            if len(loop_categories[category]) == 0:
                axes[0, col].text(0.5, 0.5, "No Data", ha="center", va="center")
                axes[1, col].text(0.5, 0.5, "No Data", ha="center", va="center")
                axes[0, col].set_title(f"{category.title()} loops")
                axes[1, col].set_title(f"{category.title()} loops")
                continue

            # Plot for cell_type1
            result_key1 = f"{category}_in_{cell_type1}"
            if result_key1 in apa_results and apa_results[result_key1]:
                matrix1 = apa_results[result_key1].apa_matrix
                self._plot_single_apa_matrix(
                    axes[0, col],
                    matrix1,
                    f"{category.title()} in {cell_type1}\n({len(loop_categories[category])} loops)",
                )

            # Plot for cell_type2
            result_key2 = f"{category}_in_{cell_type2}"
            if result_key2 in apa_results and apa_results[result_key2]:
                matrix2 = apa_results[result_key2].apa_matrix
                self._plot_single_apa_matrix(
                    axes[1, col],
                    matrix2,
                    f"{category.title()} in {cell_type2}\n({len(loop_categories[category])} loops)",
                )

        # Add colorbar
        if any(apa_results.values()):
            # Find a valid result for colorbar
            valid_result = next(
                r for r in apa_results.values() if r and r.apa_matrix is not None
            )
            im = axes[0, 0].images[0] if axes[0, 0].images else None

            if im:
                cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
                fig.colorbar(im, cax=cbar_ax)

        # Set overall title
        fig.suptitle(f"APA Comparison: {pair_name}", fontsize=16)

        # Adjust layout
        plt.tight_layout(rect=[0, 0.03, 0.9, 0.95])

        # Save figure
        output_file = self.output_dir / f"{pair_name}_apa_matrix_comparison.png"
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

        logger.info(f"APA matrix plot saved: {output_file}")
        return output_file

    def _create_apa_distribution_plot(
        self,
        pair_name: str,
        cell_type1: str,
        cell_type2: str,
        apa_results: Dict[str, APAResult],
    ) -> Path:
        """Create APA score distribution comparison plot"""

        # Extract center values for violin plot
        plot_data = []

        categories = ["common", "gained", "lost"]

        for category in categories:
            for cell_type in [cell_type1, cell_type2]:
                result_key = f"{category}_in_{cell_type}"

                if result_key in apa_results and apa_results[result_key]:
                    result = apa_results[result_key]

                    if result.apa_matrix is not None:
                        # Extract center region values
                        center_values = self._extract_center_values(
                            result.apa_matrix, size=5
                        )

                        for value in center_values:
                            plot_data.append(
                                {
                                    "Value": value,
                                    "Category": category.title(),
                                    "Cell_Type": cell_type,
                                }
                            )

        if not plot_data:
            logger.warning(f"No data for violin plot: {pair_name}")
            return None

        # Create DataFrame
        df = pd.DataFrame(plot_data)

        # Create violin plot
        plt.figure(figsize=(12, 6))

        # Create violin plot
        sns.violinplot(
            data=df,
            x="Category",
            y="Value",
            hue="Cell_Type",
            split=True,
            inner="quartile",
        )

        # Add reference line at y=1
        plt.axhline(y=1, color="black", linestyle="--", alpha=0.5)

        plt.title(f"APA Score Distribution: {pair_name}", fontsize=14)
        plt.xlabel("Loop Category", fontsize=12)
        plt.ylabel("APA Enrichment Value", fontsize=12)

        # Improve legend
        plt.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()

        # Save figure
        output_file = self.output_dir / f"{pair_name}_apa_distribution.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"APA distribution plot saved: {output_file}")
        return output_file

    def _plot_single_apa_matrix(
        self,
        ax,
        matrix: np.ndarray,
        title: str,
        vmin: float = 0,
        vmax: float = 3,
        show_scores: bool = True,
    ):
        """Plot a single APA matrix"""

        if matrix is None or matrix.size == 0:
            ax.text(0.5, 0.5, "No Data", ha="center", va="center")
            ax.set_title(title)
            return

        # Plot matrix
        im = ax.imshow(
            matrix,
            cmap="coolwarm",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            origin="lower",
        )

        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])

        if show_scores:
            # Calculate and display scores
            center_score, corner_scores = calculate_apa_scores(matrix)

            # Center score
            ax.text(
                matrix.shape[1] // 2,
                matrix.shape[0] // 2,
                f"{center_score:.2f}",
                ha="center",
                va="center",
                color="white",
                fontsize=10,
                fontweight="bold",
                bbox=dict(facecolor="black", alpha=0.5, boxstyle="round,pad=0.2"),
            )

            # Corner scores
            corner_size = 3
            positions = {
                "top_left": (corner_size // 2, corner_size // 2),
                "top_right": (matrix.shape[1] - corner_size // 2, corner_size // 2),
                "bottom_left": (corner_size // 2, matrix.shape[0] - corner_size // 2),
                "bottom_right": (
                    matrix.shape[1] - corner_size // 2,
                    matrix.shape[0] - corner_size // 2,
                ),
            }

            for corner, (x, y) in positions.items():
                if corner in corner_scores:
                    ax.text(
                        x,
                        y,
                        f"{corner_scores[corner]:.2f}",
                        ha="center",
                        va="center",
                        color="white",
                        fontsize=8,
                        fontweight="bold",
                        bbox=dict(
                            facecolor="black", alpha=0.5, boxstyle="round,pad=0.1"
                        ),
                    )

        return im

    def _extract_center_values(self, matrix: np.ndarray, size: int = 5) -> np.ndarray:
        """Extract center region values from APA matrix"""

        if matrix.size == 0:
            return np.array([])

        n_rows, n_cols = matrix.shape
        center_row = n_rows // 2
        center_col = n_cols // 2

        half = size // 2

        center_region = matrix[
            center_row - half : center_row + half + 1,
            center_col - half : center_col + half + 1,
        ]

        return center_region.flatten()

    def _create_pairwise_summary(self, results: Dict[str, Any]) -> pd.DataFrame:
        """Create summary DataFrame for pairwise comparisons"""

        summary_data = []

        for pair_name, result in results.items():
            base_data = {
                "comparison": pair_name,
                "cell_type1": result["cell_type1"],
                "cell_type2": result["cell_type2"],
                "common_loops": result["loop_counts"]["common"],
                "gained_loops": result["loop_counts"]["gained"],
                "lost_loops": result["loop_counts"]["lost"],
            }

            # Add APA scores for each category and cell type
            for category in ["common", "gained", "lost"]:
                for cell_type in [result["cell_type1"], result["cell_type2"]]:
                    result_key = f"{category}_in_{cell_type}"

                    if (
                        result_key in result["apa_results"]
                        and result["apa_results"][result_key]
                    ):
                        apa_result = result["apa_results"][result_key]
                        score_key = f"{category}_{cell_type}_apa_score"
                        base_data[score_key] = apa_result.center_score
                    else:
                        score_key = f"{category}_{cell_type}_apa_score"
                        base_data[score_key] = None

            summary_data.append(base_data)

        return pd.DataFrame(summary_data)

    def _load_loops_file(self, bedpe_file: Union[str, Path]) -> pd.DataFrame:
        """Load loops from BEDPE file"""
        return pd.read_csv(
            bedpe_file,
            sep="\t",
            header=None,
            names=[
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "score",
                "freq",
            ],
        )

    def _loops_to_set(self, loops_df: pd.DataFrame) -> Set[Tuple]:
        """Convert loops DataFrame to set of coordinate tuples"""
        return set(loops_df.iloc[:, :6].apply(lambda row: tuple(row), axis=1))

    def _set_to_dataframe(self, loop_set: Set[Tuple]) -> pd.DataFrame:
        """Convert set of loop tuples back to DataFrame"""
        return pd.DataFrame(
            list(loop_set),
            columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
        )

    def _extract_apa_matrix(self, pileup: pd.DataFrame) -> np.ndarray:
        """Extract APA matrix from pileup DataFrame"""
        try:
            mat_obj = pileup["data"].iloc[0]

            if isinstance(mat_obj, (list, np.ndarray)):
                return np.asarray(mat_obj, dtype=float)
            else:
                # Try to evaluate as string representation
                import ast

                return np.array(ast.literal_eval(str(mat_obj)), dtype=float)

        except Exception as e:
            logger.error(f"Failed to extract APA matrix: {e}")
            return np.array([])


class APAComparison:
    """Simple container for APA comparison results"""

    def __init__(
        self,
        comparison_name: str,
        cell_type1: str,
        cell_type2: str,
        common_loops: int,
        gained_loops: int,
        lost_loops: int,
    ):
        self.comparison_name = comparison_name
        self.cell_type1 = cell_type1
        self.cell_type2 = cell_type2
        self.common_loops = common_loops
        self.gained_loops = gained_loops
        self.lost_loops = lost_loops
        self.apa_results = {}
        self.visualization_files = []
