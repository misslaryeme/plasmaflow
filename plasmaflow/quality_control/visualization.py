"""
Visualization functions for APA analysis
"""

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from .apa import APAResult

logger = logging.getLogger(__name__)


class APAPlotter:
    """Class for creating APA visualization plots"""

    def __init__(
        self, output_dir: Union[str, Path], figsize_per_plot: Tuple[int, int] = (6, 5)
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.figsize_per_plot = figsize_per_plot

    def plot_single_apa(
        self,
        result: APAResult,
        vmin: float = 0,
        vmax: float = 3,
        cmap: str = "coolwarm",
        show_scores: bool = True,
        title: Optional[str] = None,
    ) -> plt.Figure:
        """
        Plot APA matrix for a single sample

        Args:
            result: APAResult object
            vmin: Minimum value for colormap
            vmax: Maximum value for colormap
            cmap: Colormap name
            show_scores: Whether to show APA scores on plot
            title: Optional custom title

        Returns:
            matplotlib Figure object
        """
        if not result.success or result.apa_matrix is None:
            raise ValueError(
                f"Cannot plot APA for unsuccessful result: {result.sample_name}"
            )

        fig, ax = plt.subplots(1, 1, figsize=self.figsize_per_plot)

        matrix = result.apa_matrix

        # Create the heatmap
        im = ax.imshow(
            matrix,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            origin="lower",
        )

        # Set title
        if title:
            ax.set_title(title)
        else:
            ax.set_title(f"APA: {result.sample_name}")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("APA Score", rotation=270, labelpad=15)

        # Remove tick labels
        ax.set_xticks([])
        ax.set_yticks([])

        if show_scores and result.center_score is not None:
            # Add center score
            ax.text(
                matrix.shape[1] // 2,
                matrix.shape[0] // 2,
                f"{result.center_score:.2f}",
                ha="center",
                va="center",
                color="white",
                fontsize=12,
                fontweight="bold",
                bbox=dict(facecolor="black", alpha=0.7, boxstyle="round,pad=0.3"),
            )

            # Add corner scores if available
            if result.corner_scores:
                corner_size = 3

                positions = {
                    "top_left": (corner_size // 2, corner_size // 2),
                    "top_right": (matrix.shape[1] - corner_size // 2, corner_size // 2),
                    "bottom_left": (
                        corner_size // 2,
                        matrix.shape[0] - corner_size // 2,
                    ),
                    "bottom_right": (
                        matrix.shape[1] - corner_size // 2,
                        matrix.shape[0] - corner_size // 2,
                    ),
                }

                for corner, (x, y) in positions.items():
                    if corner in result.corner_scores:
                        ax.text(
                            x,
                            y,
                            f"{result.corner_scores[corner]:.2f}",
                            ha="center",
                            va="center",
                            color="white",
                            fontsize=10,
                            fontweight="bold",
                            bbox=dict(
                                facecolor="black", alpha=0.5, boxstyle="round,pad=0.2"
                            ),
                        )

        plt.tight_layout()
        return fig

    def plot_apa_grid(
        self,
        results: Dict[str, APAResult],
        vmin: float = 0,
        vmax: float = 3,
        cmap: str = "coolwarm",
        show_scores: bool = True,
        ncols: Optional[int] = None,
        save_path: Optional[Union[str, Path]] = None,
    ) -> plt.Figure:
        """
        Plot APA matrices for multiple samples in a grid

        Args:
            results: Dictionary of sample_name -> APAResult
            vmin: Minimum value for colormap
            vmax: Maximum value for colormap
            cmap: Colormap name
            show_scores: Whether to show APA scores
            ncols: Number of columns (auto-calculated if None)
            save_path: Optional path to save figure

        Returns:
            matplotlib Figure object
        """
        # Filter successful results
        successful_results = {
            k: v for k, v in results.items() if v.success and v.apa_matrix is not None
        }

        if not successful_results:
            raise ValueError("No successful APA results to plot")

        n_plots = len(successful_results)

        # Calculate grid dimensions
        if ncols is None:
            ncols = int(math.ceil(math.sqrt(n_plots)))
        nrows = int(math.ceil(n_plots / ncols))

        # Create figure
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(
                self.figsize_per_plot[0] * ncols,
                self.figsize_per_plot[1] * nrows,
            ),
        )

        if n_plots == 1:
            axes = np.array([axes])
        else:
            axes = np.array(axes).flatten()

        # Plot each sample
        for i, (sample_name, result) in enumerate(successful_results.items()):
            ax = axes[i]
            matrix = result.apa_matrix

            # Create heatmap
            im = ax.imshow(
                matrix,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
                origin="lower",
            )

            ax.set_title(f"APA: {sample_name}")
            ax.set_xticks([])
            ax.set_yticks([])

            # Add colorbar to each subplot
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            if show_scores and result.center_score is not None:
                # Center score
                ax.text(
                    matrix.shape[1] // 2,
                    matrix.shape[0] // 2,
                    f"{result.center_score:.2f}",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=12,
                    fontweight="bold",
                    bbox=dict(facecolor="black", alpha=0.7, boxstyle="round,pad=0.3"),
                )

                # Corner scores
                if result.corner_scores:
                    corner_size = 3
                    positions = {
                        "top_left": (corner_size // 2, corner_size // 2),
                        "top_right": (
                            matrix.shape[1] - corner_size // 2,
                            corner_size // 2,
                        ),
                        "bottom_left": (
                            corner_size // 2,
                            matrix.shape[0] - corner_size // 2,
                        ),
                        "bottom_right": (
                            matrix.shape[1] - corner_size // 2,
                            matrix.shape[0] - corner_size // 2,
                        ),
                    }

                    for corner, (x, y) in positions.items():
                        if corner in result.corner_scores:
                            ax.text(
                                x,
                                y,
                                f"{result.corner_scores[corner]:.2f}",
                                ha="center",
                                va="center",
                                color="white",
                                fontsize=10,
                                fontweight="bold",
                                bbox=dict(
                                    facecolor="black",
                                    alpha=0.5,
                                    boxstyle="round,pad=0.2",
                                ),
                            )

        # Turn off unused subplots
        for i in range(n_plots, len(axes)):
            axes[i].axis("off")

        plt.tight_layout()

        # Save if requested
        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info(f"APA grid plot saved to {save_path}")

        return fig

    def plot_apa_comparison_bar(
        self,
        results: Dict[str, APAResult],
        metric: str = "center",
        title: Optional[str] = None,
        save_path: Optional[Union[str, Path]] = None,
    ) -> plt.Figure:
        """
        Create bar plot comparing APA scores across samples

        Args:
            results: Dictionary of sample_name -> APAResult
            metric: Which metric to plot ('center', 'top_left', etc.)
            title: Optional custom title
            save_path: Optional path to save figure

        Returns:
            matplotlib Figure object
        """
        # Extract scores
        sample_names = []
        scores = []

        for sample_name, result in results.items():
            if not result.success:
                continue

            if metric == "center" and result.center_score is not None:
                sample_names.append(sample_name)
                scores.append(result.center_score)
            elif metric in result.corner_scores:
                sample_names.append(sample_name)
                scores.append(result.corner_scores[metric])

        if not scores:
            raise ValueError(f"No valid {metric} scores found")

        # Create bar plot
        fig, ax = plt.subplots(1, 1, figsize=(max(8, len(scores) * 0.8), 6))

        bars = ax.bar(sample_names, scores, alpha=0.8)

        # Customize plot
        ax.set_ylabel(f"APA Score ({metric})")
        ax.set_xlabel("Sample")

        if title:
            ax.set_title(title)
        else:
            ax.set_title(f'APA {metric.replace("_", " ").title()} Score Comparison')

        # Add value labels on bars
        for bar, score in zip(bars, scores):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{score:.3f}",
                ha="center",
                va="bottom",
                fontweight="bold",
            )

        # Rotate x labels if many samples
        if len(sample_names) > 4:
            plt.xticks(rotation=45, ha="right")

        plt.tight_layout()

        # Save if requested
        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info(f"APA comparison plot saved to {save_path}")

        return fig


def create_apa_plots(
    results: Dict[str, APAResult],
    output_dir: Union[str, Path],
    plot_types: List[str] = ["grid", "individual", "comparison"],
    **kwargs,
) -> Dict[str, Path]:
    """
    Create comprehensive APA plots

    Args:
        results: Dictionary of APAResult objects
        output_dir: Directory to save plots
        plot_types: Types of plots to create
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        Dictionary of plot_type -> saved_file_path
    """
    plotter = APAPlotter(output_dir)
    saved_plots = {}

    # Grid plot
    if "grid" in plot_types:
        try:
            grid_path = Path(output_dir) / "apa_grid_plot.png"
            fig = plotter.plot_apa_grid(results, save_path=grid_path, **kwargs)
            plt.close(fig)
            saved_plots["grid"] = grid_path
        except Exception as e:
            logger.error(f"Failed to create grid plot: {e}")

    # Individual plots
    if "individual" in plot_types:
        individual_dir = Path(output_dir) / "individual_plots"
        individual_dir.mkdir(parents=True, exist_ok=True)

        for sample_name, result in results.items():
            if result.success:
                try:
                    fig = plotter.plot_single_apa(result, **kwargs)
                    plot_path = individual_dir / f"{sample_name}_apa.png"
                    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
                    plt.close(fig)
                except Exception as e:
                    logger.error(
                        f"Failed to create individual plot for {sample_name}: {e}"
                    )

        saved_plots["individual"] = individual_dir

    # Comparison bar plot
    if "comparison" in plot_types:
        try:
            comparison_path = Path(output_dir) / "apa_score_comparison.png"
            fig = plotter.plot_apa_comparison_bar(results, save_path=comparison_path)
            plt.close(fig)
            saved_plots["comparison"] = comparison_path
        except Exception as e:
            logger.error(f"Failed to create comparison plot: {e}")

    return saved_plots
