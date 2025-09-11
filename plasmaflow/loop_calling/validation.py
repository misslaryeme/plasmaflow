"""
Loop validation and quality control
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class LoopQualityMetrics:
    """Quality metrics for loop calling results"""

    sample_name: str
    num_loops: int

    # Distance statistics
    mean_distance: float
    median_distance: float
    min_distance: float
    max_distance: float

    # Score statistics
    mean_score: Optional[float] = None
    median_score: Optional[float] = None
    min_score: Optional[float] = None
    max_score: Optional[float] = None

    # Chromosome distribution
    intra_chromosomal: int = 0
    inter_chromosomal: int = 0

    # Quality flags
    has_scores: bool = False
    has_frequency: bool = False
    valid_coordinates: bool = True

    def __post_init__(self):
        """Calculate derived metrics"""
        if self.num_loops > 0:
            self.intra_chromosomal_fraction = self.intra_chromosomal / self.num_loops
            self.inter_chromosomal_fraction = self.inter_chromosomal / self.num_loops
        else:
            self.intra_chromosomal_fraction = 0.0
            self.inter_chromosomal_fraction = 0.0


def load_loops_file(loops_file: Union[str, Path]) -> pd.DataFrame:
    """
    Load loops from BEDPE file with flexible column handling

    Args:
        loops_file: Path to BEDPE file

    Returns:
        DataFrame with standardized column names
    """
    loops_path = Path(loops_file)

    if not loops_path.exists():
        raise FileNotFoundError(f"Loops file not found: {loops_path}")

    # Try to detect file format
    with open(loops_path, "r") as f:
        first_line = f.readline().strip()

    # Check if header exists
    has_header = not first_line.split("\t")[0].startswith("chr")

    try:
        if has_header:
            loops_df = pd.read_csv(loops_path, sep="\t")
        else:
            # Standard BEDPE format
            cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]

            # Read first line to determine number of columns
            sample_df = pd.read_csv(loops_path, sep="\t", nrows=1, header=None)
            n_cols = len(sample_df.columns)

            if n_cols > 6:
                # Add additional columns
                additional_cols = [f"col{i}" for i in range(7, n_cols + 1)]
                cols.extend(additional_cols)

                # Common additional columns
                if n_cols >= 7:
                    cols[6] = "score"
                if n_cols >= 8:
                    cols[7] = "freq"

            loops_df = pd.read_csv(loops_path, sep="\t", header=None, names=cols)

    except Exception as e:
        logger.error(f"Error loading loops file {loops_path}: {e}")
        raise

    # Standardize column names if needed
    column_mapping = {
        "chr1": "chrom1",
        "chr2": "chrom2",
        "x1": "start1",
        "x2": "start2",
        "y1": "end1",
        "y2": "end2",
    }

    loops_df = loops_df.rename(columns=column_mapping)

    # Ensure required columns exist
    required_cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    missing_cols = [col for col in required_cols if col not in loops_df.columns]

    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return loops_df


def calculate_loop_distances(loops_df: pd.DataFrame) -> pd.Series:
    """Calculate distances between loop anchors"""

    # Only calculate for intra-chromosomal loops
    intra_mask = loops_df["chrom1"] == loops_df["chrom2"]
    intra_loops = loops_df[intra_mask].copy()

    if len(intra_loops) == 0:
        return pd.Series(dtype=float)

    # Calculate distance as midpoint to midpoint
    mid1 = (intra_loops["start1"] + intra_loops["end1"]) / 2
    mid2 = (intra_loops["start2"] + intra_loops["end2"]) / 2

    distances = np.abs(mid2 - mid1)

    return distances


def validate_loop_coordinates(loops_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate loop coordinates for common issues

    Returns:
        Dictionary with validation results
    """
    results = {"valid": True, "issues": [], "warnings": []}

    # Check for negative coordinates
    coord_cols = ["start1", "end1", "start2", "end2"]
    for col in coord_cols:
        if col in loops_df.columns:
            negative_coords = (loops_df[col] < 0).sum()
            if negative_coords > 0:
                results["valid"] = False
                results["issues"].append(
                    f"Found {negative_coords} negative coordinates in {col}"
                )

    # Check for invalid intervals (start > end)
    if "start1" in loops_df.columns and "end1" in loops_df.columns:
        invalid_intervals1 = (loops_df["start1"] > loops_df["end1"]).sum()
        if invalid_intervals1 > 0:
            results["valid"] = False
            results["issues"].append(
                f"Found {invalid_intervals1} invalid intervals (start1 > end1)"
            )

    if "start2" in loops_df.columns and "end2" in loops_df.columns:
        invalid_intervals2 = (loops_df["start2"] > loops_df["end2"]).sum()
        if invalid_intervals2 > 0:
            results["valid"] = False
            results["issues"].append(
                f"Found {invalid_intervals2} invalid intervals (start2 > end2)"
            )

    # Check for very small or very large intervals
    if "start1" in loops_df.columns and "end1" in loops_df.columns:
        interval_sizes1 = loops_df["end1"] - loops_df["start1"]
        very_small = (interval_sizes1 < 1000).sum()  # < 1kb
        very_large = (interval_sizes1 > 100000).sum()  # > 100kb

        if very_small > 0:
            results["warnings"].append(
                f"Found {very_small} very small intervals (< 1kb) in anchor 1"
            )
        if very_large > 0:
            results["warnings"].append(
                f"Found {very_large} very large intervals (> 100kb) in anchor 1"
            )

    # Check chromosome naming consistency
    chromosomes = set(loops_df["chrom1"].unique()) | set(loops_df["chrom2"].unique())
    chr_formats = []
    for chrom in chromosomes:
        if str(chrom).startswith("chr"):
            chr_formats.append("chr_prefix")
        else:
            chr_formats.append("no_prefix")

    if len(set(chr_formats)) > 1:
        results["warnings"].append(
            "Inconsistent chromosome naming (mixed 'chr' prefix)"
        )

    return results


def calculate_quality_metrics(
    loops_file: Union[str, Path], sample_name: Optional[str] = None
) -> LoopQualityMetrics:
    """
    Calculate comprehensive quality metrics for loops

    Args:
        loops_file: Path to loops BEDPE file
        sample_name: Name of sample (inferred from filename if not provided)

    Returns:
        LoopQualityMetrics object
    """
    loops_path = Path(loops_file)

    if sample_name is None:
        sample_name = loops_path.stem

    # Load loops
    loops_df = load_loops_file(loops_path)
    num_loops = len(loops_df)

    if num_loops == 0:
        logger.warning(f"No loops found in {loops_path}")
        return LoopQualityMetrics(
            sample_name=sample_name,
            num_loops=0,
            mean_distance=0,
            median_distance=0,
            min_distance=0,
            max_distance=0,
        )

    # Calculate distances
    distances = calculate_loop_distances(loops_df)

    if len(distances) == 0:
        # Only inter-chromosomal loops
        distance_stats = {
            "mean_distance": np.nan,
            "median_distance": np.nan,
            "min_distance": np.nan,
            "max_distance": np.nan,
        }
    else:
        distance_stats = {
            "mean_distance": float(distances.mean()),
            "median_distance": float(distances.median()),
            "min_distance": float(distances.min()),
            "max_distance": float(distances.max()),
        }

    # Score statistics
    score_stats = {}
    if "score" in loops_df.columns:
        scores = loops_df["score"].dropna()
        if len(scores) > 0:
            score_stats = {
                "mean_score": float(scores.mean()),
                "median_score": float(scores.median()),
                "min_score": float(scores.min()),
                "max_score": float(scores.max()),
                "has_scores": True,
            }
        else:
            score_stats = {"has_scores": False}
    else:
        score_stats = {"has_scores": False}

    # Chromosome distribution
    intra_chromosomal = (loops_df["chrom1"] == loops_df["chrom2"]).sum()
    inter_chromosomal = num_loops - intra_chromosomal

    # Validate coordinates
    validation = validate_loop_coordinates(loops_df)

    return LoopQualityMetrics(
        sample_name=sample_name,
        num_loops=num_loops,
        intra_chromosomal=intra_chromosomal,
        inter_chromosomal=inter_chromosomal,
        has_frequency="freq" in loops_df.columns,
        valid_coordinates=validation["valid"],
        **distance_stats,
        **score_stats,
    )


def validate_loops(
    loops_files: Union[str, Path, List[Union[str, Path]]],
    output_file: Optional[Union[str, Path]] = None,
) -> Dict[str, LoopQualityMetrics]:
    """
    Validate multiple loops files and generate quality report

    Args:
        loops_files: Single file or list of loop files to validate
        output_file: Optional path to save quality report

    Returns:
        Dictionary mapping sample names to quality metrics
    """
    if isinstance(loops_files, (str, Path)):
        loops_files = [loops_files]

    results = {}

    for loops_file in loops_files:
        try:
            metrics = calculate_quality_metrics(loops_file)
            results[metrics.sample_name] = metrics

            logger.info(
                f"Sample {metrics.sample_name}: {metrics.num_loops} loops, "
                f"mean distance: {metrics.mean_distance:.0f} bp"
            )

        except Exception as e:
            logger.error(f"Failed to validate {loops_file}: {e}")
            continue

    # Save report if requested
    if output_file and results:
        save_quality_report(results, output_file)

    return results


def save_quality_report(
    quality_metrics: Dict[str, LoopQualityMetrics], output_file: Union[str, Path]
) -> None:
    """Save quality metrics to CSV file"""

    report_data = []
    for sample_name, metrics in quality_metrics.items():
        report_data.append(
            {
                "sample": metrics.sample_name,
                "num_loops": metrics.num_loops,
                "intra_chromosomal": metrics.intra_chromosomal,
                "inter_chromosomal": metrics.inter_chromosomal,
                "intra_chromosomal_fraction": getattr(
                    metrics, "intra_chromosomal_fraction", None
                ),
                "mean_distance": metrics.mean_distance,
                "median_distance": metrics.median_distance,
                "min_distance": metrics.min_distance,
                "max_distance": metrics.max_distance,
                "mean_score": metrics.mean_score,
                "median_score": metrics.median_score,
                "has_scores": metrics.has_scores,
                "has_frequency": metrics.has_frequency,
                "valid_coordinates": metrics.valid_coordinates,
            }
        )

    report_df = pd.DataFrame(report_data)
    report_df.to_csv(output_file, index=False)

    logger.info(f"Quality report saved to {output_file}")
