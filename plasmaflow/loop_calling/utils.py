"""
Utility functions for loop calling module
"""

import logging
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .caller import LoopCaller, LoopCallingResult
from .validation import validate_loops

logger = logging.getLogger(__name__)


def process_sample_batch(
    samples: Dict[str, str],
    config: Any,
    read_counts: Optional[Dict[str, int]] = None,
    max_workers: Optional[int] = None,
) -> Dict[str, LoopCallingResult]:
    """
    Process multiple samples in parallel

    Args:
        samples: Dictionary of sample_name -> cool_file_path
        config: PlasmaFlow configuration object
        read_counts: Optional read counts for each sample
        max_workers: Maximum number of parallel workers

    Returns:
        Dictionary of results for each sample
    """
    if max_workers is None:
        max_workers = min(len(samples), mp.cpu_count())

    logger.info(f"Processing {len(samples)} samples with {max_workers} workers")

    # For now, run sequentially since Peakachu can be memory intensive
    # TODO: Implement proper parallel processing with memory management

    loop_caller = LoopCaller(config)
    results = {}

    for sample_name, cool_file in samples.items():
        read_count = read_counts.get(sample_name) if read_counts else None
        result = loop_caller.call_loops_single_sample(
            sample_name=sample_name, cool_file=cool_file, read_count=read_count
        )
        results[sample_name] = result

    return results


def merge_loop_results(
    results_dict: Dict[str, LoopCallingResult], output_dir: Union[str, Path]
) -> pd.DataFrame:
    """
    Merge and summarize loop calling results

    Args:
        results_dict: Dictionary of sample results
        output_dir: Directory to save merged results

    Returns:
        Summary DataFrame
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Create summary
    summary_data = []
    all_loops = []

    for sample_name, result in results_dict.items():
        summary_data.append(
            {
                "sample": sample_name,
                "success": result.success,
                "num_loops": result.num_loops or 0,
                "execution_time_minutes": (
                    result.execution_time / 60 if result.execution_time else None
                ),
                "loops_file": str(result.loops_file) if result.loops_file else None,
                "error": result.error_message,
            }
        )

        # Load loops for successful results
        if result.success and result.loops_file and result.loops_file.exists():
            try:
                loops_df = pd.read_csv(result.loops_file, sep="\t", header=None)
                loops_df.columns = [
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "score",
                    "freq",
                ][: len(loops_df.columns)]
                loops_df["sample"] = sample_name
                all_loops.append(loops_df)
            except Exception as e:
                logger.warning(f"Could not load loops for {sample_name}: {e}")

    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_path / "loop_calling_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"Summary saved to {summary_file}")

    # Merge all loops
    if all_loops:
        merged_loops = pd.concat(all_loops, ignore_index=True)
        merged_file = output_path / "all_loops_merged.csv"
        merged_loops.to_csv(merged_file, index=False)
        logger.info(f"Merged loops saved to {merged_file}")

        # Create sample-specific BED files for downstream analysis
        for sample_name in merged_loops["sample"].unique():
            sample_loops = merged_loops[merged_loops["sample"] == sample_name].copy()

            # Convert to BED format (loop anchors)
            bed_data = []

            # Add first anchor
            bed1 = sample_loops[["chrom1", "start1", "end1"]].copy()
            bed1.columns = ["chrom", "start", "end"]
            bed1["type"] = "anchor1"
            bed1["loop_id"] = range(len(bed1))
            bed1["sample"] = sample_name

            # Add second anchor
            bed2 = sample_loops[["chrom2", "start2", "end2"]].copy()
            bed2.columns = ["chrom", "start", "end"]
            bed2["type"] = "anchor2"
            bed2["loop_id"] = range(len(bed2))
            bed2["sample"] = sample_name

            bed_merged = pd.concat([bed1, bed2], ignore_index=True)

            bed_file = output_path / f"{sample_name}_loop_anchors.bed"
            bed_merged[["chrom", "start", "end", "type", "sample"]].to_csv(
                bed_file, sep="\t", header=False, index=False
            )
            logger.debug(f"BED file saved for {sample_name}: {bed_file}")

    return summary_df


def create_loop_statistics(
    loops_files: List[Union[str, Path]], output_file: Union[str, Path]
) -> pd.DataFrame:
    """
    Create comprehensive statistics for loop calling results

    Args:
        loops_files: List of loop files to analyze
        output_file: Path to save statistics

    Returns:
        Statistics DataFrame
    """
    # Validate all loop files
    quality_metrics = validate_loops(loops_files)

    # Create statistics DataFrame
    stats_data = []

    for sample_name, metrics in quality_metrics.items():
        stats_data.append(
            {
                "sample": sample_name,
                "total_loops": metrics.num_loops,
                "intra_chromosomal": metrics.intra_chromosomal,
                "inter_chromosomal": metrics.inter_chromosomal,
                "intra_chromosomal_fraction": getattr(
                    metrics, "intra_chromosomal_fraction", 0
                ),
                "mean_distance_kb": (
                    metrics.mean_distance / 1000 if metrics.mean_distance else None
                ),
                "median_distance_kb": (
                    metrics.median_distance / 1000 if metrics.median_distance else None
                ),
                "min_distance_kb": (
                    metrics.min_distance / 1000 if metrics.min_distance else None
                ),
                "max_distance_kb": (
                    metrics.max_distance / 1000 if metrics.max_distance else None
                ),
                "mean_score": metrics.mean_score,
                "median_score": metrics.median_score,
                "has_scores": metrics.has_scores,
                "has_frequency": metrics.has_frequency,
                "valid_coordinates": metrics.valid_coordinates,
            }
        )

    stats_df = pd.DataFrame(stats_data)

    # Save statistics
    stats_df.to_csv(output_file, index=False)
    logger.info(f"Loop statistics saved to {output_file}")

    return stats_df


def filter_loops_by_distance(
    loops_file: Union[str, Path],
    output_file: Union[str, Path],
    min_distance: int = 20000,
    max_distance: int = 2000000,
) -> int:
    """
    Filter loops by genomic distance

    Args:
        loops_file: Input loops file
        output_file: Output filtered loops file
        min_distance: Minimum loop distance in bp
        max_distance: Maximum loop distance in bp

    Returns:
        Number of loops after filtering
    """
    # Load loops
    loops_df = pd.read_csv(loops_file, sep="\t", header=None)

    # Ensure we have enough columns
    if len(loops_df.columns) < 6:
        raise ValueError("Loops file must have at least 6 columns (BEDPE format)")

    loops_df.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"] + [
        f"col{i}" for i in range(6, len(loops_df.columns))
    ]

    # Filter intra-chromosomal loops by distance
    intra_mask = loops_df["chrom1"] == loops_df["chrom2"]

    # Calculate distances for intra-chromosomal loops
    intra_loops = loops_df[intra_mask].copy()
    if len(intra_loops) > 0:
        mid1 = (intra_loops["start1"] + intra_loops["end1"]) / 2
        mid2 = (intra_loops["start2"] + intra_loops["end2"]) / 2
        distances = abs(mid2 - mid1)

        distance_mask = (distances >= min_distance) & (distances <= max_distance)
        filtered_intra = intra_loops[distance_mask]
    else:
        filtered_intra = pd.DataFrame()

    # Keep all inter-chromosomal loops
    inter_loops = loops_df[~intra_mask]

    # Combine filtered results
    filtered_loops = pd.concat([filtered_intra, inter_loops], ignore_index=True)

    # Save filtered loops
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    filtered_loops.to_csv(output_file, sep="\t", header=False, index=False)

    logger.info(
        f"Filtered loops: {len(loops_df)} -> {len(filtered_loops)} "
        f"(distance range: {min_distance}-{max_distance} bp)"
    )

    return len(filtered_loops)
