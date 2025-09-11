"""
Matrix generation functionality for PlasmaFlow
"""

import logging
import os
import shlex
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

from ..config import Config
from ..utils import get_logger
from .bed_utils import create_anchor_beds, parse_loop_coordinates
from .deeptools import DeepToolsInterface

logger = get_logger(__name__)


class MatrixGenerator:
    """Generator for deepTools matrices from differential loop results"""

    def __init__(self, config: Config):
        """
        Initialize matrix generator

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config
        self.visualization_params = config.visualization.get("heatmaps", {})

        # Create output directories
        self.output_dir = Path(config.output_dir) / "visualization"
        self.bed_output_dir = self.output_dir / "bed_files"
        self.matrix_output_dir = self.output_dir / "matrices"

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.bed_output_dir.mkdir(parents=True, exist_ok=True)
        self.matrix_output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize deepTools interface
        self.deeptools = DeepToolsInterface(config)

        # Default parameters
        self.default_params = {
            "flanking_region": 3000,
            "bin_size": 100,
            "reference_point": "center",
            "skip_zeros": True,
            "sort_regions": "keep",
            "n_processors": config.n_threads,
        }

    def generate_matrix_from_differential_results(
        self,
        differential_results: Dict[str, Any],
        bigwig_files: Dict[str, str],
        comparison_name: str,
        method_name: str = "differential",
    ) -> Dict[str, Any]:
        """
        Generate deepTools matrix from differential analysis results

        Args:
            differential_results: Dictionary with 'up_regulated', 'down_regulated', 'stable' results
            bigwig_files: Dictionary mapping sample names to BigWig file paths
            comparison_name: Name of the comparison (e.g., "prePB_vs_memB")
            method_name: Name of the analysis method

        Returns:
            Dictionary with generated files and metadata
        """
        logger.info(f"Generating matrix for {comparison_name} using {method_name}")

        # Step 1: Create BED files from differential results
        bed_files = self._create_bed_files_from_results(
            differential_results, comparison_name, method_name
        )

        if not bed_files:
            raise ValueError("No BED files were created from differential results")

        # Step 2: Validate BigWig files
        validated_bigwigs = self.deeptools.validate_bigwig_files(bigwig_files)

        if not validated_bigwigs:
            raise FileNotFoundError("No valid BigWig files found")

        # Step 3: Generate matrix using deepTools
        matrix_result = self.deeptools.run_compute_matrix(
            bed_files=bed_files,
            bigwig_files=validated_bigwigs,
            output_prefix=f"matrix_{comparison_name}_{method_name}",
            output_dir=self.matrix_output_dir,
        )

        # Combine results
        result = {
            "bed_files": bed_files,
            "bigwig_files": validated_bigwigs,
            "matrix_file": matrix_result["matrix_file"],
            "regions_file": matrix_result["regions_file"],
            "region_labels": matrix_result["region_labels"],
            "sample_labels": list(validated_bigwigs.keys()),
            "comparison_name": comparison_name,
            "method_name": method_name,
        }

        logger.info(f"Matrix generation completed: {matrix_result['matrix_file']}")
        return result

    def _create_bed_files_from_results(
        self,
        differential_results: Dict[str, Any],
        comparison_name: str,
        method_name: str,
    ) -> Dict[str, Path]:
        """Create BED files from differential analysis results"""

        logger.info("Creating BED files from differential results")

        bed_files = {}
        categories_processed = 0

        # Process each category (up, down, stable/common)
        category_mapping = {
            "up_regulated": "up",
            "down_regulated": "down",
            "stable": "common",
        }

        for result_key, category in category_mapping.items():
            if result_key not in differential_results:
                continue

            results_df = differential_results[result_key]

            if results_df is None or len(results_df) == 0:
                logger.warning(f"No results found for {category} category")
                continue

            logger.info(
                f"Processing {category.upper()} category: {len(results_df)} interactions"
            )

            # Create BED files for this category
            category_bed_files = create_anchor_beds(
                results_df,
                category_name=category,
                comparison_name=comparison_name,
                output_dir=self.bed_output_dir,
            )

            bed_files.update(category_bed_files)
            categories_processed += 1

        logger.info(f"Created BED files for {categories_processed} categories")
        return bed_files

    def generate_method_comparison_matrix(
        self,
        method_results: Dict[str, Dict[str, Any]],
        bigwig_files: Dict[str, str],
        comparison_name: str,
    ) -> Dict[str, Any]:
        """
        Generate matrices comparing different methods (diffHic, DESeq2, edgeR)

        Args:
            method_results: Dictionary of method_name -> differential results
            bigwig_files: BigWig files for matrix generation
            comparison_name: Name of comparison

        Returns:
            Dictionary with comparison matrix results
        """
        logger.info(f"Generating method comparison matrix for {comparison_name}")

        all_bed_files = {}
        method_stats = {}

        # Create BED files for each method
        for method_name, results in method_results.items():
            logger.info(f"Processing method: {method_name}")

            method_bed_files = self._create_bed_files_from_results(
                results, comparison_name, method_name
            )

            # Add method prefix to BED file keys
            for key, bed_file in method_bed_files.items():
                prefixed_key = f"{method_name}_{key}"
                all_bed_files[prefixed_key] = bed_file

            # Calculate statistics
            method_stats[method_name] = {
                "up_regulated": len(results.get("up_regulated", [])),
                "down_regulated": len(results.get("down_regulated", [])),
                "stable": len(results.get("stable", [])),
            }

        if not all_bed_files:
            raise ValueError("No BED files created for method comparison")

        # Generate combined matrix
        validated_bigwigs = self.deeptools.validate_bigwig_files(bigwig_files)

        matrix_result = self.deeptools.run_compute_matrix(
            bed_files=all_bed_files,
            bigwig_files=validated_bigwigs,
            output_prefix=f"matrix_{comparison_name}_method_comparison",
            output_dir=self.matrix_output_dir,
        )

        result = {
            "bed_files": all_bed_files,
            "bigwig_files": validated_bigwigs,
            "matrix_file": matrix_result["matrix_file"],
            "regions_file": matrix_result["regions_file"],
            "region_labels": matrix_result["region_labels"],
            "sample_labels": list(validated_bigwigs.keys()),
            "method_stats": method_stats,
            "comparison_name": comparison_name,
        }

        logger.info(
            f"Method comparison matrix generated: {matrix_result['matrix_file']}"
        )
        return result

    def create_summary_report(self, matrix_results: Dict[str, Any]) -> pd.DataFrame:
        """Create summary report of matrix generation results"""

        summary_data = []

        for result_name, result in matrix_results.items():
            summary_data.append(
                {
                    "analysis": result_name,
                    "comparison": result.get("comparison_name", "Unknown"),
                    "method": result.get("method_name", "Unknown"),
                    "matrix_file": str(result.get("matrix_file", "")),
                    "n_bed_files": len(result.get("bed_files", {})),
                    "n_bigwig_files": len(result.get("bigwig_files", {})),
                    "region_labels": ", ".join(result.get("region_labels", [])),
                    "sample_labels": ", ".join(result.get("sample_labels", [])),
                }
            )

        summary_df = pd.DataFrame(summary_data)

        # Save summary
        summary_file = self.output_dir / "matrix_generation_summary.csv"
        summary_df.to_csv(summary_file, index=False)

        logger.info(f"Matrix generation summary saved to {summary_file}")
        return summary_df


def create_bed_files(
    differential_results: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    output_dir: Union[str, Path],
    comparison_name: str = "comparison",
) -> Dict[str, Path]:
    """
    Create BED files from differential analysis results

    Args:
        differential_results: DataFrame or dict of DataFrames with loop results
        output_dir: Output directory for BED files
        comparison_name: Name for the comparison

    Returns:
        Dictionary mapping category names to BED file paths
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    bed_files = {}

    if isinstance(differential_results, pd.DataFrame):
        # Single DataFrame - create based on regulation column
        results_dict = {}

        if "regulation" in differential_results.columns:
            for regulation_type in [
                "Up-regulated",
                "Down-regulated",
                "Not Significant",
            ]:
                subset = differential_results[
                    differential_results["regulation"] == regulation_type
                ]
                if len(subset) > 0:
                    category_name = (
                        regulation_type.lower().replace("-", "_").replace(" ", "_")
                    )
                    results_dict[category_name] = subset
        else:
            # No regulation column, treat as single category
            results_dict["all"] = differential_results

        differential_results = results_dict

    # Process each category
    for category, results_df in differential_results.items():
        if len(results_df) == 0:
            continue

        category_bed_files = create_anchor_beds(
            results_df,
            category_name=category,
            comparison_name=comparison_name,
            output_dir=output_path,
        )

        bed_files.update(category_bed_files)

    return bed_files


def generate_deeptools_matrix(
    bed_files: Dict[str, Union[str, Path]],
    bigwig_files: Dict[str, Union[str, Path]],
    output_file: Union[str, Path],
    flanking_region: int = 3000,
    reference_point: str = "center",
    n_processors: int = 4,
    **kwargs,
) -> Dict[str, Any]:
    """
    Generate deepTools matrix using computeMatrix

    Args:
        bed_files: Dictionary of region_name -> BED file path
        bigwig_files: Dictionary of sample_name -> BigWig file path
        output_file: Output matrix file path
        flanking_region: Flanking region size in bp
        reference_point: Reference point for matrix ("center", "TSS", "TES")
        n_processors: Number of processors to use
        **kwargs: Additional computeMatrix parameters

    Returns:
        Dictionary with execution results
    """
    # Create temporary DeepToolsInterface
    from ..config import get_default_config

    config = get_default_config()
    deeptools = DeepToolsInterface(config)

    # Validate files
    validated_bed_files = {}
    for name, bed_file in bed_files.items():
        bed_path = Path(bed_file)
        if bed_path.exists():
            validated_bed_files[name] = bed_path
        else:
            logger.warning(f"BED file not found: {bed_file}")

    validated_bigwig_files = deeptools.validate_bigwig_files(bigwig_files)

    if not validated_bed_files:
        raise FileNotFoundError("No valid BED files found")

    if not validated_bigwig_files:
        raise FileNotFoundError("No valid BigWig files found")

    # Generate matrix
    output_path = Path(output_file)
    output_dir = output_path.parent
    output_prefix = output_path.stem

    return deeptools.run_compute_matrix(
        bed_files=validated_bed_files,
        bigwig_files=validated_bigwig_files,
        output_prefix=output_prefix,
        output_dir=output_dir,
        flanking_region=flanking_region,
        reference_point=reference_point,
        n_processors=n_processors,
        **kwargs,
    )
