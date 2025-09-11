"""
Main differential analysis coordinator
"""

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..config import Config
from ..utils import RInterface, get_logger, setup_logging
from .comparison import ComparisonManager
from .filtering import FilterManager
from .methods import DESeq2Analyzer, DiffHicAnalyzer, EdgeRAnalyzer
from .results import ResultsProcessor

logger = get_logger(__name__)


@dataclass
class DifferentialResult:
    """Result of differential analysis"""

    comparison_name: str
    method: str
    success: bool

    # Results data
    results_table: Optional[pd.DataFrame] = None

    # Statistics
    n_tested: Optional[int] = None
    n_significant: Optional[int] = None
    n_up_regulated: Optional[int] = None
    n_down_regulated: Optional[int] = None

    # Quality metrics
    fdr_threshold: float = 0.05
    logfc_threshold: float = 1.0
    min_fdr: Optional[float] = None
    max_abs_logfc: Optional[float] = None

    # Files
    output_files: Dict[str, Path] = None

    # Execution info
    execution_time: Optional[float] = None
    error_message: Optional[str] = None

    def __post_init__(self):
        if self.output_files is None:
            self.output_files = {}


class DifferentialAnalyzer:
    """Main class for differential chromatin loop analysis"""

    def __init__(self, config: Config):
        self.config = config
        self.diff_params = config.differential

        # Create output directory
        self.output_dir = Path(config.output_dir) / "differential_analysis"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.r_interface = RInterface(config.r_config)
        self.comparison_manager = ComparisonManager(config)
        self.filter_manager = FilterManager(config)
        self.results_processor = ResultsProcessor(config)

        # Initialize method analyzers
        self.analyzers = {
            "diffHic": DiffHicAnalyzer(config, self.r_interface),
            "DESeq2": DESeq2Analyzer(config, self.r_interface),
            "edgeR": EdgeRAnalyzer(config, self.r_interface),
        }

    def run_differential_analysis(
        self,
        counts_matrix: Union[str, Path, pd.DataFrame],
        metadata: Union[str, Path, pd.DataFrame],
        comparisons: Optional[List[Dict[str, str]]] = None,
        methods: Optional[List[str]] = None,
    ) -> Dict[str, DifferentialResult]:
        """
        Run differential analysis with multiple methods

        Args:
            counts_matrix: Path to counts matrix CSV or DataFrame
            metadata: Path to metadata CSV or DataFrame
            comparisons: List of comparison dictionaries
            methods: List of methods to use

        Returns:
            Dictionary of comparison_method -> DifferentialResult
        """
        logger.info("Starting differential chromatin loop analysis")

        # Load data
        if isinstance(counts_matrix, (str, Path)):
            counts_df = pd.read_csv(counts_matrix)
        else:
            counts_df = counts_matrix.copy()

        if isinstance(metadata, (str, Path)):
            metadata_df = pd.read_csv(metadata)
        else:
            metadata_df = metadata.copy()

        # Use configured comparisons if not provided
        if comparisons is None:
            comparisons = self.diff_params.get("comparisons", [])

        if methods is None:
            methods = self.diff_params.get("methods", ["diffHic"])

        logger.info(
            f"Running {len(comparisons)} comparisons with {len(methods)} methods"
        )

        results = {}

        for comparison in comparisons:
            comp_name = comparison["name"]
            logger.info(f"Processing comparison: {comp_name}")

            # Prepare data for this comparison
            comp_data = self.comparison_manager.prepare_comparison_data(
                counts_df, metadata_df, comparison
            )

            if comp_data is None:
                logger.error(f"Failed to prepare data for comparison {comp_name}")
                continue

            # Run each method
            for method in methods:
                if method not in self.analyzers:
                    logger.warning(f"Unknown method: {method}")
                    continue

                result_key = f"{comp_name}_{method}"
                logger.info(f"Running {method} for {comp_name}")

                try:
                    result = self._run_single_analysis(method, comp_data, comparison)
                    results[result_key] = result

                    if result.success:
                        logger.info(
                            f"{method} completed: {result.n_significant} significant interactions"
                        )
                    else:
                        logger.error(f"{method} failed: {result.error_message}")

                except Exception as e:
                    logger.error(f"Error running {method} for {comp_name}: {e}")
                    results[result_key] = DifferentialResult(
                        comparison_name=comp_name,
                        method=method,
                        success=False,
                        error_message=str(e),
                    )

        # Create summary report
        self._create_analysis_summary(results)

        logger.info(
            f"Differential analysis completed. {len(results)} results generated."
        )

        return results

    def _run_single_analysis(
        self, method: str, comp_data: Dict[str, Any], comparison: Dict[str, str]
    ) -> DifferentialResult:
        """Run single differential analysis method"""

        start_time = time.time()

        try:
            # Apply filters
            filtered_data = self.filter_manager.apply_filters(comp_data)

            # Run analysis
            analyzer = self.analyzers[method]
            result = analyzer.run_analysis(filtered_data, comparison)

            if result.success:
                # Process results
                processed_result = self.results_processor.process_results(
                    result, comparison, method
                )

                # Calculate execution time
                processed_result.execution_time = time.time() - start_time

                # Save results
                output_files = self._save_results(processed_result, comparison, method)
                processed_result.output_files.update(output_files)

                return processed_result
            else:
                result.execution_time = time.time() - start_time
                return result

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            return DifferentialResult(
                comparison_name=comparison["name"],
                method=method,
                success=False,
                error_message=str(e),
                execution_time=time.time() - start_time,
            )

    def _save_results(
        self, result: DifferentialResult, comparison: Dict[str, str], method: str
    ) -> Dict[str, Path]:
        """Save analysis results to files"""

        output_files = {}
        base_filename = f"{comparison['name']}_{method}"

        # Main results file
        results_file = self.output_dir / f"{base_filename}_results.csv"
        result.results_table.to_csv(results_file, index=False)
        output_files["results"] = results_file

        # Significant results only
        if result.n_significant > 0:
            sig_results = result.results_table[
                result.results_table["regulation"] != "Not Significant"
            ]
            sig_file = self.output_dir / f"{base_filename}_significant.csv"
            sig_results.to_csv(sig_file, index=False)
            output_files["significant"] = sig_file

            # Up-regulated
            up_results = result.results_table[
                result.results_table["regulation"] == "Up-regulated"
            ]
            if len(up_results) > 0:
                up_file = self.output_dir / f"{base_filename}_up_regulated.csv"
                up_results.to_csv(up_file, index=False)
                output_files["up_regulated"] = up_file

            # Down-regulated
            down_results = result.results_table[
                result.results_table["regulation"] == "Down-regulated"
            ]
            if len(down_results) > 0:
                down_file = self.output_dir / f"{base_filename}_down_regulated.csv"
                down_results.to_csv(down_file, index=False)
                output_files["down_regulated"] = down_file

        # Stable/unchanged interactions
        stable_results = result.results_table[
            result.results_table["FDR"] > result.fdr_threshold
        ]
        if len(stable_results) > 0:
            stable_file = self.output_dir / f"{base_filename}_stable.csv"
            stable_results.to_csv(stable_file, index=False)
            output_files["stable"] = stable_file

        logger.info(f"Results saved: {len(output_files)} files for {base_filename}")

        return output_files

    def _create_analysis_summary(self, results: Dict[str, DifferentialResult]) -> None:
        """Create summary report of all analyses"""

        summary_data = []

        for result_key, result in results.items():
            summary_data.append(
                {
                    "analysis": result_key,
                    "comparison": result.comparison_name,
                    "method": result.method,
                    "success": result.success,
                    "n_tested": result.n_tested,
                    "n_significant": result.n_significant,
                    "n_up_regulated": result.n_up_regulated,
                    "n_down_regulated": result.n_down_regulated,
                    "percent_significant": (
                        round(100 * result.n_significant / result.n_tested, 2)
                        if result.n_tested and result.n_tested > 0
                        else 0
                    ),
                    "min_fdr": result.min_fdr,
                    "max_abs_logfc": result.max_abs_logfc,
                    "execution_time_min": (
                        round(result.execution_time / 60, 2)
                        if result.execution_time
                        else None
                    ),
                    "error": result.error_message,
                }
            )

        summary_df = pd.DataFrame(summary_data)

        # Save summary
        summary_file = self.output_dir / "differential_analysis_summary.csv"
        summary_df.to_csv(summary_file, index=False)

        logger.info(f"Analysis summary saved to {summary_file}")

        # Print summary to log
        logger.info("=== DIFFERENTIAL ANALYSIS SUMMARY ===")
        for _, row in summary_df.iterrows():
            if row["success"]:
                logger.info(
                    f"{row['analysis']}: {row['n_significant']} significant "
                    f"({row['percent_significant']}%) in {row['execution_time_min']:.1f}min"
                )
            else:
                logger.error(f"{row['analysis']}: FAILED - {row['error']}")

    def compare_methods(
        self, results: Dict[str, DifferentialResult], comparison_name: str
    ) -> pd.DataFrame:
        """Compare results across different methods for the same comparison"""

        # Filter results for this comparison
        comp_results = {
            key: result
            for key, result in results.items()
            if result.comparison_name == comparison_name and result.success
        }

        if len(comp_results) < 2:
            logger.warning(
                f"Need at least 2 successful methods to compare for {comparison_name}"
            )
            return pd.DataFrame()

        # Create comparison DataFrame
        comparison_data = []

        for result_key, result in comp_results.items():
            method = result.method
            comparison_data.append(
                {
                    "method": method,
                    "n_tested": result.n_tested,
                    "n_significant": result.n_significant,
                    "n_up_regulated": result.n_up_regulated,
                    "n_down_regulated": result.n_down_regulated,
                    "percent_significant": round(
                        100 * result.n_significant / result.n_tested, 2
                    ),
                    "min_fdr": result.min_fdr,
                    "max_abs_logfc": result.max_abs_logfc,
                }
            )

        comparison_df = pd.DataFrame(comparison_data)

        # Save comparison
        comp_file = self.output_dir / f"{comparison_name}_method_comparison.csv"
        comparison_df.to_csv(comp_file, index=False)

        logger.info(f"Method comparison saved to {comp_file}")

        return comparison_df

    def create_consensus_results(
        self,
        results: Dict[str, DifferentialResult],
        comparison_name: str,
        min_methods: int = 2,
    ) -> Optional[pd.DataFrame]:
        """Create consensus results from multiple methods"""

        # Get successful results for this comparison
        comp_results = {
            key: result
            for key, result in results.items()
            if result.comparison_name == comparison_name and result.success
        }

        if len(comp_results) < min_methods:
            logger.warning(
                f"Need at least {min_methods} methods for consensus, "
                f"got {len(comp_results)} for {comparison_name}"
            )
            return None

        logger.info(f"Creating consensus results for {comparison_name}")

        # This would be a complex implementation to merge results across methods
        # For now, return the most conservative (highest FDR) significant results

        all_significant = []

        for result in comp_results.values():
            if result.results_table is not None:
                sig_results = result.results_table[
                    result.results_table["regulation"] != "Not Significant"
                ].copy()
                sig_results["method"] = result.method
                all_significant.append(sig_results)

        if not all_significant:
            return pd.DataFrame()

        # Simple consensus: interactions significant in multiple methods
        combined = pd.concat(all_significant, ignore_index=True)

        # Group by interaction ID and count methods
        if "loop_id" in combined.columns:
            consensus = (
                combined.groupby("loop_id")
                .agg(
                    {
                        "method": "count",
                        "FDR": "max",  # Most conservative FDR
                        "logFC": "mean",  # Average fold change
                        "regulation": lambda x: x.mode().iloc[
                            0
                        ],  # Most common regulation
                    }
                )
                .reset_index()
            )

            # Keep only interactions found in at least min_methods
            consensus = consensus[consensus["method"] >= min_methods]
            consensus = consensus.rename(columns={"method": "n_methods"})

            # Save consensus
            consensus_file = self.output_dir / f"{comparison_name}_consensus.csv"
            consensus.to_csv(consensus_file, index=False)

            logger.info(
                f"Consensus results: {len(consensus)} interactions "
                f"found in >= {min_methods} methods"
            )

            return consensus

        return pd.DataFrame()
