"""
Core PlasmaFlow analysis orchestrator
"""

import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .config import Config, load_config
from .differential import DifferentialAnalyzer
from .enrichment import EnrichmentAnalyzer  # Will create this
from .genomics import (ExpressionAnalyzer,  # Will create these
                       GeneProximityAnalyzer)
from .loop_calling import LoopCaller
from .quality_control import APAAnalyzer
from .utils import get_logger, setup_logging, validate_environment
from .visualization import create_apa_plots, create_volcano_plot

logger = get_logger(__name__)


class PlasmaFlowAnalysis:
    """
    Main orchestrator class for PlasmaFlow Hi-C chromatin loop analysis pipeline

    This class coordinates the complete workflow from loop calling through
    pathway enrichment analysis.
    """

    def __init__(
        self,
        config: Union[str, Path, Config, Dict[str, Any]],
        log_level: str = "INFO",
        log_file: Optional[Union[str, Path]] = None,
    ):
        """
        Initialize PlasmaFlow analysis

        Args:
            config: Configuration file path, Config object, or config dict
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
            log_file: Optional log file path
        """

        # Set up logging first
        setup_logging(level=log_level, log_file=log_file)
        logger.info("Initializing PlasmaFlow analysis pipeline")

        # Load configuration
        if isinstance(config, (str, Path)):
            self.config = load_config(config)
        elif isinstance(config, dict):
            self.config = Config(**config)
        elif isinstance(config, Config):
            self.config = config
        else:
            raise ValueError(
                "Invalid config type. Expected str, Path, dict, or Config object"
            )

        # Validate environment
        self._validate_environment()

        # Initialize components
        self._initialize_components()

        # Track results
        self.results = {}
        self.execution_times = {}

        logger.info("PlasmaFlow pipeline initialized successfully")

    def _validate_environment(self) -> None:
        """Validate that the environment is properly set up"""

        logger.info("Validating environment...")

        # Validate configuration
        from .config import validate_config

        issues = validate_config(self.config)

        if issues:
            logger.warning("Configuration issues found:")
            for issue in issues:
                logger.warning(f"  - {issue}")

        # Check dependencies
        env_issues = validate_environment()
        if env_issues:
            logger.warning("Environment issues found:")
            for issue in env_issues:
                logger.warning(f"  - {issue}")

        # Create output directories
        if self.config.output_dir:
            Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)

        if self.config.temp_dir:
            Path(self.config.temp_dir).mkdir(parents=True, exist_ok=True)

    def _initialize_components(self) -> None:
        """Initialize analysis components"""

        logger.info("Initializing analysis components...")

        # Core analyzers
        self.loop_caller = LoopCaller(self.config)
        self.apa_analyzer = APAAnalyzer(self.config)
        self.differential_analyzer = DifferentialAnalyzer(self.config)

        # Optional analyzers (will be created later)
        self.gene_proximity_analyzer = None
        self.expression_analyzer = None
        self.enrichment_analyzer = None

        try:
            # These will be implemented in subsequent modules
            # For now, we'll create placeholder implementations
            pass
        except ImportError as e:
            logger.warning(f"Some analyzers not available: {e}")

    def run_full_pipeline(
        self,
        samples: Optional[Dict[str, Dict[str, str]]] = None,
        steps: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Run the complete PlasmaFlow analysis pipeline

        Args:
            samples: Dictionary of sample configurations
            steps: List of pipeline steps to run

        Returns:
            Dictionary containing all analysis results
        """

        logger.info("=" * 60)
        logger.info("Starting PlasmaFlow complete analysis pipeline")
        logger.info("=" * 60)

        start_time = time.time()

        # Default steps
        if steps is None:
            steps = [
                "loop_calling",
                "quality_control",
                "differential_analysis",
                "visualization",
            ]

        # Use configured samples if not provided
        if samples is None:
            samples = self._get_default_samples()

        # Run each step
        for step in steps:
            try:
                step_start = time.time()
                logger.info(f"\n{'='*20} STEP: {step.upper()} {'='*20}")

                if step == "loop_calling":
                    self.results["loop_calling"] = self.run_loop_calling(samples)
                elif step == "quality_control":
                    self.results["quality_control"] = self.run_quality_control()
                elif step == "differential_analysis":
                    self.results["differential_analysis"] = (
                        self.run_differential_analysis()
                    )
                elif step == "visualization":
                    self.results["visualization"] = self.run_visualization()
                elif step == "gene_proximity":
                    self.results["gene_proximity"] = self.run_gene_proximity_analysis()
                elif step == "expression_analysis":
                    self.results["expression_analysis"] = self.run_expression_analysis()
                elif step == "enrichment_analysis":
                    self.results["enrichment_analysis"] = self.run_enrichment_analysis()
                else:
                    logger.warning(f"Unknown pipeline step: {step}")
                    continue

                step_time = time.time() - step_start
                self.execution_times[step] = step_time
                logger.info(f"Step {step} completed in {step_time:.2f} seconds")

            except Exception as e:
                logger.error(f"Step {step} failed: {e}", exc_info=True)
                self.results[step] = {"success": False, "error": str(e)}

        total_time = time.time() - start_time
        self.execution_times["total"] = total_time

        # Create final summary
        self._create_pipeline_summary()

        logger.info("=" * 60)
        logger.info(f"PlasmaFlow pipeline completed in {total_time:.2f} seconds")
        logger.info("=" * 60)

        return self.results

    def run_loop_calling(self, samples: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
        """Run loop calling analysis"""

        logger.info("Running loop calling analysis")

        # Convert samples format for LoopCaller
        cool_files = {}
        read_counts = {}

        for sample_name, sample_info in samples.items():
            cool_files[sample_name] = sample_info.get("cool_file")
            read_counts[sample_name] = sample_info.get("read_count")

        # Run loop calling
        results = self.loop_caller.call_loops_batch(
            samples=cool_files, read_counts=read_counts
        )

        # Create summary
        summary = self.loop_caller.summarize_results(results)

        return {"results": results, "summary": summary, "success": True}

    def run_quality_control(self) -> Dict[str, Any]:
        """Run quality control analysis (APA)"""

        logger.info("Running quality control (APA) analysis")

        # Get loop calling results
        if "loop_calling" not in self.results:
            raise RuntimeError("Loop calling must be run before quality control")

        loop_results = self.results["loop_calling"]["results"]

        # Prepare samples for APA
        apa_samples = {}
        for sample_name, loop_result in loop_results.items():
            if loop_result.success:
                # Need to get corresponding cool file
                # This would be improved with better sample management
                apa_samples[sample_name] = {
                    "cool_file": f"/path/to/{sample_name}.cool",  # Placeholder
                    "bedpe_file": loop_result.loops_file,
                }

        # Run APA analysis
        apa_results = self.apa_analyzer.analyze_batch(apa_samples)

        # Create summary
        summary = self.apa_analyzer.create_summary_report(apa_results)

        return {"results": apa_results, "summary": summary, "success": True}

    def run_differential_analysis(self) -> Dict[str, Any]:
        """Run differential analysis"""

        logger.info("Running differential analysis")

        # This would typically use a counts matrix and metadata
        # For now, we'll use placeholder paths
        counts_matrix = (
            Path(self.config.input_dir)
            / "raw_whole_genome_counts_matrix_with_origins.csv"
        )
        metadata = Path(self.config.input_dir) / "metadata.csv"

        if not counts_matrix.exists():
            logger.warning(f"Counts matrix not found: {counts_matrix}")
            return {"success": False, "error": "Counts matrix not found"}

        if not metadata.exists():
            logger.warning(f"Metadata not found: {metadata}")
            return {"success": False, "error": "Metadata not found"}

        # Run differential analysis
        results = self.differential_analyzer.run_differential_analysis(
            counts_matrix=counts_matrix, metadata=metadata
        )

        return {"results": results, "success": True}

    def run_visualization(self) -> Dict[str, Any]:
        """Run visualization generation"""

        logger.info("Running visualization generation")

        vis_results = {}

        # APA plots
        if "quality_control" in self.results:
            apa_results = self.results["quality_control"]["results"]
            if apa_results:
                apa_plots = create_apa_plots(
                    apa_results,
                    output_dir=Path(self.config.output_dir) / "visualizations" / "apa",
                )
                vis_results["apa_plots"] = apa_plots

        # Volcano plots from differential analysis
        if "differential_analysis" in self.results:
            diff_results = self.results["differential_analysis"]["results"]
            for result_key, diff_result in diff_results.items():
                if diff_result.success and diff_result.results_table is not None:
                    volcano_plot = create_volcano_plot(
                        diff_result.results_table,
                        output_file=Path(self.config.output_dir)
                        / "visualizations"
                        / f"{result_key}_volcano.png",
                    )
                    vis_results[f"{result_key}_volcano"] = volcano_plot

        return {"plots": vis_results, "success": True}

    def run_gene_proximity_analysis(self) -> Dict[str, Any]:
        """Run gene proximity analysis (placeholder)"""

        logger.info("Gene proximity analysis not yet implemented")
        return {"success": False, "error": "Not implemented"}

    def run_expression_analysis(self) -> Dict[str, Any]:
        """Run expression analysis (placeholder)"""

        logger.info("Expression analysis not yet implemented")
        return {"success": False, "error": "Not implemented"}

    def run_enrichment_analysis(self) -> Dict[str, Any]:
        """Run pathway enrichment analysis (placeholder)"""

        logger.info("Enrichment analysis not yet implemented")
        return {"success": False, "error": "Not implemented"}

    def _get_default_samples(self) -> Dict[str, Dict[str, str]]:
        """Get default sample configuration"""

        from .config import get_default_sample_config

        sample_config = get_default_sample_config()

        samples = {}
        for cell_type_name, cell_type_config in sample_config.cell_types.items():
            samples[cell_type_name] = {
                "cool_file": cell_type_config.cool_file,
                "read_count": cell_type_config.read_count,
            }

        return samples

    def _create_pipeline_summary(self) -> None:
        """Create comprehensive pipeline summary"""

        logger.info("\n" + "=" * 50)
        logger.info("PLASMAFLOW PIPELINE SUMMARY")
        logger.info("=" * 50)

        # Execution times
        logger.info("EXECUTION TIMES:")
        for step, exec_time in self.execution_times.items():
            if step != "total":
                logger.info(f"  {step}: {exec_time:.2f} seconds")
        logger.info(f"  TOTAL: {self.execution_times.get('total', 0):.2f} seconds")

        # Results summary
        logger.info("\nRESULTS SUMMARY:")
        for step, result in self.results.items():
            if isinstance(result, dict) and "success" in result:
                status = "✓ SUCCESS" if result["success"] else "✗ FAILED"
                logger.info(f"  {step}: {status}")

                if not result["success"] and "error" in result:
                    logger.info(f"    Error: {result['error']}")

        # Save summary to file
        summary_file = Path(self.config.output_dir) / "pipeline_summary.txt"
        with open(summary_file, "w") as f:
            f.write("PlasmaFlow Pipeline Summary\\n")
            f.write("=" * 30 + "\\n\\n")

            f.write("Configuration:\\n")
            f.write(f"  Project: {self.config.project_name}\\n")
            f.write(f"  Resolution: {self.config.resolution}\\n")
            f.write(f"  Output directory: {self.config.output_dir}\\n\\n")

            f.write("Execution Times:\\n")
            for step, exec_time in self.execution_times.items():
                f.write(f"  {step}: {exec_time:.2f} seconds\\n")

            f.write("\\nResults:\\n")
            for step, result in self.results.items():
                if isinstance(result, dict) and "success" in result:
                    status = "SUCCESS" if result["success"] else "FAILED"
                    f.write(f"  {step}: {status}\\n")

        logger.info(f"\\nPipeline summary saved to: {summary_file}")

    def get_results(self) -> Dict[str, Any]:
        """Get all pipeline results"""
        return self.results

    def get_execution_times(self) -> Dict[str, float]:
        """Get execution times for all steps"""
        return self.execution_times

    def save_results(
        self, output_file: Union[str, Path], format: str = "pickle"
    ) -> None:
        """
        Save results to file

        Args:
            output_file: Path to output file
            format: Format ('pickle', 'json')
        """

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "pickle":
            import pickle

            with open(output_path, "wb") as f:
                pickle.dump(self.results, f)
        elif format == "json":
            import json

            # Convert results to JSON-serializable format
            json_results = self._results_to_json()
            with open(output_path, "w") as f:
                json.dump(json_results, f, indent=2)
        else:
            raise ValueError(f"Unsupported format: {format}")

        logger.info(f"Results saved to {output_path}")

    def _results_to_json(self) -> Dict[str, Any]:
        """Convert results to JSON-serializable format"""

        json_results = {}

        for key, value in self.results.items():
            if isinstance(value, dict):
                json_results[key] = {}
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, pd.DataFrame):
                        json_results[key][subkey] = subvalue.to_dict("records")
                    elif hasattr(subvalue, "__dict__"):
                        json_results[key][subkey] = str(subvalue)
                    else:
                        json_results[key][subkey] = subvalue
            else:
                json_results[key] = str(value)

        return json_results
