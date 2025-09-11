"""
Main loop calling functionality using Peakachu
"""

import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

from ..config import Config
from ..utils import setup_logging, validate_file_exists
from .models import ModelManager, get_recommended_model

logger = logging.getLogger(__name__)


@dataclass
class LoopCallingResult:
    """Result of loop calling for a single sample"""

    sample_name: str
    success: bool
    scores_file: Optional[Path] = None
    loops_file: Optional[Path] = None
    num_loops: Optional[int] = None
    error_message: Optional[str] = None
    execution_time: Optional[float] = None


class PeakachuRunner:
    """Runner for Peakachu commands with error handling and logging"""

    def __init__(self, config: Config):
        self.config = config
        self.model_manager = ModelManager(config.loop_calling.get("model_dir"))

    def run_score_genome(
        self,
        input_cool: Union[str, Path],
        output_file: Union[str, Path],
        model_path: Union[str, Path],
        resolution: int = 10000,
        weight_name: str = "weight",
    ) -> bool:
        """
        Run Peakachu score_genome command

        Args:
            input_cool: Path to input cool file
            output_file: Path to output BEDPE file
            model_path: Path to Peakachu model file
            resolution: Hi-C resolution
            weight_name: Weight column name in cool file

        Returns:
            bool: True if successful, False otherwise
        """
        cmd = [
            "peakachu",
            "score_genome",
            "-r",
            str(resolution),
            "--clr-weight-name",
            weight_name,
            "-p",
            str(input_cool),
            "-O",
            str(output_file),
            "-m",
            str(model_path),
        ]

        return self._execute_command(cmd, "Peakachu score_genome", str(output_file))

    def run_pool_peaks(
        self,
        scores_file: Union[str, Path],
        output_file: Union[str, Path],
        resolution: int = 10000,
        threshold: float = 0.95,
    ) -> bool:
        """
        Run Peakachu pool command to identify final loops

        Args:
            scores_file: Path to scores BEDPE file
            output_file: Path to output loops BEDPE file
            resolution: Hi-C resolution
            threshold: Score threshold for loop calling

        Returns:
            bool: True if successful, False otherwise
        """
        cmd = [
            "peakachu",
            "pool",
            "-r",
            str(resolution),
            "-i",
            str(scores_file),
            "-o",
            str(output_file),
            "-t",
            str(threshold),
        ]

        return self._execute_command(cmd, "Peakachu pool", str(output_file))

    def _execute_command(
        self, cmd: List[str], operation_name: str, expected_output_file: str
    ) -> bool:
        """Execute command with error handling and validation"""
        try:
            logger.info(f"Starting {operation_name}...")
            logger.debug(f"Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=7200,  # 2 hour timeout
            )

            if result.stdout:
                logger.debug(f"{operation_name} stdout: {result.stdout}")

            if result.stderr:
                logger.warning(f"{operation_name} stderr: {result.stderr}")

            # Validate output file was created
            if os.path.exists(expected_output_file):
                file_size = os.path.getsize(expected_output_file)
                logger.info(
                    f"{operation_name} completed successfully. "
                    f"Output file size: {file_size/1024:.2f} KB"
                )
                return True
            else:
                logger.error(
                    f"Expected output file was not created: {expected_output_file}"
                )
                return False

        except subprocess.TimeoutExpired:
            logger.error(f"{operation_name} timed out after 2 hours")
            return False
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running {operation_name}: {e}")
            if e.stderr:
                logger.error(f"Error details: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error in {operation_name}: {e}")
            return False


class LoopCaller:
    """Main class for loop calling analysis"""

    def __init__(self, config: Config):
        self.config = config
        self.peakachu = PeakachuRunner(config)
        self.model_manager = ModelManager(config.loop_calling.get("model_dir"))

        # Create output directory
        self.output_dir = Path(config.output_dir) / "loop_calling"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def call_loops_single_sample(
        self,
        sample_name: str,
        cool_file: Union[str, Path],
        read_count: Optional[int] = None,
        model_override: Optional[str] = None,
    ) -> LoopCallingResult:
        """
        Call loops for a single sample

        Args:
            sample_name: Name of the sample
            cool_file: Path to cool file
            read_count: Number of reads (for model selection)
            model_override: Override automatic model selection

        Returns:
            LoopCallingResult: Results of loop calling
        """
        import time

        start_time = time.time()

        logger.info(f"Processing sample: {sample_name}")

        # Validate input file
        cool_path = Path(cool_file)
        if not cool_path.exists():
            error_msg = f"Cool file not found: {cool_path}"
            logger.error(error_msg)
            return LoopCallingResult(
                sample_name=sample_name, success=False, error_message=error_msg
            )

        # Get model
        if model_override:
            model_name = model_override
        elif read_count:
            model_name = get_recommended_model(read_count)
        else:
            # Try to get from config
            sample_config = self.config.loop_calling.get("models", {}).get(sample_name)
            if sample_config:
                model_name = sample_config.get("model", "100million")
            else:
                model_name = "100million"  # Default

        model_path = self.model_manager.get_model_path(model_name)
        if not model_path:
            error_msg = f"Model not found: {model_name}"
            logger.error(error_msg)
            return LoopCallingResult(
                sample_name=sample_name, success=False, error_message=error_msg
            )

        logger.info(f"Using model: {model_name} ({model_path})")

        # Define output files
        scores_file = self.output_dir / f"{sample_name}_scores.bedpe"
        loops_file = self.output_dir / f"{sample_name}_loops_10kb_non_redundant.bedpe"

        # Step 1: Score genome
        logger.info("Step 1: Running score_genome")
        success = self.peakachu.run_score_genome(
            input_cool=cool_path,
            output_file=scores_file,
            model_path=model_path,
            resolution=self.config.resolution,
            weight_name=self.config.loop_calling.get("weight_name", "weight"),
        )

        if not success:
            error_msg = "Failed to run score_genome"
            logger.error(error_msg)
            return LoopCallingResult(
                sample_name=sample_name,
                success=False,
                error_message=error_msg,
                execution_time=time.time() - start_time,
            )

        # Step 2: Pool peaks
        logger.info("Step 2: Running pool")
        success = self.peakachu.run_pool_peaks(
            scores_file=scores_file,
            output_file=loops_file,
            resolution=self.config.resolution,
            threshold=self.config.loop_calling.get("threshold", 0.95),
        )

        if not success:
            error_msg = "Failed to run pool"
            logger.error(error_msg)
            return LoopCallingResult(
                sample_name=sample_name,
                success=False,
                scores_file=scores_file if scores_file.exists() else None,
                error_message=error_msg,
                execution_time=time.time() - start_time,
            )

        # Count loops
        num_loops = None
        if loops_file.exists():
            try:
                loops_df = pd.read_csv(loops_file, sep="\t", header=None)
                num_loops = len(loops_df)
                logger.info(f"Successfully called {num_loops} loops for {sample_name}")
            except Exception as e:
                logger.warning(f"Could not count loops: {e}")

        execution_time = time.time() - start_time

        return LoopCallingResult(
            sample_name=sample_name,
            success=True,
            scores_file=scores_file,
            loops_file=loops_file,
            num_loops=num_loops,
            execution_time=execution_time,
        )

    def call_loops_batch(
        self,
        samples: Dict[str, Union[str, Path]],
        read_counts: Optional[Dict[str, int]] = None,
        parallel: bool = False,
    ) -> Dict[str, LoopCallingResult]:
        """
        Call loops for multiple samples

        Args:
            samples: Dictionary of sample_name -> cool_file_path
            read_counts: Optional dictionary of sample_name -> read_count
            parallel: Whether to run samples in parallel (not implemented yet)

        Returns:
            Dictionary of sample_name -> LoopCallingResult
        """
        results = {}

        for sample_name, cool_file in samples.items():
            read_count = read_counts.get(sample_name) if read_counts else None

            result = self.call_loops_single_sample(
                sample_name=sample_name, cool_file=cool_file, read_count=read_count
            )

            results[sample_name] = result

        # Log summary
        successful = sum(1 for r in results.values() if r.success)
        total = len(results)
        logger.info(f"Loop calling completed: {successful}/{total} samples successful")

        return results

    def summarize_results(self, results: Dict[str, LoopCallingResult]) -> pd.DataFrame:
        """Create summary DataFrame of results"""

        summary_data = []
        for sample_name, result in results.items():
            summary_data.append(
                {
                    "sample": sample_name,
                    "success": result.success,
                    "num_loops": result.num_loops,
                    "execution_time_min": (
                        result.execution_time / 60 if result.execution_time else None
                    ),
                    "scores_file": (
                        str(result.scores_file) if result.scores_file else None
                    ),
                    "loops_file": str(result.loops_file) if result.loops_file else None,
                    "error": result.error_message,
                }
            )

        return pd.DataFrame(summary_data)
