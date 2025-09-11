"""
Aggregate Peak Analysis (APA) functionality
"""

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cooler
import cooltools
import numpy as np
import pandas as pd
from coolpuppy import coolpup

from ..config import Config
from ..utils import setup_logging, validate_file_exists
from .caching import load_cached_apa, save_apa_cache

logger = logging.getLogger(__name__)


@dataclass
class APAResult:
    """Result of APA analysis for a single sample"""

    sample_name: str
    success: bool

    # APA matrix and metadata
    apa_matrix: Optional[np.ndarray] = None
    pileup_data: Optional[pd.DataFrame] = None

    # APA scores
    center_score: Optional[float] = None
    corner_scores: Optional[Dict[str, float]] = None

    # Analysis parameters
    flank: Optional[int] = None
    min_diag: Optional[int] = None
    maxdist: Optional[int] = None

    # File paths
    cache_file: Optional[Path] = None

    # Execution info
    execution_time: Optional[float] = None
    error_message: Optional[str] = None


class APAAnalyzer:
    """Main class for Aggregate Peak Analysis"""

    def __init__(self, config: Config):
        self.config = config
        self.apa_params = config.quality_control.get("apa", {})

        # Create output directory
        self.output_dir = Path(config.output_dir) / "quality_control" / "apa"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Cache directory
        self.cache_dir = self.output_dir / "cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Default parameters
        self.default_params = {
            "flank": 200000,
            "min_diag": 3,
            "maxdist": 1000000,
            "nproc": 4,
        }

    def generate_expected_values(
        self, clr: cooler.Cooler, nproc: int = 4
    ) -> pd.DataFrame:
        """Generate expected interaction values for a cooler file"""
        logger.info("Generating expected values...")
        return cooltools.expected_cis(clr=clr, nproc=nproc)

    def load_loops(self, bedpe_file: Union[str, Path]) -> pd.DataFrame:
        """Load loops from a BEDPE file with standardized column names"""
        loops_df = pd.read_csv(
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

        # Handle cases with fewer columns
        actual_cols = len(loops_df.columns)
        expected_cols = [
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "score",
            "freq",
        ]
        loops_df.columns = expected_cols[:actual_cols]

        logger.info(f"Loaded {len(loops_df)} loops from {bedpe_file}")
        return loops_df

    def perform_pileup(
        self,
        clr: cooler.Cooler,
        loops: pd.DataFrame,
        expected_df: pd.DataFrame,
        flank: int = 200000,
        min_diag: int = 3,
        maxdist: Optional[int] = None,
        nproc: int = 4,
    ) -> pd.DataFrame:
        """Perform pile-up analysis on loops"""
        logger.info(
            f"Performing pileup with flank={flank}, min_diag={min_diag}, maxdist={maxdist}"
        )

        return coolpup.pileup(
            clr,
            loops,
            features_format="bedpe",
            expected_df=expected_df,
            flank=flank,
            min_diag=min_diag,
            maxdist=maxdist,
            nproc=nproc,
        )

    def analyze_single_sample(
        self,
        sample_name: str,
        cool_file: Union[str, Path],
        bedpe_file: Union[str, Path],
        flank: Optional[int] = None,
        min_diag: Optional[int] = None,
        maxdist: Optional[int] = None,
        nproc: Optional[int] = None,
        force_recompute: bool = False,
    ) -> APAResult:
        """
        Perform APA analysis on a single sample

        Args:
            sample_name: Name of the sample
            cool_file: Path to Hi-C cool file
            bedpe_file: Path to loops BEDPE file
            flank: Flanking region size (bp)
            min_diag: Minimum diagonal distance
            maxdist: Maximum interaction distance
            nproc: Number of processors
            force_recompute: Force recomputation even if cached result exists

        Returns:
            APAResult object
        """
        start_time = time.time()

        logger.info(f"Starting APA analysis for sample: {sample_name}")

        # Use provided parameters or defaults
        params = {
            "flank": flank
            or self.apa_params.get("flank", self.default_params["flank"]),
            "min_diag": min_diag
            or self.apa_params.get("min_diag", self.default_params["min_diag"]),
            "maxdist": maxdist
            or self.apa_params.get("maxdist", self.default_params["maxdist"]),
            "nproc": nproc
            or self.apa_params.get("nproc", self.default_params["nproc"]),
        }

        logger.info(f"APA parameters: {params}")

        # Check for cached result
        cache_file = self.cache_dir / f"{sample_name}_apa_cache.pkl"

        if not force_recompute and cache_file.exists():
            try:
                cached_result = load_cached_apa(cache_file, params)
                if cached_result:
                    logger.info(f"Using cached APA result for {sample_name}")
                    return cached_result
            except Exception as e:
                logger.warning(f"Could not load cached result: {e}")

        # Validate input files
        cool_path = Path(cool_file)
        bedpe_path = Path(bedpe_file)

        if not cool_path.exists():
            error_msg = f"Cool file not found: {cool_path}"
            logger.error(error_msg)
            return APAResult(
                sample_name=sample_name, success=False, error_message=error_msg
            )

        if not bedpe_path.exists():
            error_msg = f"BEDPE file not found: {bedpe_path}"
            logger.error(error_msg)
            return APAResult(
                sample_name=sample_name, success=False, error_message=error_msg
            )

        try:
            # Load Hi-C contact map
            logger.info("Loading Hi-C contact map...")
            clr = cooler.Cooler(str(cool_path))

            # Generate expected values
            cis_exp = self.generate_expected_values(clr, params["nproc"])

            # Load loops
            loops = self.load_loops(bedpe_path)

            if len(loops) == 0:
                error_msg = "No loops found in BEDPE file"
                logger.error(error_msg)
                return APAResult(
                    sample_name=sample_name, success=False, error_message=error_msg
                )

            # Perform pileup
            pileup = self.perform_pileup(
                clr=clr, loops=loops, expected_df=cis_exp, **params
            )

            # Extract APA matrix
            apa_matrix = self._extract_apa_matrix(pileup)

            # Calculate APA scores
            center_score, corner_scores = calculate_apa_scores(apa_matrix)

            execution_time = time.time() - start_time

            # Create result
            result = APAResult(
                sample_name=sample_name,
                success=True,
                apa_matrix=apa_matrix,
                pileup_data=pileup,
                center_score=center_score,
                corner_scores=corner_scores,
                cache_file=cache_file,
                execution_time=execution_time,
                **params,
            )

            # Cache the result
            try:
                save_apa_cache(result, cache_file, params)
            except Exception as e:
                logger.warning(f"Could not cache result: {e}")

            logger.info(
                f"APA analysis completed for {sample_name}. "
                f"Center score: {center_score:.3f}, "
                f"Execution time: {execution_time:.1f}s"
            )

            return result

        except Exception as e:
            error_msg = f"APA analysis failed: {e}"
            logger.error(error_msg, exc_info=True)
            return APAResult(
                sample_name=sample_name,
                success=False,
                error_message=error_msg,
                execution_time=time.time() - start_time,
            )

    def _extract_apa_matrix(self, pileup: pd.DataFrame) -> np.ndarray:
        """Extract APA matrix from pileup DataFrame"""
        try:
            # Get the matrix from the "data" column
            mat_obj = pileup["data"].iloc[0]

            if isinstance(mat_obj, (list, np.ndarray)):
                matrix = np.asarray(mat_obj, dtype=float)
            else:
                raise ValueError("Data column does not contain a list/array")

            if matrix.size == 0:
                raise ValueError("Empty matrix")

            return matrix

        except Exception as e:
            logger.error(f"Failed to extract APA matrix: {e}")
            raise

    def analyze_batch(
        self, samples: Dict[str, Dict[str, Union[str, Path]]], **kwargs
    ) -> Dict[str, APAResult]:
        """
        Analyze multiple samples

        Args:
            samples: Dict of sample_name -> {"cool_file": path, "bedpe_file": path}
            **kwargs: Additional parameters passed to analyze_single_sample

        Returns:
            Dictionary of sample_name -> APAResult
        """
        results = {}

        logger.info(f"Starting batch APA analysis for {len(samples)} samples")

        for sample_name, files in samples.items():
            result = self.analyze_single_sample(
                sample_name=sample_name,
                cool_file=files["cool_file"],
                bedpe_file=files["bedpe_file"],
                **kwargs,
            )
            results[sample_name] = result

        # Log summary
        successful = sum(1 for r in results.values() if r.success)
        logger.info(
            f"Batch APA analysis completed: {successful}/{len(results)} samples successful"
        )

        return results

    def create_summary_report(self, results: Dict[str, APAResult]) -> pd.DataFrame:
        """Create summary report of APA results"""

        summary_data = []

        for sample_name, result in results.items():
            summary_data.append(
                {
                    "sample": sample_name,
                    "success": result.success,
                    "center_score": result.center_score,
                    "top_left_score": (
                        result.corner_scores.get("top_left")
                        if result.corner_scores
                        else None
                    ),
                    "top_right_score": (
                        result.corner_scores.get("top_right")
                        if result.corner_scores
                        else None
                    ),
                    "bottom_left_score": (
                        result.corner_scores.get("bottom_left")
                        if result.corner_scores
                        else None
                    ),
                    "bottom_right_score": (
                        result.corner_scores.get("bottom_right")
                        if result.corner_scores
                        else None
                    ),
                    "flank": result.flank,
                    "min_diag": result.min_diag,
                    "maxdist": result.maxdist,
                    "execution_time_min": (
                        result.execution_time / 60 if result.execution_time else None
                    ),
                    "error": result.error_message,
                }
            )

        summary_df = pd.DataFrame(summary_data)

        # Save summary
        summary_file = self.output_dir / "apa_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"APA summary saved to {summary_file}")

        return summary_df


def calculate_apa_scores(matrix: np.ndarray) -> Tuple[float, Dict[str, float]]:
    """
    Calculate APA scores for center and corners of the matrix

    Args:
        matrix: 2D numpy array of APA values

    Returns:
        Tuple of (center_score, corner_scores_dict)
    """
    if matrix.size == 0:
        return 0.0, {}

    n_rows, n_cols = matrix.shape

    # Center score
    center_row = n_rows // 2
    center_col = n_cols // 2
    center_score = float(matrix[center_row, center_col])

    # Corner scores (use 3x3 regions)
    corner_size = 3

    corner_scores = {
        "top_left": float(np.mean(matrix[:corner_size, :corner_size])),
        "top_right": float(np.mean(matrix[:corner_size, -corner_size:])),
        "bottom_left": float(np.mean(matrix[-corner_size:, :corner_size])),
        "bottom_right": float(np.mean(matrix[-corner_size:, -corner_size:])),
    }

    return center_score, corner_scores
