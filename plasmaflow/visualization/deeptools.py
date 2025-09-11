"""
DeepTools integration for PlasmaFlow
"""

import logging
import os
import shlex
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from ..config import Config
from ..utils import get_logger

logger = get_logger(__name__)


class DeepToolsInterface:
    """Interface for deepTools computeMatrix and related operations"""

    def __init__(self, config: Config):
        """
        Initialize deepTools interface

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config
        self.visualization_params = config.visualization.get("heatmaps", {})

        # Default parameters
        self.default_params = {
            "flanking_region": self.visualization_params.get("flanking_region", 3000),
            "bin_size": self.visualization_params.get("bin_size", 100),
            "reference_point": "center",
            "skip_zeros": True,
            "sort_regions": "keep",
            "n_processors": config.n_threads or 4,
        }

    def validate_bigwig_files(
        self, bigwig_files: Dict[str, Union[str, Path]]
    ) -> Dict[str, Path]:
        """
        Validate BigWig files exist and are readable

        Args:
            bigwig_files: Dictionary mapping sample names to BigWig file paths

        Returns:
            Dictionary of validated BigWig files
        """
        validated_files = {}

        for sample_name, bigwig_path in bigwig_files.items():
            bigwig_file = Path(bigwig_path)

            if not bigwig_file.exists():
                logger.warning(
                    f"BigWig file not found for {sample_name}: {bigwig_file}"
                )
                continue

            if not bigwig_file.is_file():
                logger.warning(
                    f"BigWig path is not a file for {sample_name}: {bigwig_file}"
                )
                continue

            # Check file extension
            if not bigwig_file.suffix.lower() in [".bw", ".bigwig"]:
                logger.warning(
                    f"File does not have BigWig extension for {sample_name}: {bigwig_file}"
                )

            validated_files[sample_name] = bigwig_file
            logger.debug(f"Validated BigWig file for {sample_name}: {bigwig_file}")

        logger.info(
            f"Validated {len(validated_files)}/{len(bigwig_files)} BigWig files"
        )
        return validated_files

    def run_compute_matrix(
        self,
        bed_files: Dict[str, Path],
        bigwig_files: Dict[str, Path],
        output_prefix: str,
        output_dir: Union[str, Path],
        **kwargs,
    ) -> Dict[str, Any]:
        """
        Run computeMatrix with the given BED and BigWig files

        Args:
            bed_files: Dictionary mapping region names to BED file paths
            bigwig_files: Dictionary mapping sample names to BigWig file paths
            output_prefix: Prefix for output files
            output_dir: Output directory
            **kwargs: Additional parameters for computeMatrix

        Returns:
            Dictionary with execution results
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Get parameters (kwargs override defaults)
        params = {**self.default_params, **kwargs}

        # Define output files
        matrix_file = output_path / f"{output_prefix}.gz"
        regions_file = output_path / f"{output_prefix}_regions.bed"

        # Prepare region files and labels (maintain order)
        regions_ordered = []
        region_labels = []

        # Define preferred order for region categories
        preferred_order = [
            "up_anchor1",
            "up_anchor2",
            "down_anchor1",
            "down_anchor2",
            "common_anchor1",
            "common_anchor2",
            "stable_anchor1",
            "stable_anchor2",
        ]

        # Add regions in preferred order
        for region_key in preferred_order:
            if region_key in bed_files:
                regions_ordered.append(str(bed_files[region_key]))
                # Create readable label
                label = region_key.replace("_", " ").title()
                region_labels.append(label)

        # Add any remaining regions
        for region_key, bed_file in bed_files.items():
            if region_key not in preferred_order:
                regions_ordered.append(str(bed_file))
                label = region_key.replace("_", " ").title()
                region_labels.append(label)

        # Prepare BigWig files list (maintain sample order)
        samples_ordered = []
        sample_labels = list(bigwig_files.keys())

        for sample_name in sorted(sample_labels):  # Sort for consistency
            samples_ordered.append(str(bigwig_files[sample_name]))

        logger.info(f"Running computeMatrix with:")
        logger.info(f"  - {len(regions_ordered)} region files")
        logger.info(f"  - {len(samples_ordered)} BigWig files")
        logger.info(f"  - Output matrix: {matrix_file}")

        # Build computeMatrix command
        cmd = (
            [
                "computeMatrix",
                "reference-point",
                "--referencePoint",
                params.get("reference_point", "center"),
                "-b",
                str(params.get("flanking_region", 3000)),
                "-a",
                str(params.get("flanking_region", 3000)),
                "-R",
            ]
            + regions_ordered
            + ["-S"]
            + samples_ordered
            + [
                "--outFileName",
                str(matrix_file),
                "--outFileSortedRegions",
                str(regions_file),
                "--numberOfProcessors",
                str(params.get("n_processors", 4)),
            ]
        )

        # Add optional parameters
        if params.get("skip_zeros", True):
            cmd.append("--skipZeros")

        if params.get("sort_regions"):
            cmd.extend(["--sortRegions", params["sort_regions"]])

        if params.get("bin_size"):
            cmd.extend(["--binSize", str(params["bin_size"])])

        # Add verbosity
        cmd.append("--verbose")

        logger.debug(f"computeMatrix command: {shlex.join(cmd)}")

        # Execute command
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
            )

            logger.info("computeMatrix completed successfully")

            if result.stdout:
                logger.debug(f"computeMatrix stdout: {result.stdout}")

            if result.stderr:
                logger.debug(f"computeMatrix stderr: {result.stderr}")

            # Verify output files exist
            if not matrix_file.exists():
                raise RuntimeError(f"Matrix file was not created: {matrix_file}")

            return {
                "success": True,
                "matrix_file": matrix_file,
                "regions_file": regions_file if regions_file.exists() else None,
                "region_labels": region_labels,
                "sample_labels": sample_labels,
                "command": shlex.join(cmd),
                "stdout": result.stdout,
                "stderr": result.stderr,
            }

        except subprocess.CalledProcessError as e:
            error_msg = f"computeMatrix failed (exit code {e.returncode}): {e.stderr}"
            logger.error(error_msg)

            return {
                "success": False,
                "error": error_msg,
                "command": shlex.join(cmd),
                "stdout": e.stdout,
                "stderr": e.stderr,
            }

        except subprocess.TimeoutExpired:
            error_msg = "computeMatrix timed out after 1 hour"
            logger.error(error_msg)

            return {"success": False, "error": error_msg, "command": shlex.join(cmd)}

        except FileNotFoundError:
            error_msg = (
                "computeMatrix command not found. Make sure deepTools is installed."
            )
            logger.error(error_msg)

            return {"success": False, "error": error_msg, "command": shlex.join(cmd)}

    def run_plot_heatmap(
        self, matrix_file: Union[str, Path], output_file: Union[str, Path], **kwargs
    ) -> Dict[str, Any]:
        """
        Generate heatmap from matrix using plotHeatmap

        Args:
            matrix_file: Path to computeMatrix output file
            output_file: Path for output heatmap image
            **kwargs: Additional parameters for plotHeatmap

        Returns:
            Dictionary with execution results
        """
        matrix_path = Path(matrix_file)
        output_path = Path(output_file)

        if not matrix_path.exists():
            raise FileNotFoundError(f"Matrix file not found: {matrix_path}")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Default heatmap parameters
        heatmap_params = {
            "colorMap": self.visualization_params.get("colormap", "RdBu_r"),
            "zMin": self.visualization_params.get("vmin", -3),
            "zMax": self.visualization_params.get("vmax", 3),
            "dpi": 300,
        }

        # Update with kwargs
        heatmap_params.update(kwargs)

        # Build plotHeatmap command
        cmd = [
            "plotHeatmap",
            "-m",
            str(matrix_path),
            "-o",
            str(output_path),
            "--colorMap",
            heatmap_params["colorMap"],
            "--zMin",
            str(heatmap_params["zMin"]),
            "--zMax",
            str(heatmap_params["zMax"]),
            "--dpi",
            str(heatmap_params["dpi"]),
        ]

        # Add optional parameters
        if "plotTitle" in heatmap_params:
            cmd.extend(["--plotTitle", heatmap_params["plotTitle"]])

        logger.info(f"Running plotHeatmap: {output_path}")
        logger.debug(f"plotHeatmap command: {shlex.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout
            )

            logger.info("plotHeatmap completed successfully")

            return {
                "success": True,
                "output_file": output_path,
                "command": shlex.join(cmd),
                "stdout": result.stdout,
                "stderr": result.stderr,
            }

        except subprocess.CalledProcessError as e:
            error_msg = f"plotHeatmap failed (exit code {e.returncode}): {e.stderr}"
            logger.error(error_msg)

            return {
                "success": False,
                "error": error_msg,
                "command": shlex.join(cmd),
                "stderr": e.stderr,
            }

        except subprocess.TimeoutExpired:
            error_msg = "plotHeatmap timed out after 10 minutes"
            logger.error(error_msg)

            return {"success": False, "error": error_msg, "command": shlex.join(cmd)}

        except FileNotFoundError:
            error_msg = (
                "plotHeatmap command not found. Make sure deepTools is installed."
            )
            logger.error(error_msg)

            return {"success": False, "error": error_msg, "command": shlex.join(cmd)}


def validate_bigwig_files(bigwig_files: Dict[str, Union[str, Path]]) -> Dict[str, bool]:
    """
    Validate BigWig files (standalone function)

    Args:
        bigwig_files: Dictionary of sample_name -> file_path

    Returns:
        Dictionary of sample_name -> is_valid
    """
    validation_results = {}

    for sample_name, file_path in bigwig_files.items():
        file_path_obj = Path(file_path)

        is_valid = (
            file_path_obj.exists()
            and file_path_obj.is_file()
            and file_path_obj.suffix.lower() in [".bw", ".bigwig"]
        )

        validation_results[sample_name] = is_valid

        if not is_valid:
            logger.warning(f"Invalid BigWig file for {sample_name}: {file_path}")

    return validation_results
