"""
Path configuration and validation for PlasmaFlow
"""

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)


@dataclass
class PathConfig:
    """Configuration for file paths used in PlasmaFlow analysis"""

    # Base directories
    project_root: Optional[Path] = None
    input_dir: Optional[Path] = None
    output_dir: Optional[Path] = None
    temp_dir: Optional[Path] = None

    # Hi-C data paths
    hic_data_dir: Optional[Path] = None
    cool_files_dir: Optional[Path] = None
    loops_dir: Optional[Path] = None

    # Expression data paths
    expression_data_dir: Optional[Path] = None

    # Reference data paths
    reference_dir: Optional[Path] = None
    gtf_file: Optional[Path] = None
    genome_fasta: Optional[Path] = None

    # External tool paths
    peakachu_models_dir: Optional[Path] = None
    deeptools_bin: Optional[Path] = None

    # Bigwig files
    bigwig_dir: Optional[Path] = None
    histone_marks_dir: Optional[Path] = None

    def __post_init__(self):
        """Convert string paths to Path objects"""
        for field_name, field_value in self.__dict__.items():
            if field_value is not None and isinstance(field_value, str):
                setattr(self, field_name, Path(field_value))

    def create_output_dirs(self) -> None:
        """Create output directories if they don't exist"""
        dirs_to_create = [
            self.output_dir,
            self.temp_dir,
        ]

        for dir_path in dirs_to_create:
            if dir_path is not None:
                try:
                    dir_path.mkdir(parents=True, exist_ok=True)
                    logger.info(f"Created directory: {dir_path}")
                except Exception as e:
                    logger.error(f"Failed to create directory {dir_path}: {e}")
                    raise

    def get_sample_cool_file(self, sample_name: str) -> Optional[Path]:
        """Get cool file path for a specific sample"""
        if self.cool_files_dir is None:
            return None

        # Try different naming patterns
        patterns = [
            f"{sample_name}_merged.10000_balanced.cool",
            f"{sample_name}.10000_balanced.cool",
            f"{sample_name}_10kb.cool",
            f"{sample_name}.cool",
        ]

        for pattern in patterns:
            cool_file = self.cool_files_dir / pattern
            if cool_file.exists():
                return cool_file

        logger.warning(f"Cool file not found for sample {sample_name}")
        return None

    def get_sample_loops_file(self, sample_name: str) -> Optional[Path]:
        """Get loops file path for a specific sample"""
        if self.loops_dir is None:
            return None

        # Try different naming patterns
        patterns = [
            f"{sample_name}_loops_10kb_non_redundant.bedpe",
            f"{sample_name}_loops.bedpe",
            f"{sample_name}.bedpe",
        ]

        for pattern in patterns:
            loops_file = self.loops_dir / pattern
            if loops_file.exists():
                return loops_file

        logger.warning(f"Loops file not found for sample {sample_name}")
        return None

    def get_output_subdir(self, subdir_name: str) -> Path:
        """Get a subdirectory within the output directory"""
        if self.output_dir is None:
            raise ValueError("Output directory not set")

        subdir = self.output_dir / subdir_name
        subdir.mkdir(parents=True, exist_ok=True)
        return subdir

    def validate(self) -> List[str]:
        """Validate path configuration"""
        issues = []

        # Check input directories exist
        input_dirs = [
            ("input_dir", self.input_dir),
            ("hic_data_dir", self.hic_data_dir),
            ("cool_files_dir", self.cool_files_dir),
            ("reference_dir", self.reference_dir),
        ]

        for name, path in input_dirs:
            if path is not None and not path.exists():
                issues.append(f"{name} does not exist: {path}")

        # Check required files exist
        required_files = [
            ("GTF file", self.gtf_file),
            ("Genome FASTA", self.genome_fasta),
        ]

        for name, file_path in required_files:
            if file_path is not None and not file_path.exists():
                issues.append(f"{name} does not exist: {file_path}")

        # Check output directory can be created
        if self.output_dir is not None:
            try:
                self.output_dir.parent.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                issues.append(f"Cannot create output directory parent: {e}")

        return issues


def get_default_paths() -> PathConfig:
    """Get default path configuration"""

    # Try to detect common paths
    home = Path.home()

    # Default paths based on the original notebook structure
    default_paths = PathConfig(
        project_root=Path.cwd(),
        input_dir=Path("/zdata/HMCL/Hi-C_treated_data"),
        output_dir=home / "Documents" / "loops" / "plasmaflow_output",
        temp_dir=home / "Documents" / "loops" / "temp",
        # Hi-C data
        hic_data_dir=Path(
            "/zdata/HMCL/Hi-C_treated_data/hic_data_for_postpross/data_diff/balanced/data_10kb/merged"
        ),
        cool_files_dir=Path("/zdata/HMCL/Hi-C_treated_data"),
        loops_dir=Path("/home/miss-leriem.zellagui/loop_calling/output"),
        # Expression data
        expression_data_dir=home / "Documents" / "loops" / "expression",
        # Reference data
        reference_dir=home / "Documents" / "loops" / "reference",
        # External tools
        peakachu_models_dir=Path("/home/miss-leriem.zellagui/loop_calling/models"),
        # Bigwig files
        bigwig_dir=home / "Documents" / "loops" / "Bigwig",
    )

    return default_paths


def validate_paths(path_config: PathConfig) -> Dict[str, Any]:
    """Validate path configuration and return detailed results"""

    results = {"valid": True, "issues": [], "warnings": [], "created_dirs": []}

    # Validate configuration
    issues = path_config.validate()
    results["issues"] = issues

    if issues:
        results["valid"] = False
        return results

    # Try to create output directories
    try:
        path_config.create_output_dirs()
        if path_config.output_dir:
            results["created_dirs"].append(str(path_config.output_dir))
        if path_config.temp_dir:
            results["created_dirs"].append(str(path_config.temp_dir))
    except Exception as e:
        results["valid"] = False
        results["issues"].append(f"Failed to create output directories: {e}")

    return results
