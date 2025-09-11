"""
Validation utilities for PlasmaFlow
"""

import importlib
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


def validate_file_exists(file_path: Union[str, Path], file_type: str = "file") -> bool:
    """
    Validate that a file exists

    Args:
        file_path: Path to file
        file_type: Type description for error messages

    Returns:
        True if file exists, False otherwise
    """
    path = Path(file_path)

    if not path.exists():
        logger.error(f"{file_type} not found: {path}")
        return False

    if not path.is_file():
        logger.error(f"{file_type} is not a file: {path}")
        return False

    return True


def validate_directory_exists(
    dir_path: Union[str, Path], create_if_missing: bool = False
) -> bool:
    """
    Validate that a directory exists

    Args:
        dir_path: Path to directory
        create_if_missing: Whether to create directory if missing

    Returns:
        True if directory exists or was created, False otherwise
    """
    path = Path(dir_path)

    if not path.exists():
        if create_if_missing:
            try:
                path.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created directory: {path}")
                return True
            except Exception as e:
                logger.error(f"Could not create directory {path}: {e}")
                return False
        else:
            logger.error(f"Directory not found: {path}")
            return False

    if not path.is_dir():
        logger.error(f"Path is not a directory: {path}")
        return False

    return True


def validate_python_packages(packages: List[str]) -> Dict[str, bool]:
    """
    Check if Python packages are available

    Args:
        packages: List of package names

    Returns:
        Dictionary mapping package names to availability status
    """
    results = {}

    for package in packages:
        try:
            importlib.import_module(package)
            results[package] = True
            logger.debug(f"Package {package}: available")
        except ImportError:
            results[package] = False
            logger.debug(f"Package {package}: not available")

    return results


def validate_external_tools(tools: List[str]) -> Dict[str, bool]:
    """
    Check if external command-line tools are available

    Args:
        tools: List of tool names/commands

    Returns:
        Dictionary mapping tool names to availability status
    """
    results = {}

    for tool in tools:
        try:
            result = subprocess.run(
                [tool, "--version"], capture_output=True, text=True, timeout=10
            )
            results[tool] = result.returncode == 0
            if results[tool]:
                logger.debug(f"Tool {tool}: available")
            else:
                logger.debug(f"Tool {tool}: not working properly")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            results[tool] = False
            logger.debug(f"Tool {tool}: not found")

    return results


def validate_r_environment() -> Dict[str, Any]:
    """
    Validate R environment and packages

    Returns:
        Dictionary with R validation results
    """
    results = {"r_available": False, "r_version": None, "packages": {}}

    # Check R availability
    try:
        result = subprocess.run(
            ["R", "--version"], capture_output=True, text=True, timeout=10
        )

        if result.returncode == 0:
            results["r_available"] = True
            # Extract R version
            version_lines = result.stdout.split("\n")
            for line in version_lines:
                if "R version" in line:
                    results["r_version"] = line.strip()
                    break

    except (subprocess.TimeoutExpired, FileNotFoundError):
        logger.debug("R not found")
        return results

    # Check R packages if R is available
    if results["r_available"]:
        try:
            from .r_utils import check_r_packages

            required_packages = [
                "diffHic",
                "InteractionSet",
                "GenomicRanges",
                "DESeq2",
                "edgeR",
                "clusterProfiler",
                "msigdbr",
                "ggplot2",
                "dplyr",
            ]

            package_status = check_r_packages(required_packages)
            results["packages"] = package_status

        except Exception as e:
            logger.warning(f"Could not check R packages: {e}")

    return results


def validate_environment() -> List[str]:
    """
    Comprehensive environment validation

    Returns:
        List of validation issues found
    """
    issues = []

    logger.info("Validating PlasmaFlow environment...")

    # Check Python version
    if sys.version_info < (3, 8):
        issues.append(
            f"Python 3.8+ required, found {sys.version_info.major}.{sys.version_info.minor}"
        )

    # Check core Python packages
    core_packages = [
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "cooler",
        "cooltools",
        "coolpuppy",
        "bioframe",
    ]

    package_status = validate_python_packages(core_packages)
    missing_packages = [
        pkg for pkg, available in package_status.items() if not available
    ]

    if missing_packages:
        issues.append(f"Missing Python packages: {', '.join(missing_packages)}")

    # Check external tools
    external_tools = ["peakachu", "computeMatrix"]
    tool_status = validate_external_tools(external_tools)
    missing_tools = [tool for tool, available in tool_status.items() if not available]

    if missing_tools:
        issues.append(f"Missing external tools: {', '.join(missing_tools)}")

    # Check R environment
    r_results = validate_r_environment()

    if not r_results["r_available"]:
        issues.append("R not available")
    else:
        missing_r_packages = [
            pkg for pkg, available in r_results["packages"].items() if not available
        ]
        if missing_r_packages:
            issues.append(f"Missing R packages: {', '.join(missing_r_packages)}")

    if issues:
        logger.warning(f"Environment validation found {len(issues)} issues")
        for issue in issues:
            logger.warning(f"  - {issue}")
    else:
        logger.info("Environment validation passed")

    return issues


def validate_input_files(config: Any, check_samples: bool = True) -> List[str]:
    """
    Validate input files specified in configuration

    Args:
        config: PlasmaFlow configuration object
        check_samples: Whether to check sample-specific files

    Returns:
        List of validation issues
    """
    issues = []

    # Check input directory
    if config.input_dir and not validate_directory_exists(config.input_dir):
        issues.append(f"Input directory not found: {config.input_dir}")

    # Check reference files if specified
    if hasattr(config, "genomics") and config.genomics:
        ref_files = config.genomics.get("reference_files", {})

        for file_type, file_path in ref_files.items():
            if file_path and not validate_file_exists(file_path, file_type):
                issues.append(f"{file_type} file not found: {file_path}")

    # Check sample files if requested
    if check_samples and hasattr(config, "samples") and config.samples:
        for sample_name, sample_info in config.samples.items():
            if isinstance(sample_info, dict):
                for file_type in ["cool_file", "loops_file", "expression_file"]:
                    file_path = sample_info.get(file_type)
                    if file_path and not validate_file_exists(
                        file_path, f"{sample_name} {file_type}"
                    ):
                        issues.append(
                            f"{sample_name} {file_type} not found: {file_path}"
                        )

    return issues


def validate_output_permissions(output_dir: Union[str, Path]) -> bool:
    """
    Check if output directory is writable

    Args:
        output_dir: Output directory path

    Returns:
        True if writable, False otherwise
    """
    try:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Try to create a test file
        test_file = output_path / ".write_test"
        test_file.write_text("test")
        test_file.unlink()

        return True

    except Exception as e:
        logger.error(f"Output directory not writable: {e}")
        return False


def get_system_info() -> Dict[str, Any]:
    """
    Get system information for debugging

    Returns:
        Dictionary with system information
    """
    import platform

    info = {
        "platform": platform.platform(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "working_directory": os.getcwd(),
        "home_directory": str(Path.home()),
    }

    # Add environment variables relevant to PlasmaFlow
    relevant_env_vars = [
        "PATH",
        "PYTHONPATH",
        "R_HOME",
        "R_LIBS",
        "R_LIBS_USER",
        "CONDA_DEFAULT_ENV",
        "VIRTUAL_ENV",
    ]

    env_vars = {}
    for var in relevant_env_vars:
        if var in os.environ:
            env_vars[var] = os.environ[var]

    info["environment_variables"] = env_vars

    return info
