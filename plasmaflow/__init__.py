"""
PlasmaFlow: A comprehensive Hi-C chromatin loop analysis pipeline

PlasmaFlow is a Python package for analyzing Hi-C chromatin interaction data
through the lens of plasma cell differentiation. It provides a complete pipeline
from loop calling to advanced pathway enrichment analysis with MSigDB integration.

Main Components:
- Loop calling with Peakachu
- Quality control and APA analysis
- Differential loop analysis (diffHic, DESeq2, edgeR)
- Matrix generation and visualization
- Gene proximity mapping
- Differential expression analysis
- GSEA and over-representation analysis

Example:
    >>> from plasmaflow import PlasmaFlowAnalysis
    >>> analysis = PlasmaFlowAnalysis(config_file="config.yaml")
    >>> results = analysis.run_full_pipeline()
"""

import logging
import sys
from importlib import metadata
from typing import Any, Dict

try:
    __version__ = metadata.version("plasmaflow")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0-dev"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

# Module imports
from . import (differential, enrichment, genomics, loop_calling,
               quality_control, utils, visualization)
from .config import Config, load_config
# Main imports
from .core import PlasmaFlowAnalysis
from .utils import setup_logging, validate_environment

__all__ = [
    "__version__",
    "PlasmaFlowAnalysis",
    "Config",
    "load_config",
    "setup_logging",
    "validate_environment",
    "loop_calling",
    "quality_control",
    "differential",
    "visualization",
    "genomics",
    "enrichment",
    "utils",
]


def get_info() -> Dict[str, Any]:
    """Get package information."""
    return {
        "name": "PlasmaFlow",
        "version": __version__,
        "description": "Hi-C chromatin loop analysis pipeline for plasma cell differentiation",
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "modules": __all__[6:],  # Just the module names
    }


def check_dependencies() -> Dict[str, bool]:
    """Check if key dependencies are available."""
    dependencies = {}

    # Core packages
    try:
        import numpy

        dependencies["numpy"] = True
    except ImportError:
        dependencies["numpy"] = False

    try:
        import pandas

        dependencies["pandas"] = True
    except ImportError:
        dependencies["pandas"] = False

    try:
        import cooler

        dependencies["cooler"] = True
    except ImportError:
        dependencies["cooler"] = False

    try:
        import cooltools

        dependencies["cooltools"] = True
    except ImportError:
        dependencies["cooltools"] = False

    try:
        import rpy2

        dependencies["rpy2"] = True
    except ImportError:
        dependencies["rpy2"] = False

    return dependencies


# Initialize package
logger = logging.getLogger(__name__)
logger.info(f"PlasmaFlow v{__version__} initialized")

# Check critical dependencies
deps = check_dependencies()
missing_deps = [dep for dep, available in deps.items() if not available]
if missing_deps:
    logger.warning(f"Missing dependencies: {missing_deps}")
    logger.info("Run 'pip install plasmaflow[all]' to install all dependencies")
