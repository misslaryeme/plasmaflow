"""
Core configuration management for PlasmaFlow
"""

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

logger = logging.getLogger(__name__)


@dataclass
class Config:
    """Main configuration class for PlasmaFlow analysis"""

    # General settings
    project_name: str = "PlasmaFlow_Analysis"
    resolution: int = 10000
    random_seed: int = 42
    n_threads: int = 4

    # Input/Output paths
    input_dir: Optional[str] = None
    output_dir: Optional[str] = None
    temp_dir: Optional[str] = None

    # Sample configuration
    samples: Dict[str, Any] = field(default_factory=dict)
    cell_types: List[str] = field(default_factory=lambda: ["memB", "prePB", "PB", "PC"])

    # Analysis parameters
    loop_calling: Dict[str, Any] = field(default_factory=dict)
    quality_control: Dict[str, Any] = field(default_factory=dict)
    differential: Dict[str, Any] = field(default_factory=dict)
    visualization: Dict[str, Any] = field(default_factory=dict)
    genomics: Dict[str, Any] = field(default_factory=dict)
    enrichment: Dict[str, Any] = field(default_factory=dict)

    # R configuration
    r_config: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize default configurations"""
        if not self.loop_calling:
            self.loop_calling = self._get_default_loop_calling()
        if not self.quality_control:
            self.quality_control = self._get_default_quality_control()
        if not self.differential:
            self.differential = self._get_default_differential()
        if not self.visualization:
            self.visualization = self._get_default_visualization()
        if not self.genomics:
            self.genomics = self._get_default_genomics()
        if not self.enrichment:
            self.enrichment = self._get_default_enrichment()
        if not self.r_config:
            self.r_config = self._get_default_r_config()

    def _get_default_loop_calling(self) -> Dict[str, Any]:
        """Default loop calling configuration"""
        return {
            "threshold": 0.95,
            "weight_name": "weight",
            "model_dir": "/path/to/peakachu/models",
            "models": {
                "memB": {"model": "100million", "reads": 113922292},
                "PC": {"model": "100million", "reads": 91590221},
                "PB": {"model": "150million", "reads": 158184861},
                "prePB": {"model": "150million", "reads": 159514688},
            },
        }

    def _get_default_quality_control(self) -> Dict[str, Any]:
        """Default quality control configuration"""
        return {
            "apa": {"flank": 200000, "min_diag": 3, "maxdist": 1000000, "nproc": 4},
            "cache_results": True,
            "recompute_if_params_changed": True,
        }

    def _get_default_differential(self) -> Dict[str, Any]:
        """Default differential analysis configuration"""
        return {
            "methods": ["diffHic", "DESeq2", "edgeR"],
            "comparisons": [
                {"control": "memB", "treatment": "prePB", "name": "prePB_vs_memB"}
            ],
            "fdr_threshold": 0.05,
            "logfc_threshold": 1.0,
            "min_counts": 10,
            "filter_params": {
                "min_interaction_distance": 20000,
                "max_interaction_distance": 2000000,
            },
        }

    def _get_default_visualization(self) -> Dict[str, Any]:
        """Default visualization configuration"""
        return {
            "heatmaps": {
                "flanking_region": 3000,
                "bin_size": 100,
                "colormap": "RdBu_r",
                "vmin": -3,
                "vmax": 3,
            },
            "volcano_plots": {
                "fdr_threshold": 0.05,
                "logfc_threshold": 1.0,
                "point_size": 0.5,
                "alpha": 0.6,
            },
            "save_formats": ["png", "pdf", "svg"],
        }

    def _get_default_genomics(self) -> Dict[str, Any]:
        """Default genomics configuration"""
        return {
            "proximity": {
                "max_distance": 250000,
                "distance_bins": [
                    0,
                    500,
                    2000,
                    5000,
                    10000,
                    25000,
                    50000,
                    100000,
                    250000,
                    float("inf"),
                ],
                "distance_labels": [
                    "TSS-proximal\n(0-500bp)",
                    "Promoter\n(500bp-2kb)",
                    "Local-enhancer\n(2-5kb)",
                    "Nearby-enhancer\n(5-10kb)",
                    "Distal-enhancer\n(10-25kb)",
                    "Long-range\n(25-50kb)",
                    "Very-long-range\n(50-100kb)",
                    "Distant\n(100-250kb)",
                    "Very-distant\n(>250kb)",
                ],
                "transcript_filters": {
                    "coding_only": True,
                    "allowed_prefixes": ["NM_"],
                    "excluded_prefixes": ["NR_"],
                },
                "deduplication": {
                    "location_threshold": 100,
                    "selection_method": "longest",
                },
                "validation": {
                    "check_trans_loops": True,
                    "warn_trans_loops": True,
                    "exclude_trans_loops": True,
                    "max_trans_loops_warning": 5,
                },
            },
            "expression_analysis": {
                "fdr_threshold": 0.05,
                "logfc_threshold": 0.5,
                "min_expression": 1.0,
                "expression_methods": ["DESeq2", "edgeR"],
                "volcano_plot_params": {
                    "point_size": 1.0,
                    "alpha": 0.6,
                    "label_top_genes": 20,
                },
            },
            "visualization": {
                "category_colors": {
                    "up": "#E69F00",
                    "down": "#56B4E9",
                    "common": "#009E73",
                },
                "cluster_colors": {
                    "1": "#E69F00",
                    "2": "#56B4E9",
                    "3": "#009E73",
                    "4": "#F0E442",
                    "5": "#0072B2",
                    "6": "#D55E00",
                    "7": "#CC79A7",
                },
                "statistical": {
                    "confidence_level": 0.95,
                    "n_bootstrap": 1000,
                    "alpha": 0.05,
                    "multiple_testing_correction": "fdr_bh",
                },
            },
            "reference_files": {
                "annotation_source": "ucsc_hg19",
                "gtf": "/path/to/genes.gtf",
                "genome_fasta": "/path/to/genome.fa",
            },
        }

    def _get_default_enrichment(self) -> Dict[str, Any]:
        """Default enrichment analysis configuration"""
        return {
            "gsea": {
                "gene_sets": ["C2", "C5", "HALLMARK"],
                "permutations": 1000,
                "min_size": 15,
                "max_size": 500,
            },
            "ora": {
                "gene_sets": ["C2", "C5", "HALLMARK"],
                "background": "genome",
                "fdr_threshold": 0.05,
            },
            "msigdb_species": "Homo sapiens",
            "msigdb_version": "7.5.1",
        }

    def _get_default_r_config(self) -> Dict[str, Any]:
        """Default R configuration"""
        return {
            "r_home": None,  # Auto-detect
            "bioconductor_repos": [
                "https://bioconductor.org/packages/release/bioc",
                "https://cloud.r-project.org",
            ],
            "required_packages": [
                "diffHic",
                "InteractionSet",
                "GenomicRanges",
                "DESeq2",
                "edgeR",
                "clusterProfiler",
                "msigdbr",
                "ggplot2",
                "dplyr",
            ],
        }


def load_config(config_file: Union[str, Path]) -> Config:
    """Load configuration from YAML or JSON file"""
    config_path = Path(config_file)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    logger.info(f"Loading configuration from {config_path}")

    with open(config_path, "r") as f:
        if config_path.suffix.lower() in [".yaml", ".yml"]:
            config_dict = yaml.safe_load(f)
        elif config_path.suffix.lower() == ".json":
            config_dict = json.load(f)
        else:
            raise ValueError(f"Unsupported config file format: {config_path.suffix}")

    return Config(**config_dict)


def save_config(config: Config, output_file: Union[str, Path]) -> None:
    """Save configuration to YAML file"""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert config to dict
    config_dict = {
        "project_name": config.project_name,
        "resolution": config.resolution,
        "random_seed": config.random_seed,
        "n_threads": config.n_threads,
        "input_dir": config.input_dir,
        "output_dir": config.output_dir,
        "temp_dir": config.temp_dir,
        "samples": config.samples,
        "cell_types": config.cell_types,
        "loop_calling": config.loop_calling,
        "quality_control": config.quality_control,
        "differential": config.differential,
        "visualization": config.visualization,
        "genomics": config.genomics,
        "enrichment": config.enrichment,
        "r_config": config.r_config,
    }

    with open(output_path, "w") as f:
        yaml.dump(config_dict, f, default_flow_style=False, indent=2)

    logger.info(f"Configuration saved to {output_path}")


def validate_config(config: Config) -> List[str]:
    """Validate configuration and return list of issues"""
    issues = []

    # Check required paths exist
    if config.input_dir and not Path(config.input_dir).exists():
        issues.append(f"Input directory does not exist: {config.input_dir}")

    if config.output_dir:
        try:
            Path(config.output_dir).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            issues.append(f"Cannot create output directory {config.output_dir}: {e}")

    # Validate resolution
    if config.resolution <= 0:
        issues.append("Resolution must be positive")

    # Validate thread count
    if config.n_threads <= 0:
        issues.append("Number of threads must be positive")

    # Validate cell types
    if not config.cell_types:
        issues.append("At least one cell type must be specified")

    # Validate differential analysis comparisons
    if config.differential.get("comparisons"):
        for comp in config.differential["comparisons"]:
            if not all(k in comp for k in ["control", "treatment", "name"]):
                issues.append(
                    "Each comparison must have 'control', 'treatment', and 'name'"
                )

    return issues


def get_default_config() -> Config:
    """Get default configuration object"""
    return Config()
