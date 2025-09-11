"""
Gene proximity analysis for differential chromatin loops

This module provides functionality for analyzing genes near
differential chromatin loops with TSS-based distance calculations.
"""

import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import Config
from ..utils import get_logger
from .annotations import GeneAnnotationManager, validate_gene_annotations
from .statistical import GeneDistanceStatistics

logger = get_logger(__name__)


@dataclass
class ProximityResult:
    """Result container for gene proximity analysis"""

    file_type: str  # 'category' or 'cluster'
    category: str  # 'up', 'down', 'common'
    cluster: Optional[str]
    distances_df: pd.DataFrame
    total_loops: int
    valid_loops: int
    total_genes_found: int
    failed_loops: int
    validation_results: Dict[str, Any]


class GeneProximityAnalyzer:
    """Analyzer for gene proximity to differential chromatin loops"""

    def __init__(self, config: Config):
        """
        Initialize gene proximity analyzer

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config

        # Get genomics configuration
        self.genomics_params = config.genomics.get("proximity", {})

        # Distance and regulatory zone configuration
        self.max_distance = self.genomics_params.get("max_distance", 250000)

        self.distance_bins = self.genomics_params.get(
            "distance_bins",
            [0, 500, 2000, 5000, 10000, 25000, 50000, 100000, 250000, float("inf")],
        )

        self.distance_labels = self.genomics_params.get(
            "distance_labels",
            [
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
        )

        # Validation configuration
        self.validation_config = self.genomics_params.get(
            "validation",
            {
                "check_trans_loops": True,
                "warn_trans_loops": True,
                "exclude_trans_loops": True,
                "max_trans_loops_warning": 5,
            },
        )

        # Initialize gene annotation manager
        self.gene_manager = GeneAnnotationManager(config)

        # Create output directories
        self.output_dir = Path(config.output_dir) / "genomics" / "proximity"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.csv_output_dir = self.output_dir / "csv_exports"
        self.csv_output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize statistics calculator
        self.statistics = GeneDistanceStatistics(config)

    def analyze_differential_results(
        self,
        differential_csv_files: List[Union[str, Path]],
        comparison_name: str,
        method_name: str = "differential",
    ) -> Dict[str, Any]:
        """
        Analyze gene proximity for differential loop results

        Args:
            differential_csv_files: List of CSV files with differential results
            comparison_name: Name of comparison (e.g., "prePB_vs_memB")
            method_name: Analysis method name (e.g., "diffHic")

        Returns:
            Dictionary with analysis results
        """

        logger.info(f"Starting gene proximity analysis for {comparison_name}")
        logger.info(f"Method: {method_name}, Files: {len(differential_csv_files)}")

        # Load gene annotations
        self.gene_manager.load_gene_annotations()
        genes_by_chr = self.gene_manager.genes_by_chr

        # Validate gene annotations
        validation = validate_gene_annotations(genes_by_chr)
        if not validation["validation_passed"]:
            logger.warning(
                f"Gene annotation validation failed: {validation['missing_required_columns']}"
            )

        logger.info(
            f"Loaded annotations for {validation['total_genes']} genes across {validation['chromosomes_covered']} chromosomes"
        )

        # Process each CSV file
        all_results = []
        validation_summary = {
            "total_files": len(differential_csv_files),
            "total_loops": 0,
            "total_cis_loops": 0,
            "total_trans_loops": 0,
            "files_with_trans": 0,
        }

        for csv_file in differential_csv_files:
            logger.info(f"Processing: {Path(csv_file).name}")

            try:
                result = self._analyze_csv_file(csv_file, genes_by_chr)
                if result:
                    all_results.append(result)

                    # Update validation summary
                    val_res = result.validation_results
                    validation_summary["total_loops"] += val_res["total_loops"]
                    validation_summary["total_cis_loops"] += val_res["cis_loops"]
                    validation_summary["total_trans_loops"] += val_res["trans_loops"]
                    if val_res["trans_loops"] > 0:
                        validation_summary["files_with_trans"] += 1

            except Exception as e:
                logger.error(f"Failed to process {csv_file}: {e}")

        logger.info(
            f"Successfully processed {len(all_results)}/{len(differential_csv_files)} files"
        )

        # Log validation summary
        self._log_validation_summary(validation_summary)

        # Export results
        exported_files = self._export_results(all_results, comparison_name, method_name)

        # Create visualizations and statistical analysis
        visualization_files = self.statistics.create_comprehensive_analysis(
            all_results, comparison_name, method_name, self.output_dir
        )

        # Create summary
        summary_df = self._create_analysis_summary(all_results)
        summary_file = (
            self.output_dir / f"{comparison_name}_{method_name}_proximity_summary.csv"
        )
        summary_df.to_csv(summary_file, index=False)

        logger.info(f"Gene proximity analysis completed for {comparison_name}")

        return {
            "comparison_name": comparison_name,
            "method_name": method_name,
            "results": all_results,
            "validation_summary": validation_summary,
            "exported_files": exported_files,
            "visualization_files": visualization_files,
            "summary_file": summary_file,
            "gene_annotation_summary": validation,
        }

    def _analyze_csv_file(
        self, csv_file: Union[str, Path], genes_by_chr: Dict[str, pd.DataFrame]
    ) -> Optional[ProximityResult]:
        """Analyze a single CSV file for gene proximity"""

        csv_path = Path(csv_file)

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            logger.error(f"Error reading {csv_path}: {e}")
            return None

        # Validate required columns
        required_cols = [
            "original_loop_id",
            "anchor1_chr",
            "anchor1_start",
            "anchor1_end",
            "anchor2_chr",
            "anchor2_start",
            "anchor2_end",
        ]
        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            logger.error(f"Missing columns in {csv_path}: {missing_cols}")
            return None

        # Validate loops
        validation_results = self._validate_loops(df)

        # Apply validation mask
        df_filtered = df[validation_results["valid_loops_mask"]].copy()

        # Detect file type and category
        file_info = self._detect_file_type(csv_path.name)
        if not file_info:
            return None

        file_type, category, cluster_num = file_info

        logger.info(f"Detected: {file_type} - {category} - cluster {cluster_num}")
        logger.info(
            f"Processing {len(df_filtered)} valid loops (filtered from {len(df)} total)"
        )

        # Find genes near loops
        gene_data = []
        failed_loops = 0

        for _, row in df_filtered.iterrows():
            try:
                genes = self._find_genes_near_loop(
                    str(row["anchor1_chr"]).replace("chr", ""),
                    int(row["anchor1_start"]),
                    int(row["anchor1_end"]),
                    str(row["anchor2_chr"]).replace("chr", ""),
                    int(row["anchor2_start"]),
                    int(row["anchor2_end"]),
                    genes_by_chr,
                )

                for gene in genes:
                    gene_data.append(
                        {
                            "loop_id": row["original_loop_id"],
                            "genomic_loop_id": row.get("genomic_loop_id", ""),
                            "gene_name": gene["name"],
                            "gene_symbol": gene["name2"],
                            "distance_to_loop": gene["distance"],
                            "tss_distance": gene["tss_distance"],
                            "tss_position": gene["tss_position"],
                            "strand": gene["strand"],
                            "transcript_length": gene["transcript_length"],
                            "interaction_type": gene["interaction_type"],
                            "is_within_loop": gene["is_within_loop"],
                            "category": category,
                            "cluster": cluster_num,
                            "cluster_size": row.get("cluster_size", 0),
                        }
                    )

            except Exception as e:
                failed_loops += 1
                logger.debug(f"Failed to process loop {row['original_loop_id']}: {e}")

        # Create distances DataFrame
        distances_df = pd.DataFrame(gene_data)

        if len(distances_df) > 0:
            # Add distance binning
            distances_df["distance_bin"] = pd.cut(
                distances_df["distance_to_loop"],
                bins=self.distance_bins,
                labels=self.distance_labels[: len(self.distance_bins) - 1],
                right=False,
            )

            # Add regulatory category
            distances_df = self._add_regulatory_categories(distances_df)

        return ProximityResult(
            file_type=file_type,
            category=category,
            cluster=cluster_num,
            distances_df=distances_df,
            total_loops=len(df),
            valid_loops=len(df_filtered),
            total_genes_found=len(distances_df),
            failed_loops=failed_loops,
            validation_results=validation_results,
        )

    def _validate_loops(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Validate loop data and check for trans-chromosomal loops"""

        logger.debug("Validating loop data")

        validation_results = {
            "total_loops": len(df),
            "cis_loops": 0,
            "trans_loops": 0,
            "trans_loop_details": [],
            "valid_loops_mask": None,
        }

        if not self.validation_config.get("check_trans_loops", True):
            logger.debug("Trans loop validation disabled")
            validation_results["valid_loops_mask"] = pd.Series(
                [True] * len(df), index=df.index
            )
            validation_results["cis_loops"] = len(df)
            return validation_results

        # Clean chromosome names for comparison
        anchor1_chr_clean = (
            df["anchor1_chr"].astype(str).str.replace("chr", "", regex=False)
        )
        anchor2_chr_clean = (
            df["anchor2_chr"].astype(str).str.replace("chr", "", regex=False)
        )

        # Identify cis vs trans loops
        is_cis = anchor1_chr_clean == anchor2_chr_clean
        is_trans = ~is_cis

        validation_results["cis_loops"] = is_cis.sum()
        validation_results["trans_loops"] = is_trans.sum()

        if validation_results["trans_loops"] > 0:
            trans_loops_df = df[is_trans].copy()
            trans_loops_df["anchor1_chr_clean"] = anchor1_chr_clean[is_trans]
            trans_loops_df["anchor2_chr_clean"] = anchor2_chr_clean[is_trans]

            # Store trans loop details
            for idx, row in trans_loops_df.iterrows():
                validation_results["trans_loop_details"].append(
                    {
                        "loop_id": row.get("original_loop_id", idx),
                        "anchor1_chr": row["anchor1_chr_clean"],
                        "anchor2_chr": row["anchor2_chr_clean"],
                        "anchor1_pos": f"{row['anchor1_start']}-{row['anchor1_end']}",
                        "anchor2_pos": f"{row['anchor2_start']}-{row['anchor2_end']}",
                    }
                )

            # Log warnings
            if self.validation_config.get("warn_trans_loops", True):
                logger.warning(
                    f"Found {validation_results['trans_loops']} trans-chromosomal loops!"
                )
                logger.info(
                    f"Loop breakdown: {validation_results['cis_loops']} cis, {validation_results['trans_loops']} trans"
                )

                # Show details for first few trans loops
                max_warnings = self.validation_config.get("max_trans_loops_warning", 5)
                for i, trans_loop in enumerate(
                    validation_results["trans_loop_details"][:max_warnings]
                ):
                    logger.info(
                        f"  Trans loop {trans_loop['loop_id']}: Chr{trans_loop['anchor1_chr']} ↔ Chr{trans_loop['anchor2_chr']}"
                    )

                if validation_results["trans_loops"] > max_warnings:
                    remaining = validation_results["trans_loops"] - max_warnings
                    logger.info(f"  ... and {remaining} more trans loops")
        else:
            logger.info(
                f"All {validation_results['cis_loops']} loops are cis-chromosomal"
            )

        # Create mask for valid loops
        if (
            self.validation_config.get("exclude_trans_loops", True)
            and validation_results["trans_loops"] > 0
        ):
            validation_results["valid_loops_mask"] = is_cis
            logger.info(
                f"Excluding {validation_results['trans_loops']} trans loops from analysis"
            )
            logger.info(f"Proceeding with {validation_results['cis_loops']} cis loops")
        else:
            validation_results["valid_loops_mask"] = pd.Series(
                [True] * len(df), index=df.index
            )
            logger.info(
                f"Proceeding with all {validation_results['total_loops']} loops"
            )

        return validation_results

    def _detect_file_type(self, filename: str) -> Optional[Tuple[str, str, str]]:
        """Detect file type, category, and cluster from filename"""

        if "_up_clusters_" in filename:
            return ("category", "up", "all")
        elif "_down_clusters_" in filename:
            return ("category", "down", "all")
        elif "_common_clusters_" in filename:
            return ("category", "common", "all")
        elif "cluster" in filename:
            parts = filename.replace("_ENHANCED.csv", "").split("_")
            category = parts[-1] if parts[-1] in ["up", "down", "common"] else "unknown"

            # Find cluster number
            cluster_num = "all"
            for part in parts:
                if part.startswith("cluster"):
                    cluster_num = part.replace("cluster", "")
                    break

            return ("cluster", category, cluster_num)
        else:
            logger.warning(f"Cannot determine file type for {filename}")
            return None

    def _find_genes_near_loop(
        self,
        anchor1_chr: str,
        anchor1_start: int,
        anchor1_end: int,
        anchor2_chr: str,
        anchor2_start: int,
        anchor2_end: int,
        genes_by_chr: Dict[str, pd.DataFrame],
    ) -> List[Dict[str, Any]]:
        """Find genes near a chromatin loop using TSS-based distances"""

        try:
            anchor1_center = (anchor1_start + anchor1_end) // 2
            anchor2_center = (anchor2_start + anchor2_end) // 2
        except (ValueError, TypeError):
            return []

        all_genes = []

        # Same chromosome analysis
        if anchor1_chr == anchor2_chr and anchor1_chr in genes_by_chr:
            chr_genes = genes_by_chr[anchor1_chr]

            # Calculate loop boundaries
            loop_start = min(anchor1_center, anchor2_center)
            loop_end = max(anchor1_center, anchor2_center)

            for _, gene in chr_genes.iterrows():
                tss_pos = gene["TSS"]

                # Calculate TSS distances to both anchors
                dist_to_anchor1 = abs(tss_pos - anchor1_center)
                dist_to_anchor2 = abs(tss_pos - anchor2_center)
                min_tss_distance = min(dist_to_anchor1, dist_to_anchor2)

                # Check if TSS is within loop span
                is_within_loop = loop_start <= tss_pos <= loop_end

                # Determine regulatory distance and type
                if is_within_loop:
                    regulatory_distance = 0
                    interaction_type = "within_loop"
                    include_gene = True
                else:
                    regulatory_distance = min_tss_distance
                    interaction_type = "flanking"
                    include_gene = regulatory_distance <= self.max_distance

                if include_gene:
                    all_genes.append(
                        {
                            "name": gene["name"],
                            "name2": gene["name2"],
                            "distance": regulatory_distance,
                            "tss_distance": min_tss_distance,
                            "tss_position": tss_pos,
                            "strand": gene["strand"],
                            "transcript_length": gene["transcript_length"],
                            "interaction_type": interaction_type,
                            "is_within_loop": is_within_loop,
                            "loop_span": loop_end - loop_start,
                        }
                    )

        else:
            # Trans-chromosomal analysis
            for anchor_chr, anchor_center in [
                (anchor1_chr, anchor1_center),
                (anchor2_chr, anchor2_center),
            ]:
                if anchor_chr in genes_by_chr:
                    chr_genes = genes_by_chr[anchor_chr]

                    for _, gene in chr_genes.iterrows():
                        tss_distance = abs(gene["TSS"] - anchor_center)

                        if tss_distance <= self.max_distance:
                            all_genes.append(
                                {
                                    "name": gene["name"],
                                    "name2": gene["name2"],
                                    "distance": tss_distance,
                                    "tss_distance": tss_distance,
                                    "tss_position": gene["TSS"],
                                    "strand": gene["strand"],
                                    "transcript_length": gene["transcript_length"],
                                    "interaction_type": "trans_flanking",
                                    "is_within_loop": False,
                                    "anchor_chr": anchor_chr,
                                }
                            )

        # Remove duplicates, keeping closest TSS distance
        gene_dict = {}
        for gene in all_genes:
            name = gene["name2"]
            if name not in gene_dict or gene["distance"] < gene_dict[name]["distance"]:
                gene_dict[name] = gene

        return list(gene_dict.values())

    def _add_regulatory_categories(self, distances_df: pd.DataFrame) -> pd.DataFrame:
        """Add regulatory category classifications"""

        distances_df["regulatory_category"] = "Unknown"

        for idx, row in distances_df.iterrows():
            if row.get("is_within_loop", False):
                distances_df.loc[idx, "regulatory_category"] = "Within-loop"
            else:
                dist = row["distance_to_loop"]
                if dist <= 500:
                    distances_df.loc[idx, "regulatory_category"] = "TSS-proximal"
                elif dist <= 2000:
                    distances_df.loc[idx, "regulatory_category"] = "Promoter"
                elif dist <= 5000:
                    distances_df.loc[idx, "regulatory_category"] = "Local-enhancer"
                elif dist <= 10000:
                    distances_df.loc[idx, "regulatory_category"] = "Nearby-enhancer"
                elif dist <= 25000:
                    distances_df.loc[idx, "regulatory_category"] = "Distal-enhancer"
                elif dist <= 50000:
                    distances_df.loc[idx, "regulatory_category"] = "Long-range"
                elif dist <= 100000:
                    distances_df.loc[idx, "regulatory_category"] = "Very-long-range"
                elif dist <= 250000:
                    distances_df.loc[idx, "regulatory_category"] = "Distant"
                else:
                    distances_df.loc[idx, "regulatory_category"] = "Very-distant"

        return distances_df

    def _log_validation_summary(self, validation_summary: Dict[str, Any]):
        """Log validation summary"""

        logger.info("Loop Validation Summary:")
        logger.info(f"  Total loops: {validation_summary['total_loops']}")
        logger.info(f"  Cis loops: {validation_summary['total_cis_loops']}")
        logger.info(f"  Trans loops: {validation_summary['total_trans_loops']}")
        logger.info(
            f"  Files with trans loops: {validation_summary['files_with_trans']}/{validation_summary['total_files']}"
        )

        if validation_summary["total_trans_loops"] > 0:
            trans_pct = (
                validation_summary["total_trans_loops"]
                / validation_summary["total_loops"]
            ) * 100
            logger.info(f"  Trans loop percentage: {trans_pct:.2f}%")
        else:
            logger.info("  All loops are cis-chromosomal!")

    def _export_results(
        self, results: List[ProximityResult], comparison_name: str, method_name: str
    ) -> Dict[str, Path]:
        """Export analysis results to CSV files"""

        logger.info("Exporting gene proximity results")

        exported_files = {}

        for result in results:
            if len(result.distances_df) > 0:
                if result.file_type == "category":
                    filename = f"{comparison_name}_{method_name}_{result.category}_genes_enhanced.csv"
                else:
                    filename = f"{comparison_name}_{method_name}_{result.category}_cluster{result.cluster}_genes_enhanced.csv"

                filepath = self.csv_output_dir / filename
                result.distances_df.to_csv(filepath, index=False)
                exported_files[f"{result.category}_{result.cluster}"] = filepath

                logger.debug(f"Exported: {filename}")

        logger.info(f"Exported {len(exported_files)} CSV files")
        return exported_files

    def _create_analysis_summary(self, results: List[ProximityResult]) -> pd.DataFrame:
        """Create summary DataFrame of analysis results"""

        summary_data = []

        for result in results:
            summary_data.append(
                {
                    "file_type": result.file_type,
                    "category": result.category,
                    "cluster": result.cluster,
                    "total_loops": result.total_loops,
                    "valid_loops": result.valid_loops,
                    "total_genes_found": result.total_genes_found,
                    "failed_loops": result.failed_loops,
                    "cis_loops": result.validation_results.get("cis_loops", 0),
                    "trans_loops": result.validation_results.get("trans_loops", 0),
                    "genes_within_loop": (
                        len(result.distances_df[result.distances_df["is_within_loop"]])
                        if len(result.distances_df) > 0
                        else 0
                    ),
                    "median_distance": (
                        result.distances_df["distance_to_loop"].median()
                        if len(result.distances_df) > 0
                        else np.nan
                    ),
                    "mean_distance": (
                        result.distances_df["distance_to_loop"].mean()
                        if len(result.distances_df) > 0
                        else np.nan
                    ),
                }
            )

        return pd.DataFrame(summary_data)


def find_genes_near_loops(
    loop_coordinates: pd.DataFrame,
    genes_by_chr: Dict[str, pd.DataFrame],
    max_distance: int = 250000,
) -> pd.DataFrame:
    """
    Convenience function to find genes near loops

    Args:
        loop_coordinates: DataFrame with loop coordinates
        genes_by_chr: Dictionary of chromosome -> genes DataFrame
        max_distance: Maximum distance to consider (bp)

    Returns:
        DataFrame with gene-loop proximity results
    """

    required_cols = [
        "anchor1_chr",
        "anchor1_start",
        "anchor1_end",
        "anchor2_chr",
        "anchor2_start",
        "anchor2_end",
    ]

    if not all(col in loop_coordinates.columns for col in required_cols):
        raise ValueError(f"Loop coordinates missing required columns: {required_cols}")

    # Initialize analyzer with default config
    from ..config import get_default_config

    config = get_default_config()
    analyzer = GeneProximityAnalyzer(config)

    # Process each loop
    all_gene_data = []

    for idx, loop in loop_coordinates.iterrows():
        genes = analyzer._find_genes_near_loop(
            str(loop["anchor1_chr"]).replace("chr", ""),
            int(loop["anchor1_start"]),
            int(loop["anchor1_end"]),
            str(loop["anchor2_chr"]).replace("chr", ""),
            int(loop["anchor2_start"]),
            int(loop["anchor2_end"]),
            genes_by_chr,
        )

        for gene in genes:
            gene_data = gene.copy()
            gene_data["loop_index"] = idx

            # Add loop information if available
            for col in loop_coordinates.columns:
                if col not in gene_data:
                    gene_data[f"loop_{col}"] = loop[col]

            all_gene_data.append(gene_data)

    return pd.DataFrame(all_gene_data)
