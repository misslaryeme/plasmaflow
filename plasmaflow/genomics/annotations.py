"""
Gene annotation management for PlasmaFlow

This module handles gene annotation loading, filtering, and processing
for gene proximity analysis.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..config import Config
from ..utils import get_logger

logger = get_logger(__name__)


class GeneAnnotationManager:
    """Manager for gene annotations with configurable filtering and processing"""

    def __init__(self, config: Config):
        """
        Initialize gene annotation manager

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config

        # Get genomics configuration
        self.genomics_params = config.genomics.get("proximity", {})

        # Default transcript filters
        self.transcript_filters = self.genomics_params.get(
            "transcript_filters",
            {
                "coding_only": True,
                "allowed_prefixes": ["NM_"],
                "excluded_prefixes": ["NR_"],
            },
        )

        # Default deduplication config
        self.dedup_config = self.genomics_params.get(
            "deduplication", {"location_threshold": 100, "selection_method": "longest"}
        )

        # Cache for loaded annotations
        self.genes_by_chr = {}
        self._annotation_loaded = False

        # Output directory for annotations
        self.output_dir = Path(config.output_dir) / "genomics" / "annotations"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def load_gene_annotations(
        self, annotation_source: str = "ucsc_hg19", force_reload: bool = False
    ) -> Dict[str, pd.DataFrame]:
        """
        Load and process gene annotations

        Args:
            annotation_source: Source of annotations ("ucsc_hg19", "ensembl", or file path)
            force_reload: Force reload even if already cached

        Returns:
            Dictionary mapping chromosome names to gene DataFrames
        """

        if self._annotation_loaded and not force_reload:
            logger.info("Using cached gene annotations")
            return self.genes_by_chr

        logger.info(f"Loading gene annotations from {annotation_source}")

        if annotation_source == "ucsc_hg19":
            self.genes_by_chr = self._load_ucsc_hg19()
        elif annotation_source.startswith(("http://", "https://", "ftp://")):
            self.genes_by_chr = self._load_from_url(annotation_source)
        elif Path(annotation_source).exists():
            self.genes_by_chr = self._load_from_file(annotation_source)
        else:
            raise ValueError(f"Unknown annotation source: {annotation_source}")

        self._annotation_loaded = True

        # Save processed annotations
        self._save_processed_annotations()

        logger.info(f"Loaded annotations for {len(self.genes_by_chr)} chromosomes")
        return self.genes_by_chr

    def _load_ucsc_hg19(self) -> Dict[str, pd.DataFrame]:
        """Load UCSC hg19 gene annotations"""
        logger.info("Loading UCSC hg19 gene annotations")

        url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz"

        try:
            # Load required columns
            genes_df = pd.read_csv(
                url,
                sep="\t",
                compression="gzip",
                header=None,
                usecols=[1, 2, 3, 4, 5, 12],
                names=["name", "chrom", "strand", "txStart", "txEnd", "name2"],
            )

            # Save raw annotations
            raw_file = self.output_dir / "refGene_hg19_raw.csv"
            genes_df.to_csv(raw_file, index=False)
            logger.debug(f"Saved raw annotations to {raw_file}")

            logger.info(f"Initial annotations loaded: {len(genes_df)} entries")

            # Process annotations
            return self._process_annotations(genes_df)

        except Exception as e:
            logger.error(f"Failed to load UCSC hg19 annotations: {e}")
            raise

    def _load_from_url(self, url: str) -> Dict[str, pd.DataFrame]:
        """Load annotations from URL"""
        logger.info(f"Loading annotations from URL: {url}")

        try:
            genes_df = pd.read_csv(url)
            return self._process_annotations(genes_df)
        except Exception as e:
            logger.error(f"Failed to load annotations from URL: {e}")
            raise

    def _load_from_file(self, file_path: str) -> Dict[str, pd.DataFrame]:
        """Load annotations from file"""
        logger.info(f"Loading annotations from file: {file_path}")

        try:
            genes_df = pd.read_csv(file_path)
            return self._process_annotations(genes_df)
        except Exception as e:
            logger.error(f"Failed to load annotations from file: {e}")
            raise

    def _process_annotations(self, genes_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Process raw gene annotations with filtering and deduplication"""

        logger.info("Processing gene annotations")

        # Step 1: Transcript filtering
        genes_df = self._apply_transcript_filters(genes_df)

        # Step 2: Chromosome filtering
        genes_df = self._filter_chromosomes(genes_df)

        # Step 3: Calculate genomic metrics
        genes_df = self._calculate_genomic_metrics(genes_df)

        # Step 4: Deduplication
        genes_df = self._deduplicate_genes(genes_df)

        # Step 5: Location-based deduplication
        genes_df = self._remove_redundant_locations(genes_df)

        # Step 6: Organize by chromosome
        genes_by_chr = self._organize_by_chromosome(genes_df)

        # Log summary statistics
        self._log_processing_summary(genes_df, genes_by_chr)

        return genes_by_chr

    def _apply_transcript_filters(self, genes_df: pd.DataFrame) -> pd.DataFrame:
        """Apply transcript type filtering"""

        logger.info("Applying transcript filters")

        if self.transcript_filters.get("coding_only", True):
            allowed_prefixes = self.transcript_filters.get("allowed_prefixes", ["NM_"])
            logger.info(f"Filtering for coding transcripts: {allowed_prefixes}")

            allowed_mask = genes_df["name"].str.startswith(tuple(allowed_prefixes))
            genes_df = genes_df[allowed_mask].copy()

        # Apply exclusions
        excluded_prefixes = self.transcript_filters.get("excluded_prefixes", [])
        if excluded_prefixes:
            logger.info(f"Excluding prefixes: {excluded_prefixes}")
            excluded_mask = ~genes_df["name"].str.startswith(tuple(excluded_prefixes))
            genes_df = genes_df[excluded_mask].copy()

        logger.info(f"After transcript filtering: {len(genes_df)} entries")
        return genes_df

    def _filter_chromosomes(self, genes_df: pd.DataFrame) -> pd.DataFrame:
        """Filter for standard chromosomes"""

        logger.info("Filtering chromosomes")

        # Clean chromosome names
        genes_df["chrom"] = genes_df["chrom"].str.replace("chr", "", regex=False)

        # Standard chromosomes
        standard_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
        genes_df = genes_df[genes_df["chrom"].isin(standard_chroms)].copy()

        logger.info(f"After chromosome filtering: {len(genes_df)} entries")
        return genes_df

    def _calculate_genomic_metrics(self, genes_df: pd.DataFrame) -> pd.DataFrame:
        """Calculate genomic metrics including TSS positions"""

        logger.info("Calculating genomic metrics")

        # Basic metrics
        genes_df["transcript_length"] = genes_df["txEnd"] - genes_df["txStart"]
        genes_df["gene_center"] = (genes_df["txStart"] + genes_df["txEnd"]) // 2

        # Strand-aware TSS calculation
        genes_df["TSS"] = np.where(
            genes_df["strand"] == "+",
            genes_df["txStart"],  # + strand: TSS at start
            genes_df["txEnd"],  # - strand: TSS at end
        )

        # Analysis of TSS vs gene center
        genes_df["tss_center_diff"] = np.abs(genes_df["TSS"] - genes_df["gene_center"])

        logger.info("TSS vs Gene Center Analysis:")
        logger.info(f"  Mean difference: {genes_df['tss_center_diff'].mean():.0f} bp")
        logger.info(
            f"  Median difference: {genes_df['tss_center_diff'].median():.0f} bp"
        )
        logger.info(f"  Max difference: {genes_df['tss_center_diff'].max():.0f} bp")
        logger.info(
            f"  Genes with >10kb difference: {(genes_df['tss_center_diff'] > 10000).sum()}"
        )

        # Strand distribution
        strand_counts = genes_df["strand"].value_counts()
        logger.info(f"Strand distribution: {strand_counts.to_dict()}")

        return genes_df

    def _deduplicate_genes(self, genes_df: pd.DataFrame) -> pd.DataFrame:
        """Deduplicate genes by symbol using configured method"""

        method = self.dedup_config.get("selection_method", "longest")
        logger.info(f"Deduplicating by gene symbol using '{method}' method")

        if method == "longest":
            # Select transcript with maximum length
            idx_selected = genes_df.groupby("name2")["transcript_length"].idxmax()

        elif method == "first":
            # Select first occurrence
            idx_selected = genes_df.groupby("name2").first().index

        elif method == "canonical":
            # Prefer canonical transcripts (NM_ prefixed)
            genes_df["is_canonical"] = genes_df["name"].str.startswith("NM_")
            genes_df = genes_df.sort_values(
                ["name2", "is_canonical", "transcript_length"],
                ascending=[True, False, False],
            )
            idx_selected = genes_df.groupby("name2").first().index

        else:
            raise ValueError(f"Unknown selection method: {method}")

        genes_df_selected = genes_df.loc[idx_selected].copy()
        logger.info(
            f"After gene symbol deduplication: {len(genes_df_selected)} entries"
        )

        return genes_df_selected

    def _remove_redundant_locations(self, genes_df: pd.DataFrame) -> pd.DataFrame:
        """Remove redundant locations using TSS-based distance threshold"""

        threshold = self.dedup_config.get("location_threshold", 100)
        logger.info(f"Removing redundant locations (TSS threshold: {threshold}bp)")

        final_entries = []

        for chrom in genes_df["chrom"].unique():
            chr_df = genes_df[genes_df["chrom"] == chrom].copy()
            chr_df = chr_df.sort_values("TSS")

            for _, row in chr_df.iterrows():
                gene_symbol = row["name2"]
                tss_pos = row["TSS"]

                # Check for redundant entries
                is_redundant = False
                for processed_entry in final_entries:
                    if (
                        processed_entry["name2"] == gene_symbol
                        and processed_entry["chrom"] == chrom
                        and abs(processed_entry["TSS"] - tss_pos) <= threshold
                    ):
                        is_redundant = True
                        break

                if not is_redundant:
                    final_entries.append(row.to_dict())

        final_df = pd.DataFrame(final_entries)
        logger.info(f"After TSS-based location deduplication: {len(final_df)} entries")

        return final_df

    def _organize_by_chromosome(
        self, genes_df: pd.DataFrame
    ) -> Dict[str, pd.DataFrame]:
        """Organize genes by chromosome"""

        logger.info("Organizing genes by chromosome")

        genes_by_chr = {}
        standard_chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]

        for chrom in standard_chroms:
            chr_genes = genes_df[genes_df["chrom"] == chrom]
            if len(chr_genes) > 0:
                genes_by_chr[chrom] = chr_genes.reset_index(drop=True)

        return genes_by_chr

    def _log_processing_summary(
        self, genes_df: pd.DataFrame, genes_by_chr: Dict[str, pd.DataFrame]
    ):
        """Log processing summary statistics"""

        logger.info("Gene Annotation Processing Summary:")
        logger.info(f"  Configuration:")
        logger.info(
            f"    Coding only: {self.transcript_filters.get('coding_only', True)}"
        )
        logger.info(
            f"    Allowed prefixes: {self.transcript_filters.get('allowed_prefixes', [])}"
        )
        logger.info(
            f"    Excluded prefixes: {self.transcript_filters.get('excluded_prefixes', [])}"
        )
        logger.info(
            f"    Selection method: {self.dedup_config.get('selection_method', 'longest')}"
        )
        logger.info(
            f"    Location threshold: {self.dedup_config.get('location_threshold', 100)}bp"
        )

        logger.info(f"  Results:")
        logger.info(f"    Total unique genes: {len(genes_df)}")
        logger.info(f"    Chromosomes covered: {len(genes_by_chr)}")
        logger.info(
            f"    Average transcript length: {genes_df['transcript_length'].mean():.0f}bp"
        )
        logger.info(
            f"    Median transcript length: {genes_df['transcript_length'].median():.0f}bp"
        )

        logger.info(f"  Genes per chromosome:")
        for chrom in sorted(
            genes_by_chr.keys(), key=lambda x: int(x) if x.isdigit() else 100
        ):
            logger.info(f"    Chr {chrom}: {len(genes_by_chr[chrom])} genes")

    def _save_processed_annotations(self):
        """Save processed annotations for future use"""

        output_file = self.output_dir / "processed_gene_annotations.csv"

        # Combine all chromosome data
        all_genes = []
        for chrom, genes_df in self.genes_by_chr.items():
            all_genes.append(genes_df)

        if all_genes:
            combined_df = pd.concat(all_genes, ignore_index=True)
            combined_df.to_csv(output_file, index=False)
            logger.info(f"Saved processed annotations to {output_file}")

    def get_genes_for_chromosome(self, chromosome: str) -> Optional[pd.DataFrame]:
        """
        Get genes for a specific chromosome

        Args:
            chromosome: Chromosome name (e.g., '1', '2', 'X', 'Y')

        Returns:
            DataFrame with genes for the chromosome, or None if not found
        """

        # Clean chromosome name
        chrom_clean = str(chromosome).replace("chr", "")
        return self.genes_by_chr.get(chrom_clean)

    def get_gene_info(self, gene_symbol: str) -> Optional[Dict]:
        """
        Get information for a specific gene symbol

        Args:
            gene_symbol: Gene symbol to search for

        Returns:
            Dictionary with gene information, or None if not found
        """

        for chrom, genes_df in self.genes_by_chr.items():
            gene_matches = genes_df[genes_df["name2"] == gene_symbol]
            if len(gene_matches) > 0:
                gene_info = gene_matches.iloc[0].to_dict()
                gene_info["chromosome"] = chrom
                return gene_info

        return None

    def get_total_gene_count(self) -> int:
        """Get total number of unique genes"""
        return sum(len(genes_df) for genes_df in self.genes_by_chr.values())

    def get_chromosome_coverage(self) -> List[str]:
        """Get list of chromosomes with gene annotations"""
        return list(self.genes_by_chr.keys())


def validate_gene_annotations(genes_by_chr: Dict[str, pd.DataFrame]) -> Dict[str, any]:
    """
    Validate gene annotation data

    Args:
        genes_by_chr: Dictionary of chromosome -> genes DataFrame

    Returns:
        Dictionary with validation results
    """

    validation_results = {
        "total_genes": 0,
        "chromosomes_covered": len(genes_by_chr),
        "missing_required_columns": [],
        "genes_per_chromosome": {},
        "strand_distribution": {},
        "length_statistics": {},
        "validation_passed": True,
    }

    required_columns = ["name", "name2", "chrom", "strand", "txStart", "txEnd", "TSS"]

    for chrom, genes_df in genes_by_chr.items():
        # Check required columns
        missing_cols = [col for col in required_columns if col not in genes_df.columns]
        if missing_cols:
            validation_results["missing_required_columns"].extend(missing_cols)
            validation_results["validation_passed"] = False

        # Count genes
        n_genes = len(genes_df)
        validation_results["total_genes"] += n_genes
        validation_results["genes_per_chromosome"][chrom] = n_genes

        # Strand distribution
        if "strand" in genes_df.columns:
            strand_counts = genes_df["strand"].value_counts()
            validation_results["strand_distribution"][chrom] = strand_counts.to_dict()

        # Length statistics
        if "transcript_length" in genes_df.columns:
            validation_results["length_statistics"][chrom] = {
                "mean": genes_df["transcript_length"].mean(),
                "median": genes_df["transcript_length"].median(),
                "min": genes_df["transcript_length"].min(),
                "max": genes_df["transcript_length"].max(),
            }

    # Remove duplicates from missing columns
    validation_results["missing_required_columns"] = list(
        set(validation_results["missing_required_columns"])
    )

    return validation_results
