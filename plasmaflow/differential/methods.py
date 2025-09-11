"""
Individual differential analysis methods (diffHic, DESeq2, edgeR)
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..utils import RInterface, get_logger
from .analyzer import DifferentialResult

logger = get_logger(__name__)


class BaseAnalyzer(ABC):
    """Base class for differential analysis methods"""

    def __init__(self, config, r_interface: RInterface):
        self.config = config
        self.r_interface = r_interface
        self.diff_params = config.differential

    @abstractmethod
    def run_analysis(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> DifferentialResult:
        """Run the differential analysis"""
        pass

    def _classify_results(
        self,
        results_df: pd.DataFrame,
        fdr_threshold: float = 0.05,
        logfc_threshold: float = 1.0,
    ) -> pd.DataFrame:
        """Classify results into up/down/not significant"""

        results_df = results_df.copy()

        # Initialize regulation column
        results_df["regulation"] = "Not Significant"

        # Classify based on FDR and log fold change
        up_mask = (results_df["FDR"] <= fdr_threshold) & (
            results_df["logFC"] > logfc_threshold
        )
        down_mask = (results_df["FDR"] <= fdr_threshold) & (
            results_df["logFC"] < -logfc_threshold
        )

        results_df.loc[up_mask, "regulation"] = "Up-regulated"
        results_df.loc[down_mask, "regulation"] = "Down-regulated"

        # Add confidence levels
        results_df["confidence"] = "Low"
        results_df.loc[results_df["FDR"] <= 0.05, "confidence"] = "Medium"
        results_df.loc[results_df["FDR"] <= 0.01, "confidence"] = "High"
        results_df.loc[results_df["FDR"] <= 0.001, "confidence"] = "Very High"

        # Add effect size categories
        results_df["effect_size"] = "Small"
        results_df.loc[np.abs(results_df["logFC"]) > 0.58, "effect_size"] = (
            "Medium"  # ~1.5-fold
        )
        results_df.loc[np.abs(results_df["logFC"]) > 1.0, "effect_size"] = (
            "Large"  # 2-fold
        )
        results_df.loc[np.abs(results_df["logFC"]) > 1.58, "effect_size"] = (
            "Very Large"  # 3-fold
        )

        return results_df

    def _calculate_statistics(self, results_df: pd.DataFrame) -> Dict[str, Any]:
        """Calculate summary statistics"""

        stats = {
            "n_tested": len(results_df),
            "n_up_regulated": sum(results_df["regulation"] == "Up-regulated"),
            "n_down_regulated": sum(results_df["regulation"] == "Down-regulated"),
            "min_fdr": results_df["FDR"].min() if len(results_df) > 0 else None,
            "max_abs_logfc": (
                np.abs(results_df["logFC"]).max() if len(results_df) > 0 else None
            ),
        }

        stats["n_significant"] = stats["n_up_regulated"] + stats["n_down_regulated"]

        return stats


class DiffHicAnalyzer(BaseAnalyzer):
    """diffHic-based differential analysis using R"""

    def __init__(self, config, r_interface: RInterface):
        super().__init__(config, r_interface)
        self.r_packages = [
            "diffHic",
            "InteractionSet",
            "GenomicRanges",
            "edgeR",
            "ggplot2",
            "dplyr",
        ]

    def run_analysis(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> DifferentialResult:
        """Run diffHic analysis"""

        logger.info("Starting diffHic analysis")

        try:
            # Check R packages
            if not self.r_interface.check_packages(self.r_packages):
                missing = self.r_interface.install_packages(self.r_packages)
                if missing:
                    raise RuntimeError(f"Failed to install R packages: {missing}")

            # Prepare R script
            r_script = self._create_diffhic_script(data, comparison)

            # Execute R analysis
            result = self.r_interface.run_script(r_script)

            if not result["success"]:
                return DifferentialResult(
                    comparison_name=comparison["name"],
                    method="diffHic",
                    success=False,
                    error_message=result.get("error", "Unknown R error"),
                )

            # Load results
            results_df = self._parse_diffhic_results(result)

            # Classify results
            results_df = self._classify_results(
                results_df,
                fdr_threshold=self.diff_params.get("fdr_threshold", 0.05),
                logfc_threshold=self.diff_params.get("logfc_threshold", 1.0),
            )

            # Calculate statistics
            stats = self._calculate_statistics(results_df)

            return DifferentialResult(
                comparison_name=comparison["name"],
                method="diffHic",
                success=True,
                results_table=results_df,
                **stats,
            )

        except Exception as e:
            logger.error(f"diffHic analysis failed: {e}")
            return DifferentialResult(
                comparison_name=comparison["name"],
                method="diffHic",
                success=False,
                error_message=str(e),
            )

    def _create_diffhic_script(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> str:
        """Create R script for diffHic analysis"""

        counts_file = data["counts_file"]
        metadata_file = data["metadata_file"]
        control_group = comparison["control"]
        treatment_group = comparison["treatment"]

        r_script = f"""
# Load required libraries
library(diffHic)
library(InteractionSet)
library(GenomicRanges)
library(edgeR)
library(dplyr)

# Set working directory
setwd("{data.get('output_dir', '.')}")

# Load data
counts_df <- read.csv("{counts_file}")
metadata <- read.csv("{metadata_file}")

# Parse loop coordinates to GRanges
parse_to_granges <- function(loop_ids) {{
    coords_split <- strsplit(loop_ids, "_")
    
    # Parse anchor 1
    anchor1_str <- sapply(coords_split, function(x) x[1])
    chr1 <- gsub("(chr[^:]+):.*", "\\\\1", anchor1_str)
    pos1_str <- gsub("chr[^:]+:(\\\\d+)-(\\\\d+)", "\\\\1-\\\\2", anchor1_str)
    start1 <- as.numeric(gsub("(\\\\d+)-\\\\d+", "\\\\1", pos1_str))
    end1 <- as.numeric(gsub("\\\\d+-(\\\\d+)", "\\\\1", pos1_str))
    
    # Parse anchor 2
    anchor2_str <- sapply(coords_split, function(x) x[2])
    chr2 <- gsub("(chr[^:]+):.*", "\\\\1", anchor2_str)
    pos2_str <- gsub("chr[^:]+:(\\\\d+)-(\\\\d+)", "\\\\1-\\\\2", anchor2_str)
    start2 <- as.numeric(gsub("(\\\\d+)-\\\\d+", "\\\\1", pos2_str))
    end2 <- as.numeric(gsub("\\\\d+-(\\\\d+)", "\\\\1", pos2_str))
    
    # Create GRanges
    anchor1 <- GRanges(seqnames = chr1, ranges = IRanges(start1, end1))
    anchor2 <- GRanges(seqnames = chr2, ranges = IRanges(start2, end2))
    
    return(list(anchor1 = anchor1, anchor2 = anchor2))
}}

# Prepare sample information
treatment_samples <- metadata$sample[metadata$cell_type == "{treatment_group}"]
control_samples <- metadata$sample[metadata$cell_type == "{control_group}"]
all_samples <- c(treatment_samples, control_samples)

# Prepare count matrix
count_matrix <- as.matrix(counts_df[, all_samples])
original_rownames <- counts_df$loop_id
rownames(count_matrix) <- NULL  # Remove for InteractionSet

# Parse coordinates
cat("Creating GInteractions objects...\\n")
loop_coords <- parse_to_granges(counts_df$loop_id)

interactions <- GInteractions(
    anchor1 = loop_coords$anchor1,
    anchor2 = loop_coords$anchor2
)

# Create InteractionSet
condition <- factor(c(rep("{treatment_group}", length(treatment_samples)), 
                     rep("{control_group}", length(control_samples))))

col_data <- DataFrame(
    condition = condition,
    sample_id = all_samples,
    totals = colSums(count_matrix),
    row.names = all_samples
)

hic_data <- InteractionSet(
    assays = list(counts = count_matrix),
    interactions = interactions,
    colData = col_data
)

# Store metadata
metadata(hic_data)$loop_ids <- original_rownames
metadata(hic_data)$source_info <- counts_df$source

# Normalization (simplified for moderate counts)
cat("Applying normalization...\\n")
distances <- pairdist(interactions)
valid_distances <- !is.na(distances) & distances > 0

if (any(valid_distances)) {{
    # Light distance correction for short-range data
    offsets <- matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix))
    log_dist <- log10(distances[valid_distances])
    dist_effect <- (log_dist - mean(log_dist)) * 0.05
    
    for (j in 1:ncol(count_matrix)) {{
        offsets[valid_distances, j] <- dist_effect
    }}
    
    assays(hic_data)$offset <- offsets
}}

# Statistical analysis with edgeR
cat("Running statistical analysis...\\n")
design <- model.matrix(~ 0 + condition, colData(hic_data))
colnames(design) <- levels(colData(hic_data)$condition)

# Convert to DGEList
y <- asDGEList(hic_data)

# Filter low counts
keep <- filterByExpr(y, design, min.count = {self.diff_params.get('min_counts', 1)})
y <- y[keep, ]
hic_data_filtered <- hic_data[keep, ]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")

# Dispersion estimation
y <- estimateDisp(y, design, robust = TRUE)

# GLM fitting
fit <- glmQLFit(y, design, robust = TRUE)

# Differential testing
contrast <- makeContrasts({treatment_group} - {control_group}, levels = design)
result <- glmQLFTest(fit, contrast = contrast)

# Extract results
results_table <- topTags(result, n = Inf, sort.by = "none")$table

# Add metadata back
filtered_positions <- which(keep)
if (length(filtered_positions) == nrow(results_table)) {{
    results_table$loop_id <- original_rownames[filtered_positions]
    if ("source" %in% colnames(counts_df)) {{
        results_table$source <- counts_df$source[filtered_positions]
    }}
}}

# Add genomic coordinates
coords_df <- as.data.frame(interactions(hic_data_filtered))
results_table <- cbind(results_table, coords_df)

# Calculate distances
if ("seqnames1" %in% colnames(results_table)) {{
    results_table$distance <- ifelse(
        results_table$seqnames1 == results_table$seqnames2,
        abs(results_table$start1 - results_table$start2),
        NA
    )
    results_table$interaction_type <- ifelse(
        results_table$seqnames1 == results_table$seqnames2,
        "cis", "trans"
    )
}}

# Save results
write.csv(results_table, "diffhic_results.csv", row.names = FALSE)

cat("diffHic analysis completed successfully\\n")
"""

        return r_script

    def _parse_diffhic_results(self, r_result: Dict[str, Any]) -> pd.DataFrame:
        """Parse diffHic results from R output"""

        # Look for results file
        results_file = Path(r_result.get("working_dir", ".")) / "diffhic_results.csv"

        if not results_file.exists():
            raise FileNotFoundError("diffHic results file not found")

        results_df = pd.read_csv(results_file)

        # Ensure required columns exist
        required_cols = ["logFC", "FDR", "PValue"]
        missing_cols = [col for col in required_cols if col not in results_df.columns]

        if missing_cols:
            raise ValueError(
                f"Missing required columns in diffHic results: {missing_cols}"
            )

        return results_df


class DESeq2Analyzer(BaseAnalyzer):
    """DESeq2-based differential analysis"""

    def __init__(self, config, r_interface: RInterface):
        super().__init__(config, r_interface)
        self.r_packages = ["DESeq2", "dplyr"]

    def run_analysis(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> DifferentialResult:
        """Run DESeq2 analysis"""

        logger.info("Starting DESeq2 analysis")

        try:
            # Check R packages
            if not self.r_interface.check_packages(self.r_packages):
                missing = self.r_interface.install_packages(self.r_packages)
                if missing:
                    raise RuntimeError(f"Failed to install R packages: {missing}")

            # Create R script
            r_script = self._create_deseq2_script(data, comparison)

            # Execute R analysis
            result = self.r_interface.run_script(r_script)

            if not result["success"]:
                return DifferentialResult(
                    comparison_name=comparison["name"],
                    method="DESeq2",
                    success=False,
                    error_message=result.get("error", "Unknown R error"),
                )

            # Parse results
            results_df = self._parse_deseq2_results(result)

            # Classify results
            results_df = self._classify_results(
                results_df,
                fdr_threshold=self.diff_params.get("fdr_threshold", 0.05),
                logfc_threshold=self.diff_params.get("logfc_threshold", 1.0),
            )

            # Calculate statistics
            stats = self._calculate_statistics(results_df)

            return DifferentialResult(
                comparison_name=comparison["name"],
                method="DESeq2",
                success=True,
                results_table=results_df,
                **stats,
            )

        except Exception as e:
            logger.error(f"DESeq2 analysis failed: {e}")
            return DifferentialResult(
                comparison_name=comparison["name"],
                method="DESeq2",
                success=False,
                error_message=str(e),
            )

    def _create_deseq2_script(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> str:
        """Create R script for DESeq2 analysis"""

        counts_file = data["counts_file"]
        metadata_file = data["metadata_file"]
        control_group = comparison["control"]
        treatment_group = comparison["treatment"]

        r_script = f"""
# Load required libraries
library(DESeq2)
library(dplyr)

# Set working directory
setwd("{data.get('output_dir', '.')}")

# Load data
counts_df <- read.csv("{counts_file}")
metadata <- read.csv("{metadata_file}")

# Prepare sample information
treatment_samples <- metadata$sample[metadata$cell_type == "{treatment_group}"]
control_samples <- metadata$sample[metadata$cell_type == "{control_group}"]
all_samples <- c(treatment_samples, control_samples)

# Prepare count matrix (integers for DESeq2)
count_matrix <- round(as.matrix(counts_df[, all_samples]))
rownames(count_matrix) <- counts_df$loop_id

# Create sample metadata
sample_data <- data.frame(
    sample = all_samples,
    condition = factor(c(rep("{treatment_group}", length(treatment_samples)), 
                        rep("{control_group}", length(control_samples)))),
    stringsAsFactors = FALSE
)
rownames(sample_data) <- all_samples

# Filter low counts
keep <- rowSums(count_matrix >= {self.diff_params.get('min_counts', 1)}) >= 2
count_matrix <- count_matrix[keep, ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_data,
    design = ~ condition
)

# Set reference level
dds$condition <- relevel(dds$condition, ref = "{control_group}")

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\\n")
dds <- DESeq(dds)

# Extract results
results_deseq2 <- results(dds, 
                         contrast = c("condition", "{treatment_group}", "{control_group}"),
                         alpha = {self.diff_params.get('fdr_threshold', 0.05)})

# Convert to data frame
results_df <- as.data.frame(results_deseq2)
results_df$loop_id <- rownames(results_df)

# Rename columns to match standard format
colnames(results_df)[colnames(results_df) == "log2FoldChange"] <- "logFC"
colnames(results_df)[colnames(results_df) == "padj"] <- "FDR" 
colnames(results_df)[colnames(results_df) == "pvalue"] <- "PValue"

# Remove rows with NA FDR
results_df <- results_df[!is.na(results_df$FDR), ]

# Add source information if available
if ("source" %in% colnames(counts_df)) {{
    loop_to_source <- setNames(counts_df$source, counts_df$loop_id)
    results_df$source <- loop_to_source[results_df$loop_id]
}}

# Save results
write.csv(results_df, "deseq2_results.csv", row.names = FALSE)

cat("DESeq2 analysis completed successfully\\n")
"""

        return r_script

    def _parse_deseq2_results(self, r_result: Dict[str, Any]) -> pd.DataFrame:
        """Parse DESeq2 results from R output"""

        results_file = Path(r_result.get("working_dir", ".")) / "deseq2_results.csv"

        if not results_file.exists():
            raise FileNotFoundError("DESeq2 results file not found")

        results_df = pd.read_csv(results_file)

        # Ensure required columns exist
        required_cols = ["logFC", "FDR", "PValue"]
        missing_cols = [col for col in required_cols if col not in results_df.columns]

        if missing_cols:
            raise ValueError(
                f"Missing required columns in DESeq2 results: {missing_cols}"
            )

        return results_df


class EdgeRAnalyzer(BaseAnalyzer):
    """edgeR-based differential analysis"""

    def __init__(self, config, r_interface: RInterface):
        super().__init__(config, r_interface)
        self.r_packages = ["edgeR", "dplyr"]

    def run_analysis(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> DifferentialResult:
        """Run edgeR analysis"""

        logger.info("Starting edgeR analysis")

        try:
            # Check R packages
            if not self.r_interface.check_packages(self.r_packages):
                missing = self.r_interface.install_packages(self.r_packages)
                if missing:
                    raise RuntimeError(f"Failed to install R packages: {missing}")

            # Create R script
            r_script = self._create_edger_script(data, comparison)

            # Execute R analysis
            result = self.r_interface.run_script(r_script)

            if not result["success"]:
                return DifferentialResult(
                    comparison_name=comparison["name"],
                    method="edgeR",
                    success=False,
                    error_message=result.get("error", "Unknown R error"),
                )

            # Parse results
            results_df = self._parse_edger_results(result)

            # Classify results
            results_df = self._classify_results(
                results_df,
                fdr_threshold=self.diff_params.get("fdr_threshold", 0.05),
                logfc_threshold=self.diff_params.get("logfc_threshold", 1.0),
            )

            # Calculate statistics
            stats = self._calculate_statistics(results_df)

            return DifferentialResult(
                comparison_name=comparison["name"],
                method="edgeR",
                success=True,
                results_table=results_df,
                **stats,
            )

        except Exception as e:
            logger.error(f"edgeR analysis failed: {e}")
            return DifferentialResult(
                comparison_name=comparison["name"],
                method="edgeR",
                success=False,
                error_message=str(e),
            )

    def _create_edger_script(
        self, data: Dict[str, Any], comparison: Dict[str, str]
    ) -> str:
        """Create R script for edgeR analysis"""

        counts_file = data["counts_file"]
        metadata_file = data["metadata_file"]
        control_group = comparison["control"]
        treatment_group = comparison["treatment"]

        r_script = f"""
# Load required libraries
library(edgeR)
library(dplyr)

# Set working directory
setwd("{data.get('output_dir', '.')}")

# Load data
counts_df <- read.csv("{counts_file}")
metadata <- read.csv("{metadata_file}")

# Prepare sample information
treatment_samples <- metadata$sample[metadata$cell_type == "{treatment_group}"]
control_samples <- metadata$sample[metadata$cell_type == "{control_group}"]
all_samples <- c(treatment_samples, control_samples)

# Prepare count matrix
count_matrix <- as.matrix(counts_df[, all_samples])
rownames(count_matrix) <- counts_df$loop_id

# Create sample groups
group <- factor(c(rep("{treatment_group}", length(treatment_samples)), 
                 rep("{control_group}", length(control_samples))))

# Create DGEList
y <- DGEList(counts = count_matrix, group = group)

# Filter low counts
keep <- filterByExpr(y, min.count = {self.diff_params.get('min_counts', 1)})
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y)

# Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Estimate dispersions
y <- estimateDisp(y, design, robust = TRUE)

# Quasi-likelihood fitting
fit <- glmQLFit(y, design, robust = TRUE)

# Make contrast
contrast <- makeContrasts({treatment_group} - {control_group}, levels = design)

# Quasi-likelihood test
qlf <- glmQLFTest(fit, contrast = contrast)

# Extract results
results_table <- topTags(qlf, n = Inf, sort.by = "none")$table
results_table$loop_id <- rownames(results_table)

# Rename columns to match standard format
if ("FDR" %in% colnames(results_table)) {{
    # Column already correctly named
}} else if ("adj.P.Val" %in% colnames(results_table)) {{
    colnames(results_table)[colnames(results_table) == "adj.P.Val"] <- "FDR"
}}

# Add source information if available
if ("source" %in% colnames(counts_df)) {{
    # Map results back to original data
    original_indices <- match(rownames(results_table), counts_df$loop_id)
    results_table$source <- counts_df$source[original_indices[!is.na(original_indices)]]
}}

# Remove rows with missing values
results_table <- results_table[!is.na(results_table$FDR), ]

# Save results
write.csv(results_table, "edger_results.csv", row.names = FALSE)

cat("edgeR analysis completed successfully\\n")
"""

        return r_script

    def _parse_edger_results(self, r_result: Dict[str, Any]) -> pd.DataFrame:
        """Parse edgeR results from R output"""

        results_file = Path(r_result.get("working_dir", ".")) / "edger_results.csv"

        if not results_file.exists():
            raise FileNotFoundError("edgeR results file not found")

        results_df = pd.read_csv(results_file)

        # Ensure required columns exist
        required_cols = ["logFC", "FDR", "PValue"]
        missing_cols = [col for col in required_cols if col not in results_df.columns]

        if missing_cols:
            raise ValueError(
                f"Missing required columns in edgeR results: {missing_cols}"
            )

        return results_df
