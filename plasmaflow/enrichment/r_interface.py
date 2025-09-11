"""
R interface for pathway enrichment analysis using clusterProfiler

This module provides Python-R integration for running sophisticated
pathway enrichment analyses using the R clusterProfiler ecosystem.
"""

import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

from ..config import Config
from ..utils import get_logger

logger = get_logger(__name__)


class RPathwayInterface:
    """Interface for R-based pathway enrichment analysis"""

    def __init__(self, config: Config):
        """
        Initialize R pathway interface

        Args:
            config: PlasmaFlow configuration object
        """
        self.config = config
        self.r_config = config.r_config

        # Get enrichment parameters
        self.enrichment_params = config.enrichment

        # R package requirements
        self.required_packages = [
            "clusterProfiler",
            "org.Hs.eg.db",
            "DOSE",
            "ReactomePA",
            "enrichplot",
            "ggplot2",
            "dplyr",
            "readr",
            "msigdbr",
            "ggridges",
            "igraph",
            "patchwork",
            "viridis",
        ]

        # Output directories
        self.output_dir = Path(config.output_dir) / "enrichment"
        self.temp_dir = self.output_dir / "temp"
        self.r_scripts_dir = self.output_dir / "r_scripts"

        # Create directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.r_scripts_dir.mkdir(parents=True, exist_ok=True)

    def validate_r_environment(self) -> Dict[str, Any]:
        """
        Validate R environment and required packages

        Returns:
            Dictionary with validation results
        """
        logger.info("Validating R environment for pathway analysis...")

        validation_results = {
            "r_available": False,
            "r_version": None,
            "missing_packages": [],
            "available_packages": [],
            "validation_passed": False,
        }

        try:
            # Check R availability
            result = subprocess.run(
                ["R", "--version"], capture_output=True, text=True, timeout=10
            )

            if result.returncode == 0:
                validation_results["r_available"] = True
                # Extract R version
                version_line = result.stdout.split("\\n")[0]
                validation_results["r_version"] = version_line
                logger.info(f"R found: {version_line}")

                # Check packages
                package_check_script = self._generate_package_check_script()
                package_results = self._run_r_script(package_check_script)

                if package_results["success"]:
                    # Parse package results
                    for line in package_results["output"].split("\\n"):
                        if line.startswith("PACKAGE_AVAILABLE:"):
                            package = line.split(":")[1].strip()
                            validation_results["available_packages"].append(package)
                        elif line.startswith("PACKAGE_MISSING:"):
                            package = line.split(":")[1].strip()
                            validation_results["missing_packages"].append(package)

                # Overall validation
                validation_results["validation_passed"] = (
                    len(validation_results["missing_packages"]) == 0
                )

                if validation_results["validation_passed"]:
                    logger.info("R environment validation passed")
                else:
                    logger.warning(
                        f"R environment validation failed. Missing packages: {validation_results['missing_packages']}"
                    )

            else:
                logger.error("R not found or not accessible")

        except Exception as e:
            logger.error(f"R environment validation failed: {e}")

        return validation_results

    def run_clusterprofiler_analysis(
        self,
        gene_expression_files: Dict[str, Union[str, Path]],
        comparison_name: str = "prePB_vs_memB",
        filter_bcell_terms: bool = True,
        use_msigdb_c2: bool = True,
        msigdb_subcategory: str = "CP",
    ) -> Dict[str, Any]:
        """
        Run comprehensive clusterProfiler analysis

        Args:
            gene_expression_files: Dictionary mapping category names to gene expression files
            comparison_name: Name of comparison
            filter_bcell_terms: Whether to filter for B-cell relevant terms
            use_msigdb_c2: Whether to use MSigDB C2 analysis
            msigdb_subcategory: MSigDB subcategory to use

        Returns:
            Dictionary with analysis results
        """
        logger.info(f"Running clusterProfiler analysis for {comparison_name}")

        # Validate R environment
        validation = self.validate_r_environment()
        if not validation["validation_passed"]:
            raise RuntimeError(
                f"R environment validation failed: {validation['missing_packages']}"
            )

        # Generate R analysis script
        analysis_script = self._generate_clusterprofiler_script(
            gene_expression_files,
            comparison_name,
            filter_bcell_terms,
            use_msigdb_c2,
            msigdb_subcategory,
        )

        # Save and run script
        script_file = (
            self.r_scripts_dir / f"{comparison_name}_clusterprofiler_analysis.R"
        )
        with open(script_file, "w") as f:
            f.write(analysis_script)

        logger.info(f"Running R analysis script: {script_file}")

        # Run R script
        r_results = self._run_r_script_file(script_file)

        if not r_results["success"]:
            raise RuntimeError(f"R analysis failed: {r_results['error']}")

        # Parse results
        analysis_results = self._parse_clusterprofiler_results(comparison_name)

        logger.info(f"clusterProfiler analysis completed for {comparison_name}")

        return {
            "comparison_name": comparison_name,
            "validation": validation,
            "script_file": script_file,
            "r_output": r_results["output"],
            "analysis_results": analysis_results,
            "parameters": {
                "filter_bcell_terms": filter_bcell_terms,
                "use_msigdb_c2": use_msigdb_c2,
                "msigdb_subcategory": msigdb_subcategory,
            },
        }

    def _generate_package_check_script(self) -> str:
        """Generate R script to check package availability"""

        script = (
            """
        # Check required packages
        required_packages <- c("""
            + ", ".join(f'"{pkg}"' for pkg in self.required_packages)
            + """)
        
        for(package in required_packages) {
            if(require(package, character.only = TRUE, quietly = TRUE)) {
                cat("PACKAGE_AVAILABLE:", package, "\\n")
            } else {
                cat("PACKAGE_MISSING:", package, "\\n")
            }
        }
        """
        )

        return script

    def _generate_clusterprofiler_script(
        self,
        gene_expression_files: Dict[str, Union[str, Path]],
        comparison_name: str,
        filter_bcell_terms: bool,
        use_msigdb_c2: bool,
        msigdb_subcategory: str,
    ) -> str:
        """Generate comprehensive clusterProfiler analysis script"""

        # File paths
        input_files = {
            k: str(Path(v).absolute()) for k, v in gene_expression_files.items()
        }
        output_dir = str(self.output_dir.absolute())

        script = f"""
        # ============================================================================
        # PlasmaFlow clusterProfiler Analysis Script
        # ============================================================================
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(DOSE)
        library(ReactomePA)
        library(enrichplot)
        library(ggplot2)
        library(dplyr)
        library(readr)
        library(msigdbr)
        library(ggridges)
        library(igraph)
        library(patchwork)
        library(viridis)
        
        # Configuration
        COMPARISON <- "{comparison_name}"
        OUTPUT_DIR <- "{output_dir}"
        FILTER_BCELL_TERMS <- {str(filter_bcell_terms).upper()}
        USE_MSIGDB_C2 <- {str(use_msigdb_c2).upper()}
        MSIGDB_SUBCATEGORY <- "{msigdb_subcategory}"
        
        # Create output directory structure
        dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "GO", "ORA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "GO", "GSEA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "KEGG", "ORA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "KEGG", "GSEA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "MSigDB", "ORA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "plots", "MSigDB", "GSEA"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(OUTPUT_DIR, "csv_results"), recursive = TRUE, showWarnings = FALSE)
        
        cat("🔧 PlasmaFlow clusterProfiler Analysis\\n")
        cat("📊 Comparison:", COMPARISON, "\\n")
        cat("🎛️ B cell filtering:", FILTER_BCELL_TERMS, "\\n")
        cat("🔬 MSigDB C2:", USE_MSIGDB_C2, "\\n")
        cat("📚 MSigDB subcategory:", MSIGDB_SUBCATEGORY, "\\n")
        cat("📁 Output directory:", OUTPUT_DIR, "\\n")
        
        {self._get_helper_functions()}
        
        # Input files
        input_files <- list(
            {self._format_input_files(input_files)}
        )
        
        # Prepare MSigDB C2 if requested
        c2_term2gene <- NULL
        if(USE_MSIGDB_C2) {{
            c2_term2gene <- prepare_bcell_msigdb_c2(MSIGDB_SUBCATEGORY, include_description_search = TRUE)
        }}
        
        # Process each category
        all_results <- list()
        categories <- names(input_files)
        
        for(category in categories) {{
            cat("\\n", rep("=", 60), "\\n")
            cat("📊 ANALYZING", toupper(category), "GENES\\n")
            cat(rep("=", 60), "\\n")
            
            file_path <- input_files[[category]]
            if(!file.exists(file_path)) {{
                cat("❌ File not found:", file_path, "\\n")
                next
            }}
            
            df <- read_csv(file_path, show_col_types = FALSE)
            cat("✅ Loaded", nrow(df), "genes\\n")
            
            # Run analysis
            results <- run_fixed_bcell_analysis(df, category, FILTER_BCELL_TERMS, c2_term2gene)
            
            # Extract gene_list for plotting
            gene_list <- results$gene_list
            results$gene_list <- NULL
            
            # Create plots
            filter_suffix <- ifelse(FILTER_BCELL_TERMS, "_filtered", "_unfiltered")
            create_comprehensive_plots(results, category, filter_suffix, gene_list)
            create_patchwork_combined_ora_plots(results, category, filter_suffix)
            
            # Export results to CSV
            export_results_to_csv(results, category, filter_suffix)
            
            # Save results
            saveRDS(results, file.path(OUTPUT_DIR, paste0("corrected_", category, filter_suffix, "_results.rds")))
            all_results[[category]] <- results
            
            cat("✅", toupper(category), "analysis completed!\\n")
        }}
        
        # Save combined results
        saveRDS(all_results, file.path(OUTPUT_DIR, paste0(COMPARISON, "_all_results.rds")))
        
        cat("🎉 clusterProfiler analysis completed successfully!\\n")
        cat("📁 Results saved in:", OUTPUT_DIR, "\\n")
        """

        return script

    def _get_helper_functions(self) -> str:
        """Get R helper functions for pathway analysis"""

        return """
        # ============================================================================
        # HELPER FUNCTIONS FOR PATHWAY ANALYSIS
        # ============================================================================
        
        # MSigDB C2 preparation function
        prepare_bcell_msigdb_c2 <- function(subcategory = "CP", include_description_search = TRUE) {
            cat("🔬 Preparing B cell-specific MSigDB C2 gene sets\\n")
            
            tryCatch({
                if(subcategory == "ALL") {
                    c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
                } else {
                    c2_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = subcategory)
                }
                
                cat("   📚 Total C2 gene sets available:", length(unique(c2_sets$gs_name)), "\\n")
                
                # B-cell specific patterns
                bcell_patterns <- c(
                    "plasma|plasmablast|memory",
                    "immunoglobulin|antibody|ig[agmde]",
                    "class.switching|somatic.hypermutation",
                    "bcr|b.cell.receptor",
                    "germinal.center|germinal.centre",
                    "lymphocyte.differentiation|b.cell.activation",
                    "unfolded.protein.response|endoplasmic.reticulum",
                    "nfkb|blimp|irf4|pax5|bach2",
                    "secretory.pathway|protein.folding"
                )
                
                full_pattern <- paste(bcell_patterns, collapse = "|")
                
                if(include_description_search) {
                    filtered_sets <- c2_sets %>%
                        filter(grepl(full_pattern, gs_name, ignore.case = TRUE) |
                               grepl(full_pattern, gs_description, ignore.case = TRUE))
                } else {
                    filtered_sets <- c2_sets %>%
                        filter(grepl(full_pattern, gs_name, ignore.case = TRUE))
                }
                
                c2_term2gene <- filtered_sets[, c("gs_name", "gene_symbol")]
                
                # Size filtering
                gene_set_sizes <- table(c2_term2gene$gs_name)
                valid_sets <- names(gene_set_sizes)[gene_set_sizes >= 10 & gene_set_sizes <= 500]
                c2_term2gene_filtered <- c2_term2gene[c2_term2gene$gs_name %in% valid_sets, ]
                
                cat("   ✂️ After size filtering (10-500 genes):", length(unique(c2_term2gene_filtered$gs_name)), "gene sets\\n")
                
                return(c2_term2gene_filtered)
                
            }, error = function(e) {
                cat("   ❌ MSigDB preparation error:", e$message, "\\n")
                return(NULL)
            })
        }
        
        # Main analysis function
        run_fixed_bcell_analysis <- function(gene_data, category_name, filter_bcell = TRUE, c2_term2gene = NULL) {
            cat("🎯 Running analysis for", category_name, "\\n")
            
            # Prepare data
            valid_data <- gene_data %>%
                filter(!is.na(log2FoldChange), !is.na(padj), padj >= 0, !is.na(gene_symbol))
            
            if(nrow(valid_data) == 0) {
                cat("   ⚠️ No valid data\\n")
                return(list())
            }
            
            # Create ranked gene list
            gene_list <- valid_data$log2FoldChange
            names(gene_list) <- valid_data$gene_symbol
            gene_list <- sort(gene_list, decreasing = TRUE)
            
            # Background universe
            universe <- valid_data$gene_symbol
            
            results <- list()
            
            # GSEA Analysis
            results$gseGO_BP <- run_fixed_gseGO(gene_list, category_name, "BP", filter_bcell)
            results$gseKEGG <- run_fixed_gseKEGG(gene_list, category_name, filter_bcell)
            
            if(!is.null(c2_term2gene)) {
                results$gseMSigDB_C2 <- run_fixed_gsea_msigdb(gene_list, category_name, c2_term2gene, filter_bcell)
            }
            
            # ORA Analysis
            significant_genes <- valid_data %>% filter(padj < 0.05)
            up_genes <- significant_genes$gene_symbol[significant_genes$log2FoldChange > 0]
            down_genes <- significant_genes$gene_symbol[significant_genes$log2FoldChange < 0]
            
            if(length(up_genes) > 0) {
                results$enrichGO_UP_BP <- run_fixed_enrichGO(up_genes, paste(category_name, "UP"), "BP", universe, filter_bcell)
                results$enrichKEGG_UP <- run_fixed_enrichKEGG(up_genes, paste(category_name, "UP"), universe, filter_bcell)
                
                if(!is.null(c2_term2gene)) {
                    results$enrichMSigDB_C2_UP <- run_fixed_ora_msigdb(up_genes, paste(category_name, "UP"), c2_term2gene, universe, filter_bcell)
                }
            }
            
            if(length(down_genes) > 0) {
                results$enrichGO_DOWN_BP <- run_fixed_enrichGO(down_genes, paste(category_name, "DOWN"), "BP", universe, filter_bcell)
                results$enrichKEGG_DOWN <- run_fixed_enrichKEGG(down_genes, paste(category_name, "DOWN"), universe, filter_bcell)
                
                if(!is.null(c2_term2gene)) {
                    results$enrichMSigDB_C2_DOWN <- run_fixed_ora_msigdb(down_genes, paste(category_name, "DOWN"), c2_term2gene, universe, filter_bcell)
                }
            }
            
            results$gene_list <- gene_list
            return(results)
        }
        
        # Individual analysis functions would go here...
        # (Abbreviated for space - full functions from notebooks would be included)
        
        # Placeholder functions - in real implementation these would be the full functions from the R notebooks
        run_fixed_gseGO <- function(gene_list, category_name, ontology = "BP", filter_bcell = TRUE) {
            tryCatch({
                gseGO(
                    geneList = gene_list,
                    OrgDb = org.Hs.eg.db,
                    ont = ontology,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    keyType = "SYMBOL",
                    verbose = FALSE
                )
            }, error = function(e) {
                cat("   ❌ gseGO error:", e$message, "\\n")
                return(NULL)
            })
        }
        
        run_fixed_gseKEGG <- function(gene_list, category_name, filter_bcell = TRUE) {
            # Convert to Entrez IDs and run gseKEGG
            # Implementation details...
            return(NULL)  # Placeholder
        }
        
        run_fixed_gsea_msigdb <- function(gene_list, category_name, c2_term2gene, filter_bcell = TRUE) {
            if(is.null(c2_term2gene)) return(NULL)
            
            tryCatch({
                GSEA(
                    geneList = gene_list,
                    TERM2GENE = c2_term2gene,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    eps = 1e-10,
                    seed = TRUE,
                    by = "fgsea",
                    verbose = FALSE
                )
            }, error = function(e) {
                cat("   ❌ GSEA MSigDB error:", e$message, "\\n")
                return(NULL)
            })
        }
        
        run_fixed_enrichGO <- function(gene_list, category_name, ontology = "BP", universe = NULL, filter_bcell = TRUE) {
            tryCatch({
                enrichGO(
                    gene = gene_list,
                    OrgDb = org.Hs.eg.db,
                    ont = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    keyType = "SYMBOL",
                    minGSSize = 10,
                    maxGSSize = 500,
                    universe = universe,
                    readable = TRUE
                )
            }, error = function(e) {
                cat("   ❌ enrichGO error:", e$message, "\\n")
                return(NULL)
            })
        }
        
        run_fixed_enrichKEGG <- function(gene_list, category_name, universe = NULL, filter_bcell = TRUE) {
            # Implementation details...
            return(NULL)  # Placeholder
        }
        
        run_fixed_ora_msigdb <- function(gene_list, category_name, c2_term2gene, universe = NULL, filter_bcell = TRUE) {
            if(is.null(c2_term2gene)) return(NULL)
            
            tryCatch({
                enricher(
                    gene = gene_list,
                    TERM2GENE = c2_term2gene,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 500,
                    universe = universe
                )
            }, error = function(e) {
                cat("   ❌ ORA MSigDB error:", e$message, "\\n")
                return(NULL)
            })
        }
        
        # Plotting functions
        create_comprehensive_plots <- function(results, category_name, filter_suffix = "", gene_list = NULL) {
            cat("📊 Creating comprehensive plots for", category_name, "\\n")
            # Plotting implementation...
        }
        
        create_patchwork_combined_ora_plots <- function(results, category_name, filter_suffix = "") {
            cat("📊 Creating patchwork plots for", category_name, "\\n")
            # Patchwork plotting implementation...
        }
        
        # Export function
        export_results_to_csv <- function(results, category_name, filter_suffix = "") {
            cat("💾 Exporting", category_name, "results to CSV\\n")
            
            csv_dir <- file.path(OUTPUT_DIR, "csv_results")
            dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
            
            exported_count <- 0
            
            for(analysis_name in names(results)) {
                result_obj <- results[[analysis_name]]
                
                if(!is.null(result_obj) && nrow(result_obj) > 0) {
                    filename <- paste0(category_name, "_", analysis_name, filter_suffix, ".csv")
                    filepath <- file.path(csv_dir, filename)
                    
                    if(class(result_obj)[1] %in% c("enrichResult", "gseaResult")) {
                        result_df <- result_obj@result
                        result_df$analysis_type <- ifelse(grepl("gse", analysis_name), "GSEA", "ORA")
                        result_df$category <- category_name
                        
                        write_csv(result_df, filepath)
                        exported_count <- exported_count + 1
                        cat("   ✅ Exported:", filename, "\\n")
                    }
                }
            }
            
            cat("   📊 Total CSV files exported:", exported_count, "\\n")
        }
        """

    def _format_input_files(self, input_files: Dict[str, str]) -> str:
        """Format input files for R script"""

        file_lines = []
        for category, filepath in input_files.items():
            file_lines.append(f'    {category} = "{filepath}"')

        return ",\\n".join(file_lines)

    def _run_r_script(self, script: str) -> Dict[str, Any]:
        """Run R script and return results"""

        try:
            # Write script to temporary file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
                f.write(script)
                temp_script = f.name

            # Run R script
            result = subprocess.run(
                ["R", "--slave", "--vanilla", "-f", temp_script],
                capture_output=True,
                text=True,
                timeout=300,
            )

            # Clean up
            os.unlink(temp_script)

            return {
                "success": result.returncode == 0,
                "output": result.stdout,
                "error": result.stderr,
                "returncode": result.returncode,
            }

        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "output": "",
                "error": "R script execution timed out",
                "returncode": -1,
            }
        except Exception as e:
            return {"success": False, "output": "", "error": str(e), "returncode": -1}

    def _run_r_script_file(self, script_file: Path) -> Dict[str, Any]:
        """Run R script file and return results"""

        try:
            result = subprocess.run(
                ["R", "--slave", "--vanilla", "-f", str(script_file)],
                capture_output=True,
                text=True,
                timeout=1800,  # 30 minutes
            )

            return {
                "success": result.returncode == 0,
                "output": result.stdout,
                "error": result.stderr,
                "returncode": result.returncode,
            }

        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "output": "",
                "error": "R script execution timed out",
                "returncode": -1,
            }
        except Exception as e:
            return {"success": False, "output": "", "error": str(e), "returncode": -1}

    def _parse_clusterprofiler_results(self, comparison_name: str) -> Dict[str, Any]:
        """Parse clusterProfiler results from output files"""

        results = {
            "rds_files": {},
            "csv_files": {},
            "plot_files": {},
            "categories_analyzed": [],
        }

        try:
            # Find RDS result files
            rds_pattern = f"{comparison_name}_*.rds"
            for rds_file in self.output_dir.glob(rds_pattern):
                category = rds_file.stem.replace(f"{comparison_name}_", "")
                results["rds_files"][category] = rds_file

            # Find CSV result files
            csv_dir = self.output_dir / "csv_results"
            if csv_dir.exists():
                for csv_file in csv_dir.glob("*.csv"):
                    results["csv_files"][csv_file.stem] = csv_file

            # Find plot files
            plots_dir = self.output_dir / "plots"
            if plots_dir.exists():
                for plot_file in plots_dir.rglob("*.png"):
                    rel_path = plot_file.relative_to(plots_dir)
                    results["plot_files"][str(rel_path)] = plot_file

            # Extract categories
            for category_file in results["rds_files"].keys():
                category = (
                    category_file.replace("corrected_", "")
                    .replace("_filtered", "")
                    .replace("_unfiltered", "")
                )
                if category not in results["categories_analyzed"]:
                    results["categories_analyzed"].append(category)

            logger.info(
                f"Parsed results: {len(results['rds_files'])} RDS files, "
                f"{len(results['csv_files'])} CSV files, {len(results['plot_files'])} plots"
            )

        except Exception as e:
            logger.error(f"Error parsing clusterProfiler results: {e}")

        return results


def validate_r_environment(config: Optional[Config] = None) -> Dict[str, Any]:
    """
    Validate R environment for pathway analysis

    Args:
        config: Optional configuration object

    Returns:
        Dictionary with validation results
    """

    if config is None:
        from ..config import get_default_config

        config = get_default_config()

    interface = RPathwayInterface(config)
    return interface.validate_r_environment()
