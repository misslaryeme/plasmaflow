Changelog
=========

All notable changes to PlasmaFlow will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[Unreleased]
------------

[1.0.0] - 2024-XX-XX
--------------------

Initial release of PlasmaFlow - A comprehensive Python package for Hi-C chromatin loop analysis.

Added
~~~~~

**Core Infrastructure:**

- Complete package structure with modular design
- YAML-based configuration system with validation
- Comprehensive error handling and logging
- Command-line interface using Click
- Automated testing framework with pytest
- Sphinx documentation with RTD theme

**Loop Calling Module:**

- Integration with Peakachu framework for loop identification
- Support for multiple Hi-C resolutions and file formats
- Batch processing capabilities for multiple samples
- Model management system for different cell types
- Comprehensive result validation and QC metrics

**Quality Control Module:**

- Aggregate Peak Analysis (APA) implementation using coolpuppy
- Comparative APA analysis across conditions
- Loop validation with statistical significance testing
- Publication-ready APA visualizations
- Quality metrics calculation and reporting

**Differential Analysis Module:**

- Integration with R/Bioconductor statistical methods
- Support for diffHic, DESeq2, and edgeR approaches
- Multi-method comparative analysis
- FDR correction and effect size filtering
- Comprehensive result export and visualization

**Visualization Module:**

- Hi-C contact matrix heatmap generation
- Differential loop intensity visualization
- APA plot creation with customizable styling
- Matrix generation using deepTools integration
- Multi-panel figure composition for publications

**Genomics Module:**

- TSS-based gene proximity calculations
- Regulatory region classification (promoter, enhancer, gene body)
- Differential gene expression integration
- Enhanced volcano plots with loop proximity annotation
- Statistical testing for gene expression enrichment near loops

**Enrichment Module:**

- R-based pathway enrichment using clusterProfiler
- MSigDB C2 gene set integration with B-cell focus
- Support for GO, KEGG, and Reactome databases
- Both ORA (Over-Representation Analysis) and GSEA methods
- Custom B-cell differentiation gene set curation
- Publication-ready pathway visualization plots

**Workflow Integration:**

- Complete end-to-end analysis orchestration
- Configurable pipeline execution
- Result caching and intermediate file management
- Progress tracking and error recovery
- Multi-sample batch processing

**Development Tools:**

- Pre-commit hooks for code quality
- GitHub Actions CI/CD pipeline
- Docker containerization support
- Comprehensive test coverage (>90%)
- Type hints throughout codebase

Technical Details
~~~~~~~~~~~~~~~~

**Supported File Formats:**

- Hi-C: .cool, .mcool formats via cooler
- Loops: .bedpe format 
- Annotations: GTF/GFF3 gene annotation files
- Expression: CSV/TSV count matrices
- Output: PDF, PNG, SVG visualizations

**External Dependencies:**

- Peakachu: Loop calling framework
- deepTools: Matrix computation and visualization
- R/Bioconductor: Statistical analysis packages
- cooler/cooltools: Hi-C data manipulation

**Statistical Methods:**

- diffHic: Hi-C specific differential analysis
- DESeq2: RNA-seq method adapted for Hi-C
- edgeR: Count-based statistical testing
- clusterProfiler: Comprehensive pathway analysis

**Visualization Capabilities:**

- Contact matrix heatmaps with loop annotations
- APA plots with statistical significance
- Violin plots for comparative analysis
- Network plots for pathway relationships
- Volcano plots for expression analysis

Documentation
~~~~~~~~~~~~

- Complete API reference with examples
- Step-by-step tutorials for common workflows
- Configuration guide with templates
- Installation instructions for all platforms
- Contributing guidelines for developers

**Analysis Examples:**

- B-cell to plasma cell differentiation
- Comparative Hi-C analysis workflows
- Integration with RNA-seq data
- Pathway enrichment interpretation

Performance
~~~~~~~~~~

- Multi-core processing support
- Memory-efficient large dataset handling
- Result caching for improved performance
- Progress tracking for long-running analyses

Security
~~~~~~~

- Input validation for all file operations
- Safe R script execution with subprocess isolation
- Temporary file cleanup
- Error handling prevents data corruption