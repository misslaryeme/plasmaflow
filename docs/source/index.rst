PlasmaFlow Documentation
========================

Welcome to PlasmaFlow
---------------------

PlasmaFlow is a comprehensive Python package for Hi-C chromatin loop analysis, specifically designed for studying plasma cell differentiation. This package provides a complete workflow from raw Hi-C data processing to advanced pathway enrichment analysis.

Key Features
-----------

* **Complete Hi-C Analysis Pipeline**: From loop calling to differential analysis
* **Quality Control**: Comprehensive APA (Aggregate Peak Analysis) validation
* **Statistical Analysis**: Integration with R/Bioconductor packages (diffHic, DESeq2, edgeR)
* **Visualization**: Publication-ready plots and heatmaps
* **Genomics Integration**: Gene proximity and expression analysis
* **Pathway Enrichment**: Advanced analysis with clusterProfiler and MSigDB
* **Command Line Interface**: Easy-to-use CLI for complete workflows
* **Extensible Architecture**: Modular design for custom analysis

Quick Start
-----------

Installation::

    pip install plasmaflow

Basic usage::

    from plasmaflow import PlasmaFlowAnalysis
    
    # Initialize analysis
    analyzer = PlasmaFlowAnalysis(config_file="config.yaml")
    
    # Run complete workflow
    results = analyzer.run_complete_workflow()

For detailed examples, see the :doc:`tutorials/index` section.

Package Overview
----------------

PlasmaFlow is organized into several key modules:

* :doc:`api/loop_calling`: Peakachu-based loop calling functionality
* :doc:`api/quality_control`: APA analysis and quality assessment
* :doc:`api/differential`: Statistical differential loop analysis
* :doc:`api/visualization`: Comprehensive plotting and visualization
* :doc:`api/genomics`: Gene proximity and expression analysis  
* :doc:`api/enrichment`: Pathway enrichment with R integration

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   tutorials/index
   api/index
   configuration
   examples/index
   contributing
   changelog

