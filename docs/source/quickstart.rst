Quick Start Guide
=================

This guide will help you get started with PlasmaFlow for Hi-C chromatin loop analysis.

Basic Workflow
--------------

PlasmaFlow follows a sequential analysis workflow:

1. **Loop Calling**: Identify chromatin loops from Hi-C data using Peakachu
2. **Quality Control**: Validate loops using Aggregate Peak Analysis (APA)
3. **Differential Analysis**: Compare loop changes between conditions
4. **Visualization**: Generate publication-ready plots and heatmaps
5. **Gene Analysis**: Map genes near differential loops
6. **Pathway Enrichment**: Analyze biological pathways

Setup
-----

First, create a configuration file:

.. code-block:: bash

    plasmaflow init-config

This creates ``plasmaflow_config.yaml``. Edit this file to specify your data paths and analysis parameters.

Running the Complete Workflow
-----------------------------

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~

Run the complete analysis pipeline:

.. code-block:: bash

    # Run complete workflow
    plasmaflow run-workflow --config plasmaflow_config.yaml
    
    # Or run individual steps
    plasmaflow loop-calling --config plasmaflow_config.yaml
    plasmaflow quality-control --config plasmaflow_config.yaml
    plasmaflow differential-analysis --config plasmaflow_config.yaml

Python API
~~~~~~~~~~

.. code-block:: python

    from plasmaflow import PlasmaFlowAnalysis
    
    # Initialize analysis
    analyzer = PlasmaFlowAnalysis(config_file="plasmaflow_config.yaml")
    
    # Run complete workflow
    results = analyzer.run_complete_workflow()
    
    # Or run individual steps
    loop_results = analyzer.run_loop_calling()
    qc_results = analyzer.run_quality_control()
    diff_results = analyzer.run_differential_analysis()

Step-by-Step Analysis
--------------------

1. Loop Calling
~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.loop_calling import PeakachuAnalyzer
    from plasmaflow.core.config import Config
    
    # Load configuration
    config = Config.from_yaml("plasmaflow_config.yaml")
    
    # Initialize loop caller
    analyzer = PeakachuAnalyzer(config)
    
    # Run loop calling for all samples
    results = analyzer.run_batch_analysis()

2. Quality Control
~~~~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.quality_control import APAAnalyzer
    
    # Initialize APA analyzer
    apa_analyzer = APAAnalyzer(config)
    
    # Run APA analysis
    apa_results = apa_analyzer.run_apa_analysis(
        loops_file="output/sample1_loops.bedpe",
        cool_file="data/sample1.cool"
    )
    
    # Generate QC plots
    apa_analyzer.plot_apa_results(apa_results)

3. Differential Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.differential import DifferentialAnalyzer
    
    # Initialize differential analyzer
    diff_analyzer = DifferentialAnalyzer(config)
    
    # Compare conditions
    diff_results = diff_analyzer.run_differential_analysis(
        sample_info="samples.csv",
        comparison="condition_A_vs_condition_B"
    )
    
    # Get significant loops
    significant_loops = diff_results.get_significant_loops(
        fdr_threshold=0.05,
        logfc_threshold=1.0
    )

4. Visualization
~~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.visualization import HeatmapPlotter
    
    # Create heatmaps
    plotter = HeatmapPlotter(config)
    
    # Generate differential loop heatmap
    heatmap_path = plotter.create_differential_heatmap(
        differential_results=diff_results,
        output_file="differential_loops_heatmap.pdf"
    )

5. Gene Proximity Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.genomics import GeneProximityAnalyzer
    
    # Analyze genes near differential loops
    gene_analyzer = GeneProximityAnalyzer(config)
    
    proximity_results = gene_analyzer.analyze_gene_proximity(
        loops_file="significant_loops.bedpe",
        distance_threshold=100000  # 100kb
    )

6. Pathway Enrichment
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from plasmaflow.enrichment import ClusterProfilerAnalyzer
    
    # Run pathway enrichment
    enrichment_analyzer = ClusterProfilerAnalyzer(config)
    
    # Get differentially expressed genes near loops
    de_genes = proximity_results.get_de_genes()
    
    # Run comprehensive enrichment analysis
    enrichment_results = enrichment_analyzer.run_comprehensive_analysis(
        gene_data=de_genes,
        include_msigdb=True
    )

Configuration Example
--------------------

Here's an example configuration file:

.. code-block:: yaml

    # PlasmaFlow Configuration
    analysis:
      output_dir: "/path/to/output"
      temp_dir: "/tmp/plasmaflow"
      n_cores: 8
      
    samples:
      - name: "memB_rep1"
        cool_file: "/data/memB_rep1.cool"
        condition: "memB"
      - name: "prePB_rep1" 
        cool_file: "/data/prePB_rep1.cool"
        condition: "prePB"
        
    loop_calling:
      peakachu_model: "GM12878"
      resolution: 10000
      balance: true
      
    differential:
      methods: ["diffHic", "DESeq2"]
      comparisons:
        - name: "prePB_vs_memB"
          condition1: "prePB"
          condition2: "memB"
          
    visualization:
      flanking_bp: 100000
      colormap: "RdBu_r"
      
    r_config:
      r_executable: "/usr/bin/R"
      library_path: "/usr/local/lib/R/site-library"

Data Requirements
----------------

Input Files
~~~~~~~~~~

* **Hi-C Contact Matrices**: `.cool` format files (10kb resolution recommended)
* **Sample Information**: CSV file with sample metadata
* **Genome Annotation**: GTF file for gene proximity analysis (optional)
* **Gene Expression**: Count matrix for pathway analysis (optional)

File Organization
~~~~~~~~~~~~~~~~

Recommended directory structure:

.. code-block:: text

    project/
    ├── data/
    │   ├── sample1.cool
    │   ├── sample2.cool
    │   └── samples.csv
    ├── config/
    │   └── plasmaflow_config.yaml
    ├── output/
    │   ├── loop_calling/
    │   ├── quality_control/
    │   ├── differential/
    │   └── plots/
    └── scripts/
        └── run_analysis.py

Output Files
-----------

PlasmaFlow generates the following outputs:

* **Loop calling**: `.bedpe` files with called loops
* **Quality control**: APA plots and QC metrics
* **Differential analysis**: CSV files with statistical results
* **Visualization**: PDF plots and heatmaps
* **Gene analysis**: Gene lists and proximity statistics
* **Pathway enrichment**: Enrichment tables and network plots

Next Steps
----------

* Check the :doc:`tutorials/index` for detailed examples
* See :doc:`api/index` for complete API documentation
* Visit :doc:`examples/index` for real-world analysis examples