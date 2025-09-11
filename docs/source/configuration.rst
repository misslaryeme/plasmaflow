Configuration
=============

PlasmaFlow uses YAML configuration files to manage analysis parameters, input data paths, and output settings. This section describes the complete configuration system.

Configuration File Structure
----------------------------

A typical PlasmaFlow configuration file has the following structure:

.. code-block:: yaml

    # PlasmaFlow Configuration File
    analysis:
      output_dir: "/path/to/output"
      temp_dir: "/tmp/plasmaflow"
      n_cores: 8
      log_level: "INFO"
      
    samples:
      - name: "sample1"
        cool_file: "/data/sample1.cool"
        condition: "control"
        replicate: 1
      - name: "sample2"
        cool_file: "/data/sample2.cool" 
        condition: "treatment"
        replicate: 1
        
    loop_calling:
      peakachu_executable: "peakachu"
      model: "GM12878"
      resolution: 10000
      balance: true
      output_format: "bedpe"
      
    quality_control:
      flanking_bp: 100000
      min_separation: 20000
      expected_file: null
      
    differential:
      methods: ["diffHic", "DESeq2"]
      comparisons:
        - name: "treatment_vs_control"
          condition1: "treatment"
          condition2: "control"
          
    visualization:
      flanking_bp: 100000
      bin_size: 3000
      colormap: "RdBu_r" 
      output_format: "pdf"
      
    genomics:
      gtf_file: "/data/genes.gtf"
      distance_threshold: 100000
      tss_flank: 2000
      
    enrichment:
      databases: ["GO_BP", "GO_MF", "GO_CC", "KEGG"]
      include_msigdb: true
      organism: "org.Hs.eg.db"
      
    r_config:
      r_executable: "/usr/bin/R"
      library_path: "/usr/local/lib/R/site-library"

Configuration Sections
----------------------

Analysis Section
~~~~~~~~~~~~~~~

Controls general analysis parameters:

.. code-block:: yaml

    analysis:
      output_dir: "/path/to/output"      # Main output directory
      temp_dir: "/tmp/plasmaflow"        # Temporary files location
      n_cores: 8                         # Number of CPU cores to use
      log_level: "INFO"                  # Logging level (DEBUG, INFO, WARNING, ERROR)
      cache_enabled: true                # Enable result caching
      overwrite: false                   # Overwrite existing results

Samples Section
~~~~~~~~~~~~~~

Defines input samples and metadata:

.. code-block:: yaml

    samples:
      - name: "memB_rep1"               # Unique sample identifier
        cool_file: "/data/memB_rep1.cool"  # Hi-C contact matrix file
        condition: "memB"                # Experimental condition
        replicate: 1                     # Replicate number
        batch: "batch1"                  # Optional batch information
      - name: "prePB_rep1"
        cool_file: "/data/prePB_rep1.cool"
        condition: "prePB"
        replicate: 1

Loop Calling Section
~~~~~~~~~~~~~~~~~~~

Parameters for Peakachu loop calling:

.. code-block:: yaml

    loop_calling:
      peakachu_executable: "peakachu"   # Path to Peakachu executable
      model: "GM12878"                  # Pre-trained model to use
      resolution: 10000                 # Hi-C resolution (bp)
      balance: true                     # Use balanced matrices
      output_format: "bedpe"            # Output format
      custom_model: null                # Path to custom model (optional)
      additional_args: []               # Extra Peakachu arguments

Quality Control Section
~~~~~~~~~~~~~~~~~~~~~~

APA analysis parameters:

.. code-block:: yaml

    quality_control:
      flanking_bp: 100000               # Flanking region size for APA
      min_separation: 20000             # Minimum loop separation
      expected_file: null               # Expected interactions file (optional)
      normalization: "expected"         # Normalization method
      aggregation: "mean"               # Aggregation method

Differential Analysis Section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statistical comparison settings:

.. code-block:: yaml

    differential:
      methods: ["diffHic", "DESeq2", "edgeR"]  # Statistical methods
      comparisons:                             # Pairwise comparisons
        - name: "prePB_vs_memB"
          condition1: "prePB"
          condition2: "memB"
          fdr_threshold: 0.05
          logfc_threshold: 1.0
      filter_low_counts: true                  # Filter low count interactions
      min_count: 5                            # Minimum interaction count

Visualization Section
~~~~~~~~~~~~~~~~~~~~

Plotting and figure generation:

.. code-block:: yaml

    visualization:
      flanking_bp: 100000               # Region size for heatmaps
      bin_size: 3000                    # Bin size for matrices
      colormap: "RdBu_r"                # Matplotlib colormap
      output_format: "pdf"              # Figure format (pdf, png, svg)
      dpi: 300                          # Resolution for raster formats
      figure_size: [12, 8]              # Figure dimensions (inches)

Genomics Section
~~~~~~~~~~~~~~~

Gene analysis parameters:

.. code-block:: yaml

    genomics:
      gtf_file: "/data/gencode.gtf"     # Gene annotation file
      distance_threshold: 100000        # Maximum gene-loop distance (bp)
      tss_flank: 2000                   # TSS flanking region (bp)
      feature_types: ["gene", "transcript"]  # GTF feature types
      gene_biotypes: ["protein_coding"] # Gene biotypes to include

Enrichment Section
~~~~~~~~~~~~~~~~~

Pathway enrichment analysis:

.. code-block:: yaml

    enrichment:
      databases: ["GO_BP", "GO_MF", "GO_CC", "KEGG"]  # Databases to query
      include_msigdb: true                             # Include MSigDB C2
      organism: "org.Hs.eg.db"                         # Organism database
      pvalue_cutoff: 0.05                             # P-value threshold
      qvalue_cutoff: 0.2                              # Q-value threshold

R Configuration Section
~~~~~~~~~~~~~~~~~~~~~~

R environment setup:

.. code-block:: yaml

    r_config:
      r_executable: "/usr/bin/R"                       # R executable path
      library_path: "/usr/local/lib/R/site-library"    # R library path
      packages:                                        # Required packages
        - "diffHic"
        - "DESeq2"
        - "edgeR"
        - "clusterProfiler"

Configuration Validation
-----------------------

PlasmaFlow automatically validates configuration files and provides helpful error messages:

.. code-block:: python

    from plasmaflow.core.config import Config
    
    # Load and validate configuration
    try:
        config = Config.from_yaml("plasmaflow_config.yaml")
        print("Configuration valid!")
    except ConfigurationError as e:
        print(f"Configuration error: {e}")

Environment Variables
--------------------

Some configuration values can be overridden with environment variables:

.. code-block:: bash

    # Override output directory
    export PLASMAFLOW_OUTPUT_DIR="/custom/output/path"
    
    # Override R executable
    export PLASMAFLOW_R_EXECUTABLE="/usr/local/bin/R"
    
    # Override number of cores
    export PLASMAFLOW_N_CORES=16

Configuration Templates
----------------------

Generate configuration templates for different use cases:

.. code-block:: bash

    # Generate basic configuration
    plasmaflow init-config
    
    # Generate configuration with examples
    plasmaflow init-config --with-examples
    
    # Generate minimal configuration  
    plasmaflow init-config --minimal

Best Practices
--------------

1. **Use absolute paths** for all file locations
2. **Validate configuration** before running analysis
3. **Version control** your configuration files
4. **Document custom settings** with comments
5. **Test with small datasets** before full analysis
6. **Backup important configurations** for reproducibility

Example Configurations
---------------------

Complete B-cell Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # B-cell differentiation Hi-C analysis
    analysis:
      output_dir: "/project/bcell_hic/output"
      n_cores: 16
      
    samples:
      - {name: "memB_rep1", cool_file: "/data/memB_rep1.cool", condition: "memB"}
      - {name: "memB_rep2", cool_file: "/data/memB_rep2.cool", condition: "memB"}
      - {name: "prePB_rep1", cool_file: "/data/prePB_rep1.cool", condition: "prePB"}
      - {name: "prePB_rep2", cool_file: "/data/prePB_rep2.cool", condition: "prePB"}
      
    differential:
      comparisons:
        - name: "prePB_vs_memB"
          condition1: "prePB"
          condition2: "memB"
          
    enrichment:
      include_msigdb: true
      databases: ["GO_BP", "KEGG", "MSigDB_C2"]

Minimal Configuration
~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # Minimal configuration for testing
    analysis:
      output_dir: "/tmp/test_output"
      
    samples:
      - name: "test_sample"
        cool_file: "/data/test.cool"
        condition: "test"
        
    loop_calling:
      model: "GM12878"