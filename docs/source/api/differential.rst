Differential Analysis Module
============================

The differential module provides statistical analysis of chromatin loop changes between conditions using multiple methods including diffHic, DESeq2, and edgeR.

Main Classes
-----------

.. automodule:: plasmaflow.differential
   :members:
   :undoc-members:
   :show-inheritance:

Differential Analyzer
~~~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.differential.DifferentialAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

R Integration
~~~~~~~~~~~~

.. autoclass:: plasmaflow.differential.RDifferentialInterface
   :members:
   :undoc-members:
   :show-inheritance:

Result Classes
~~~~~~~~~~~~~

.. autoclass:: plasmaflow.differential.DifferentialResult
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: plasmaflow.differential.ComparisonResult
   :members:
   :undoc-members:
   :show-inheritance:

Statistical Methods
~~~~~~~~~~~~~~~~~

The module supports multiple statistical approaches:

* **diffHic**: Specialized Hi-C differential analysis from Bioconductor
* **DESeq2**: Popular RNA-seq method adapted for Hi-C data
* **edgeR**: Another RNA-seq method for count-based analysis

Each method has different strengths:

- diffHic: Designed specifically for Hi-C data, handles technical biases
- DESeq2: Robust normalization and dispersion estimation
- edgeR: Fast computation, good for large datasets

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.differential import DifferentialAnalyzer
   from plasmaflow.core.config import Config
   
   # Initialize analyzer
   config = Config.from_yaml("config.yaml")
   analyzer = DifferentialAnalyzer(config)
   
   # Run differential analysis
   results = analyzer.run_differential_analysis(
       sample_info="samples.csv",
       comparison="treatment_vs_control"
   )
   
   # Get significant loops
   significant_loops = results.get_significant_loops(
       fdr_threshold=0.05,
       logfc_threshold=1.0
   )
   
   # Export results
   results.save_results("output/differential_results.csv")