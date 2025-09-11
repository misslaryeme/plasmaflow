Genomics Module
===============

The genomics module provides functionality for gene proximity analysis, expression analysis, and regulatory region classification in relation to chromatin loops.

Main Classes
-----------

.. automodule:: plasmaflow.genomics
   :members:
   :undoc-members:
   :show-inheritance:

Gene Proximity Analyzer
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.genomics.GeneProximityAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Expression Analyzer
~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.genomics.GeneExpressionAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Annotation Manager
~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.genomics.GeneAnnotationManager
   :members:
   :undoc-members:
   :show-inheritance:

Statistical Analysis
~~~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.genomics.GeneDistanceStatistics
   :members:
   :undoc-members:
   :show-inheritance:

Result Classes
~~~~~~~~~~~~~

.. autoclass:: plasmaflow.genomics.ProximityResult
   :members:
   :undoc-members:
   :show-inheritance:

Analysis Features
~~~~~~~~~~~~~~~~

The genomics module provides:

* **TSS-based Distance Calculations**: Accurate gene-loop distance measurements
* **Regulatory Region Classification**: Promoter, enhancer, and gene body annotations
* **Expression Integration**: Combine loop changes with gene expression data
* **Statistical Testing**: Test for enrichment of expression changes near loops
* **Volcano Plots**: Visualize differential expression with loop proximity

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.genomics import GeneProximityAnalyzer, GeneExpressionAnalyzer
   from plasmaflow.core.config import Config
   
   # Initialize analyzers
   config = Config.from_yaml("config.yaml")
   proximity_analyzer = GeneProximityAnalyzer(config)
   expression_analyzer = GeneExpressionAnalyzer(config)
   
   # Analyze gene proximity to differential loops
   proximity_results = proximity_analyzer.analyze_gene_proximity(
       loops_file="differential_loops.bedpe",
       distance_threshold=100000  # 100kb
   )
   
   # Integrate with expression data
   expression_results = expression_analyzer.analyze_expression_overlap(
       proximity_results=proximity_results,
       expression_file="gene_expression.csv"
   )
   
   # Create volcano plot
   volcano_path = expression_analyzer.create_volcano_plot(
       expression_data=expression_results,
       highlight_near_loops=True
   )