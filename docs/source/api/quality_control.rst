Quality Control Module
======================

The quality_control module provides comprehensive quality assessment for chromatin loops using Aggregate Peak Analysis (APA) and comparative analysis.

Main Classes
-----------

.. automodule:: plasmaflow.quality_control
   :members:
   :undoc-members:
   :show-inheritance:

APA Analyzer
~~~~~~~~~~~

.. autoclass:: plasmaflow.quality_control.APAAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Comparative APA
~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.quality_control.ComparativeAPAAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Result Classes
~~~~~~~~~~~~~

.. autoclass:: plasmaflow.quality_control.APAResult
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: plasmaflow.quality_control.ComparativeAPAResult
   :members:
   :undoc-members:
   :show-inheritance:

Visualization
~~~~~~~~~~~~

.. autoclass:: plasmaflow.quality_control.APAPlotter
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.quality_control import APAAnalyzer
   from plasmaflow.core.config import Config
   
   # Initialize analyzer
   config = Config.from_yaml("config.yaml")
   apa_analyzer = APAAnalyzer(config)
   
   # Run APA analysis
   apa_result = apa_analyzer.run_apa_analysis(
       loops_file="loops.bedpe",
       cool_file="sample.cool"
   )
   
   # Generate plots
   plot_path = apa_analyzer.plot_apa_results(apa_result)
   
   # Comparative analysis
   comp_analyzer = ComparativeAPAAnalyzer(config)
   comp_results = comp_analyzer.run_pairwise_comparison(
       sample_pairs=[("sample1", "sample2")]
   )