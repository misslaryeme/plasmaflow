Loop Calling Module
===================

The loop_calling module provides functionality for identifying chromatin loops from Hi-C contact matrices using the Peakachu framework.

Main Classes
-----------

.. automodule:: plasmaflow.loop_calling
   :members:
   :undoc-members:
   :show-inheritance:

Peakachu Analyzer
~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.loop_calling.PeakachuAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

Result Classes
~~~~~~~~~~~~~

.. autoclass:: plasmaflow.loop_calling.LoopCallingResult
   :members:
   :undoc-members:
   :show-inheritance:

Model Management
~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.loop_calling.PeakachuModelManager
   :members:
   :undoc-members:
   :show-inheritance:

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.loop_calling import PeakachuAnalyzer
   from plasmaflow.core.config import Config
   
   # Load configuration
   config = Config.from_yaml("config.yaml")
   
   # Initialize analyzer
   analyzer = PeakachuAnalyzer(config)
   
   # Run loop calling for single sample
   result = analyzer.run_single_analysis(
       cool_file="sample.cool",
       output_prefix="sample_loops"
   )
   
   # Run batch analysis
   batch_results = analyzer.run_batch_analysis()