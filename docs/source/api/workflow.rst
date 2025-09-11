Workflow Module
===============

The workflow module contains the main orchestrator class that coordinates the complete PlasmaFlow analysis pipeline.

Main Workflow Class
------------------

.. automodule:: plasmaflow.workflow
   :members:
   :undoc-members:
   :show-inheritance:

PlasmaFlowAnalysis
~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.workflow.PlasmaFlowAnalysis
   :members:
   :undoc-members:
   :show-inheritance:

Workflow Steps
~~~~~~~~~~~~~

The PlasmaFlowAnalysis class coordinates the following analysis steps:

1. **Loop Calling**: Identify chromatin loops using Peakachu
2. **Quality Control**: Validate loops with APA analysis
3. **Differential Analysis**: Compare loops between conditions
4. **Visualization**: Generate plots and heatmaps
5. **Gene Analysis**: Map genes near differential loops
6. **Pathway Enrichment**: Analyze biological pathways

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow import PlasmaFlowAnalysis
   
   # Initialize workflow
   analyzer = PlasmaFlowAnalysis(config_file="config.yaml")
   
   # Run complete analysis
   results = analyzer.run_complete_workflow()
   
   # Or run individual steps
   loop_results = analyzer.run_loop_calling()
   qc_results = analyzer.run_quality_control()
   diff_results = analyzer.run_differential_analysis()