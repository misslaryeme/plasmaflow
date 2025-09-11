API Reference
=============

This section contains the complete API documentation for all PlasmaFlow modules.

Core Modules
-----------

.. toctree::
   :maxdepth: 2

   core
   workflow

Analysis Modules  
----------------

.. toctree::
   :maxdepth: 2

   loop_calling
   quality_control
   differential
   visualization
   genomics
   enrichment

Module Overview
--------------

Core Infrastructure
~~~~~~~~~~~~~~~~~~

* :doc:`core`: Configuration, exceptions, and logging utilities
* :doc:`workflow`: Main analysis orchestrator class

Analysis Pipeline
~~~~~~~~~~~~~~~~

* :doc:`loop_calling`: Peakachu-based chromatin loop identification
* :doc:`quality_control`: APA validation and quality metrics
* :doc:`differential`: Statistical comparison of loop changes
* :doc:`visualization`: Publication-ready plots and heatmaps
* :doc:`genomics`: Gene proximity and expression analysis
* :doc:`enrichment`: Pathway enrichment with R integration

Quick Reference
--------------

Main Classes
~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   plasmaflow.PlasmaFlowAnalysis
   plasmaflow.loop_calling.PeakachuAnalyzer
   plasmaflow.quality_control.APAAnalyzer
   plasmaflow.differential.DifferentialAnalyzer
   plasmaflow.visualization.HeatmapPlotter
   plasmaflow.genomics.GeneProximityAnalyzer
   plasmaflow.enrichment.ClusterProfilerAnalyzer

Configuration
~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   plasmaflow.core.config.Config
   plasmaflow.core.config.AnalysisConfig
   plasmaflow.core.config.SampleConfig

Exceptions
~~~~~~~~~

.. autosummary::
   :nosignatures:

   plasmaflow.core.exceptions.PlasmaFlowError
   plasmaflow.core.exceptions.ConfigurationError
   plasmaflow.core.exceptions.AnalysisError