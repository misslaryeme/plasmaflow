Visualization Module
====================

The visualization module provides comprehensive plotting capabilities for Hi-C analysis results, including heatmaps, APA plots, and publication-ready figures.

Main Classes
-----------

.. automodule:: plasmaflow.visualization
   :members:
   :undoc-members:
   :show-inheritance:

Heatmap Plotter
~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.visualization.HeatmapPlotter
   :members:
   :undoc-members:
   :show-inheritance:

APA Plotter
~~~~~~~~~~

.. autoclass:: plasmaflow.visualization.APAPlotter
   :members:
   :undoc-members:
   :show-inheritance:

Comparative Plotter
~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.visualization.ComparativePlotter
   :members:
   :undoc-members:
   :show-inheritance:

Matrix Generator
~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.visualization.MatrixGenerator
   :members:
   :undoc-members:
   :show-inheritance:

Plot Types
~~~~~~~~~

The visualization module supports various plot types:

* **Heatmaps**: Differential loop intensity visualization
* **APA Plots**: Aggregate peak analysis validation
* **Comparative Plots**: Side-by-side condition comparisons
* **Violin Plots**: Distribution comparisons
* **Scatter Plots**: Correlation analysis
* **Matrix Plots**: Hi-C contact matrix visualization

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.visualization import HeatmapPlotter
   from plasmaflow.core.config import Config
   
   # Initialize plotter
   config = Config.from_yaml("config.yaml")
   plotter = HeatmapPlotter(config)
   
   # Create differential heatmap
   heatmap_path = plotter.create_differential_heatmap(
       differential_results=diff_results,
       output_file="differential_heatmap.pdf"
   )
   
   # Generate comprehensive plot set
   plot_paths = plotter.create_comprehensive_plots(
       results_dict=all_results,
       output_prefix="comprehensive"
   )
   
   # Create custom matrix plot
   matrix_plot = plotter.plot_contact_matrix(
       cool_file="sample.cool",
       region="chr1:1000000-2000000",
       loops_file="loops.bedpe"
   )