Enrichment Module
=================

The enrichment module provides comprehensive pathway enrichment analysis using R's clusterProfiler package and MSigDB gene sets, specifically tailored for B-cell and plasma cell biology.

Main Classes
-----------

.. automodule:: plasmaflow.enrichment
   :members:
   :undoc-members:
   :show-inheritance:

ClusterProfiler Analyzer
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.enrichment.ClusterProfilerAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

MSigDB Analyzer
~~~~~~~~~~~~~~

.. autoclass:: plasmaflow.enrichment.MSigDBAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

R Interface
~~~~~~~~~~

.. autoclass:: plasmaflow.enrichment.RPathwayInterface
   :members:
   :undoc-members:
   :show-inheritance:

Visualization
~~~~~~~~~~~~

.. autoclass:: plasmaflow.enrichment.PathwayPlotter
   :members:
   :undoc-members:
   :show-inheritance:

Result Classes
~~~~~~~~~~~~~

.. autoclass:: plasmaflow.enrichment.PathwayEnrichmentResult
   :members:
   :undoc-members:
   :show-inheritance:

Enrichment Methods
~~~~~~~~~~~~~~~~~

The module supports two main enrichment approaches:

* **Over-Representation Analysis (ORA)**: Tests for enrichment of significant genes in pathways
* **Gene Set Enrichment Analysis (GSEA)**: Uses ranked gene lists to identify enriched pathways

Supported Databases
~~~~~~~~~~~~~~~~~~

* **Gene Ontology (GO)**: Biological Process, Molecular Function, Cellular Component
* **KEGG**: Kyoto Encyclopedia of Genes and Genomes pathways
* **MSigDB C2**: Curated gene sets including Reactome and BioCarta
* **Custom B-cell Sets**: Plasma cell differentiation specific pathways

Example Usage
~~~~~~~~~~~~

.. code-block:: python

   from plasmaflow.enrichment import ClusterProfilerAnalyzer
   from plasmaflow.core.config import Config
   
   # Initialize analyzer
   config = Config.from_yaml("config.yaml")
   enrichment_analyzer = ClusterProfilerAnalyzer(config)
   
   # Run ORA analysis
   significant_genes = ["PRDM1", "XBP1", "IRF4", "CD38", "CD138"]
   ora_results = enrichment_analyzer.run_ora_analysis(
       gene_list=significant_genes,
       databases=["GO_BP", "KEGG", "MSigDB_C2"]
   )
   
   # Run GSEA analysis
   ranked_genes = {"PRDM1": 3.2, "XBP1": 2.8, "IRF4": 2.1}
   gsea_results = enrichment_analyzer.run_gsea_analysis(
       ranked_genes=ranked_genes,
       databases=["GO_BP", "KEGG"]
   )
   
   # Comprehensive analysis
   all_results = enrichment_analyzer.run_comprehensive_analysis(
       gene_data=expression_df,
       include_msigdb=True
   )
   
   # Create visualizations
   from plasmaflow.enrichment import create_comprehensive_pathway_plots
   plot_paths = create_comprehensive_pathway_plots(
       results=all_results,
       config=config
   )