Tutorials
=========

This section provides step-by-step tutorials for using PlasmaFlow to analyze Hi-C chromatin loop data.

.. toctree::
   :maxdepth: 2

   basic_workflow
   advanced_analysis
   pathway_enrichment
   custom_analysis

Tutorial Overview
----------------

**Basic Tutorials:**

* :doc:`basic_workflow`: Complete workflow from loop calling to visualization
* :doc:`advanced_analysis`: Advanced statistical analysis and interpretation

**Specialized Tutorials:**

* :doc:`pathway_enrichment`: Comprehensive pathway enrichment analysis
* :doc:`custom_analysis`: Customizing analysis for specific research questions

Before You Start
----------------

Prerequisites
~~~~~~~~~~~~

- Basic understanding of Hi-C data analysis
- Familiarity with Python programming
- Understanding of statistical concepts (p-values, FDR correction)
- Knowledge of pathway analysis concepts

Required Data
~~~~~~~~~~~~

For the tutorials, you'll need:

- Hi-C contact matrices in .cool format
- Sample metadata file
- Gene annotation file (GTF format)
- (Optional) Gene expression data

Getting Tutorial Data
~~~~~~~~~~~~~~~~~~~~

Download example datasets:

.. code-block:: bash

    # Download tutorial data
    plasmaflow download-examples --output-dir tutorial_data/
    
    # Or use your own data following the format specifications

Tutorial Structure
-----------------

Each tutorial includes:

- **Learning objectives** 
- **Required input data**
- **Step-by-step instructions**
- **Expected outputs**
- **Interpretation guidelines**
- **Common troubleshooting**

Difficulty Levels
----------------

Tutorials are marked with difficulty levels:

- 🟢 **Beginner**: Basic workflow execution
- 🟡 **Intermediate**: Parameter customization and interpretation
- 🔴 **Advanced**: Method comparison and custom analysis

Next Steps
----------

Start with the :doc:`basic_workflow` tutorial to learn the fundamentals of PlasmaFlow analysis.