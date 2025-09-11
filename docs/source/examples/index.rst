Examples
========

This section contains practical examples and case studies demonstrating PlasmaFlow usage for real-world Hi-C analysis scenarios.

.. toctree::
   :maxdepth: 2

   bcell_differentiation
   comparative_analysis  
   integration_examples
   custom_workflows

Example Categories
-----------------

**Biological Applications:**

* :doc:`bcell_differentiation`: B-cell to plasma cell differentiation analysis
* :doc:`comparative_analysis`: Multi-condition Hi-C comparison studies

**Technical Examples:**

* :doc:`integration_examples`: Integrating with other genomics data types
* :doc:`custom_workflows`: Building custom analysis pipelines

Example Data
-----------

All examples use publicly available datasets or simulated data that demonstrates key analysis concepts.

Download Example Data
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Download all example datasets
    plasmaflow download-examples --all
    
    # Download specific example
    plasmaflow download-examples --example bcell_differentiation

Running Examples
---------------

Each example can be run as:

1. **Jupyter Notebooks**: Interactive analysis with explanations
2. **Python Scripts**: Automated execution for batch processing  
3. **Command Line**: Using PlasmaFlow CLI interface

Example Structure
----------------

Each example includes:

- **Background**: Biological context and research questions
- **Data Description**: Input datasets and experimental design
- **Analysis Steps**: Complete workflow with code
- **Results**: Expected outputs and interpretation
- **Extensions**: Ideas for further analysis

Performance Notes
----------------

Examples are designed to run on standard hardware. For large datasets:

- Use the provided configuration templates
- Consider cloud computing resources
- Enable result caching for iterative analysis

Getting Started
--------------

Start with the :doc:`bcell_differentiation` example to see a complete PlasmaFlow workflow applied to a real biological question.