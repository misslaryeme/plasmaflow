Installation
============

Requirements
------------

PlasmaFlow requires Python 3.8 or higher and has the following dependencies:

**Core Python packages:**

* numpy >= 1.20.0
* pandas >= 1.3.0
* scipy >= 1.7.0
* matplotlib >= 3.4.0
* seaborn >= 0.11.0
* click >= 8.0.0
* pyyaml >= 5.4.0
* tqdm >= 4.60.0

**Hi-C analysis packages:**

* cooler >= 0.8.11
* cooltools >= 0.4.1
* coolpuppy >= 0.9.2

**R integration:**

* rpy2 >= 3.4.5 (optional, for advanced statistical analysis)

**External dependencies:**

* **Peakachu**: Loop calling framework (https://github.com/tariks/peakachu)
* **deepTools**: For matrix computations (https://deeptools.readthedocs.io/)
* **R**: For Bioconductor packages (diffHic, DESeq2, edgeR, clusterProfiler)

Installation Methods
-------------------

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    pip install plasmaflow

From Source
~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/your-org/plasmaflow.git
    cd plasmaflow
    pip install -e .

Development Installation
~~~~~~~~~~~~~~~~~~~~~~~

For development and contributing:

.. code-block:: bash

    git clone https://github.com/your-org/plasmaflow.git
    cd plasmaflow
    pip install -e ".[dev]"
    
    # Install pre-commit hooks
    pre-commit install

External Dependencies
--------------------

Peakachu Installation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Using conda (recommended)
    conda install -c bioconda peakachu
    
    # Or from source
    git clone https://github.com/tariks/peakachu.git
    cd peakachu
    pip install .

deepTools Installation
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Using conda
    conda install -c bioconda deeptools
    
    # Using pip
    pip install deeptools

R and Bioconductor Setup
~~~~~~~~~~~~~~~~~~~~~~~

Install R (version 4.0 or higher) and required Bioconductor packages:

.. code-block:: r

    # In R console
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install(c(
        "diffHic",
        "DESeq2", 
        "edgeR",
        "clusterProfiler",
        "org.Hs.eg.db",
        "DOSE",
        "ReactomePA",
        "enrichplot",
        "msigdbr"
    ))

Docker Installation
------------------

For a complete environment with all dependencies:

.. code-block:: bash

    # Pull the Docker image
    docker pull plasmaflow/plasmaflow:latest
    
    # Run with mounted data directory
    docker run -v /path/to/data:/data plasmaflow/plasmaflow:latest

Verification
-----------

Test your installation:

.. code-block:: bash

    # Test CLI
    plasmaflow --help
    
    # Test Python import
    python -c "import plasmaflow; print(plasmaflow.__version__)"
    
    # Test R integration (optional)
    plasmaflow validate-r

Configuration
------------

Create a basic configuration file:

.. code-block:: bash

    plasmaflow init-config

This creates a ``plasmaflow_config.yaml`` file that you can customize for your analysis.

Troubleshooting
--------------

Common Issues
~~~~~~~~~~~~

**ImportError for cooler/cooltools:**

.. code-block:: bash

    # Install Hi-C packages
    conda install -c bioconda cooler cooltools coolpuppy

**R integration issues:**

.. code-block:: bash

    # Verify R installation
    which R
    R --version
    
    # Test rpy2 installation
    python -c "import rpy2.robjects as ro; print(ro.r('R.version.string'))"

**Peakachu not found:**

.. code-block:: bash

    # Ensure Peakachu is in PATH
    which peakachu
    
    # If not, add to PATH or specify full path in config

Getting Help
-----------

If you encounter installation issues:

1. Check the `GitHub Issues <https://github.com/your-org/plasmaflow/issues>`_
2. Create a new issue with your system information and error message
3. Join our `Discussion Forum <https://github.com/your-org/plasmaflow/discussions>`_