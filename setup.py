#!/usr/bin/env python3

import pathlib

from setuptools import find_packages, setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# Version
__version__ = "0.1.0"

setup(
    name="plasmaflow",
    version=__version__,
    author="PlasmaFlow Development Team",
    author_email="plasmaflow@example.com",
    description="A comprehensive Hi-C chromatin loop analysis pipeline for plasma cell differentiation",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/misslaryeme/plasmaflow",
    project_urls={
        "Bug Reports": "https://github.com/misslaryeme/plasmaflow/issues",
        "Source": "https://github.com/misslaryeme/plasmaflow",
        "Documentation": "https://plasmaflow.readthedocs.io/",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        # Core scientific computing
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        # Visualization
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "plotly>=5.0.0",
        # Hi-C specific
        "cooler>=0.9.0",
        "cooltools>=0.5.0",
        "coolpuppy>=1.1.0",
        "bioframe>=0.4.0",
        "pairix>=0.3.7",
        # R integration
        "rpy2>=3.5.0",
        # Genomics
        "pybedtools>=0.9.0",
        "pyBigWig>=0.3.18",
        "pysam>=0.19.0",
        # Configuration and utilities
        "pyyaml>=6.0",
        "tqdm>=4.60.0",
        "click>=8.0.0",
        "colorlog>=6.6.0",
        "importlib-resources>=5.0.0;python_version<'3.9'",
        # Data handling
        "h5py>=3.6.0",
        "tables>=3.7.0",
        "openpyxl>=3.0.9",
        "xlrd>=2.0.1",
        # Parallel processing
        "joblib>=1.1.0",
        "multiprocessing-logging>=0.3.4",
        # Statistical analysis
        "scikit-learn>=1.0.0",
        "statsmodels>=0.13.0",
        # Network analysis
        "networkx>=2.6.0",
    ],
    extras_require={
        "dev": [
            "black>=22.0.0",
            "isort>=5.10.0",
            "flake8>=5.0.0",
            "mypy>=0.991",
            "pre-commit>=2.20.0",
        ],
        "docs": [
            "sphinx>=5.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "sphinx-autodoc-typehints>=1.19.0",
            "nbsphinx>=0.8.9",
            "pandoc>=2.0.0",
        ],
        "notebook": [
            "jupyter>=1.0.0",
            "ipykernel>=6.15.0",
            "ipywidgets>=8.0.0",
        ],
        "gpu": [
            "cupy>=11.0.0",  # GPU-accelerated NumPy
        ],
    },
    entry_points={
        "console_scripts": [
            "plasmaflow=plasmaflow.cli:main",
            "pf=plasmaflow.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "plasmaflow": [
            "config/*.yaml",
            "config/*.json",
            "data/*.bed",
            "data/*.gtf",
        ],
    },
    zip_safe=False,
    keywords=[
        "Hi-C",
        "chromatin",
        "loops",
        "3D-genomics",
        "bioinformatics",
        "plasma-cells",
        "B-cells",
        "transcriptomics",
        "epigenomics",
        "differential-analysis",
        "pathway-enrichment",
        "MSigDB",
    ],
)
