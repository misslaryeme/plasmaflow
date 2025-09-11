# PlasmaFlow

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation Status](https://readthedocs.org/projects/plasmaflow/badge/?version=latest)](https://plasmaflow.readthedocs.io/en/latest/?badge=latest)

**PlasmaFlow** is a comprehensive Python package for analyzing Hi-C chromatin interaction data through the lens of plasma cell differentiation. It provides a complete pipeline from loop calling to advanced pathway enrichment analysis with MSigDB integration.

## 🌟 Features

- **Loop Calling**: Automated chromatin loop identification using Peakachu
- **Quality Control**: Comprehensive APA (Aggregate Peak Analysis) for loop validation
- **Differential Analysis**: Multiple statistical methods (diffHic, DESeq2, edgeR) for comparative analysis
- **Matrix Generation**: Integration with deepTools for visualization matrices
- **Gene Proximity Analysis**: TSS-based distance calculations and clustering
- **Expression Analysis**: Differential gene expression with enhanced visualizations
- **Pathway Enrichment**: GSEA and ORA analysis with MSigDB C2 integration for B-cell pathways
- **Seamless R Integration**: Native R package integration through rpy2
- **Comprehensive CLI**: Full command-line interface for pipeline automation
- **Extensible Architecture**: Modular design for easy customization and extension

## 🔬 Scientific Background

PlasmaFlow was designed specifically for analyzing chromatin architecture changes during B-cell to plasma cell differentiation. The pipeline supports analysis of four key cell types:

- **memB**: Memory B cells
- **prePB**: Pre-plasmablasts  
- **PB**: Plasmablasts
- **PC**: Plasma cells

## 🚀 Quick Start

### Installation

```bash
# Install from PyPI (when released)
pip install plasmaflow

# Or install with all dependencies
pip install plasmaflow[all]

# Development installation
git clone https://github.com/misslaryeme/plasmaflow.git
cd plasmaflow
pip install -e .[dev]
```

### Environment Setup

PlasmaFlow requires both Python and R dependencies:

```bash
# Check environment
plasmaflow check-env

# Create conda environment with all dependencies
conda env create -f environment.yml
conda activate plasmaflow
```

### Basic Usage

```python
from plasmaflow import PlasmaFlowAnalysis

# Initialize analysis
analysis = PlasmaFlowAnalysis(config_file="config.yaml")

# Run complete pipeline
results = analysis.run_full_pipeline()

# Access results
loop_calling_results = results['loop_calling']
apa_results = results['quality_control']
differential_results = results['differential_analysis']
```

### Command Line Interface

```bash
# Initialize configuration
plasmaflow init-config my_config.yaml

# Run complete pipeline
plasmaflow --config my_config.yaml run

# Run individual steps
plasmaflow --config my_config.yaml loop-calling --samples samples.csv
plasmaflow differential --counts counts.csv --metadata metadata.csv
```

## 📖 Documentation

Comprehensive documentation is available at [plasmaflow.readthedocs.io](https://plasmaflow.readthedocs.io/)

### Key Documentation Sections

- [**Installation Guide**](https://plasmaflow.readthedocs.io/en/latest/installation.html): Detailed setup instructions
- [**Tutorial**](https://plasmaflow.readthedocs.io/en/latest/tutorials/index.html): Step-by-step analysis walkthrough
- [**API Reference**](https://plasmaflow.readthedocs.io/en/latest/api/index.html): Complete function documentation
- [**Configuration**](https://plasmaflow.readthedocs.io/en/latest/configuration.html): Parameter customization guide
- [**Examples**](https://plasmaflow.readthedocs.io/en/latest/examples/index.html): Jupyter notebook examples

## 🧬 Pipeline Overview

PlasmaFlow implements an 8-step analysis workflow:

### 1. Loop Calling (`plasmaflow.loop_calling`)
- **Method**: Peakachu framework
- **Input**: Hi-C contact matrices (.cool format)
- **Output**: Loop coordinates (.bedpe format)
- **Features**: Automatic model selection, quality validation, batch processing

### 2. Quality Control (`plasmaflow.quality_control`)
- **Method**: Aggregate Peak Analysis (APA)
- **Input**: Hi-C matrices + loop coordinates
- **Output**: Quality metrics, APA visualizations
- **Features**: Caching, comparative analysis, statistical scoring

### 3. Differential Analysis (`plasmaflow.differential`)
- **Methods**: diffHic, DESeq2, edgeR
- **Input**: Count matrices + sample metadata
- **Output**: Differential loop lists, statistical results
- **Features**: Multiple method comparison, consensus calling

### 4. Matrix Generation (`plasmaflow.visualization`)
- **Method**: deepTools integration
- **Input**: Differential results + reference files
- **Output**: Visualization matrices, heatmaps
- **Features**: Multiple file formats, customizable parameters

### 5. Comparative APA (`plasmaflow.quality_control`)
- **Method**: Pairwise comparative analysis
- **Input**: Multiple sample APA results
- **Output**: Statistical comparisons, plots
- **Features**: Cross-sample validation, effect size analysis

### 6. Gene Proximity (`plasmaflow.genomics`)
- **Method**: TSS-based distance calculations
- **Input**: Differential loops + gene annotations
- **Output**: Gene-loop associations, clustering
- **Features**: Multiple distance metrics, cluster analysis

### 7. Expression Analysis (`plasmaflow.genomics`)
- **Method**: Differential gene expression
- **Input**: Expression data + loop-gene associations
- **Output**: Enhanced volcano plots, categorized results
- **Features**: Comprehensive subgroup analysis

### 8. Pathway Enrichment (`plasmaflow.enrichment`)
- **Method**: clusterProfiler + MSigDB
- **Input**: Differential genes from loop-associated regions
- **Output**: GSEA/ORA results, pathway visualizations
- **Features**: B-cell specific pathways, multiple gene set databases

## 🛠️ Configuration

PlasmaFlow uses YAML configuration files for parameter management:

```yaml
# Example configuration
project_name: "My_PlasmaFlow_Analysis"
resolution: 10000
n_threads: 8

# Input/Output paths
input_dir: "/path/to/hic/data"
output_dir: "/path/to/results"

# Sample configuration
cell_types: ["memB", "prePB", "PB", "PC"]

# Analysis parameters
loop_calling:
  threshold: 0.95
  models:
    memB: {model: "100million", reads: 113922292}
    prePB: {model: "150million", reads: 159514688}

quality_control:
  apa:
    flank: 200000
    min_diag: 3
    maxdist: 1000000

differential:
  methods: ["diffHic", "DESeq2", "edgeR"]
  fdr_threshold: 0.05
  logfc_threshold: 1.0
  comparisons:
    - {control: "memB", treatment: "prePB", name: "prePB_vs_memB"}
```

## 🔧 Requirements

### Python Dependencies
- Python 3.8+
- Core: numpy, pandas, scipy, matplotlib, seaborn
- Hi-C: cooler, cooltools, coolpuppy, bioframe
- Statistics: scikit-learn, statsmodels
- R integration: rpy2

### R Dependencies
- Bioconductor: diffHic, InteractionSet, GenomicRanges
- Statistics: DESeq2, edgeR
- Pathways: clusterProfiler, msigdbr
- Visualization: ggplot2

### External Tools
- **Peakachu**: Loop calling (installable via conda)
- **deepTools**: Matrix computation (conda/pip)
- **R**: Statistical computing environment

### Development Setup

```bash
# Clone repository
git clone https://github.com/misslaryeme/plasmaflow.git
cd plasmaflow

# Install development dependencies
pip install -e .[dev]

# Install pre-commit hooks
pre-commit install

# Validate setup
plasmaflow check-env
```

## 📚 Citation

If you use PlasmaFlow in your research, please cite:

```bibtex
@software{plasmaflow2025,
  title={PlasmaFlow: Hi-C Chromatin Loop Analysis for Plasma Cell Differentiation},
  author={Miss Leriem ZELLAGUI},
  year={2025},
  url={https://github.com/misslaryeme/plasmaflow},
  version={0.1.0}
}
```

## 🆘 Support

### Getting Help
- **Documentation**: [plasmaflow.readthedocs.io](https://plasmaflow.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/misslaryeme/plasmaflow/issues)
- **Discussions**: [GitHub Discussions](https://github.com/misslaryeme/plasmaflow/discussions)

## 📄 License

PlasmaFlow is released under the MIT License. See [LICENSE](LICENSE) for details.
