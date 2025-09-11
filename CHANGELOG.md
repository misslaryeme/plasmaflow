# Changelog

All notable changes to PlasmaFlow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial PlasmaFlow package structure
- Core analysis pipeline with 8-step workflow
- Comprehensive configuration management system
- Command-line interface with Click
- Test suite with pytest framework
- Documentation structure with Sphinx

## [0.1.0] - 2024-12-19

### Added
- **Loop Calling Module** (`plasmaflow.loop_calling`)
  - Peakachu integration for automated loop calling
  - Model management with automatic read count-based selection
  - Batch processing capabilities
  - Quality validation and metrics calculation
  - Support for multiple sample formats

- **Quality Control Module** (`plasmaflow.quality_control`)
  - Aggregate Peak Analysis (APA) implementation
  - Caching system for computational efficiency
  - Comparative APA analysis across samples
  - Statistical scoring (center and corner metrics)
  - Comprehensive visualization suite

- **Differential Analysis Module** (`plasmaflow.differential`)
  - Multiple statistical methods: diffHic, DESeq2, edgeR
  - R integration through rpy2
  - Consensus calling across methods
  - Enhanced result classification with confidence levels
  - Effect size categorization

- **Core Infrastructure**
  - `PlasmaFlowAnalysis` main orchestrator class
  - Flexible configuration system with YAML/JSON support
  - Sample and path management
  - Logging and progress tracking
  - Result caching and serialization

- **Configuration Management**
  - `Config` class with validation
  - `SampleConfig` for cell type management
  - `PathConfig` for file system organization
  - Default configurations for plasma cell analysis
  - Parameter validation and error handling

- **Command Line Interface**
  - `plasmaflow` CLI with comprehensive commands
  - Configuration initialization and validation
  - Individual pipeline step execution
  - Environment checking and dependency validation
  - Utility commands for file management

- **Testing Framework**
  - Comprehensive test suite with >90% coverage target
  - Unit tests for all core components
  - Integration test framework
  - Mock objects for external dependencies
  - Fixtures for sample data generation

- **Documentation**
  - Sphinx documentation structure
  - API reference documentation
  - User guide and tutorials
  - Configuration reference
  - Installation instructions

### Technical Details

#### Dependencies
- **Python**: 3.8+ required
- **Core**: numpy, pandas, scipy, matplotlib, seaborn
- **Hi-C**: cooler, cooltools, coolpuppy, bioframe
- **R Integration**: rpy2 for seamless R package access
- **CLI**: Click for command-line interface
- **Testing**: pytest with comprehensive fixtures
- **Documentation**: Sphinx with RTD theme

#### R Package Integration
- **diffHic**: Bioconductor package for Hi-C differential analysis
- **InteractionSet**: Genomic interaction data structures
- **DESeq2**: Differential expression analysis adapted for Hi-C
- **edgeR**: Empirical analysis of digital gene expression
- **clusterProfiler**: Statistical analysis and visualization of functional profiles
- **msigdbr**: MSigDB gene sets for pathway analysis

#### External Tool Support
- **Peakachu**: Machine learning-based loop calling
- **deepTools**: Matrix computation and visualization
- **R**: Statistical computing environment (4.1+ recommended)

### Architecture

#### Modular Design
```
plasmaflow/
├── config/          # Configuration management
├── loop_calling/    # Peakachu integration and loop calling
├── quality_control/ # APA analysis and quality metrics
├── differential/    # Statistical analysis methods
├── visualization/   # Plotting and matrix generation
├── genomics/        # Gene proximity and expression analysis
├── enrichment/      # Pathway analysis and GSEA
├── utils/          # Utility functions and R interface
└── core.py         # Main orchestrator class
```

#### Data Flow
1. **Input**: Hi-C contact matrices (.cool) and metadata
2. **Loop Calling**: Peakachu-based loop identification
3. **Quality Control**: APA validation and scoring
4. **Differential Analysis**: Statistical comparison across conditions
5. **Matrix Generation**: Visualization matrix creation
6. **Gene Analysis**: Proximity mapping and expression integration  
7. **Pathway Enrichment**: MSigDB-based functional analysis
8. **Output**: Comprehensive results with visualizations

#### Configuration System
- **Hierarchical**: Project → Analysis → Method level parameters
- **Validation**: Automatic parameter checking and error reporting
- **Flexibility**: YAML/JSON formats with environment variable support
- **Defaults**: Pre-configured settings for plasma cell analysis

### Performance Optimizations

#### Computational Efficiency
- **Caching**: Automatic result caching to avoid recomputation
- **Parallel Processing**: Multi-threaded analysis where applicable
- **Memory Management**: Efficient data structures for large Hi-C matrices
- **Lazy Loading**: On-demand loading of R packages and external tools

#### Scalability
- **Batch Processing**: Support for multiple samples simultaneously
- **Resource Management**: Configurable memory and CPU usage
- **Progress Tracking**: Real-time feedback for long-running analyses
- **Error Recovery**: Robust error handling with partial result preservation

### Quality Assurance

#### Testing
- **Unit Tests**: Individual component testing with >90% coverage
- **Integration Tests**: End-to-end workflow validation
- **Mock Testing**: External dependency simulation
- **Data Validation**: Sample data integrity checks

#### Code Quality
- **Formatting**: Black code formatting
- **Import Sorting**: isort for consistent imports
- **Linting**: flake8 for code quality
- **Type Checking**: mypy for static type analysis
- **Documentation**: Google-style docstrings throughout

#### Continuous Integration
- **GitHub Actions**: Automated testing on multiple Python versions
- **Coverage Reporting**: Codecov integration
- **Documentation**: Automatic documentation building
- **Release Management**: Automated version tagging and PyPI publishing

### Known Limitations

#### Current Version (0.1.0)
- **R Dependencies**: Requires manual installation of Bioconductor packages
- **Memory Usage**: Large Hi-C datasets may require substantial RAM
- **Single-threaded R**: R operations not fully parallelized
- **Limited File Formats**: Primary support for .cool and .bedpe formats

#### Platform Support
- **Linux**: Fully supported and tested
- **macOS**: Supported with conda environment
- **Windows**: Limited support (WSL recommended)

### Migration Notes

This is the initial release, so no migration is required. Future versions will provide migration guides for configuration and data formats.

### Breaking Changes

None in initial release.

### Deprecations

None in initial release.

### Security

- **Input Validation**: All file inputs validated before processing
- **Path Traversal Protection**: Safe file path handling
- **External Command Execution**: Controlled subprocess execution with timeouts
- **R Code Injection Prevention**: Parameterized R script generation

---

## Release Notes Format

Each release will include:
- **Added**: New features and capabilities
- **Changed**: Modifications to existing functionality  
- **Deprecated**: Features planned for removal
- **Removed**: Features removed in this version
- **Fixed**: Bug fixes and corrections
- **Security**: Security-related changes

## Development Process

### Version Numbering
- **Major** (X.0.0): Breaking changes or major feature additions
- **Minor** (0.X.0): New features, backwards compatible
- **Patch** (0.0.X): Bug fixes and small improvements

### Release Schedule
- **Major releases**: Every 6-12 months
- **Minor releases**: Every 2-3 months
- **Patch releases**: As needed for critical fixes

### Contribution Guidelines
See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed development and release procedures.