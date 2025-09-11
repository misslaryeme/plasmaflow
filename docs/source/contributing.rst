Contributing to PlasmaFlow
=========================

We welcome contributions to PlasmaFlow! This guide will help you get started with contributing to the project.

Getting Started
--------------

Setting Up Development Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Fork and clone the repository:**

   .. code-block:: bash

       git clone https://github.com/misslaryeme/plasmaflow.git
       cd plasmaflow

2. **Create a development environment:**

   .. code-block:: bash

       # Using conda
       conda create -n plasmaflow-dev python=3.8
       conda activate plasmaflow-dev
       
       # Or using virtualenv
       python -m venv plasmaflow-dev
       source plasmaflow-dev/bin/activate

3. **Install in development mode:**

   .. code-block:: bash

       pip install -e ".[dev]"

4. **Install pre-commit hooks:**

   .. code-block:: bash

       pre-commit install

Development Workflow
-------------------

1. **Create a feature branch:**

   .. code-block:: bash

       git checkout -b feature/your-feature-name

2. **Make your changes and test them:**

   .. code-block:: bash

       # Run tests
       pytest
       
       # Run type checking
       mypy plasmaflow
       
       # Run linting
       flake8 plasmaflow
       black plasmaflow

3. **Commit your changes:**

   .. code-block:: bash

       git add .
       git commit -m "Add your descriptive commit message"

4. **Push and create a pull request:**

   .. code-block:: bash

       git push origin feature/your-feature-name

Code Style Guidelines
--------------------

Python Code Style
~~~~~~~~~~~~~~~~~

We follow PEP 8 with some modifications:

- **Line length**: 88 characters (Black default)
- **Import sorting**: Use isort
- **Type hints**: Use type hints for all public functions
- **Docstrings**: Use Google-style docstrings

Example:

.. code-block:: python

    from typing import List, Optional, Dict, Any
    import pandas as pd
    
    
    def analyze_loops(
        loops_file: str,
        cool_file: str,
        output_dir: Optional[str] = None
    ) -> Dict[str, Any]:
        """Analyze chromatin loops from Hi-C data.
        
        Args:
            loops_file: Path to BEDPE file with loops
            cool_file: Path to cooler file with Hi-C data
            output_dir: Output directory for results
            
        Returns:
            Dictionary containing analysis results
            
        Raises:
            FileNotFoundError: If input files don't exist
            ValueError: If file formats are invalid
        """
        # Implementation here
        pass

Documentation Style
~~~~~~~~~~~~~~~~~~

- Use **reStructuredText** for documentation
- Include **code examples** in docstrings and documentation
- Keep **API documentation** up to date
- Add **type information** in docstrings

Testing Guidelines
-----------------

Test Organization
~~~~~~~~~~~~~~~~

Tests are organized in the ``tests/`` directory:

.. code-block:: text

    tests/
    ├── unit/                 # Unit tests
    │   ├── test_loop_calling.py
    │   ├── test_quality_control.py
    │   └── ...
    ├── integration/          # Integration tests
    │   ├── test_workflows.py
    │   └── ...
    ├── fixtures/            # Test data
    │   ├── sample.cool
    │   ├── loops.bedpe
    │   └── ...
    └── conftest.py          # Pytest configuration

Writing Tests
~~~~~~~~~~~~

Use pytest for all tests:

.. code-block:: python

    import pytest
    from plasmaflow.loop_calling import PeakachuAnalyzer
    from plasmaflow.core.config import Config


    class TestPeakachuAnalyzer:
        
        def test_init(self, sample_config):
            """Test analyzer initialization."""
            analyzer = PeakachuAnalyzer(sample_config)
            assert analyzer.config == sample_config
            
        def test_run_analysis(self, sample_config, sample_cool_file):
            """Test loop calling analysis."""
            analyzer = PeakachuAnalyzer(sample_config)
            result = analyzer.run_single_analysis(
                cool_file=sample_cool_file,
                output_prefix="test"
            )
            
            assert result.loops_file.exists()
            assert len(result.loops_df) > 0

Test Coverage
~~~~~~~~~~~~

Maintain high test coverage:

.. code-block:: bash

    # Run tests with coverage
    pytest --cov=plasmaflow --cov-report=html
    
    # View coverage report
    open htmlcov/index.html

Documentation
------------

Building Documentation
~~~~~~~~~~~~~~~~~~~~~

Documentation is built with Sphinx:

.. code-block:: bash

    # Build HTML documentation
    cd docs
    make html
    
    # View documentation
    open build/html/index.html

Adding Documentation
~~~~~~~~~~~~~~~~~~~

1. **API documentation**: Auto-generated from docstrings
2. **Tutorials**: Add to ``docs/source/tutorials/``
3. **Examples**: Add to ``docs/source/examples/``
4. **Configuration**: Update ``docs/source/configuration.rst``

Contribution Types
-----------------

Bug Reports
~~~~~~~~~~

When reporting bugs, please include:

- PlasmaFlow version
- Python version
- Operating system
- Complete error traceback
- Minimal reproduction example
- Expected vs actual behavior

Feature Requests
~~~~~~~~~~~~~~~

For feature requests, please provide:

- Clear description of the feature
- Use case and motivation
- Proposed API (if applicable)
- Willingness to implement

Code Contributions
~~~~~~~~~~~~~~~~~

Types of contributions we welcome:

- **Bug fixes**
- **New features**
- **Performance improvements**
- **Documentation improvements**
- **Test coverage improvements**
- **Example notebooks**

Pull Request Process
-------------------

1. **Check existing issues** and PRs to avoid duplication
2. **Discuss major changes** in an issue first
3. **Write tests** for new functionality
4. **Update documentation** as needed
5. **Ensure CI passes** before requesting review
6. **Respond to review feedback** promptly

PR Requirements
~~~~~~~~~~~~~~

All pull requests must:

- Pass all existing tests
- Include tests for new functionality
- Maintain or improve test coverage
- Follow code style guidelines
- Update relevant documentation
- Include a clear description

Release Process
--------------

Versioning
~~~~~~~~~

We use semantic versioning (SemVer):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

Release Checklist
~~~~~~~~~~~~~~~~

1. Update version in ``setup.py``
2. Update ``CHANGELOG.md``
3. Run full test suite
4. Build and test documentation
5. Create release PR
6. Tag release after merge
7. Deploy to PyPI

Community Guidelines
-------------------

Code of Conduct
~~~~~~~~~~~~~~

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and constructive
- Focus on the technical aspects
- Help others learn and grow
- Report unacceptable behavior

Getting Help
~~~~~~~~~~~

If you need help contributing:

- Check existing documentation
- Search closed issues and PRs
- Ask questions in discussions
- Contact maintainers directly

Recognition
----------

All contributors are recognized in:

- ``AUTHORS.md`` file
- Release notes
- Documentation acknowledgments

Thank you for contributing to PlasmaFlow!