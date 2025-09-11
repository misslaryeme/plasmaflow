# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath("../../"))

# Mock heavy dependencies for Read the Docs
class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

MOCK_MODULES = [
    'scipy',
    'scipy.stats',
    'scipy.sparse',
    'scipy.ndimage',
    'matplotlib',
    'matplotlib.pyplot',
    'seaborn',
    'plotly',
    'plotly.graph_objects',
    'plotly.express',
    'cooler',
    'cooltools',
    'coolpuppy',
    'bioframe',
    'pairix',
    'rpy2',
    'rpy2.robjects',
    'pybedtools',
    'pyBigWig',
    'pysam',
    'h5py',
    'tables',
    'openpyxl',
    'xlrd',
    'joblib',
    'multiprocessing_logging',
    'sklearn',
    'scikit-learn',
    'statsmodels',
    'networkx',
    'click',
    'colorlog',
    'tqdm'
]

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = Mock()

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "PlasmaFlow"
copyright = "2025, PlasmaFlow Contributors"
author = "PlasmaFlow Contributors"

version = "1.0.0"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Extension configuration -------------------------------------------------

# Autodoc configuration
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# Todo extension
todo_include_todos = True

# Typehints
typehints_fully_qualified = False
always_document_param_types = True
