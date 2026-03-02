# -- Path setup --------------------------------------------------------------

import os
import sys

# Add the src directory to sys.path
sys.path.insert(0, os.path.abspath('../../src'))


# -- Project information -----------------------------------------------------

project = 'Pysotope'
copyright = '2026, Dr. Gerard Otiniano and Dr. Elizabeth Thomas'
author = 'Dr. Gerard Otiniano and Dr. Elizabeth Thomas'
release = '1.12.3'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
]

templates_path = ['_templates']
exclude_patterns = []

# Automatically include class __init__ docstrings
autoclass_content = 'both'

# Show type hints in description
autodoc_typehints = 'description'


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']