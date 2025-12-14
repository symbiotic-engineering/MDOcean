import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join("..", "mdocean"))) # for python autosummary to find the package
sys.path.insert(0, os.path.abspath('_ext')) # custom matlab_autosummary extension

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MDOcean'
copyright = '2025, Rebecca McCabe, Madison Dietrich, Olivia Murphy, Iris Ren, Maha Haji'
author = 'Rebecca McCabe, Madison Dietrich, Olivia Murphy, Iris Ren, Maha Haji'
release = 'v2.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
 #   'sphinxcontrib.bibtex',
    'sphinxcontrib.matlab',
    'sphinx.ext.autosummary', # must come before matlab_autosummary to avoid ref warnings
    'matlab_autosummary',
    'myst_parser'
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
primary_domain = "mat"
matlab_src_dir = os.path.abspath(os.path.join("..", "mdocean"))
matlab_auto_link = "all"
nitpicky = True


autosummary_generate = True
autosummary_output_dir = "generated"
autosummary_context = {
    'inputs': {'matmodules': ['inputs.validation', 'inputs.wave_conditions']},
    'optimization': {'matmodules': ['optimization.sensitivities', 'optimization.multiobjective']},
    'simulation': {'matmodules': ['simulation.modules', 'simulation.modules.dynamics', 'simulation.modules.econ']},
    'plots': {'matmodules': ['plots.util']}
}

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
