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
    'myst_parser',
    'sphinx_design',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    #'sphinx_last_updated_by_git'
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
    'simulation': {'matmodules': ['simulation.modules', 'simulation.run']},
    'simulation.modules': {'matmodules': ['simulation.modules.dynamics', 'simulation.modules.econ']},
    'plots': {'matmodules': ['plots.util']}
}

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}

globaltoc_maxdepth = 2

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/SEALab_Logo_Light_202101_120ht.png'

html_context = {
    "display_github": True,
    "github_user": "symbiotic-engineering",
    "github_repo": "MDOcean",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

html_theme_options = {
    'version_selector': True,
}