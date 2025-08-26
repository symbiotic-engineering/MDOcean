import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join("..", "mdocean")))
sys.path.insert(0, os.path.abspath('_ext'))

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
    'sphinx.ext.autosummary',
    'matlab_autosummary'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
primary_domain = "mat"
matlab_src_dir = "../mdocean"
matlab_auto_link = "all"
nitpicky = True


autosummary_generate = True

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
