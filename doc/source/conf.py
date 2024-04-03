# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os

project = 'casos'
copyright = '2024, Torbjørn Cunis, Jan Olucak, Renato Loureiro, Fabian Geyer'
author = 'Torbjørn Cunis, Jan Olucak, Renato Loureiro, Fabian Geyer'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))
matlab_src_dir = os.path.abspath("../..")

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.matlab',
    'sphinx.ext.napoleon',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

primary_domain = "mat"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
