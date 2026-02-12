import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "trails"
copyright = "2025"
author = "Paul Scherrer Institute"

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode"]
templates_path = ["_templates"]
exclude_patterns = []

# docs/conf.py
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "bw2calc",
    "bw2data",
    "bw2io",
    "bw_processing",
    "bw2analyzer",
    "bw2parameters",
    "bw_migrations",
    "pypardiso",
    "scikits",
    "scikits.umfpack",
    "fair",
    "plotly",
    "pyvis",
    "networkx",
    "datapackage",
    "sparse",
    "xarray",
    "prettytable",
    "tqdm",
    "constructive_geometries",
]


html_theme = "alabaster"
html_static_path = ["_static"]
html_logo = "../assets/permanent/trails_logo_grey_on_white.png"

import os
import sys

sys.path.insert(0, os.path.abspath("../"))  # or '../src' if your code is in src/
