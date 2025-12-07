"""
edges: A Python package for calculating the environmental impact of products by
applying characterization factors conditioned by the context of exchanges.
"""

__all__ = (
    "A3",
    "lca",
    "get_lcia_method_names",
    "plot_temporal_scores"
)

__version__ = "1.0.0"

from .a3 import A3
from .lca import lca
from .lcia import get_lcia_method_names
from .plotting import plot_temporal_scores
