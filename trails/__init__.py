"""
TRAILS: Temporal Routing And Aggregation of Impacts across Life-cycle Systems.
"""

__all__ = (
    "Trails",
    "lca",
    "lci",
    "lcia",
    "get_lcia_method_names",
    "get_edges_lcia_method_names",
    "plot_temporal_scores",
    "plot_rf",
    "plot_temp",
    "plot_adaptive_sankey",
    "clear_cache",
    "search_activity",
)

__version__ = "1.0.1"

from .trails import Trails
from .lca import lca, lci
from .lcia import get_lcia_method_names, lcia
from .edges_matrix import get_edges_lcia_method_names
from .plotting import (
    plot_temporal_scores,
    plot_rf,
    plot_temp,
    plot_adaptive_sankey,
)
from .cache import clear_cache, _warn_if_cache_large
from .search import search_activity

_warn_if_cache_large()
