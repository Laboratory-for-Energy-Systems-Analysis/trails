"""
TRAILS: Temporal Routing And Aggregation of Impacts across Life-cycle Systems.
"""

__all__ = (
    "Trails",
    "lca",
    "get_lcia_method_names",
    "plot_temporal_scores",
    "clear_cache",
    "search_activity",
)

__version__ = "1.0.0"

from .trails import Trails
from .lca import lca
from .lcia import get_lcia_method_names
from .plotting import plot_temporal_scores
from .cache import clear_cache
from .search import search_activity
