import logging
from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores
from datapackage import Package
from trails.logging import configure_trails_logging, trails_log_context

configure_trails_logging(file_level=logging.DEBUG)

dp = Package("/Users/romain/GitHub/premise/dev/trails_2026-01-25.zip")
trails = Trails(dp, interpolate_annual=True, debug=False)
methods = [
    "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - global warming potential (GWP100)",
]


trails.temporal_routing(
    start_year=2050,
    start_act_idx=28673,
    max_depth=4,
    show_progress=True,
    attribute_to_roots=True,
)
