import logging
from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores
from datapackage import Package
from trails.logging import configure_trails_logging, trails_log_context

configure_trails_logging(file_level=logging.DEBUG)

# dp = Package("datapackage.json")
# dp = Package("/Users/romain/GitHub/premise/dev/trails_2026-01-07.zip")
dp = Package("/Users/romain/GitHub/premise/dev/trails_2025-12-31 3.zip")
trails = Trails(dp)

methods = [
    "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - global warming potential (GWP100)",
]

ref_year = 2050

trails.lca(
    start_year=ref_year,
    start_act_idx=27827,
    methods=methods,
    max_depth=1,
    #max_depth=4,
    show_progress=True,
    min_amount=1e-22,
    attribute_to_roots=True,
    debug=False
)
