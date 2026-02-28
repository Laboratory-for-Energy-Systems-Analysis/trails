import logging
from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores
from datapackage import Package
from trails.logging import configure_trails_logging, trails_log_context

configure_trails_logging(file_level=logging.DEBUG)

dp = Package("/Users/romain/GitHub/premise/dev/trails_2026-02-22.zip")
trails = Trails(
    dp,
    interpolate_annual=True,
    debug=False,
    interpolation_start_year_offset=-20,
    interpolation_end_year_offset=20
)
methods = [
    "IPCC 2021 - climate change: total (excl. biogenic CO2) - global warming potential (GWP100)",
]

trails.import_excel_inventory("../lci-case-study-daccs_storage_risk.xlsx")

ref_year = 2025
idx = 41792
trails.temporal_routing(
    start_year=ref_year,
    start_act_idx=idx,
    amount=1,
    max_depth=4,
    show_progress=True,
    attribute_to_roots=True,
)

trails.lca(
    methods=methods,
    show_progress=True,
    compute_score=True,
    store_inventory=True,
    solver_mode="iterative",
    iterative_rtol=1e-4,
)
