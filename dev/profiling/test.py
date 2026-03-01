import logging
import cProfile
import pstats
import time
from pathlib import Path

from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores
from datapackage import Package
from trails.logging import configure_trails_logging, trails_log_context
from trails.fair_rf import run_fair_delta_rf

configure_trails_logging(file_level=logging.DEBUG)

HERE = Path(__file__).resolve().parent

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

trails.import_excel_inventory(str(HERE.parent / "lci-case-study-daccs_storage_risk.xlsx"))

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

t0 = time.perf_counter()
trails.lca(
    methods=methods,
    show_progress=True,
    compute_score=True,
    store_inventory=True,
    solver_mode="iterative",
    iterative_rtol=1e-4,
)
print(f"lca_seconds={time.perf_counter() - t0:.3f}")

prof = cProfile.Profile()
t1 = time.perf_counter()
prof.enable()
rf = run_fair_delta_rf(
    trails,
    scenario="REMIND|SSP2-PkBudg650",
)
prof.disable()
rf_seconds = time.perf_counter() - t1
print(f"run_fair_delta_rf_seconds={rf_seconds:.3f}")
print(
    f"rf_shape={rf.shape} rf_nnz={int(getattr(rf.data, 'nnz', 0))} "
    f"rf_sum={float(rf.sum().item()):.12f}"
)

profile_out = Path(__file__).resolve().parent / "fair_rf_profile.txt"
with profile_out.open("w") as handle:
    stats = pstats.Stats(prof, stream=handle).sort_stats("cumtime")
    stats.print_stats(120)
    stats.print_callees("trails/fair_rf.py")
print(f"wrote profile to {profile_out}")
