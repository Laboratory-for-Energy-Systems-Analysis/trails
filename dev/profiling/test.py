import logging

from datapackage import Package

from trails import Trails, lca
from trails.logging import configure_trails_logging

def main() -> None:
    configure_trails_logging(file_level=logging.DEBUG)

    # dp = Package("../example/datapackage.json")
    dp = Package("/Users/romain/GitHub/premise/dev/trails_2025-12-16.zip")
    trails = Trails(dp, interpolate_annual=True)

    methods = [
        "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - global warming potential (GWP100)",
    ]

    ref_year = 2050

    results_by_year_ICEV, provenance = lca(
        trails,
        start_year=ref_year,
        start_act_idx=13,
        methods=methods,
        max_depth=4,
        return_provenance=True,
        min_amount=1e-18,
    )


if __name__ == "__main__":
    main()
