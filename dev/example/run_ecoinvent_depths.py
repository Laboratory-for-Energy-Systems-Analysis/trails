import os
from pathlib import Path

import logging
from datapackage import Package

from trails import Trails
from trails.logging import configure_trails_logging, trails_log_context


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    os.chdir(base_dir)

    # Match the example notebook defaults
    dp_path = Path("/Users/romain/GitHub/premise/dev/trails_2026-01-25.zip")
    methods = [
        "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - global warming potential (GWP100)",
    ]
    ref_year = 2050
    start_act_idx = 28673

    for depth in (0, 1, 2):
        log_name = f"trails_depth_{depth}.log"
        configure_trails_logging(file_level=logging.DEBUG, filename=log_name)

        with trails_log_context(run_id="ecoinvent_example", year=ref_year, depth=depth):
            dp = Package(str(dp_path))
            trails = Trails(dp, interpolate_annual=True, debug=True)

            trails.temporal_routing(
                start_year=ref_year,
                start_act_idx=start_act_idx,
                amount=1.0,
                max_depth=depth,
                show_progress=True,
                attribute_to_roots=True,
                debug=True,
            )

            trails.lca(
                methods=methods,
                show_progress=True,
                attribute_to_roots=True,
                compute_score=True,
                store_inventory=False,
                debug=True,
            )

            trails.static_lca(
                year=ref_year,
                act_idx=start_act_idx,
                methods=methods,
            )

            print(f"Depth {depth} complete. Log: {log_name}")


if __name__ == "__main__":
    main()
