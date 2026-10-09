"""Check the synthetic package produced by premise/tests/test_stock_vintage.py.

This is a cross-environment integration check: premise and TRAILS currently
require different Python versions. It does not use an empirical stock profile.
"""

import argparse
import json
from pathlib import Path

from datapackage import Package
import numpy as np

from trails import Trails


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    package = Package(str(args.package))
    results = []
    with Trails(
        package, stock_vintage=True, cache_interpolation=False, value_dtype=np.float64
    ) as model:
        indices = model.activity_indices[model.template_labels[0]]
        caller = next(
            i for i, d in indices.items() if d["name"] == "synthetic fleet service"
        )
        supplier = next(
            i for i, d in indices.items() if d["name"] == "synthetic new asset"
        )
        expected = {
            2020: {2010: 0.75, 2015: 0.25},
            2021: {2010: 0.4, 2021: 0.6},
            2022: {2010: 0.2, 2021: 0.8},
        }
        for year in (2022, 2020, 2021, 2020):
            for amount in (1.0, -2.0):
                demand = model.expand_temporal_exchanges(
                    year=year, act_idx=caller, amount=amount
                )
                actual = {
                    y: values[supplier]
                    for y, values in demand.items()
                    if supplier in values
                }
                target = {
                    y: weight * amount * 0.25 for y, weight in expected[year].items()
                }
                assert actual.keys() == target.keys()
                assert all(
                    np.isclose(actual[y], target[y], rtol=1e-12, atol=1e-12)
                    for y in target
                )
                model.temporal_routing(
                    start_year=year,
                    start_act_idx=caller,
                    amount=amount,
                    max_depth=1,
                    adaptive_relative_score_cutoff=None,
                    min_amount=0,
                    show_progress=False,
                )
                inventory = model.lci(show_progress=False, solver_mode="direct")
                total = float(inventory.sum().item())
                assert np.isclose(total, 0.5 * amount, rtol=1e-12, atol=1e-12)
                results.append(
                    {
                        "service_year": year,
                        "functional_unit": amount,
                        "capital_pulses": actual,
                        "co2_total": total,
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({"evidence_tier": "synthetic", "checks": results}, indent=2) + "\n"
    )
    print(f"Validated {len(results)} signed annual routing/LCI cases: {args.output}")


if __name__ == "__main__":
    main()
