"""Validate the premise-generated lifecycle fixture through actual TRAILS LCI."""

import argparse
from collections import defaultdict
import json
from pathlib import Path

from datapackage import Package
import numpy as np

from trails import Trails


def identity(metadata):
    return tuple(
        metadata[key] for key in ("name", "reference product", "unit", "location")
    )


def validate(
    directory,
    *,
    solver_mode="direct",
    attribute_to_roots=True,
    interpolate_annual=True,
    adaptive_cutoff=None,
):
    directory = Path(directory)
    expected = json.loads((directory / "expectations.json").read_text())
    package = Package(str(directory / "synthetic_lifecycle.zip"))
    checks = []
    with Trails(
        package,
        stock_vintage=True,
        cache_interpolation=False,
        value_dtype=np.float64,
        interpolate_annual=interpolate_annual,
    ) as model:
        lookup = {
            identity(d): index
            for index, d in model.activity_indices[model.template_labels[0]].items()
        }
        caller = lookup[identity(expected["caller"])]
        capital = lookup[identity(expected["rewrite"]["capital_root"])]
        lifecycle = lookup[
            identity(expected["rewrite"]["lifted_exchanges"][0]["supplier"])
        ]
        for year in (2022, 2020, 2021, 2020):
            manufacture = next(
                r for r in expected["manufacture"] if r["service_year"] == year
            )
            retirement = next(
                r for r in expected["retirement"] if r["service_year"] == year
            )
            for amount in (1.0, -2.0):
                direct = model.expand_temporal_exchanges(
                    year=year, act_idx=caller, amount=amount
                )
                for supplier, record, coefficient in (
                    (capital, manufacture, 0.1),
                    (lifecycle, retirement, -0.1),
                ):
                    actual = {
                        y: values[supplier]
                        for y, values in direct.items()
                        if supplier in values
                    }
                    target = dict(
                        zip(
                            record["event_years"],
                            [amount * coefficient * w for w in record["weights"]],
                        )
                    )
                    assert actual.keys() == target.keys()
                    assert all(
                        np.isclose(actual[y], target[y], atol=1e-12, rtol=1e-12)
                        for y in actual
                    )
                routing = {
                    "max_depth": 4,
                    "adaptive_relative_score_cutoff": None,
                }
                if adaptive_cutoff is not None:
                    routing = {
                        "max_depth": None,
                        "adaptive_relative_score_cutoff": adaptive_cutoff,
                        "adaptive_methods": [
                            "CML v4.8 2016 no LT - climate change no LT - "
                            "global warming potential (GWP100) no LT"
                        ],
                    }
                model.temporal_routing(
                    start_year=year,
                    start_act_idx=caller,
                    amount=amount,
                    **routing,
                    min_amount=0,
                    show_progress=False,
                    attribute_to_roots=attribute_to_roots,
                )
                actual_capital_years = {n[0] for n in model.graph if n[-1] == capital}
                assert actual_capital_years == set(manufacture["event_years"])
                assert {n[0] for n in model.graph if n[-1] == lifecycle} == set(
                    retirement["event_years"]
                )
                assert min(retirement["event_years"]) > year
                inventory = model.lci(show_progress=False, solver_mode=solver_mode)
                total = float(inventory.sum().item())
                assert np.isclose(
                    total, amount * expected["co2_per_service"], atol=1e-12, rtol=1e-12
                )
                # With constant technology, the exact biosphere calendar series
                # has operating CO2 now, manufacture CO2 at cohort years and a
                # signed recycling credit after service. This detects double
                # shifts inside wrappers, not just direct capital-port timing.
                target_bio = defaultdict(float, {year: amount * 0.3})
                for event, weight in zip(
                    manufacture["event_years"], manufacture["weights"]
                ):
                    target_bio[event] += amount * 0.2 * weight
                for event, weight in zip(
                    retirement["event_years"], retirement["weights"]
                ):
                    target_bio[event] -= amount * 0.05 * weight
                annual = inventory.sum(dim=[d for d in inventory.dims if d != "year"])
                values = (
                    annual.data.todense()
                    if hasattr(annual.data, "todense")
                    else np.asarray(annual.data)
                )
                actual_bio = {
                    int(y): float(v)
                    for y, v in zip(annual.coords["year"].values, values)
                    if abs(v) > 1e-14
                }
                assert actual_bio.keys() == target_bio.keys()
                assert all(
                    np.isclose(actual_bio[y], target_bio[y], atol=1e-12, rtol=1e-12)
                    for y in actual_bio
                )
                checks.append(
                    {
                        "service_year": year,
                        "functional_unit": amount,
                        "biosphere_by_year": actual_bio,
                        "summed_co2": total,
                        "all_retirements_after_service": True,
                    }
                )
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    checks = validate(args.directory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({"evidence_tier": "synthetic", "checks": checks}, indent=2) + "\n"
    )
    print(
        f"Validated {len(checks)} manufacture/retirement routing and signed LCI calendars: {args.output}"
    )


if __name__ == "__main__":
    main()
