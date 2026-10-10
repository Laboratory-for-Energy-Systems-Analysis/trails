"""Verify actual REMIND cohort profiles through premise export and TRAILS.

The inventory is deliberately synthetic. This validates the new data route,
signs, calendar years and varying-background lookup; it is not a regional LCA.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

from datapackage import Package
import numpy as np

from trails import Trails


def identity(record):
    return tuple(record[k] for k in ("name", "reference product", "unit", "location"))


def verify(directory):
    source = json.loads((directory / "fixture-evidence.json").read_text())
    package = Path(source["package"])
    if hashlib.sha256(package.read_bytes()).hexdigest() != source["sha256"]:
        raise ValueError("Fixture package hash mismatch")
    checks = []
    with Trails(
        Package(str(package)),
        stock_vintage=True,
        cache_interpolation=False,
        interpolate_annual=True,
        value_dtype=np.float64,
    ) as model:
        lookup = {
            identity(metadata): index
            for index, metadata in model.activity_indices[
                model.template_labels[0]
            ].items()
        }
        profiles = {
            p["id"]: {r["service_year"]: r for r in p["years"]}
            for p in source["profiles"]
        }
        for asset in source["assets"]:
            caller, supplier = (
                lookup[identity(asset["caller"])],
                lookup[identity(asset["supplier"])],
            )
            for year in (2050, 2023, 2030, 2023):
                record = profiles[asset["asset"]][year]
                for amount in (1.0, -2.0):
                    pulses = model.expand_temporal_exchanges(
                        year=year, act_idx=caller, amount=amount
                    )
                    actual = {
                        y: row[supplier]
                        for y, row in pulses.items()
                        if supplier in row and row[supplier] != 0
                    }
                    expected = {
                        y: amount * asset["coefficient"] * w
                        for y, w in zip(record["event_years"], record["weights"])
                    }
                    if actual.keys() != expected.keys() or any(
                        not math.isclose(
                            actual[y], expected[y], rel_tol=1e-10, abs_tol=1e-30
                        )
                        for y in expected
                    ):
                        raise AssertionError(
                            "TRAILS changed the prepared cohort calendar or weights"
                        )
                    quantity = math.fsum(actual.values())
                    expected_quantity = amount * asset["coefficient"]
                    if not math.isclose(
                        quantity, expected_quantity, rel_tol=1e-10, abs_tol=1e-30
                    ):
                        raise AssertionError(
                            "Signed inventory coefficient was not conserved"
                        )
                    # Independent analytical sum: factor(year) is linearly interpolated
                    # between fixture anchors and clamped only for technology lookup.
                    expected_lci = math.fsum(
                        quantity
                        * asset["biosphere_factor"]
                        * float(np.interp(y, source["anchors"], source["multipliers"]))
                        for y, quantity in expected.items()
                    )
                    model.temporal_routing(
                        start_year=year,
                        start_act_idx=caller,
                        amount=amount,
                        max_depth=1,
                        adaptive_relative_score_cutoff=None,
                        min_amount=0.0,
                        show_progress=False,
                    )
                    graph_years = {
                        node[0] for node in model.graph if node[-1] == supplier
                    }
                    if graph_years != expected.keys():
                        raise AssertionError(
                            "Physical event dates were clamped or changed"
                        )
                    inventory = model.lci(show_progress=False, solver_mode="direct")
                    value = inventory.sum().data
                    actual_lci = float(
                        value.todense() if hasattr(value, "todense") else value
                    )
                    if not math.isclose(
                        actual_lci, expected_lci, rel_tol=1e-9, abs_tol=1e-15
                    ):
                        raise AssertionError(
                            f"Dated LCI mismatch: {actual_lci} vs {expected_lci}"
                        )
                    checks.append(
                        {
                            "asset": asset["asset"],
                            "year": year,
                            "amount": amount,
                            "cohorts": len(actual),
                            "quantity": quantity,
                            "expected_quantity": expected_quantity,
                            "lci": actual_lci,
                            "expected_lci": expected_lci,
                            "passed": True,
                        }
                    )
        try:
            model.expand_temporal_exchanges(year=2051, act_idx=caller, amount=1)
        except ValueError as exc:
            if "service year" not in str(exc):
                raise
        else:
            raise AssertionError("Missing cohort year did not fail closed")
    result = {
        "schema_version": 1,
        "scope": source["kind"],
        "package_sha256": source["sha256"],
        "varying_background": source["varying_background"],
        "checks": checks,
        "missing_year_rejected": True,
        "all_passed": True,
    }
    (directory / "trails-validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(
        json.dumps(
            {"directory": str(directory), "checks": len(checks), "all_passed": True}
        )
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("integration_directory", type=Path)
    args = parser.parse_args()
    for variant in ("constant", "varying"):
        verify(args.integration_directory / variant)
