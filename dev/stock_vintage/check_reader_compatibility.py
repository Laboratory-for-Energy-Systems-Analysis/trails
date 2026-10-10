"""Probe an actual TRAILS revision with the synthetic premise stock package.

Run once with the preserved baseline on PYTHONPATH and once with the corrected
reader. This demonstrates behaviour; metadata alone cannot protect old readers.
"""

import argparse
import hashlib
import json
from pathlib import Path

from datapackage import Package
import numpy as np

import trails
from trails import Trails


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("--reader-revision", required=True)
    parser.add_argument("--expect", choices=["legacy", "corrected"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    corrected = args.expect == "corrected"
    options = {"stock_vintage": True} if corrected else {}
    with Trails(
        Package(str(args.package)),
        cache_interpolation=False,
        value_dtype=np.float64,
        **options,
    ) as model:
        indices = model.activity_indices[model.template_labels[0]]
        caller = next(
            i for i, d in indices.items() if d["name"] == "synthetic fleet service"
        )
        supplier = next(
            i for i, d in indices.items() if d["name"] == "synthetic new asset"
        )
        demand = model.expand_temporal_exchanges(year=2021, act_idx=caller, amount=1)
        pulses = {
            y: values[supplier] for y, values in demand.items() if supplier in values
        }
        target = {2010: 0.1, 2021: 0.15}
        respected = pulses.keys() == target.keys() and all(
            np.isclose(pulses[y], target[y], rtol=1e-12, atol=1e-12) for y in target
        )
        if respected != corrected:
            raise AssertionError(
                "Reader behaviour differs from compatibility expectation"
            )
        result = {
            "reader_revision": args.reader_revision,
            "reader_file": trails.__file__,
            "package_sha256": hashlib.sha256(args.package.read_bytes()).hexdigest(),
            "loaded": True,
            "stock_flag_supplied": corrected,
            "annual_resource_respected": bool(respected),
            "service_year": 2021,
            "actual_capital_pulses": pulses,
            "expected_capital_pulses": target,
            "conclusion": (
                "Corrected opt-in reader uses actual-year resource"
                if corrected
                else "Unsupported old reader loads package but ignores required annual profile"
            ),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
