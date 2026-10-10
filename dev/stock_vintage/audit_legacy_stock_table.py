"""Freeze the complete legacy stock table and reproduce its consumer weights.

This is a baseline diagnostic, not scientific approval or parameter migration.
Unscoped exchanges retain their legacy status; only explicit pilot bindings are
revised. The asymmetric triangular CDF comparison is an independent diagnostic.
"""

import argparse
from collections import Counter
import csv
import hashlib
import inspect
import json
import math
from pathlib import Path

from trails.datapackage import _parse_temporal_exchange_row
from trails.temporal_distributions import TemporalDistribution


def triangular_cdf(x, lower, mode, upper):
    if x <= lower:
        return 0.0
    if x >= upper:
        return 1.0
    if x <= mode and mode > lower:
        return (x - lower) ** 2 / ((upper - lower) * (mode - lower))
    return 1 - (upper - x) ** 2 / ((upper - lower) * (upper - mode))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with args.table.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    identities = Counter(
        (r["name"], r["reference product"])
        for r in source
        if r["temporal_tag"] == "stock_asset"
    )
    rows, profiles = [], []
    for index, row in enumerate(source, start=2):
        if row["temporal_tag"] != "stock_asset":
            continue
        fields = {
            target: float(row[key]) if row[key].strip() else None
            for target, key in {
                "temporal_distribution": "age distribution type",
                "temporal_loc": "loc",
                "temporal_scale": "scale",
                "temporal_min": "minimum",
                "temporal_max": "maximum",
            }.items()
        }
        exchange = _parse_temporal_exchange_row(fields)
        if exchange is None:
            raise ValueError(f"Legacy stock row does not parse: {index}")
        pulses = list(TemporalDistribution(exchange).iter_offsets_and_weights())
        offsets, weights = zip(*pulses)
        if any(not math.isfinite(w) or w < 0 for w in weights):
            raise ValueError(f"Invalid legacy weights: {index}")
        if not math.isclose(math.fsum(weights), 1, abs_tol=1e-12):
            raise ValueError(f"Legacy mass does not close: {index}")
        uniform = all(math.isclose(w, 1 / len(weights), abs_tol=1e-12) for w in weights)
        lower, mode, upper = (
            fields[k] for k in ("temporal_min", "temporal_loc", "temporal_max")
        )
        reference_mean, reference_mass, continuous_mean = None, None, None
        if exchange.distribution == 5 and lower < upper and lower <= mode <= upper:
            independent = [
                triangular_cdf(o + 0.5, lower, mode, upper)
                - triangular_cdf(o - 0.5, lower, mode, upper)
                for o in offsets
            ]
            mass = math.fsum(independent)
            reference_mass = mass
            continuous_mean = -(lower + mode + upper) / 3
            reference_mean = (
                -math.fsum(o * w for o, w in zip(offsets, independent)) / mass
            )
        record = {
            "csv_row": index,
            "name": row["name"],
            "reference_product": row["reference product"],
            "identity_occurrences": identities[(row["name"], row["reference product"])],
            "distribution": exchange.distribution,
            "source_minimum": lower,
            "source_maximum": upper,
            "source_loc": mode,
            "consumed_minimum": min(offsets),
            "consumed_maximum": max(offsets),
            "consumed_mean_age": -math.fsum(o * w for o, w in pulses),
            "uniform_weights": uniform,
            "asymmetric_triangle_binned_mean_age": reference_mean,
            "asymmetric_triangle_retained_mass": reference_mass,
            "asymmetric_triangle_continuous_mean_age": continuous_mean,
            "disposition": "Legacy retained outside exact opt-in bindings; no blanket scientific approval",
        }
        rows.append(record)
        profiles.append({"csv_row": index, "offsets": offsets, "weights": weights})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "legacy-stock-row-audit.csv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    runtime = Path(inspect.getfile(TemporalDistribution))
    report = {
        "source_table": str(args.table.resolve()),
        "source_sha256": hashlib.sha256(args.table.read_bytes()).hexdigest(),
        "runtime_file": str(runtime),
        "runtime_sha256": hashlib.sha256(runtime.read_bytes()).hexdigest(),
        "stock_rows": len(rows),
        "unique_identities": len(identities),
        "duplicate_identities": sum(n > 1 for n in identities.values()),
        "distribution_counts": dict(Counter(r["distribution"] for r in rows)),
        "lognormal_uniform_count": sum(
            r["distribution"] == 2 and r["uniform_weights"] for r in rows
        ),
        "lognormal_maximum_zero_count": sum(
            r["distribution"] == 2 and r["source_maximum"] == 0 for r in rows
        ),
        "triangle_cdf_comparisons": sum(
            r["asymmetric_triangle_binned_mean_age"] is not None for r in rows
        ),
        "row_audit_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
        "interpretation": "Complete legacy numerical/disposition audit, not full-database exchange coverage or approval of the other asset groups",
        "profiles": profiles,
    }
    (args.output_dir / "legacy-stock-baseline.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({k: v for k, v in report.items() if k != "profiles"}))


if __name__ == "__main__":
    main()
