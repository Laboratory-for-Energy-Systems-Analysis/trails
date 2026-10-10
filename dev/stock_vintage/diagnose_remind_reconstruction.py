"""Summarise constraints and sensitivity in the archived REMIND reconstructions.

This reads eight existing cohort files, verifies their archive hashes and
compares stock-target with reported-sales/additions cases. No IAM or LCI is
rerun. Keep output in ignored local research storage.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def percent_difference(value, reference):
    return 100 * (value - reference) / reference if reference else None


def summarise(path):
    data = json.loads(path.read_text())
    primary = data["cases"]["primary"]
    final = primary["annual"][-1]
    vehicle = "total_stock_vehicles" in final
    alternate_name = "reported_sales" if vehicle else "reported_additions"
    alternate = data["cases"][alternate_name]["annual"][-1]
    reported_field = (
        "iam_scaled_sales_vehicles" if vehicle
        else "iam_scaled_reported_additions_mw_ac"
    )
    stock_field = "total_stock_vehicles" if vehicle else "total_capacity_mw_ac"
    weighting = "declining_utilisation" if vehicle else "capacity_weights"
    assert alternate["service_year"] == final["service_year"]
    reported = math.fsum(b[reported_field] for b in primary["balances"])
    inferred = math.fsum(b["gross_additions"] for b in primary["balances"])
    ages = [case["annual"][-1]["mean_age_years"]
            for case in data["cases"].values()]
    return {
        "file": path.name,
        "scenario": data["iam_source"]["scenario"],
        "reference_year": data["reference_year"],
        "final_year": final["service_year"],
        "stock_unit": "vehicle" if vehicle else "MW AC",
        "inferred_additions_sum": inferred,
        "reported_additions_sum": reported,
        "addition_difference_absolute": inferred - reported,
        "addition_difference_pct": percent_difference(inferred, reported),
        "primary_final_stock": final[stock_field],
        "reported_additions_final_stock": alternate[stock_field],
        "reported_case_stock_difference_pct": percent_difference(
            alternate[stock_field], final[stock_field]
        ),
        "primary_mean_age": final["mean_age_years"],
        "reported_case_mean_age_difference": (
            alternate["mean_age_years"] - final["mean_age_years"]
        ),
        "weighting_case": weighting,
        "weighting_case_mean_age_difference": (
            data["cases"][weighting]["annual"][-1]["mean_age_years"]
            - final["mean_age_years"]
        ),
        "all_cases_min_age": min(ages),
        "all_cases_max_age": max(ages),
        "primary_early_exits_sum": math.fsum(
            b["early_exits"] for b in primary["balances"]
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=Path.home() / (
        "Library/Application Support/trails/stock-vintage-evidence/"
        "2026-10-10-pilots"
    ))
    parser.add_argument("--output-dir", type=Path, default=(
        Path(__file__).resolve().parents[2] / "build/reports/remind-robustness"
    ))
    args = parser.parse_args()
    manifest_path = args.archive / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    files = {r["path"]: r for r in manifest["files"]}
    used, rows = {}, []
    for path in sorted((args.archive / "cohorts").glob("*.json")):
        if path.name == "water-cohorts.json":
            continue
        relative = str(path.relative_to(args.archive))
        actual = digest(path)
        assert actual == files[relative]["sha256"]
        assert path.stat().st_size == files[relative]["bytes"]
        used[relative] = actual
        rows.append(summarise(path))
    assert len(rows) == 8
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "reconstruction-diagnostics.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "diagnostics-provenance.json").write_text(json.dumps({
        "archive": str(args.archive),
        "manifest_sha256": digest(manifest_path),
        "script_sha256": digest(Path(__file__)),
        "output_sha256": digest(output),
        "sources": used,
        "definitions": {
            "addition_difference": "Primary inferred minus reported, 2023–2030",
            "stock_difference": "Reported-additions case minus primary at 2030",
            "percentages": "Addition denominator: reported sum; stock: primary",
            "interpretation": (
                "Differences between labelled reconstructions, not errors proven "
                "in REMIND. Mean ages alone cannot validate quantity consistency."
            ),
        },
        "rights": "Restricted local IAM-derived evidence; keep outside Git",
    }, indent=2) + "\n")
    print(json.dumps({"output": str(output), "verified_cohort_files": len(rows)}))


if __name__ == "__main__":
    main()
