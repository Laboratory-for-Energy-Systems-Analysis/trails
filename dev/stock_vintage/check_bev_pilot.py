"""Verify the real passenger-BEV exporter with the shared scoped-pilot matrix/LCI checks.

This retains legacy direct battery and maintenance roles in the comparison.
Full source inventories and derived temporal series remain restricted local data.
"""

import argparse
import json
from pathlib import Path

try:
    from .check_water_pilot import validate
except ImportError:
    from check_water_pilot import validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--variant", choices=["legacy", "corrected"], default="corrected"
    )
    parser.add_argument("--no-interpolation", action="store_true")
    parser.add_argument("--cache", action="store_true")
    parser.add_argument("--pulses-only", action="store_true")
    parser.add_argument("--lci-years", type=int, nargs="+", default=[2023])
    args = parser.parse_args()
    evidence = json.loads((args.directory / "export-evidence.json").read_text())
    if evidence.get("pilot") != "uk-passenger-bev":
        raise ValueError("Require the reviewed passenger-BEV exporter evidence")
    report = validate(
        args.directory,
        variant=args.variant,
        interpolate=not args.no_interpolation,
        cache=args.cache,
        lci_years=args.lci_years,
        pulses_only=args.pulses_only,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if args.variant == "corrected" and not all(
        r["static_conservation_passed"] for r in report["lci_checks"]
    ):
        raise AssertionError(
            f"BEV conservation failed; diagnostic report: {args.output}"
        )
    print(f"Recorded BEV pilot evidence: {args.output}")


if __name__ == "__main__":
    main()
