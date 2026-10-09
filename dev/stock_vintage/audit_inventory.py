"""Inspect candidate pilot activities without exporting licensed inventories.

Run in an environment with the local Brightway project. Outputs contain activity
identity and selected capital/maintenance exchange metadata only; no biosphere
or complete inventory datasets are exported. This script does not modify a DB.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import re

CANDIDATES = {
    "passenger_cars": r"^transport, passenger,? car,.*(medium|diesel|petrol|electric)",
    "heavy_trucks": r"^transport, freight, lorry",
    "pv": r"^electricity production, photovoltaic",
    "gas_power": r"^electricity production, natural gas, combined cycle",
    "water": r"^(tap water production|water supply network|market for tap water)",
}
CAPITAL = re.compile(
    r"passenger car|lorry|photovoltaic|power plant|water supply network|"
    r"pipe|maintenance|treatment of|construction|installation|building|pump|"
    r"inverter|battery|dismantl|disposal|renewal",
    re.IGNORECASE,
)


def identity(activity) -> dict:
    return {
        "database": activity.key[0],
        "code": activity.key[1],
        "name": activity.get("name"),
        "reference_product": activity.get("reference product"),
        "location": activity.get("location"),
        "unit": activity.get("unit"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--database", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--locations", nargs="+", default=["GB", "RER", "CA", "CA-QC"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inspect-codes", nargs="*", default=[])
    args = parser.parse_args()

    import bw2data as bd

    existing_projects = {project.name for project in bd.projects}
    if args.project not in existing_projects:
        raise ValueError(f"Existing project required: {args.project!r}")
    previous = bd.projects.current
    try:
        bd.projects.set_current(args.project)
        if args.database not in bd.databases:
            raise ValueError(f"Database unavailable: {args.database!r}")
        patterns = {key: re.compile(value, re.I) for key, value in CANDIDATES.items()}
        candidates = {key: [] for key in patterns}
        inspections = []
        requested = set(args.inspect_codes)
        found = set()
        for activity in bd.Database(args.database):
            name = activity.get("name", "")
            location = activity.get("location")
            if "*" in args.locations or location in args.locations:
                for key, pattern in patterns.items():
                    if pattern.search(name):
                        candidates[key].append(identity(activity))
            if activity.key[1] not in requested:
                continue
            found.add(activity.key[1])
            exchanges = []
            for exchange in activity.technosphere():
                supplier = exchange.input
                if not CAPITAL.search(supplier.get("name", "")):
                    continue
                row = identity(supplier)
                row["amount"] = float(exchange["amount"])
                # Comments can contain complete licensed inventory descriptions;
                # inspect those in the source DB, not in this metadata artifact.
                row["exchange_comment_present"] = bool(exchange.get("comment"))
                exchanges.append(row)
            inspections.append(
                {
                    "activity": identity(activity),
                    "selected_capital_exchanges": exchanges,
                    "total_technosphere_exchanges": len(list(activity.technosphere())),
                }
            )
        if requested - found:
            raise ValueError(f"Unknown activity codes: {sorted(requested - found)}")
        for rows in candidates.values():
            rows.sort(key=lambda row: (row["name"], row["location"]))
        report = {
            "scope": "candidate identities and selected capital exchange metadata",
            "project": args.project,
            "database": args.database,
            "bw2data_version": importlib.metadata.version("bw2data"),
            "locations": args.locations,
            "candidates": candidates,
            "inspections": inspections,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(
            json.dumps(
                {
                    "candidates": {k: len(v) for k, v in candidates.items()},
                    "inspected": len(inspections),
                    "output": str(args.output),
                }
            )
        )
    finally:
        bd.projects.set_current(previous)


if __name__ == "__main__":
    main()
