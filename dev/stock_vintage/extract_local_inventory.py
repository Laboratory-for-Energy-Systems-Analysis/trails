"""Read an existing Brightway inventory into a restricted local export input.

Run in the source Brightway environment. No database is modified. The compressed
JSON is licensed inventory data: keep it outside repositories and distribution
packages. It is an input to the companion premise pilot exporter, not a public
test fixture. Descriptive comments are deliberately omitted.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--database", default="ecoinvent-3.12-cutoff")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    import bw2data as bd

    if args.project not in {p.name for p in bd.projects}:
        raise ValueError("An existing Brightway project is required")
    previous = bd.projects.current
    try:
        bd.projects.set_current(args.project)
        if args.database not in bd.databases:
            raise ValueError("Existing source database is required")
        source = bd.Database(args.database).load()
        bio_databases = {
            e["input"][0]
            for d in source.values()
            for e in d["exchanges"]
            if e["type"] == "biosphere"
        }
        suppliers = dict(source)
        for name in sorted(bio_databases):
            suppliers.update(bd.Database(name).load())
        result = []
        fields = ("name", "reference product", "unit", "location")
        seen = set()
        for key, dataset in sorted(source.items()):
            identity = tuple(dataset[f] for f in fields)
            if identity in seen:
                raise ValueError(f"Ambiguous exporter identity: {identity}")
            seen.add(identity)
            record = {f: dataset[f] for f in fields}
            record.update(database=key[0], code=key[1], exchanges=[])
            for original in dataset["exchanges"]:
                kind, supplier_key = original["type"], tuple(original["input"])
                if kind not in {"production", "technosphere", "biosphere"}:
                    raise ValueError(f"Unsupported exchange type: {kind}")
                if kind != "biosphere" and supplier_key not in source:
                    raise ValueError("The selected inventory is not closed")
                supplier = suppliers[supplier_key]
                amount = float(original["amount"])
                if not math.isfinite(amount):
                    raise ValueError("Nonfinite inventory coefficient")
                exchange = dict(
                    name=supplier["name"],
                    unit=supplier["unit"],
                    input=supplier_key,
                    output=key,
                    type=kind,
                    amount=amount,
                )
                if kind == "biosphere":
                    exchange["categories"] = supplier["categories"]
                else:
                    exchange.update(
                        product=supplier["reference product"],
                        location=supplier["location"],
                    )
                for field in (
                    "uncertainty type",
                    "loc",
                    "scale",
                    "shape",
                    "minimum",
                    "maximum",
                    "negative",
                ):
                    if field in original:
                        value = original[field]
                        exchange[field] = (
                            None
                            if isinstance(value, float) and not math.isfinite(value)
                            else value
                        )
                record["exchanges"].append(exchange)
            result.append(record)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(args.output, "wt", encoding="utf-8") as stream:
            json.dump(
                {
                    "database": args.database,
                    "project": args.project,
                    "source_version": "3.12",
                    "system_model": "cutoff",
                    "coefficient_precision": "source metadata, not processed float32 matrix",
                    "rights": "restricted local ecoinvent inventory; do not redistribute",
                    "inventory": result,
                },
                stream,
                allow_nan=False,
            )
        manifest = {
            "file": args.output.name,
            "sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
            "activities": len(result),
            "exchanges": sum(len(d["exchanges"]) for d in result),
        }
        args.output.with_suffix(args.output.suffix + ".manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        print(json.dumps(manifest))
    finally:
        bd.projects.set_current(previous)


if __name__ == "__main__":
    main()
