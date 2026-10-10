"""Compare raw and actually routed parent manufacture calendars for five pilots.

This graph-only diagnostic separates the legacy distribution from its inventory
year clamping. It does not replace full LCI or validate lifetime assumptions.
Keep restricted package-derived output outside version control.
"""

import argparse
from collections import defaultdict
import csv
import gc
import hashlib
import json
import math
import os
from pathlib import Path

from datapackage import Package
import numpy as np

from trails import Trails

PILOTS = {
    "Passenger BEV chassis": "bev-real-primary",
    "Heavy truck": "trucks-real-primary",
    "PV parent": "pv-real-primary",
    "Gas CCGT": "ccgt-real-npi",
    "Water pipe/network": "water-real-primary",
}


def identity(row):
    return tuple(row[k] for k in ("name", "reference product", "unit", "location"))


def mean_age(pulses, year):
    return math.fsum((year - y) * v for y, v in pulses.items()) / math.fsum(
        pulses.values()
    )


def compare(directory, year):
    evidence = json.loads((directory / "export-evidence.json").read_text())
    audit, reports = evidence["rewrite"], {}
    for variant in ("legacy", "corrected"):
        resource = evidence["packages"][variant]
        path = Path(resource["path"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != resource["sha256"]:
            raise ValueError("Package differs from recorded export")
        corrected = variant == "corrected"
        with Trails(
            Package(str(path)),
            stock_vintage=corrected,
            cache_interpolation=False,
            interpolate_annual=False,
            value_dtype=np.float64,
        ) as model:
            indices = model.activity_indices[model.template_labels[0]]
            lookup = {identity(d): i for i, d in indices.items()}
            caller = lookup[
                identity(audit["caller" if corrected else "original_caller"])
            ]
            capital = lookup[
                identity(audit["capital_root" if corrected else "original_root"])
            ]
            constructor = dict(
                next(
                    row["supplier"]
                    for row in audit["zero_shift_bindings"]
                    if identity(row["caller"]) == identity(audit["capital_root"])
                )
            )
            if not corrected:
                constructor["name"] = constructor["name"].removesuffix(
                    f" [stock context {audit['context_id']}: capital]"
                )
            construction = lookup[identity(constructor)]
            net = float(model.A[0, caller, caller])
            target = audit["preserved_capital_amount"] / net
            raw = defaultdict(float)
            direct = model.expand_temporal_exchanges(
                year=year, act_idx=caller, amount=1 / net
            )
            for date, inputs in direct.items():
                if not inputs.get(capital, 0):
                    continue
                for event, nested in model.expand_temporal_exchanges(
                    year=date, act_idx=capital, amount=inputs[capital]
                ).items():
                    if nested.get(construction, 0):
                        raw[event] += nested[construction]
            model.temporal_routing(
                start_year=year,
                start_act_idx=caller,
                amount=1,
                max_depth=2,
                adaptive_relative_score_cutoff=None,
                min_amount=0,
                show_progress=False,
                attribute_to_roots=True,
            )
            routed = defaultdict(float)
            for node, data in model.graph.nodes(data=True):
                if node[-1] == construction:
                    routed[node[0]] += data["amount"]
            for values in (raw, routed):
                if not math.isclose(math.fsum(values.values()), target, rel_tol=1e-10):
                    raise AssertionError("Parent calendar changed its signed total")
            if corrected and (
                raw.keys() != routed.keys()
                or any(
                    not math.isclose(raw[y], routed[y], rel_tol=1e-10, abs_tol=1e-30)
                    for y in raw
                )
            ):
                raise AssertionError("Corrected routing changed physical parent dates")
            reports[variant] = {
                "package_sha256": resource["sha256"],
                "parent_supplier": constructor,
                "quantity_per_net_service_unit": target,
                "raw_construction_pulses": dict(sorted(raw.items())),
                "routed_construction_pulses": dict(sorted(routed.items())),
                "raw_mean_age_years": mean_age(raw, year),
                "routed_mean_age_years": mean_age(routed, year),
            }
        del model
        gc.collect()
    return {
        "service_year": year,
        "amount": 1,
        "interpolate_annual": False,
        "scope": "Parent manufacture only; component, lifecycle and sign-error evidence is separate",
        "rights": "Restricted local ecoinvent/IAM-derived diagnostic",
        "reports": reports,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("research_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--service-year", type=int, default=2030)
    args = parser.parse_args()
    rows, checks = [], {}
    for label, name in PILOTS.items():
        directory = args.research_dir / name
        report = compare(directory, args.service_year)
        path = directory / f"checks-parent-calendar-{args.service_year}.json"
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        checks[label] = {
            "path": str(path.resolve()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        rows.append(
            {
                "pilot": label,
                "service_year": args.service_year,
                "legacy_raw_mean_age": report["reports"]["legacy"][
                    "raw_mean_age_years"
                ],
                "legacy_routed_mean_age": report["reports"]["legacy"][
                    "routed_mean_age_years"
                ],
                "corrected_mean_age": report["reports"]["corrected"][
                    "routed_mean_age_years"
                ],
            }
        )
        print(json.dumps(rows[-1]), flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "parent-calendar-comparison.csv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "checks": checks,
                "table_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                "rights": "Restricted local evidence; do not redistribute",
            },
            indent=2,
        )
        + "\n"
    )
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(rows))
    for offset, key, label, colour in (
        (-0.25, "legacy_raw_mean_age", "Legacy distribution", "#999999"),
        (0, "legacy_routed_mean_age", "Legacy routed calendar", "#E69F00"),
        (0.25, "corrected_mean_age", "Corrected stock composition", "#0072B2"),
    ):
        ax.bar(x + offset, [r[key] for r in rows], 0.25, label=label, color=colour)
    ax.set_xticks(x, [r["pilot"].replace(" ", "\n", 1) for r in rows])
    ax.set(
        ylabel="Mean parent age (years)",
        title=f"Parent manufacture for service in {args.service_year}",
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    fig.text(
        0.5,
        0.01,
        "Common amortisation retained. Parent timing only; separate component/lifecycle and signed-waste checks remain necessary.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(args.output_dir / f"parent-calendar-comparison.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
