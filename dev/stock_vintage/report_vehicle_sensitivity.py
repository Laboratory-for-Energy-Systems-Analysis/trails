"""Summarise restricted local car/truck cohort sensitivities and provenance."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

CASES = {
    "primary": ("Primary", "#0072B2", "-"),
    "half_year_entrants": ("Half-year entrants", "#009E73", "--"),
    "declining_utilisation": ("Use declines with age", "#D55E00", ":"),
    "reported_sales": ("Reported-sales constraint", "#CC79A7", "-."),
    "shorter_survival": ("Shorter survival support", "#E69F00", "--"),
    "longer_survival": ("Longer survival support", "#555555", ":"),
}


def mean_offset(record):
    return math.fsum(
        (y - record["service_year"]) * w
        for y, w in zip(record["event_years"], record["weights"])
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohorts", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reports, sources, rows = {}, [], []
    for path in args.cohorts:
        report = json.loads(path.read_text())
        key = (report["group"], report["iam_source"]["scenario"])
        if key in reports:
            raise ValueError("Duplicate group/scenario report")
        reports[key] = report
        sources.append(
            {
                "file": str(path.resolve()),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
        for name, case in report["cases"].items():
            ends = {r["service_year"]: r for r in case["parent_retirement"]}
            batteries = {
                r["service_year"]: r for r in case.get("battery_manufacture", [])
            }
            battery_ends = {
                r["service_year"]: r for r in case.get("battery_retirement", [])
            }
            balances = {r["year"]: r for r in case["balances"]}
            for annual in case["annual"]:
                year = annual["service_year"]
                balance = balances.get(year, {})
                rows.append(
                    {
                        "group": key[0],
                        "scenario": key[1],
                        "case": name,
                        "service_year": year,
                        "stock_vehicles": annual["total_stock_vehicles"],
                        "service_index": annual["total_service_index"],
                        "parent_service_weighted_age_years": annual["mean_age_years"],
                        "assumed_parent_mean_years_to_disposal": mean_offset(
                            ends[year]
                        ),
                        "assumed_battery_service_weighted_age_years": (
                            -mean_offset(batteries[year]) if batteries else None
                        ),
                        "assumed_battery_mean_years_to_disposal": (
                            mean_offset(battery_ends[year]) if battery_ends else None
                        ),
                        "stock_balance_residual_vehicles": balance.get(
                            "balance_residual", 0
                        ),
                        "stock_target_residual_vehicles": balance.get(
                            "stock_target_residual_vehicles", 0
                        ),
                        "sales_residual_vehicles": balance.get(
                            "sales_residual_vehicles", 0
                        ),
                        "territorial_exits_vehicles": balance.get("early_exits", 0),
                        "selected_reference_vehicles": case["observation_selection"][
                            "selected_vehicles"
                        ],
                        "exchange_amount_changed": False,
                    }
                )
    groups = ["passenger_bev", "heavy_trucks_32_40t"]
    scenarios = sorted({key[1] for key in reports})
    if set(reports) != {(g, s) for g in groups for s in scenarios}:
        raise ValueError("Require both vehicle groups for every supplied scenario")
    with (args.output_dir / "vehicle-sensitivities.csv").open(
        "w", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(
        2, len(scenarios), figsize=(6 * len(scenarios), 8), squeeze=False
    )
    for i, group in enumerate(groups):
        for j, scenario in enumerate(scenarios):
            ax = axes[i, j]
            report = reports[group, scenario]
            for name, (label, colour, style) in CASES.items():
                annual = report["cases"][name]["annual"]
                ax.plot(
                    [r["service_year"] for r in annual],
                    [r["mean_age_years"] for r in annual],
                    label=label,
                    color=colour,
                    linestyle=style,
                    linewidth=2 if name == "primary" else 1.4,
                )
            title = (
                "Passenger BEVs"
                if group == "passenger_bev"
                else "Articulated diesel trucks, 32–40 t"
            )
            ax.set(
                title=f"{title}\n{scenario}",
                xlabel="Service year",
                ylabel="Service-weighted parent age (years)",
            )
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(axis="y", alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.025),
    )
    fig.suptitle(
        "Vehicle timing: observed UK first-use cohorts and declared future assumptions",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.012,
        "End-year stock snapshots; EUR scenario proxy. Original common-amortised quantities stay unchanged.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.10, 1, 0.95))
    for extension in ("png", "svg", "pdf"):
        fig.savefig(args.output_dir / f"vehicle-sensitivities.{extension}", dpi=180)
    plt.close(fig)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "sources": sources,
                "rows": len(rows),
                "rights": "Restricted local IAM-derived evidence; do not redistribute without source authorisation",
                "interpretation": "Timing under unchanged common amortisation; service, survival and battery histories are assumptions, not exact native IAM vintages",
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"output_dir": str(args.output_dir), "rows": len(rows)}))


if __name__ == "__main__":
    main()
