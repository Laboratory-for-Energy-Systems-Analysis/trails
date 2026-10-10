"""Plot local PV parent/component timing sensitivities with source hashes.

Keep IAM-derived tables and figures local under source rights. These are timing
changes under the original common amortisation, not new material-demand totals.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

CASES = {
    "primary": ("30-year parent; 15-year inverter", "#0072B2", "-"),
    "capacity_weights": ("Equal service per MW AC", "#009E73", "--"),
    "shorter_parent": ("Parent: 25 years", "#D55E00", "-"),
    "longer_parent": ("Parent: 40 years", "#CC79A7", "-"),
    "quartic_parent": ("Quartic parent survival", "#555555", ":"),
    "shorter_inverter": ("Inverter: 10 years", "#56B4E9", "--"),
    "longer_inverter": ("Inverter: 20 years", "#882255", "-."),
    "reported_additions": ("Reported additions", "#E69F00", "-."),
}


def mean_offset(record):
    return math.fsum(
        (year - record["service_year"]) * weight
        for year, weight in zip(record["event_years"], record["weights"])
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohorts", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources, reports, rows = [], [], []
    for path in args.cohorts:
        report = json.loads(path.read_text())
        reports.append(report)
        sources.append(
            {
                "file": str(path.resolve()),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
        for name, case in report["cases"].items():
            balances = {r["year"]: r for r in case["balances"]}
            for annual, parent_end, inverter, inverter_end in zip(
                case["annual"],
                case["parent_retirement"],
                case["inverter_manufacture"],
                case["inverter_retirement"],
            ):
                year = annual["service_year"]
                balance = balances.get(year, {})
                rows.append(
                    {
                        "scenario": report["iam_source"]["scenario"],
                        "case": name,
                        "service_year": year,
                        "parent_service_weighted_age_years": annual["mean_age_years"],
                        "inverter_service_weighted_age_years": -mean_offset(inverter),
                        "parent_mean_remaining_years": mean_offset(parent_end),
                        "inverter_mean_remaining_years": mean_offset(inverter_end),
                        "capacity_mw_ac": annual["total_capacity_mw_ac"],
                        "generation_mwh": annual["total_generation_mwh"],
                        "stock_balance_residual_mw_ac": balance.get(
                            "balance_residual", 0
                        ),
                        "capacity_target_residual_mw_ac": balance.get(
                            "capacity_target_residual_mw_ac", 0
                        ),
                        "additions_residual_mw_ac": balance.get(
                            "additions_residual_mw_ac", 0
                        ),
                        "exchange_amount_changed": False,
                    }
                )
    with (args.output_dir / "pv-sensitivities.csv").open("w", newline="") as stream:
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
        2, 1 + len(reports), figsize=(6 * (1 + len(reports)), 9), squeeze=False
    )
    initial = reports[0]
    stock = {int(k): v for k, v in initial["initial_capacity_mw_ac"].items()}
    output = {int(k): v for k, v in initial["initial_generation_mwh"].items()}
    years = sorted(stock)
    axes[0, 0].bar(
        [y - 0.2 for y in years],
        [100 * stock[y] / math.fsum(stock.values()) for y in years],
        width=0.4,
        color="#56B4E9",
        label="AC stock share",
    )
    axes[0, 0].bar(
        [y + 0.2 for y in years],
        [100 * output[y] / math.fsum(output.values()) for y in years],
        width=0.4,
        color="#0072B2",
        label="Electricity service share",
    )
    axes[0, 0].set(
        title="Observed 2022 parent portfolio",
        xlabel="Plant commissioning year",
        ylabel="Share (%)",
    )
    axes[0, 0].legend(frameon=False)
    inverter = initial["cases"]["primary"]["inverter_manufacture"][0]
    axes[1, 0].bar(
        inverter["event_years"], [100 * w for w in inverter["weights"]], color="#0072B2"
    )
    axes[1, 0].set(
        title="Modelled active inverters in 2022",
        xlabel="Inverter manufacture year",
        ylabel="Service share (%)",
    )
    axes[1, 0].text(
        0.03,
        0.97,
        "Regular 15-year replacement assumed\nComponent dates are not observed",
        transform=axes[1, 0].transAxes,
        va="top",
        fontsize=9,
    )
    for column, report in enumerate(reports, start=1):
        for name, (label, colour, style) in CASES.items():
            case = report["cases"][name]
            dates = [r["service_year"] for r in case["annual"]]
            for row, values in enumerate(
                [
                    [r["mean_age_years"] for r in case["annual"]],
                    [-mean_offset(r) for r in case["inverter_manufacture"]],
                ]
            ):
                axes[row, column].plot(
                    dates,
                    values,
                    label=label,
                    color=colour,
                    linestyle=style,
                    linewidth=2 if name == "primary" else 1.4,
                )
                axes[row, column].set(
                    xlabel="Service year", ylabel="Service-weighted age (years)"
                )
                axes[row, column].grid(axis="y", alpha=0.2)
        axes[0, column].set_title(report["iam_source"]["scenario"] + ": parent")
        axes[1, column].set_title(report["iam_source"]["scenario"] + ": inverter")
    for ax in axes.flat:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.025),
    )
    fig.suptitle(
        "PV timing: observed parent stock and declared component/scenario assumptions",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.012,
        "Bounded fixed-crystalline WECC population; USA relative IAM proxy. All cases preserve original inventory coefficients.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.095, 1, 0.95))
    for extension in ("png", "svg", "pdf"):
        fig.savefig(args.output_dir / f"pv-sensitivities.{extension}", dpi=180)
    plt.close(fig)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "sources": sources,
                "rows": len(rows),
                "rights": "Local IAM-derived evidence; do not redistribute without source authorization",
                "interpretation": "Timing sensitivity under unchanged common amortisation; not observed component history or exact native IAM vintages",
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"output_dir": str(args.output_dir), "rows": len(rows)}))


if __name__ == "__main__":
    main()
