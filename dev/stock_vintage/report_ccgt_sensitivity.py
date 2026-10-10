"""Write local CCGT cohort sensitivity tables and a publication-style figure.

These are timing sensitivities under unchanged common amortisation. They are
not native IAM vintages or evidence of changed total production burdens. Keep
IAM-derived artifacts local under the input source rights.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

CASES = {
    "primary": ("35-year capacity law", "#0072B2", "-"),
    "capacity_weights": ("Equal service per MW", "#009E73", "--"),
    "short_old_tail": ("Old-stock residual life: 3 y", "#999999", ":"),
    "long_old_tail": ("Old-stock residual life: 10 y", "#555555", "--"),
    "shorter_depreciation": ("25-year capacity law", "#D55E00", "-"),
    "longer_depreciation": ("40-year capacity law", "#CC79A7", "-"),
    "reported_additions": ("Reported-additions constraint", "#E69F00", "-."),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohorts", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources, data, rows = [], [], []
    for path in args.cohorts:
        report = json.loads(path.read_text())
        scenario = report["iam_source"]["scenario"]
        sources.append(
            {
                "file": str(path.resolve()),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
        data.append(report)
        for name, case in report["cases"].items():
            balances = {r["year"]: r for r in case["balances"]}
            for row in case["annual"]:
                year = row["service_year"]
                balance = balances.get(year, {})
                rows.append(
                    {
                        "scenario": scenario,
                        "case": name,
                        "service_year": year,
                        "service_weighted_age_years": row["mean_age_years"],
                        "capacity_mw_ac": row["total_capacity_mw_ac"],
                        "generation_mwh": row["total_generation_mwh"],
                        "maximum_cohort_full_load_hours": max(
                            row["cohort_full_load_hours"].values()
                        ),
                        "capacity_target_residual_mw_ac": balance.get(
                            "capacity_target_residual_mw_ac", 0
                        ),
                        "additions_residual_mw_ac": balance.get(
                            "additions_residual_mw_ac", 0
                        ),
                        "stock_balance_residual_mw_ac": balance.get(
                            "balance_residual", 0
                        ),
                        "construction_exchange_amount_changed": False,
                    }
                )
    with (args.output_dir / "ccgt-sensitivities.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    figure, axes = plt.subplots(
        1, len(data) + 1, figsize=(6 * (len(data) + 1), 5), squeeze=False
    )
    ax = axes[0, 0]
    initial = data[0]
    block_count = len(initial["observation_audit"]["candidate_blocks"])
    stock = {int(k): v for k, v in initial["initial_capacity_mw_ac"].items()}
    output = {int(k): v for k, v in initial["initial_generation_mwh"].items()}
    years = sorted(stock)
    total_stock, total_output = math.fsum(stock.values()), math.fsum(output.values())
    ax.bar(
        [y - 0.2 for y in years],
        [100 * stock[y] / total_stock for y in years],
        width=0.4,
        label="Capacity shares",
        color="#56B4E9",
    )
    ax.bar(
        [y + 0.2 for y in years],
        [100 * output[y] / total_output for y in years],
        width=0.4,
        label="Electricity service shares",
        color="#0072B2",
    )
    ax.set(
        xlabel="Initial commercial-operation year",
        ylabel="Share of the selected portfolio (%)",
        title="Observed 2022 CCGT portfolio",
    )
    ax.legend(frameon=False)
    ax.text(
        0.02,
        0.98,
        f"{block_count} complete blocks\n{total_stock:,.1f} MW AC\nEIA-860 / EIA-923",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    for ax, report in zip(axes[0, 1:], data):
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
        ax.set(
            xlabel="Service year",
            ylabel="Service-weighted stock age (years)",
            title=report["iam_source"]["scenario"] + " reconstruction",
        )
        ax.grid(axis="y", alpha=0.2)
        ax.legend(frameon=False, fontsize=8, loc="best")
    figure.suptitle(
        "Stock-vintage timing: observed reference and declared reconstruction sensitivities",
        fontsize=14,
    )
    figure.text(
        0.5,
        0.01,
        "Selected US-WECC non-CHP CCGT blocks; USA relative IAM trajectories are a proxy. Every case preserves the original capital coefficient.",
        ha="center",
        fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.035, 1, 0.95))
    for extension in ("png", "svg", "pdf"):
        figure.savefig(args.output_dir / f"ccgt-sensitivities.{extension}", dpi=180)
    plt.close(figure)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "sources": sources,
                "rows": len(rows),
                "interpretation": "Cohort timing sensitivity; neither changed common amortisation nor native IAM reproduction",
                "rights": "Local IAM-derived scenario evidence; original source rights apply",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        f"Wrote {len(rows)} scenario/case/year records and three figure formats to {args.output_dir}"
    )


if __name__ == "__main__":
    main()
