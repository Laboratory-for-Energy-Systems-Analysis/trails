"""Export auditable water-cohort sensitivities and a publication-ready figure."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohorts", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = json.loads(args.cohorts.read_text())
    rows = []
    for case, result in data["cases"].items():
        for asset, record in result["assets"].items():
            for annual, end in zip(record["annual"], record["retirement"]):
                year = annual["service_year"]
                if end["service_year"] != year:
                    raise ValueError("Service and retirement records are misaligned")
                rows.append(
                    {
                        "case": case,
                        "asset": asset,
                        "service_year": year,
                        "stock": annual["total_stock"],
                        "mean_age": annual["mean_age_years"],
                        "mean_remaining_years": math.fsum(
                            (y - year) * w
                            for y, w in zip(end["event_years"], end["weights"])
                        ),
                        "unknown_share_2022": record["initial_reconciliation"][
                            "unknown_share"
                        ],
                        "terminal_tail_mass": end["tail_mass_folded_into_last_year"],
                    }
                )
    with (args.output_dir / "water-sensitivities.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    figure, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    cases = list(data["cases"])
    labels = [c.replace("_", " ") for c in cases]
    for column, asset in enumerate(("pipes", "storage")):
        observed = [
            r for r in rows if r["asset"] == asset and r["service_year"] == 2022
        ]
        future = [r for r in rows if r["asset"] == asset and r["service_year"] == 2030]
        axes[0, column].barh(labels, [r["mean_age"] for r in observed], color="#377eb8")
        axes[0, column].invert_yaxis()
        axes[0, column].set_xlabel("Service-weighted mean age in 2022 (years)")
        axes[0, column].set_title(
            f"Quebec {asset}: {100*observed[0]['unknown_share_2022']:.1f}% unknown dates"
        )
        axes[1, column].plot(
            [r["mean_age"] for r in future],
            labels,
            "o",
            label="Service-weighted age in 2030",
            color="#377eb8",
        )
        axes[1, column].plot(
            [r["mean_remaining_years"] for r in future],
            labels,
            "s",
            label="Mean remaining life in 2030",
            color="#e66101",
        )
        axes[1, column].invert_yaxis()
        axes[1, column].set_xlabel("Years")
        axes[1, column].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(axis="x", alpha=0.2)
        axis.set_axisbelow(True)
        axis.tick_params(labelsize=9)
    figure.suptitle(
        "Water-network timing: explicit assumptions, preserved stock quantities\nStatistics Canada 2022 bins; sensitivity scenarios are not confidence intervals",
        fontsize=13,
    )
    for extension in ("png", "svg", "pdf"):
        figure.savefig(args.output_dir / f"water-sensitivities.{extension}", dpi=180)
    plt.close(figure)
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "cohorts_path": str(args.cohorts.resolve()),
                "cohorts_sha256": hashlib.sha256(args.cohorts.read_bytes()).hexdigest(),
                "rows": len(rows),
                "matplotlib_version": matplotlib.__version__,
                "scope": "Stock and conditional-retirement sensitivities; not temporal LCI results",
                "source": data["sources"],
                "rights": data["rights"],
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Wrote {len(rows)} sensitivity records and PNG/SVG/PDF to {args.output_dir}")


if __name__ == "__main__":
    main()
