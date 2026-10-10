"""Compare verified PV waste-timing endpoints without changing inventory totals.

Input/output series contain restricted ecoinvent/IAM-derived results. Keep them
outside version control. Emission dates include upstream processes and are not
physical disposal dates or dynamic climate-impact scores.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

CASES = {
    "primary": ("pv-real-primary", None),
    "polymer_at_manufacture": ("pv-polymer-manufacture", "panel_mixed_polymer"),
    "mounting_at_retirement": ("pv-mounting-eol", "mounting_ambiguous_packaging"),
}


def read_case(directory):
    path = directory / "checks-2030.json"
    report = json.loads(path.read_text())
    export = json.loads((directory / "export-evidence.json").read_text())
    if report["package_sha256"] != export["packages"]["corrected"]["sha256"]:
        raise ValueError("Check and export refer to different packages")
    if {(r["service_year"], r["amount"]) for r in report["lci_checks"]} != {
        (2030, 1.0),
        (2030, -2.0),
    }:
        raise ValueError("Require both signed full-LCI checks for 2030")
    if not all(
        r["static_conservation_passed"]
        and r["flow_accumulation_passed"]
        and r["operator_conservation"]["operator_conservation_passed"]
        for r in report["lci_checks"]
    ):
        raise ValueError("Full-LCI conservation failed")
    row = next(r for r in report["lci_checks"] if r["amount"] == 1)
    series = {
        role: {int(y): v for y, v in values.items()}
        for role, values in row["fossil_co2_kg_by_role_and_year"].items()
    }
    return series, {
        "check_file": str(path.resolve()),
        "check_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "package_sha256": report["package_sha256"],
        "max_absolute_flow_error": max(
            r["max_absolute_flow_error"] for r in report["lci_checks"]
        ),
        "max_operator_relative_residual": max(
            r["operator_conservation"]["operator_max_componentwise_relative_residual"]
            for r in report["lci_checks"]
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("research_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sources, data, rows, summaries = {}, {}, [], []
    for case, (directory, _) in CASES.items():
        data[case], sources[case] = read_case(args.research_dir / directory)
    baseline = data["primary"]
    for case, series in data.items():
        if series.keys() != baseline.keys():
            raise ValueError("Endpoint role boundaries differ")
        changed_role = CASES[case][1]
        for role, values in series.items():
            total = math.fsum(values.values())
            original = math.fsum(baseline[role].values())
            if not math.isclose(total, original, rel_tol=1e-10, abs_tol=1e-14):
                raise ValueError(f"Endpoint changed total burden: {case}/{role}")
            difference = math.fsum(
                abs(values.get(y, 0) - baseline[role].get(y, 0))
                for y in values.keys() | baseline[role].keys()
            )
            if role != changed_role and difference > 1e-12:
                raise ValueError(f"Unrelated role changed dates: {case}/{role}")
            summaries.append(
                {
                    "case": case,
                    "role": role,
                    "total_fossil_co2_kg_per_kwh": total,
                    "signed_burden_weighted_mean_emission_year": (
                        math.fsum(y * v for y, v in values.items()) / total
                        if total
                        else None
                    ),
                    "calendar_l1_difference_from_primary_kg_per_kwh": difference,
                }
            )
            rows.extend(
                {
                    "case": case,
                    "role": role,
                    "emission_year": y,
                    "fossil_co2_kg_per_kwh": v,
                }
                for y, v in sorted(values.items())
            )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "pv-endpoint-emissions.csv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "pv-endpoint-summary.json").write_text(
        json.dumps(
            {
                "service_year": 2030,
                "functional_unit": "1 kWh",
                "rights": "Restricted local ecoinvent/IAM-derived evidence",
                "interpretation": "Signed upstream fossil-CO2 dates, not physical event dates or climate scores; endpoint totals and unaffected role calendars conserved",
                "sources": sources,
                "roles": summaries,
                "table_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mplconfig"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, role, endpoint, title in zip(
        axes,
        ["panel_mixed_polymer", "mounting_ambiguous_packaging"],
        ["polymer_at_manufacture", "mounting_at_retirement"],
        ["Mixed panel-polymer waste", "Ambiguous mounting packaging"],
    ):
        for case, label, colour in (
            ("primary", "Primary convention", "#0072B2"),
            (endpoint, "Alternative endpoint", "#D55E00"),
        ):
            values = data[case][role]
            total, cumulative = math.fsum(values.values()), 0
            if total <= 0 or any(v < -1e-18 for v in values.values()):
                raise ValueError("Cumulative plot requires nonnegative burdens")
            years, shares = [], []
            for year, value in sorted(values.items()):
                cumulative += value
                years.append(year)
                shares.append(cumulative / total * 100)
            ax.step(years, shares, where="post", label=label, color=colour)
        ax.axvline(2030, color="#777777", linestyle=":", linewidth=1)
        ax.set(title=title, xlabel="Emission year", ylabel="Cumulative fossil CO₂ (%)")
        ax.set_xlim(1980, 2070)
        ax.set_ylim(-2, 102)
        ax.legend(frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("PV waste timing: same total burdens for one 2030 kWh")
    fig.text(
        0.5,
        0.01,
        "Upstream emission dates include background processes; plotted window 1980–2070. Full calendars are retained in the CSV.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(args.output_dir / f"pv-endpoint-emissions.{suffix}", dpi=180)
    plt.close(fig)
    print(json.dumps({"rows": len(rows), "roles": len(summaries)}))


if __name__ == "__main__":
    main()
