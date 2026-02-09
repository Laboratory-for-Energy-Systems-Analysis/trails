"""Merge FAIR historical extensions with REMIND emissions scenarios."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REMIND_PATH = Path("../trails/data/scenarios/remind_emissions.csv")
FAIR_EXT_PATH = Path("../trails/data/scenarios/extensions_1750-2500.csv")
OUTPUT_PATH = Path("../trails/data/scenarios/remind_fair_1750-2500.csv")

EXTENSION_SCENARIO = "medium-extension"

REMIND_TO_FAIR = {
    "Emi|BC": "BC",
    "Emi|OC": "OC",
    "Emi|NH3": "NH3",
    "Emi|NOX": "NOx",
    "Emi|VOC": "VOC",
    "Emi|CO": "CO",
    "Emi|SO2": "Sulfur",
    "Emi|CH4": "CH4",
    "Emi|N2O": "N2O",
    "Emi|CO2": "CO2 FFI",
    "Emi|CO2|+|Land-Use Change": "CO2 AFOLU",
    "Emi|CF4": "CF4",
    "Emi|C2F6": "C2F6",
    "Emi|C6F14": "C6F14",
    "Emi|SF6": "SF6",
    "Emi|HFC|HFC125": "HFC-125",
    "Emi|HFC|HFC134a": "HFC-134a",
    "Emi|HFC|HFC143a": "HFC-143a",
    "Emi|HFC|HFC227ea": "HFC-227ea",
    "Emi|HFC|HFC23": "HFC-23",
    "Emi|HFC|HFC245fa": "HFC-245fa",
    "Emi|HFC|HFC32": "HFC-32",
    "Emi|HFC|HFC43-10": "HFC-4310mee",
}


def _year_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.replace(".", "", 1).isdigit()]


def _fair_year_map(cols: list[str]) -> dict[str, int]:
    mapping = {}
    for c in cols:
        if c.isdigit():
            mapping[c] = int(c)
        else:
            v = float(c)
            mapping[c] = int(v - 0.5)
    return mapping


def _load_fair_extensions(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.lower() for c in df.columns]
    df = df[(df["scenario"] == EXTENSION_SCENARIO) & (df["region"] == "World")].copy()

    cols = _year_cols(df)
    year_map = _fair_year_map(cols)
    df = df[["scenario", "region", "variable", "unit"] + cols].copy()
    df = df.rename(columns=year_map)
    return df


def _load_remind(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.lower() for c in df.columns]
    df = df[df["region"] == "World"].copy()
    df["variable"] = df["variable"].map(REMIND_TO_FAIR).fillna(df["variable"])
    df["scenario"] = df["model"].astype(str) + " | " + df["scenario"].astype(str)
    df = df.drop(columns=["model"])
    return df


def merge_emissions() -> Path:
    fair_df = _load_fair_extensions(FAIR_EXT_PATH)
    remind_df = _load_remind(REMIND_PATH)

    year_cols_fair = [c for c in fair_df.columns if isinstance(c, int)]
    year_cols_remind = [c for c in remind_df.columns if c.isdigit()]
    remind_df = remind_df.rename(columns={c: int(c) for c in year_cols_remind})
    remind_years = sorted(int(c) for c in year_cols_remind) if year_cols_remind else []
    remind_end_year = max(remind_years) if remind_years else 2100

    years = list(range(1750, 2501))

    fair_by_var = {v: row for _, row in fair_df.iterrows() for v in [row["variable"]]}

    rows = []
    for scenario in sorted(remind_df["scenario"].unique()):
        rem_s = remind_df[remind_df["scenario"] == scenario]
        rem_by_var = {row["variable"]: row for _, row in rem_s.iterrows()}

        variables = sorted(set(fair_by_var.keys()) | set(rem_by_var.keys()))
        for var in variables:
            fair_row = fair_by_var.get(var)
            rem_row = rem_by_var.get(var)

            unit = (
                rem_row["unit"]
                if rem_row is not None
                else (fair_row["unit"] if fair_row is not None else "")
            )

            data = {y: float("nan") for y in years}

            if fair_row is not None:
                for y in year_cols_fair:
                    if 1750 <= y <= 2500:
                        data[y] = fair_row[y]

            if rem_row is not None:
                for y in years:
                    if 2005 <= y <= remind_end_year and y in rem_row:
                        data[y] = rem_row[y]

            # Hold at last REMIND year through 2500
            hold = data.get(remind_end_year)
            if hold is not None and pd.notna(hold):
                for y in range(remind_end_year + 1, 2501):
                    data[y] = hold

            row = {
                "scenario": scenario,
                "region": "World",
                "variable": var,
                "unit": unit,
            }
            row.update(data)
            rows.append(row)

    out = pd.DataFrame(rows)
    out = out[["scenario", "region", "variable", "unit"] + years]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_PATH, index=False)
    return OUTPUT_PATH


if __name__ == "__main__":
    out = merge_emissions()
    print(f"Wrote {out}")
