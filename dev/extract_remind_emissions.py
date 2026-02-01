"""Extract REMIND emissions variables from .mif files into a single CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

INPUT_DIR = Path("/Users/romain/Desktop/remind scenarios")
OUTPUT_PATH = Path("../trails/data/scenarios/remind_emissions.csv")

VARIABLES = [
    "Emi|BC",
    "Emi|C2F6",
    "Emi|C6F14",
    "Emi|CF4",
    "Emi|CH4",
    "Emi|CO",
    "Emi|CO2",
    "Emi|CO2|+|Land-Use Change",
    "Emi|HFC|HFC125",
    "Emi|HFC|HFC134a",
    "Emi|HFC|HFC143a",
    "Emi|HFC|HFC227ea",
    "Emi|HFC|HFC23",
    "Emi|HFC|HFC245fa",
    "Emi|HFC|HFC32",
    "Emi|HFC|HFC43-10",
    "Emi|N2O",
    "Emi|NH3",
    "Emi|NOX",
    "Emi|OC",
    "Emi|SF6",
    "Emi|SO2",
    "Emi|VOC",
]


def _year_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.isdigit()]


def extract_remind_emissions(
    input_dir: Path = INPUT_DIR,
    output_path: Path = OUTPUT_PATH,
) -> Path:
    frames: list[pd.DataFrame] = []
    for path in sorted(input_dir.glob("*.mif")):
        df = pd.read_csv(path, sep=";")
        df.columns = [c.lower() for c in df.columns]
        df = df[(df["region"] == "World") & (df["variable"].isin(VARIABLES))].copy()
        if df.empty:
            continue
        frames.append(df)

    if not frames:
        raise RuntimeError(f"No matching variables found in {input_dir}.")

    out = pd.concat(frames, ignore_index=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


if __name__ == "__main__":
    out = extract_remind_emissions()
    print(f"Wrote {out}")
