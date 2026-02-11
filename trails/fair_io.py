"""FaIR IO helpers and utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

import fair
from fair.io import DEFAULT_PROPERTIES_FILE

from .filesystem_constants import DATA_DIR

DEFAULT_EMISSIONS_CSV = (
    DATA_DIR / "scenarios" / "extensions_1750-2500_with_scenarios.csv"
)
DEFAULT_MAPPING_YAML = DATA_DIR / "scenarios" / "fair_species_map.yaml"
DEFAULT_CONFIGS_CSV: Path | None = (
    DATA_DIR / "scenarios" / "calibrated_constrained_parameters_calibration1.4.1.csv"
)
DEFAULT_PROPERTIES_CSV: Path | None = (
    DATA_DIR / "scenarios" / "species_configs_properties_calibration1.4.1.csv"
)


def _extend_years_freeze_last(
    df: pd.DataFrame, *, target_year: int = 2500
) -> pd.DataFrame:
    """Extend year columns to target_year by freezing last available values."""
    year_cols, year_vals = _extract_year_columns(df)
    if not year_cols:
        return df
    max_year = max(year_vals)
    has_half = any(abs(v - round(v)) > 1e-9 for v in year_vals)
    target = float(target_year) + (0.5 if has_half else 0.0)
    if max_year >= target:
        return df
    last_col = year_cols[-1]
    # Step by 1.0 year (preserving .5 if present)
    y = max_year + 1.0
    while y <= target + 1e-9:
        col = f"{y:.1f}" if has_half else str(int(y))
        df[col] = df[last_col]
        y += 1.0
    return df


def load_emissions_csv(path: str | Path = DEFAULT_EMISSIONS_CSV) -> pd.DataFrame:
    """Load the merged REMIND/FAIR emissions CSV (IAMC-style)."""
    path = Path(path)
    df = pd.read_csv(path)
    df = _normalize_emissions_columns(df)
    df = _extend_years_freeze_last(df, target_year=2500)
    return df


def load_species_mapping(
    path: str | Path = DEFAULT_MAPPING_YAML,
) -> tuple[dict[object, str], dict[object, float]]:
    """Load biosphere flow -> FaIR species mapping and sign overrides."""
    data = yaml.safe_load(Path(path).read_text()) or {}
    species_map = data.get("species_map", {}) or {}
    signs = data.get("signs", {}) or {}
    parsed_map: dict[object, str] = {}
    for k, v in species_map.items():
        if isinstance(k, (list, tuple)) and len(k) >= 1:
            parsed_map[tuple(k)] = str(v)
        else:
            parsed_map[str(k)] = str(v)
    parsed_signs: dict[object, float] = {}
    for k, v in signs.items():
        if isinstance(k, (list, tuple)) and len(k) >= 1:
            parsed_signs[tuple(k)] = float(v)
        else:
            parsed_signs[str(k)] = float(v)
    return parsed_map, parsed_signs


def _candidate_paths(filename: str) -> list[Path]:
    """Return candidate paths in the FAIR repo tree for a given filename."""
    root = Path(fair.__file__).resolve().parent
    targets = [
        Path("examples") / "data" / "importing-data" / filename,
        Path("examples") / "data" / "calibrated_constrained_ensemble" / filename,
    ]
    return [
        parent / target for parent in [root, *root.parents[:6]] for target in targets
    ]


def _find_fair_repo_file(filename: str) -> Path | None:
    """Try to locate a file in the FAIR repo tree near the installed package."""
    for candidate in _candidate_paths(filename):
        if candidate.exists():
            return candidate
    return None


def _find_any_fair_repo_file(filenames: list[str]) -> Path | None:
    """Return the first matching FAIR repo file found from a list of names."""
    for name in filenames:
        path = _find_fair_repo_file(name)
        if path is not None:
            return path
    return None


def _convert_kg_to_unit(values_kg: np.ndarray, unit: str) -> np.ndarray:
    unit = unit.strip()
    if unit.startswith("Gt "):
        return values_kg / 1e12
    if unit.startswith("Mt "):
        return values_kg / 1e9
    if unit.startswith("kt "):
        return values_kg / 1e6
    if unit.startswith("t "):
        return values_kg / 1e3
    raise ValueError(f"Unsupported unit for emissions conversion: {unit}")


def _convert_unit_to_kg(values: np.ndarray, unit: str) -> np.ndarray:
    """Convert from IAMC mass units back to kg for emissions."""
    factor = _convert_kg_to_unit(np.array([1.0], dtype=float), unit)[0]
    if factor == 0:
        raise ValueError(f"Unsupported unit for emissions conversion: {unit}")
    return values / factor


def _normalize_emissions_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize emissions DataFrame columns for FaIR fill_from_pandas."""
    cols = []
    for c in df.columns:
        if isinstance(c, str):
            c_str = c.strip()
        else:
            c_str = str(c)
        if c_str.lower() in {"scenario", "region", "variable", "unit"}:
            cols.append(c_str.lower())
            continue
        try:
            c_val = float(c_str)
            if c_val.is_integer():
                cols.append(str(int(c_val)))
                continue
        except ValueError:
            pass
        cols.append(c_str.lower())
    df = df.copy()
    df.columns = cols
    return df


def _extract_year_columns(df: pd.DataFrame) -> tuple[list[str], list[float]]:
    """Return ordered IAMC year columns and their numeric values."""
    meta_cols = {"scenario", "region", "variable", "unit"}
    year_cols: list[str] = []
    year_vals: list[float] = []
    for col in df.columns:
        if col in meta_cols:
            continue
        try:
            val = float(col)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(val):
            continue
        year_cols.append(col)
        year_vals.append(float(val))

    if not year_cols:
        return [], []

    order = np.argsort(np.array(year_vals, dtype=float))
    year_cols = [year_cols[i] for i in order]
    year_vals = [year_vals[i] for i in order]
    return year_cols, year_vals


def _default_properties_file() -> Path:
    """Return the default properties file path from fair.io."""
    return Path(DEFAULT_PROPERTIES_FILE)
