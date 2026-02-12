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
    """extend years freeze last.

    :param df: Value for `df`.
    :type df: pd.DataFrame
    :param target_year: Value for `target_year`.
    :type target_year: int
    :returns: Return value.
    :rtype: pd.DataFrame"""
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
    """Load emissions csv.

    :param path: Value for `path`.
    :type path: str | Path
    :returns: Return value.
    :rtype: pd.DataFrame"""
    path = Path(path)
    df = pd.read_csv(path)
    df = _normalize_emissions_columns(df)
    df = _extend_years_freeze_last(df, target_year=2500)
    return df


def load_species_mapping(
    path: str | Path = DEFAULT_MAPPING_YAML,
) -> tuple[dict[object, str], dict[object, float]]:
    """Load species mapping.

    :param path: Value for `path`.
    :type path: str | Path
    :returns: Return value.
    :rtype: tuple[dict[object, str], dict[object, float]]"""
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
    """candidate paths.

    :param filename: Value for `filename`.
    :type filename: str
    :returns: Return value.
    :rtype: list[Path]"""
    root = Path(fair.__file__).resolve().parent
    targets = [
        Path("examples") / "data" / "importing-data" / filename,
        Path("examples") / "data" / "calibrated_constrained_ensemble" / filename,
    ]
    return [
        parent / target for parent in [root, *root.parents[:6]] for target in targets
    ]


def _find_fair_repo_file(filename: str) -> Path | None:
    """find fair repo file.

    :param filename: Value for `filename`.
    :type filename: str
    :returns: Return value.
    :rtype: Path | None"""
    for candidate in _candidate_paths(filename):
        if candidate.exists():
            return candidate
    return None


def _find_any_fair_repo_file(filenames: list[str]) -> Path | None:
    """find any fair repo file.

    :param filenames: Value for `filenames`.
    :type filenames: list[str]
    :returns: Return value.
    :rtype: Path | None"""
    for name in filenames:
        path = _find_fair_repo_file(name)
        if path is not None:
            return path
    return None


def _convert_kg_to_unit(values_kg: np.ndarray, unit: str) -> np.ndarray:
    """convert kg to unit.

    :param values_kg: Value for `values_kg`.
    :type values_kg: np.ndarray
    :param unit: Value for `unit`.
    :type unit: str
    :returns: Return value.
    :rtype: np.ndarray
    :raises ValueError: If an error occurs."""
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
    """convert unit to kg.

    :param values: Value for `values`.
    :type values: np.ndarray
    :param unit: Value for `unit`.
    :type unit: str
    :returns: Return value.
    :rtype: np.ndarray
    :raises ValueError: If an error occurs."""
    factor = _convert_kg_to_unit(np.array([1.0], dtype=float), unit)[0]
    if factor == 0:
        raise ValueError(f"Unsupported unit for emissions conversion: {unit}")
    return values / factor


def _normalize_emissions_columns(df: pd.DataFrame) -> pd.DataFrame:
    """normalize emissions columns.

    :param df: Value for `df`.
    :type df: pd.DataFrame
    :returns: Return value.
    :rtype: pd.DataFrame"""
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
    """extract year columns.

    :param df: Value for `df`.
    :type df: pd.DataFrame
    :returns: Return value.
    :rtype: tuple[list[str], list[float]]"""
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
    """default properties file.

    :returns: Return value.
    :rtype: Path"""
    return Path(DEFAULT_PROPERTIES_FILE)
