"""FaIR integration helpers."""

from __future__ import annotations

import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import os
from pathlib import Path
from threading import Lock
from typing import Any, Dict, Tuple
import warnings

import pandas as pd
import numpy as np
import sparse
import xarray as xr

from tqdm import tqdm

import fair  # hard dependency
from fair.io import read_properties
from fair.interface import initialise

from .fair_io import (
    DEFAULT_CONFIGS_CSV,
    DEFAULT_EMISSIONS_CSV,
    DEFAULT_MAPPING_YAML,
    DEFAULT_PROPERTIES_CSV,
    _convert_kg_to_unit,
    _convert_unit_to_kg,
    _default_properties_file,
    _extract_year_columns,
    _find_any_fair_repo_file,
    _normalize_emissions_columns,
    load_emissions_csv,
    load_species_mapping,
)
from .chunked_inventory import is_chunked_sparse, iter_sparse_blocks

# FaIR's setup phase is not reliably thread-safe with all SciPy/FaIR
# combinations. Unprepared/fallback runs remain serialized; prepared runs reuse
# setup generated once in the cached template and operate on independent copies.
_FAIR_RUN_LOCK = Lock()


def _resolve_flow_mapping(mapping: dict[object, Any], flow_key: object) -> Any:
    """Resolve mappings that may use either a full flow key or only its name."""
    name_key = (
        str(flow_key[0])
        if isinstance(flow_key, tuple) and len(flow_key) >= 1
        else str(flow_key)
    )
    if name_key in mapping:
        return mapping[name_key]
    return mapping.get(flow_key)


def _inventory_flow_year_root(
    inv: xr.DataArray,
    *,
    flow_positions: list[int] | None = None,
    builder: Any | None = None,
    show_progress: bool = False,
) -> sparse.COO:
    """Reduce an eager or chunked inventory to flow/year/root sequentially."""
    canonical = inv.transpose("activity", "flow", "year", "root activity")
    data = canonical.data
    if not is_chunked_sparse(data):
        if not isinstance(data, sparse.COO):
            data = sparse.COO.from_numpy(np.asarray(data, dtype=float))
        return data.sum(axis=0)

    if (
        builder is not None
        and getattr(builder, "_finalized", False)
        and flow_positions is not None
    ):
        return builder.reduce_activity_for_flows(
            flow_positions,
            show_progress=show_progress,
        )

    coords_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []
    for slices, block in iter_sparse_blocks(data, primary_axis=2):
        flow_offset = int(slices[1].start or 0)
        if flow_positions is not None:
            local_flows = np.asarray(
                [
                    int(flow) - flow_offset
                    for flow in flow_positions
                    if flow_offset <= int(flow) < int(slices[1].stop)
                ],
                dtype=np.int64,
            )
            if not local_flows.size:
                continue
            reduced = block[:, local_flows, :, :].sum(axis=0)
        else:
            local_flows = None
            reduced = block.sum(axis=0)
        if not isinstance(reduced, sparse.COO) or not reduced.nnz:
            continue
        coords = reduced.coords.astype(np.int64, copy=True)
        if local_flows is None:
            coords[0] += flow_offset
        else:
            coords[0] = local_flows[coords[0]] + flow_offset
        coords[1] += int(slices[2].start or 0)
        coords[2] += int(slices[3].start or 0)
        coords_parts.append(coords)
        data_parts.append(reduced.data)
    shape = (
        int(canonical.sizes["flow"]),
        int(canonical.sizes["year"]),
        int(canonical.sizes["root activity"]),
    )
    if not coords_parts:
        return sparse.zeros(shape, dtype=canonical.dtype)
    return sparse.COO(
        np.concatenate(coords_parts, axis=1),
        np.concatenate(data_parts),
        shape=shape,
    )


def _sanitize_emissions_year_values(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce emissions year columns to numeric and replace missing with zero."""
    year_cols, _ = _extract_year_columns(df)
    if not year_cols:
        return df
    out = df.copy()
    # Assign each year column explicitly so pandas>=3 doesn't try in-place casts
    # on existing StringDtype blocks during a multi-column .loc assignment.
    for col in year_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    return out


def _safe_nanpercentile(values: np.ndarray, quantiles: list[float]) -> np.ndarray:
    """Percentiles that avoid warnings on all-NaN slices by treating them as zero."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.size == 0:
        return np.zeros((len(quantiles), 0), dtype=float)
    all_nan_by_time = np.all(np.isnan(arr), axis=0)
    if np.any(all_nan_by_time):
        arr = arr.copy()
        arr[:, all_nan_by_time] = 0.0
    out = np.nanpercentile(arr, quantiles, axis=0)
    return np.nan_to_num(out, nan=0.0)


def _array_to_builtin(value: Any) -> float | list[Any]:
    """Convert numpy scalars/arrays to JSON-friendly values."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return float(arr)
    return arr.tolist()


def _trapezoid(values: np.ndarray, years: np.ndarray, *, axis: int) -> np.ndarray:
    """Compatibility wrapper for NumPy trapezoidal integration."""
    if hasattr(np, "trapezoid"):
        return np.trapezoid(values, x=years, axis=axis)
    return np.trapz(values, x=years, axis=axis)


def integrate_window(
    values: Any,
    years: Any,
    start_year: int = 2050,
    end_year: int = 2100,
    *,
    time_axis: int = -1,
) -> float | np.ndarray:
    """Integrate a time series over a closed year window.

    :param values: Time series values. Multidimensional arrays are supported.
    :type values: Any
    :param years: Year coordinates corresponding to the time axis.
    :type years: Any
    :param start_year: First year included in the assessment window.
    :type start_year: int
    :param end_year: Last year included in the assessment window.
    :type end_year: int
    :param time_axis: Axis in ``values`` corresponding to ``years``.
    :type time_axis: int
    :returns: Trapezoidal integral in value unit times year.
    :rtype: float | np.ndarray
    :raises ValueError: If fewer than two points are available in the window.
    """
    years_arr = np.asarray(years, dtype=float)
    values_arr = np.asarray(values, dtype=float)
    if values_arr.ndim == 0:
        raise ValueError("values must include a time axis.")
    axis = int(time_axis)
    if axis < 0:
        axis += values_arr.ndim
    if axis < 0 or axis >= values_arr.ndim:
        raise ValueError("time_axis is out of range for values.")
    if values_arr.shape[axis] != years_arr.shape[0]:
        raise ValueError("years length must match values along time_axis.")

    mask = (years_arr >= float(start_year)) & (years_arr <= float(end_year))
    if int(mask.sum()) < 2:
        raise ValueError(
            "At least two time points are required for trapezoidal integration."
        )

    indices = np.flatnonzero(mask)
    window_values = np.take(values_arr, indices, axis=axis)
    window_years = years_arr[indices]
    result = _trapezoid(window_values, window_years, axis=axis)
    if np.asarray(result).ndim == 0:
        return float(result)
    return result


def calculate_co2_pulse_equivalents(
    dRF_lca: Any,
    dT_lca: Any,
    dRF_ref: Any,
    dT_ref: Any,
    years: Any,
    reference_pulse_mass: float,
    start_year: int = 2050,
    end_year: int = 2100,
    *,
    time_axis: int = -1,
) -> dict[str, Any]:
    """Calculate fixed-window CO2 pulse-equivalent indicators.

    The ratios are computed before any uncertainty summary, so callers can pass
    arrays shaped like ``config, time`` and then summarize the returned values.

    :param dRF_lca: Incremental RF response of the LCA perturbation.
    :type dRF_lca: Any
    :param dT_lca: Incremental temperature response of the LCA perturbation.
    :type dT_lca: Any
    :param dRF_ref: Incremental RF response of the reference CO2 pulse.
    :type dRF_ref: Any
    :param dT_ref: Incremental temperature response of the reference CO2 pulse.
    :type dT_ref: Any
    :param years: Time coordinate in years.
    :type years: Any
    :param reference_pulse_mass: Size of the reference CO2 pulse.
    :type reference_pulse_mass: float
    :param start_year: First year included in the assessment window.
    :type start_year: int
    :param end_year: Last year included in the assessment window.
    :type end_year: int
    :param time_axis: Axis corresponding to ``years`` in each response array.
    :type time_axis: int
    :returns: Integrated responses and pulse-equivalent values.
    :rtype: dict[str, Any]
    :raises ZeroDivisionError: If a reference response has a zero integral.
    """
    int_rf_lca = integrate_window(
        dRF_lca, years, start_year, end_year, time_axis=time_axis
    )
    int_t_lca = integrate_window(
        dT_lca, years, start_year, end_year, time_axis=time_axis
    )
    int_rf_ref = integrate_window(
        dRF_ref, years, start_year, end_year, time_axis=time_axis
    )
    int_t_ref = integrate_window(
        dT_ref, years, start_year, end_year, time_axis=time_axis
    )

    rf_ref_arr = np.asarray(int_rf_ref, dtype=float)
    t_ref_arr = np.asarray(int_t_ref, dtype=float)
    if np.any(~np.isfinite(rf_ref_arr)) or np.any(rf_ref_arr == 0.0):
        raise ZeroDivisionError(
            "Integrated RF response of reference CO2 pulse is zero or non-finite."
        )
    if np.any(~np.isfinite(t_ref_arr)) or np.any(t_ref_arr == 0.0):
        raise ZeroDivisionError(
            "Integrated temperature response of reference CO2 pulse is zero "
            "or non-finite."
        )

    rf_equivalent = (
        float(reference_pulse_mass) * np.asarray(int_rf_lca, dtype=float) / rf_ref_arr
    )
    temperature_equivalent = (
        float(reference_pulse_mass) * np.asarray(int_t_lca, dtype=float) / t_ref_arr
    )

    return {
        "window": (int(start_year), int(end_year)),
        "reference_pulse_mass": float(reference_pulse_mass),
        "integrated_rf_lca": int_rf_lca,
        "integrated_rf_ref": int_rf_ref,
        "integrated_temperature_lca": int_t_lca,
        "integrated_temperature_ref": int_t_ref,
        "co2_pulse_equivalent_integrated_rf": rf_equivalent,
        "co2_pulse_equivalent_integrated_temperature": temperature_equivalent,
    }


def _summarize_indicator_values(values: Any) -> dict[str, float]:
    """Summarize scalar or per-configuration indicator values."""
    arr = np.asarray(values, dtype=float).reshape(-1)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return {"value": np.nan, "median": np.nan, "p025": np.nan, "p975": np.nan}
    median = float(np.nanmedian(finite))
    return {
        "value": median,
        "median": median,
        "p025": float(np.nanquantile(finite, 0.025)),
        "p975": float(np.nanquantile(finite, 0.975)),
    }


def make_reference_co2_pulse_emissions(
    emissions_df: pd.DataFrame,
    *,
    scenario: str,
    pulse_year: int = 2050,
    pulse_mass_kg: float = 1.0e9,
    co2_species_name: str = "CO2 FFI",
) -> pd.DataFrame:
    """Return baseline emissions with a one-year reference CO2 pulse added.

    The existing Trails-FaIR interface uses IAMC-style annual emission-rate
    tables. For a one-year FaIR time step, adding the pulse mass to the selected
    year column in the row's native unit represents a pulse over that interval.

    :param emissions_df: Baseline emissions table.
    :type emissions_df: pd.DataFrame
    :param scenario: Baseline scenario to perturb.
    :type scenario: str
    :param pulse_year: Calendar year receiving the reference pulse.
    :type pulse_year: int
    :param pulse_mass_kg: Pulse mass in kg CO2.
    :type pulse_mass_kg: float
    :param co2_species_name: Preferred CO2 emissions row to perturb.
    :type co2_species_name: str
    :returns: Perturbed emissions table.
    :rtype: pd.DataFrame
    :raises ValueError: If scenario, year, species, or units cannot be resolved.
    """
    out = _normalize_emissions_columns(emissions_df)
    scenario_mask = (out["scenario"] == scenario) & (
        out["region"].str.lower() == "world"
    )
    if not bool(np.any(scenario_mask)):
        raise ValueError(f"Scenario '{scenario}' not found in emissions data.")

    year_cols, year_vals = _extract_year_columns(out)
    if not year_cols:
        raise ValueError("No year columns found in emissions data.")
    has_half_years = any(abs(v - round(v)) > 1e-9 for v in year_vals)
    pulse_coord = float(pulse_year) + (0.5 if has_half_years else 0.0)
    year_lookup = {float(v): col for col, v in zip(year_cols, year_vals)}
    ycol = year_lookup.get(pulse_coord)
    if ycol is None:
        raise ValueError(
            f"Reference pulse year {pulse_year} is not available in emissions data."
        )

    candidates: list[str] = []
    for candidate in [co2_species_name, "CO2 FFI", "CO2", "CO2 AFOLU"]:
        candidate_str = str(candidate)
        if candidate_str not in candidates:
            candidates.append(candidate_str)

    row_idx: Any | None = None
    for candidate in candidates:
        matches = out.index[scenario_mask & (out["variable"] == candidate)].tolist()
        if matches:
            row_idx = matches[0]
            break
    if row_idx is None:
        raise ValueError(
            "No CO2 emissions row found for reference pulse. Tried: "
            + ", ".join(candidates)
        )

    unit = str(out.loc[row_idx, "unit"])
    add_value = _convert_kg_to_unit(np.array([pulse_mass_kg], dtype=float), unit)[0]
    current = pd.to_numeric(out.loc[row_idx, ycol], errors="coerce")
    if pd.isna(current):
        current = 0.0
    out.loc[row_idx, ycol] = float(current) + float(add_value)
    return out


# Emissions precursor species whose forcing mostly appears in derived FaIR species.
_PRECURSOR_RESPONSE_SPECIES: dict[str, tuple[str, ...]] = {
    "CO": ("Ozone",),
    "NH3": ("Aerosol-radiation interactions", "Ozone"),
    "NOx": ("Ozone", "Aerosol-radiation interactions"),
    "Sulfur": (
        "Aerosol-radiation interactions",
        "Aerosol-cloud interactions",
        "Ozone",
    ),
    "VOC": ("Ozone", "Aerosol-radiation interactions"),
    "BC": (
        "Aerosol-radiation interactions",
        "Aerosol-cloud interactions",
        "Light absorbing particles on snow and ice",
    ),
    "OC": ("Aerosol-radiation interactions", "Aerosol-cloud interactions"),
}


_AEROSOL_OZONE_PRECURSOR_SPECIES = frozenset(_PRECURSOR_RESPONSE_SPECIES)


def _exclude_aerosol_ozone_precursor_mapping(
    species_map: dict[object, str],
    signs: dict[object, float],
) -> tuple[dict[object, str], dict[object, float]]:
    """Remove aerosol/ozone precursor species from a FaIR flow mapping."""
    filtered_map = {
        flow_key: specie
        for flow_key, specie in species_map.items()
        if str(specie) not in _AEROSOL_OZONE_PRECURSOR_SPECIES
    }
    filtered_signs = {
        flow_key: sign for flow_key, sign in signs.items() if flow_key in filtered_map
    }
    return filtered_map, filtered_signs


def _ensure_response_species_rows(
    emissions_df: pd.DataFrame,
    *,
    scenario: str,
    drivers: list[str],
    debug: bool = False,
) -> pd.DataFrame:
    """Ensure required calculated response species exist as zero rows in emissions."""
    needed: set[str] = set()
    for driver in drivers:
        needed.update(_PRECURSOR_RESPONSE_SPECIES.get(str(driver), ()))
    if not needed:
        return emissions_df

    out = _normalize_emissions_columns(emissions_df)
    scen_mask = (out["scenario"] == scenario) & (out["region"].str.lower() == "world")
    if not bool(np.any(scen_mask)):
        return out

    year_cols, _ = _extract_year_columns(out)
    if not year_cols:
        return out

    existing = set(out.loc[scen_mask, "variable"].astype(str).tolist())
    add_rows: list[dict[str, object]] = []
    for specie in sorted(needed):
        if specie in existing:
            continue
        row: dict[str, object] = {
            "scenario": scenario,
            "region": "World",
            "variable": specie,
            # Calculated forcing channels are represented in forcing units.
            "unit": "W m-2",
        }
        for col in year_cols:
            row[col] = 0.0
        add_rows.append(row)

    if not add_rows:
        return out

    out = pd.concat([out, pd.DataFrame(add_rows)], ignore_index=True)
    if debug:
        print(
            "FAIR debug: added response species rows",
            [str(r["variable"]) for r in add_rows],
        )
    return out


@lru_cache(maxsize=2)
def _build_fair_template_cached(
    *,
    scenario: str,
    start_year: float,
    end_year: float,
    config_csv: str,
    properties_csv: str,
    species_key: tuple[str, ...],
    configs_key: tuple[object, ...],
    ghg_method: str | None,
    temperature_prescribed: bool | None,
) -> fair.FAIR:
    """Build and cache a configured FaIR template for repeated perturbation runs."""
    species, properties = read_properties(
        filename=properties_csv, species=list(species_key)
    )

    if temperature_prescribed is None:
        f = fair.FAIR()
    else:
        f = fair.FAIR(temperature_prescribed=bool(temperature_prescribed))
    if ghg_method is not None:
        f.ghg_method = ghg_method
    f.define_time(float(start_year), float(end_year), 1)
    f.define_scenarios([str(scenario)])
    f.define_configs(list(configs_key))
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(filename=properties_csv)
    f.override_defaults(config_csv)

    # Initialize arrays to stable values before each deep-copied run.
    f.emissions.data[...] = 0
    initialise(f.temperature, 0)
    initialise(f.forcing, 0)

    # These setup products only depend on the time axis, species, and calibrated
    # configurations, all of which are part of this function's cache key.  FaIR
    # normally rebuilds them at the start of every run; doing that for every
    # per-species perturbation is particularly costly for large ensembles.
    f._check_properties()
    f._make_indices()
    if f._routine_flags["temperature"]:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=RuntimeWarning,
                module="scipy.stats._multivariate",
            )
            f._make_ebms()
    f._trails_prepared_setup = True
    return f


def _run_prepared_fair(f: fair.FAIR, *, progress: bool) -> None:
    """Run FaIR while reusing immutable setup cached in a copied template."""
    if not getattr(f, "_trails_prepared_setup", False):
        f.run(progress=progress)
        return

    # ``FAIR.run`` unconditionally rebuilds these objects.  Shadow the setup
    # methods only for this call; the cached indices and EBM matrices copied
    # from the template remain available to the numerical integration.
    method_names = ("_check_properties", "_make_indices", "_make_ebms")
    for name in method_names:
        setattr(f, name, lambda: None)
    try:
        f.run(progress=progress)
    finally:
        for name in method_names:
            delattr(f, name)


def _fill_emissions_from_df_fast(
    f: fair.FAIR,
    emissions_df: pd.DataFrame,
    *,
    scenario: str,
    year_cols: list[str],
    year_vals: list[float],
) -> None:
    """Fast path for filling FaIR emissions without per-species xarray indexing."""
    from fair.exceptions import DuplicateScenarioError
    from fair.io.fill_from import _emissions_unit_convert

    if emissions_df.empty:
        return

    dup_counts = emissions_df.groupby("variable", sort=False).size()
    duplicates = dup_counts[dup_counts > 1]
    if not duplicates.empty:
        specie = str(duplicates.index[0])
        raise DuplicateScenarioError(
            "Input data for emissions contains duplicate "
            f"rows for variable='{specie}, scenario='{scenario}'."
        )

    df_by_var = emissions_df.set_index("variable", drop=False)
    times = np.asarray(year_vals, dtype=float)
    target_times = np.asarray(f.timepoints, dtype=float)
    same_grid = times.shape == target_times.shape and np.allclose(times, target_times)

    emissions_data = np.asarray(f.emissions.values, dtype=float)
    for specie_idx, specie in enumerate(f.species):
        if f.properties_df.loc[specie, "input_mode"] != "emissions":
            continue
        if specie not in df_by_var.index:
            continue

        row = df_by_var.loc[specie]
        data_in = row[year_cols].to_numpy(dtype=float, copy=False)
        if same_grid:
            data = data_in
        else:
            data = np.interp(target_times, times, data_in, left=np.nan, right=np.nan)

        unit = str(row["unit"])
        is_ghg = bool(f.properties_df.loc[specie, "greenhouse_gas"])
        converted = _emissions_unit_convert(data, unit, specie, is_ghg)
        emissions_data[:, 0, :, specie_idx] = converted[:, None]

    f.emissions.data = emissions_data


def _inventory_emissions_by_fair_species(
    inv_data: sparse.COO,
    inv_years: list[int],
    n_flow: int,
    flow_pos_to_key: dict[int, tuple[str, str, str] | str],
    species_map: dict[object, str],
    signs: dict[object, float],
) -> pd.DataFrame:
    """inventory emissions by fair species.

    :param inv_data: Activity-reduced inventory as flow/year/root sparse COO.
    :type inv_data: sparse.COO
    :param inv_years: Inventory years (matching inv_data year axis).
    :type inv_years: list[int]
    :param n_flow: Number of flows.
    :type n_flow: int
    :param flow_pos_to_key: Flow-position mapping to biosphere flow keys.
    :type flow_pos_to_key: dict[int, tuple[str, str, str] | str]
    :param species_map: Value for `species_map`.
    :type species_map: dict[object, str]
    :param signs: Value for `signs`.
    :type signs: dict[object, float]
    :returns: Return value.
    :rtype: pd.DataFrame
    :raises ValueError: If an error occurs."""
    if not flow_pos_to_key:
        return pd.DataFrame(index=inv_years)

    def _flow_name(flow_key: object) -> str:
        if isinstance(flow_key, tuple) and len(flow_key) >= 1:
            return str(flow_key[0])
        return str(flow_key)

    def _resolve_mapping(mapping: dict[object, Any], flow_key: object) -> Any:
        name_key = _flow_name(flow_key)
        if name_key in mapping:
            return mapping[name_key]
        return mapping.get(flow_key)

    n_year = len(inv_years)
    species_order: list[str] = []
    species_to_idx: dict[str, int] = {}
    flow_to_species = np.full(n_flow, -1, dtype=int)
    flow_to_sign = np.ones(n_flow, dtype=float)

    for pos, flow_key in flow_pos_to_key.items():
        fair_species = _resolve_mapping(species_map, flow_key)
        if fair_species is None:
            continue
        specie = str(fair_species)
        if specie not in species_to_idx:
            species_to_idx[specie] = len(species_order)
            species_order.append(specie)
        flow_to_species[int(pos)] = species_to_idx[specie]
        sign_value = _resolve_mapping(signs, flow_key)
        flow_to_sign[int(pos)] = 1.0 if sign_value is None else float(sign_value)

    if not species_order:
        return pd.DataFrame(index=inv_years)

    coo = inv_data

    flow_idx = coo.coords[0].astype(int, copy=False)
    year_idx = coo.coords[1].astype(int, copy=False)
    vals = coo.data.astype(float, copy=False)

    specie_idx = flow_to_species[flow_idx]
    valid = specie_idx >= 0
    agg = np.zeros((len(species_order), n_year), dtype=float)
    if np.any(valid):
        specie_valid = specie_idx[valid]
        year_valid = year_idx[valid]
        vals_valid = vals[valid]
        signs_valid = flow_to_sign[flow_idx[valid]]
        signed_vals = np.where(
            signs_valid == 1.0, vals_valid, np.abs(vals_valid) * signs_valid
        )
        flat_idx = specie_valid * n_year + year_valid
        flat = agg.reshape(-1)
        np.add.at(flat, flat_idx, signed_vals)

    out = pd.DataFrame(agg.T, index=inv_years, columns=species_order)
    out.index.name = "year"
    return out


def _run_fair_emissions(
    emissions_df: pd.DataFrame,
    scenario: str,
    *,
    config_csv: str | Path | None = DEFAULT_CONFIGS_CSV,
    properties_csv: str | Path | None = DEFAULT_PROPERTIES_CSV,
    config_name: str | None = None,
    config_names: list[str] | None = None,
    ghg_method: str | None = "myhre1998",
    temperature_prescribed: bool | None = None,
    debug: bool = False,
    progress: bool = False,
) -> fair.FAIR:
    """run fair emissions.

    :param emissions_df: Value for `emissions_df`.
    :type emissions_df: pd.DataFrame
    :param scenario: Value for `scenario`.
    :type scenario: str
    :param config_csv: Value for `config_csv`.
    :type config_csv: str | Path | None
    :param properties_csv: Value for `properties_csv`.
    :type properties_csv: str | Path | None
    :param config_name: Value for `config_name`.
    :type config_name: str | None
    :param config_names: Value for `config_names`.
    :type config_names: list[str] | None
    :param ghg_method: Value for `ghg_method`.
    :type ghg_method: str | None
    :param temperature_prescribed: Value for `temperature_prescribed`.
    :type temperature_prescribed: bool | None
    :param debug: Value for `debug`.
    :type debug: bool
    :param progress: Value for `progress`.
    :type progress: bool
    :returns: Return value.
    :rtype: fair.FAIR
    :raises ValueError: If an error occurs."""
    df = _normalize_emissions_columns(emissions_df)
    df = df[(df["scenario"] == scenario) & (df["region"].str.lower() == "world")].copy()
    if df.empty:
        raise ValueError(f"Scenario '{scenario}' not found in emissions data.")

    year_cols, year_vals = _extract_year_columns(df)
    if not year_cols:
        raise ValueError("No year columns found in emissions data.")

    meta_cols = ["scenario", "region", "variable", "unit"]
    df = df[list(meta_cols) + year_cols]
    df = _sanitize_emissions_year_values(df)
    years = year_vals
    start_year = min(years)
    end_year = max(years)

    species = sorted(df["variable"].unique())
    if ("CO2 FFI" in species or "CO2 AFOLU" in species) and "CO2" not in species:
        species.append("CO2")

        if properties_csv is None:
            properties_csv = _find_any_fair_repo_file(
                [
                    "species_configs_properties_calibration1.4.1.csv",
                    "species_configs_properties.csv",
                ]
            )
            if properties_csv is None:
                properties_csv = _default_properties_file()
    properties_csv = Path(properties_csv)

    if config_csv is None:
        config_csv = _find_any_fair_repo_file(
            [
                "calibrated_constrained_parameters_calibration1.4.1.csv",
                "calibrated_constrained_parameters.csv",
            ]
        )
        if config_csv is None:
            raise ValueError(
                "config_csv not provided and calibration CSV not found.\n"
                "Paths tried:\n"
                + "\n".join(
                    str(p)
                    for name in [
                        "calibrated_constrained_parameters_calibration1.4.1.csv",
                        "calibrated_constrained_parameters.csv",
                    ]
                    for p in _candidate_paths(name)
                )
            )
    config_csv = Path(config_csv)
    cfg: pd.DataFrame | None = None

    def _load_cfg() -> pd.DataFrame:
        """Load calibration configs only when name normalization requires them."""
        nonlocal cfg
        if cfg is None:
            cfg = pd.read_csv(config_csv, index_col=0)
            if cfg.empty:
                raise ValueError(f"No configs found in {config_csv}.")
        return cfg

    def _normalize_config_name(name: object) -> object:
        """normalize config name.

        :param name: Value for `name`.
        :type name: object
        :returns: Return value.
        :rtype: object"""
        cfg_index = _load_cfg().index
        if name in cfg_index:
            return name
        try:
            name_str = str(name)
            if name_str in cfg_index:
                return name_str
        except Exception:
            pass
        try:
            name_int = int(name)
            if name_int in cfg_index:
                return name_int
        except Exception:
            pass
        return name

    if config_names is not None:
        if config_csv.exists():
            configs = [_normalize_config_name(c) for c in config_names]
        else:
            configs = list(config_names)
    else:
        cfg = _load_cfg()
        if config_name is None:
            config_name = cfg.index[0]
        configs = [_normalize_config_name(config_name)]

    template = _build_fair_template_cached(
        scenario=str(scenario),
        start_year=float(start_year),
        end_year=float(end_year),
        config_csv=str(config_csv.resolve()),
        properties_csv=str(properties_csv.resolve()),
        species_key=tuple(str(s) for s in species),
        configs_key=tuple(configs),
        ghg_method=ghg_method,
        temperature_prescribed=temperature_prescribed,
    )
    f = copy.deepcopy(template)
    _fill_emissions_from_df_fast(
        f,
        df,
        scenario=str(scenario),
        year_cols=year_cols,
        year_vals=year_vals,
    )
    if getattr(f, "_trails_prepared_setup", False):
        _run_prepared_fair(f, progress=progress)
    else:
        with _FAIR_RUN_LOCK:
            _run_prepared_fair(f, progress=progress)
    if not np.isfinite(f.forcing.values).any():
        forcing = _compute_ghg_forcing_from_concentration(f)
        if forcing is not None:
            forcing_array = f.forcing.data
            forcing_array[..., f._ghg_indices] = forcing[..., f._ghg_indices]
            f.forcing.data = forcing_array
    return f


def _compute_ghg_forcing_from_concentration(f: fair.FAIR) -> np.ndarray | None:
    """compute ghg forcing from concentration.

    :param f: Value for `f`.
    :type f: fair.FAIR
    :returns: Return value.
    :rtype: np.ndarray | None"""
    try:
        from fair.forcing import ghg as ghg_mod
    except Exception:
        return None

    conc = np.asarray(f.concentration.values, dtype=float)
    base_arr = np.asarray(
        f.species_configs["baseline_concentration"].values, dtype=float
    )
    ref_arr = np.asarray(
        f.species_configs["forcing_reference_concentration"].values, dtype=float
    )
    scale_arr = np.asarray(f.species_configs["forcing_scale"].values, dtype=float)
    rad_eff = np.asarray(
        f.species_configs["greenhouse_gas_radiative_efficiency"].values, dtype=float
    )
    base_full = base_arr[None, None, ...] * np.ones(
        (1, f._n_scenarios, f._n_configs, f._n_species)
    )
    ref_full = ref_arr[None, None, ...] * np.ones(
        (1, f._n_scenarios, f._n_configs, f._n_species)
    )
    scale_full = scale_arr[None, None, ...] * np.ones(
        (1, f._n_scenarios, f._n_configs, f._n_species)
    )
    rad_full = rad_eff[None, None, ...] * np.ones(
        (1, f._n_scenarios, f._n_configs, f._n_species)
    )

    method = getattr(f, "ghg_method", "meinshausen2020")
    if method == "leach2021":
        return ghg_mod.leach2021ghg(
            conc,
            base_full,
            scale_full,
            rad_full,
            f._co2_indices,
            f._ch4_indices,
            f._n2o_indices,
            f._minor_ghg_indices,
        )
    if method == "meinshausen2020":
        return ghg_mod.meinshausen2020(
            conc,
            ref_full,
            scale_full,
            rad_full,
            f._co2_indices,
            f._ch4_indices,
            f._n2o_indices,
            f._minor_ghg_indices,
        )
    if method == "etminan2016":
        return ghg_mod.etminan2016(
            conc,
            base_full,
            scale_full,
            rad_full,
            f._co2_indices,
            f._ch4_indices,
            f._n2o_indices,
            f._minor_ghg_indices,
        )
    if method == "myhre1998":
        return ghg_mod.myhre1998(
            conc,
            base_full,
            scale_full,
            rad_full,
            f._co2_indices,
            f._ch4_indices,
            f._n2o_indices,
            f._minor_ghg_indices,
        )
    return None


def _extract_fair_timeseries(da: xr.DataArray) -> np.ndarray:
    """extract fair timeseries.

    :param da: Value for `da`.
    :type da: xr.DataArray
    :returns: Return value.
    :rtype: np.ndarray"""
    if "timebounds" in da.dims:
        time_dim = "timebounds"
    elif "time" in da.dims:
        time_dim = "time"
    else:
        time_dim = da.dims[0]

    if "layer" in da.dims and da.dims != (time_dim,):
        layer_vals = da.coords.get("layer", None)
        if layer_vals is not None:
            layers = [str(x) for x in layer_vals.values.tolist()]
            if "surface" in layers:
                da = da.sel(layer="surface")
            else:
                da = da.isel(layer=0)
        else:
            da = da.isel(layer=0)

    # Drop any remaining non-time dims by selecting first index
    for d in list(da.dims):
        if d == time_dim:
            continue
        da = da.isel({d: 0})

    if da.dims != (time_dim,):
        da = da.transpose(time_dim)

    return np.asarray(da.values, dtype=float)


def _extract_fair_timeseries_by_config(da: xr.DataArray) -> np.ndarray:
    """extract fair timeseries by config.

    :param da: Value for `da`.
    :type da: xr.DataArray
    :returns: Return value.
    :rtype: np.ndarray"""
    if "timebounds" in da.dims:
        time_dim = "timebounds"
    elif "time" in da.dims:
        time_dim = "time"
    else:
        time_dim = da.dims[0]

    if "layer" in da.dims and da.dims != (time_dim,):
        layer_vals = da.coords.get("layer", None)
        if layer_vals is not None:
            layers = [str(x) for x in layer_vals.values.tolist()]
            if "surface" in layers:
                da = da.sel(layer="surface")
            else:
                da = da.isel(layer=0)
        else:
            da = da.isel(layer=0)

    for d in list(da.dims):
        if d in (time_dim, "config"):
            continue
        da = da.isel({d: 0})

    if "config" not in da.dims:
        da = da.expand_dims({"config": [0]})

    if da.dims != ("config", time_dim):
        da = da.transpose("config", time_dim)

    return np.asarray(da.values, dtype=float)


def _total_response_by_config(
    da: xr.DataArray,
    *,
    scenario: str,
    sum_species: bool,
) -> xr.DataArray:
    """Return a FaIR response as config x time DataArray."""
    if "scenario" in da.dims:
        da = da.sel(scenario=scenario)

    if "timebounds" in da.dims:
        time_dim = "timebounds"
    elif "time" in da.dims:
        time_dim = "time"
    else:
        time_dim = da.dims[0]

    if "layer" in da.dims and da.dims != (time_dim,):
        layer_vals = da.coords.get("layer", None)
        if layer_vals is not None:
            layers = [str(x) for x in layer_vals.values.tolist()]
            if "surface" in layers:
                da = da.sel(layer="surface")
            else:
                da = da.isel(layer=0)
        else:
            da = da.isel(layer=0)

    if sum_species and "specie" in da.dims:
        da = da.fillna(0.0).sum(dim="specie")

    for dim in list(da.dims):
        if dim in ("config", time_dim):
            continue
        da = da.isel({dim: 0})

    if "config" not in da.dims:
        da = da.expand_dims({"config": [0]})

    if da.dims != ("config", time_dim):
        da = da.transpose("config", time_dim)
    return da


def _all_config_names_from_csv(
    config_csv: str | Path | None,
    config_name: str | None,
    config_names: list[str] | None,
) -> list[str] | None:
    """Resolve default FaIR ensemble configs when no specific config is requested."""
    if config_names is not None or config_name is not None:
        return config_names

    cfg_path = Path(config_csv) if config_csv is not None else None
    if cfg_path is None:
        cfg_path = _find_any_fair_repo_file(
            [
                "calibrated_constrained_parameters_calibration1.4.1.csv",
                "calibrated_constrained_parameters.csv",
            ]
        )
    if cfg_path is None:
        return None
    cfg = pd.read_csv(cfg_path, index_col=0)
    return [str(x) for x in cfg.index.tolist()]


def _inventory_delta_by_fair_species_for_trails(
    trails: Any,
    species_map: dict[object, str],
    signs: dict[object, float],
) -> tuple[pd.DataFrame, list[int]]:
    """Map a Trails inventory to annual kg perturbations by FaIR species."""
    inv = trails.inventory
    if inv is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    if "root activity" not in inv.dims:
        raise ValueError(
            "Trails.inventory must include 'root activity' for CO2 pulse equivalents."
        )

    dims = list(inv.dims)
    if "activity" not in dims or "flow" not in dims or "year" not in dims:
        raise ValueError("Trails.inventory must include activity/flow/year dimensions.")
    inv_data = _inventory_flow_year_root(inv)

    flow_coord = inv.coords["flow"].values
    coord_value_to_pos = {int(v): i for i, v in enumerate(flow_coord)}

    def _is_flow_category_mappable(compartment: str, subcompartment: str) -> bool:
        comp = str(compartment).strip().lower()
        sub = str(subcompartment).strip().lower()
        if comp == "air":
            return True
        if comp == "natural resource" and sub in {"air", "in air"}:
            return True
        return False

    flow_pos_to_key: dict[int, tuple[str, str, str] | str] = {}
    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        for fid, md in meta.items():
            if not isinstance(md, dict):
                continue
            name = md.get("name")
            if not name:
                continue
            compartment = (md.get("compartment") or "").strip()
            subcompartment = (md.get("subcompartment") or "").strip()
            if not _is_flow_category_mappable(compartment, subcompartment):
                continue
            flow_key = (str(name), compartment, subcompartment)
            key = int(fid)
            if key in coord_value_to_pos:
                pos = coord_value_to_pos[key]
            elif 0 <= key < len(flow_coord):
                pos = key
            else:
                continue
            flow_pos_to_key.setdefault(pos, flow_key)

    inv_years = [int(y) for y in inv.coords["year"].values.tolist()]
    delta_by_species = _inventory_emissions_by_fair_species(
        inv_data,
        inv_years,
        int(inv.sizes["flow"]),
        flow_pos_to_key,
        species_map,
        signs,
    )
    return delta_by_species, inv_years


def run_fair_co2_pulse_equivalents(
    trails: Any,
    *,
    scenario: str,
    emissions_csv: str | Path = DEFAULT_EMISSIONS_CSV,
    mapping_yaml: str | Path = DEFAULT_MAPPING_YAML,
    exclude_aerosol_ozone_precursors: bool = False,
    config_csv: str | Path | None = DEFAULT_CONFIGS_CSV,
    properties_csv: str | Path | None = DEFAULT_PROPERTIES_CSV,
    config_name: str | None = None,
    config_names: list[str] | None = None,
    ghg_method: str | None = "myhre1998",
    temperature_prescribed: bool | None = False,
    scale_factor: float | None = None,
    scale_target_fraction: float = 0.01,
    reference_pulse_year: int = 2050,
    window_start: int = 2050,
    window_end: int = 2100,
    reference_pulse_mass_kg: float = 1.0e9,
    co2_species_name: str = "CO2 FFI",
) -> dict[str, Any]:
    """Run FaIR and calculate fixed-window CO2 pulse-equivalent indicators.

    This performs three scenario-consistent FaIR runs: baseline, baseline plus
    the Trails inventory perturbation, and baseline plus a reference CO2 pulse.
    Equivalents are calculated per FaIR configuration first, then summarized.

    :param trails: Trails object with a stored inventory.
    :type trails: Any
    :param scenario: FaIR background scenario.
    :type scenario: str
    :param emissions_csv: Baseline emissions CSV.
    :type emissions_csv: str | Path
    :param mapping_yaml: Trails-flow to FaIR-species mapping.
    :type mapping_yaml: str | Path
    :param exclude_aerosol_ozone_precursors: Exclude mapped aerosol and ozone
        precursor species, such as ``Sulfur``, ``NOx``, ``CO``, ``VOC``,
        ``NH3``, ``BC``, and ``OC``, from the Trails perturbation.
    :type exclude_aerosol_ozone_precursors: bool
    :param config_csv: FaIR calibration config CSV.
    :type config_csv: str | Path | None
    :param properties_csv: FaIR species properties CSV.
    :type properties_csv: str | Path | None
    :param config_name: Optional single FaIR config name.
    :type config_name: str | None
    :param config_names: Optional FaIR config names.
    :type config_names: list[str] | None
    :param ghg_method: FaIR greenhouse-gas forcing method.
    :type ghg_method: str | None
    :param temperature_prescribed: Passed to ``fair.FAIR``.
    :type temperature_prescribed: bool | None
    :param scale_factor: Optional perturbation scaling factor.
    :type scale_factor: float | None
    :param scale_target_fraction: Baseline-relative automatic scaling target.
    :type scale_target_fraction: float
    :param reference_pulse_year: Calendar year of the reference CO2 pulse.
    :type reference_pulse_year: int
    :param window_start: First year included in the assessment window.
    :type window_start: int
    :param window_end: Last year included in the assessment window.
    :type window_end: int
    :param reference_pulse_mass_kg: Size of reference pulse in kg CO2.
    :type reference_pulse_mass_kg: float
    :param co2_species_name: Preferred CO2 emissions row for the reference pulse.
    :type co2_species_name: str
    :returns: Structured CO2 pulse-equivalent result.
    :rtype: dict[str, Any]
    """
    debug = bool(getattr(trails, "debug", False))
    if scale_factor is None:
        if scale_target_fraction <= 0:
            raise ValueError("scale_target_fraction must be > 0.")
    elif scale_factor <= 0:
        raise ValueError("scale_factor must be > 0.")

    df = load_emissions_csv(emissions_csv)
    species_map, signs = load_species_mapping(mapping_yaml)
    if exclude_aerosol_ozone_precursors:
        species_map, signs = _exclude_aerosol_ozone_precursor_mapping(
            species_map,
            signs,
        )
    config_names = _all_config_names_from_csv(config_csv, config_name, config_names)

    delta_by_species, inv_years = _inventory_delta_by_fair_species_for_trails(
        trails,
        species_map,
        signs,
    )
    no_perturbation = delta_by_species.empty
    if no_perturbation and scale_factor is None:
        scale_factor = 1.0

    if not no_perturbation:
        df = _ensure_response_species_rows(
            df,
            scenario=scenario,
            drivers=[str(s) for s in delta_by_species.columns.tolist()],
            debug=debug,
        )

    f_base = _run_fair_emissions(
        df,
        scenario,
        config_csv=config_csv,
        properties_csv=properties_csv,
        config_name=config_name,
        config_names=config_names,
        ghg_method=ghg_method,
        temperature_prescribed=temperature_prescribed,
        debug=debug,
        progress=False,
    )

    base_year_cols, base_year_vals = _extract_year_columns(df)
    has_half_years = any(abs(v - round(v)) > 1e-9 for v in base_year_vals)

    def _year_col_name(year: int) -> str:
        if has_half_years:
            return f"{year + 0.5:.1f}"
        return str(int(year))

    if scale_factor is None and not delta_by_species.empty:
        df_base = df[
            (df["scenario"] == scenario) & (df["region"].str.lower() == "world")
        ].copy()
        df_base = _normalize_emissions_columns(df_base)
        candidates = []
        for specie in delta_by_species.columns:
            rows = df_base[df_base["variable"] == specie]
            if rows.empty:
                continue
            unit = rows["unit"].iloc[0]
            base_row = rows.iloc[0]
            for year, val in delta_by_species[specie].items():
                ycol = _year_col_name(int(year))
                if ycol not in df_base.columns:
                    continue
                base_val = base_row[ycol]
                if base_val == 0 or pd.isna(base_val):
                    continue
                delta_unit = _convert_kg_to_unit(np.array([val]), unit)[0]
                if delta_unit == 0:
                    continue
                candidates.append(
                    scale_target_fraction * abs(base_val) / abs(delta_unit)
                )
        if candidates:
            scale_factor = float(min(candidates))
        else:
            scale_factor = 1.0

    if scale_factor != 1.0 and not delta_by_species.empty:
        delta_by_species = delta_by_species * float(scale_factor)

    def _build_lca_perturbed_df(base_df: pd.DataFrame) -> pd.DataFrame:
        df_pert_local = _normalize_emissions_columns(base_df)
        df_pert_local = df_pert_local[
            (df_pert_local["scenario"] == scenario)
            & (df_pert_local["region"].str.lower() == "world")
        ].copy()
        if no_perturbation:
            return df_pert_local
        for specie_name in delta_by_species.columns:
            rows = df_pert_local[df_pert_local["variable"] == specie_name]
            if rows.empty:
                continue
            unit = rows["unit"].iloc[0]
            idx = rows.index[0]
            for year, val in delta_by_species[specie_name].items():
                ycol = _year_col_name(int(year))
                if ycol not in df_pert_local.columns:
                    continue
                add = _convert_kg_to_unit(np.array([val]), unit)[0]
                df_pert_local.loc[idx, ycol] = df_pert_local.loc[idx, ycol] + add
        return df_pert_local

    f_lca = _run_fair_emissions(
        _build_lca_perturbed_df(df),
        scenario,
        config_csv=config_csv,
        properties_csv=properties_csv,
        config_name=config_name,
        config_names=config_names,
        ghg_method=ghg_method,
        temperature_prescribed=temperature_prescribed,
        debug=debug,
        progress=False,
    )
    f_ref = _run_fair_emissions(
        make_reference_co2_pulse_emissions(
            df,
            scenario=scenario,
            pulse_year=reference_pulse_year,
            pulse_mass_kg=reference_pulse_mass_kg,
            co2_species_name=co2_species_name,
        ),
        scenario,
        config_csv=config_csv,
        properties_csv=properties_csv,
        config_name=config_name,
        config_names=config_names,
        ghg_method=ghg_method,
        temperature_prescribed=temperature_prescribed,
        debug=debug,
        progress=False,
    )

    rf_base = _total_response_by_config(
        f_base.forcing,
        scenario=scenario,
        sum_species=True,
    )
    rf_lca = _total_response_by_config(
        f_lca.forcing,
        scenario=scenario,
        sum_species=True,
    )
    rf_ref = _total_response_by_config(
        f_ref.forcing,
        scenario=scenario,
        sum_species=True,
    )
    temp_base = _total_response_by_config(
        f_base.temperature,
        scenario=scenario,
        sum_species=False,
    )
    temp_lca = _total_response_by_config(
        f_lca.temperature,
        scenario=scenario,
        sum_species=False,
    )
    temp_ref = _total_response_by_config(
        f_ref.temperature,
        scenario=scenario,
        sum_species=False,
    )

    time_dim = "timebounds" if "timebounds" in rf_base.dims else "time"
    years = np.asarray(rf_base.coords[time_dim].values, dtype=float)
    dRF_lca = np.asarray((rf_lca - rf_base).values, dtype=float)
    dT_lca = np.asarray((temp_lca - temp_base).values, dtype=float)
    dRF_ref = np.asarray((rf_ref - rf_base).values, dtype=float)
    dT_ref = np.asarray((temp_ref - temp_base).values, dtype=float)

    effective_scale = 1.0 if scale_factor is None else float(scale_factor)
    if effective_scale != 1.0:
        dRF_lca = dRF_lca / effective_scale
        dT_lca = dT_lca / effective_scale

    calculated = calculate_co2_pulse_equivalents(
        dRF_lca,
        dT_lca,
        dRF_ref,
        dT_ref,
        years,
        reference_pulse_mass_kg,
        start_year=window_start,
        end_year=window_end,
        time_axis=-1,
    )
    rf_values = calculated["co2_pulse_equivalent_integrated_rf"]
    temp_values = calculated["co2_pulse_equivalent_integrated_temperature"]
    config_labels = [str(x) for x in rf_base.coords["config"].values.tolist()]
    indicator_unit = "kg CO2 pulse-equivalent"

    indicator = {
        "reference_year": int(reference_pulse_year),
        "window_start": int(window_start),
        "window_end": int(window_end),
        "reference_pulse_mass": float(reference_pulse_mass_kg),
        "mass_unit": "kg CO2",
        "integrated_rf": {
            **_summarize_indicator_values(rf_values),
            "label": (
                "Integrated RF CO2 pulse equivalent, "
                f"{int(window_start)}-{int(window_end)}"
            ),
            "unit": indicator_unit,
            "by_config": _array_to_builtin(rf_values),
        },
        "integrated_temperature": {
            **_summarize_indicator_values(temp_values),
            "label": (
                "Integrated temperature CO2 pulse equivalent, "
                f"{int(window_start)}-{int(window_end)}"
            ),
            "unit": indicator_unit,
            "by_config": _array_to_builtin(temp_values),
        },
        "diagnostics": {
            "config": config_labels,
            "reference_pulse_mass_kg": float(reference_pulse_mass_kg),
            "scaling_factor": effective_scale,
            "integrated_rf_lca": _array_to_builtin(calculated["integrated_rf_lca"]),
            "integrated_rf_ref": _array_to_builtin(calculated["integrated_rf_ref"]),
            "integrated_temperature_lca": _array_to_builtin(
                calculated["integrated_temperature_lca"]
            ),
            "integrated_temperature_ref": _array_to_builtin(
                calculated["integrated_temperature_ref"]
            ),
            "rf_unit": "W m-2 yr",
            "temperature_unit": "K yr",
        },
    }
    result = {"co2_pulse_equivalent": indicator}
    setattr(trails, "co2_pulse_equivalent", indicator)
    return result


def run_fair_delta_rf(
    trails: Any,
    *,
    scenario: str,
    emissions_csv: str | Path = DEFAULT_EMISSIONS_CSV,
    mapping_yaml: str | Path = DEFAULT_MAPPING_YAML,
    exclude_aerosol_ozone_precursors: bool = False,
    config_csv: str | Path | None = DEFAULT_CONFIGS_CSV,
    properties_csv: str | Path | None = DEFAULT_PROPERTIES_CSV,
    config_name: str | None = None,
    config_names: list[str] | None = None,
    ghg_method: str | None = "myhre1998",
    temperature_prescribed: bool | None = False,
    scale_factor: float | None = None,
    scale_target_fraction: float = 0.01,
    scaling_factor: float | None = None,
    validate_emissions_delta: bool = True,
    validate_atol: float = 1e-6,
    validate_rtol: float = 1e-6,
    validate_raise: bool = False,
    per_species_runs: bool = True,
    per_species_workers: int | None = None,
    quantiles: list[float] | None = None,
    show_progress: bool = True,
) -> xr.DataArray:
    """Run fair delta rf.

    :param trails: Value for `trails`.
    :type trails: Any
    :param scenario: Value for `scenario`.
    :type scenario: str
    :param emissions_csv: Value for `emissions_csv`.
    :type emissions_csv: str | Path
    :param mapping_yaml: Value for `mapping_yaml`.
    :type mapping_yaml: str | Path
    :param exclude_aerosol_ozone_precursors: Exclude mapped aerosol and ozone
        precursor species, such as ``Sulfur``, ``NOx``, ``CO``, ``VOC``,
        ``NH3``, ``BC``, and ``OC``, from the Trails perturbation.
    :type exclude_aerosol_ozone_precursors: bool
    :param config_csv: Value for `config_csv`.
    :type config_csv: str | Path | None
    :param properties_csv: Value for `properties_csv`.
    :type properties_csv: str | Path | None
    :param config_name: Value for `config_name`.
    :type config_name: str | None
    :param config_names: Value for `config_names`.
    :type config_names: list[str] | None
    :param ghg_method: Value for `ghg_method`.
    :type ghg_method: str | None
    :param temperature_prescribed: Value for `temperature_prescribed`.
    :type temperature_prescribed: bool | None
    :param scale_factor: Value for `scale_factor`.
    :type scale_factor: float | None
    :param scale_target_fraction: Value for `scale_target_fraction`.
    :type scale_target_fraction: float
    :param scaling_factor: Value for `scaling_factor`.
    :type scaling_factor: float | None
    :param validate_emissions_delta: Value for `validate_emissions_delta`.
    :type validate_emissions_delta: bool
    :param validate_atol: Value for `validate_atol`.
    :type validate_atol: float
    :param validate_rtol: Value for `validate_rtol`.
    :type validate_rtol: float
    :param validate_raise: Value for `validate_raise`.
    :type validate_raise: bool
    :param per_species_runs: Value for `per_species_runs`.
    :type per_species_runs: bool
    :param per_species_workers: Number of workers for per-species FaIR runs.
        If ``None``, an automatic worker count is used.
    :type per_species_workers: int | None
    :param quantiles: Value for `quantiles`.
    :type quantiles: list[float] | None
    :param show_progress: Show inventory-preparation progress for disk-backed
        inventories.
    :type show_progress: bool
    :returns: Return value.
    :rtype: xr.DataArray
    :raises ValueError: If an error occurs."""
    debug = bool(getattr(trails, "debug", False))
    if scaling_factor is not None:
        scale_factor = float(scaling_factor)
    if scale_factor is None:
        if scale_target_fraction <= 0:
            raise ValueError("scale_target_fraction must be > 0.")
    elif scale_factor <= 0:
        raise ValueError("scale_factor must be > 0.")
    if per_species_workers is not None and int(per_species_workers) < 1:
        raise ValueError("per_species_workers must be >= 1 when provided.")
    if per_species_workers is not None:
        per_species_workers = int(per_species_workers)
    df = load_emissions_csv(emissions_csv)
    species_map, signs = load_species_mapping(mapping_yaml)
    if exclude_aerosol_ozone_precursors:
        species_map, signs = _exclude_aerosol_ozone_precursor_mapping(
            species_map,
            signs,
        )
    if quantiles is None:
        quantiles = [2.5, 25, 50, 75, 97.5]
    quantiles = [float(q) for q in quantiles]
    if config_names is None and config_name is None:
        cfg_path = Path(config_csv) if config_csv is not None else None
        if cfg_path is None:
            cfg_path = _find_any_fair_repo_file(
                [
                    "calibrated_constrained_parameters_calibration1.4.1.csv",
                    "calibrated_constrained_parameters.csv",
                ]
            )
        if cfg_path is not None:
            cfg = pd.read_csv(cfg_path, index_col=0)
            config_names = [str(x) for x in cfg.index.tolist()]

    inv = trails.inventory
    if inv is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    if "root activity" not in inv.dims:
        raise ValueError("Trails.inventory must include 'root activity' for delta RF.")

    dims = list(inv.dims)
    if "activity" not in dims or "flow" not in dims or "year" not in dims:
        raise ValueError("Trails.inventory must include activity/flow/year dimensions.")
    n_flow = int(inv.sizes["flow"])
    n_root = int(inv.sizes["root activity"])

    flow_coord = inv.coords["flow"].values
    coord_value_to_pos = {int(v): i for i, v in enumerate(flow_coord)}

    def _is_flow_category_mappable(compartment: str, subcompartment: str) -> bool:
        """Restrict inventory-to-FaIR mapping to atmospheric exchanges."""
        comp = str(compartment).strip().lower()
        sub = str(subcompartment).strip().lower()
        if comp == "air":
            return True
        if comp == "natural resource" and sub in {"air", "in air"}:
            return True
        return False

    flow_pos_to_key: dict[int, tuple[str, str, str] | str] = {}
    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        for fid, md in meta.items():
            if not isinstance(md, dict):
                continue
            name = md.get("name")
            if not name:
                continue
            compartment = (md.get("compartment") or "").strip()
            subcompartment = (md.get("subcompartment") or "").strip()
            if not _is_flow_category_mappable(compartment, subcompartment):
                continue
            flow_key = (str(name), compartment, subcompartment)
            key = int(fid)
            if key in coord_value_to_pos:
                pos = coord_value_to_pos[key]
            elif 0 <= key < len(flow_coord):
                pos = key
            else:
                continue
            flow_pos_to_key.setdefault(pos, flow_key)

    inv_years = [int(y) for y in inv.coords["year"].values.tolist()]

    mapped_flow_positions = sorted(
        int(position)
        for position, flow_key in flow_pos_to_key.items()
        if _resolve_flow_mapping(species_map, flow_key) is not None
    )
    # Reduce only FaIR-mapped atmospheric flows and reuse this compact COO for
    # species deltas and root allocations. The chunked builder streams directly
    # from its finalized runs and caches the result for subsequent FaIR calls.
    inv_data = _inventory_flow_year_root(
        inv,
        flow_positions=mapped_flow_positions,
        builder=getattr(trails, "_inventory_builder", None),
        show_progress=bool(show_progress),
    )

    # Build perturbations from activity-reduced inventory (single sparse
    # reduction path).
    delta_by_species = _inventory_emissions_by_fair_species(
        inv_data,
        inv_years,
        n_flow,
        flow_pos_to_key,
        species_map,
        signs,
    )
    no_perturbation = delta_by_species.empty
    if no_perturbation and scale_factor is None:
        scale_factor = 1.0

    # Precursor species (NOx, VOC, etc.) force climate through calculated
    # response channels (Ozone, ARI/ACI, ...). Add zero rows so FaIR includes
    # these channels in the run when needed.
    if not no_perturbation:
        df = _ensure_response_species_rows(
            df,
            scenario=scenario,
            drivers=[str(s) for s in delta_by_species.columns.tolist()],
            debug=debug,
        )

    # Baseline run
    f_base = _run_fair_emissions(
        df,
        scenario,
        config_csv=config_csv,
        properties_csv=properties_csv,
        config_name=config_name,
        config_names=config_names,
        ghg_method=ghg_method,
        temperature_prescribed=temperature_prescribed,
        debug=debug,
        progress=False,
    )

    # Determine if the emissions baseline uses half-year columns (e.g., 1750.5)
    base_year_cols, base_year_vals = _extract_year_columns(df)
    has_half_years = any(abs(v - round(v)) > 1e-9 for v in base_year_vals)

    def _year_col_name(year: int) -> str:
        """year col name.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: str"""
        if has_half_years:
            return f"{year + 0.5:.1f}"
        return str(int(year))

    if scale_factor is None and not delta_by_species.empty:
        df_base = df[
            (df["scenario"] == scenario) & (df["region"].str.lower() == "world")
        ].copy()
        df_base = _normalize_emissions_columns(df_base)
        candidates = []
        for specie in delta_by_species.columns:
            rows = df_base[df_base["variable"] == specie]
            if rows.empty:
                continue
            unit = rows["unit"].iloc[0]
            base_row = rows.iloc[0]
            for year, val in delta_by_species[specie].items():
                ycol = _year_col_name(int(year))
                if ycol not in df_base.columns:
                    continue
                base_val = base_row[ycol]
                if base_val == 0 or pd.isna(base_val):
                    continue
                delta_unit = _convert_kg_to_unit(np.array([val]), unit)[0]
                if delta_unit == 0:
                    continue
                candidates.append(
                    scale_target_fraction * abs(base_val) / abs(delta_unit)
                )
        if candidates:
            scale_factor = float(min(candidates))
        else:
            scale_factor = 1.0

    if scale_factor != 1.0 and not delta_by_species.empty:
        delta_by_species = delta_by_species * float(scale_factor)

    def _build_perturbed_df(
        base_df: pd.DataFrame,
        specie: str | None,
        delta_series: pd.Series | None = None,
    ) -> pd.DataFrame:
        """build perturbed df.

        :param base_df: Value for `base_df`.
        :type base_df: pd.DataFrame
        :param specie: Value for `specie`.
        :type specie: str | None
        :param delta_series: Value for `delta_series`.
        :type delta_series: pd.Series | None
        :returns: Return value.
        :rtype: pd.DataFrame"""
        df_pert_local = _normalize_emissions_columns(base_df)
        df_pert_local = df_pert_local[
            (df_pert_local["scenario"] == scenario)
            & (df_pert_local["region"].str.lower() == "world")
        ].copy()
        if specie is None:
            for specie_name in delta_by_species.columns:
                rows = df_pert_local[df_pert_local["variable"] == specie_name]
                if rows.empty:
                    continue
                unit = rows["unit"].iloc[0]
                series = delta_by_species[specie_name]
                for year, val in series.items():
                    ycol = _year_col_name(int(year))
                    if ycol not in df_pert_local.columns:
                        continue
                    idx = rows.index[0]
                    add = _convert_kg_to_unit(np.array([val]), unit)[0]
                    df_pert_local.loc[idx, ycol] = df_pert_local.loc[idx, ycol] + add
            return df_pert_local
        rows = df_pert_local[df_pert_local["variable"] == specie]
        if rows.empty:
            return df_pert_local
        unit = rows["unit"].iloc[0]
        series = delta_series if delta_series is not None else delta_by_species[specie]
        for year, val in series.items():
            ycol = _year_col_name(int(year))
            if ycol not in df_pert_local.columns:
                continue
            idx = rows.index[0]
            add = _convert_kg_to_unit(np.array([val]), unit)[0]
            df_pert_local.loc[idx, ycol] = df_pert_local.loc[idx, ycol] + add
        return df_pert_local

    if validate_emissions_delta:
        df_base_chk = _normalize_emissions_columns(df)
        df_base_chk = df_base_chk[
            (df_base_chk["scenario"] == scenario)
            & (df_base_chk["region"].str.lower() == "world")
        ].copy()
        df_pert_chk = df_base_chk if no_perturbation else _build_perturbed_df(df, None)
        year_cols, year_vals = _extract_year_columns(df_base_chk)
        if year_cols:
            issues: list[str] = []
            for specie in delta_by_species.columns:
                base_rows = df_base_chk[df_base_chk["variable"] == specie]
                pert_rows = df_pert_chk[df_pert_chk["variable"] == specie]
                if base_rows.empty or pert_rows.empty:
                    issues.append(f"{specie}: missing in IAMC data")
                    continue
                unit = base_rows["unit"].iloc[0]
                base_vals = base_rows[year_cols].astype(float).sum(axis=0).to_numpy()
                pert_vals = pert_rows[year_cols].astype(float).sum(axis=0).to_numpy()
                diff_unit = pert_vals - base_vals
                try:
                    diff_kg = _convert_unit_to_kg(diff_unit, unit)
                except ValueError:
                    issues.append(f"{specie}: unsupported unit '{unit}'")
                    continue
                expected = delta_by_species[specie].copy()
                if has_half_years:
                    expected.index = expected.index + 0.5
                expected = expected.reindex(year_vals, fill_value=0.0)
                expected_kg = expected.to_numpy(dtype=float)
                abs_diff = np.abs(diff_kg - expected_kg)
                denom = np.maximum(np.abs(expected_kg), 1e-12)
                rel_diff = abs_diff / denom
                if np.any((abs_diff > validate_atol) & (rel_diff > validate_rtol)):
                    max_abs = float(np.nanmax(abs_diff))
                    max_rel = float(np.nanmax(rel_diff))
                    issues.append(
                        f"{specie}: max_abs={max_abs:.3e} kg/yr max_rel={max_rel:.3e}"
                    )
                    if debug:
                        idx = int(np.nanargmax(abs_diff))
                        year = int(year_vals[idx]) if idx < len(year_vals) else idx
                        exp_val = float(expected_kg[idx])
                        got_val = float(diff_kg[idx])
                        print(
                            "FAIR debug: worst mismatch",
                            specie,
                            "year",
                            year,
                            "expected",
                            f"{exp_val:.3e}",
                            "got",
                            f"{got_val:.3e}",
                            "abs",
                            f"{abs_diff[idx]:.3e}",
                            "rel",
                            f"{rel_diff[idx]:.3e}",
                        )
            if issues:
                msg = "FAIR emissions delta check failed: " + "; ".join(issues)
                if validate_raise:
                    raise ValueError(msg)
                if debug:
                    print(msg)

    if config_names is not None and len(config_names) > 1:
        forcing_base = f_base.forcing.sel(scenario=scenario)
        temp_base = f_base.temperature.sel(scenario=scenario)
    else:
        forcing_base = f_base.forcing.sel(scenario=scenario, config=f_base.configs[0])
        temp_base = f_base.temperature.sel(scenario=scenario, config=f_base.configs[0])
    fair_years = [float(y) for y in forcing_base.coords["timebounds"].values.tolist()]
    year_to_fair_idx = {y: i for i, y in enumerate(fair_years)}

    # Build per-species sparse entries
    flow_idx = inv_data.coords[0]
    year_idx = inv_data.coords[1]
    root_idx = inv_data.coords[2]
    data = inv_data.data.astype(float, copy=False)

    coords_out = []
    data_out = []
    temp_coords_out = []
    temp_data_out = []
    n_quant = len(quantiles)

    def _flow_name(flow_key: object) -> str:
        if isinstance(flow_key, tuple) and len(flow_key) >= 1:
            return str(flow_key[0])
        return str(flow_key)

    def _resolve_mapping(mapping: dict[object, Any], flow_key: object) -> Any:
        name_key = _flow_name(flow_key)
        if name_key in mapping:
            return mapping[name_key]
        return mapping.get(flow_key)

    species_to_positions: Dict[str, list[int]] = {}
    for pos, flow_key in flow_pos_to_key.items():
        specie = _resolve_mapping(species_map, flow_key)
        if specie is None:
            continue
        species_to_positions.setdefault(specie, []).append(pos)

    rf_alias = {"CO2 FFI": "CO2", "CO2 AFOLU": "CO2"}

    def _append_allocated_rf(
        specie: str,
        rf_series: np.ndarray,
        sign_mode: str,
        *,
        quantile_idx: int,
    ) -> None:
        """append allocated rf.

        :param specie: Value for `specie`.
        :type specie: str
        :param rf_series: Value for `rf_series`.
        :type rf_series: np.ndarray
        :param sign_mode: Value for `sign_mode`.
        :type sign_mode: str
        :param quantile_idx: Value for `quantile_idx`.
        :type quantile_idx: int
        :raises ValueError: If an error occurs."""
        positions = species_to_positions.get(specie, [])
        if not positions:
            return
        positions = [int(p) for p in positions]
        pos_set = set(positions)
        mask = np.isin(flow_idx, list(pos_set))
        if not np.any(mask):
            return

        f_idx = flow_idx[mask]
        y_idx = year_idx[mask]
        r_idx = root_idx[mask]
        vals = data[mask]

        # Apply sign overrides
        if signs:
            sign_arr = np.ones_like(vals)
            for i, p in enumerate(f_idx):
                flow_key = flow_pos_to_key.get(int(p))
                sign_value = _resolve_mapping(signs, flow_key)
                sign_val = 1.0 if sign_value is None else float(sign_value)
                if sign_val != 1.0:
                    sign_arr[i] = sign_val
            vals = np.abs(vals) * sign_arr

        # Map inventory year index -> fair year index
        fair_idx = []
        for yi in y_idx:
            y = inv_years[int(yi)]
            y_key = float(y) + (0.5 if has_half_years else 0.0)
            fi = year_to_fair_idx.get(y_key)
            if fi is None:
                fair_idx.append(-1)
            else:
                fair_idx.append(fi)
        fair_idx = np.asarray(fair_idx, dtype=int)
        valid = fair_idx >= 0
        if not np.any(valid):
            return

        f_idx = f_idx[valid]
        r_idx = r_idx[valid]
        vals = vals[valid]
        fair_idx = fair_idx[valid]

        # Build per-year emissions by (flow, root)
        pair_keys = list({(int(f), int(r)) for f, r in zip(f_idx, r_idx)})
        pair_index = {p: i for i, p in enumerate(pair_keys)}
        n_pairs = len(pair_keys)
        E = np.zeros((n_pairs, len(fair_years)), dtype=float)
        for f, r, fi, v in zip(f_idx, r_idx, fair_idx, vals):
            E[pair_index[(int(f), int(r))], int(fi)] += float(v)

        if sign_mode == "pos":
            E_use = np.maximum(E, 0.0)
        elif sign_mode == "neg":
            E_use = np.minimum(E, 0.0)
        else:
            raise ValueError("sign_mode must be 'pos' or 'neg'.")

        # Allocate by cumulative signed emissions so uptake yields negative tail
        cumE = np.cumsum(E_use, axis=1)
        total = np.sum(cumE, axis=0)
        if not np.any(np.isfinite(total)):
            return
        per_kg = np.zeros_like(total, dtype=float)
        prev = None
        eps = 1e-12
        for yi, denom in enumerate(total):
            if not np.isfinite(denom) or abs(denom) <= eps:
                if prev is not None:
                    per_kg[yi] = prev
                continue
            per_kg[yi] = rf_series[yi] / denom
            prev = per_kg[yi]
        RF_alloc = cumE * per_kg[None, :]

        row_idx, col_idx = np.nonzero(RF_alloc)
        if row_idx.size == 0:
            return
        pairs_arr = np.asarray(pair_keys, dtype=int)
        f_out = pairs_arr[row_idx, 0]
        r_out = pairs_arr[row_idx, 1]
        y_out = col_idx.astype(int)
        q_out = np.full_like(y_out, int(quantile_idx))
        coords_out.append(np.vstack([q_out, y_out, f_out, r_out]))
        data_out.append(RF_alloc[row_idx, col_idx].astype(float, copy=False))

    def _append_allocated_temp(
        specie: str,
        temp_series: np.ndarray,
        sign_mode: str,
        *,
        quantile_idx: int,
    ) -> None:
        """append allocated temp.

        :param specie: Value for `specie`.
        :type specie: str
        :param temp_series: Value for `temp_series`.
        :type temp_series: np.ndarray
        :param sign_mode: Value for `sign_mode`.
        :type sign_mode: str
        :param quantile_idx: Value for `quantile_idx`.
        :type quantile_idx: int
        :raises ValueError: If an error occurs."""
        positions = species_to_positions.get(specie, [])
        if not positions:
            return
        positions = [int(p) for p in positions]
        pos_set = set(positions)
        mask = np.isin(flow_idx, list(pos_set))
        if not np.any(mask):
            return

        f_idx = flow_idx[mask]
        y_idx = year_idx[mask]
        r_idx = root_idx[mask]
        vals = data[mask]

        if signs:
            sign_arr = np.ones_like(vals)
            for i, p in enumerate(f_idx):
                flow_key = flow_pos_to_key.get(int(p))
                sign_value = _resolve_mapping(signs, flow_key)
                sign_val = 1.0 if sign_value is None else float(sign_value)
                if sign_val != 1.0:
                    sign_arr[i] = sign_val
            vals = np.abs(vals) * sign_arr

        fair_idx = []
        for yi in y_idx:
            y = inv_years[int(yi)]
            y_key = float(y) + (0.5 if has_half_years else 0.0)
            fi = year_to_fair_idx.get(y_key)
            if fi is None:
                fair_idx.append(-1)
            else:
                fair_idx.append(fi)
        fair_idx = np.asarray(fair_idx, dtype=int)
        valid = fair_idx >= 0
        if not np.any(valid):
            return

        f_idx = f_idx[valid]
        r_idx = r_idx[valid]
        vals = vals[valid]
        fair_idx = fair_idx[valid]

        pair_keys = list({(int(f), int(r)) for f, r in zip(f_idx, r_idx)})
        pair_index = {p: i for i, p in enumerate(pair_keys)}
        n_pairs = len(pair_keys)
        E = np.zeros((n_pairs, len(fair_years)), dtype=float)
        for f, r, fi, v in zip(f_idx, r_idx, fair_idx, vals):
            E[pair_index[(int(f), int(r))], int(fi)] += float(v)

        if sign_mode == "pos":
            E_use = np.maximum(E, 0.0)
        elif sign_mode == "neg":
            E_use = np.minimum(E, 0.0)
        else:
            raise ValueError("sign_mode must be 'pos' or 'neg'.")

        cumE = np.cumsum(E_use, axis=1)
        total = np.sum(cumE, axis=0)
        if not np.any(np.isfinite(total)):
            return
        per_kg = np.zeros_like(total, dtype=float)
        prev = None
        eps = 1e-12
        for yi, denom in enumerate(total):
            if not np.isfinite(denom) or abs(denom) <= eps:
                if prev is not None:
                    per_kg[yi] = prev
                continue
            per_kg[yi] = temp_series[yi] / denom
            prev = per_kg[yi]
        TEMP_alloc = cumE * per_kg[None, :]

        row_idx, col_idx = np.nonzero(TEMP_alloc)
        if row_idx.size == 0:
            return
        pairs_arr = np.asarray(pair_keys, dtype=int)
        f_out = pairs_arr[row_idx, 0]
        r_out = pairs_arr[row_idx, 1]
        y_out = col_idx.astype(int)
        q_out = np.full_like(y_out, int(quantile_idx))
        temp_coords_out.append(np.vstack([q_out, y_out, f_out, r_out]))
        temp_data_out.append(TEMP_alloc[row_idx, col_idx].astype(float, copy=False))

    def _compute_species_quantiles(
        specie: str, delta_series: pd.Series
    ) -> tuple[np.ndarray, np.ndarray] | None:
        df_pert = _build_perturbed_df(df, specie, delta_series)
        f_pert = _run_fair_emissions(
            df_pert,
            scenario,
            config_csv=config_csv,
            properties_csv=properties_csv,
            config_name=config_name,
            config_names=config_names,
            ghg_method=ghg_method,
            temperature_prescribed=temperature_prescribed,
            debug=debug,
            progress=False,
        )
        if config_names is not None and len(config_names) > 1:
            forcing_pert = f_pert.forcing.sel(scenario=scenario)
            temp_pert = f_pert.temperature.sel(scenario=scenario)
        else:
            forcing_pert = f_pert.forcing.sel(
                scenario=scenario, config=f_pert.configs[0]
            )
            temp_pert = f_pert.temperature.sel(
                scenario=scenario, config=f_pert.configs[0]
            )
        delta_forcing = forcing_pert - forcing_base
        delta_temp = temp_pert - temp_base
        if debug:
            try:
                arr = np.asarray(delta_temp.values, dtype=float)
                print(
                    "FAIR debug: delta_temp min/max/finite",
                    float(np.nanmin(arr)),
                    float(np.nanmax(arr)),
                    int(np.isfinite(arr).sum()),
                )
            except Exception:
                pass

        effective_scale = 1.0 if scale_factor is None else float(scale_factor)
        if effective_scale != 1.0:
            delta_forcing = delta_forcing / effective_scale
            delta_temp = delta_temp / effective_scale

        available_species = {
            str(s) for s in delta_forcing.coords["specie"].values.tolist()
        }
        primary = rf_alias.get(specie, specie)
        target_species: list[str] = []
        for candidate in [primary, *_PRECURSOR_RESPONSE_SPECIES.get(primary, ())]:
            cand = str(candidate)
            if cand in available_species and cand not in target_species:
                target_species.append(cand)
        if not target_species:
            return None

        if config_names is not None and len(config_names) > 1:
            rf_parts = [
                _extract_fair_timeseries_by_config(delta_forcing.sel(specie=target))
                for target in target_species
            ]
            if len(rf_parts) == 1:
                rf_series = rf_parts[0]
            else:
                rf_series = np.sum(np.stack(rf_parts, axis=0), axis=0)
            rf_quant = _safe_nanpercentile(rf_series, quantiles)
            temp_series = _extract_fair_timeseries_by_config(delta_temp)
            temp_quant = _safe_nanpercentile(temp_series, quantiles)
        else:
            rf_series = np.zeros(len(fair_years), dtype=float)
            for target in target_species:
                part = np.asarray(delta_forcing.sel(specie=target).values, dtype=float)
                part = np.nan_to_num(part, nan=0.0)
                part = np.ravel(part)
                rf_series += part
            rf_quant = np.tile(rf_series[None, :], (n_quant, 1))
            temp_series = _extract_fair_timeseries(delta_temp)
            temp_series = np.nan_to_num(temp_series, nan=0.0)
            temp_quant = np.tile(temp_series[None, :], (n_quant, 1))
        return rf_quant, temp_quant

    if no_perturbation:
        if debug:
            print("FAIR debug: no perturbations; delta RF will be zero.")
    elif per_species_runs:
        delta_by_species_pos = delta_by_species.clip(lower=0.0)
        delta_by_species_neg = delta_by_species.clip(upper=0.0)

        work_items: list[tuple[str, str, pd.Series]] = []
        for specie in delta_by_species.columns.tolist():
            specie_str = str(specie)
            series_pos = delta_by_species_pos[specie]
            if np.any(series_pos.values != 0):
                work_items.append((specie_str, "pos", series_pos.copy()))
            series_neg = delta_by_species_neg[specie]
            if np.any(series_neg.values != 0):
                work_items.append((specie_str, "neg", series_neg.copy()))

        results: dict[tuple[str, str], tuple[np.ndarray, np.ndarray] | None] = {}
        if work_items:
            if per_species_workers is None:
                auto_workers = os.cpu_count() or 1
                max_workers = min(4, auto_workers, len(work_items))
            else:
                max_workers = min(per_species_workers, len(work_items))
            max_workers = max(1, int(max_workers))
            if debug:
                print(
                    "FAIR debug: per-species workers",
                    max_workers,
                    "items",
                    len(work_items),
                )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_key = {
                    executor.submit(
                        _compute_species_quantiles,
                        specie,
                        series,
                    ): (specie, sign_mode)
                    for specie, sign_mode, series in work_items
                }
                with tqdm(
                    total=len(work_items),
                    desc="FaIR per-species",
                    unit="run",
                ) as pbar:
                    for future in as_completed(future_to_key):
                        key = future_to_key[future]
                        results[key] = future.result()
                        pbar.set_postfix_str(f"{key[0]}:{key[1]}")
                        pbar.update(1)

        for specie in delta_by_species.columns.tolist():
            specie_str = str(specie)
            for sign_mode in ("pos", "neg"):
                quantiles_pair = results.get((specie_str, sign_mode))
                if quantiles_pair is None:
                    continue
                rf_quant, temp_quant = quantiles_pair
                for qi in range(n_quant):
                    _append_allocated_rf(
                        specie_str, rf_quant[qi], sign_mode, quantile_idx=qi
                    )
                    _append_allocated_temp(
                        specie_str, temp_quant[qi], sign_mode, quantile_idx=qi
                    )
    else:
        df_pert = _build_perturbed_df(df, None)
        f_pert = _run_fair_emissions(
            df_pert,
            scenario,
            config_csv=config_csv,
            properties_csv=properties_csv,
            config_name=config_name,
            config_names=config_names,
            ghg_method=ghg_method,
            temperature_prescribed=temperature_prescribed,
            debug=debug,
            progress=False,
        )
        if config_names is not None and len(config_names) > 1:
            forcing_pert = f_pert.forcing.sel(scenario=scenario)
            temp_pert = f_pert.temperature.sel(scenario=scenario)
        else:
            forcing_pert = f_pert.forcing.sel(
                scenario=scenario, config=f_pert.configs[0]
            )
            temp_pert = f_pert.temperature.sel(
                scenario=scenario, config=f_pert.configs[0]
            )
        delta_forcing = forcing_pert - forcing_base
        delta_temp = temp_pert - temp_base
        if debug:
            try:
                arr = np.asarray(delta_temp.values, dtype=float)
                print(
                    "FAIR debug: delta_temp min/max/finite",
                    float(np.nanmin(arr)),
                    float(np.nanmax(arr)),
                    int(np.isfinite(arr).sum()),
                )
            except Exception:
                pass
        if scale_factor is None:
            scale_factor = 1.0
        if scale_factor != 1.0:
            delta_forcing = delta_forcing / float(scale_factor)
            delta_temp = delta_temp / float(scale_factor)

        # Fast (single-run) mode can't isolate species-specific temperature from
        # FaIR directly. Allocate total delta temperature to each specie by its
        # yearly RF share so species contributions sum back to total temperature.
        if "specie" in delta_forcing.dims:
            total_forcing = delta_forcing.sum(dim="specie")
        else:
            total_forcing = delta_forcing

        def _rf_share(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
            numer_arr = np.nan_to_num(np.asarray(numer, dtype=float), nan=0.0)
            denom_arr = np.nan_to_num(np.asarray(denom, dtype=float), nan=0.0)
            # Use exact non-zero denominator masking. A fixed absolute cutoff can
            # zero-out physically meaningful tiny signals.
            return np.divide(
                numer_arr,
                denom_arr,
                out=np.zeros_like(numer_arr, dtype=float),
                where=np.isfinite(denom_arr) & (denom_arr != 0.0),
            )

        available_species = set()
        if "specie" in delta_forcing.dims:
            available_species = {
                str(s) for s in delta_forcing.coords["specie"].values.tolist()
            }

        source_species = [str(s) for s in delta_by_species.columns.tolist()]
        if not source_species and "specie" in delta_forcing.dims:
            source_species = [
                str(s) for s in delta_forcing.coords["specie"].values.tolist()
            ]

        source_to_primary = {
            str(specie): str(rf_alias.get(str(specie), str(specie)))
            for specie in source_species
        }
        primary_to_sources: dict[str, list[str]] = {}
        for source_specie, primary in source_to_primary.items():
            primary_to_sources.setdefault(primary, []).append(source_specie)

        def _cumulative_source_emissions(source_specie: str) -> np.ndarray:
            out = np.zeros(len(fair_years), dtype=float)
            if source_specie not in delta_by_species.columns:
                return out
            series = delta_by_species[source_specie]
            for year, val in series.items():
                y_key = float(int(year)) + (0.5 if has_half_years else 0.0)
                fi = year_to_fair_idx.get(y_key)
                if fi is None:
                    continue
                out[int(fi)] += float(val)
            return np.cumsum(out)

        cumulative_by_source = {
            source_specie: _cumulative_source_emissions(source_specie)
            for source_specie in source_species
            if source_specie in delta_by_species.columns
        }

        if config_names is not None and len(config_names) > 1:
            rf_total_series_cfg = _extract_fair_timeseries_by_config(total_forcing)
            temp_total_series_cfg = _extract_fair_timeseries_by_config(delta_temp)
        else:
            rf_total_series_1d = _extract_fair_timeseries(total_forcing)
            rf_total_series_1d = np.nan_to_num(rf_total_series_1d, nan=0.0)
            rf_total_series_1d = np.ravel(rf_total_series_1d)
            temp_total_series_1d = _extract_fair_timeseries(delta_temp)
            temp_total_series_1d = np.nan_to_num(temp_total_series_1d, nan=0.0)
            temp_total_series_1d = np.ravel(temp_total_series_1d)

        for source_specie in source_species:
            primary = source_to_primary.get(source_specie, source_specie)
            target_species: list[str] = []
            for candidate in [primary, *_PRECURSOR_RESPONSE_SPECIES.get(primary, ())]:
                cand = str(candidate)
                if cand in available_species and cand not in target_species:
                    target_species.append(cand)
            if not target_species:
                if source_specie in available_species:
                    target_species = [source_specie]
                else:
                    continue

            if config_names is not None and len(config_names) > 1:
                rf_parts = [
                    _extract_fair_timeseries_by_config(delta_forcing.sel(specie=target))
                    for target in target_species
                ]
                if len(rf_parts) == 1:
                    rf_target = rf_parts[0]
                else:
                    rf_target = np.sum(np.stack(rf_parts, axis=0), axis=0)

                same_sources = [
                    s
                    for s in primary_to_sources.get(primary, [source_specie])
                    if s in cumulative_by_source
                ]
                if len(same_sources) > 1 and source_specie in cumulative_by_source:
                    total_cum = np.sum(
                        np.stack(
                            [cumulative_by_source[s] for s in same_sources], axis=0
                        ),
                        axis=0,
                    )
                    frac = _rf_share(cumulative_by_source[source_specie], total_cum)
                    rf_series = rf_target * frac[None, :]
                else:
                    rf_series = rf_target

                share = _rf_share(rf_series, rf_total_series_cfg)
                temp_series = np.asarray(
                    temp_total_series_cfg, dtype=float
                ) * np.nan_to_num(share, nan=0.0)
                rf_quant = _safe_nanpercentile(rf_series, quantiles)
                temp_quant = _safe_nanpercentile(temp_series, quantiles)
            else:
                rf_target = np.zeros(len(fair_years), dtype=float)
                for target in target_species:
                    part = np.asarray(
                        delta_forcing.sel(specie=target).values, dtype=float
                    )
                    part = np.nan_to_num(part, nan=0.0)
                    rf_target += np.ravel(part)

                same_sources = [
                    s
                    for s in primary_to_sources.get(primary, [source_specie])
                    if s in cumulative_by_source
                ]
                if len(same_sources) > 1 and source_specie in cumulative_by_source:
                    total_cum = np.sum(
                        np.stack(
                            [cumulative_by_source[s] for s in same_sources], axis=0
                        ),
                        axis=0,
                    )
                    frac = _rf_share(cumulative_by_source[source_specie], total_cum)
                    rf_series = rf_target * frac
                else:
                    rf_series = rf_target

                share = _rf_share(rf_series, rf_total_series_1d)
                temp_series = temp_total_series_1d * np.nan_to_num(share, nan=0.0)
                rf_quant = np.tile(rf_series[None, :], (n_quant, 1))
                temp_quant = np.tile(temp_series[None, :], (n_quant, 1))

            has_pos = True
            has_neg = False
            if source_specie in delta_by_species.columns:
                vals = np.asarray(delta_by_species[source_specie].values, dtype=float)
                has_pos = bool(np.any(vals > 0.0))
                has_neg = bool(np.any(vals < 0.0))

            for qi in range(n_quant):
                rf_q = np.asarray(rf_quant[qi], dtype=float)
                temp_q = np.asarray(temp_quant[qi], dtype=float)

                if has_pos and has_neg:
                    rf_pos = np.maximum(rf_q, 0.0)
                    rf_neg = np.minimum(rf_q, 0.0)
                    temp_pos = np.maximum(temp_q, 0.0)
                    temp_neg = np.minimum(temp_q, 0.0)
                    if np.any(rf_pos) or np.any(temp_pos):
                        _append_allocated_rf(
                            source_specie,
                            rf_pos,
                            "pos",
                            quantile_idx=qi,
                        )
                        _append_allocated_temp(
                            source_specie,
                            temp_pos,
                            "pos",
                            quantile_idx=qi,
                        )
                    if np.any(rf_neg) or np.any(temp_neg):
                        _append_allocated_rf(
                            source_specie,
                            rf_neg,
                            "neg",
                            quantile_idx=qi,
                        )
                        _append_allocated_temp(
                            source_specie,
                            temp_neg,
                            "neg",
                            quantile_idx=qi,
                        )
                elif has_neg and not has_pos:
                    _append_allocated_rf(
                        source_specie,
                        rf_q,
                        "neg",
                        quantile_idx=qi,
                    )
                    _append_allocated_temp(
                        source_specie,
                        temp_q,
                        "neg",
                        quantile_idx=qi,
                    )
                else:
                    _append_allocated_rf(
                        source_specie,
                        rf_q,
                        "pos",
                        quantile_idx=qi,
                    )
                    _append_allocated_temp(
                        source_specie,
                        temp_q,
                        "pos",
                        quantile_idx=qi,
                    )

    if coords_out:
        coords = np.hstack(coords_out)
        data = np.concatenate(data_out)
        rf_sparse = sparse.COO(
            coords=coords,
            data=data,
            shape=(n_quant, len(fair_years), n_flow, n_root),
        )
    else:
        rf_sparse = sparse.COO(
            coords=np.zeros((4, 0), dtype=int),
            data=np.array([], dtype=float),
            shape=(n_quant, len(fair_years), n_flow, n_root),
        )

    year_coord = np.array(fair_years, dtype=float)
    trails.instant_radiative_forcing = xr.DataArray(
        rf_sparse,
        dims=("quantile", "year", "flow", "root activity"),
        coords={
            "quantile": np.array(quantiles, dtype=float),
            "year": year_coord,
            "flow": inv.coords["flow"],
            "root activity": inv.coords["root activity"],
        },
    )
    if temp_coords_out:
        tcoords = np.hstack(temp_coords_out)
        tdata = np.concatenate(temp_data_out)
        temp_sparse = sparse.COO(
            coords=tcoords,
            data=tdata,
            shape=(n_quant, len(fair_years), n_flow, n_root),
        )
    else:
        temp_sparse = sparse.COO(
            coords=np.zeros((4, 0), dtype=int),
            data=np.array([], dtype=float),
            shape=(n_quant, len(fair_years), n_flow, n_root),
        )
    trails.delta_temperature = xr.DataArray(
        temp_sparse,
        dims=("quantile", "year", "flow", "root activity"),
        coords={
            "quantile": np.array(quantiles, dtype=float),
            "year": year_coord,
            "flow": inv.coords["flow"],
            "root activity": inv.coords["root activity"],
        },
    )
    return trails.instant_radiative_forcing
