"""FaIR integration helpers."""

from __future__ import annotations

import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import os
from pathlib import Path
from typing import Any, Dict, Tuple

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


def _sanitize_emissions_year_values(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce emissions year columns to numeric and replace missing with zero."""
    year_cols, _ = _extract_year_columns(df)
    if not year_cols:
        return df
    out = df.copy()
    out.loc[:, year_cols] = (
        out.loc[:, year_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    )
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
    initialise(f.temperature, 0)
    initialise(f.forcing, 0)
    return f


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

    n_year = len(inv_years)
    species_order: list[str] = []
    species_to_idx: dict[str, int] = {}
    flow_to_species = np.full(n_flow, -1, dtype=int)
    flow_to_sign = np.ones(n_flow, dtype=float)

    for pos, flow_key in flow_pos_to_key.items():
        fair_species = species_map.get(flow_key)
        if fair_species is None and isinstance(flow_key, tuple):
            fair_species = species_map.get(flow_key[0])
        if fair_species is None:
            continue
        specie = str(fair_species)
        if specie not in species_to_idx:
            species_to_idx[specie] = len(species_order)
            species_order.append(specie)
        flow_to_species[int(pos)] = species_to_idx[specie]
        flow_to_sign[int(pos)] = float(
            signs.get(
                flow_key,
                signs.get(
                    flow_key[0] if isinstance(flow_key, tuple) else flow_key, 1.0
                ),
            )
        )

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
    cfg = pd.read_csv(config_csv, index_col=0)
    if cfg.empty:
        raise ValueError(f"No configs found in {config_csv}.")

    def _normalize_config_name(name: object) -> object:
        """normalize config name.

        :param name: Value for `name`.
        :type name: object
        :returns: Return value.
        :rtype: object"""
        if name in cfg.index:
            return name
        try:
            name_str = str(name)
            if name_str in cfg.index:
                return name_str
        except Exception:
            pass
        try:
            name_int = int(name)
            if name_int in cfg.index:
                return name_int
        except Exception:
            pass
        return name

    if config_names is not None:
        configs = [_normalize_config_name(c) for c in config_names]
    else:
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
    f.run(progress=progress)
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


def run_fair_delta_rf(
    trails: Any,
    *,
    scenario: str,
    emissions_csv: str | Path = DEFAULT_EMISSIONS_CSV,
    mapping_yaml: str | Path = DEFAULT_MAPPING_YAML,
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
    activity_axis = dims.index("activity")
    flow_axis = dims.index("flow")
    year_axis = dims.index("year")
    root_axis = dims.index("root activity")

    inv_raw = inv.data
    if isinstance(inv_raw, sparse.COO):
        coo_full = inv_raw
    else:
        coo_full = sparse.COO.from_numpy(np.asarray(inv_raw, dtype=float))

    # Reduce once over activity and reuse this COO for both species deltas and
    # per-flow/root allocations to avoid duplicate sparse reductions.
    inv_data = coo_full.sum(axis=activity_axis)
    remaining_axes = [i for i in range(coo_full.ndim) if i != activity_axis]
    flow_new_axis = remaining_axes.index(flow_axis)
    year_new_axis = remaining_axes.index(year_axis)
    root_new_axis = remaining_axes.index(root_axis)
    if inv_data.ndim != 3:
        raise ValueError(
            "Unexpected inventory shape after aggregation to flow/year/root."
        )
    if (flow_new_axis, year_new_axis, root_new_axis) != (0, 1, 2):
        inv_data = inv_data.transpose((flow_new_axis, year_new_axis, root_new_axis))

    n_flow = int(inv.sizes["flow"])
    n_root = int(inv.sizes["root activity"])

    flow_coord = inv.coords["flow"].values
    coord_value_to_pos = {int(v): i for i, v in enumerate(flow_coord)}

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

    # Build perturbations from activity-reduced inventory (single sparse reduction path).
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

    species_to_positions: Dict[str, list[int]] = {}
    for pos, flow_key in flow_pos_to_key.items():
        specie = species_map.get(flow_key)
        if specie is None and isinstance(flow_key, tuple):
            specie = species_map.get(flow_key[0])
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
                sign_val = float(
                    signs.get(
                        flow_key,
                        signs.get(
                            flow_key[0] if isinstance(flow_key, tuple) else flow_key,
                            1.0,
                        ),
                    )
                )
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
                sign_val = float(
                    signs.get(
                        flow_key,
                        signs.get(
                            flow_key[0] if isinstance(flow_key, tuple) else flow_key,
                            1.0,
                        ),
                    )
                )
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

        alias = rf_alias.get(specie, specie)
        if alias not in delta_forcing.coords["specie"].values:
            return None

        if config_names is not None and len(config_names) > 1:
            rf_series = _extract_fair_timeseries_by_config(
                delta_forcing.sel(specie=alias)
            )
            rf_quant = _safe_nanpercentile(rf_series, quantiles)
            temp_series = _extract_fair_timeseries_by_config(delta_temp)
            temp_quant = _safe_nanpercentile(temp_series, quantiles)
        else:
            rf_series = np.asarray(delta_forcing.sel(specie=alias).values, dtype=float)
            rf_series = np.nan_to_num(rf_series, nan=0.0)
            rf_series = np.ravel(rf_series)
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

        for specie in delta_forcing.coords["specie"].values:
            if config_names is not None and len(config_names) > 1:
                rf_series = _extract_fair_timeseries_by_config(
                    delta_forcing.sel(specie=specie)
                )
                rf_quant = _safe_nanpercentile(rf_series, quantiles)
                temp_series = _extract_fair_timeseries_by_config(delta_temp)
                temp_quant = _safe_nanpercentile(temp_series, quantiles)
            else:
                rf_series = np.asarray(
                    delta_forcing.sel(specie=specie).values, dtype=float
                )
                rf_series = np.nan_to_num(rf_series, nan=0.0)
                rf_series = np.ravel(rf_series)
                rf_quant = np.tile(rf_series[None, :], (n_quant, 1))
                temp_series = _extract_fair_timeseries(delta_temp)
                temp_series = np.nan_to_num(temp_series, nan=0.0)
                temp_quant = np.tile(temp_series[None, :], (n_quant, 1))
            for qi in range(n_quant):
                _append_allocated_rf(str(specie), rf_quant[qi], "pos", quantile_idx=qi)
                _append_allocated_temp(
                    str(specie), temp_quant[qi], "pos", quantile_idx=qi
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
