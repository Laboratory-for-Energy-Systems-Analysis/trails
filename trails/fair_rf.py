"""FaIR integration helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

import pandas as pd
import numpy as np
import sparse
import xarray as xr

from tqdm import tqdm

import fair  # hard dependency
from fair.io import read_properties

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




def _inventory_emissions_by_fair_species(
    trails: Any,
    species_map: dict[object, str],
    signs: dict[object, float],
) -> pd.DataFrame:
    """Aggregate trails.inventory into FaIR species emissions (kg/yr)."""
    inv = trails.inventory
    if inv is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    if "root activity" in inv.dims:
        inv_sum = inv.sum(dim=["activity", "root activity"])
    else:
        inv_sum = inv.sum(dim=["activity"])

    inv_sum = inv_sum.transpose("flow", "year")
    years = [int(y) for y in inv_sum.coords["year"].values.tolist()]

    # flow position -> name
    flow_coord = inv_sum.coords["flow"].values
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
            fid_int = int(fid)
            if fid_int in coord_value_to_pos:
                pos = coord_value_to_pos[fid_int]
            elif 0 <= fid_int < len(flow_coord):
                pos = fid_int
            else:
                continue
            flow_pos_to_key.setdefault(pos, flow_key)

    data = inv_sum.data
    if hasattr(data, "todense"):
        dense = pd.DataFrame(
            data.todense(), index=range(len(flow_coord)), columns=years
        )
    else:
        dense = pd.DataFrame(data, index=range(len(flow_coord)), columns=years)

    agg: Dict[str, pd.Series] = {}
    for pos, flow_key in flow_pos_to_key.items():
        fair_species = species_map.get(flow_key)
        if fair_species is None and isinstance(flow_key, tuple):
            fair_species = species_map.get(flow_key[0])
        if fair_species is None:
            continue
        sign = float(
            signs.get(
                flow_key,
                signs.get(flow_key[0] if isinstance(flow_key, tuple) else flow_key, 1.0),
            )
        )
        series = dense.loc[pos]
        if sign != 1.0:
            series = series.abs() * sign
        agg[fair_species] = agg.get(fair_species, 0.0) + series

    if not agg:
        return pd.DataFrame(index=years)

    out = pd.DataFrame(agg)
    out.index.name = "year"
    return out


def _run_fair_emissions(
    emissions_df: pd.DataFrame,
    scenario: str,
    *,
    config_csv: str | Path | None = DEFAULT_CONFIGS_CSV,
    properties_csv: str | Path | None = DEFAULT_PROPERTIES_CSV,
    config_name: str | None = None,
    ghg_method: str | None = "myhre1998",
    progress: bool = False,
) -> fair.FAIR:
    df = _normalize_emissions_columns(emissions_df)
    df = df[(df["scenario"] == scenario) & (df["region"].str.lower() == "world")].copy()
    if df.empty:
        raise ValueError(f"Scenario '{scenario}' not found in emissions data.")

    year_cols, year_vals = _extract_year_columns(df)
    if not year_cols:
        raise ValueError("No year columns found in emissions data.")

    meta_cols = ["scenario", "region", "variable", "unit"]
    df = df[list(meta_cols) + year_cols]
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
    species, properties = read_properties(filename=str(properties_csv), species=species)

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
    if config_name is None:
        config_name = cfg.index[0]
    configs = [config_name]

    f = fair.FAIR()
    if ghg_method is not None:
        f.ghg_method = ghg_method
    f.define_time(start_year, end_year, 1)
    f.define_scenarios([scenario])
    f.define_configs(configs)
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(filename=str(properties_csv))
    f.override_defaults(str(config_csv))
    f.fill_from_pandas(mode="emissions", df=df)
    f.run(progress=progress)
    if not np.isfinite(f.forcing.values).any():
        forcing = _compute_ghg_forcing_from_concentration(f)
        if forcing is not None:
            forcing_array = f.forcing.data
            forcing_array[..., f._ghg_indices] = forcing[..., f._ghg_indices]
            f.forcing.data = forcing_array
    return f


def _compute_ghg_forcing_from_concentration(f: fair.FAIR) -> np.ndarray | None:
    """Compute GHG forcing from concentration using FAIR's own forcing functions."""
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


def run_fair_delta_rf(
    trails: Any,
    *,
    scenario: str,
    emissions_csv: str | Path = DEFAULT_EMISSIONS_CSV,
    mapping_yaml: str | Path = DEFAULT_MAPPING_YAML,
    config_csv: str | Path | None = DEFAULT_CONFIGS_CSV,
    properties_csv: str | Path | None = DEFAULT_PROPERTIES_CSV,
    config_name: str | None = None,
    ghg_method: str | None = "myhre1998",
    scale_factor: float | None = None,
    scale_target_fraction: float = 0.01,
    scaling_factor: float | None = None,
    validate_emissions_delta: bool = True,
    validate_atol: float = 1e-6,
    validate_rtol: float = 1e-6,
    validate_raise: bool = False,
    per_species_runs: bool = True,
) -> xr.DataArray:
    """Run FaIR baseline/perturbed and store delta RF on trails."""
    debug = bool(getattr(trails, "debug", False))
    if scaling_factor is not None:
        scale_factor = float(scaling_factor)
    if scale_factor is None:
        if scale_target_fraction <= 0:
            raise ValueError("scale_target_fraction must be > 0.")
    elif scale_factor <= 0:
        raise ValueError("scale_factor must be > 0.")
    df = load_emissions_csv(emissions_csv)
    species_map, signs = load_species_mapping(mapping_yaml)

    # Baseline run
    f_base = _run_fair_emissions(
        df,
        scenario,
        config_csv=config_csv,
        properties_csv=properties_csv,
        config_name=config_name,
        ghg_method=ghg_method,
        progress=False,
    )

    # Build perturbations from Trails inventory
    delta_by_species = _inventory_emissions_by_fair_species(
        trails, species_map, signs
    )
    no_perturbation = delta_by_species.empty
    if no_perturbation:
        if debug:
            print("FAIR debug: no inventory emissions matched; returning zero delta RF.")
        if scale_factor is None:
            scale_factor = 1.0

    if scale_factor is None and not delta_by_species.empty:
        df_base = df[(df["scenario"] == scenario) & (df["region"].str.lower() == "world")].copy()
        df_base = _normalize_emissions_columns(df_base)
        candidates = []
        for specie in delta_by_species.columns:
            rows = df_base[df_base["variable"] == specie]
            if rows.empty:
                continue
            unit = rows["unit"].iloc[0]
            base_row = rows.iloc[0]
            for year, val in delta_by_species[specie].items():
                ycol = str(int(year))
                if ycol not in df_base.columns:
                    continue
                base_val = base_row[ycol]
                if base_val == 0 or pd.isna(base_val):
                    continue
                delta_unit = _convert_kg_to_unit(np.array([val]), unit)[0]
                if delta_unit == 0:
                    continue
                candidates.append(scale_target_fraction * abs(base_val) / abs(delta_unit))
        if candidates:
            scale_factor = float(min(candidates))
        else:
            scale_factor = 1.0
        if debug:
            print("FAIR debug: auto scale_factor", scale_factor)

    if scale_factor != 1.0 and not delta_by_species.empty:
        delta_by_species = delta_by_species * float(scale_factor)

    def _build_perturbed_df(
        base_df: pd.DataFrame,
        specie: str | None,
        delta_series: pd.Series | None = None,
    ) -> pd.DataFrame:
        df_pert_local = _normalize_emissions_columns(base_df)
        df_pert_local = df_pert_local[
            (df_pert_local["scenario"] == scenario)
            & (df_pert_local["region"].str.lower() == "world")
        ].copy()
        if specie is None:
            return df_pert_local
        rows = df_pert_local[df_pert_local["variable"] == specie]
        if rows.empty:
            return df_pert_local
        unit = rows["unit"].iloc[0]
        series = delta_series if delta_series is not None else delta_by_species[specie]
        for year, val in series.items():
            ycol = str(int(year))
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
                expected = delta_by_species[specie].reindex(year_vals, fill_value=0.0)
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
            if issues:
                msg = "FAIR emissions delta check failed: " + "; ".join(issues)
                if validate_raise:
                    raise ValueError(msg)
                if debug:
                    print(msg)

    forcing_base = f_base.forcing.sel(scenario=scenario, config=f_base.configs[0])
    fair_years = [int(y) for y in forcing_base.coords["timebounds"].values.tolist()]
    year_to_fair_idx = {y: i for i, y in enumerate(fair_years)}

    inv = trails.inventory
    if inv is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    if "root activity" not in inv.dims:
        raise ValueError("Trails.inventory must include 'root activity' for delta RF.")

    inv_sum = inv.sum(dim=["activity"]).transpose("flow", "year", "root activity")
    inv_data = inv_sum.data
    if not isinstance(inv_data, sparse.COO):
        inv_data = sparse.COO.from_numpy(np.asarray(inv_data))

    n_flow = int(inv_sum.sizes["flow"])
    n_root = int(inv_sum.sizes["root activity"])

    flow_coord = inv_sum.coords["flow"].values
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

    # Build per-species sparse entries
    flow_idx = inv_data.coords[0]
    year_idx = inv_data.coords[1]
    root_idx = inv_data.coords[2]
    data = inv_data.data.astype(float, copy=False)

    inv_years = [int(y) for y in inv_sum.coords["year"].values.tolist()]

    coords_out = []
    data_out = []

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
        specie: str, rf_series: np.ndarray, sign_mode: str
    ) -> None:
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
            fi = year_to_fair_idx.get(int(y))
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
        coords_out.append(np.vstack([y_out, f_out, r_out]))
        data_out.append(RF_alloc[row_idx, col_idx].astype(float, copy=False))

    def _debug_nan_forcing_check(forcing_pert: xr.DataArray) -> None:
        if not debug:
            return
        if np.isfinite(forcing_base.values).any() and np.isfinite(forcing_pert.values).any():
            return
        print("FAIR debug: forcing all-NaN detected.")
        print("FAIR debug: ghg_method", getattr(f_base, "ghg_method", None))
        try:
            print("FAIR debug: ghg flag", bool(getattr(f_base, "_routine_flags", {}).get("ghg", False)))
        except Exception:
            print("FAIR debug: ghg flag unavailable")
        print("FAIR debug: species count", len(f_base.species))
        print(
            "FAIR debug: co2/ch4/n2o indices",
            int(np.sum(getattr(f_base, "_co2_indices", []))),
            int(np.sum(getattr(f_base, "_ch4_indices", []))),
            int(np.sum(getattr(f_base, "_n2o_indices", []))),
        )
        for s in ["CO2", "CH4", "N2O"]:
            if s in f_base.concentration.coords["specie"].values:
                arr = f_base.concentration.sel(specie=s).values
                print(
                    f"FAIR debug: {s} conc min/max/finite",
                    float(np.nanmin(arr)),
                    float(np.nanmax(arr)),
                    int(np.isfinite(arr).sum()),
                )
            if s in f_base.forcing.coords["specie"].values:
                farr = f_base.forcing.sel(specie=s).values
                print(
                    f"FAIR debug: {s} forcing finite",
                    int(np.isfinite(farr).sum()),
                )
        try:
            scale = f_base.species_configs["forcing_scale"].values
            basec = f_base.species_configs["baseline_concentration"].values
            print(
                "FAIR debug: forcing_scale nan count",
                int(np.isnan(scale).sum()),
            )
            print(
                "FAIR debug: baseline_concentration nan count",
                int(np.isnan(basec).sum()),
            )
            import fair.forcing.ghg as _ghg
            conc = f_base.concentration.values
            base_arr = f_base.species_configs["baseline_concentration"].values
            scale_arr = f_base.species_configs["forcing_scale"].values
            rad_eff = f_base.species_configs[
                "greenhouse_gas_radiative_efficiency"
            ].values
            base_full = base_arr[None, None, ...] * np.ones(
                (1, f_base._n_scenarios, f_base._n_configs, f_base._n_species)
            )
            scale_full = scale_arr[None, None, ...] * np.ones(
                (1, f_base._n_scenarios, f_base._n_configs, f_base._n_species)
            )
            rad_full = rad_eff[None, None, ...] * np.ones(
                (1, f_base._n_scenarios, f_base._n_configs, f_base._n_species)
            )
            test = _ghg.myhre1998(
                conc,
                base_full,
                scale_full,
                rad_full,
                f_base._co2_indices,
                f_base._ch4_indices,
                f_base._n2o_indices,
                f_base._minor_ghg_indices,
            )
            print(
                "FAIR debug: myhre1998 finite",
                int(np.isfinite(test).sum()),
            )
        except Exception as exc:
            print("FAIR debug: species_configs check failed:", exc)

    if no_perturbation:
        if debug:
            print("FAIR debug: no perturbations; delta RF will be zero.")
    elif per_species_runs:
        delta_by_species_pos = delta_by_species.clip(lower=0.0)
        delta_by_species_neg = delta_by_species.clip(upper=0.0)
        pbar = tqdm(
            delta_by_species.columns.tolist(),
            desc="FaIR per-species",
            unit="specie",
        )
        for specie in pbar:
            pbar.set_postfix_str(str(specie))
            series_pos = delta_by_species_pos[specie]
            if np.any(series_pos.values != 0):
                df_pert = _build_perturbed_df(df, specie, series_pos)
                f_pert = _run_fair_emissions(
                    df_pert,
                    scenario,
                    config_csv=config_csv,
                    properties_csv=properties_csv,
                    config_name=config_name,
                    ghg_method=ghg_method,
                    progress=False,
                )
                forcing_pert = f_pert.forcing.sel(
                    scenario=scenario, config=f_pert.configs[0]
                )
                _debug_nan_forcing_check(forcing_pert)
                delta_forcing = forcing_pert - forcing_base
                if scale_factor is None:
                    scale_factor = 1.0
                if scale_factor != 1.0:
                    delta_forcing = delta_forcing / float(scale_factor)

                alias = rf_alias.get(specie, specie)
                if alias in delta_forcing.coords["specie"].values:
                    rf_series = np.asarray(
                        delta_forcing.sel(specie=alias).values, dtype=float
                    )
                    rf_series = np.nan_to_num(rf_series, nan=0.0)
                    _append_allocated_rf(specie, rf_series, "pos")

            series_neg = delta_by_species_neg[specie]
            if np.any(series_neg.values != 0):
                df_pert = _build_perturbed_df(df, specie, series_neg)
                f_pert = _run_fair_emissions(
                    df_pert,
                    scenario,
                    config_csv=config_csv,
                    properties_csv=properties_csv,
                    config_name=config_name,
                    ghg_method=ghg_method,
                    progress=False,
                )
                forcing_pert = f_pert.forcing.sel(
                    scenario=scenario, config=f_pert.configs[0]
                )
                _debug_nan_forcing_check(forcing_pert)
                delta_forcing = forcing_pert - forcing_base
                if scale_factor is None:
                    scale_factor = 1.0
                if scale_factor != 1.0:
                    delta_forcing = delta_forcing / float(scale_factor)

                alias = rf_alias.get(specie, specie)
                if alias in delta_forcing.coords["specie"].values:
                    rf_series = np.asarray(
                        delta_forcing.sel(specie=alias).values, dtype=float
                    )
                    rf_series = np.nan_to_num(rf_series, nan=0.0)
                    _append_allocated_rf(specie, rf_series, "neg")
    else:
        df_pert = _build_perturbed_df(df, None)
        f_pert = _run_fair_emissions(
            df_pert,
            scenario,
            config_csv=config_csv,
            properties_csv=properties_csv,
            config_name=config_name,
            ghg_method=ghg_method,
            progress=False,
        )
        forcing_pert = f_pert.forcing.sel(
            scenario=scenario, config=f_pert.configs[0]
        )
        _debug_nan_forcing_check(forcing_pert)
        delta_forcing = forcing_pert - forcing_base
        if scale_factor is None:
            scale_factor = 1.0
        if scale_factor != 1.0:
            delta_forcing = delta_forcing / float(scale_factor)

        for specie in delta_forcing.coords["specie"].values:
            rf_series = np.asarray(delta_forcing.sel(specie=specie).values, dtype=float)
            rf_series = np.nan_to_num(rf_series, nan=0.0)
            _append_allocated_rf(str(specie), rf_series, "pos")

    if coords_out:
        coords = np.hstack(coords_out)
        data = np.concatenate(data_out)
        rf_sparse = sparse.COO(
            coords=coords,
            data=data,
            shape=(len(fair_years), n_flow, n_root),
        )
    else:
        rf_sparse = sparse.COO(
            coords=np.zeros((3, 0), dtype=int),
            data=np.array([], dtype=float),
            shape=(len(fair_years), n_flow, n_root),
        )

    trails.instant_radiative_forcing = xr.DataArray(
        rf_sparse,
        dims=("year", "flow", "root activity"),
        coords={
            "year": np.array(fair_years, dtype=int),
            "flow": inv_sum.coords["flow"],
            "root activity": inv_sum.coords["root activity"],
        },
    )
    return trails.instant_radiative_forcing
