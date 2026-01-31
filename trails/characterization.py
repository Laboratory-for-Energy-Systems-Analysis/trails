from __future__ import annotations

from pathlib import Path
from typing import Any, List, TYPE_CHECKING

import yaml

import numpy as np
import sparse
import xarray as xr

from .dynamic_gwp import CO2Params, WMGHGParams, rf_co2_from_annual_emissions, rf_wmghg_from_annual_emissions
from .lcia import get_lcia_methods

if TYPE_CHECKING:
    from .trails import Trails


def _build_flowkey_to_flowindex(trails: Trails) -> dict[tuple, int]:
    """
    Build (name, compartment, subcompartment) -> flow_index mapping,
    where flow_index is the position along the Trails/B/inventory flow dimension.
    """
    # Prefer explicit flow coordinate if available
    flow_coord = None
    if (
        getattr(trails, "inventory", None) is not None
        and "flow" in trails.inventory.coords
    ):
        flow_coord = trails.inventory.coords["flow"].values
    elif getattr(trails, "B", None) is not None:
        # If Trails has a flow coordinate elsewhere, use it; otherwise fall back to range
        flow_coord = np.arange(int(trails.B.shape[2]), dtype=int)

    # Build mapping from whatever the coordinate values are -> position
    coord_value_to_pos = {int(v): i for i, v in enumerate(flow_coord)}

    out: dict[tuple, int] = {}
    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        if not meta:
            continue

        k0 = next(iter(meta.keys()))
        if not isinstance(k0, (int, np.integer)):
            raise TypeError(
                f"Unexpected biosphere_indices structure for label {_label}: "
                f"expected int keys, got {type(k0)}"
            )

        for key_int, md in meta.items():
            if not isinstance(md, dict):
                continue
            name = md.get("name")
            comp = md.get("compartment")
            sub = md.get("subcompartment")
            if name is None or comp is None or sub is None:
                continue

            # key_int might be a coordinate value (often it is), not necessarily a positional index
            # Convert coordinate value -> positional index where possible; otherwise assume it is already positional
            k_int = int(key_int)
            flow_index = coord_value_to_pos.get(k_int, None)

            if flow_index is None:
                # Fallback: if key_int looks like a position and is in-range, accept it
                if 0 <= k_int < len(flow_coord):
                    flow_index = k_int
                else:
                    continue  # cannot align this entry safely

            out.setdefault((name, comp, sub), int(flow_index))

    return out


def _build_flowkey_to_flowid(trails: Trails) -> dict[tuple, int]:
    """Build a flow-key to flow-id mapping across labels.

    :param trails: Trails instance with biosphere metadata.
    :type trails: Trails
    :returns: Mapping of ``(name, compartment, subcompartment)`` to flow id.
    :rtype: dict[tuple, int]
    """
    out: dict[tuple, int] = {}

    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        if not meta:
            continue

        k0 = next(iter(meta.keys()))

        # Expected: {flow_id: {"name":..., "compartment":..., "subcompartment":...}}
        if isinstance(k0, (int, np.integer)):
            for flow_id, md in meta.items():
                if not isinstance(md, dict):
                    continue
                name = md.get("name")
                comp = md.get("compartment")
                sub = md.get("subcompartment")
                if name is None or comp is None or sub is None:
                    continue
                out.setdefault((name, comp, sub), int(flow_id))
        else:
            raise TypeError(
                f"Unexpected biosphere_indices structure for label {_label}: "
                f"expected int keys, got {type(k0)}"
            )

    return out


def _build_cf_vector_flowid_space(
    trails: Trails,
    methods: List[str],
    ei_version: str,
    char_cache: dict[tuple, np.ndarray],
    debug: bool = False,
) -> np.ndarray:
    """Build a CF vector aligned with Trails flow-id space.

    :param trails: Trails instance with biosphere metadata.
    :type trails: Trails
    :param methods: LCIA methods to include.
    :type methods: list[str]
    :param ei_version: Ecoinvent release identifier.
    :type ei_version: str
    :param char_cache: Cache mapping for characterization vectors.
    :type char_cache: dict
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: CF vector in flow-id space.
    :rtype: numpy.ndarray
    """
    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        return np.zeros(0, dtype=np.float64)

    cache_key = ("cf_vector_flowid_space", ei_version, tuple(methods))
    if cache_key in char_cache:
        return char_cache[cache_key]

    flowkey_to_flowid = _build_flowkey_to_flowindex(trails)

    methods_dict = get_lcia_methods(methods=methods, ei_version=ei_version)

    cf = np.zeros(n_flows, dtype=np.float64)

    for _mname, exc in methods_dict.items():
        for flow_key, val in exc.items():
            fid = flowkey_to_flowid.get(flow_key)
            if fid is None:
                continue
            cf[fid] += float(val)

    char_cache[cache_key] = cf
    return cf


def build_characterized_inventory(
    trails: Trails,
    methods: List[str],
    char_cache: dict[tuple, np.ndarray],
    *,
    debug: bool = False,
    ei_version: str = "3.11",
) -> xr.DataArray:
    """Characterize a Trails inventory into a sparse array."""
    if trails.inventory is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    cf = _build_cf_vector_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )

    inventory = trails.inventory

    # Enforce canonical dim order for safe broadcasting
    if "root activity" in inventory.dims:
        inventory = inventory.transpose("activity", "flow", "year", "root activity")
    else:
        inventory = inventory.transpose("activity", "flow", "year")

    inv_data = inventory.data
    if not isinstance(inv_data, sparse.COO):
        inv_data = sparse.COO.from_numpy(np.asarray(inv_data))

    has_root = "root activity" in inventory.dims
    if has_root:
        characterized = inv_data * cf.astype(np.float64)[None, :, None, None]
    else:
        characterized = inv_data * cf.astype(np.float64)[None, :, None]

    dims = ("activity", "flow", "year")
    coords = {
        "activity": inventory.coords["activity"],
        "flow": inventory.coords["flow"],
        "year": inventory.coords["year"],
    }
    if has_root:
        dims = ("activity", "flow", "year", "root activity")
        coords["root activity"] = inventory.coords["root activity"]

    trails.characterized_inventory = xr.DataArray(
        characterized,
        dims=dims,
        coords=coords,
    )
    return trails.characterized_inventory


def _load_ipcc_ar6_biosphere_params(
    params_path: str | Path | None,
) -> dict[str, dict[str, Any]]:
    """
    Returns: biosphere flow name -> spec dict
      spec["model"] in {"wmgHg", "co2"}
      spec["params"] is WMGHGParams or CO2Params
      spec["molar_mass_g_per_mol"] present for wmgHg
    """
    if params_path is None:
        params_path = Path(__file__).resolve().parent / "data" / "ipcc_ar6_ghg_params.yaml"
    else:
        params_path = Path(params_path)

    data = yaml.safe_load(params_path.read_text()) or {}
    mapping: dict[str, dict[str, Any]] = {}

    for gas_name, p in data.items():
        biosphere_names = p.get("biosphere_name") or []
        if not biosphere_names:
            continue

        model = (p.get("model") or "wmgHg").strip()
        if p.get("concentration_from_emissions"):
            model = "co2"

        # Prefer explicit molar mass from YAML (recommended)
        molar_mass = p.get("molar_mass_g_per_mol")
        if molar_mass is not None:
            molar_mass = float(molar_mass)

        if model == "co2":
            # Build CO2Params from YAML
            cfe = p.get("concentration_from_emissions") or {}
            forcing = p.get("forcing") or {}

            a0 = float(cfe["a0"])
            a = tuple(float(x) for x in cfe["a"])
            tau = tuple(float(x) for x in cfe["tau"])

            params = CO2Params(
                a0=a0, a=a, tau=tau,
                forcing_formula=str(forcing.get("formula", "myhre1998")),
                C0_ppm=float(forcing.get("C0_ppm", 278.3)),
            )

            spec = {"model": "co2", "params": params, "molar_mass_g_per_mol": molar_mass}

        elif model == "wmgHg":
            lifetime = p.get("lifetime_yr")
            rad_eff = p.get("rad_eff_Wm2_per_ppb")

            if lifetime is None or rad_eff is None:
                # Don't guess tau=inf. Just skip entries that aren't properly specified.
                continue

            params = WMGHGParams(
                lifetime_yr=float(lifetime),
                rad_eff_Wm2_per_ppb=float(rad_eff),
            )

            if molar_mass is None:
                continue

            spec = {"model": "wmgHg", "params": params, "molar_mass_g_per_mol": molar_mass}

        else:
            # Unknown model tag
            continue

        biosphere_signs = p.get("biosphere_signs") or {}

        for name in biosphere_names:
            if name:
                spec_with_sign = dict(spec)
                sign = biosphere_signs.get(name, 1)
                try:
                    spec_with_sign["sign"] = float(sign)
                except (TypeError, ValueError):
                    spec_with_sign["sign"] = 1.0
                mapping[str(name)] = spec_with_sign

    return mapping

def build_instant_radiative_forcing(
    trails: Trails,
    *,
    params_path: str | Path | None = None,
    pulse_placement: str = "midyear",
) -> xr.DataArray:
    if trails.inventory is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    debug = bool(getattr(trails, "debug", False))
    name_to_spec = _load_ipcc_ar6_biosphere_params(params_path)

    inventory = trails.inventory
    has_root = "root activity" in inventory.dims
    if has_root:
        inv_ordered = inventory.transpose("activity", "flow", "year", "root activity")
    else:
        inv_ordered = inventory.transpose("activity", "flow", "year")

    years = np.asarray(inv_ordered.coords["year"].values, dtype=int)
    flow_coord = inv_ordered.coords["flow"].values
    coord_value_to_pos = {int(v): i for i, v in enumerate(flow_coord)}

    # flow position -> biosphere name (align flow ids with coordinate values)
    flow_pos_to_name: dict[int, str] = {}
    total_meta = 0
    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        for flow_id, md in meta.items():
            total_meta += 1
            name = md.get("name") if isinstance(md, dict) else None
            if not name:
                continue
            fid = int(flow_id)
            if fid in coord_value_to_pos:
                pos = coord_value_to_pos[fid]
            elif 0 <= fid < len(flow_coord):
                pos = fid
            else:
                continue
            if pos not in flow_pos_to_name:
                flow_pos_to_name[pos] = str(name)

    inv_data = inv_ordered.data
    n_flow = int(inv_ordered.sizes["flow"])
    n_year = int(inv_ordered.sizes["year"])

    if isinstance(inv_data, sparse.COO):
        coords = inv_data.coords
        flow_idx = coords[1]
        year_idx = coords[2]
        if has_root:
            root_idx = coords[3]
            summed = sparse.COO(
                coords=[flow_idx, year_idx, root_idx],
                data=inv_data.data,
                shape=(n_flow, n_year, int(inv_ordered.sizes["root activity"])),
            )
        else:
            summed = sparse.COO(
                coords=[flow_idx, year_idx],
                data=inv_data.data,
                shape=(n_flow, n_year),
            )
    else:
        if has_root:
            summed = (
                inventory.sum(dim=["activity"])
                .transpose("flow", "year", "root activity")
                .data
            )
        else:
            summed = inventory.sum(dim=["activity"]).transpose("flow", "year").data

    if not has_root:
        inv_dense = (
            np.asarray(summed.todense()) if isinstance(summed, sparse.COO) else np.asarray(summed)
        )
        rf = np.zeros_like(inv_dense, dtype=float)
        matched = 0
        emitted = 0
        for pos, name in flow_pos_to_name.items():
            spec = name_to_spec.get(name)
            if spec is None:
                continue
            matched += 1

            sign = float(spec.get("sign", 1.0))
            emissions = inv_dense[pos, :] * sign
            if np.all(emissions == 0):
                continue
            emitted += 1
            if spec["model"] == "wmgHg":
                p = spec["params"]
                mm = spec["molar_mass_g_per_mol"]
                rf[pos, :] = rf_wmghg_from_annual_emissions(
                    years=years,
                    emissions_kg_per_yr=emissions,
                    params=p,
                    molar_mass_g_per_mol=mm,
                    pulse_placement=pulse_placement,
                )
            else:
                p = spec["params"]
                rf[pos, :] = rf_co2_from_annual_emissions(
                    years=years,
                    emissions_kg_per_yr=emissions,
                    params=p,
                    pulse_placement=pulse_placement,
                )

        trails.instant_radiative_forcing = xr.DataArray(
            rf,
            dims=("flow", "year"),
            coords={"flow": inv_ordered.coords["flow"], "year": inv_ordered.coords["year"]},
        )
        return trails.instant_radiative_forcing

    # Root activity preserved: build sparse RF (flow, year, root activity)
    root_coord = inv_ordered.coords["root activity"]
    n_root = int(inv_ordered.sizes["root activity"])

    rf_coords: list[np.ndarray] = []
    rf_data: list[np.ndarray] = []

    if isinstance(summed, sparse.COO):
        summed_flow = summed
    else:
        summed_flow = sparse.COO.from_numpy(np.asarray(summed))

    matched = 0
    emitted = 0
    for pos, name in flow_pos_to_name.items():
        spec = name_to_spec.get(name)
        if spec is None:
            continue
        matched += 1

        flow_slice = summed_flow[pos, :, :]  # (year, root)
        if flow_slice.nnz == 0:
            continue
        emitted += 1
        year_idx = flow_slice.coords[0]
        root_idx = flow_slice.coords[1]
        sign = float(spec.get("sign", 1.0))
        data = flow_slice.data * sign

        # Build emissions matrix for this flow (year x root_nonzero)
        root_ids = np.unique(root_idx)
        root_pos = {int(r): i for i, r in enumerate(root_ids)}
        E = np.zeros((n_year, len(root_ids)), dtype=float)
        np.add.at(
            E,
            (year_idx, np.array([root_pos[int(r)] for r in root_idx], dtype=int)),
            data,
        )

        if spec["model"] == "wmgHg":
            p = spec["params"]
            mm = spec["molar_mass_g_per_mol"]
            RF = rf_wmghg_from_annual_emissions(
                years=years,
                emissions_kg_per_yr=E,
                params=p,
                molar_mass_g_per_mol=mm,
                pulse_placement=pulse_placement,
            )
        else:
            p = spec["params"]
            RF = rf_co2_from_annual_emissions(
                years=years,
                emissions_kg_per_yr=E,
                params=p,
                pulse_placement=pulse_placement,
            )

        # Store dense RF for these roots into sparse coord lists
        yy = np.repeat(np.arange(n_year, dtype=int), len(root_ids))
        rr = np.tile(root_ids.astype(int), n_year)
        ff = np.full(yy.shape, int(pos), dtype=int)
        rf_coords.append(np.vstack([ff, yy, rr]))
        rf_data.append(RF.reshape(-1))

    if rf_coords:
        coords = np.concatenate(rf_coords, axis=1)
        data = np.concatenate(rf_data)
        rf_sparse = sparse.COO(
            coords=coords,
            data=data,
            shape=(n_flow, n_year, n_root),
        )
    else:
        rf_sparse = sparse.COO(
            coords=np.zeros((3, 0), dtype=int),
            data=np.array([], dtype=float),
            shape=(n_flow, n_year, n_root),
        )


    trails.instant_radiative_forcing = xr.DataArray(
        rf_sparse,
        dims=("flow", "year", "root activity"),
        coords={
            "flow": inv_ordered.coords["flow"],
            "year": inv_ordered.coords["year"],
            "root activity": root_coord,
        },
    )


def get_cf_vector(
    trails: Trails,
    methods: List[str],
    char_cache: dict[tuple, np.ndarray],
    *,
    debug: bool = False,
    ei_version: str = "3.11",
) -> np.ndarray:
    """Return a dense CF vector aligned to trails.B flow dimension."""
    return _build_cf_vector_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )
