from __future__ import annotations

from typing import List, TYPE_CHECKING

import numpy as np
import sparse
import xarray as xr

from .lcia import get_lcia_methods

if TYPE_CHECKING:
    from .trails import Trails


def _build_flowkey_to_flowindex(trails: Trails) -> dict[tuple, int]:
    """build flowkey to flowindex.

    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: dict[tuple, int]
    :raises TypeError: If an error occurs."""
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
    """build flowkey to flowid.

    :param trails: Value for `trails`.
    :type trails: Trails
    :returns: Return value.
    :rtype: dict[tuple, int]
    :raises TypeError: If an error occurs."""
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


def _build_cf_matrix_flowid_space(
    trails: Trails,
    methods: List[str],
    ei_version: str,
    char_cache: dict[tuple, np.ndarray],
    debug: bool = False,
) -> np.ndarray:
    """build cf matrix flowid space.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param methods: Value for `methods`.
    :type methods: List[str]
    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :param char_cache: Value for `char_cache`.
    :type char_cache: dict[tuple, np.ndarray]
    :param debug: Value for `debug`.
    :type debug: bool
    :returns: Return value.
    :rtype: np.ndarray
    :raises ValueError: If an error occurs."""
    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        return np.zeros((0, 0), dtype=np.float64)

    cache_key = ("cf_matrix_flowid_space", ei_version, tuple(methods))
    if cache_key in char_cache:
        return char_cache[cache_key]

    flowkey_to_flowid = _build_flowkey_to_flowindex(trails)

    methods_dict = get_lcia_methods(methods=methods, ei_version=ei_version)

    cf = np.zeros((len(methods), n_flows), dtype=np.float64)

    for m, mname in enumerate(methods):
        exc = methods_dict.get(mname)
        if exc is None:
            raise ValueError(f"LCIA method not found: {mname}")
        for flow_key, val in exc.items():
            fid = flowkey_to_flowid.get(flow_key)
            if fid is None:
                continue
            cf[m, fid] += float(val)

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
    """Build characterized inventory.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param methods: Value for `methods`.
    :type methods: List[str]
    :param char_cache: Value for `char_cache`.
    :type char_cache: dict[tuple, np.ndarray]
    :param debug: Value for `debug`.
    :type debug: bool
    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: xr.DataArray
    :raises ValueError: If an error occurs."""
    if trails.inventory is None:
        raise ValueError("Trails.inventory is empty; run LCA first.")

    cf = _build_cf_matrix_flowid_space(
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
    cf = cf.astype(np.float64, copy=False)
    if len(methods) == 1:
        if has_root:
            cf_b = cf[0][None, :, None, None]
            characterized = inv_data * cf_b
        else:
            cf_b = cf[0][None, :, None]
            characterized = inv_data * cf_b
        characterized = sparse.stack([characterized], axis=0)
    else:
        blocks = []
        if has_root:
            for m in range(cf.shape[0]):
                blocks.append(inv_data * cf[m][None, :, None, None])
        else:
            for m in range(cf.shape[0]):
                blocks.append(inv_data * cf[m][None, :, None])
        characterized = sparse.stack(blocks, axis=0)

    dims = ("method", "activity", "flow", "year")
    coords = {
        "method": np.asarray(methods, dtype=object),
        "activity": inventory.coords["activity"],
        "flow": inventory.coords["flow"],
        "year": inventory.coords["year"],
    }
    if has_root:
        dims = ("method", "activity", "flow", "year", "root activity")
        coords["root activity"] = inventory.coords["root activity"]

    trails.characterized_inventory = xr.DataArray(
        characterized,
        dims=dims,
        coords=coords,
    )
    return trails.characterized_inventory


def get_cf_vector(
    trails: Trails,
    methods: List[str],
    char_cache: dict[tuple, np.ndarray],
    *,
    debug: bool = False,
    ei_version: str = "3.11",
) -> np.ndarray:
    """Get cf vector.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param methods: Value for `methods`.
    :type methods: List[str]
    :param char_cache: Value for `char_cache`.
    :type char_cache: dict[tuple, np.ndarray]
    :param debug: Value for `debug`.
    :type debug: bool
    :param ei_version: Value for `ei_version`.
    :type ei_version: str
    :returns: Return value.
    :rtype: np.ndarray"""
    cf = _build_cf_matrix_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )
    return cf.sum(axis=0)
