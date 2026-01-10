from __future__ import annotations

from typing import List, TYPE_CHECKING

import numpy as np
import sparse
import xarray as xr

from .lcia import get_lcia_methods

if TYPE_CHECKING:
    from .trails import Trails


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

    flowkey_to_flowid = _build_flowkey_to_flowid(trails)

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
