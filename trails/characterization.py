from __future__ import annotations

import hashlib
import json
from typing import List, TYPE_CHECKING

import numpy as np
import sparse
import xarray as xr

from .lcia import get_lcia_methods
from .chunked_inventory import is_chunked_sparse, iter_sparse_blocks

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

    flowkey_to_flowid = _build_flowkey_to_flowindex(trails)
    hasher = hashlib.blake2b(digest_size=20)
    hasher.update(str(n_flows).encode("ascii"))
    for flow_key, flow_position in sorted(
        flowkey_to_flowid.items(),
        key=lambda item: (
            str(item[0][0]),
            str(item[0][1]),
            str(item[0][2]),
            int(item[1]),
        ),
    ):
        record = [*(str(part) for part in flow_key), int(flow_position)]
        hasher.update(
            json.dumps(record, ensure_ascii=False, separators=(",", ":")).encode(
                "utf-8"
            )
        )
        hasher.update(b"\x00")
    flow_fingerprint = hasher.hexdigest()

    cache_key = (
        "cf_matrix_flowid_space",
        ei_version,
        tuple(methods),
        flow_fingerprint,
    )
    if cache_key in char_cache:
        return char_cache[cache_key]

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
        raise ValueError("Trails.inventory is empty; run lci() first.")

    cf = _build_cf_matrix_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )

    inventory = trails.inventory

    if is_chunked_sparse(inventory.data):
        cf_da = xr.DataArray(
            cf,
            dims=("method", "flow"),
            coords={
                "method": np.asarray(methods, dtype=object),
                "flow": inventory.coords["flow"],
            },
        )
        characterized = (
            inventory.expand_dims(method=np.asarray(methods, dtype=object)) * cf_da
        )
        dims = ["method", "activity", "flow", "year"]
        if "root activity" in inventory.dims:
            dims.append("root activity")
        trails.characterized_inventory = characterized.transpose(*dims)
        return trails.characterized_inventory

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


def score_inventory_with_regular_methods(
    trails: Trails,
    methods: List[str],
    char_cache: dict[tuple, np.ndarray],
    *,
    ei_version: str = "3.11",
    show_progress: bool = True,
) -> xr.DataArray:
    """Incrementally reduce a finalized inventory with regular CF vectors."""
    if trails.inventory is None:
        raise ValueError("Regular LCIA requires a finalized Trails inventory.")
    if not methods:
        raise ValueError("methods must contain at least one regular LCIA method.")

    cf = _build_cf_matrix_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
    )
    inventory = trails.inventory.transpose(
        "activity",
        "flow",
        "year",
        *[
            dim
            for dim in trails.inventory.dims
            if dim not in {"activity", "flow", "year"}
        ],
    )
    has_root = "root activity" in inventory.dims
    n_methods = len(methods)
    n_activities = int(inventory.sizes["activity"])
    n_years = int(inventory.sizes["year"])
    n_roots = int(inventory.sizes["root activity"]) if has_root else 0
    selected_flows = np.flatnonzero(np.any(cf != 0.0, axis=0)).astype(np.int64)

    output_coords: list[np.ndarray] = []
    output_data: list[np.ndarray] = []

    def append_entries(
        activities: np.ndarray,
        flows: np.ndarray,
        years: np.ndarray,
        values: np.ndarray,
        roots: np.ndarray | None,
    ) -> None:
        for method_idx in range(n_methods):
            scored = values.astype(np.float64, copy=False) * cf[method_idx, flows]
            keep = scored != 0.0
            count = int(np.count_nonzero(keep))
            if not count:
                continue
            axes: list[np.ndarray] = []
            if n_methods > 1:
                axes.append(np.full(count, method_idx, dtype=np.int64))
            axes.extend([activities[keep], years[keep]])
            if roots is not None:
                axes.append(roots[keep])
            output_coords.append(np.vstack(axes))
            output_data.append(scored[keep])

    builder = getattr(trails, "_inventory_builder", None)
    if (
        has_root
        and builder is not None
        and getattr(builder, "_finalized", False)
        and hasattr(builder, "iter_entries_for_flows")
    ):
        for (
            _year_block,
            activities,
            flows,
            years,
            roots,
            values,
        ) in builder.iter_entries_for_flows(
            selected_flows,
            show_progress=show_progress,
            progress_desc="Regular LCIA inventory",
        ):
            append_entries(activities, flows, years, values, roots)
    elif is_chunked_sparse(inventory.data):
        selected_mask = np.zeros(int(inventory.sizes["flow"]), dtype=bool)
        selected_mask[selected_flows] = True
        for slices, block in iter_sparse_blocks(inventory.data, primary_axis=2):
            coords = block.coords.astype(np.int64, copy=False)
            activities = coords[0] + int(slices[0].start or 0)
            flows = coords[1] + int(slices[1].start or 0)
            years = coords[2] + int(slices[2].start or 0)
            roots = coords[3] + int(slices[3].start or 0) if has_root else None
            keep = selected_mask[flows]
            if not np.any(keep):
                continue
            append_entries(
                activities[keep],
                flows[keep],
                years[keep],
                block.data[keep],
                roots[keep] if roots is not None else None,
            )
    else:
        data = inventory.data
        if not isinstance(data, sparse.COO):
            data = sparse.COO.from_numpy(np.asarray(data))
        coords = data.coords.astype(np.int64, copy=False)
        append_entries(
            coords[0],
            coords[1],
            coords[2],
            data.data,
            coords[3] if has_root else None,
        )

    if has_root and n_methods > 1:
        shape = (n_methods, n_activities, n_years, n_roots)
        dims = ("method", "activity", "year", "root activity")
    elif has_root:
        shape = (n_activities, n_years, n_roots)
        dims = ("activity", "year", "root activity")
    elif n_methods > 1:
        shape = (n_methods, n_activities, n_years)
        dims = ("method", "activity", "year")
    else:
        shape = (n_activities, n_years)
        dims = ("activity", "year")

    coords_xr: dict[str, np.ndarray] = {
        "activity": inventory.coords["activity"].values,
        "year": inventory.coords["year"].values,
    }
    if n_methods > 1:
        coords_xr["method"] = np.asarray(methods, dtype=object)
    if has_root:
        coords_xr["root activity"] = inventory.coords["root activity"].values

    if output_coords:
        scores_data = sparse.COO(
            np.concatenate(output_coords, axis=1),
            np.concatenate(output_data).astype(np.float64, copy=False),
            shape=shape,
        )
    else:
        scores_data = sparse.zeros(shape, dtype=np.float64)
    trails.scores = xr.DataArray(scores_data, dims=dims, coords=coords_xr)
    return trails.scores


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


def get_cf_matrix(
    trails: Trails,
    methods: List[str],
    char_cache: dict[tuple, np.ndarray],
    *,
    debug: bool = False,
    ei_version: str = "3.11",
) -> np.ndarray:
    """Get characterization factor matrix aligned to Trails biosphere flows.

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
    :returns: Characterization factors with shape ``(methods, flows)``.
    :rtype: np.ndarray"""
    return _build_cf_matrix_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )
