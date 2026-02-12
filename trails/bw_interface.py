from __future__ import annotations

from typing import Any, Dict, Tuple, TYPE_CHECKING

import bw_processing as bwp
import numpy as np

if TYPE_CHECKING:
    from .trails import Trails


def _ij_from_coords(X_t: Any) -> tuple[np.ndarray, np.ndarray]:
    """ij from coords.

    :param X_t: Value for `X_t`.
    :type X_t: Any
    :returns: Return value.
    :rtype: tuple[np.ndarray, np.ndarray]
    :raises ValueError: If an error occurs."""
    coords = X_t.coords
    if coords.shape[0] == 3:
        i_idx = np.asarray(coords[1], dtype=np.int64)
        j_idx = np.asarray(coords[2], dtype=np.int64)
    elif coords.shape[0] == 2:
        i_idx = np.asarray(coords[0], dtype=np.int64)
        j_idx = np.asarray(coords[1], dtype=np.int64)
    else:
        raise ValueError(f"Unsupported coords ndim={coords.shape[0]} for sparse array")
    return i_idx, j_idx


def _resolve_matrix_label(trails: Trails, year: int) -> tuple[str, int]:
    """resolve matrix label.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :returns: Return value.
    :rtype: tuple[str, int]"""
    label_for_matrix = str(year)
    if label_for_matrix not in trails.scenario_index:
        years = np.array([int(lbl) for lbl in trails.scenario_labels])
        idx = int(np.argmin(np.abs(years - year)))
        label_for_matrix = trails.scenario_labels[idx]

    t = trails.scenario_index[label_for_matrix]
    return label_for_matrix, t


def _select_metadata_label(trails: Trails, label_for_matrix: str) -> str:
    """select metadata label.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param label_for_matrix: Value for `label_for_matrix`.
    :type label_for_matrix: str
    :returns: Return value.
    :rtype: str"""
    if label_for_matrix in trails.activity_indices:
        return label_for_matrix
    return _nearest_metadata_label_for_year(trails, int(label_for_matrix))


def _make_bw_indices_rowcol(row_idx: np.ndarray, col_idx: np.ndarray) -> np.ndarray:
    """make bw indices rowcol.

    :param row_idx: Value for `row_idx`.
    :type row_idx: np.ndarray
    :param col_idx: Value for `col_idx`.
    :type col_idx: np.ndarray
    :returns: Return value.
    :rtype: np.ndarray"""
    idx = np.empty(len(row_idx), dtype=bwp.INDICES_DTYPE)
    idx["row"] = row_idx.astype(np.uint32, copy=False)
    idx["col"] = col_idx.astype(np.uint32, copy=False)
    return idx


def _warn_on_missing_metadata(
    label_for_matrix: str,
    meta_label: str,
    A_act_idx: np.ndarray,
    A_prod_idx: np.ndarray,
    B_flow_idx: np.ndarray,
    act_meta: dict[int, dict],
    bio_meta: dict[int, dict],
) -> None:
    """warn on missing metadata.

    :param label_for_matrix: Value for `label_for_matrix`.
    :type label_for_matrix: str
    :param meta_label: Value for `meta_label`.
    :type meta_label: str
    :param A_act_idx: Value for `A_act_idx`.
    :type A_act_idx: np.ndarray
    :param A_prod_idx: Value for `A_prod_idx`.
    :type A_prod_idx: np.ndarray
    :param B_flow_idx: Value for `B_flow_idx`.
    :type B_flow_idx: np.ndarray
    :param act_meta: Value for `act_meta`.
    :type act_meta: dict[int, dict]
    :param bio_meta: Value for `bio_meta`.
    :type bio_meta: dict[int, dict]"""
    meta_act_indices = set(act_meta.keys())
    meta_bio_indices = set(bio_meta.keys())

    A_act_indices = set(A_act_idx.tolist())
    A_prod_indices = set(A_prod_idx.tolist())
    B_flow_indices = set(B_flow_idx.tolist())

    A_all_activities = A_act_indices | A_prod_indices

    missing_activities = A_all_activities - meta_act_indices
    if missing_activities:
        print(
            f"⚠️ WARNING: A[{label_for_matrix}] uses {len(missing_activities)} "
            f"activity indices not present in activity_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_activities))[:10]}"
            f"{' ...' if len(missing_activities) > 10 else ''}"
        )

    missing_flows = B_flow_indices - meta_bio_indices
    if missing_flows:
        print(
            f"⚠️ WARNING: B[{label_for_matrix}] uses {len(missing_flows)} "
            f"biosphere flow indices not present in biosphere_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_flows))[:10]}"
            f"{' ...' if len(missing_flows) > 10 else ''}"
        )


def _build_metadata_indices(
    act_meta: dict[int, dict], bio_meta: dict[int, dict]
) -> tuple[dict[tuple, int], dict[tuple, int]]:
    """build metadata indices.

    :param act_meta: Value for `act_meta`.
    :type act_meta: dict[int, dict]
    :param bio_meta: Value for `bio_meta`.
    :type bio_meta: dict[int, dict]
    :returns: Return value.
    :rtype: tuple[dict[tuple, int], dict[tuple, int]]"""
    technosphere_indices: Dict[tuple, int] = {}
    for idx, meta in act_meta.items():
        key = (
            meta["name"],
            meta["reference product"],
            meta["unit"],
            meta["location"],
        )
        technosphere_indices[key] = idx

    biosphere_indices: Dict[tuple, int] = {}
    for idx, meta in bio_meta.items():
        key = (
            meta["name"],
            meta["compartment"],
            meta["subcompartment"],
            meta["unit"],
        )
        biosphere_indices[key] = idx

    return technosphere_indices, biosphere_indices


def _nearest_metadata_label_for_year(trails: Trails, year: int) -> str:
    """nearest metadata label for year.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :returns: Return value.
    :rtype: str
    :raises ValueError: If an error occurs."""
    if not trails.activity_indices:
        raise ValueError("Trails.activity_indices is empty – no metadata available.")

    labels = list(trails.activity_indices.keys())
    years = np.array([int(lbl) for lbl in labels])
    idx = int(np.argmin(np.abs(years - year)))
    return labels[idx]


def build_datapackage_for_year_from_trails(
    trails: Trails,
    year: int,
    zero_biosphere: bool = False,
    debug: bool = False,
) -> tuple[Any, dict[tuple, int], dict[tuple, int], list[tuple[int, int]]]:
    """Build datapackage for year from trails.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :param zero_biosphere: Value for `zero_biosphere`.
    :type zero_biosphere: bool
    :param debug: Value for `debug`.
    :type debug: bool
    :returns: Return value.
    :rtype: tuple[Any, dict[tuple, int], dict[tuple, int], list[tuple[int, int]]]"""

    label_for_matrix, t = _resolve_matrix_label(trails, int(year))

    A_t = trails.A[t, :, :]  # (activity, product)
    B_t = trails.B[t, :, :]  # (activity, flow)

    A_signed = A_t.data.astype(np.float64, copy=False)
    B_signed = B_t.data.astype(np.float64, copy=False)

    meta_label = _select_metadata_label(trails, label_for_matrix)
    act_meta = trails.activity_indices.get(meta_label, {})
    bio_meta = trails.biosphere_indices.get(meta_label, {})

    # A -> bw_processing: non-negative data + flip array
    flip_A = A_signed < 0
    data_A = np.abs(A_signed)

    # B -> bw_processing: non-negative data + flip array
    if zero_biosphere:
        data_B = np.zeros_like(B_signed, dtype=np.float64)
        flip_B = np.zeros_like(B_signed, dtype=bool)
    else:
        flip_B = B_signed < 0
        data_B = np.abs(B_signed)

    # Indices
    A_act_idx, A_prod_idx = _ij_from_coords(A_t)
    B_act_idx, B_flow_idx = _ij_from_coords(B_t)

    # Brightway convention:
    # technosphere_matrix is (product, activity)
    indices_A = _make_bw_indices_rowcol(A_prod_idx, A_act_idx)
    # biosphere_matrix is (flow, activity)
    indices_B = _make_bw_indices_rowcol(B_flow_idx, B_act_idx)

    _warn_on_missing_metadata(
        label_for_matrix=label_for_matrix,
        meta_label=meta_label,
        A_act_idx=A_act_idx,
        A_prod_idx=A_prod_idx,
        B_flow_idx=B_flow_idx,
        act_meta=act_meta,
        bio_meta=bio_meta,
    )

    dp = bwp.create_datapackage()

    # Safety checks
    assert indices_A.shape == data_A.shape, (indices_A.shape, data_A.shape)
    assert indices_A.dtype == bwp.INDICES_DTYPE, indices_A.dtype
    assert np.all(data_A >= 0), "A data must be non-negative for flip convention"
    assert flip_A.shape == data_A.shape

    assert indices_B.shape == data_B.shape, (indices_B.shape, data_B.shape)
    assert indices_B.dtype == bwp.INDICES_DTYPE, indices_B.dtype
    assert np.all(data_B >= 0), "B data must be non-negative for flip convention"
    assert flip_B.shape == data_B.shape

    dp.add_persistent_vector(
        matrix="technosphere_matrix",
        indices_array=indices_A,
        data_array=data_A,
        flip_array=flip_A,
    )
    dp.add_persistent_vector(
        matrix="biosphere_matrix",
        indices_array=indices_B,
        data_array=data_B,
        flip_array=flip_B,
    )

    technosphere_indices, biosphere_indices = _build_metadata_indices(
        act_meta, bio_meta
    )
    uncertain_parameters: list[tuple[int, int]] = []

    return dp, technosphere_indices, biosphere_indices, uncertain_parameters


def _reference_product_from_activity_id(
    lca_obj: Any, activity_id: int
) -> tuple[int, float]:
    """reference product from activity id.

    :param lca_obj: Value for `lca_obj`.
    :type lca_obj: Any
    :param activity_id: Value for `activity_id`.
    :type activity_id: int
    :returns: Return value.
    :rtype: tuple[int, float]
    :raises KeyError: If an error occurs.
    :raises ValueError: If an error occurs."""
    act_map = getattr(lca_obj.dicts, "activity", None)
    prod_map = getattr(lca_obj.dicts, "product", None)
    if not act_map or not prod_map:
        raise ValueError("LCA object missing dicts.activity or dicts.product")

    activity_id = int(activity_id)

    ref_cache = getattr(lca_obj, "_ref_product_cache", None)
    if ref_cache is None:
        ref_cache = {}
        setattr(lca_obj, "_ref_product_cache", ref_cache)
    cached = ref_cache.get(activity_id)
    if cached is not None:
        return cached

    if activity_id not in act_map:
        raise KeyError(f"activity_id={activity_id} not in lca_obj.dicts.activity")

    act_pos = int(act_map[activity_id])
    pos_to_prod_id = getattr(lca_obj, "_pos_to_prod_id_cache", None)
    if pos_to_prod_id is None:
        pos_to_prod_id = {int(pos): int(pid) for pid, pos in prod_map.items()}
        setattr(lca_obj, "_pos_to_prod_id_cache", pos_to_prod_id)

    col = lca_obj.technosphere_matrix[:, act_pos]
    if hasattr(col, "tocoo"):
        coo = col.tocoo()
        rows = coo.row
        vals = coo.data
    else:
        vals = np.asarray(col).ravel()
        rows = np.where(vals != 0)[0]
        vals = vals[rows]

    if rows.size == 0:
        raise ValueError(f"No technosphere entries found for activity_id={activity_id}")

    if activity_id in prod_map:
        prod_row_pos = int(prod_map[activity_id])
        mask = rows == prod_row_pos
        if np.any(mask):
            prod_value = float(vals[mask][0])
        else:
            prod_value = 0.0
    else:
        prod_row_pos = None
        prod_value = 0.0

    if prod_row_pos is None or prod_value == 0.0:
        # Prefer a production exchange closest to +1 in magnitude
        k = int(np.argmin(np.abs(np.abs(vals) - 1.0)))
        prod_row_pos = int(rows[k])
        prod_value = float(vals[k])

    if prod_row_pos not in pos_to_prod_id:
        raise KeyError(
            f"Product row position {prod_row_pos} not found in dicts.product "
            f"(activity_id={activity_id}, act_pos={act_pos})"
        )

    result = (int(pos_to_prod_id[prod_row_pos]), prod_value)
    ref_cache[activity_id] = result
    return result


def _reference_product_id_from_activity_id(lca_obj: Any, activity_id: int) -> int:
    """reference product id from activity id.

    :param lca_obj: Value for `lca_obj`.
    :type lca_obj: Any
    :param activity_id: Value for `activity_id`.
    :type activity_id: int
    :returns: Return value.
    :rtype: int"""
    prod_id, _ = _reference_product_from_activity_id(lca_obj, activity_id)
    return prod_id


def _extract_supply_fast_cached(
    supply_array: np.ndarray,
    act_ids: np.ndarray,
    positions: np.ndarray,
    min_amount: float,
) -> Dict[int, float]:
    """extract supply fast cached.

    :param supply_array: Value for `supply_array`.
    :type supply_array: np.ndarray
    :param act_ids: Value for `act_ids`.
    :type act_ids: np.ndarray
    :param positions: Value for `positions`.
    :type positions: np.ndarray
    :param min_amount: Value for `min_amount`.
    :type min_amount: float
    :returns: Return value.
    :rtype: Dict[int, float]"""
    supply = np.asarray(supply_array, dtype=np.float64)

    vals = supply[positions]
    m = vals != 0.0
    if not m.any():
        return {}

    a = act_ids[m]
    v = vals[m]
    return {int(ai): float(vi) for ai, vi in zip(a, v)}


def _extract_supply_fast(lca_obj: Any, min_amount: float) -> Dict[int, float]:
    """extract supply fast.

    :param lca_obj: Value for `lca_obj`.
    :type lca_obj: Any
    :param min_amount: Value for `min_amount`.
    :type min_amount: float
    :returns: Return value.
    :rtype: Dict[int, float]"""
    supply = np.asarray(lca_obj.supply_array).astype(np.float64, copy=False)

    # dicts.activity maps Trails act_idx -> position in supply_array
    act_map = getattr(lca_obj.dicts, "activity", None)
    if not act_map:
        return {}

    # Build arrays once
    act_ids = np.fromiter(act_map.keys(), dtype=np.int64, count=len(act_map))
    positions = np.fromiter(act_map.values(), dtype=np.int64, count=len(act_map))

    vals = supply[positions]
    m = vals != 0.0
    if not m.any():
        return {}

    act_ids = act_ids[m]
    vals = vals[m]

    # Convert back to python dict
    return {int(a): float(v) for a, v in zip(act_ids, vals)}


def _get_datapackage(
    dp_cache: dict[tuple[int, bool], Any],
    trails: Trails,
    year: int,
    zero_bio: bool,
    debug: bool,
) -> tuple[Any, dict[tuple, int], dict[tuple, int], list[tuple[int, int]]]:
    """get datapackage.

    :param dp_cache: Value for `dp_cache`.
    :type dp_cache: dict[tuple[int, bool], Any]
    :param trails: Value for `trails`.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :param zero_bio: Value for `zero_bio`.
    :type zero_bio: bool
    :param debug: Value for `debug`.
    :type debug: bool
    :returns: Return value.
    :rtype: tuple[Any, dict[tuple, int], dict[tuple, int], list[tuple[int, int]]]"""
    dp_key = (year, zero_bio)
    if dp_key not in dp_cache:
        dp, tech_idx, bio_idx, uncertain_params = (
            build_datapackage_for_year_from_trails(
                trails=trails,
                year=year,
                zero_biosphere=zero_bio,
                debug=debug,
            )
        )
        dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
    return dp_cache[dp_key]
