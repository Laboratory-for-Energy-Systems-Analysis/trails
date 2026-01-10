from __future__ import annotations

from typing import Any, Dict, Tuple, TYPE_CHECKING

import bw_processing as bwp
import numpy as np

if TYPE_CHECKING:
    from .trails import Trails


def _ij_from_coords(X_t: Any) -> tuple[np.ndarray, np.ndarray]:
    """Extract row/column indices from a sparse COO coordinate array.

    :param X_t: Sparse array with ``coords`` attribute.
    :type X_t: sparse.COO
    :returns: Tuple of ``(row_indices, col_indices)``.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
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
    """Resolve a calendar year to the closest available matrix label.

    :param trails: Trails instance with scenario labels and index.
    :type trails: Trails
    :param year: Calendar year to resolve.
    :type year: int
    :returns: Tuple of ``(label_for_matrix, scenario_index)``.
    :rtype: tuple[str, int]
    """
    label_for_matrix = str(year)
    if label_for_matrix not in trails.scenario_index:
        years = np.array([int(lbl) for lbl in trails.scenario_labels])
        idx = int(np.argmin(np.abs(years - year)))
        label_for_matrix = trails.scenario_labels[idx]


    t = trails.scenario_index[label_for_matrix]
    return label_for_matrix, t


def _select_metadata_label(trails: Trails, label_for_matrix: str) -> str:
    """Select the metadata label matching or nearest to the matrix label.

    :param trails: Trails instance with metadata labels.
    :type trails: Trails
    :param label_for_matrix: Scenario label used for the matrix slice.
    :type label_for_matrix: str
    :returns: Metadata label to use for indices.
    :rtype: str
    """
    if label_for_matrix in trails.activity_indices:
        return label_for_matrix
    return _nearest_metadata_label_for_year(trails, int(label_for_matrix))


def _make_bw_indices_rowcol(row_idx: np.ndarray, col_idx: np.ndarray) -> np.ndarray:
    """Build a bw_processing indices array from row and column indices.

    :param row_idx: Row indices array.
    :type row_idx: numpy.ndarray
    :param col_idx: Column indices array.
    :type col_idx: numpy.ndarray
    :returns: Structured indices array compatible with bw_processing.
    :rtype: numpy.ndarray
    """
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
    """Emit warnings for matrix indices missing from metadata.

    :param label_for_matrix: Scenario label for matrix slices.
    :param meta_label: Metadata label used for indices.
    :param A_act_idx: Activity indices from A.
    :param A_prod_idx: Product indices from A.
    :param B_flow_idx: Flow indices from B.
    :param act_meta: Activity metadata mapping.
    :param bio_meta: Biosphere metadata mapping.
    :type label_for_matrix: str
    :type meta_label: str
    :type A_act_idx: numpy.ndarray
    :type A_prod_idx: numpy.ndarray
    :type B_flow_idx: numpy.ndarray
    :type act_meta: dict
    :type bio_meta: dict
    """
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
    """Build technosphere and biosphere indices from metadata dictionaries.

    :param act_meta: Activity metadata mapping.
    :type act_meta: dict
    :param bio_meta: Biosphere metadata mapping.
    :type bio_meta: dict
    :returns: Tuple of ``(technosphere_indices, biosphere_indices)``.
    :rtype: tuple[dict[tuple, int], dict[tuple, int]]
    """
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
    """
    Pick the metadata scenario_label whose numeric year is closest to `year`.
    Used only as a fallback if the exact label is not present.
    """
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


def _reference_product_id_from_activity_id(lca_obj, activity_id: int) -> int:
    act_map = getattr(lca_obj.dicts, "activity", None)
    prod_map = getattr(lca_obj.dicts, "product", None)
    if not act_map or not prod_map:
        raise ValueError("LCA object missing dicts.activity or dicts.product")

    activity_id = int(activity_id)

    # Common case: activity ids are also valid product ids (diagonal identity)
    if activity_id in prod_map:
        return activity_id

    if activity_id not in act_map:
        raise KeyError(f"activity_id={activity_id} not in lca_obj.dicts.activity")

    act_pos = int(act_map[activity_id])
    pos_to_prod_id = {int(pos): int(pid) for pid, pos in prod_map.items()}

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

    # YOUR convention: production is positive (typically +1)
    prod_mask = vals > 0

    # Prefer +1-ish production
    if np.any(prod_mask):
        pos_rows = rows[prod_mask]
        pos_vals = vals[prod_mask]
        # closest to +1
        k = int(np.argmin(np.abs(pos_vals - 1.0)))
        prod_row_pos = int(pos_rows[k])
    else:
        # No positive entries found: fall back to largest magnitude entry for diagnostics
        k = int(np.argmax(np.abs(vals)))
        prod_row_pos = int(rows[k])

    if prod_row_pos not in pos_to_prod_id:
        raise KeyError(
            f"Product row position {prod_row_pos} not found in dicts.product "
            f"(activity_id={activity_id}, act_pos={act_pos})"
        )

    return int(pos_to_prod_id[prod_row_pos])


def _extract_supply_fast_cached(
    supply_array: np.ndarray,
    act_ids: np.ndarray,
    positions: np.ndarray,
    min_amount: float,
) -> Dict[int, float]:
    """Extract supply using precomputed (act_ids, positions) mapping for this solve_year."""
    supply = np.asarray(supply_array, dtype=np.float64)

    vals = supply[positions]
    m = np.abs(vals) > float(min_amount)
    if not m.any():
        return {}

    a = act_ids[m]
    v = vals[m]
    return {int(ai): float(vi) for ai, vi in zip(a, v)}


def _extract_supply_fast(lca_obj: Any, min_amount: float) -> Dict[int, float]:
    """Extract supply in Trails activity-index space using lca_obj.dicts.activity.

    Uses vectorized thresholding on the BW-ordered supply array.
    """
    supply = np.asarray(lca_obj.supply_array).astype(np.float64, copy=False)

    # dicts.activity maps Trails act_idx -> position in supply_array
    act_map = getattr(lca_obj.dicts, "activity", None)
    if not act_map:
        return {}

    # Build arrays once
    act_ids = np.fromiter(act_map.keys(), dtype=np.int64, count=len(act_map))
    positions = np.fromiter(act_map.values(), dtype=np.int64, count=len(act_map))

    vals = supply[positions]
    m = np.abs(vals) > float(min_amount)
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
    """Fetch or build a datapackage for a given year and biosphere setting.

    :param dp_cache: Cache mapping ``(year, zero_bio)`` to datapackage tuples.
    :type dp_cache: dict[tuple[int, bool], tuple]
    :param trails: Trails instance to build datapackages from.
    :type trails: Trails
    :param year: Calendar year to load.
    :type year: int
    :param zero_bio: Whether to zero biosphere emissions.
    :type zero_bio: bool
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: Datapackage tuple ``(dp, tech_idx, bio_idx, uncertain_params)``.
    :rtype: tuple
    """
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
