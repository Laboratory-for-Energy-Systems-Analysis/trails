# datapackage.py
from __future__ import annotations

from pathlib import PurePosixPath
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import numpy as np
import sparse

from .utils import _parse_float_or_none, _parse_int_or_none
from .temporal_distributions import TemporalExchange

import logging

logger = logging.getLogger(__name__)


def _parse_intish_or_none(value: object) -> int | None:
    """Parse an integer from values that may be formatted as strings or floats.

    :param value: Input value to parse.
    :type value: object
    :returns: Parsed integer or ``None`` for empty/invalid values.
    :rtype: int | None
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


# ----------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------
def _iter_inventory_resources(
    package: Any, filename: str
) -> Iterator[tuple[str, Any]]:
    """
    Yield (year_label, resource) for resources whose path matches:

        inventories/<model>/<pathway>/<year>/<filename>

    Because each datapackage contains exactly one model and one pathway,
    we ignore them and use <year> as the scenario_label.
    :param package: Frictionless data package to scan.
    :type package: Package
    :param filename: Inventory filename to match.
    :type filename: str
    :yields: Tuple of ``(year_label, resource)``.
    :rtype: tuple[str, Resource]
    """
    for res in package.resources:
        desc = res.descriptor or {}
        path_str = desc.get("path")
        if not path_str:
            continue

        path = PurePosixPath(path_str)
        parts = path.parts

        if (
            len(parts) >= 5
            and parts[0] == "inventories"
            and parts[-1].lower() == filename.lower()
        ):
            year_label = parts[3]  # inventories/<model>/<pathway>/<year>/...
            yield year_label, res


def _parse_temporal_exchange_row(
    row: Mapping[str, Any],
) -> TemporalExchange | None:
    """Parse a temporal exchange row into a TemporalExchange object.

    :param row: Row mapping from the inventory CSV.
    :type row: dict
    :returns: Parsed TemporalExchange or ``None`` if not specified.
    :rtype: TemporalExchange | None
    """
    dist_code = _parse_intish_or_none(row.get("temporal_distribution"))
    if dist_code is None:
        return None

    loc = _parse_float_or_none(row.get("temporal_loc"))
    scale = _parse_float_or_none(row.get("temporal_scale"))
    off_min = _parse_intish_or_none(row.get("temporal_min")) or 0
    off_max = _parse_intish_or_none(row.get("temporal_max")) or 0

    scale_mode = row.get("temporal_scale_mode")
    scale_base = _parse_float_or_none(row.get("temporal_scale_base", 1.0))
    scale_rate = _parse_float_or_none(row.get("temporal_scale_rate", 0.0))

    return TemporalExchange(
        distribution=dist_code,
        loc=loc,
        scale=scale,
        offset_min=int(off_min),
        offset_max=int(off_max),
        scale_mode=scale_mode,
        scale_base=scale_base,
        scale_rate=scale_rate,
    )


def _build_sparse_matrix(
    coords: Sequence[Sequence[int]],
    data: Sequence[float],
    shape: tuple[int, ...],
    idx_dtype: np.dtype,
    val_dtype: np.dtype,
) -> sparse.COO:
    """Build a sparse.COO matrix from coordinate lists and data.

    :param coords: Coordinate lists for each axis.
    :type coords: list[list[int]]
    :param data: Data values for nonzero entries.
    :type data: list[float]
    :param shape: Shape of the sparse matrix.
    :type shape: tuple[int, ...]
    :param idx_dtype: Numpy dtype for indices.
    :type idx_dtype: numpy.dtype
    :param val_dtype: Numpy dtype for values.
    :type val_dtype: numpy.dtype
    :returns: Sparse COO matrix.
    :rtype: sparse.COO
    """
    if data:
        coords_arrays = [np.array(axis, dtype=idx_dtype) for axis in coords]
        return sparse.COO(
            coords=coords_arrays,
            data=np.array(data, dtype=val_dtype),
            shape=shape,
        )

    empty_coords = [np.array([], dtype=idx_dtype) for _ in coords]
    return sparse.COO(
        coords=empty_coords,
        data=np.array([], dtype=val_dtype),
        shape=shape,
    )


def _parse_a_exchange_row(
    row: Mapping[str, Any],
    scenario_label: str,
    t: int,
    temporal_exchanges: Dict[tuple[str, int, int], TemporalExchange],
    A_coords: Dict[str, List[int]],
    max_indices: Dict[str, int],
) -> None:
    """Parse a technosphere exchange row and append to coordinates.

    :param row: Row mapping from the technosphere inventory.
    :type row: dict
    :param scenario_label: Scenario label for the row.
    :type scenario_label: str
    :param t: Scenario index to record.
    :type t: int
    :param temporal_exchanges: Mapping for temporal exchanges to update.
    :type temporal_exchanges: dict[tuple[str, int, int], TemporalExchange]
    :param A_coords: Coordinate accumulator for A.
    :type A_coords: dict[str, list]
    :param max_indices: Tracker for maximum indices.
    :type max_indices: dict[str, int]
    :returns: None.
    :rtype: None
    """
    act_idx = int(row["index of activity"])
    prod_idx = int(row["index of product"])
    value = float(row["value"])

    flip_flag = _parse_intish_or_none(row.get("flip")) or 0
    if flip_flag == 1:
        value = -value

    A_coords["t"].append(t)
    A_coords["i"].append(act_idx)
    A_coords["j"].append(prod_idx)
    A_coords["data"].append(value)

    max_indices["max_activity_idx_for_A"] = max(
        max_indices["max_activity_idx_for_A"], act_idx
    )
    max_indices["max_product_idx"] = max(max_indices["max_product_idx"], prod_idx)

    tex = _parse_temporal_exchange_row(row)
    if tex is not None:
        temporal_exchanges[(scenario_label, act_idx, prod_idx)] = tex


def _parse_b_exchange_row(
    row: Mapping[str, Any],
    scenario_label: str,
    t: int,
    temporal_biosphere_exchanges: Dict[tuple[str, int, int], TemporalExchange],
    B_coords: Dict[str, List[int]],
    max_indices: Dict[str, int],
) -> None:
    """Parse a biosphere exchange row and append to coordinates.

    :param row: Row mapping from the biosphere inventory.
    :type row: dict
    :param scenario_label: Scenario label for the row.
    :type scenario_label: str
    :param t: Scenario index to record.
    :type t: int
    :param temporal_biosphere_exchanges: Mapping for temporal exchanges to update.
    :type temporal_biosphere_exchanges: dict[tuple[str, int, int], TemporalExchange]
    :param B_coords: Coordinate accumulator for B.
    :type B_coords: dict[str, list]
    :param max_indices: Tracker for maximum indices.
    :type max_indices: dict[str, int]
    :returns: None.
    :rtype: None
    """
    act_idx = int(row["index of activity"])
    flow_idx = int(row["index of biosphere flow"])
    value = float(row["value"])

    B_coords["t"].append(t)
    B_coords["i"].append(act_idx)
    B_coords["j"].append(flow_idx)
    B_coords["data"].append(value)

    max_indices["max_activity_idx_for_B"] = max(
        max_indices["max_activity_idx_for_B"], act_idx
    )
    max_indices["max_flow_idx"] = max(max_indices["max_flow_idx"], flow_idx)

    tex = _parse_temporal_exchange_row(row)
    if tex is not None:
        temporal_biosphere_exchanges[(scenario_label, act_idx, flow_idx)] = tex


def _label_to_year(label: str) -> int:
    """
    Scenario labels may be:
      - "2050"
      - "model/pathway/2050"

    This extracts the trailing year.
    :param label: Scenario label string.
    :type label: str
    :returns: Extracted year.
    :rtype: int
    """
    tail = str(label).split("/")[-1]
    return int(tail)


def _years_and_sorted_indices(
    scenario_labels: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Return numeric years and indices that sort scenario labels by year.

    :param scenario_labels: Scenario labels to parse.
    :type scenario_labels: list[str]
    :returns: Tuple of years array and sorting indices.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    years = np.array([_label_to_year(lbl) for lbl in scenario_labels], dtype=int)
    order = np.argsort(years)
    return years[order], order


# ----------------------------------------------------------------------
# Public loading functions
# ----------------------------------------------------------------------
def load_matrices_from_package(
    package: Any,
    value_dtype: np.dtype = np.float32,
    index_dtype: np.dtype = np.int32,
    debug: bool = False,
) -> tuple[
    sparse.COO,
    sparse.COO,
    list[str],
    dict[str, int],
    dict[tuple[str, int, int], TemporalExchange],
    dict[tuple[str, int, int], TemporalExchange],
]:
    """Load technosphere and biosphere matrices from a data package.

    :param package: Frictionless data package to load.
    :type package: object
    :param value_dtype: Data type for matrix values.
    :type value_dtype: numpy.dtype
    :param index_dtype: Data type for matrix indices.
    :type index_dtype: numpy.dtype
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: Tuple of A matrix, B matrix, scenario labels, scenario index,
        temporal technosphere exchanges, and temporal biosphere exchanges.
    :rtype: tuple
    """
    scenario_labels: List[str] = []
    scenario_index: Dict[str, int] = {}

    def get_scenario_idx(label: str) -> int:
        """Get or assign an index for a scenario label.

        :param label: Scenario label to register.
        :type label: str
        :returns: Scenario index.
        :rtype: int
        """
        if label not in scenario_index:
            scenario_index[label] = len(scenario_labels)
            scenario_labels.append(label)
        return scenario_index[label]

    temporal_exchanges: Dict[Tuple[str, int, int], TemporalExchange] = {}
    temporal_biosphere_exchanges: Dict[Tuple[str, int, int], TemporalExchange] = {}

    A_coords = {"t": [], "i": [], "j": [], "data": []}
    B_coords = {"t": [], "i": [], "j": [], "data": []}

    max_indices = {
        "max_activity_idx_for_A": -1,
        "max_product_idx": -1,
        "max_activity_idx_for_B": -1,
        "max_flow_idx": -1,
    }

    for scenario_label, res in _iter_inventory_resources(package, "A_matrix.csv"):

        t = get_scenario_idx(scenario_label)
        for row in res.iter(keyed=True):
            _parse_a_exchange_row(
                row=row,
                scenario_label=scenario_label,
                t=t,
                temporal_exchanges=temporal_exchanges,
                A_coords=A_coords,
                max_indices=max_indices,
            )

    # ---------- Load all B_matrix.csv ----------
    for scenario_label, res in _iter_inventory_resources(package, "B_matrix.csv"):
        t = get_scenario_idx(scenario_label)
        for row in res.iter(keyed=True):
            _parse_b_exchange_row(
                row=row,
                scenario_label=scenario_label,
                t=t,
                temporal_biosphere_exchanges=temporal_biosphere_exchanges,
                B_coords=B_coords,
                max_indices=max_indices,
            )

    # ---------- Deduce shapes ----------
    n_scenarios = len(scenario_labels)
    n_activities = (
        max(
            max_indices["max_activity_idx_for_A"],
            max_indices["max_activity_idx_for_B"],
        )
        + 1
    )
    n_products = (
        max_indices["max_product_idx"] + 1 if max_indices["max_product_idx"] >= 0 else 0
    )
    n_flows = max_indices["max_flow_idx"] + 1 if max_indices["max_flow_idx"] >= 0 else 0

    idx_dtype = index_dtype
    val_dtype = value_dtype

    # ---------- Build sparse 3D matrices ----------
    A = _build_sparse_matrix(
        coords=[A_coords["t"], A_coords["i"], A_coords["j"]],
        data=A_coords["data"],
        shape=(n_scenarios, n_activities, n_products),
        idx_dtype=idx_dtype,
        val_dtype=val_dtype,
    )

    B = _build_sparse_matrix(
        coords=[B_coords["t"], B_coords["i"], B_coords["j"]],
        data=B_coords["data"],
        shape=(n_scenarios, n_activities, n_flows),
        idx_dtype=idx_dtype,
        val_dtype=val_dtype,
    )

    if debug:
        logger.info(
            "Datapackage: loaded A shape=%s nnz=%d | B shape=%s nnz=%d | scenarios=%d | temporal_exchanges=%d",
            getattr(A, "shape", None),
            int(getattr(A, "nnz", 0)),
            getattr(B, "shape", None),
            int(getattr(B, "nnz", 0)),
            len(scenario_labels),
            len(temporal_exchanges),
        )

    return (
        A,
        B,
        scenario_labels,
        scenario_index,
        temporal_exchanges,
        temporal_biosphere_exchanges,
    )


def interpolate_to_annual(
    A: sparse.COO,
    B: sparse.COO,
    scenario_labels: Sequence[str],
    value_dtype: np.dtype = np.float32,
    debug: bool = False,
) -> tuple[sparse.COO, sparse.COO, list[str], dict[str, int]]:
    """Interpolate A and B to annual slices between min and max year.

    :param A: 3D technosphere matrix with leading dimension = scenario/time.
    :type A: sparse.COO
    :param B: 3D biosphere matrix with leading dimension = scenario/time.
    :type B: sparse.COO
    :param scenario_labels: Scenario labels that are parseable as years.
    :type scenario_labels: Sequence[str]
    :param value_dtype: Target floating dtype.
    :type value_dtype: numpy.dtype
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: Interpolated A matrix, interpolated B matrix, annual labels,
        and label-to-index mapping.
    :rtype: tuple[sparse.COO, sparse.COO, list[str], dict[str, int]]
    """
    years_sorted, order = _years_and_sorted_indices(scenario_labels)

    if debug:
        logger.info("Datapackage: discovered inventory years=%s", years_sorted)

    A_sorted = A[order, :, :]
    B_sorted = B[order, :, :]

    new_As, new_Bs, new_labels = _interpolate_annual_slices(
        years_sorted=years_sorted,
        A_sorted=A_sorted,
        B_sorted=B_sorted,
        val_dtype=value_dtype,
    )

    A_interp = sparse.stack(new_As, axis=0).astype(value_dtype)
    B_interp = sparse.stack(new_Bs, axis=0).astype(value_dtype)
    new_index = {label: i for i, label in enumerate(new_labels)}

    return A_interp, B_interp, new_labels, new_index


def _interpolate_annual_slices(
    years_sorted: np.ndarray,
    A_sorted: sparse.COO,
    B_sorted: sparse.COO,
    val_dtype: np.dtype,
) -> tuple[list[sparse.COO], list[sparse.COO], list[str]]:
    """Interpolate sparse matrix slices to annual resolution.

    :param years_sorted: Sorted scenario years.
    :type years_sorted: numpy.ndarray
    :param A_sorted: Sorted technosphere slices.
    :type A_sorted: sparse.COO
    :param B_sorted: Sorted biosphere slices.
    :type B_sorted: sparse.COO
    :param val_dtype: Target value dtype.
    :type val_dtype: numpy.dtype
    :returns: Tuple of interpolated A slices, B slices, and labels.
    :rtype: tuple[list[sparse.COO], list[sparse.COO], list[str]]
    """
    new_As = []
    new_Bs = []
    new_labels: List[str] = []

    y0 = years_sorted[0]
    A0 = A_sorted[0].astype(val_dtype)
    B0 = B_sorted[0].astype(val_dtype)

    new_As.append(A0)
    new_Bs.append(B0)
    new_labels.append(str(y0))

    for k in range(len(years_sorted) - 1):
        y0 = years_sorted[k]
        y1 = years_sorted[k + 1]
        A0 = A_sorted[k]
        A1 = A_sorted[k + 1]
        B0 = B_sorted[k]
        B1 = B_sorted[k + 1]

        dt = y1 - y0
        if dt <= 0:
            continue

        for y in range(y0 + 1, y1 + 1):
            w = val_dtype((y - y0) / dt)
            one_minus_w = val_dtype(1.0) - w

            A_y = (one_minus_w * A0 + w * A1).astype(val_dtype)
            B_y = (one_minus_w * B0 + w * B1).astype(val_dtype)

            new_As.append(A_y)
            new_Bs.append(B_y)
            new_labels.append(str(y))

    return new_As, new_Bs, new_labels


def _load_activity_indices(package: Any) -> Dict[str, Dict[int, dict]]:
    """Load activity index metadata by scenario label.

    :param package: Frictionless data package.
    :type package: Package
    :returns: Mapping of scenario label to activity metadata.
    :rtype: dict[str, dict[int, dict]]
    """
    activity_indices: Dict[str, Dict[int, dict]] = {}

    for scenario_label, res in _iter_inventory_resources(package, "A_matrix_index.csv"):
        rows = res.read(keyed=True)

        mapping = {}
        for row in rows:
            try:
                idx = int(row["index"])
            except (KeyError, TypeError, ValueError):
                continue

            mapping[idx] = {
                "name": row.get("name"),
                "reference product": row.get("reference product"),
                "unit": row.get("unit"),
                "location": row.get("location"),
            }

        activity_indices[scenario_label] = mapping

    return activity_indices


def _load_biosphere_indices(package: Any) -> Dict[str, Dict[int, dict]]:
    """Load biosphere index metadata by scenario label.

    :param package: Frictionless data package.
    :type package: Package
    :returns: Mapping of scenario label to biosphere metadata.
    :rtype: dict[str, dict[int, dict]]
    """
    biosphere_indices: Dict[str, Dict[int, dict]] = {}

    for scenario_label, res in _iter_inventory_resources(package, "B_matrix_index.csv"):
        rows = res.read(keyed=True)

        mapping = {}
        for row in rows:
            idx = int(row["index"])

            mapping[idx] = {
                "name": row.get("name"),
                "compartment": row.get("compartment"),
                "subcompartment": row.get("subcompartment"),
                "unit": row.get("unit"),
            }

        biosphere_indices[scenario_label] = mapping

    return biosphere_indices


def load_indices_from_package(
    package: Any,
) -> tuple[Dict[str, Dict[int, dict]], Dict[str, Dict[int, dict]]]:
    """Load scenario-specific activity and biosphere index metadata.

    :param package: Frictionless data package.
    :type package: object
    :returns: Activity and biosphere index mappings keyed by scenario label.
    :rtype: tuple[dict[str, dict[int, dict]], dict[str, dict[int, dict]]]
    """
    activity_indices = _load_activity_indices(package)
    biosphere_indices = _load_biosphere_indices(package)
    return activity_indices, biosphere_indices
