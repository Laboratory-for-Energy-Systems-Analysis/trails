# datapackage.py

from pathlib import PurePosixPath
from typing import Dict, List, Tuple

import numpy as np
import sparse

from .utils import _parse_float_or_none, _parse_int_or_none
from .temporal_distributions import TemporalExchange

def _parse_intish_or_none(value):
    """Parse an integer from values that may be given as '3', '3.0', 3.0, etc."""
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
def _iter_inventory_resources(package, filename: str):
    """
    Yield (scenario_label, resource) for resources whose path matches:

        inventories/<scenario_label>/<filename>

    `filename` is typically "A.csv", "B.csv", "A_matrix_index.csv",
    or "B_matrix_index.csv".
    """
    for res in package.resources:
        desc = res.descriptor or {}
        path_str = desc.get("path")
        if not path_str:
            # inline or remote resources, skip
            continue

        path = PurePosixPath(path_str)
        parts = path.parts

        if (
            len(parts) >= 3
            and parts[0] == "inventories"
            and parts[-1].lower() == filename.lower()
        ):
            scenario_label = parts[1]  # e.g. "2005", "2050", "2100"
            yield scenario_label, res


def _years_and_sorted_indices(scenario_labels: List[str]):
    """Return numeric years and indices that sort scenario_labels by year."""
    years = np.array([int(lbl) for lbl in scenario_labels])
    order = np.argsort(years)
    return years[order], order


# ----------------------------------------------------------------------
# Public loading functions
# ----------------------------------------------------------------------
def load_matrices_from_package(
    package,
    value_dtype=np.float32,
    index_dtype=np.int32,
) -> Tuple[sparse.COO, sparse.COO, List[str], Dict[str, int], Dict[Tuple[str, int, int], TemporalExchange], Dict[Tuple[str, int, int], TemporalExchange]]:
    """
    Collect all technosphere and biosphere exchanges across scenarios and
    build sparse 3D matrices A and B.

    Returns
    -------
    A : sparse.COO
    B : sparse.COO
    scenario_labels : List[str]
    scenario_index : Dict[str, int]
    temporal_exchanges : Dict[(scenario_label, act_idx, prod_idx), TemporalExchange]
    temporal_biosphere_exchanges : Dict[(scenario_label, act_idx, bio_idx), TemporalExchange]
    """
    scenario_labels: List[str] = []
    scenario_index: Dict[str, int] = {}

    def get_scenario_idx(label: str) -> int:
        if label not in scenario_index:
            scenario_index[label] = len(scenario_labels)
            scenario_labels.append(label)
        return scenario_index[label]

    temporal_exchanges: Dict[Tuple[str, int, int], TemporalExchange] = {}
    temporal_biosphere_exchanges: Dict[Tuple[str, int, int], TemporalExchange] = {}

    # Collect triplets for A and B: (scenario, row, col, value)
    A_t, A_i, A_j, A_data = [], [], [], []
    B_t, B_i, B_j, B_data = [], [], [], []

    max_activity_idx_for_A = -1
    max_product_idx = -1
    max_activity_idx_for_B = -1
    max_flow_idx = -1

    for scenario_label, res in _iter_inventory_resources(package, "A.csv"):
        t = get_scenario_idx(scenario_label)
        for row in res.iter(keyed=True):
            act_idx = int(row["index of activity"])
            prod_idx = int(row["index of product"])
            value = float(row["value"])

            A_t.append(t)
            A_i.append(act_idx)
            A_j.append(prod_idx)
            A_data.append(value)

            if act_idx > max_activity_idx_for_A:
                max_activity_idx_for_A = act_idx
            if prod_idx > max_product_idx:
                max_product_idx = prod_idx

            td_raw = row.get("temporal_distribution")
            if td_raw is not None and str(td_raw).strip() != "":
                try:
                    dist_code = int(float(td_raw))
                except ValueError:
                    dist_code = None
            else:
                dist_code = None

            if dist_code is not None:
                loc = _parse_float_or_none(row.get("temporal_loc"))
                scale = _parse_float_or_none(row.get("temporal_scale"))
                off_min = _parse_intish_or_none(row.get("temporal_min")) or 0
                off_max = _parse_intish_or_none(row.get("temporal_max")) or 0

                temporal_exchanges[(scenario_label, act_idx, prod_idx)] = TemporalExchange(
                    distribution=dist_code,
                    loc=loc,
                    scale=scale,
                    offset_min=int(off_min),
                    offset_max=int(off_max),
                )

    # ---------- Load all B.csv ----------
    for scenario_label, res in _iter_inventory_resources(package, "B.csv"):
        t = get_scenario_idx(scenario_label)
        for row in res.iter(keyed=True):
            act_idx = int(row["index of activity"])
            flow_idx = int(row["index of biosphere flow"])
            value = float(row["value"])

            B_t.append(t)
            B_i.append(act_idx)
            B_j.append(flow_idx)
            B_data.append(value)

            if act_idx > max_activity_idx_for_B:
                max_activity_idx_for_B = act_idx
            if flow_idx > max_flow_idx:
                max_flow_idx = flow_idx

            # Optional temporal metadata for biosphere exchanges (same columns as A.csv)
            td_raw = row.get("temporal_distribution")
            if td_raw is not None and str(td_raw).strip() != "":
                dist_code = _parse_intish_or_none(td_raw)
            else:
                dist_code = None

            if dist_code is not None:
                loc = _parse_float_or_none(row.get("temporal_loc"))
                scale = _parse_float_or_none(row.get("temporal_scale"))
                off_min = _parse_intish_or_none(row.get("temporal_min")) or 0
                off_max = _parse_intish_or_none(row.get("temporal_max")) or 0

                temporal_biosphere_exchanges[(scenario_label, act_idx, flow_idx)] = TemporalExchange(
                    distribution=dist_code,
                    loc=loc,
                    scale=scale,
                    offset_min=int(off_min),
                    offset_max=int(off_max),
                )

    # ---------- Deduce shapes ----------
    n_scenarios = len(scenario_labels)
    n_activities = max(max_activity_idx_for_A, max_activity_idx_for_B) + 1
    n_products = max_product_idx + 1 if max_product_idx >= 0 else 0
    n_flows = max_flow_idx + 1 if max_flow_idx >= 0 else 0

    idx_dtype = index_dtype
    val_dtype = value_dtype

    # ---------- Build sparse 3D matrices ----------
    if A_data:
        A = sparse.COO(
            coords=[
                np.array(A_t, dtype=idx_dtype),
                np.array(A_i, dtype=idx_dtype),
                np.array(A_j, dtype=idx_dtype),
            ],
            data=np.array(A_data, dtype=val_dtype),
            shape=(n_scenarios, n_activities, n_products),
        )
    else:
        A = sparse.COO(
            coords=[
                np.array([], dtype=idx_dtype),
                np.array([], dtype=idx_dtype),
                np.array([], dtype=idx_dtype),
            ],
            data=np.array([], dtype=val_dtype),
            shape=(n_scenarios, n_activities, n_products),
        )

    if B_data:
        B = sparse.COO(
            coords=[
                np.array(B_t, dtype=idx_dtype),
                np.array(B_i, dtype=idx_dtype),
                np.array(B_j, dtype=idx_dtype),
            ],
            data=np.array(B_data, dtype=val_dtype),
            shape=(n_scenarios, n_activities, n_flows),
        )
    else:
        B = sparse.COO(
            coords=[
                np.array([], dtype=idx_dtype),
                np.array([], dtype=idx_dtype),
                np.array([], dtype=idx_dtype),
            ],
            data=np.array([], dtype=val_dtype),
            shape=(n_scenarios, n_activities, n_flows),
        )

    return A, B, scenario_labels, scenario_index, temporal_exchanges, temporal_biosphere_exchanges



def interpolate_to_annual(
    A: sparse.COO,
    B: sparse.COO,
    scenario_labels: List[str],
    value_dtype=np.float32,
):
    """
    Linearly interpolate A and B to annual slices between min and max year.

    Parameters
    ----------
    A, B : sparse.COO
        3D matrices with leading dimension = scenario/time.
    scenario_labels : list of str
        Scenario labels that are parseable as years.
    value_dtype :
        Target floating dtype.

    Returns
    -------
    A_interp, B_interp : sparse.COO
        Interpolated 3D matrices (one slice per year).
    new_labels : List[str]
        List of yearly labels (str).
    new_index : Dict[str, int]
        Mapping from year label to index.
    """
    years_sorted, order = _years_and_sorted_indices(scenario_labels)

    # Reorder A and B to chronological order
    A_sorted = A[order, :, :]
    B_sorted = B[order, :, :]

    new_As = []
    new_Bs = []
    new_labels: List[str] = []

    val_dtype = value_dtype

    # Start with the first year exactly as in the data
    y0 = years_sorted[0]
    A0 = A_sorted[0].astype(val_dtype)
    B0 = B_sorted[0].astype(val_dtype)

    new_As.append(A0)
    new_Bs.append(B0)
    new_labels.append(str(y0))

    # Loop over each subsequent *provided* year
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

        # Fill intermediate years y0+1, ..., y1 (inclusive)
        for y in range(y0 + 1, y1 + 1):
            w = val_dtype((y - y0) / dt)
            one_minus_w = val_dtype(1.0) - w

            A_y = one_minus_w * A0 + w * A1
            B_y = one_minus_w * B0 + w * B1

            A_y = A_y.astype(val_dtype)
            B_y = B_y.astype(val_dtype)

            new_As.append(A_y)
            new_Bs.append(B_y)
            new_labels.append(str(y))

    A_interp = sparse.stack(new_As, axis=0).astype(val_dtype)
    B_interp = sparse.stack(new_Bs, axis=0).astype(val_dtype)
    new_index = {label: i for i, label in enumerate(new_labels)}

    return A_interp, B_interp, new_labels, new_index


def load_indices_from_package(package):
    """
    Load scenario-specific dictionaries mapping index -> activity metadata
    and index -> biosphere-flow metadata.

    Returns
    -------
    activity_indices : {scenario_label: {index: row_dict}}
    biosphere_indices : {scenario_label: {index: row_dict}}
    """
    activity_indices: Dict[str, Dict[int, dict]] = {}
    biosphere_indices: Dict[str, Dict[int, dict]] = {}

    # ---- Activities (A_matrix_index.csv) ----
    for scenario_label, res in _iter_inventory_resources(package, "A_matrix_index.csv"):
        rows = res.read(keyed=True)

        mapping = {}
        for row in rows:
            # headers: name;reference product;unit;location;index
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

    # ---- Biosphere flows (B_matrix_index.csv) ----
    for scenario_label, res in _iter_inventory_resources(package, "B_matrix_index.csv"):
        rows = res.read(keyed=True)

        mapping = {}
        for row in rows:
            # headers: name;compartment;subcompartment;unit;index
            try:
                idx = int(row["index"])
            except (KeyError, TypeError, ValueError):
                continue

            mapping[idx] = {
                "name": row.get("name"),
                "compartment": row.get("compartment"),
                "subcompartment": row.get("subcompartment"),
                "unit": row.get("unit"),
            }

        biosphere_indices[scenario_label] = mapping

    return activity_indices, biosphere_indices


