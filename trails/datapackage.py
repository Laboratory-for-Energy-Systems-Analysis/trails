# datapackage.py
from __future__ import annotations

from pathlib import PurePosixPath
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

import numpy as np
import sparse

from .utils import _parse_float_or_none, _parse_int_or_none
from .temporal_distributions import TemporalExchange

import logging

import os
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)


from pathlib import Path
import os


def _resource_abspath(package: Any, res: Any) -> Path:
    """
    Resolve a resource to an absolute filesystem path, using Frictionless'
    resolved source if available; otherwise fall back to package.basepath.
    """
    # 1) Best: frictionless-resolved source (often absolute or correctly rooted)
    src = getattr(res, "source", None)
    if src:
        try:
            p = Path(str(src))
            if p.exists():
                return p.resolve()
        except Exception:
            pass

    # 2) Next: package.basepath is the canonical base for relative resource paths
    base = getattr(package, "basepath", None)
    if base:
        try:
            basep = Path(str(base))
            # basepath is usually a directory; sometimes a file path
            if basep.is_file():
                basep = basep.parent
            desc = res.descriptor or {}
            rel = desc.get("path")
            if not rel:
                raise ValueError("Resource has no descriptor['path']")
            p = basep / rel
            if p.exists():
                return p.resolve()
        except Exception:
            pass

    # 3) Last resort: package.path (path to datapackage.json) if present
    pkg_path = getattr(package, "path", None)
    if pkg_path:
        pkgp = Path(str(pkg_path))
        basep = pkgp.parent if pkgp.is_file() else pkgp
        desc = res.descriptor or {}
        rel = desc.get("path")
        if not rel:
            raise ValueError("Resource has no descriptor['path']")
        p = basep / rel
        return p.resolve()

    # 4) Fail loudly with useful diagnostics
    desc = res.descriptor or {}
    raise FileNotFoundError(
        f"Could not resolve resource path. "
        f"descriptor.path={desc.get('path')!r}, resource.source={getattr(res, 'source', None)!r}, "
        f"package.basepath={getattr(package, 'basepath', None)!r}, package.path={getattr(package, 'path', None)!r}"
    )


import numpy as np
import pandas as pd

# Columns we expect in any matrix file (A)
BASE_COLS_A = [
    "index of activity",
    "index of product",
    "value",
    "uncertainty type",
    "loc",
    "scale",
    "shape",
    "minimum",
    "maximum",
    "negative",
    "flip",
]

# Optional temporal columns (may or may not be present; may include temporal_amount_source)
TEMPORAL_COLS = [
    "temporal_distribution",
    "temporal_loc",
    "temporal_scale",
    "temporal_min",
    "temporal_max",
    "temporal_amount_source",
]

# For B you likely have "index of biosphere flow" instead of "index of product"
BASE_COLS_B = [
    "index of activity",
    "index of biosphere flow",
    "value",
    "uncertainty type",
    "loc",
    "scale",
    "shape",
    "minimum",
    "maximum",
    "negative",
    "flip",
]


def _read_matrix_csv_fast(csv_path, kind="A"):
    """
    Robust CSV reader that keeps temporal columns and safely casts index columns to int.

    - Reads everything as string first to prevent pandas from shifting/casting surprises.
    - Then explicitly casts indices to int64 and numeric columns to float64.
    - Blanks are allowed in numeric columns; they become NaN (or 0 for flags).
    """
    if kind == "A":
        base = BASE_COLS_A
    elif kind == "B":
        base = BASE_COLS_B
    else:
        raise ValueError(f"Unknown kind={kind}")

    # Read header first to see what we have
    header = pd.read_csv(csv_path, sep=";", nrows=0)
    cols_present = list(header.columns)

    # Build the list of columns to read: required base + any temporal columns that exist
    cols_to_read = [c for c in base if c in cols_present]
    if len(cols_to_read) != len(base):
        missing = [c for c in base if c not in cols_present]
        raise ValueError(f"Missing required columns in {csv_path}: {missing}")

    cols_to_read += [c for c in TEMPORAL_COLS if c in cols_present]

    # Read selected columns, but as strings (object) so we can control casting
    df = pd.read_csv(
        csv_path,
        sep=";",
        usecols=cols_to_read,
        dtype=str,
        keep_default_na=False,  # keep blanks as "" so we can decide how to treat them
        na_filter=False,
        engine="c",
    )

    # Strip whitespace everywhere (cheap insurance)
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    # ---- Cast index columns (must be integer-like) ----
    index_cols = ["index of activity"]
    if kind == "A":
        index_cols.append("index of product")
    else:
        index_cols.append("index of biosphere flow")

    for col in index_cols:
        # Blank indices are not allowed; treat "" as NaN then fail with a useful error
        s = df[col].replace("", np.nan)
        arr = pd.to_numeric(s, errors="coerce").to_numpy(dtype=np.float64)

        if np.isnan(arr).any():
            bad_rows = np.where(np.isnan(arr))[0][:10]
            examples = df.loc[bad_rows, col].tolist()
            raise ValueError(
                f"Invalid/missing integer in '{col}' in {csv_path}. "
                f"Row examples (raw): {examples}"
            )

        rounded = np.rint(arr)
        if not np.allclose(arr, rounded, atol=0.0, rtol=0.0):
            bad = arr[arr != rounded][:10]
            raise ValueError(
                f"Non-integer values found in '{col}' in {csv_path}. Examples: {bad}"
            )

        df[col] = rounded.astype(np.int64)

    # ---- Cast numeric float columns (allow blanks) ----
    float_cols = [
        "value",
        "loc",
        "scale",
        "shape",
        "minimum",
        "maximum",
        "temporal_loc",
        "temporal_scale",
        "temporal_min",
        "temporal_max",
    ]
    for col in float_cols:
        if col not in df.columns:
            continue
        s = df[col]
        s = s.where(s != "", np.nan)
        df[col] = pd.to_numeric(s, errors="coerce").astype(np.float64)

    # ---- Cast int-like parameter columns where appropriate ----
    int_cols = ["uncertainty type"]
    # temporal_distribution might be int-coded; if you use strings, remove this
    if "temporal_distribution" in df.columns:
        int_cols.append("temporal_distribution")

    for col in int_cols:
        if col not in df.columns:
            continue

        s = df[col]
        s = s.where(s != "", 0)

        x = pd.to_numeric(s, errors="coerce")
        if x.isna().any():
            # keep your current permissive behavior
            x = x.fillna(0)

        df[col] = x.astype(np.int64)

    # ---- Cast flag columns (allow blanks -> 0) ----
    for col in ("negative", "flip"):
        if col not in df.columns:
            continue

        s = df[col]
        s = s.where(s != "", 0)

        x = pd.to_numeric(s, errors="coerce").fillna(0).astype(np.int64)
        df[col] = x

    return df


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
def _iter_inventory_resources(package: Any, filename: str) -> Iterator[tuple[str, Any]]:
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

    amount_source = _parse_amount_source(row.get("temporal_amount_source"))

    return TemporalExchange(
        distribution=dist_code,
        loc=loc,
        scale=scale,
        offset_min=int(off_min),
        offset_max=int(off_max),
        amount_source=amount_source,
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
    """Load technosphere (A) and biosphere (B) matrices from a data package.

    Optimizations vs. baseline:
      - Accumulate coords/data as numpy blocks, concatenate once (avoids millions of Python objects)
      - TD parsing avoids dict(zip(...)) per row
      - Still uses your existing _read_matrix_csv_fast for robustness
    """
    scenario_labels: list[str] = []
    scenario_index: dict[str, int] = {}

    def get_scenario_idx(label: str) -> int:
        idx = scenario_index.get(label)
        if idx is None:
            idx = len(scenario_labels)
            scenario_index[label] = idx
            scenario_labels.append(label)
        return idx

    temporal_exchanges: dict[tuple[str, int, int], TemporalExchange] = {}
    temporal_biosphere_exchanges: dict[tuple[str, int, int], TemporalExchange] = {}

    # Store blocks instead of Python lists of scalars
    A_t_blocks: list[np.ndarray] = []
    A_i_blocks: list[np.ndarray] = []
    A_j_blocks: list[np.ndarray] = []
    A_v_blocks: list[np.ndarray] = []

    B_t_blocks: list[np.ndarray] = []
    B_i_blocks: list[np.ndarray] = []
    B_j_blocks: list[np.ndarray] = []
    B_v_blocks: list[np.ndarray] = []

    max_activity_A = -1
    max_product = -1
    max_activity_B = -1
    max_flow = -1

    # ---- Small helper: parse TD row fields without dict construction ----
    def _parse_td_fields(
        dist,
        loc,
        scale,
        off_min,
        off_max,
        amount_source,
    ) -> TemporalExchange | None:
        # distribution required
        if dist is None:
            return None
        if isinstance(dist, str):
            s = dist.strip()
            if s == "":
                return None
            try:
                dist_code = int(float(s))
            except Exception:
                return None
        else:
            # numeric (int/float)
            try:
                if np.isnan(dist):  # type: ignore[arg-type]
                    return None
            except Exception:
                pass
            try:
                dist_code = int(dist)
            except Exception:
                try:
                    dist_code = int(float(dist))
                except Exception:
                    return None

        # amount_source optional
        src = "port"
        if amount_source is not None:
            try:
                if np.isnan(amount_source):  # type: ignore[arg-type]
                    amount_source = None
            except Exception:
                pass
            ssrc = str(amount_source).strip().lower() if amount_source is not None else ""
            if ssrc:
                if ssrc not in {"port", "matrix"}:
                    raise ValueError(f"Unknown temporal_amount_source: {ssrc}")
                src = ssrc

        def _to_int(x, default=0) -> int:
            if x is None:
                return default
            if isinstance(x, str):
                sx = x.strip()
                if sx == "":
                    return default
                try:
                    return int(float(sx))
                except Exception:
                    return default
            try:
                if np.isnan(x):  # type: ignore[arg-type]
                    return default
            except Exception:
                pass
            try:
                return int(x)
            except Exception:
                try:
                    return int(float(x))
                except Exception:
                    return default

        # loc/scale can be NaN
        def _to_float_or_none(x):
            if x is None:
                return None
            if isinstance(x, str):
                sx = x.strip()
                if sx == "":
                    return None
                try:
                    return float(sx)
                except Exception:
                    return None
            try:
                if np.isnan(x):  # type: ignore[arg-type]
                    return None
            except Exception:
                pass
            try:
                return float(x)
            except Exception:
                return None

        return TemporalExchange(
            distribution=dist_code,
            loc=_to_float_or_none(loc),
            scale=_to_float_or_none(scale),
            offset_min=_to_int(off_min, 0),
            offset_max=_to_int(off_max, 0),
            amount_source=src,
        )

    # ---------- Load all A_matrix.csv ----------
    for scenario_label, res in _iter_inventory_resources(package, "A_matrix.csv"):
        t = get_scenario_idx(scenario_label)

        csv_path = _resource_abspath(package, res)
        df = _read_matrix_csv_faster(csv_path, kind="A")

        n = int(len(df))
        if n:
            act = df["index of activity"].to_numpy(dtype=np.int64, copy=False)
            prod = df["index of product"].to_numpy(dtype=np.int64, copy=False)
            val = df["value"].to_numpy(dtype=np.float64, copy=False)

            # Apply flip if present
            if "flip" in df.columns:
                flip = df["flip"].to_numpy(dtype=np.int64, copy=False)
                m = flip == 1
                if m.any():
                    val = val.copy()
                    val[m] *= -1.0

            # Block append (avoid Python scalar lists)
            A_t_blocks.append(np.full(n, t, dtype=index_dtype))
            A_i_blocks.append(act.astype(index_dtype, copy=False))
            A_j_blocks.append(prod.astype(index_dtype, copy=False))
            A_v_blocks.append(val.astype(value_dtype, copy=False))

            max_activity_A = max(max_activity_A, int(act.max()))
            max_product = max(max_product, int(prod.max()))

        # Temporal exchanges (A)
        if "temporal_distribution" in df.columns:
            td = df["temporal_distribution"].to_numpy(copy=False)

            # Works with your robust reader: td is int64 (blanks -> 0), so filter on != 0.
            # If you ever switch to string TDs, this still behaves well for object dtype.
            if td.dtype.kind in ("U", "S", "O"):
                mask_td = td != ""
            else:
                mask_td = td != 0

            if mask_td.any():
                act_td = df.loc[mask_td, "index of activity"].to_numpy(
                    np.int64, copy=False
                )
                prod_td = df.loc[mask_td, "index of product"].to_numpy(
                    np.int64, copy=False
                )
                dist_td = df.loc[mask_td, "temporal_distribution"].to_numpy(copy=False)

                loc_td = (
                    df.loc[mask_td, "temporal_loc"].to_numpy(copy=False)
                    if "temporal_loc" in df.columns
                    else None
                )
                scale_td = (
                    df.loc[mask_td, "temporal_scale"].to_numpy(copy=False)
                    if "temporal_scale" in df.columns
                    else None
                )
                min_td = (
                    df.loc[mask_td, "temporal_min"].to_numpy(copy=False)
                    if "temporal_min" in df.columns
                    else None
                )
                max_td = (
                    df.loc[mask_td, "temporal_max"].to_numpy(copy=False)
                    if "temporal_max" in df.columns
                    else None
                )
                src_td = (
                    df.loc[mask_td, "temporal_amount_source"].to_numpy(copy=False)
                    if "temporal_amount_source" in df.columns
                    else None
                )

                for irow in range(act_td.size):
                    tex = _parse_td_fields(
                        dist_td[irow],
                        None if loc_td is None else loc_td[irow],
                        None if scale_td is None else scale_td[irow],
                        None if min_td is None else min_td[irow],
                        None if max_td is None else max_td[irow],
                        None if src_td is None else src_td[irow],
                    )
                    if tex is not None:
                        temporal_exchanges[
                            (scenario_label, int(act_td[irow]), int(prod_td[irow]))
                        ] = tex

    # ---------- Load all B_matrix.csv ----------
    for scenario_label, res in _iter_inventory_resources(package, "B_matrix.csv"):
        t = get_scenario_idx(scenario_label)

        csv_path = _resource_abspath(package, res)
        df = _read_matrix_csv_fast(csv_path, kind="B")

        n = int(len(df))
        if n:
            act = df["index of activity"].to_numpy(dtype=np.int64, copy=False)
            flow = df["index of biosphere flow"].to_numpy(dtype=np.int64, copy=False)
            val = df["value"].to_numpy(dtype=np.float64, copy=False)

            B_t_blocks.append(np.full(n, t, dtype=index_dtype))
            B_i_blocks.append(act.astype(index_dtype, copy=False))
            B_j_blocks.append(flow.astype(index_dtype, copy=False))
            B_v_blocks.append(val.astype(value_dtype, copy=False))

            max_activity_B = max(max_activity_B, int(act.max()))
            max_flow = max(max_flow, int(flow.max()))

        # Temporal exchanges (B)
        if "temporal_distribution" in df.columns:
            td = df["temporal_distribution"].to_numpy(copy=False)
            if td.dtype.kind in ("U", "S", "O"):
                mask_td = td != ""
            else:
                mask_td = td != 0

            if mask_td.any():
                act_td = df.loc[mask_td, "index of activity"].to_numpy(
                    np.int64, copy=False
                )
                flow_td = df.loc[mask_td, "index of biosphere flow"].to_numpy(
                    np.int64, copy=False
                )
                dist_td = df.loc[mask_td, "temporal_distribution"].to_numpy(copy=False)

                loc_td = (
                    df.loc[mask_td, "temporal_loc"].to_numpy(copy=False)
                    if "temporal_loc" in df.columns
                    else None
                )
                scale_td = (
                    df.loc[mask_td, "temporal_scale"].to_numpy(copy=False)
                    if "temporal_scale" in df.columns
                    else None
                )
                min_td = (
                    df.loc[mask_td, "temporal_min"].to_numpy(copy=False)
                    if "temporal_min" in df.columns
                    else None
                )
                max_td = (
                    df.loc[mask_td, "temporal_max"].to_numpy(copy=False)
                    if "temporal_max" in df.columns
                    else None
                )
                src_td = (
                    df.loc[mask_td, "temporal_amount_source"].to_numpy(copy=False)
                    if "temporal_amount_source" in df.columns
                    else None
                )

                for irow in range(act_td.size):
                    tex = _parse_td_fields(
                        dist_td[irow],
                        None if loc_td is None else loc_td[irow],
                        None if scale_td is None else scale_td[irow],
                        None if min_td is None else min_td[irow],
                        None if max_td is None else max_td[irow],
                        None if src_td is None else src_td[irow],
                    )
                    if tex is not None:
                        temporal_biosphere_exchanges[
                            (scenario_label, int(act_td[irow]), int(flow_td[irow]))
                        ] = tex

    # ---------- Deduce shapes ----------
    n_scenarios = len(scenario_labels)
    n_activities = max(max_activity_A, max_activity_B) + 1
    n_products = (max_product + 1) if max_product >= 0 else 0
    n_flows = (max_flow + 1) if max_flow >= 0 else 0

    # ---------- Concatenate blocks ----------
    def _cat(blocks: list[np.ndarray], dtype: np.dtype) -> np.ndarray:
        if not blocks:
            return np.array([], dtype=dtype)
        out = np.concatenate(blocks)
        return out.astype(dtype, copy=False)

    A_t = _cat(A_t_blocks, index_dtype)
    A_i = _cat(A_i_blocks, index_dtype)
    A_j = _cat(A_j_blocks, index_dtype)
    A_v = _cat(A_v_blocks, value_dtype)

    B_t = _cat(B_t_blocks, index_dtype)
    B_i = _cat(B_i_blocks, index_dtype)
    B_j = _cat(B_j_blocks, index_dtype)
    B_v = _cat(B_v_blocks, value_dtype)

    # ---------- Build sparse 3D matrices ----------
    A = sparse.COO(
        coords=[A_t, A_i, A_j],
        data=A_v,
        shape=(n_scenarios, n_activities, n_products),
    )

    B = sparse.COO(
        coords=[B_t, B_i, B_j],
        data=B_v,
        shape=(n_scenarios, n_activities, n_flows),
    )

    if debug:
        logger.info(
            "Datapackage: loaded A shape=%s nnz=%d | B shape=%s nnz=%d | scenarios=%d | temporal_exchanges=%d | temporal_biosphere_exchanges=%d",
            getattr(A, "shape", None),
            int(getattr(A, "nnz", 0)),
            getattr(B, "shape", None),
            int(getattr(B, "nnz", 0)),
            len(scenario_labels),
            len(temporal_exchanges),
            len(temporal_biosphere_exchanges),
        )

    return (
        A,
        B,
        scenario_labels,
        scenario_index,
        temporal_exchanges,
        temporal_biosphere_exchanges,
    )


def _read_matrix_csv_faster(csv_path: str, kind="A") -> pd.DataFrame:
    if kind == "A":
        base = BASE_COLS_A
        idx2 = "index of product"
    else:
        base = BASE_COLS_B
        idx2 = "index of biosphere flow"

    header = pd.read_csv(csv_path, sep=";", nrows=0)
    cols_present = list(header.columns)

    cols_to_read = [c for c in base if c in cols_present]
    cols_to_read += [c for c in TEMPORAL_COLS if c in cols_present]

    dtype_map = {
        "index of activity": "int64",
        idx2: "int64",
        "uncertainty type": "int64",
        "negative": "int64",
        "flip": "int64",
        # floats
        "value": "float64",
        "loc": "float64",
        "scale": "float64",
        "shape": "float64",
        "minimum": "float64",
        "maximum": "float64",
        "temporal_loc": "float64",
        "temporal_scale": "float64",
        "temporal_min": "float64",  # if these are really integers, keep float64 then cast later
        "temporal_max": "float64",
        "temporal_distribution": "float64",  # blanks -> NaN
        # temporal_amount_source left as object
    }

    # only include dtype entries for existing columns
    dtype_map = {
        k: v
        for k, v in dtype_map.items()
        if k in cols_to_read and k != "temporal_amount_source"
    }

    df = pd.read_csv(
        csv_path,
        sep=";",
        usecols=cols_to_read,
        dtype=dtype_map,
        na_values=["", " "],
        keep_default_na=True,
        engine="c",
    )

    # Fill flags if missing
    for col in ("negative", "flip", "uncertainty type"):
        if col in df.columns:
            df[col] = df[col].fillna(0).astype(np.int64)

    return df


def _parse_temporal_exchange_fields(
    dist, loc, scale, off_min, off_max, amount_source
) -> TemporalExchange | None:
    # dist can be "" or 0 or NaN depending on your read path
    if dist is None:
        return None
    if isinstance(dist, str) and dist.strip() == "":
        return None
    try:
        dist_code = int(dist)
    except Exception:
        try:
            dist_code = int(float(dist))
        except Exception:
            return None

    # amount_source: only parse if column exists, else default
    src = "port"
    if amount_source is not None:
        s = str(amount_source).strip().lower()
        if s:
            if s not in {"port", "matrix"}:
                raise ValueError(f"Unknown temporal_amount_source: {s}")
            src = s

    # offsets
    def _to_int(x, default=0):
        if x is None:
            return default
        if isinstance(x, str) and x.strip() == "":
            return default
        try:
            return int(x)
        except Exception:
            try:
                return int(float(x))
            except Exception:
                return default

    return TemporalExchange(
        distribution=dist_code,
        loc=(
            None
            if (loc is None or (isinstance(loc, float) and np.isnan(loc)))
            else float(loc)
        ),
        scale=(
            None
            if (scale is None or (isinstance(scale, float) and np.isnan(scale)))
            else float(scale)
        ),
        offset_min=_to_int(off_min, 0),
        offset_max=_to_int(off_max, 0),
        amount_source=src,
    )


def interpolate_to_annual(
    A: sparse.COO,
    B: sparse.COO,
    scenario_labels: Sequence[str],
    value_dtype: np.dtype = np.float32,
    debug: bool = False,
) -> tuple[sparse.COO, sparse.COO, list[str], dict[str, int]]:
    years_sorted, order = _years_and_sorted_indices(scenario_labels)

    if debug:
        logger.info("Datapackage: discovered inventory years=%s", years_sorted)

    idx_dtype = A.coords[0].dtype  # typically int32

    new_As: list[sparse.COO] = []
    new_Bs: list[sparse.COO] = []
    new_labels: list[str] = []

    # First anchor
    i0 = int(order[0])
    new_As.append(A[i0].astype(value_dtype))
    new_Bs.append(B[i0].astype(value_dtype))
    new_labels.append(str(int(years_sorted[0])))

    for k in range(len(years_sorted) - 1):
        y0 = int(years_sorted[k])
        y1 = int(years_sorted[k + 1])
        dt = y1 - y0
        if dt <= 0:
            continue

        i0 = int(order[k])
        i1 = int(order[k + 1])

        A0 = A[i0]
        A1 = A[i1]
        B0 = B[i0]
        B1 = B[i1]

        for y in range(y0 + 1, y1 + 1):
            w = float(y - y0) / float(dt)
            new_As.append(
                _interp_slice_union_vectorized(A0, A1, w, idx_dtype, value_dtype)
            )
            new_Bs.append(
                _interp_slice_union_vectorized(B0, B1, w, idx_dtype, value_dtype)
            )
            new_labels.append(str(y))

    A_interp = sparse.stack(new_As, axis=0).astype(value_dtype)
    B_interp = sparse.stack(new_Bs, axis=0).astype(value_dtype)
    new_index = {label: i for i, label in enumerate(new_labels)}
    return A_interp, B_interp, new_labels, new_index


def _concat_or_empty(arrs: list[np.ndarray], dtype: np.dtype) -> np.ndarray:
    if not arrs:
        return np.array([], dtype=dtype)
    return np.concatenate(arrs).astype(dtype, copy=False)


def _interpolate_annual_slices(
    years_sorted: np.ndarray,
    A_sorted: sparse.COO,
    B_sorted: sparse.COO,
    val_dtype: np.dtype,
) -> tuple[list[sparse.COO], list[sparse.COO], list[str]]:
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


def _interp_slice_union_vectorized(
    M0: sparse.COO,
    M1: sparse.COO,
    w: float,
    idx_dtype: np.dtype,
    val_dtype: np.dtype,
) -> sparse.COO:
    """
    Vectorized interpolation between two 2D COO slices.
    Computes union of nonzero coordinates and interpolates values on that union.

    Assumes M0, M1 have same shape.
    """
    ncols = int(M0.shape[1])

    i0 = M0.coords[0].astype(np.int64, copy=False)
    j0 = M0.coords[1].astype(np.int64, copy=False)
    v0 = M0.data.astype(np.float64, copy=False)

    i1 = M1.coords[0].astype(np.int64, copy=False)
    j1 = M1.coords[1].astype(np.int64, copy=False)
    v1 = M1.data.astype(np.float64, copy=False)

    # Encode (i,j) -> key so union/intersection are cheap
    k0 = i0 * ncols + j0
    k1 = i1 * ncols + j1

    # Sort once
    o0 = np.argsort(k0)
    o1 = np.argsort(k1)
    k0s, v0s = k0[o0], v0[o0]
    k1s, v1s = k1[o1], v1[o1]

    # Union of keys
    ku = np.union1d(k0s, k1s)

    # Align v0 and v1 on union positions (missing -> 0)
    pos0 = np.searchsorted(k0s, ku)
    pos1 = np.searchsorted(k1s, ku)

    v0u = np.zeros(ku.shape[0], dtype=np.float64)
    v1u = np.zeros(ku.shape[0], dtype=np.float64)

    m0 = (pos0 < k0s.size) & (k0s[pos0] == ku)
    m1 = (pos1 < k1s.size) & (k1s[pos1] == ku)

    v0u[m0] = v0s[pos0[m0]]
    v1u[m1] = v1s[pos1[m1]]

    # Interpolate
    vv = (1.0 - w) * v0u + w * v1u

    # Prune exact zeros to keep sparsity
    nz = vv != 0.0
    ku = ku[nz]
    vv = vv[nz].astype(val_dtype, copy=False)

    # Decode keys back to (i,j)
    ii = (ku // ncols).astype(idx_dtype, copy=False)
    jj = (ku % ncols).astype(idx_dtype, copy=False)

    return sparse.COO(coords=[ii, jj], data=vv, shape=M0.shape)


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


def _parse_amount_source(value: object) -> str:
    """Parse temporal_amount_source; defaults to 'port'."""
    if value is None:
        return "port"
    s = str(value).strip().lower()
    if s == "":
        return "port"
    if s in {"port", "matrix"}:
        return s
    # Be strict (recommended) or soft fallback. I recommend strict for debugging:
    raise ValueError(f"Unknown temporal_amount_source: {s}")
