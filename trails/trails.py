# trails.py

from dataclasses import dataclass
import time
from typing import Any, Dict, List, Optional, Callable
import importlib
import importlib.util
from collections import defaultdict, deque

import numpy as np
import sparse
import xarray as xr

from tqdm import tqdm

from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package,
)

from .temporal_distributions import TemporalDistribution, TemporalExchange

import logging

logger = logging.getLogger(__name__)


def _log_every(n: int, i: int) -> bool:
    return n > 0 and (i % n == 0)


@dataclass
class _BioAccumulationContext:
    base_year: int
    scenario_year: int
    t: int
    row_ptr: np.ndarray
    flow_sorted: np.ndarray
    data_sorted: np.ndarray
    act_coords: np.ndarray
    flow_coords: np.ndarray
    data: np.ndarray
    n_acts: int
    value_dtype: np.dtype
    scenario_index_get: Callable[[str], int | None]
    map_year_to_scenario: Callable[[int], int]
    tpl_label: str | None
    bio_td_get: Callable[[tuple[str, int, int]], TemporalExchange | None] | None
    td_expanded_cache: dict
    pulse_cache: dict
    year_map_cache: dict[int, int]
    t_eff_cache: dict[int, int | None]
    B_row_cache_local: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]


class Trails:
    """
    Wrapper around 3D sparse matrices A and B loaded from a Frictionless
    data package, with optional temporal interpolation.

    Dimensions:
        A: (scenario, activity, product)
        B: (scenario, activity, biosphere_flow)
    """

    def __init__(
        self,
        package: Any,
        interpolate_annual: bool = True,
        value_dtype: np.dtype = np.float32,
        index_dtype: np.dtype = np.int32,
        debug: bool = False,
    ) -> None:
        """Initialize Trails by loading matrices, indices, and metadata.

        :param package: Frictionless data package to load.
        :type package: Package
        :param interpolate_annual: Whether to interpolate matrices to annual resolution.
        :type interpolate_annual: bool
        :param value_dtype: Data type for matrix values.
        :type value_dtype: numpy.dtype
        :param index_dtype: Data type for matrix indices.
        :type index_dtype: numpy.dtype
        :param debug: Whether to emit debug logging.
        :type debug: bool
        """
        self.package = package
        self.value_dtype = value_dtype
        self.index_dtype = index_dtype
        self.debug = debug

        self.scenario_labels: List[str] = []
        self.scenario_index: Dict[str, int] = {}

        self.A: Optional[sparse.COO] = None
        self.B: Optional[sparse.COO] = None
        self.inventory: Optional[xr.DataArray] = None
        self.characterized_inventory: Optional[xr.DataArray] = None
        self.static_score: Optional[float] = None
        self._inventory_years: Optional[np.ndarray] = None
        self._inventory_year_index: dict[int, int] = {}
        self._inventory_coords: Optional[list[list[np.ndarray]]] = None
        self._inventory_data: Optional[list[np.ndarray]] = None
        self.provenance: Optional[dict] = None

        print("Loading matrices from data package          [1/3]")
        (
            self.A,
            self.B,
            self.scenario_labels,
            self.scenario_index,
            self.temporal_technosphere_exchanges,
            self.temporal_biosphere_exchanges,
        ) = load_matrices_from_package(
            package=self.package,
            value_dtype=self.value_dtype,
            index_dtype=self.index_dtype,
            debug=debug,
        )

        if debug:
            logger.info(
                "Trails init: scenarios=%d year_range=[%s..%s] A=%s nnz=%s B=%s nnz=%s temporal_exchanges=%d",
                len(self.scenario_labels),
                getattr(self, "min_year", None),
                getattr(self, "max_year", None),
                None if self.A is None else self.A.shape,
                None if self.A is None else int(self.A.nnz),
                None if self.B is None else self.B.shape,
                None if self.B is None else int(self.B.nnz),
                len(getattr(self, "temporal_exchanges", {})),
            )

        self.template_labels = list(self.scenario_labels)
        self.template_years_int = np.array(
            [int(lbl) for lbl in self.template_labels], dtype=int
        )

        self.years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
        self.min_year = int(self.years_int.min())
        self.max_year = int(self.years_int.max())

        # Load indices/metadata
        print("Loading indices from data package           [2/3]")
        (
            self.activity_indices,
            self.biosphere_indices,
        ) = load_indices_from_package(self.package)

        # Optional temporal interpolation to annual resolution
        if interpolate_annual and self.scenario_labels:
            print("Interpolating matrices to annual resolution [3/3]")
            (
                self.A,
                self.B,
                self.scenario_labels,
                self.scenario_index,
            ) = interpolate_to_annual(
                self.A,
                self.B,
                self.scenario_labels,
                value_dtype=self.value_dtype,
                debug=debug,
            )

            self.years_int = np.array(
                [int(lbl) for lbl in self.scenario_labels], dtype=int
            )
            self.min_year = int(self.years_int.min())
            self.max_year = int(self.years_int.max())

        self.scores: Optional[xr.DataArray] = (
            None  # dims: (activity, year) or (activity, year, root activity) or (+method)
        )
        self._score_years: Optional[np.ndarray] = None
        self._score_year_index: dict[int, int] = {}

    def reset_scores(
        self,
        *,
        attribute_to_roots: bool = False,
    ) -> None:
        # Always initialize score years, independent of inventory
        min_offset, max_offset = self._inventory_offset_bounds()
        years = np.arange(
            int(self.min_year + min_offset),
            int(self.max_year + max_offset) + 1,
            dtype=int,
        )

        self._score_years = years
        self._score_year_index = {int(y): int(i) for i, y in enumerate(years)}
        self._scores_has_root = bool(attribute_to_roots)

        self._score_chunk_act = []
        self._score_chunk_year = []
        self._score_chunk_root = []
        self._score_chunk_value = []

        # Bulk score builder (vectorized appends)
        self._score_bulk_act = []
        self._score_bulk_year = []
        self._score_bulk_root = []
        self._score_bulk_value = []

        self.scores = None

    def _clamp_year_to_inventory(self, year: int) -> int:
        years = self._inventory_years
        if years is None or years.size == 0:
            raise RuntimeError(
                "Inventory years not initialized; call reset_inventory() first."
            )
        y = int(year)
        if y < int(years[0]):
            return int(years[0])
        if y > int(years[-1]):
            return int(years[-1])
        return y

    def _clamp_year_to_scores(self, year: int) -> int:
        years = self._score_years
        if years is None or years.size == 0:
            raise RuntimeError(
                "Score years not initialized; call reset_scores() or reset_inventory() first."
            )
        y = int(year)
        if y < int(years[0]):
            return int(years[0])
        if y > int(years[-1]):
            return int(years[-1])
        return y

    def _append_scores_from_yearidx_map(
        self,
        act_idx: int,
        yearidx_to_value: dict[int, float],
        *,
        root_activity: int | None = None,
    ) -> None:
        """Append many score entries at once from a map {year_idx: value}."""
        if not yearidx_to_value:
            return

        # Filter zeros fast (and avoid building arrays if empty)
        items = [(yidx, v) for yidx, v in yearidx_to_value.items() if v != 0.0]
        if not items:
            return

        n = len(items)
        self._score_chunk_act.extend([int(act_idx)] * n)
        self._score_chunk_year.extend([int(yidx) for yidx, _ in items])
        self._score_chunk_value.extend([float(v) for _, v in items])

        if getattr(self, "_scores_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            self._score_chunk_root.extend([int(root_activity)] * n)

    def reset_inventory(
        self,
        *,
        attribute_to_roots: bool = False,
        reset_scores: bool = True,
    ) -> None:
        min_offset, max_offset = self._inventory_offset_bounds()
        years = np.arange(
            int(self.min_year + min_offset),
            int(self.max_year + max_offset) + 1,
            dtype=int,
        )

        self._inventory_years = years
        self._inventory_year_index = {int(y): int(i) for i, y in enumerate(years)}
        self._inventory_has_root = bool(attribute_to_roots)

        # Initialize inventory chunk builders (block-based appends)
        self._inv_chunk_act = []
        self._inv_chunk_year = []
        self._inv_chunk_root = []
        self._inv_chunk_flows = []
        self._inv_chunk_values = []
        self._inv_chunk_len = []

        # Initialize inventory bulk builders (vectorized appends)
        self._inv_bulk_act = []
        self._inv_bulk_year = []
        self._inv_bulk_flow = []
        self._inv_bulk_value = []
        self._inv_bulk_root = []
        # Reset outputs
        self.inventory = None
        self.characterized_inventory = None
        self.static_score = None
        self.provenance = None

        if reset_scores:
            # Reset scores builder/output (direct LCIA scoring)
            self._score_years = self._inventory_years
            self._score_year_index = self._inventory_year_index
            self._scores_has_root = bool(attribute_to_roots)

            self._score_chunk_act = []
            self._score_chunk_year = []
            self._score_chunk_root = []
            self._score_chunk_value = []

            self._score_bulk_act = []
            self._score_bulk_year = []
            self._score_bulk_root = []
            self._score_bulk_value = []

            self.scores = None

    def _append_score_entry(
        self,
        act_idx: int,
        year: int,
        value: float,
        *,
        root_activity: int | None = None,
    ) -> None:
        """Append a scalar score entry to the chunked score builder."""
        if not hasattr(self, "_score_chunk_value"):
            raise RuntimeError(
                "Score builders not initialized. Call reset_scores() or reset_inventory() first."
            )

        y = self._clamp_year_to_scores(int(year))
        year_idx = self._score_year_index.get(y)
        if year_idx is None:
            return

        v = float(value)
        if v == 0.0:
            return

        self._score_chunk_act.append(int(act_idx))
        self._score_chunk_year.append(int(year_idx))
        self._score_chunk_value.append(v)

        if getattr(self, "_scores_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            self._score_chunk_root.append(int(root_activity))

    def finalize_scores(self) -> xr.DataArray:
        """Finalize and store sparse scalar scores as an xarray."""
        years = self._score_years
        if years is None:
            raise RuntimeError("Scores years not initialized. Call reset_inventory().")
        if self.A is None:
            raise RuntimeError("A is None")

        n_activities = int(self.A.shape[1])
        has_root = bool(self._scores_has_root)

        has_any = bool(self._score_chunk_act) or bool(
            getattr(self, "_score_bulk_act", [])
        )
        if not has_any:
            if has_root:
                arr = sparse.COO.zeros(
                    (n_activities, len(years), n_activities), dtype=self.value_dtype
                )
                self.scores = xr.DataArray(
                    arr,
                    dims=("activity", "year", "root activity"),
                    coords={
                        "activity": np.arange(n_activities, dtype=int),
                        "year": years,
                        "root activity": np.arange(n_activities, dtype=int),
                    },
                )
            else:
                arr = sparse.COO.zeros(
                    (n_activities, len(years)), dtype=self.value_dtype
                )
                self.scores = xr.DataArray(
                    arr,
                    dims=("activity", "year"),
                    coords={
                        "activity": np.arange(n_activities, dtype=int),
                        "year": years,
                    },
                )
            return self.scores

        coords_parts = []
        data_parts = []
        root_parts = []

        # Chunk parts
        if self._score_chunk_act:
            act_c = np.asarray(self._score_chunk_act, dtype=np.int64)
            yr_c = np.asarray(self._score_chunk_year, dtype=np.int64)
            data_c = np.asarray(self._score_chunk_value, dtype=self.value_dtype)
            coords_parts.append((act_c, yr_c))
            data_parts.append(data_c)
            if has_root:
                root_c = np.asarray(self._score_chunk_root, dtype=np.int64)
                root_parts.append(root_c)

        # Bulk parts
        if getattr(self, "_score_bulk_act", []):
            act_b = np.concatenate(self._score_bulk_act).astype(np.int64, copy=False)
            yr_b = np.concatenate(self._score_bulk_year).astype(np.int64, copy=False)
            data_b = np.concatenate(self._score_bulk_value).astype(
                self.value_dtype, copy=False
            )

            coords_parts.append((act_b, yr_b))
            data_parts.append(data_b)
            if has_root:
                root_b = np.concatenate(self._score_bulk_root).astype(
                    np.int64, copy=False
                )
                root_parts.append(root_b)

        act = np.concatenate([p[0] for p in coords_parts])
        yr = np.concatenate([p[1] for p in coords_parts])
        data = np.concatenate(data_parts).astype(self.value_dtype, copy=False)

        if has_root:
            root = np.concatenate(root_parts).astype(np.int64, copy=False)
            coords = np.vstack([act, yr, root])
            arr = sparse.COO(
                coords, data, shape=(n_activities, len(years), n_activities)
            )

            self.scores = xr.DataArray(
                arr,
                dims=("activity", "year", "root activity"),
                coords={
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                    "root activity": np.arange(n_activities, dtype=int),
                },
            )
        else:
            coords = np.vstack([act, yr])
            arr = sparse.COO(coords, data, shape=(n_activities, len(years)))
            self.scores = xr.DataArray(
                arr,
                dims=("activity", "year"),
                coords={"activity": np.arange(n_activities, dtype=int), "year": years},
            )

        return self.scores

    def _inventory_offset_bounds(self) -> tuple[int, int]:
        """Return min/max year offsets implied by temporal exchange metadata."""
        min_offset = 0
        max_offset = 0
        for exchanges in (
            self.temporal_technosphere_exchanges,
            self.temporal_biosphere_exchanges,
        ):
            if not exchanges:
                continue
            for tex in exchanges.values():
                offset_min = getattr(tex, "offset_min", None)
                offset_max = getattr(tex, "offset_max", None)
                if offset_min is not None:
                    min_offset = min(min_offset, int(offset_min))
                if offset_max is not None:
                    max_offset = max(max_offset, int(offset_max))
        return min_offset, max_offset

    def _append_inventory_entries(
        self,
        act_idx: int,
        year: int,
        flows: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: int | None = None,
    ) -> None:
        """
        Append one (activity, year, [root]) block of biosphere contributions to the
        chunked inventory builder.

        Contract:
          - `year` is a *calendar year*; we map to `year_idx` using `_inventory_year_index`.
          - `flows` are flow indices; `values` are amounts for those flows.
          - Stores compact chunks; `finalize_inventory()` expands coords via np.repeat.

        Performance notes:
          - Avoids unnecessary copies.
          - Filters zeros quickly.
          - Converts dtypes only at the end, and with copy=False when possible.
        """
        # Chunk builders must exist (single source of truth)
        if not hasattr(self, "_inv_chunk_len"):
            raise RuntimeError(
                "Inventory chunk builders not initialized. Call reset_inventory() first."
            )

        y = self._clamp_year_to_inventory(int(year))
        year_idx = self._inventory_year_index.get(y)
        if year_idx is None:
            return  # should not happen unless year index map is inconsistent

        # Fast-path: accept ndarray-like without forcing copies
        flows_arr = np.asarray(flows)
        vals_arr = np.asarray(values)

        # Basic validation
        if flows_arr.size == 0 or vals_arr.size == 0:
            return
        if flows_arr.shape[0] != vals_arr.shape[0]:
            raise ValueError(
                f"flows and values must have same length. Got {flows_arr.shape[0]} and {vals_arr.shape[0]}."
            )

        # Filter exact zeros (cheap and avoids storing empty chunks)
        # (If you want tolerance-based filtering, keep that at the caller level.)
        mask = vals_arr != 0.0
        if not np.any(mask):
            return
        if not np.all(mask):
            flows_arr = flows_arr[mask]
            vals_arr = vals_arr[mask]
            if flows_arr.size == 0:
                return

        # Enforce dtypes late, avoid copies when possible
        flows_i64 = flows_arr.astype(np.int64, copy=False)

        # Values: ensure numeric dtype + match configured value_dtype
        # (This keeps memory small and makes concatenation cheaper.)
        if vals_arr.dtype != self.value_dtype:
            vals_out = vals_arr.astype(self.value_dtype, copy=False)
        else:
            vals_out = vals_arr

        n = int(flows_i64.size)
        if n == 0:
            return

        # Append one chunk
        self._inv_chunk_act.append(int(act_idx))
        self._inv_chunk_year.append(int(year_idx))

        if getattr(self, "_inventory_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            self._inv_chunk_root.append(int(root_activity))

        self._inv_chunk_flows.append(flows_i64)
        self._inv_chunk_values.append(vals_out)
        self._inv_chunk_len.append(n)

    def _append_scores_bulk(
        self,
        act_idx: np.ndarray,
        year_idx: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: np.ndarray | None = None,
    ) -> None:
        """Append aligned arrays of (act, year_idx, value[, root]) into score bulk builder."""
        if not hasattr(self, "_score_bulk_value"):
            raise RuntimeError(
                "Score bulk builders not initialized. Call reset_scores() or reset_inventory() first."
            )

        a = np.asarray(act_idx)
        y = np.asarray(year_idx)
        v = np.asarray(values)

        if a.size == 0 or y.size == 0 or v.size == 0:
            return
        if not (a.size == y.size == v.size):
            raise ValueError("act_idx, year_idx, values must have same length")

        # Filter exact zeros early
        m = v != 0.0
        if not np.any(m):
            return

        a = a[m].astype(np.int64, copy=False)
        y = y[m].astype(np.int64, copy=False)

        if v.dtype != self.value_dtype:
            v = v[m].astype(self.value_dtype, copy=False)
        else:
            v = v[m]

        self._score_bulk_act.append(a)
        self._score_bulk_year.append(y)
        self._score_bulk_value.append(v)

        if getattr(self, "_scores_has_root", False):
            if root_activity is None:
                raise ValueError(
                    "root_activity must be provided when scores have root dimension"
                )
            r = np.asarray(root_activity)
            if r.shape != a.shape:
                raise ValueError("root_activity must match act_idx shape")
            self._score_bulk_root.append(r[m].astype(np.int64, copy=False))

    def _append_inventory_entries_bulk(
        self,
        act_idx: np.ndarray,
        year: int | np.ndarray,
        flows: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: int | np.ndarray | None = None,
    ) -> None:
        """Append inventory entries for aligned arrays of (act, flow, year, value)."""
        if not hasattr(self, "_inv_bulk_flow"):
            raise RuntimeError(
                "Inventory bulk builders not initialized. Call reset_inventory() first."
            )

        acts_arr = np.asarray(act_idx)
        flows_arr = np.asarray(flows)
        vals_arr = np.asarray(values)

        if flows_arr.size == 0 or vals_arr.size == 0 or acts_arr.size == 0:
            return
        if (
            flows_arr.shape[0] != vals_arr.shape[0]
            or flows_arr.shape[0] != acts_arr.shape[0]
        ):
            raise ValueError(
                "act, flow, and value arrays must have the same length for bulk append."
            )

        if isinstance(year, np.ndarray):
            years_arr = np.asarray(year)
            if years_arr.shape[0] != flows_arr.shape[0]:
                raise ValueError(
                    "year array must match act/flow/value length for bulk append."
                )
            years_arr = np.asarray(year)
            y0 = int(self._inventory_years[0])
            y1 = int(self._inventory_years[-1])
            years_clamped = np.clip(years_arr.astype(np.int64, copy=False), y0, y1)

            year_idx = np.fromiter(
                (self._inventory_year_index[int(y)] for y in years_clamped),
                dtype=np.int64,
                count=years_clamped.size,
            )
        else:
            y = self._clamp_year_to_inventory(int(year))
            year_idx_val = self._inventory_year_index.get(y)
            if year_idx_val is None:
                return

            year_idx = np.full(flows_arr.shape[0], int(year_idx_val), dtype=np.int64)

        mask = vals_arr != 0.0
        if not np.any(mask):
            return

        acts_out = acts_arr[mask].astype(np.int64, copy=False)
        flows_out = flows_arr[mask].astype(np.int64, copy=False)
        years_out = year_idx[mask]

        if vals_arr.dtype != self.value_dtype:
            vals_out = vals_arr[mask].astype(self.value_dtype, copy=False)
        else:
            vals_out = vals_arr[mask]

        self._inv_bulk_act.append(acts_out)
        self._inv_bulk_year.append(years_out)
        self._inv_bulk_flow.append(flows_out)
        self._inv_bulk_value.append(vals_out)

        if getattr(self, "_inventory_has_root", False):
            if root_activity is None:
                root_arr = acts_out
            else:
                root_arr = np.asarray(root_activity)
                if root_arr.shape == ():
                    root_arr = np.full_like(acts_out, int(root_arr))
                elif root_arr.shape[0] != acts_out.shape[0]:
                    raise ValueError(
                        "root_activity array must match act/flow/value length for bulk append."
                    )
                else:
                    root_arr = root_arr[mask]
            self._inv_bulk_root.append(np.asarray(root_arr, dtype=np.int64))

    def finalize_inventory(self) -> xr.DataArray:
        """Finalize and store sparse inventory as an xarray."""
        if self.A is None or self.B is None:
            raise ValueError("Cannot finalize inventory: A or B is None.")

        years = self._inventory_years
        if years is None:
            raise RuntimeError(
                "Inventory years not initialized. Call reset_inventory()."
            )

        if not hasattr(self, "_inv_chunk_len"):
            raise RuntimeError(
                "Inventory chunk builders not initialized. Call reset_inventory()."
            )

        n_activities = int(self.A.shape[1])
        n_flows = int(self.B.shape[2])
        has_root = bool(getattr(self, "_inventory_has_root", False))

        coords_parts: list[np.ndarray] = []
        data_parts: list[np.ndarray] = []
        root_parts: list[np.ndarray] = []

        if self._inv_chunk_len:
            flows_chunk = np.concatenate(self._inv_chunk_flows).astype(
                np.int64, copy=False
            )
            data_chunk = np.concatenate(self._inv_chunk_values)

            lens = np.asarray(self._inv_chunk_len, dtype=np.int64)

            act_chunk = np.repeat(np.asarray(self._inv_chunk_act, dtype=np.int64), lens)
            year_chunk = np.repeat(
                np.asarray(self._inv_chunk_year, dtype=np.int64), lens
            )

            coords_parts.append(np.vstack([act_chunk, flows_chunk, year_chunk]))
            data_parts.append(data_chunk)

            if has_root:
                root_chunk = np.repeat(
                    np.asarray(self._inv_chunk_root, dtype=np.int64), lens
                )
                root_parts.append(root_chunk)

        if self._inv_bulk_flow:
            act_bulk = np.concatenate(self._inv_bulk_act).astype(np.int64, copy=False)
            flow_bulk = np.concatenate(self._inv_bulk_flow).astype(np.int64, copy=False)
            year_bulk = np.concatenate(self._inv_bulk_year).astype(np.int64, copy=False)
            data_bulk = np.concatenate(self._inv_bulk_value)
            coords_parts.append(np.vstack([act_bulk, flow_bulk, year_bulk]))
            data_parts.append(data_bulk)

            if has_root:
                root_bulk = np.concatenate(self._inv_bulk_root).astype(
                    np.int64, copy=False
                )
                root_parts.append(root_bulk)

        if data_parts:
            coords_base = np.hstack(coords_parts)
            data = np.concatenate(data_parts)
            if has_root:
                root_all = np.concatenate(root_parts)
                coords = np.vstack([coords_base, root_all])
                inv = sparse.COO(
                    coords,
                    data,
                    shape=(n_activities, n_flows, len(years), n_activities),
                )
            else:
                inv = sparse.COO(
                    coords_base,
                    data,
                    shape=(n_activities, n_flows, len(years)),
                )
        else:
            if has_root:
                inv = sparse.COO.zeros(
                    (n_activities, n_flows, len(years), n_activities),
                    dtype=self.value_dtype,
                )
            else:
                inv = sparse.COO.zeros(
                    (n_activities, n_flows, len(years)),
                    dtype=self.value_dtype,
                )

        dims = ("activity", "flow", "year")
        coords_xr = {
            "activity": np.arange(n_activities, dtype=int),
            "flow": np.arange(n_flows, dtype=int),
            "year": years,
        }
        if has_root:
            dims = ("activity", "flow", "year", "root activity")
            coords_xr["root activity"] = np.arange(n_activities, dtype=int)

        self.inventory = xr.DataArray(inv, dims=dims, coords=coords_xr)
        return self.inventory

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------
    def _interpolate_temporal_exchange(
        self,
        year: int,
        act_idx: int,
        other_idx: int,
        exchanges: Dict[tuple, TemporalExchange],
    ) -> Optional[TemporalExchange]:
        """Interpolate temporal exchange metadata for a given year and indices.

        :param year: Calendar year to resolve.
        :type year: int
        :param act_idx: Activity index.
        :type act_idx: int
        :param other_idx: Product or flow index.
        :type other_idx: int
        :param exchanges: Mapping of exchange metadata keyed by scenario label.
        :type exchanges: dict[tuple, TemporalExchange]
        :returns: Interpolated exchange metadata, if available.
        :rtype: TemporalExchange | None
        """
        if not exchanges:
            return None

        label = str(year)
        direct = exchanges.get((label, act_idx, other_idx))
        if direct is not None:
            return direct

        entries = [
            (int(lbl), tex)
            for (lbl, a_idx, o_idx), tex in exchanges.items()
            if a_idx == act_idx and o_idx == other_idx
        ]
        if not entries:
            return None

        entries.sort(key=lambda pair: pair[0])
        years = [y for y, _ in entries]

        if year <= years[0]:
            return entries[0][1]
        if year >= years[-1]:
            return entries[-1][1]

        for (y0, tex0), (y1, tex1) in zip(entries, entries[1:]):
            if y0 <= year <= y1:
                if y1 == y0:
                    return tex0
                if (
                    tex0.distribution != tex1.distribution
                    or tex0.offset_min != tex1.offset_min
                    or tex0.offset_max != tex1.offset_max
                    or getattr(tex0, "amount_source", "port")
                    != getattr(tex1, "amount_source", "port")
                ):
                    return tex0 if (year - y0) <= (y1 - year) else tex1

                w = (year - y0) / (y1 - y0)

                def interp_optional(v0: float | None, v1: float | None) -> float | None:
                    """Interpolate optional numeric values with nearest fallback.

                    :param v0: First value to interpolate.
                    :type v0: float | None
                    :param v1: Second value to interpolate.
                    :type v1: float | None
                    :returns: Interpolated or nearest value.
                    :rtype: float | None
                    """
                    if v0 is None or v1 is None:
                        return v0 if (year - y0) <= (y1 - year) else v1
                    return float(v0) + (float(v1) - float(v0)) * w

                return TemporalExchange(
                    distribution=tex0.distribution,
                    loc=interp_optional(tex0.loc, tex1.loc),
                    scale=interp_optional(tex0.scale, tex1.scale),
                    offset_min=tex0.offset_min,
                    offset_max=tex0.offset_max,
                    amount_source=getattr(tex0, "amount_source", "port"),
                )

        return None

    def _map_year_to_scenario_year(self, year: int) -> int:
        """
        Map an arbitrary year to the closest available scenario year (where A/B exist),
        clipped to [min_year, max_year].
        """
        y = max(self.min_year, min(self.max_year, int(year)))

        # If we have a full annual grid, this is effectively identity after clipping
        if len(self.years_int) == (self.max_year - self.min_year + 1):
            return y

        # Otherwise: snap to nearest scenario year
        idx = int(np.abs(self.years_int - y).argmin())
        return int(self.years_int[idx])

    def _map_year_to_template_year(self, year: int) -> int:
        """
        Map an arbitrary year to the nearest original year that has temporal metadata.
        """
        y = int(year)
        idx = int(np.abs(self.template_years_int - y).argmin())
        return int(self.template_years_int[idx])

    def _get_scenario_context(self, year: int) -> tuple[int, str, int] | None:
        """Return the scenario tuple for a given year if available.

        :param year: Calendar year to map.
        :type year: int
        :returns: ``(scenario_year, scenario_label, scenario_index)`` or ``None``.
        :rtype: tuple[int, str, int] | None
        """
        scenario_year = self._map_year_to_scenario_year(year)
        scenario_label = str(scenario_year)
        if scenario_label not in self.scenario_index:
            return None
        t = self.scenario_index[scenario_label]
        return scenario_year, scenario_label, t

    @staticmethod
    def _add_demand_entry(
        demand: dict[int, dict[int, float]],
        target_year: int,
        product_index: int,
        exchange_amount: float,
    ) -> None:
        """Accumulate a demand amount for a given year and product index.

        :param demand: Nested demand mapping to update.
        :type demand: dict[int, dict[int, float]]
        :param target_year: Scenario year to update.
        :type target_year: int
        :param product_index: Product index to add demand for.
        :type product_index: int
        :param exchange_amount: Amount to add.
        :type exchange_amount: float
        """
        demand.setdefault(target_year, {})
        demand[target_year][product_index] = (
            demand[target_year].get(product_index, 0.0) + exchange_amount
        )

    @staticmethod
    def _child_amount(parent_amount: float, exchange_value: float) -> float:
        """Convert an exchange value into a child demand amount.

        :param parent_amount: Upstream demand amount.
        :type parent_amount: float
        :param exchange_value: Exchange coefficient from the A matrix.
        :type exchange_value: float
        :returns: Child demand amount with sign handled.
        :rtype: float
        """
        if exchange_value < 0.0:
            return parent_amount * (-exchange_value)
        return parent_amount * exchange_value

    def _apply_temporal_distribution_to_demand(
        self,
        *,
        year: int,
        product_index: int,
        child_amount: float,
        tex: TemporalExchange,
        demand: dict[int, dict[int, float]],
        debug: bool,
    ) -> None:
        """Apply temporal distribution metadata to a demand entry.

        :param year: Base year for the demand.
        :type year: int
        :param scenario_year: Scenario year used for A/B lookup.
        :type scenario_year: int
        :param product_index: Product index receiving the demand.
        :type product_index: int
        :param child_amount: Amount to distribute across offsets.
        :type child_amount: float
        :param tex: Temporal exchange metadata to apply.
        :type tex: TemporalExchange
        :param demand: Demand mapping to update.
        :type demand: dict[int, dict[int, float]]
        :param debug: Whether to emit debug logging.
        :type debug: bool
        """
        td = TemporalDistribution(tex)

        offsets_and_weights = list(td.iter_offsets_and_weights(debug=debug))
        if not offsets_and_weights:
            if debug:
                logger.warning(
                    "expand_temporal_exchanges: TD produced no offsets/weights for (year=%d prod=%d) -> dropping exchange",
                    year,
                    product_index,
                )
            return

        if debug:
            logger.debug(
                "expand_temporal_exchanges: TD offsets=%s (sum_w=%g)",
                [p[0] for p in offsets_and_weights],
                float(sum(p[1] for p in offsets_and_weights)),
            )

        for offset, weight in offsets_and_weights:
            raw_year = year + offset
            self._add_demand_entry(
                demand,
                int(raw_year),
                product_index,
                child_amount * float(weight),
            )

    def get_A_for_scenario(self, label: str) -> sparse.COO:
        """Return the 2D A matrix (activity x product) for a given scenario label."""
        t = self.scenario_index[label]
        return self.A[t, :, :]

    def get_B_for_scenario(self, label: str) -> sparse.COO:
        """Return the 2D B matrix (activity x flow) for a given scenario label."""
        t = self.scenario_index[label]
        return self.B[t, :, :]

    def get_temporal_exchange(
        self, year: int, act_idx: int, prod_idx: int
    ) -> TemporalExchange | None:
        """Return temporal exchange metadata for a given activity/product.

        :param year: Calendar year to query.
        :type year: int
        :param act_idx: Activity index.
        :type act_idx: int
        :param prod_idx: Product index.
        :type prod_idx: int
        :returns: Temporal exchange metadata if available.
        :rtype: TemporalExchange | None
        """
        return self._interpolate_temporal_exchange(
            year,
            act_idx,
            prod_idx,
            self.temporal_technosphere_exchanges,
        )

    def get_temporal_distribution(
        self, year: int, act_idx: int, prod_idx: int
    ) -> TemporalDistribution | None:
        """Return temporal distribution metadata for a given activity/product.

        :param year: Calendar year to query.
        :type year: int
        :param act_idx: Activity index.
        :type act_idx: int
        :param prod_idx: Product index.
        :type prod_idx: int
        :returns: Temporal distribution if metadata exists.
        :rtype: TemporalDistribution | None
        """
        tex = self._interpolate_temporal_exchange(
            year,
            act_idx,
            prod_idx,
            self.temporal_technosphere_exchanges,
        )
        if tex is None:
            return None
        return TemporalDistribution(tex)

    def lca(self, *args: Any, **kwargs: Any) -> None:
        from .lca import lca as lca_fn

        if "debug" in kwargs:
            self.debug = bool(kwargs.pop("debug"))

        return lca_fn(self, *args, **kwargs)

    def static_lca(
        self,
        year: int,
        act_idx: int,
        methods: list[str],
        amount: float = 1.0,
        debug: bool = False,
    ) -> None:
        """Run static LCA using this Trails instance."""
        from .lca import lca_static_simple

        return lca_static_simple(
            trails=self,
            year=int(year),
            fu_act_idx=int(act_idx),
            methods=methods,
            amount=float(amount),
            debug=debug,
        )

    def expand_temporal_exchanges(
        self,
        year: int,
        act_idx: int,
        amount: float = 1.0,
        *,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> dict[int, dict[int, float]]:
        """Expand activity-year demand into temporally distributed demands."""
        demand: dict[int, dict[int, float]] = {}
        # --- summary counters for low-noise debugging ---
        n_exchanges = 0
        n_skipped_prod = 0
        n_no_td = 0
        n_td_ported = 0
        n_td_matrix = 0
        min_raw_year = None
        max_raw_year = None

        context = self._get_scenario_context(year)
        if context is None:
            return demand
        scenario_year, scenario_label, t = context

        if debug:
            logger.info(
                "expand_tech: base_year=%d scen_year=%d t=%d act=%d amount=%g",
                int(year),
                int(scenario_year),
                int(t),
                int(act_idx),
                float(amount),
            )

        A_row = self.A[t, act_idx, :]
        if A_row.nnz == 0:
            return demand

        product_indices = A_row.coords[0]
        values = A_row.data

        for product_index, exchange_value in zip(product_indices, values):
            n_exchanges += 1
            product_index = int(product_index)
            exchange_value = float(exchange_value)

            if exchange_value == 0.0:
                continue

            # Skip canonical production exchange (A[act, act] = 1)
            if product_index == act_idx:
                n_skipped_prod += 1
                continue

            # Fetch TD metadata (template-year lookup; stable across interpolation)
            tex = self._get_tech_temporal_exchange(year, act_idx, product_index)

            # ------------------------------------------------------------------
            # No temporal distribution (or disabled): status quo
            # ------------------------------------------------------------------
            if (tex is None) or (not use_temporal_distributions):
                n_no_td += 1
                child_amount = self._child_amount(amount, exchange_value)
                if child_amount != 0.0:
                    self._add_demand_entry(
                        demand, int(year), product_index, float(child_amount)
                    )
                continue

            amount_source = getattr(tex, "amount_source", "port")

            # ------------------------------------------------------------------
            # TD + matrix-sourced magnitude: read A at each pulse year
            # ------------------------------------------------------------------
            if amount_source == "matrix":
                n_td_matrix += 1
                self._apply_temporal_distribution_matrix_sourced_to_demand(
                    year=year,
                    act_idx=act_idx,
                    product_index=product_index,
                    parent_amount=amount,
                    tex=tex,
                    demand=demand,
                    debug=debug,
                )
                continue

            # ------------------------------------------------------------------
            # TD + ported magnitude (default): distribute anchor-year child amount
            # ------------------------------------------------------------------
            n_td_ported += 1
            child_amount = self._child_amount(amount, exchange_value)
            if child_amount == 0.0:
                continue

            self._apply_temporal_distribution_to_demand(
                year=year,
                product_index=product_index,
                child_amount=child_amount,
                tex=tex,
                demand=demand,
                debug=debug,
            )

        if demand:
            years_out = list(demand.keys())
            min_raw_year = int(min(years_out))
            max_raw_year = int(max(years_out))

        if debug:
            out_years = int(len(demand))
            out_edges = int(sum(len(v) for v in demand.values()))
            logger.info(
                "expand_tech_done: base_year=%d scen_year=%d act=%d amount=%g "
                "exchanges=%d skipped_prod=%d no_td=%d td_ported=%d td_matrix=%d "
                "years=[%s..%s] out_years=%d out_edges=%d",
                int(year),
                int(scenario_year),
                int(act_idx),
                float(amount),
                int(n_exchanges),
                int(n_skipped_prod),
                int(n_no_td),
                int(n_td_ported),
                int(n_td_matrix),
                str(min_raw_year) if min_raw_year is not None else "NA",
                str(max_raw_year) if max_raw_year is not None else "NA",
                out_years,
                out_edges,
            )

        return demand

    def _get_tech_temporal_exchange(
        self, year: int, act_idx: int, prod_idx: int
    ) -> Optional[TemporalExchange]:
        """
        For technosphere TD metadata, do NOT interpolate across years.
        Instead, map to the nearest template year and do a direct lookup.
        This prevents TD metadata from 'bleeding' into years where it wasn't specified
        and makes TD availability stable across interpolated scenario years.
        """
        if not self.temporal_technosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        return self.temporal_technosphere_exchanges.get(
            (str(y_tpl), int(act_idx), int(prod_idx))
        )

    def _get_bio_temporal_exchange(
        self, year: int, act_idx: int, flow_idx: int
    ) -> Optional[TemporalExchange]:
        """
        For biosphere TD metadata, do NOT interpolate across years.
        Instead, map to the nearest template year and do a direct lookup.
        This prevents TD metadata from 'bleeding' into years where it wasn't specified.
        """
        if not self.temporal_biosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        return self.temporal_biosphere_exchanges.get(
            (str(y_tpl), int(act_idx), int(flow_idx))
        )

    def _get_biosphere_slice(
        self, base_year: int, debug: bool
    ) -> tuple[int, int, sparse.COO, int] | None:
        """Return biosphere slice metadata for a base year.

        :param base_year: Calendar year used to select the scenario slice.
        :type base_year: int
        :param debug: Whether to emit debug logging.
        :type debug: bool
        :returns: Tuple ``(scenario_year, t_index, B_slice, n_flows)`` or ``None``.
        :rtype: tuple[int, int, sparse.COO, int] | None
        """
        if self.B is None:
            if debug:
                logger.warning("accumulate_bio: B is None -> nothing to accumulate")
            return None

        context = self._get_scenario_context(base_year)
        if context is None:
            if debug:
                logger.error(
                    "accumulate_bio: scenario_label not in scenario_index (base_year=%d) -> abort",
                    int(base_year),
                )
            return None

        scenario_year, scenario_label, t = context
        B_t = self.B[t, :, :]
        B_t_nnz = int(getattr(B_t, "nnz", 0))
        if debug:
            logger.info(
                "accumulate_bio: B slice t=%d nnz=%d shape=%s",
                int(t),
                B_t_nnz,
                getattr(B_t, "shape", None),
            )
        if B_t_nnz == 0:
            if debug:
                logger.warning("accumulate_bio: B_t.nnz==0 -> nothing to accumulate")
            return None

        n_flows = int(self.B.shape[2])
        return scenario_year, t, B_t, n_flows

    def _get_B_csr_for_t(self, t: int) -> sparse.GCXS:
        """Return (and cache) a CSR slice for B[t,:,:]."""
        if not hasattr(self, "_B_csr_cache"):
            self._B_csr_cache = {}  # type: ignore[attr-defined]
        cache = self._B_csr_cache  # type: ignore[attr-defined]
        cached = cache.get(int(t))
        if cached is not None:
            return cached

        if self.B is None:
            raise RuntimeError("B is None")

        B_t = self.B[int(t), :, :]
        csr = B_t.tocsr()
        cache[int(t)] = csr
        return csr

    def _get_B_row_cache_for_t(
        self, t: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return CSR-like row cache for B[t,:,:] as (row_ptr, flow_sorted, data_sorted),
        where entries are sorted by (act, flow). This enables fast per-row lookup.
        """
        if not hasattr(self, "_B_row_cache"):
            self._B_row_cache = {}  # type: ignore[attr-defined]

        cache = self._B_row_cache  # type: ignore[attr-defined]
        cached = cache.get(int(t))
        if cached is not None:
            return cached

        if self.B is None:
            raise RuntimeError("B is None")

        csr = self._get_B_csr_for_t(int(t))
        n_acts = int(csr.shape[0])
        nnz = int(getattr(csr, "nnz", 0))

        if nnz == 0:
            row_ptr = np.zeros(n_acts + 1, dtype=np.int64)
            flow_sorted = np.zeros(0, dtype=np.int32)
            data_sorted = np.zeros(0, dtype=self.value_dtype)
            cached = (row_ptr, flow_sorted, data_sorted)
            cache[int(t)] = cached
            return cached

        row_ptr = csr.indptr.astype(np.int64, copy=False)
        flow_sorted = csr.indices.astype(np.int32, copy=False)
        data = csr.data
        if self.value_dtype == np.float32:
            data_sorted = data.astype(np.float32, copy=False)
        else:
            data_sorted = data

        cached = (row_ptr, flow_sorted, data_sorted)
        cache[int(t)] = cached
        return cached

    def _get_B_row_index_map_for_t_act(self, t: int, act: int) -> np.ndarray:
        """Return a cached flow->position map for a given (t, act) row."""
        if not hasattr(self, "_B_row_index_map"):
            self._B_row_index_map = {}  # type: ignore[attr-defined]
        cache = self._B_row_index_map  # type: ignore[attr-defined]
        key = (int(t), int(act))
        cached = cache.get(key)
        if cached is not None:
            return cached

        if self.B is None:
            raise RuntimeError("B is None")

        row_ptr, flow_sorted, _ = self._get_B_row_cache_for_t(int(t))
        start = int(row_ptr[int(act)])
        end = int(row_ptr[int(act) + 1])
        n_flows = int(self.B.shape[2])
        index_map = np.full(n_flows, -1, dtype=np.int64)
        if end > start:
            row_flows = flow_sorted[start:end].astype(np.int64, copy=False)
            index_map[row_flows] = np.arange(end - start, dtype=np.int64)
        cache[key] = index_map
        return index_map

    @staticmethod
    def _row_values_for_flows_sorted(
        row_flows_sorted: np.ndarray,
        row_vals_sorted: np.ndarray,
        query_flows_sorted: np.ndarray,
    ) -> np.ndarray:
        """
        Given a row with flows sorted ascending, return values for query_flows_sorted
        using searchsorted. Flows not present get value 0.
        """
        # positions where each query flow would be inserted
        pos = np.searchsorted(row_flows_sorted, query_flows_sorted)
        out = np.zeros(query_flows_sorted.size, dtype=np.float64)
        m = (pos < row_flows_sorted.size) & (
            row_flows_sorted[pos] == query_flows_sorted
        )
        if np.any(m):
            out[m] = row_vals_sorted[pos[m]].astype(np.float64, copy=False)
        return out

    def _build_bio_accumulation_context(
        self,
        base_year: int,
        *,
        use_temporal_distributions: bool,
        debug: bool,
    ) -> _BioAccumulationContext | None:
        """Build and return shared context for biosphere accumulation."""
        biosphere_slice = self._get_biosphere_slice(base_year, debug)
        if biosphere_slice is None:
            return None
        scenario_year, t, B_t, _ = biosphere_slice

        value_dtype = self.value_dtype
        scenario_index_get = self.scenario_index.get
        map_year_to_scenario = self._map_year_to_scenario_year

        bio_td = (
            self.temporal_biosphere_exchanges if use_temporal_distributions else None
        )
        if bio_td:
            tpl_label = str(self._map_year_to_template_year(base_year))
            bio_td_get = bio_td.get
        else:
            tpl_label = None
            bio_td_get = None  # type: ignore[assignment]

        # Global caches on self (persist across calls)
        if not hasattr(self, "_bio_td_expanded_cache"):
            self._bio_td_expanded_cache = {}  # type: ignore[attr-defined]
        td_expanded_cache = self._bio_td_expanded_cache  # type: ignore[attr-defined]

        if not hasattr(self, "_td_pulse_cache"):
            self._td_pulse_cache = {}  # type: ignore[attr-defined]
        pulse_cache = self._td_pulse_cache  # type: ignore[attr-defined]

        # Fast row-structure cache (per t)
        if not hasattr(self, "_B_row_cache_actsorted"):
            self._B_row_cache_actsorted = {}  # type: ignore[attr-defined]
        row_cache = self._B_row_cache_actsorted  # type: ignore[attr-defined]

        act_coords = B_t.coords[0].astype(np.int32, copy=False)
        flow_coords = B_t.coords[1].astype(np.int32, copy=False)
        data = (
            B_t.data.astype(np.float32, copy=False)
            if value_dtype == np.float32
            else B_t.data
        )

        cached = row_cache.get(int(t))
        if cached is None:
            n_acts = int(B_t.shape[0])
            nnz = int(getattr(B_t, "nnz", 0))
            if nnz == 0:
                row_cache[int(t)] = (
                    np.zeros(n_acts + 1, dtype=np.int64),
                    flow_coords,
                    data,
                )
                return None

            order = np.argsort(act_coords, kind="mergesort")
            act_sorted = act_coords[order]
            flow_sorted = flow_coords[order]
            data_sorted = data[order]

            row_ptr = np.zeros(n_acts + 1, dtype=np.int64)
            counts = np.bincount(act_sorted, minlength=n_acts)
            np.cumsum(counts, out=row_ptr[1:])

            cached = (row_ptr, flow_sorted, data_sorted)
            row_cache[int(t)] = cached

        row_ptr, flow_sorted, data_sorted = cached

        return _BioAccumulationContext(
            base_year=int(base_year),
            scenario_year=int(scenario_year),
            t=int(t),
            row_ptr=row_ptr,
            flow_sorted=flow_sorted,
            data_sorted=data_sorted,
            act_coords=act_coords,
            flow_coords=flow_coords,
            data=data,
            n_acts=int(B_t.shape[0]),
            value_dtype=value_dtype,
            scenario_index_get=scenario_index_get,
            map_year_to_scenario=map_year_to_scenario,
            tpl_label=tpl_label,
            bio_td_get=bio_td_get,
            td_expanded_cache=td_expanded_cache,
            pulse_cache=pulse_cache,
            year_map_cache={},
            t_eff_cache={},
            B_row_cache_local={},
        )

    def _filter_idx_with_keep(
        self, idx_full: np.ndarray, keep_full: np.ndarray
    ) -> np.ndarray:
        """
        Safe filter: return idx_full restricted to positions where keep_full[pos] is True.

        Defensive against stale/out-of-bounds cached indices (observed in your IndexError).
        """
        if idx_full is None or idx_full.size == 0:
            return idx_full

        # Ensure integer array
        idx_full = np.asarray(idx_full, dtype=np.intp)

        # Guard against out-of-bounds indices (stale cache or unexpected row changes)
        n = int(keep_full.size)
        if n <= 0:
            return idx_full[:0]

        valid = (idx_full >= 0) & (idx_full < n)
        if not np.any(valid):
            return idx_full[:0]

        idxv = idx_full[valid]
        return idxv[keep_full[idxv]]

    def _get_B_cf_activity_vector(self, t: int, cf: np.ndarray) -> np.ndarray:
        """
        Compute (and cache) v = B_t @ cf as float64 activity vector.

        Cache key uses (t, id(cf), cf.shape, cf.dtype) so repeated solve-years reuse the matvec.
        """
        if not hasattr(self, "_B_cf_actvec_cache"):
            self._B_cf_actvec_cache = {}  # type: ignore[attr-defined]

        key = (int(t), int(id(cf)), tuple(cf.shape), str(cf.dtype))
        cached = self._B_cf_actvec_cache.get(key)  # type: ignore[attr-defined]
        if cached is not None:
            return cached

        # Use CSR for fast matvec
        csr = self._get_B_csr_for_t(int(t))

        # csr is (activity x flow). cf is (flow,). result is (activity,)
        v = csr @ np.asarray(cf, dtype=np.float64)

        # Ensure plain ndarray float64
        v = np.asarray(v, dtype=np.float64)
        self._B_cf_actvec_cache[key] = v  # type: ignore[attr-defined]
        return v

    def accumulate_temporalized_biosphere_score_matrix(
        self,
        base_year: int,
        supply_matrix: np.ndarray,  # shape (n_acts, n_roots)
        root_activities: np.ndarray,  # shape (n_roots,)
        cf: np.ndarray,  # shape (n_flows,)
        *,
        min_amount: float = 0.0,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """
        Score many roots for one base_year in one pass, using dense supply_matrix.

        This keeps TD semantics and supports future cf(year) by letting caller pass cf per year.
        """
        if self.B is None:
            return
        if supply_matrix.size == 0:
            return

        base_year = int(base_year)

        # Ensure score builders exist
        if not hasattr(self, "_score_year_index") or not hasattr(
            self, "_score_bulk_value"
        ):
            self.reset_scores(attribute_to_roots=True)

        # Resolve B slice
        biosphere_slice = self._get_biosphere_slice(base_year, debug)
        if biosphere_slice is None:
            return
        _scenario_year, t, _B_t, _ = biosphere_slice

        # Map calendar year -> year_idx
        year_to_idx = self._score_year_index
        base_year_idx = year_to_idx.get(base_year)
        if base_year_idx is None:
            return

        # Validate CF
        cf = np.asarray(cf, dtype=np.float64)
        if cf.ndim != 1 or cf.size != int(self.B.shape[2]):
            raise ValueError("cf must be 1D and aligned to B flow dimension")

        # Validate shapes
        X = np.asarray(supply_matrix, dtype=np.float64)
        if X.ndim != 2:
            raise ValueError("supply_matrix must be 2D (n_acts, n_roots)")
        roots = np.asarray(root_activities, dtype=np.int64)
        if roots.ndim != 1 or roots.size != X.shape[1]:
            raise ValueError("root_activities must be 1D length n_roots")

        n_acts = int(self.B.shape[1])
        if X.shape[0] != n_acts:
            raise ValueError(f"supply_matrix has {X.shape[0]} acts; expected {n_acts}")

        # ---------- No-TD fast path (avoid dense S = v[:, None] * X) ----------
        # ---------- No-TD fast path (chunked; avoids full S allocation) ----------
        if not use_temporal_distributions or not self.temporal_biosphere_exchanges:
            # v[a] = sum_f B[a,f] * cf[f]  (cached per t/cf)
            v = self._get_B_cf_activity_vector(int(t), cf)  # float64 (n_acts,)

            # We will append sparse triplets (act, year_idx, root) in blocks
            year_idx_scalar = int(base_year_idx)

            # Choose a block size that fits cache; tune if needed
            block = 2048

            n_acts, n_roots = X.shape

            # If X is float64 already, keep it; otherwise convert once
            # (UMFPACK returns float64; so this usually does nothing)
            if X.dtype != np.float64:
                X = X.astype(np.float64, copy=False)

            out_act_parts: list[np.ndarray] = []
            out_year_parts: list[np.ndarray] = []
            out_root_parts: list[np.ndarray] = []
            out_val_parts: list[np.ndarray] = []

            for a0 in range(0, n_acts, block):
                a1 = min(a0 + block, n_acts)

                v_blk = v[a0:a1]  # (blk,)
                if not np.any(v_blk):  # all zeros
                    continue

                X_blk = X[a0:a1, :]  # (blk, n_roots)

                # Quick reject: if X block has no nonzeros, skip
                # (works even if dense; cheap check)
                if not np.any(X_blk):
                    continue

                # No threshold: only skip exact zeros
                S_blk = v_blk[:, None] * X_blk
                M = S_blk != 0.0
                if not np.any(M):
                    continue

                rr, cc = np.nonzero(M)
                vals = S_blk[rr, cc]

                acts = (a0 + rr).astype(np.int64, copy=False)
                years = np.full(acts.size, year_idx_scalar, dtype=np.int64)
                roots_out = roots[cc].astype(np.int64, copy=False)

                out_act_parts.append(acts)
                out_year_parts.append(years)
                out_root_parts.append(roots_out)
                out_val_parts.append(np.asarray(vals, dtype=np.float64))

            if not out_val_parts:
                return

            act_all = np.concatenate(out_act_parts)
            year_all = np.concatenate(out_year_parts)
            root_all = np.concatenate(out_root_parts)
            val_all = np.concatenate(out_val_parts)

            self._append_scores_bulk(
                act_idx=act_all,
                year_idx=year_all,
                values=val_all,
                root_activity=root_all,
            )
            return

        # ---------- TD enabled (keep semantics) ----------
        # For TD, we reuse your existing per-activity row-char cache logic,
        # but multiply scalars by X[a, :] across roots instead of by one supply_amt.

        # Pull row cache for anchor-year
        row_ptr, flow_sorted, data_sorted = self._get_B_row_cache_for_t(int(t))

        tpl_label = str(self._map_year_to_template_year(base_year))
        bio_td_get = self.temporal_biosphere_exchanges.get

        # Caches (reuse the same ones as scalar method)
        if not hasattr(self, "_td_pulse_cache"):
            self._td_pulse_cache = {}
        pulse_cache = self._td_pulse_cache

        if not hasattr(self, "_bio_score_row_char_cache"):
            self._bio_score_row_char_cache = {}
        row_char_cache = self._bio_score_row_char_cache

        def td_key(tex: TemporalExchange) -> tuple:
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
            )

        def pulses_from_key(k: tuple) -> list[tuple[int, float]]:
            dist, loc, scale, off_min, off_max, amt_src = k
            tex = TemporalExchange(
                distribution=dist,
                loc=loc,
                scale=scale,
                offset_min=off_min,
                offset_max=off_max,
                amount_source=amt_src,
            )
            return [
                (int(o), float(w))
                for o, w in TemporalDistribution(tex).iter_offsets_and_weights(
                    debug=False
                )
            ]

        # We will append contributions for (act, year_idx, root) in bulk by accumulating
        # sparse triplets in local python lists then flush via _append_scores_bulk per group.
        out_act = []
        out_year = []
        out_root = []
        out_val = []

        # Iterate only activities that actually contribute to any root
        active_acts = np.where(np.any(X != 0.0, axis=1))[0]
        for a in active_acts:
            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            flows_full = flow_sorted[start:end].astype(np.intp, copy=False)
            vals_full = data_sorted[start:end].astype(np.float64, copy=False)

            cf_key = int(id(cf))
            cache_key = (tpl_label, int(t), int(a), cf_key)
            cached = row_char_cache.get(cache_key)

            if cached is None:
                no_td_pos = []
                port_groups_pos: dict[tuple, list[int]] = {}
                matrix_entries_pos: list[tuple[int, TemporalExchange]] = []

                for p, f in enumerate(flows_full):
                    tex = bio_td_get((tpl_label, int(a), int(f)))
                    if tex is None:
                        no_td_pos.append(p)
                        continue
                    if getattr(tex, "amount_source", "port") == "matrix":
                        matrix_entries_pos.append((p, tex))
                    else:
                        k = td_key(tex)
                        port_groups_pos.setdefault(k, []).append(p)

                # Pre-characterize anchor-year coefficients
                no_td_coeff = (
                    float(np.dot(vals_full[no_td_pos], cf[flows_full[no_td_pos]]))
                    if no_td_pos
                    else 0.0
                )

                ported_coeffs = {}
                for k, plist in port_groups_pos.items():
                    pos = np.asarray(plist, dtype=np.intp)
                    ported_coeffs[k] = (
                        float(np.dot(vals_full[pos], cf[flows_full[pos]]))
                        if pos.size
                        else 0.0
                    )

                grouped: dict[tuple, list[int]] = {}
                for p, tex in matrix_entries_pos:
                    grouped.setdefault(td_key(tex), []).append(int(p))
                matrix_groups = {
                    k: np.asarray(v, dtype=np.intp) for k, v in grouped.items()
                }

                cached = (no_td_coeff, ported_coeffs, matrix_groups)
                row_char_cache[cache_key] = cached

            no_td_coeff, ported_coeffs, matrix_groups = cached

            # X row across roots
            x_row = X[a, :]  # (n_roots,)

            # 1) No TD at base year
            if no_td_coeff != 0.0:
                vals = no_td_coeff * x_row
                m = vals != 0.0
                if np.any(m):
                    r_idx = np.where(m)[0]
                    out_act.append(np.full(r_idx.size, int(a), dtype=np.int64))
                    out_year.append(
                        np.full(r_idx.size, int(base_year_idx), dtype=np.int64)
                    )
                    out_root.append(roots[r_idx])
                    out_val.append(vals[r_idx])

            # 2) Ported TD groups
            for k, coeff_k in ported_coeffs.items():
                if coeff_k == 0.0:
                    continue
                pulses = pulse_cache.get(k)
                if pulses is None:
                    pulses = pulses_from_key(k)
                    pulse_cache[k] = pulses
                if not pulses:
                    continue

                vals_anchor = coeff_k * x_row  # per-root
                for offset, weight in pulses:
                    if weight == 0.0:
                        continue
                    raw = int(base_year + offset)
                    y_clamped = self._clamp_year_to_scores(raw)
                    yidx = year_to_idx[int(y_clamped)]

                    vals = vals_anchor * float(weight)
                    m = vals != 0.0
                    if np.any(m):
                        r_idx = np.where(m)[0]
                        out_act.append(np.full(r_idx.size, int(a), dtype=np.int64))
                        out_year.append(np.full(r_idx.size, int(yidx), dtype=np.int64))
                        out_root.append(roots[r_idx])
                        out_val.append(vals[r_idx])

            # 3) Matrix-sourced groups
            # (kept exactly correct; still not cheap, but now shared across roots)
            if matrix_groups:
                # minimal implementation: compute score_per_supply for each pulse year as you do,
                # then multiply by x_row.
                scenario_index_get = self.scenario_index.get
                year_map_cache: dict[int, int] = {}
                t_eff_cache: dict[int, int | None] = {}
                B_row_cache_local: dict[
                    int, tuple[np.ndarray, np.ndarray, np.ndarray]
                ] = {}

                def map_year_cached(raw_year: int) -> int:
                    y = year_map_cache.get(raw_year)
                    if y is None:
                        y = int(self._map_year_to_scenario_year(raw_year))
                        year_map_cache[raw_year] = y
                    return y

                for k, idx_full in matrix_groups.items():
                    if idx_full.size == 0:
                        continue
                    f_arr = flows_full[idx_full]
                    if f_arr.size == 0:
                        continue
                    ord_f = np.argsort(f_arr, kind="mergesort")
                    f_sorted = f_arr[ord_f].astype(np.intp, copy=False)
                    cf_sorted = cf[f_sorted]

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        pulses = pulses_from_key(k)
                        pulse_cache[k] = pulses
                    if not pulses:
                        continue

                    for offset, weight in pulses:
                        if weight == 0.0:
                            continue

                        raw = int(base_year + offset)
                        y_clamped = self._clamp_year_to_scores(raw)
                        yidx = year_to_idx[int(y_clamped)]

                        y_eff = map_year_cached(raw)
                        t_eff = t_eff_cache.get(y_eff)
                        if t_eff is None and y_eff not in t_eff_cache:
                            t_eff = scenario_index_get(str(y_eff))
                            t_eff_cache[y_eff] = t_eff
                        if t_eff is None:
                            continue

                        t_eff_i = int(t_eff)
                        cached_eff = B_row_cache_local.get(t_eff_i)
                        if cached_eff is None:
                            cached_eff = self._get_B_row_cache_for_t(t_eff_i)
                            B_row_cache_local[t_eff_i] = cached_eff
                        row_ptr_eff, flow_sorted_eff, data_sorted_eff = cached_eff

                        start_eff = int(row_ptr_eff[a])
                        end_eff = int(row_ptr_eff[a + 1])
                        if start_eff == end_eff:
                            continue

                        row_vals_eff = data_sorted_eff[start_eff:end_eff].astype(
                            np.float64, copy=False
                        )

                        # Map flows -> position in this row (cached)
                        row_index_map = self._get_B_row_index_map_for_t_act(t_eff_i, a)

                        # Extract the B-values aligned to f_sorted (safe even if row flows are unsorted)
                        matrix_kernel = self._get_numba_matrix_kernel()
                        if matrix_kernel is not None:
                            vals_eff, valid = matrix_kernel(
                                f_sorted.astype(np.int64, copy=False),
                                row_index_map,
                                row_vals_eff,
                            )
                            if not np.any(valid):
                                continue
                            vals_eff = vals_eff[valid]
                            cf_use = cf_sorted[valid]
                        else:
                            pos = row_index_map[f_sorted]
                            valid = pos >= 0
                            if not np.any(valid):
                                continue
                            vals_eff = row_vals_eff[pos[valid]]
                            cf_use = cf_sorted[valid]

                        score_per_supply = float(np.dot(vals_eff, cf_use)) * float(
                            weight
                        )

                        if score_per_supply == 0.0:
                            continue

                        vals = score_per_supply * x_row
                        m = vals != 0.0
                        if np.any(m):
                            r_idx = np.where(m)[0]
                            out_act.append(np.full(r_idx.size, int(a), dtype=np.int64))
                            out_year.append(
                                np.full(r_idx.size, int(yidx), dtype=np.int64)
                            )
                            out_root.append(roots[r_idx])
                            out_val.append(vals[r_idx])

        # Flush
        if out_val:
            act = np.concatenate(out_act)
            yr = np.concatenate(out_year)
            root = np.concatenate(out_root)
            val = np.concatenate(out_val).astype(np.float64, copy=False)

            self._append_scores_bulk(act, yr, val, root_activity=root)

    def accumulate_temporalized_biosphere_score(
        self,
        base_year: int,
        supply_by_activity: Dict[int, float],
        cf: np.ndarray,
        *,
        min_amount: float = 0.0,
        store_activity: int | None = None,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """
        Accumulate *characterized* biosphere impacts directly into Trails.scores.

        Radical speed approach:
          - Cache per-(t, act) *characterized* row coefficients:
              no-TD coeff:     sum_f B[a,f] * cf[f] over flows with no TD
              ported group:    sum_{f in group(k)} B[a,f] * cf[f]
            so runtime becomes mostly scalar arithmetic per activity.
          - Matrix-sourced TD stays year-dependent (must read B at pulse years),
            but overhead is minimized.

        Also fixes the IndexError by defensively validating any cached position arrays
        against the current row length.
        """
        # Ensure score builders exist
        if not hasattr(self, "_score_year_index") or not hasattr(
            self, "_score_chunk_value"
        ):
            self.reset_scores(
                attribute_to_roots=bool(getattr(self, "_scores_has_root", False))
            )

        if not supply_by_activity:
            return
        if self.B is None:
            return

        base_year = int(base_year)

        biosphere_slice = self._get_biosphere_slice(base_year, debug)
        if biosphere_slice is None:
            return
        scenario_year, t, B_t, _n_flows = biosphere_slice

        # Calendar year -> score year index
        year_to_idx = self._score_year_index
        base_year_idx = year_to_idx.get(base_year)
        if base_year_idx is None:
            return

        # Validate CF
        cf = np.asarray(cf, dtype=np.float64)
        if cf.ndim != 1:
            raise ValueError(
                "cf must be a 1D vector aligned to trails.B flow dimension."
            )
        if cf.size != int(self.B.shape[2]):
            raise ValueError(
                f"cf length {cf.size} does not match B flows {int(self.B.shape[2])}"
            )

        scenario_index_get = self.scenario_index.get
        map_year_to_scenario = self._map_year_to_scenario_year

        # TD metadata lookup (template-year stable)
        bio_td = (
            self.temporal_biosphere_exchanges if use_temporal_distributions else None
        )
        if bio_td:
            tpl_label = str(self._map_year_to_template_year(base_year))
            bio_td_get = bio_td.get
        else:
            tpl_label = None
            bio_td_get = None  # type: ignore[assignment]

        # --- persistent caches on self ---
        # Cache: pulses for TD parameter key
        if not hasattr(self, "_td_pulse_cache"):
            self._td_pulse_cache = {}  # type: ignore[attr-defined]
        pulse_cache = self._td_pulse_cache  # type: ignore[attr-defined]

        # Cache: for each (tpl_label, t, act) store pre-characterized coefficients:
        #   (no_td_coeff: float,
        #    ported_coeffs: dict[k, float],
        #    matrix_groups: dict[k, np.ndarray[pos]]   # positions within the row
        #   )
        if not hasattr(self, "_bio_score_row_char_cache"):
            self._bio_score_row_char_cache = {}  # type: ignore[attr-defined]
        row_char_cache = self._bio_score_row_char_cache  # type: ignore[attr-defined]

        # Per-call caches
        year_map_cache: dict[int, int] = {}
        t_eff_cache: dict[int, int | None] = {}
        B_row_cache_local: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        def map_year_cached(raw_year: int) -> int:
            y = year_map_cache.get(raw_year)
            if y is None:
                y = int(map_year_to_scenario(raw_year))
                year_map_cache[raw_year] = y
            return y

        def td_key(tex: TemporalExchange) -> tuple:
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
            )

        def pulses_from_key(k: tuple) -> list[tuple[int, float]]:
            dist, loc, scale, off_min, off_max, amt_src = k
            tex = TemporalExchange(
                distribution=dist,
                loc=loc,
                scale=scale,
                offset_min=off_min,
                offset_max=off_max,
                amount_source=amt_src,
            )
            return [
                (int(o), float(w))
                for o, w in TemporalDistribution(tex).iter_offsets_and_weights(
                    debug=False
                )
            ]

        # Anchor-year CSR-like row cache (fast row slicing)
        row_ptr, flow_sorted, data_sorted = self._get_B_row_cache_for_t(int(t))

        def _safe_filter_positions(
            pos: np.ndarray, keep: np.ndarray | None, row_len: int
        ) -> np.ndarray:
            """
            pos: positions within [0, row_len)
            keep: boolean mask of length row_len or None
            Returns filtered positions, safely (never IndexError).
            """
            if pos.size == 0:
                return pos
            # guard against stale/corrupt cached positions
            if pos.max(initial=-1) >= row_len or pos.min(initial=0) < 0:
                pos = pos[(pos >= 0) & (pos < row_len)]
                if pos.size == 0:
                    return pos
            if keep is None:
                return pos
            # keep[pos] is now safe
            return pos[keep[pos]]

        # Main loop
        for act_idx, supply_amt in supply_by_activity.items():
            supply_amt = float(supply_amt)
            if supply_amt == 0.0:
                continue

            a = int(act_idx)
            if a < 0 or a + 1 >= len(row_ptr):
                continue

            # Store into which activity index?
            has_root = bool(getattr(self, "_scores_has_root", False))
            if has_root:
                score_act = a
                root_activity = int(store_activity) if store_activity is not None else a
            else:
                score_act = int(store_activity) if store_activity is not None else a
                root_activity = None

            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            row_len = end - start
            flows_full = flow_sorted[start:end].astype(np.intp, copy=False)
            vals_full = data_sorted[start:end].astype(np.float64, copy=False)

            # Filter zeros at the contribution level to avoid useless work for huge rows.

            # ---- NO TD fast path: score = supply_amt * sum_f B[a,f]*cf[f] ----
            if not bio_td:
                score = supply_amt * float(np.dot(vals_full, cf[flows_full]))

                if score != 0.0:
                    self._append_scores_from_yearidx_map(
                        score_act, {base_year_idx: score}, root_activity=root_activity
                    )
                continue

            # ---- TD enabled: use cached pre-characterized row structure ----
            cf_key = int(id(cf))
            cache_key = (tpl_label, int(t), int(a), cf_key)

            cached = row_char_cache.get(cache_key)

            if cached is None:
                # Build once per (tpl_label, t, act)
                no_td_pos: list[int] = []
                port_groups_pos: dict[tuple, list[int]] = {}
                matrix_entries_pos: list[tuple[int, TemporalExchange]] = []

                # classify flows by TD
                for p, f in enumerate(flows_full):
                    tex = bio_td_get((tpl_label, a, int(f)))  # type: ignore[misc]
                    if tex is None:
                        no_td_pos.append(p)
                        continue
                    amount_source = getattr(tex, "amount_source", "port")
                    if amount_source == "matrix":
                        matrix_entries_pos.append((p, tex))
                    else:
                        k = td_key(tex)
                        port_groups_pos.setdefault(k, []).append(p)

                # Pre-characterize coefficients (UNSCALED by supply):
                #   coeff = sum_{pos} B_val[pos] * cf[flow[pos]]
                no_td_coeff = 0.0
                if no_td_pos:
                    pos = np.asarray(no_td_pos, dtype=np.intp)
                    # safe by construction here, but keep it clean:
                    if pos.size:
                        no_td_coeff = float(np.dot(vals_full[pos], cf[flows_full[pos]]))

                ported_coeffs: dict[tuple, float] = {}
                for k, plist in port_groups_pos.items():
                    pos = np.asarray(plist, dtype=np.intp)
                    if pos.size:
                        ported_coeffs[k] = float(
                            np.dot(vals_full[pos], cf[flows_full[pos]])
                        )
                    else:
                        ported_coeffs[k] = 0.0

                # Matrix groups store positions (year-dependent values, so cannot pre-characterize)
                grouped: dict[tuple, list[int]] = {}
                for p, tex in matrix_entries_pos:
                    k = td_key(tex)
                    grouped.setdefault(k, []).append(int(p))
                matrix_groups = {
                    k: np.asarray(v, dtype=np.intp) for k, v in grouped.items()
                }

                cached = (no_td_coeff, ported_coeffs, matrix_groups)
                row_char_cache[cache_key] = cached

            no_td_coeff, ported_coeffs, matrix_groups = cached

            # Per-activity accumulator: year_idx -> score
            acc_yearidx: dict[int, float] = {}

            # 1) No TD subset (anchor year): supply * no_td_coeff, with optional flow filtering
            # If keep_full exists, we must recompute the no-TD coefficient on the kept subset.
            if no_td_coeff != 0.0:
                score = supply_amt * float(no_td_coeff)
                if score != 0.0:
                    acc_yearidx[base_year_idx] = (
                        acc_yearidx.get(base_year_idx, 0.0) + score
                    )

            # 2) Ported TD groups: distribute scalar (supply * coeff_k) across pulse years
            if ported_coeffs:
                for k, coeff_k in ported_coeffs.items():
                    if coeff_k == 0.0:
                        continue
                    score_anchor = supply_amt * float(coeff_k)

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        pulses = pulses_from_key(k)
                        pulse_cache[k] = pulses

                    if not pulses:
                        continue

                    for offset, weight in pulses:
                        if weight == 0.0:
                            continue
                        raw = int(base_year + offset)
                        y_clamped = self._clamp_year_to_scores(raw)
                        yidx = year_to_idx[int(y_clamped)]

                        acc_yearidx[yidx] = acc_yearidx.get(
                            yidx, 0.0
                        ) + score_anchor * float(weight)

            # 3) Matrix-sourced TD groups: year-dependent row lookup + dot on that subset
            if matrix_groups:
                for k, idx_full in matrix_groups.items():
                    if idx_full.size == 0:
                        continue

                    # Apply keep_full safely (never IndexError)
                    idx = _safe_filter_positions(idx_full, None, row_len)

                    if idx.size == 0:
                        continue

                    # subset flows for this matrix group (positions within row)
                    f_arr = flows_full[idx]
                    if f_arr.size == 0:
                        continue

                    ord_f = np.argsort(f_arr, kind="mergesort")
                    f_sorted = f_arr[ord_f].astype(np.intp, copy=False)
                    cf_sorted = cf[f_sorted]

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        pulses = pulses_from_key(k)
                        pulse_cache[k] = pulses
                    if not pulses:
                        continue

                    for offset, weight in pulses:
                        if weight == 0.0:
                            continue

                        raw = int(base_year + offset)
                        y_clamped = self._clamp_year_to_scores(raw)
                        yidx = year_to_idx[int(y_clamped)]

                        y_eff = map_year_cached(raw)
                        t_eff = t_eff_cache.get(y_eff)
                        if t_eff is None and y_eff not in t_eff_cache:
                            t_eff = scenario_index_get(str(y_eff))
                            t_eff_cache[y_eff] = t_eff
                        if t_eff is None:
                            continue

                        t_eff_i = int(t_eff)
                        cached_eff = B_row_cache_local.get(t_eff_i)
                        if cached_eff is None:
                            cached_eff = self._get_B_row_cache_for_t(t_eff_i)
                            B_row_cache_local[t_eff_i] = cached_eff
                        row_ptr_eff, flow_sorted_eff, data_sorted_eff = cached_eff

                        start_eff = int(row_ptr_eff[a])
                        end_eff = int(row_ptr_eff[a + 1])
                        if start_eff == end_eff:
                            continue

                        row_flows_eff = flow_sorted_eff[start_eff:end_eff].astype(
                            np.intp, copy=False
                        )
                        row_vals_eff = data_sorted_eff[start_eff:end_eff].astype(
                            np.float64, copy=False
                        )

                        if row_flows_eff.size > 1 and np.any(
                            row_flows_eff[1:] < row_flows_eff[:-1]
                        ):
                            order = np.argsort(row_flows_eff, kind="mergesort")
                            row_flows_eff = row_flows_eff[order]
                            row_vals_eff = row_vals_eff[order]

                        vals_sorted = self._row_values_for_flows_sorted(
                            row_flows_eff, row_vals_eff, f_sorted
                        )
                        # score per unit supply:
                        score_per_supply = float(np.dot(vals_sorted, cf_sorted))
                        score = supply_amt * float(weight) * score_per_supply

                        if score != 0.0:
                            acc_yearidx[yidx] = acc_yearidx.get(yidx, 0.0) + score

            if acc_yearidx:
                self._append_scores_from_yearidx_map(
                    score_act, acc_yearidx, root_activity=root_activity
                )

    def _accumulate_temporalized_biosphere_inventory_core(
        self,
        ctx: _BioAccumulationContext,
        supply_by_activity: Dict[int, float],
        *,
        min_amount: float,
        store_activity: int | None,
        use_temporal_distributions: bool,
        debug: bool,
    ) -> None:
        """Core accumulation routine reused across batch calls."""
        if not supply_by_activity:
            return

        scenario_index_get = ctx.scenario_index_get
        map_year_to_scenario = ctx.map_year_to_scenario
        row_ptr = ctx.row_ptr
        flow_sorted = ctx.flow_sorted
        data_sorted = ctx.data_sorted

        base_year = int(ctx.base_year)

        bio_td_get = ctx.bio_td_get
        tpl_label = ctx.tpl_label
        td_expanded_cache = ctx.td_expanded_cache
        pulse_cache = ctx.pulse_cache

        year_map_cache = ctx.year_map_cache
        t_eff_cache = ctx.t_eff_cache
        B_row_cache_local = ctx.B_row_cache_local

        if use_temporal_distributions and ctx.bio_td_get is None:
            use_temporal_distributions = False

        def map_year_cached(raw_year: int) -> int:
            y = year_map_cache.get(raw_year)
            if y is None:
                y = int(map_year_to_scenario(raw_year))
                year_map_cache[raw_year] = y
            return y

        def td_key(tex: TemporalExchange) -> tuple:
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
            )

        if not use_temporal_distributions:
            supply_vec = np.zeros(ctx.n_acts, dtype=np.float64)
            for act_idx, supply_amt in supply_by_activity.items():
                a_idx = int(act_idx)
                if 0 <= a_idx < ctx.n_acts:
                    supply_vec[a_idx] = float(supply_amt)

            scaled = (
                ctx.data.astype(np.float64, copy=False) * supply_vec[ctx.act_coords]
            )

            act_coords = ctx.act_coords
            flow_coords = ctx.flow_coords

            has_root = bool(getattr(self, "_inventory_has_root", False))
            if has_root:
                inventory_act = act_coords
                root_activity = (
                    np.full_like(act_coords, int(store_activity))
                    if store_activity is not None
                    else act_coords
                )
            else:
                if store_activity is not None:
                    inventory_act = np.full_like(act_coords, int(store_activity))
                else:
                    inventory_act = act_coords
                root_activity = None

            self._append_inventory_entries_bulk(
                inventory_act,
                base_year,
                flow_coords,
                scaled,
                root_activity=root_activity,
            )
            return

        for act_idx, supply_amt in supply_by_activity.items():
            supply_amt = float(supply_amt)
            if supply_amt == 0.0:
                continue

            a = int(act_idx)
            has_root = bool(getattr(self, "_inventory_has_root", False))
            if has_root:
                inventory_act = a
                root_activity = int(store_activity) if store_activity is not None else a
            else:
                inventory_act = int(store_activity) if store_activity is not None else a
                root_activity = None
            if a < 0 or a + 1 >= len(row_ptr):
                continue

            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            flows_full = flow_sorted[start:end].astype(np.intp, copy=False)
            vals_full = data_sorted[start:end]

            scaled_full = supply_amt * vals_full.astype(np.float64, copy=False)

            thr = float(min_amount)
            if thr > 0.0:
                temporalize = np.abs(scaled_full) >= thr
            else:
                temporalize = None

            # TD metadata is keyed by (tpl_label, act, flow) and does not depend on scenario slice t.
            cache_key = (tpl_label, int(a))

            td_struct = td_expanded_cache.get(cache_key)

            if td_struct is None:
                no_td_pos: list[int] = []
                port_groups_pos: dict[tuple, list[int]] = {}
                matrix_entries: list[tuple[int, TemporalExchange]] = []

                for p, f in enumerate(flows_full):
                    tex = bio_td_get((tpl_label, a, int(f)))  # type: ignore[misc]
                    if tex is None:
                        no_td_pos.append(p)
                        continue

                    amount_source = getattr(tex, "amount_source", "port")
                    if amount_source == "matrix":
                        matrix_entries.append((p, tex))
                    else:
                        k = td_key(tex)
                        port_groups_pos.setdefault(k, []).append(p)

                no_td_idx = np.array(no_td_pos, dtype=np.intp) if no_td_pos else None

                # Cache base-row positions by TD key (no pulse expansion here)
                ported_groups = {
                    k: np.asarray(v, dtype=np.intp)
                    for k, v in port_groups_pos.items()
                    if v
                }

                td_struct = (
                    no_td_idx,
                    ported_groups,  # <-- dict[k, np.ndarray[pos]]
                    matrix_entries,
                )
                td_expanded_cache[cache_key] = td_struct

            no_td_idx, ported_groups, matrix_entries = td_struct

            if no_td_idx is not None:
                idx = no_td_idx
                if idx.size:
                    self._append_inventory_entries(
                        inventory_act,
                        base_year,
                        flows_full[idx],
                        scaled_full[idx],
                        root_activity=root_activity,
                    )

            # -------------------------
            # PORTED TD groups (min_amount controls temporalization only)
            # -------------------------
            if ported_groups:
                row_len = int(end - start)
                thr = float(min_amount)

                for k, idx_full in ported_groups.items():
                    if idx_full is None or idx_full.size == 0:
                        continue

                    # Defensive: ensure indices are within current row bounds (stale cache protection)
                    if (
                        idx_full.max(initial=-1) >= row_len
                        or idx_full.min(initial=0) < 0
                    ):
                        # Cache is stale -> drop it and rebuild next time
                        td_expanded_cache.pop(cache_key, None)
                        idx_full = idx_full[(idx_full >= 0) & (idx_full < row_len)]
                        if idx_full.size == 0:
                            continue

                    if temporalize is None:
                        # No thresholding requested: everything temporalized
                        idx_td = idx_full
                        idx_anchor = None
                    else:
                        # Below threshold -> anchor to base_year (NOT omitted)
                        idx_anchor = idx_full[~temporalize[idx_full]]
                        idx_td = idx_full[temporalize[idx_full]]

                    # 1) Anchor below-threshold contributions at base_year
                    if idx_anchor is not None and idx_anchor.size:
                        self._append_inventory_entries(
                            inventory_act,
                            base_year,
                            flows_full[idx_anchor],
                            scaled_full[idx_anchor],
                            root_activity=root_activity,
                        )

                    # 2) Temporalize above-threshold contributions
                    if idx_td.size == 0:
                        continue

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        # Reconstruct a TemporalExchange from td_key tuple (your existing convention)
                        dist, loc, scale, off_min, off_max, amt_src = k
                        tex0 = TemporalExchange(
                            distribution=dist,
                            loc=loc,
                            scale=scale,
                            offset_min=off_min,
                            offset_max=off_max,
                            amount_source=amt_src,
                        )
                        pulses = [
                            (int(o), float(w))
                            for o, w in TemporalDistribution(
                                tex0
                            ).iter_offsets_and_weights(debug=False)
                        ]
                        pulse_cache[k] = pulses

                    if not pulses:
                        # If TD produces nothing, treat as anchor (still not omitted)
                        self._append_inventory_entries(
                            inventory_act,
                            base_year,
                            flows_full[idx_td],
                            scaled_full[idx_td],
                            root_activity=root_activity,
                        )
                        continue

                    # Expand pulses in a vectorized way:
                    idx_rep = np.repeat(idx_td, len(pulses))
                    offsets_arr = np.fromiter(
                        (o for o, _ in pulses), dtype=np.int64, count=len(pulses)
                    )
                    weights_arr = np.fromiter(
                        (w for _, w in pulses), dtype=np.float64, count=len(pulses)
                    )

                    offsets_rep = np.tile(offsets_arr, idx_td.size)
                    weights_rep = np.tile(weights_arr, idx_td.size)

                    # Build contributions
                    flows_use = flows_full[idx_rep]
                    years_use = base_year + offsets_rep
                    contrib = scaled_full[idx_rep] * weights_rep

                    acts_use = np.full_like(flows_use, inventory_act, dtype=np.int64)
                    self._append_inventory_entries_bulk(
                        acts_use,
                        years_use,
                        flows_use,
                        contrib,
                        root_activity=root_activity,
                    )

            if matrix_entries:
                if isinstance(matrix_entries, list):
                    grouped: dict[tuple, tuple[list[int], TemporalExchange]] = {}
                    for p, tex in matrix_entries:
                        k = td_key(tex)
                        if k in grouped:
                            grouped[k][0].append(int(p))
                        else:
                            grouped[k] = ([int(p)], tex)
                    matrix_groups = {}
                    for k, (pos_list, tex0) in grouped.items():
                        idx_full = np.array(pos_list, dtype=np.intp)
                        pulses = pulse_cache.get(k)
                        if pulses is None:
                            pulses = [
                                (int(o), float(w))
                                for o, w in TemporalDistribution(
                                    tex0
                                ).iter_offsets_and_weights(debug=False)
                            ]
                            pulse_cache[k] = pulses
                        offsets_arr = np.array([o for o, _ in pulses], dtype=np.int64)
                        weights_arr = np.array([w for _, w in pulses], dtype=np.float64)
                        matrix_groups[k] = (idx_full, offsets_arr, weights_arr)

                    td_struct = (
                        no_td_idx,
                        ported_groups,
                        matrix_groups,
                    )
                    td_expanded_cache[cache_key] = td_struct
                    matrix_entries = matrix_groups

                for k, (idx_full, offsets_arr, weights_arr) in matrix_entries.items():  # type: ignore[union-attr]
                    idx = idx_full
                    if temporalize is not None:
                        # anchor below-threshold
                        idx_anchor = idx[~temporalize[idx]]
                        if idx_anchor.size:
                            self._append_inventory_entries(
                                inventory_act,
                                base_year,
                                flows_full[idx_anchor],
                                scaled_full[idx_anchor],
                                root_activity=root_activity,
                            )
                        # TD only for above-threshold
                        idx = idx[temporalize[idx]]
                        if idx.size == 0:
                            continue

                    f_arr = flows_full[idx]

                    year_groups: dict[int, list[tuple[int, float]]] = {}
                    for offset, weight in zip(offsets_arr, weights_arr):
                        if weight == 0.0:
                            continue
                        raw_year = base_year + int(offset)
                        y_eff = map_year_cached(raw_year)
                        year_groups.setdefault(int(y_eff), []).append(
                            (int(raw_year), float(weight))
                        )

                    for y_eff, year_weights in year_groups.items():
                        t_eff = t_eff_cache.get(y_eff)
                        if t_eff is None and y_eff not in t_eff_cache:
                            t_eff = scenario_index_get(str(y_eff))
                            t_eff_cache[y_eff] = t_eff
                        if t_eff is None:
                            continue

                        t_eff_i = int(t_eff)
                        cached_eff = B_row_cache_local.get(t_eff_i)
                        if cached_eff is None:
                            cached_eff = self._get_B_row_cache_for_t(t_eff_i)
                            B_row_cache_local[t_eff_i] = cached_eff
                        row_ptr_eff, flow_sorted_eff, data_sorted_eff = cached_eff

                        start_eff = int(row_ptr_eff[a])
                        end_eff = int(row_ptr_eff[a + 1])
                        if start_eff == end_eff:
                            continue

                        row_vals_eff = data_sorted_eff[start_eff:end_eff]
                        if f_arr.size == 0:
                            continue
                        row_index_map = self._get_B_row_index_map_for_t_act(t_eff_i, a)
                        matrix_kernel = self._get_numba_matrix_kernel()
                        if matrix_kernel is not None:
                            vals_eff, valid = matrix_kernel(
                                f_arr.astype(np.int64, copy=False),
                                row_index_map,
                                row_vals_eff,
                            )
                            if not np.any(valid):
                                continue
                            vals_eff = vals_eff[valid]
                            f_use = f_arr[valid]
                        else:
                            pos = row_index_map[f_arr]
                            valid = pos >= 0
                            if not np.any(valid):
                                continue
                            vals_eff = row_vals_eff[pos[valid]]
                            f_use = f_arr[valid]

                        years_vec = np.array(
                            [yw[0] for yw in year_weights], dtype=np.int64
                        )
                        weights_vec = np.array(
                            [yw[1] for yw in year_weights], dtype=np.float64
                        )

                        flows_rep = np.repeat(f_use, weights_vec.size)
                        years_rep = np.tile(years_vec, f_use.size)
                        contrib = (
                            supply_amt
                            * np.repeat(
                                vals_eff.astype(np.float64, copy=False),
                                weights_vec.size,
                            )
                            * np.tile(weights_vec, f_use.size)
                        )

                        acts_use = np.full_like(
                            flows_rep, inventory_act, dtype=np.int64
                        )
                        self._append_inventory_entries_bulk(
                            acts_use,
                            years_rep,
                            flows_rep,
                            contrib,
                            root_activity=root_activity,
                        )

    def _get_numba_ported_kernel(self) -> Callable | None:
        """Return a cached Numba kernel for ported TD accumulation, if available."""
        if hasattr(self, "_numba_ported_kernel_checked"):
            return getattr(self, "_numba_ported_kernel", None)

        setattr(self, "_numba_ported_kernel_checked", True)
        spec = importlib.util.find_spec("numba")
        if spec is None:
            setattr(self, "_numba_ported_kernel", None)
            return None

        numba = importlib.import_module("numba")
        np_local = np

        @numba.njit(cache=True)
        def _ported_kernel(
            flows_full: np_local.ndarray,
            scaled_full: np_local.ndarray,
            ported_flow_idx: np_local.ndarray,
            ported_offsets: np_local.ndarray,
            ported_weights: np_local.ndarray,
            base_year: int,
        ) -> tuple[np_local.ndarray, np_local.ndarray, np_local.ndarray]:
            n = ported_flow_idx.size
            flows_out = np_local.empty(n, dtype=np_local.int64)
            years_out = np_local.empty(n, dtype=np_local.int64)
            contrib_out = np_local.empty(n, dtype=np_local.float64)
            for i in range(n):
                idx = ported_flow_idx[i]
                flows_out[i] = flows_full[idx]
                years_out[i] = base_year + ported_offsets[i]
                contrib_out[i] = scaled_full[idx] * ported_weights[i]
            return flows_out, years_out, contrib_out

        setattr(self, "_numba_ported_kernel", _ported_kernel)
        return _ported_kernel

    def _get_numba_matrix_kernel(self) -> Callable | None:
        """Return a cached Numba kernel for matrix TD lookups, if available."""
        if hasattr(self, "_numba_matrix_kernel_checked"):
            return getattr(self, "_numba_matrix_kernel", None)

        setattr(self, "_numba_matrix_kernel_checked", True)
        spec = importlib.util.find_spec("numba")
        if spec is None:
            setattr(self, "_numba_matrix_kernel", None)
            return None

        numba = importlib.import_module("numba")
        np_local = np

        @numba.njit(cache=True)
        def _matrix_kernel(
            flows: np_local.ndarray,
            row_index_map: np_local.ndarray,
            row_vals: np_local.ndarray,
        ) -> tuple[np_local.ndarray, np_local.ndarray]:
            n = flows.size
            vals = np_local.empty(n, dtype=np_local.float64)
            valid = np_local.zeros(n, dtype=np_local.bool_)
            for i in range(n):
                pos = row_index_map[flows[i]]
                if pos >= 0:
                    vals[i] = row_vals[pos]
                    valid[i] = True
            return vals, valid

        setattr(self, "_numba_matrix_kernel", _matrix_kernel)
        return _matrix_kernel

    def _accumulate_no_td_batch(
        self,
        ctx: _BioAccumulationContext,
        supplies: List[tuple[Dict[int, float], int | None]],
        *,
        min_amount: float,
        debug: bool,
    ) -> bool:
        """Fast-path bulk no-TD accumulation across multiple supplies."""
        if not supplies:
            return True

        has_root = bool(getattr(self, "_inventory_has_root", False))
        if not has_root:
            return False

        if any(store_activity is None for _, store_activity in supplies):
            return False

        base_year = int(ctx.base_year)

        root_ids = np.array([int(store) for _, store in supplies], dtype=np.int64)
        n_roots = int(root_ids.size)
        if n_roots == 0:
            return True

        supply_matrix = np.zeros((ctx.n_acts, n_roots), dtype=np.float64)
        for col, (supply_by_activity, _) in enumerate(supplies):
            for act_idx, supply_amt in supply_by_activity.items():
                a_idx = int(act_idx)
                if 0 <= a_idx < ctx.n_acts:
                    supply_matrix[a_idx, col] += float(supply_amt)

        nnz = int(ctx.data.shape[0])
        if nnz == 0:
            return True

        chunk_size = 200_000
        for start in range(0, nnz, chunk_size):
            end = min(start + chunk_size, nnz)
            act_chunk = ctx.act_coords[start:end].astype(np.int64, copy=False)
            flow_chunk = ctx.flow_coords[start:end].astype(np.int64, copy=False)
            data_chunk = ctx.data[start:end].astype(np.float64, copy=False)

            supply_chunk = supply_matrix[act_chunk, :]
            values = data_chunk[:, None] * supply_chunk

            acts_rep = np.repeat(act_chunk, n_roots)
            flows_rep = np.repeat(flow_chunk, n_roots)
            roots_rep = np.tile(root_ids, int(act_chunk.size))
            values_flat = values.ravel()

            self._append_inventory_entries_bulk(
                acts_rep,
                base_year,
                flows_rep,
                values_flat,
                root_activity=roots_rep,
            )

        return True

    def accumulate_temporalized_biosphere_inventory(
        self,
        base_year: int,
        supply_by_activity: Dict[int, float],
        *,
        min_amount: float = 0.0,
        store_activity: int | None = None,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """
        Accumulate temporally shifted biosphere emissions for a solved supply vector,
        storing results in the Trails inventory builder.

        Performance strategy:
          - Iterate only over supplied activities.
          - Use cached CSR-like row structure for B_t.
          - Vectorize no-TD adds with np.add.at.
          - For PORTED TD, group row flows by TD parameter key and apply pulses to
            whole vectors (scaled values) at once, using np.add.at per (key, pulse, year).
          - Keep MATRIX TD semantics as scalar per (flow, pulse) because values depend on year.
          - NEW: Cache TD classification per (tpl_label, t, act) so we don't redo per-flow TD lookups.

        Semantics preserved:
          - No TD: anchor to scenario_year of B slice.
          - TD + ported: distribute anchor-year scaled amount across pulse years.
          - TD + matrix: read B at each pulse-year, multiply by supply and weight.
          - Optional store_activity: attribute biosphere flows to a root activity index.
        """
        ctx = self._build_bio_accumulation_context(
            base_year,
            use_temporal_distributions=use_temporal_distributions,
            debug=debug,
        )
        if ctx is None:
            return
        if debug or self.debug:
            logger.debug(
                "accumulate_bio: base_year=%d scenario_year=%d t=%d inv_years=[%s..%s]",
                int(ctx.base_year),
                int(ctx.scenario_year),
                int(ctx.t),
                (
                    int(self._inventory_years[0])
                    if self._inventory_years is not None
                    else -1
                ),
                (
                    int(self._inventory_years[-1])
                    if self._inventory_years is not None
                    else -1
                ),
            )
        self._accumulate_temporalized_biosphere_inventory_core(
            ctx,
            supply_by_activity,
            min_amount=min_amount,
            store_activity=store_activity,
            use_temporal_distributions=use_temporal_distributions,
            debug=debug,
        )

    def accumulate_temporalized_biosphere_inventory_batch(
        self,
        base_year: int,
        supplies: List[tuple[Dict[int, float], int | None]],
        *,
        min_amount: float = 0.0,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """Accumulate multiple supply vectors for the same base year in one pass."""
        if not supplies:
            return
        ctx = self._build_bio_accumulation_context(
            base_year,
            use_temporal_distributions=use_temporal_distributions,
            debug=debug,
        )
        if ctx is None:
            return
        if debug or self.debug:
            logger.debug(
                "accumulate_bio_batch: base_year=%d scenario_year=%d t=%d count=%d",
                int(ctx.base_year),
                int(ctx.scenario_year),
                int(ctx.t),
                len(supplies),
            )
        if not use_temporal_distributions:
            if self._accumulate_no_td_batch(
                ctx, supplies, min_amount=min_amount, debug=debug
            ):
                return

            merged_supply: Dict[int, float] = {}
            merged_store_activity: int | None = None
            for supply_by_activity, store_activity in supplies:
                if store_activity is not None and merged_store_activity is None:
                    merged_store_activity = int(store_activity)
                for act_idx, supply_amt in supply_by_activity.items():
                    merged_supply[int(act_idx)] = merged_supply.get(
                        int(act_idx), 0.0
                    ) + float(supply_amt)
            self._accumulate_temporalized_biosphere_inventory_core(
                ctx,
                merged_supply,
                min_amount=min_amount,
                store_activity=merged_store_activity,
                use_temporal_distributions=False,
                debug=debug,
            )
            return
        for supply_by_activity, store_activity in supplies:
            if not supply_by_activity:
                continue
            self._accumulate_temporalized_biosphere_inventory_core(
                ctx,
                supply_by_activity,
                min_amount=min_amount,
                store_activity=store_activity,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
            )

    def _map_year_to_available(self, year: int) -> int:
        """
        Backwards-compatible alias for _map_year_to_scenario_year.
        """
        return self._map_year_to_scenario_year(year)

    @staticmethod
    def _estimate_total_from_depth(max_depth: int) -> int | None:
        """Estimate a traversal size for a given maximum depth.

        :param max_depth: Maximum traversal depth.
        :type max_depth: int
        :returns: Estimated total size or ``None`` if not available.
        :rtype: int | None
        """
        DEPTH_TOTALS = {
            1: 10,
            2: 50,
            3: 400,
            4: 4000,
            5: 40000,
            6: 400000,
            7: 4000000,
            8: 40000000,
        }

        EMPIRICAL_SAFETY_FACTOR = 1.05

        if max_depth in DEPTH_TOTALS:
            return int(max(1, DEPTH_TOTALS[max_depth] * EMPIRICAL_SAFETY_FACTOR))

        depths = sorted(DEPTH_TOTALS.keys())
        if len(depths) >= 2:
            lo = max([d for d in depths if d < max_depth], default=None)
            hi = min([d for d in depths if d > max_depth], default=None)
            if (
                lo is not None
                and hi is not None
                and DEPTH_TOTALS[lo] > 0
                and DEPTH_TOTALS[hi] > 0
            ):
                import math

                y0 = math.log(DEPTH_TOTALS[lo])
                y1 = math.log(DEPTH_TOTALS[hi])
                t = (max_depth - lo) / (hi - lo)
                est = math.exp(y0 + t * (y1 - y0))
                return int(max(1, est * EMPIRICAL_SAFETY_FACTOR))

        return None

    @staticmethod
    def _record_frontier(
        frontier_total: dict[tuple[int, int], float],
        provenance_roots: dict[tuple[int, int], dict[int, float]],
        y: int,
        a: int,
        x: float,
        r: Optional[int],
        return_provenance: bool,
    ) -> None:
        """Record frontier totals and optional provenance entries.

        :param frontier_total: Mapping of (year, activity) to totals.
        :type frontier_total: dict[tuple[int, int], float]
        :param provenance_roots: Mapping of (year, activity) to root totals.
        :type provenance_roots: dict[tuple[int, int], dict[int, float]]
        :param y: Year to record.
        :type y: int
        :param a: Activity index to record.
        :type a: int
        :param x: Amount to add.
        :type x: float
        :param r: Root activity index, if any.
        :type r: int | None
        :param return_provenance: Whether to populate provenance data.
        :type return_provenance: bool
        """
        frontier_total[(int(y), int(a))] += float(x)
        if return_provenance and (r is not None):
            provenance_roots[(int(y), int(a))][int(r)] += float(x)

    @staticmethod
    def _record_direct_bio(
        direct_bio_total: dict[tuple[int, int], float],
        direct_bio_roots: dict[tuple[int, int], dict[int, float]],
        y: int,
        a: int,
        x: float,
        r: Optional[int],
        return_provenance: bool,
    ) -> None:
        """Record direct biosphere totals and optional provenance entries.

        :param direct_bio_total: Mapping of (year, activity) to totals.
        :type direct_bio_total: dict[tuple[int, int], float]
        :param direct_bio_roots: Mapping of (year, activity) to root totals.
        :type direct_bio_roots: dict[tuple[int, int], dict[int, float]]
        :param y: Year to record.
        :type y: int
        :param a: Activity index to record.
        :type a: int
        :param x: Amount to add.
        :type x: float
        :param r: Root activity index, if any.
        :type r: int | None
        :param return_provenance: Whether to populate provenance data.
        :type return_provenance: bool
        """
        direct_bio_total[(int(y), int(a))] += float(x)
        if return_provenance and (r is not None):
            direct_bio_roots[(int(y), int(a))][int(r)] += float(x)

    def _has_direct_biosphere(
        self, scenario_year: int, act: int, bio_cache: dict
    ) -> bool:
        """Check whether an activity has direct biosphere exchanges.

        :param scenario_year: Scenario year to query.
        :type scenario_year: int
        :param act: Activity index to query.
        :type act: int
        :param bio_cache: Cache mapping ``(year, activity)`` to bool.
        :type bio_cache: dict[tuple[int, int], bool]
        :returns: True if the activity has direct biosphere exchanges.
        :rtype: bool
        """
        label = str(scenario_year)
        if label in self.scenario_index and (self.B is not None):
            t = self.scenario_index[label]
            key = (scenario_year, act)
            if key in bio_cache:
                return bio_cache[key]
            has_direct_bio = self.B[t, act, :].nnz > 0
            bio_cache[key] = has_direct_bio
            return has_direct_bio
        return False

    def temporal_traversal(
        self,
        start_year: int,
        start_act_idx: int,
        amount: float = 1.0,
        max_depth: int = 3,
        min_amount: float = 1e-12,
        return_provenance: bool = False,
        show_progress: bool = False,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> tuple[dict, dict] | tuple[dict, dict, dict, dict]:
        """
        Traverse the temporal-technosphere graph starting from (start_year, start_act_idx).

        Progress:
          - Prefer an empirical total estimate based on max_depth (DEPTH_TOTALS).
          - If not available, fall back to a short warm-up branching estimate.
        """

        if debug:
            logger.info(
                "temporal_traversal start: start_year=%d start_act=%d amount=%g max_depth=%d min_amount=%g use_td=%s",
                start_year,
                start_act_idx,
                amount,
                max_depth,
                min_amount,
                use_temporal_distributions,
            )
            run_tag = f"{int(start_year)}:{int(start_act_idx)}:{int(max_depth)}:{time.time_ns()}"
            t0 = time.perf_counter()
            LOG_EVERY = 5000  # adjust as you like

        # ------------------------------------------------------------------
        # Progress bar setup
        # ------------------------------------------------------------------

        # ------------------------------------------------------------------
        # Progress estimation params (self-contained)
        # ------------------------------------------------------------------
        WARMUP_LIMIT = 1000
        BRANCHING_PERCENTILE = 95.0
        BRANCHING_SAFETY_FACTOR = 1.2
        EMPIRICAL_SAFETY_FACTOR = 1.05  # also used as a conservative headroom

        def estimate_total_from_branching(branching_samples: list[int]) -> int:
            """Estimate total nodes from observed branching samples.

            :param branching_samples: List of branching counts observed in warm-up.
            :type branching_samples: list[int]
            :returns: Estimated total node count.
            :rtype: int
            """
            if not branching_samples:
                return 1
            s = sorted(branching_samples)
            k = int((BRANCHING_PERCENTILE / 100.0) * (len(s) - 1))
            b = max(1.0, float(s[k]) * BRANCHING_SAFETY_FACTOR)
            # geometric series sum up to depth max_depth
            if abs(b - 1.0) < 1e-9:
                return max_depth + 1
            return int((b ** (max_depth + 1) - 1.0) / (b - 1.0))

        # Trackers needed by pbar helpers
        nodes_processed = 0
        branching_samples = []

        pbar = None
        total_est = self._estimate_total_from_depth(max_depth)
        try:
            if not show_progress:
                pbar = None
            elif total_est is None:
                # Indeterminate until warm-up can estimate
                pbar = tqdm(
                    total=None,
                    desc="Temporal traversal",
                    unit="node",
                    dynamic_ncols=True,
                )
            else:
                pbar = tqdm(
                    total=total_est,
                    desc="Temporal traversal",
                    unit="node",
                    dynamic_ncols=True,
                )
        except Exception:
            pbar = None

        # Track actual processed count and keep tqdm sane even if total was misestimated.
        # Policy:
        #  - If we exceed total, expand total so the bar never runs beyond 100%.
        #  - At the end, snap total to exactly nodes_processed so the bar finishes at 100%.
        def _pbar_step() -> None:
            """Advance the progress bar, expanding total if needed."""
            nonlocal nodes_processed, pbar
            nodes_processed += 1
            if pbar is None:
                return

            # If total is unknown, just advance.
            if pbar.total is None:
                pbar.update(1)
                return

            # If we're about to exceed the total, expand it with some headroom.
            if pbar.n + 1 > pbar.total:
                new_total = int(max(pbar.n + 1, pbar.total * 1.2, pbar.total + 100))
                pbar.total = new_total
                pbar.refresh()

            pbar.update(1)

        def _pbar_finalize() -> None:
            """Finalize the progress bar, snapping total to actual count."""
            nonlocal pbar
            if pbar is None:
                return
            try:
                # Snap to actual processed count so we end at 100%
                pbar.total = int(pbar.n)
                pbar.refresh()
                pbar.close()
            except Exception:
                pass

        # ------------------------------------------------------------------
        # Traversal state
        # ------------------------------------------------------------------
        queue = deque()
        queue.append((int(start_year), int(start_act_idx), float(amount), 0, (), None))

        frontier_total = defaultdict(float)  # (year, act) -> amt
        provenance_roots = defaultdict(
            lambda: defaultdict(float)
        )  # (year, act) -> {root_act: amt}
        direct_bio_total = defaultdict(
            float
        )  # (year, act) -> amt (nodes with direct biosphere we do NOT want to solve)
        direct_bio_roots = defaultdict(
            lambda: defaultdict(float)
        )  # same but by root, if provenance requested

        bio_cache: dict[tuple[int, int], bool] = {}

        while queue:
            year, act, amt, depth, path, root_act = queue.popleft()

            # ------------------------------------------------------------------
            # Ensure we always have a valid root for any non-root node.
            # If root_act is missing, recover it from the traversal path.
            # This prevents fallback attribution to the FU in lca.py.
            # ------------------------------------------------------------------
            if root_act is None and depth > 0:
                if path:
                    first = path[0]
                    # path stores ((year, act), ...) tuples
                    if isinstance(first, (tuple, list)) and len(first) >= 2:
                        root_act = int(first[1])
                    else:
                        root_act = int(act)
                else:
                    root_act = int(act)

            if amt == 0.0:
                continue

            _pbar_step()

            if debug and _log_every(LOG_EVERY, nodes_processed):
                logger.info(
                    "traversal_progress run=%s nodes=%d queue=%d frontier=%d direct_bio=%d elapsed_s=%.1f",
                    run_tag,
                    int(nodes_processed),
                    int(len(queue)),
                    int(len(frontier_total)),
                    int(len(direct_bio_total)),
                    float(time.perf_counter() - t0),
                )

            # Map to scenario year for "has direct biosphere" test (fast cutoff logic)
            scenario_year = self._map_year_to_scenario_year(year)
            has_direct_bio = self._has_direct_biosphere(scenario_year, act, bio_cache)

            # Stop expanding at max_depth
            if depth >= max_depth:
                self._record_frontier(
                    frontier_total,
                    provenance_roots,
                    year,
                    act,
                    amt,
                    root_act,
                    return_provenance,
                )
                continue

            # Expand this node (only if not cut)
            child_demands = self.expand_temporal_exchanges(
                year=year,
                act_idx=act,
                amount=amt,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
            )

            # Leaf: record it
            if not child_demands:
                self._record_frontier(
                    frontier_total,
                    provenance_roots,
                    year,
                    act,
                    amt,
                    root_act,
                    return_provenance,
                )
                continue

            # --------------------------------------------------------------
            # Warm-up: if tqdm started indeterminate (total=None),
            # estimate a total after WARMUP_LIMIT processed nodes using
            # observed branching and then set pbar.total.
            # --------------------------------------------------------------
            if show_progress and pbar is not None and pbar.total is None:
                # Branching sample = number of children edges we would enqueue
                # for this node (consistent with traversal).
                if nodes_processed <= WARMUP_LIMIT:
                    # Count children that would actually be enqueued
                    cnt = 0
                    if child_demands:
                        for _cy, _mapping in child_demands.items():
                            for _ca, _camt in _mapping.items():
                                cnt += 1
                    branching_samples.append(cnt)

                # Once warm-up complete, set a total estimate
                if nodes_processed == WARMUP_LIMIT:
                    est = estimate_total_from_branching(branching_samples)
                    # Keep a bit conservative so it doesn't finish early
                    est = int(max(est, pbar.n + 1) * EMPIRICAL_SAFETY_FACTOR)
                    pbar.total = est
                    pbar.refresh()

            # Enqueue children
            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    child_amt = float(child_amt)

                    # Preserve remainder monotonically: do not expand tiny contributions,
                    # keep them at the frontier so depth increases cannot reduce totals.
                    if abs(child_amt) < float(min_amount):
                        child_year = int(child_year)
                        child_act = int(child_act)
                        child_root = child_act if depth == 0 else root_act
                        self._record_frontier(
                            frontier_total,
                            provenance_roots,
                            child_year,
                            child_act,
                            child_amt,
                            child_root,
                            return_provenance,
                        )
                        continue

                    child_year = int(child_year)
                    child_act = int(child_act)

                    # Root propagation (Option A)
                    # - If parent is the start node (depth == 0): root becomes this child activity
                    # - Otherwise: propagate existing root_act
                    if depth == 0:
                        child_root = child_act
                        child_path = ((child_year, child_act),)
                    else:
                        child_root = root_act
                        child_path = path + ((child_year, child_act),)

                    queue.append(
                        (
                            child_year,
                            child_act,
                            child_amt,
                            depth + 1,
                            child_path,
                            child_root,
                        )
                    )

        _pbar_finalize()

        if debug:
            logger.info(
                "traversal_done run=%s nodes=%d frontier=%d direct_bio=%d elapsed_s=%.2f",
                run_tag,
                int(nodes_processed),
                int(len(frontier_total)),
                int(len(direct_bio_total)),
                float(time.perf_counter() - t0),
            )

        # Normalize provenance to plain dicts
        if return_provenance:
            provenance_roots = {k: dict(v) for k, v in provenance_roots.items()}
            direct_bio_roots = {k: dict(v) for k, v in direct_bio_roots.items()}
            return (
                dict(frontier_total),
                provenance_roots,
                dict(direct_bio_total),
                direct_bio_roots,
            )

        return dict(frontier_total), dict(direct_bio_total)

    def frontier_to_demand_vectors(self, frontier: dict) -> dict[int, np.ndarray]:
        """
        Convert a (year, activity) -> amount frontier into per-year demand vectors.

        Calendar years are preserved (no mapping to scenario years here).
        """
        if self.A is None:
            raise ValueError("Cannot build demand vectors: A is None")

        n_activities = int(self.A.shape[1])
        dtype = self.value_dtype

        f_by_year: dict[int, np.ndarray] = {}

        for key, amt in frontier.items():
            if not isinstance(key, tuple):
                raise ValueError(
                    f"Frontier key must be a tuple (year, act). Got {type(key)}: {key}"
                )

            if len(key) != 2:
                raise ValueError(
                    f"Frontier key must be (year, act). Got len={len(key)}: {key}"
                )

            year, act_idx = key
            y = int(year)
            if self._inventory_years is not None and self._inventory_years.size:
                y = max(
                    int(self._inventory_years[0]),
                    min(int(self._inventory_years[-1]), y),
                )
            else:
                # fallback: clamp to scenario range
                y = max(int(self.min_year), min(int(self.max_year), y))

            a = int(act_idx)

            if y not in f_by_year:
                f_by_year[y] = np.zeros(n_activities, dtype=dtype)

            f_by_year[y][a] += dtype(amt)

        return f_by_year

    def collect_traversal_edges(
        self,
        start_year: int,
        start_act_idx: int,
        amount: float = 1.0,
        max_depth: int = 3,
        min_amount: float = 1e-12,
    ) -> dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]]:
        """Traverse the temporal graph and record edges by depth.

        :param start_year: Start year for traversal.
        :type start_year: int
        :param start_act_idx: Start activity index.
        :type start_act_idx: int
        :param amount: Functional unit amount.
        :type amount: float
        :param max_depth: Maximum traversal depth.
        :type max_depth: int
        :param min_amount: Ignored (filtering disabled).
        :type min_amount: float
        :returns: Mapping of depth to edge amounts.
        :rtype: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]]
        """

        queue = deque()
        queue.append((int(start_year), int(start_act_idx), float(amount), 0))

        edges_by_depth: dict[
            int, dict[tuple[tuple[int, int], tuple[int, int]], float]
        ] = defaultdict(lambda: defaultdict(float))

        while queue:
            year, act, amt, depth = queue.popleft()

            if amt == 0.0:
                continue

            if depth >= max_depth:
                continue

            child_demands = self.expand_temporal_exchanges(
                year=int(year), act_idx=int(act), amount=float(amt)
            )
            if not child_demands:
                continue

            parent_node = (int(year), int(act))

            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    child_amt = float(child_amt)
                    if child_amt == 0.0:
                        continue

                    child_node = (int(child_year), int(child_act))
                    edges_by_depth[int(depth)][(parent_node, child_node)] += child_amt

                    queue.append(
                        (int(child_year), int(child_act), child_amt, int(depth) + 1)
                    )

        return {d: dict(edges) for d, edges in edges_by_depth.items()}

    def _apply_temporal_distribution_matrix_sourced_to_demand(
        self,
        *,
        year: int,
        act_idx: int,
        product_index: int,
        parent_amount: float,
        tex: TemporalExchange,
        demand: dict[int, dict[int, float]],
        debug: bool,
    ) -> None:
        if self.A is None:
            return

        if int(product_index) == int(act_idx):
            return

        td = TemporalDistribution(tex)
        offsets_and_weights = list(td.iter_offsets_and_weights(debug=debug))
        if not offsets_and_weights:
            return

        for offset, weight in offsets_and_weights:
            raw_year = int(year + int(offset))
            y_eff = int(self._map_year_to_scenario_year(raw_year))
            t_eff = self.scenario_index.get(str(y_eff))
            if t_eff is None:
                continue

            exchange_value = float(self.A[t_eff, int(act_idx), int(product_index)])
            if exchange_value == 0.0:
                continue

            # requirement magnitude at that pulse-year
            child_amount = self._child_amount(float(parent_amount), exchange_value)
            if child_amount == 0.0:
                continue

            # IMPORTANT: distribute mass using TD weights
            weighted_child_amount = float(child_amount) * float(weight)
            if weighted_child_amount == 0.0:
                continue

            if debug:
                logger.debug(
                    "expand_temporal_exchanges: matrix pulse raw_year=%d mapped_year=%d weight=%g",
                    int(raw_year),
                    int(y_eff),
                    float(weight),
                )
            self._add_demand_entry(
                demand, int(raw_year), int(product_index), weighted_child_amount
            )
