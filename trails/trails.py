# trails.py

from dataclasses import dataclass
import time
import os
from typing import Any, Dict, List, Optional, Callable
import importlib
import importlib.util
import warnings
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import sparse
import xarray as xr

from tqdm import tqdm

from .cache_interpolation import (
    load_cached_interpolation,
    save_cached_interpolation,
)
from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package,
)

from .temporal_distributions import TemporalDistribution, TemporalExchange
from .lca import lca_static
from .static_activity_scores import (
    _activity_score_potential,
    _ensure_static_activity_scores,
)
from .chunked_inventory import (
    ChunkedInventoryBuilder,
    DEFAULT_INVENTORY_MEMORY_BUDGET,
    estimate_flush_peak_bytes,
    estimate_materialization_peak_bytes,
    is_chunked_sparse,
)
from .factorized_inventory import FactorizedInventoryBuilder

import logging

logger = logging.getLogger(__name__)


DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF = 1e-4


class _DefaultAdaptiveRelativeScoreCutoff:
    def __repr__(self) -> str:
        """Represent the adaptive routing default in public signatures."""
        return repr(DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF)


_DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF = _DefaultAdaptiveRelativeScoreCutoff()


def _log_every(n: int, i: int) -> bool:
    """log every.

    :param n: Value for `n`.
    :type n: int
    :param i: Value for `i`.
    :type i: int
    :returns: Return value.
    :rtype: bool"""
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


@dataclass(frozen=True)
class _TechExpansionEntry:
    """Cached technosphere child expansion for one parent activity row."""

    child_act_idx: int
    scale: float
    amount_source: str
    offsets: tuple[int, ...] = ()
    weights: tuple[float, ...] = ()


class Trails:
    """Wrapper around time-indexed technosphere and biosphere matrices.

    ``Trails`` loads 3D sparse matrices from a Frictionless data package, can
    interpolate scenario slices to annual resolution, and stores optional
    default LCIA configuration. The recommended workflow is
    ``temporal_routing()``, ``lci()``, then one or more ``lcia()`` calls.
    Constructor ``methods`` are optional characterization defaults; routing
    screening methods are configured independently.

    Dimensions: ``A`` is ``(scenario, activity, product)`` and ``B`` is
    ``(scenario, activity, biosphere_flow)``.
    """

    def __init__(
        self,
        package: Any,
        interpolate_annual: bool = True,
        cache_interpolation: bool = True,
        interpolation_start_year_offset: int = -1,
        interpolation_end_year_offset: int = 1,
        value_dtype: np.dtype = np.float32,
        index_dtype: np.dtype = np.int32,
        methods: list[Any] | None = None,
        method_backend: str = "auto",
        edges_methods: list[Any] | None = None,
        ei_version: str = "3.11",
        debug: bool = False,
    ) -> None:
        """Initialize a Trails data package wrapper.

        :param package: Frictionless data package or compatible package object.
        :type package: Any
        :param interpolate_annual: Interpolate scenario matrices to annual
            resolution.
        :type interpolate_annual: bool
        :param cache_interpolation: Load/write annual interpolation caches.
        :type cache_interpolation: bool
        :param interpolation_start_year_offset: Offset from min inventory year.
        :type interpolation_start_year_offset: int
        :param interpolation_end_year_offset: Offset from max inventory year.
        :type interpolation_end_year_offset: int
        :param value_dtype: Sparse matrix value dtype.
        :type value_dtype: np.dtype
        :param index_dtype: Sparse matrix coordinate dtype.
        :type index_dtype: np.dtype
        :param methods: Optional default regular or EDGES methods for
            ``lcia()``. Methods can instead be supplied to each ``lcia()`` call.
        :type methods: list[Any] | None
        :param method_backend: Default characterization backend. ``"auto"``
            infers mapping-based EDGES methods and treats string names as
            regular methods.
        :type method_backend: str
        :param edges_methods: Deprecated alias for ``methods`` with
            ``method_backend="edges"``.
        :type edges_methods: list[Any] | None
        :param ei_version: Default ecoinvent LCIA data version for regular
            methods and adaptive routing.
        :type ei_version: str
        :param debug: Enable diagnostic logging.
        :type debug: bool"""
        self.package = package
        self.value_dtype = value_dtype
        self.index_dtype = index_dtype
        self.debug = debug
        if method_backend not in {"auto", "regular", "edges"}:
            raise ValueError(
                "method_backend must be one of {'auto', 'regular', 'edges'}"
            )
        if methods and edges_methods:
            raise ValueError("methods and edges_methods are mutually exclusive.")
        if edges_methods:
            warnings.warn(
                "Trails(..., edges_methods=...) is deprecated; pass the EDGES "
                "method through methods=... with method_backend='edges'.",
                FutureWarning,
                stacklevel=2,
            )
            methods = list(edges_methods)
            method_backend = "edges"
        self.methods = list(methods) if methods else None
        self.method_backend = str(method_backend)
        self.edges_methods = list(edges_methods) if edges_methods else None
        self.ei_version = str(ei_version)
        self.default_methods = self.methods
        self.default_edges_methods = self.edges_methods
        self.default_method_backend = self.method_backend
        self.default_ei_version = self.ei_version
        self.interpolation_start_year_offset = int(interpolation_start_year_offset)
        self.interpolation_end_year_offset = int(interpolation_end_year_offset)

        # If a zip archive is provided, Frictionless unpacks it to a temp basepath.
        pkg_path = getattr(self.package, "path", None)
        pkg_base = getattr(self.package, "basepath", None)
        if pkg_path and str(pkg_path).lower().endswith(".zip"):
            base_str = f" -> {pkg_base}" if pkg_base else ""
            print(f"Data package unarchived from: {pkg_path}{base_str}")

        self.scenario_labels: List[str] = []
        self.scenario_index: Dict[str, int] = {}

        self.A: Optional[sparse.COO] = None
        self.B: Optional[sparse.COO] = None
        self.inventory: Optional[xr.DataArray] = None
        self.characterized_inventory: Optional[xr.DataArray] = None
        self.static_score: Optional[float | list[float]] = None
        self._inventory_years: Optional[np.ndarray] = None
        self._inventory_year_index: dict[int, int] = {}
        self._inventory_coords: Optional[list[list[np.ndarray]]] = None
        self._inventory_data: Optional[list[np.ndarray]] = None
        self.provenance: Optional[dict] = None

        cache_loaded = False
        if interpolate_annual and cache_interpolation:
            (
                A_cached,
                B_cached,
                labels,
                template_labels_cached,
                temporal_tech,
                temporal_bio,
                indices_cached,
                cache_dir,
            ) = load_cached_interpolation(
                self.package,
                value_dtype=str(self.value_dtype),
                index_dtype=str(self.index_dtype),
                interpolation_start_year_offset=self.interpolation_start_year_offset,
                interpolation_end_year_offset=self.interpolation_end_year_offset,
            )
            if (
                A_cached is not None
                and B_cached is not None
                and labels
                and temporal_tech is not None
                and temporal_bio is not None
                and indices_cached is not None
            ):
                self.A = A_cached
                self.B = B_cached
                self.scenario_labels = list(labels)
                self.scenario_index = {
                    lbl: i for i, lbl in enumerate(self.scenario_labels)
                }
                self.template_labels = (
                    list(template_labels_cached)
                    if template_labels_cached
                    else list(labels)
                )
                self.temporal_technosphere_exchanges = temporal_tech
                self.temporal_biosphere_exchanges = temporal_bio
                self.activity_indices = indices_cached.get("activity_indices", {})
                self.biosphere_indices = indices_cached.get("biosphere_indices", {})
                cache_loaded = True
                print(f"Loaded interpolated matrices from cache: {cache_dir}")
                if debug:
                    logger.info(
                        "Trails init: loaded interpolated matrices from cache %s",
                        str(cache_dir),
                    )

        if not cache_loaded:
            print("Loading matrices from data package          [1/4]")
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

        if not hasattr(self, "template_labels"):
            self.template_labels = list(self.scenario_labels)
        self.template_years_int = np.array(
            [int(lbl) for lbl in self.template_labels], dtype=int
        )

        self.years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
        self.min_year = int(self.years_int.min())
        self.max_year = int(self.years_int.max())

        if not cache_loaded:
            # Load indices/metadata
            print("Loading indices from data package           [2/4]")
            (
                self.activity_indices,
                self.biosphere_indices,
            ) = load_indices_from_package(self.package)

        # Optional temporal interpolation to annual resolution
        if interpolate_annual and self.scenario_labels and not cache_loaded:
            print("Interpolating matrices to annual resolution [3/4]")
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
                start_year_offset=self.interpolation_start_year_offset,
                end_year_offset=self.interpolation_end_year_offset,
                debug=debug,
            )

            self.years_int = np.array(
                [int(lbl) for lbl in self.scenario_labels], dtype=int
            )
            self.min_year = int(self.years_int.min())
            self.max_year = int(self.years_int.max())

            if cache_interpolation:
                try:
                    print("Building cache                              [4/4]")
                    cache_dir = save_cached_interpolation(
                        self.package,
                        value_dtype=str(self.value_dtype),
                        index_dtype=str(self.index_dtype),
                        interpolation_start_year_offset=(
                            self.interpolation_start_year_offset
                        ),
                        interpolation_end_year_offset=self.interpolation_end_year_offset,
                        A=self.A,
                        B=self.B,
                        scenario_labels=self.scenario_labels,
                        template_labels=self.template_labels,
                        temporal_technosphere_exchanges=self.temporal_technosphere_exchanges,
                        temporal_biosphere_exchanges=self.temporal_biosphere_exchanges,
                        activity_indices=self.activity_indices,
                        biosphere_indices=self.biosphere_indices,
                    )
                    print(f"Data package cached at: {cache_dir}")
                    if debug:
                        logger.info(
                            "Trails init: cached interpolated matrices at %s",
                            str(cache_dir),
                        )
                except Exception:
                    pass

        self.scores: Optional[xr.DataArray] = (
            None  # dims: (activity, year) or (activity, year, root activity) or (+method)
        )
        self._score_years: Optional[np.ndarray] = None
        self._score_year_index: dict[int, int] = {}
        self.graph = None
        self._routing_attribute_to_roots: Optional[bool] = None
        self._routing_params: Optional[dict[str, object]] = None
        self._td_offsets_cache: dict[tuple, list[tuple[int, float]]] = {}
        self._tech_td_cache: dict[tuple[int, int, int], Optional[TemporalExchange]] = {}
        self._A_row_cache: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
        self._production_amount_cache: dict[tuple[int, int], float] = {}
        self._production_amount_vector_cache: dict[int, np.ndarray] = {}
        self._direct_bio_cache_by_year: dict[int, np.ndarray] = {}
        self._tech_td_expanded_cache: dict[
            tuple[int, int, int],
            tuple[Optional[TemporalExchange], list[tuple[int, float]]],
        ] = {}
        self._tech_expansion_template_cache: dict[
            tuple[int, int, int, bool], tuple[_TechExpansionEntry, ...]
        ] = {}
        self._scenario_year_map_cache: dict[int, int] = {}
        self._template_year_map_cache: dict[int, int] = {}
        self._static_activity_score_cache: dict[str, object] = {}
        self._static_activity_score_fingerprint: tuple[tuple, str] | None = None
        self._inventory_backend_requested = "coo"
        self._inventory_memory_budget = DEFAULT_INVENTORY_MEMORY_BUDGET
        self._inventory_store: str | Path | None = None
        self._inventory_builder: (
            ChunkedInventoryBuilder | FactorizedInventoryBuilder | None
        ) = None
        self._inventory_builders: list[
            ChunkedInventoryBuilder | FactorizedInventoryBuilder
        ] = []
        self.inventory_diagnostics: dict[str, object] = {}
        self.lca_diagnostics: dict[str, object] = {}
        self.lci_diagnostics: dict[str, object] = {}
        self.lcia_diagnostics: dict[str, object] = {}
        self.lcia_results: dict[str, dict[str, Any]] = {}
        self.current_lcia_result: str | None = None

    def configure_inventory_storage(
        self,
        *,
        backend: str = "auto",
        memory_budget: int = DEFAULT_INVENTORY_MEMORY_BUDGET,
        store: str | Path | None = None,
    ) -> None:
        """Configure storage used by the next ``reset_inventory`` call."""
        if backend not in {"auto", "coo", "chunked", "factorized"}:
            raise ValueError(
                "inventory_backend must be one of "
                "{'auto', 'coo', 'chunked', 'factorized'}"
            )
        budget = int(memory_budget)
        if budget <= 0:
            raise ValueError("inventory_memory_budget must be a positive integer")
        self._inventory_backend_requested = str(backend)
        self._inventory_memory_budget = budget
        self._inventory_store = store

    def _new_chunked_inventory_builder(self) -> ChunkedInventoryBuilder:
        if self.A is None or self.B is None or self._inventory_years is None:
            raise RuntimeError("Inventory dimensions are not initialized")
        builder = ChunkedInventoryBuilder(
            n_activities=int(self.A.shape[1]),
            n_flows=int(self.B.shape[2]),
            n_years=int(self._inventory_years.size),
            has_root=bool(self._inventory_has_root),
            value_dtype=self.value_dtype,
            memory_budget=int(self._inventory_memory_budget),
            store=self._inventory_store,
        )
        self._inventory_builders.append(builder)
        return builder

    def _new_factorized_inventory_builder(self) -> FactorizedInventoryBuilder:
        if self.A is None or self.B is None or self._inventory_years is None:
            raise RuntimeError("Inventory dimensions are not initialized")
        builder = FactorizedInventoryBuilder(
            n_activities=int(self.A.shape[1]),
            n_flows=int(self.B.shape[2]),
            n_years=int(self._inventory_years.size),
            has_root=bool(self._inventory_has_root),
            value_dtype=self.value_dtype,
            memory_budget=int(self._inventory_memory_budget),
            store=self._inventory_store,
        )
        self._inventory_builders.append(builder)
        return builder

    def _promote_inventory_to_chunked(self) -> None:
        if self._inventory_builder is not None:
            return
        builder = self._new_chunked_inventory_builder()
        for keys, values in zip(self._inv_key_parts, self._inv_value_parts):
            builder.append_linear_global(keys, values)
        self._inv_key_parts.clear()
        self._inv_value_parts.clear()
        self._inventory_builder = builder

    def _maybe_promote_inventory_to_chunked(self) -> None:
        if self._inventory_builder is not None:
            return
        if getattr(self, "_inventory_backend_requested", "coo") != "auto":
            return
        entries = sum(int(part.size) for part in self._inv_key_parts)
        ndim = 4 if bool(getattr(self, "_inventory_has_root", False)) else 3
        predicted = max(
            estimate_flush_peak_bytes(entries, value_dtype=self.value_dtype),
            estimate_materialization_peak_bytes(
                entries,
                ndim=ndim,
                value_dtype=self.value_dtype,
            ),
        )
        if predicted > int(self._inventory_memory_budget):
            self._promote_inventory_to_chunked()

    def reset_scores(
        self,
        *,
        attribute_to_roots: bool = False,
        methods: list[str] | None = None,
    ) -> None:
        """Reset scores.

        :param attribute_to_roots: Value for `attribute_to_roots`.
        :type attribute_to_roots: bool"""
        min_offset, max_offset = self._inventory_offset_bounds()
        # Extend the inventory time axis to allow biosphere emissions to spread/decay
        # beyond the matrix years (e.g., for long-lived gases).
        tail_years = 500
        years = np.arange(
            int(self.min_year) + int(min_offset),
            int(self.max_year) + int(max_offset) + int(tail_years) + 1,
            dtype=int,
        )

        self._score_years = years
        self._score_year_index = {int(y): int(i) for i, y in enumerate(years)}
        self._scores_has_root = bool(attribute_to_roots)
        self._score_methods = list(methods) if methods else None

        self._score_chunk_act = []
        self._score_chunk_year = []
        self._score_chunk_root = []
        self._score_chunk_method = []
        self._score_chunk_value = []

        # Bulk score builder (vectorized appends)
        self._score_bulk_act = []
        self._score_bulk_year = []
        self._score_bulk_root = []
        self._score_bulk_method = []
        self._score_bulk_value = []

        self.scores = None

    def _clamp_year_to_inventory(self, year: int) -> int:
        """clamp year to inventory.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: int
        :raises RuntimeError: If an error occurs."""
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
        """clamp year to scores.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: int
        :raises RuntimeError: If an error occurs."""
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

    def _get_debug_flow_filters(self, *, debug: bool) -> dict | None:
        """get debug flow filters.

        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: dict | None"""
        if not (debug or self.debug):
            return None
        if getattr(self, "_debug_flow_filters", None) is not None:
            return self._debug_flow_filters
        if not hasattr(self, "_debug_flow_filters_logged"):
            self._debug_flow_filters_logged = False

        flow_id = os.getenv("TRAILS_DEBUG_FLOW_ID")
        year = os.getenv("TRAILS_DEBUG_YEAR")
        act = os.getenv("TRAILS_DEBUG_ACTIVITY")
        max_pulses = os.getenv("TRAILS_DEBUG_MAX_PULSES")
        max_matches = os.getenv("TRAILS_DEBUG_MAX_MATCHES")

        self._debug_flow_filters = {
            "flow_id": int(flow_id) if flow_id not in (None, "") else None,
            "year": int(year) if year not in (None, "") else None,
            "act": int(act) if act not in (None, "") else None,
            "max_pulses": int(max_pulses) if max_pulses not in (None, "") else 12,
            "max_matches": int(max_matches) if max_matches not in (None, "") else 50,
            "matches": 0,
        }

        if not self._debug_flow_filters_logged:
            logger.debug(
                "debug_flow_filters: %s",
                {
                    "flow_id": self._debug_flow_filters["flow_id"],
                    "year": self._debug_flow_filters["year"],
                    "act": self._debug_flow_filters["act"],
                    "max_pulses": self._debug_flow_filters["max_pulses"],
                    "max_matches": self._debug_flow_filters["max_matches"],
                },
            )
            self._debug_flow_filters_logged = True

        return self._debug_flow_filters

    def _append_scores_from_yearidx_map(
        self,
        act_idx: int,
        yearidx_to_value: dict[int, float],
        *,
        root_activity: int | None = None,
        method_idx: int | None = None,
    ) -> None:
        """append scores from yearidx map.

        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param yearidx_to_value: Value for `yearidx_to_value`.
        :type yearidx_to_value: dict[int, float]
        :param root_activity: Value for `root_activity`.
        :type root_activity: int | None"""
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

        methods = getattr(self, "_score_methods", None)
        if methods is not None:
            if method_idx is None:
                raise ValueError(
                    "method_idx must be provided when scores have method dimension."
                )
            self._score_chunk_method.extend([int(method_idx)] * n)

        if getattr(self, "_scores_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            self._score_chunk_root.extend([int(root_activity)] * n)

    def reset_inventory(
        self,
        *,
        attribute_to_roots: bool = False,
        reset_scores: bool = True,
        score_methods: list[str] | None = None,
    ) -> None:
        """Reset inventory.

        :param attribute_to_roots: Value for `attribute_to_roots`.
        :type attribute_to_roots: bool
        :param reset_scores: Value for `reset_scores`.
        :type reset_scores: bool
        :raises RuntimeError: If an error occurs."""
        if self.min_year is None or self.max_year is None:
            if not getattr(self, "scenario_labels", None):
                raise RuntimeError(
                    "Trails scenario years not initialized; cannot reset inventory."
                )
            years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
            if years_int.size == 0:
                raise RuntimeError(
                    "Trails scenario years not initialized; cannot reset inventory."
                )
            self.min_year = int(years_int.min())
            self.max_year = int(years_int.max())
        min_offset, max_offset = self._inventory_offset_bounds()
        tail_years = 500
        years = np.arange(
            int(self.min_year) + int(min_offset),
            int(self.max_year) + int(max_offset) + int(tail_years) + 1,
            dtype=int,
        )

        self._inventory_years = years
        self._inventory_year_index = {int(y): int(i) for i, y in enumerate(years)}
        self._inventory_has_root = bool(attribute_to_roots)

        # Linearized inventory builders (fast append + finalize dedup)
        self._inv_key_parts = []
        self._inv_value_parts = []
        self._inventory_builder = None
        backend = getattr(self, "_inventory_backend_requested", "coo")
        if backend == "chunked":
            self._inventory_builder = self._new_chunked_inventory_builder()
        elif backend == "factorized":
            self._inventory_builder = self._new_factorized_inventory_builder()

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
            self._score_methods = list(score_methods) if score_methods else None

            self._score_chunk_act = []
            self._score_chunk_year = []
            self._score_chunk_root = []
            self._score_chunk_method = []
            self._score_chunk_value = []

            self._score_bulk_act = []
            self._score_bulk_year = []
            self._score_bulk_root = []
            self._score_bulk_method = []
            self._score_bulk_value = []

            self.scores = None

    def import_excel_inventory(
        self,
        path: str | Path | list[str | Path],
        *,
        year: int | None = None,
        scenario_label: str | None = None,
        cache_import: bool = False,
    ) -> dict[str, int]:
        """Import excel inventory.

        :param path: Value for `path`.
        :type path: str | Path | list[str | Path]
        :param year: Value for `year`.
        :type year: int | None
        :param scenario_label: Value for `scenario_label`.
        :type scenario_label: str | None
        :param cache_import: Value for `cache_import`.
        :type cache_import: bool
        :returns: Return value.
        :rtype: dict[str, int]"""
        from .importer import import_excel_inventory

        self._invalidate_calculation_results(close_inventory=True)
        self.graph = None

        return import_excel_inventory(
            self,
            path,
            year=year,
            scenario_label=scenario_label,
            cache_import=cache_import,
        )

    def _invalidate_calculation_results(self, *, close_inventory: bool) -> None:
        """Invalidate inventory and characterization state after model changes."""
        if close_inventory:
            for builder in getattr(self, "_inventory_builders", []):
                builder.close()
            self._inventory_builders = []
            self._inventory_builder = None
        self.inventory = None
        self.characterized_inventory = None
        self.scores = None
        self.lci_diagnostics = {}
        self.lcia_diagnostics = {}
        self.lcia_results = {}
        self.current_lcia_result = None

    def _append_score_entry(
        self,
        act_idx: int,
        year: int,
        value: float,
        *,
        root_activity: int | None = None,
        method_idx: int | None = None,
    ) -> None:
        """append score entry.

        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param year: Value for `year`.
        :type year: int
        :param value: Value for `value`.
        :type value: float
        :param root_activity: Value for `root_activity`.
        :type root_activity: int | None
        :raises RuntimeError: If an error occurs."""
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

        methods = getattr(self, "_score_methods", None)
        if methods is not None:
            if method_idx is None:
                raise ValueError(
                    "method_idx must be provided when scores have method dimension."
                )
            self._score_chunk_method.append(int(method_idx))

        if getattr(self, "_scores_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            self._score_chunk_root.append(int(root_activity))

    def finalize_scores(self) -> xr.DataArray:
        """Finalize scores.

        :returns: Return value.
        :rtype: xr.DataArray
        :raises RuntimeError: If an error occurs."""
        years = self._score_years
        if years is None:
            raise RuntimeError("Scores years not initialized. Call reset_inventory().")
        if self.A is None:
            raise RuntimeError("A is None")

        n_activities = int(self.A.shape[1])
        has_root = bool(self._scores_has_root)
        methods = getattr(self, "_score_methods", None)
        has_method = methods is not None

        has_any = bool(self._score_chunk_act) or bool(
            getattr(self, "_score_bulk_act", [])
        )
        if not has_any:
            if has_root and has_method:
                arr = sparse.zeros(
                    (len(methods), n_activities, len(years), n_activities),
                    dtype=self.value_dtype,
                )
                self.scores = xr.DataArray(
                    arr,
                    dims=("method", "activity", "year", "root activity"),
                    coords={
                        "method": np.asarray(methods, dtype=object),
                        "activity": np.arange(n_activities, dtype=int),
                        "year": years,
                        "root activity": np.arange(n_activities, dtype=int),
                    },
                )
            elif has_root:
                arr = sparse.zeros(
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
            elif has_method:
                arr = sparse.zeros(
                    (len(methods), n_activities, len(years)), dtype=self.value_dtype
                )
                self.scores = xr.DataArray(
                    arr,
                    dims=("method", "activity", "year"),
                    coords={
                        "method": np.asarray(methods, dtype=object),
                        "activity": np.arange(n_activities, dtype=int),
                        "year": years,
                    },
                )
            else:
                arr = sparse.zeros((n_activities, len(years)), dtype=self.value_dtype)
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
        method_parts = []

        # Chunk parts
        if self._score_chunk_act:
            act_c = np.asarray(self._score_chunk_act, dtype=np.int64)
            yr_c = np.asarray(self._score_chunk_year, dtype=np.int64)
            data_c = np.asarray(self._score_chunk_value, dtype=self.value_dtype)
            coords_parts.append((act_c, yr_c))
            data_parts.append(data_c)
            if has_method:
                method_c = np.asarray(self._score_chunk_method, dtype=np.int64)
                method_parts.append(method_c)
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
            if has_method:
                method_b = np.concatenate(self._score_bulk_method).astype(
                    np.int64, copy=False
                )
                method_parts.append(method_b)
            if has_root:
                root_b = np.concatenate(self._score_bulk_root).astype(
                    np.int64, copy=False
                )
                root_parts.append(root_b)

        act = np.concatenate([p[0] for p in coords_parts])
        yr = np.concatenate([p[1] for p in coords_parts])
        data = np.concatenate(data_parts).astype(self.value_dtype, copy=False)

        if has_root and has_method:
            method_idx = np.concatenate(method_parts).astype(np.int64, copy=False)
            root = np.concatenate(root_parts).astype(np.int64, copy=False)
            coords = np.vstack([method_idx, act, yr, root])
            arr = sparse.COO(
                coords,
                data,
                shape=(len(methods), n_activities, len(years), n_activities),
            )
            self.scores = xr.DataArray(
                arr,
                dims=("method", "activity", "year", "root activity"),
                coords={
                    "method": np.asarray(methods, dtype=object),
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
                    "root activity": np.arange(n_activities, dtype=int),
                },
            )
        elif has_root:
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
        elif has_method:
            method_idx = np.concatenate(method_parts).astype(np.int64, copy=False)
            coords = np.vstack([method_idx, act, yr])
            arr = sparse.COO(
                coords, data, shape=(len(methods), n_activities, len(years))
            )
            self.scores = xr.DataArray(
                arr,
                dims=("method", "activity", "year"),
                coords={
                    "method": np.asarray(methods, dtype=object),
                    "activity": np.arange(n_activities, dtype=int),
                    "year": years,
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
        """inventory offset bounds.

        :returns: Return value.
        :rtype: tuple[int, int]"""
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

    def _biosphere_offset_bounds(self) -> tuple[int, int]:
        """biosphere offset bounds.

        :returns: Return value.
        :rtype: tuple[int, int]"""
        min_offset = 0
        max_offset = 0
        exchanges = self.temporal_biosphere_exchanges
        if exchanges:
            for tex in exchanges.values():
                offset_min = getattr(tex, "offset_min", None)
                offset_max = getattr(tex, "offset_max", None)
                if offset_min is not None:
                    min_offset = min(min_offset, int(offset_min))
                if offset_max is not None:
                    max_offset = max(max_offset, int(offset_max))
        return min_offset, max_offset

    def _linearize_inventory_entries(
        self,
        act_idx: np.ndarray,
        flow_idx: np.ndarray,
        year_idx: np.ndarray,
        *,
        root_idx: np.ndarray | None = None,
    ) -> np.ndarray:
        """Linearize inventory coordinates into 1D integer keys."""
        if self.A is None or self.B is None or self._inventory_years is None:
            raise RuntimeError(
                "Inventory dimensions not initialized. Call reset_inventory() first."
            )

        n_acts = int(self.A.shape[1])
        n_flows = int(self.B.shape[2])
        n_years = int(self._inventory_years.size)

        acts = np.asarray(act_idx, dtype=np.int64)
        flows = np.asarray(flow_idx, dtype=np.int64)
        years = np.asarray(year_idx, dtype=np.int64)
        keys = ((acts * n_flows) + flows) * n_years + years
        if root_idx is not None:
            roots = np.asarray(root_idx, dtype=np.int64)
            keys = keys * n_acts + roots
        return keys

    def _append_inventory_entries(
        self,
        act_idx: int,
        year: int,
        flows: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: int | None = None,
    ) -> None:
        """append inventory entries.

        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param year: Value for `year`.
        :type year: int
        :param flows: Value for `flows`.
        :type flows: np.ndarray
        :param values: Value for `values`.
        :type values: np.ndarray
        :param root_activity: Value for `root_activity`.
        :type root_activity: int | None
        :raises RuntimeError: If an error occurs.
        :raises ValueError: If an error occurs."""
        # Inventory builders must exist
        if not hasattr(self, "_inv_key_parts"):
            raise RuntimeError(
                "Inventory builders not initialized. Call reset_inventory() first."
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

        acts = np.full(n, int(act_idx), dtype=np.int64)
        years = np.full(n, int(year_idx), dtype=np.int64)
        if getattr(self, "_inventory_has_root", False):
            if root_activity is None:
                root_activity = int(act_idx)
            roots = np.full(n, int(root_activity), dtype=np.int64)
        else:
            roots = None

        keys = self._linearize_inventory_entries(
            acts,
            flows_i64,
            years,
            root_idx=roots,
        )
        if self._inventory_builder is not None:
            self._inventory_builder.append(
                acts,
                flows_i64,
                years,
                vals_out,
                roots=roots,
            )
        else:
            self._inv_key_parts.append(keys)
            self._inv_value_parts.append(vals_out)
            self._maybe_promote_inventory_to_chunked()

    def _append_scores_bulk(
        self,
        act_idx: np.ndarray,
        year_idx: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: np.ndarray | None = None,
        method_idx: int | np.ndarray | None = None,
    ) -> None:
        """append scores bulk.

        :param act_idx: Value for `act_idx`.
        :type act_idx: np.ndarray
        :param year_idx: Value for `year_idx`.
        :type year_idx: np.ndarray
        :param values: Value for `values`.
        :type values: np.ndarray
        :param root_activity: Value for `root_activity`.
        :type root_activity: np.ndarray | None
        :raises RuntimeError: If an error occurs.
        :raises ValueError: If an error occurs."""
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

        methods = getattr(self, "_score_methods", None)
        if methods is not None:
            if method_idx is None:
                raise ValueError(
                    "method_idx must be provided when scores have method dimension."
                )
            if np.isscalar(method_idx):
                method_arr = np.full(a.shape, int(method_idx), dtype=np.int64)
            else:
                method_arr = np.asarray(method_idx)
                if method_arr.shape != a.shape:
                    raise ValueError("method_idx must match act_idx shape")
                method_arr = method_arr[m].astype(np.int64, copy=False)
            self._score_bulk_method.append(method_arr)

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
        """append inventory entries bulk.

        :param act_idx: Value for `act_idx`.
        :type act_idx: np.ndarray
        :param year: Value for `year`.
        :type year: int | np.ndarray
        :param flows: Value for `flows`.
        :type flows: np.ndarray
        :param values: Value for `values`.
        :type values: np.ndarray
        :param root_activity: Value for `root_activity`.
        :type root_activity: int | np.ndarray | None
        :raises RuntimeError: If an error occurs.
        :raises ValueError: If an error occurs."""
        if not hasattr(self, "_inv_key_parts"):
            raise RuntimeError(
                "Inventory builders not initialized. Call reset_inventory() first."
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

        years_axis = self._inventory_years
        if years_axis is None or years_axis.size == 0:
            raise RuntimeError(
                "Inventory years not initialized. Call reset_inventory() first."
            )
        y0 = int(years_axis[0])
        y1 = int(years_axis[-1])

        if isinstance(year, np.ndarray):
            years_arr = np.asarray(year)
            if years_arr.shape[0] != flows_arr.shape[0]:
                raise ValueError(
                    "year array must match act/flow/value length for bulk append."
                )
            years_clamped = np.clip(years_arr.astype(np.int64, copy=False), y0, y1)
            year_idx = years_clamped - y0
        else:
            y = self._clamp_year_to_inventory(int(year))
            year_idx = np.full(flows_arr.shape[0], int(y - y0), dtype=np.int64)

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

        if getattr(self, "_inventory_has_root", False):
            if root_activity is None:
                root_arr = acts_out
            else:
                root_arr = np.asarray(root_activity)
                if root_arr.shape == ():
                    root_arr = np.full_like(acts_out, int(root_arr))
                elif root_arr.shape[0] != acts_arr.shape[0]:
                    raise ValueError(
                        "root_activity array must match act/flow/value length for bulk append."
                    )
                else:
                    root_arr = root_arr[mask]
            root_out = np.asarray(root_arr, dtype=np.int64)
        else:
            root_out = None

        if self._inventory_builder is not None:
            self._inventory_builder.append(
                acts_out,
                flows_out,
                years_out,
                vals_out,
                roots=root_out,
            )
        else:
            keys = self._linearize_inventory_entries(
                acts_out,
                flows_out,
                years_out,
                root_idx=root_out,
            )
            self._inv_key_parts.append(keys)
            self._inv_value_parts.append(vals_out)
            self._maybe_promote_inventory_to_chunked()

    def finalize_inventory(self, *, show_progress: bool = False) -> xr.DataArray:
        """Finalize inventory.

        :returns: Return value.
        :rtype: xr.DataArray
        :raises RuntimeError: If an error occurs.
        :raises ValueError: If an error occurs."""
        if self.A is None or self.B is None:
            raise ValueError("Cannot finalize inventory: A or B is None.")

        years = self._inventory_years
        if years is None:
            raise RuntimeError(
                "Inventory years not initialized. Call reset_inventory()."
            )

        if not hasattr(self, "_inv_key_parts"):
            raise RuntimeError(
                "Inventory builders not initialized. Call reset_inventory()."
            )

        n_activities = int(self.A.shape[1])
        n_flows = int(self.B.shape[2])
        has_root = bool(getattr(self, "_inventory_has_root", False))

        if self._inventory_builder is not None:
            inv = self._inventory_builder.finalize(show_progress=show_progress)
            self.inventory_diagnostics = self._inventory_builder.diagnostics()
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
            self._inv_key_parts = []
            self._inv_value_parts = []
            return self.inventory

        eager_entries = sum(int(part.size) for part in self._inv_key_parts)
        eager_peak = max(
            estimate_flush_peak_bytes(eager_entries, value_dtype=self.value_dtype),
            estimate_materialization_peak_bytes(
                eager_entries,
                ndim=4 if has_root else 3,
                value_dtype=self.value_dtype,
            ),
        )
        if eager_peak > int(getattr(self, "_inventory_memory_budget", 0)):
            backend = getattr(self, "_inventory_backend_requested", "coo")
            raise MemoryError(
                f"Eager inventory finalization requires an estimated "
                f"{eager_peak / 2**30:.2f} GiB, above the configured "
                f"{self._inventory_memory_budget / 2**30:.2f} GiB budget "
                f"for inventory_backend={backend!r}. Use backend='auto' or "
                "'chunked', or explicitly raise the memory budget."
            )

        if self._inv_key_parts:
            if len(self._inv_key_parts) != len(self._inv_value_parts):
                raise RuntimeError(
                    "Inventory key/value builders are inconsistent. "
                    "Call reset_inventory() and rerun lci()."
                )

            keys = np.concatenate(self._inv_key_parts).astype(np.int64, copy=False)
            data = np.concatenate(self._inv_value_parts).astype(
                self.value_dtype, copy=False
            )

            # Free builder references once contiguous arrays are materialized.
            self._inv_key_parts.clear()
            self._inv_value_parts.clear()

            if keys.size:
                # Reassign sorted arrays to drop pre-sort buffers earlier and
                # keep the fast argsort path.
                order = np.argsort(keys, kind="quicksort")
                keys = keys[order]
                data = data[order]
                del order

                first = np.empty(keys.size, dtype=bool)
                first[0] = True
                first[1:] = keys[1:] != keys[:-1]
                group_starts = np.flatnonzero(first)
                data_agg = np.add.reduceat(data, group_starts).astype(
                    self.value_dtype, copy=False
                )
                keys_agg = keys[group_starts].copy()

                keep = data_agg != 0.0
                keys_agg = keys_agg[keep]
                data_agg = data_agg[keep]

                del first
                del group_starts
            else:
                keys_agg = np.empty(0, dtype=np.int64)
                data_agg = np.empty(0, dtype=self.value_dtype)

            n_years = int(len(years))
            if has_root:
                shape = (n_activities, n_flows, n_years, n_activities)
            else:
                shape = (n_activities, n_flows, n_years)

            # Keep inventory coordinates int64 to avoid uint64 reduce index
            # promotion in sparse reductions on large reshaped arrays.
            coord_dtype = np.int64

            if keys_agg.size:
                # Decode linearized keys in-place to avoid allocating several
                # large temporary coordinate arrays.
                if has_root:
                    coords = np.empty((4, keys_agg.size), dtype=coord_dtype)
                    q = keys_agg
                    np.remainder(q, n_activities, out=coords[3])
                    np.floor_divide(q, n_activities, out=q)
                    np.remainder(q, n_years, out=coords[2])
                    np.floor_divide(q, n_years, out=q)
                    np.remainder(q, n_flows, out=coords[1])
                    np.floor_divide(q, n_flows, out=coords[0])
                else:
                    coords = np.empty((3, keys_agg.size), dtype=coord_dtype)
                    q = keys_agg
                    np.remainder(q, n_years, out=coords[2])
                    np.floor_divide(q, n_years, out=q)
                    np.remainder(q, n_flows, out=coords[1])
                    np.floor_divide(q, n_flows, out=coords[0])

                inv = sparse.COO(
                    coords,
                    data_agg,
                    shape=shape,
                    has_duplicates=False,
                    sorted=True,
                    idx_dtype=coord_dtype,
                )
            else:
                inv = sparse.zeros(shape, dtype=self.value_dtype)
        else:
            if has_root:
                inv = sparse.zeros(
                    (n_activities, n_flows, len(years), n_activities),
                    dtype=self.value_dtype,
                )
            else:
                inv = sparse.zeros(
                    (n_activities, n_flows, len(years)), dtype=self.value_dtype
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
        # Builders are no longer needed once inventory is finalized.
        self._inv_key_parts = []
        self._inv_value_parts = []
        self.inventory_diagnostics = {
            "backend": "coo",
            "raw_entries": eager_entries,
            "canonical_entries": int(inv.nnz),
            "estimated_peak_bytes": eager_peak,
        }
        return self.inventory

    def _materialize_dataarray(
        self,
        value: xr.DataArray | None,
        *,
        memory_budget: int | None,
        label: str,
    ) -> xr.DataArray:
        if value is None:
            raise ValueError(f"Trails.{label} is empty")
        if not is_chunked_sparse(value.data):
            return value
        builder = self._inventory_builder
        if builder is None:
            raise RuntimeError("Chunked inventory metadata is unavailable")
        methods = int(value.sizes.get("method", 1))
        predicted = estimate_materialization_peak_bytes(
            int(builder.nnz) * methods,
            ndim=int(value.ndim),
            value_dtype=value.dtype,
        )
        budget = (
            int(self._inventory_memory_budget)
            if memory_budget is None
            else int(memory_budget)
        )
        if predicted > budget:
            raise MemoryError(
                f"Materializing {label} requires an estimated "
                f"{predicted / 2**30:.2f} GiB, above the "
                f"{budget / 2**30:.2f} GiB budget."
            )
        data = value.data.compute(scheduler="synchronous")
        return value.copy(data=data)

    def materialize_inventory(
        self, *, memory_budget: int | None = None
    ) -> xr.DataArray:
        """Safely replace a lazy inventory with one eager sparse COO."""
        result = self._materialize_dataarray(
            self.inventory,
            memory_budget=memory_budget,
            label="inventory",
        )
        self.inventory = result
        return result

    def materialize_characterized_inventory(
        self, *, memory_budget: int | None = None
    ) -> xr.DataArray:
        """Safely replace a lazy characterized inventory with eager COO."""
        result = self._materialize_dataarray(
            self.characterized_inventory,
            memory_budget=memory_budget,
            label="characterized_inventory",
        )
        self.characterized_inventory = result
        return result

    def close(self) -> None:
        """Release managed disk-backed inventory stores."""
        for builder in getattr(self, "_inventory_builders", []):
            builder.close()
        self._inventory_builders = []

    def __enter__(self) -> "Trails":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

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
        """interpolate temporal exchange.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param other_idx: Value for `other_idx`.
        :type other_idx: int
        :param exchanges: Value for `exchanges`.
        :type exchanges: Dict[tuple, TemporalExchange]
        :returns: Return value.
        :rtype: Optional[TemporalExchange]"""
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
                    or tuple(getattr(tex0, "offsets", ()) or ())
                    != tuple(getattr(tex1, "offsets", ()) or ())
                    or tuple(getattr(tex0, "weights", ()) or ())
                    != tuple(getattr(tex1, "weights", ()) or ())
                ):
                    return tex0 if (year - y0) <= (y1 - year) else tex1

                w = (year - y0) / (y1 - y0)

                def interp_optional(v0: float | None, v1: float | None) -> float | None:
                    """Interp optional.

                    :param v0: Value for `v0`.
                    :type v0: float | None
                    :param v1: Value for `v1`.
                    :type v1: float | None
                    :returns: Return value.
                    :rtype: float | None"""
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
                    offsets=getattr(tex0, "offsets", None),
                    weights=getattr(tex0, "weights", None),
                )

        return None

    def _map_year_to_scenario_year(self, year: int) -> int:
        """map year to scenario year.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: int"""
        year_int = int(year)
        cached = self._scenario_year_map_cache.get(year_int)
        if cached is not None:
            return cached

        y = max(self.min_year, min(self.max_year, year_int))

        # If we have a full annual grid, this is effectively identity after clipping
        if len(self.years_int) == (self.max_year - self.min_year + 1):
            self._scenario_year_map_cache[year_int] = int(y)
            return int(y)

        # Otherwise: snap to nearest scenario year
        idx = int(np.abs(self.years_int - y).argmin())
        mapped = int(self.years_int[idx])
        self._scenario_year_map_cache[year_int] = mapped
        return mapped

    def _map_year_to_template_year(self, year: int) -> int:
        """map year to template year.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: int"""
        y = int(year)
        cached = self._template_year_map_cache.get(y)
        if cached is not None:
            return cached
        idx = int(np.abs(self.template_years_int - y).argmin())
        mapped = int(self.template_years_int[idx])
        self._template_year_map_cache[y] = mapped
        return mapped

    def _get_scenario_context(self, year: int) -> tuple[int, str, int] | None:
        """get scenario context.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: tuple[int, str, int] | None"""
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
        """add demand entry.

        :param demand: Value for `demand`.
        :type demand: dict[int, dict[int, float]]
        :param target_year: Value for `target_year`.
        :type target_year: int
        :param product_index: Value for `product_index`.
        :type product_index: int
        :param exchange_amount: Value for `exchange_amount`.
        :type exchange_amount: float"""
        demand.setdefault(target_year, {})
        demand[target_year][product_index] = (
            demand[target_year].get(product_index, 0.0) + exchange_amount
        )

    @staticmethod
    def _child_amount(parent_amount: float, exchange_value: float) -> float:
        """child amount.

        :param parent_amount: Value for `parent_amount`.
        :type parent_amount: float
        :param exchange_value: Value for `exchange_value`.
        :type exchange_value: float
        :returns: Return value.
        :rtype: float"""
        if exchange_value < 0.0:
            return parent_amount * (-exchange_value)
        return parent_amount * exchange_value

    def _production_amount(self, t: int, act_idx: int) -> float:
        """Return the absolute production amount for an activity in a scenario."""
        if self.A is None:
            return 1.0
        act = int(act_idx)
        vector = self._production_amount_vector(int(t))
        if act < 0 or act >= int(vector.size):
            return 1.0
        return float(vector[act])

    def _production_amount_vector(self, t: int) -> np.ndarray:
        """Return absolute production amounts for all activities in a scenario."""
        if self.A is None:
            return np.ones(0, dtype=np.float64)
        t_int = int(t)
        vector_cache = getattr(self, "_production_amount_vector_cache", None)
        if vector_cache is None:
            vector_cache = self._production_amount_vector_cache = {}
        cached_vector = vector_cache.get(t_int)
        if cached_vector is not None:
            return cached_vector

        n_activities = int(self.A.shape[1])
        production = np.ones(n_activities, dtype=np.float64)
        try:
            A_t = self.A[t_int, :, :]
            coords = np.asarray(A_t.coords, dtype=np.int64)
            data = np.asarray(A_t.data, dtype=np.float64)
            if coords.ndim == 2 and coords.shape[0] == 2 and data.size:
                diag_mask = coords[0] == coords[1]
                if np.any(diag_mask):
                    indices = coords[0, diag_mask].astype(np.intp, copy=False)
                    values = np.abs(data[diag_mask])
                    valid = (
                        (indices >= 0)
                        & (indices < n_activities)
                        & np.isfinite(values)
                        & (values >= 1e-30)
                    )
                    if np.any(valid):
                        production[indices[valid]] = values[valid]
        except Exception:
            pass

        vector_cache[t_int] = production

        # Keep the legacy scalar cache coherent for callers that inspect it.
        cache = getattr(self, "_production_amount_cache", None)
        if cache is None:
            cache = self._production_amount_cache = {}
        for act, amount in enumerate(production):
            cache[(t_int, int(act))] = float(amount)
        return production

    def _activity_amount_from_product_demand(
        self,
        t: int,
        act_idx: int,
        product_amount: float,
    ) -> float:
        """Convert product demand into an activity scaling amount."""
        return float(product_amount) / self._production_amount(t, act_idx)

    def _child_activity_amount(
        self,
        *,
        t: int,
        product_index: int,
        parent_amount: float,
        exchange_value: float,
    ) -> float:
        """Convert a technosphere product exchange into child activity scaling."""
        product_amount = self._child_amount(parent_amount, exchange_value)
        if product_amount == 0.0:
            return 0.0
        return float(product_amount) / self._production_amount(t, product_index)

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
        """apply temporal distribution to demand.

        :param year: Value for `year`.
        :type year: int
        :param product_index: Value for `product_index`.
        :type product_index: int
        :param child_amount: Value for `child_amount`.
        :type child_amount: float
        :param tex: Value for `tex`.
        :type tex: TemporalExchange
        :param demand: Value for `demand`.
        :type demand: dict[int, dict[int, float]]
        :param debug: Value for `debug`.
        :type debug: bool"""
        offsets_and_weights = self._get_td_offsets(tex=tex, debug=debug)
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

    def _get_tech_td_expanded(
        self, *, year: int, act_idx: int, prod_idx: int, debug: bool
    ) -> tuple[Optional[TemporalExchange], list[tuple[int, float]]]:
        """get tech td expanded.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param prod_idx: Value for `prod_idx`.
        :type prod_idx: int
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: tuple[Optional[TemporalExchange], list[tuple[int, float]]]"""
        y_tpl = self._map_year_to_template_year(year)
        key = (int(y_tpl), int(act_idx), int(prod_idx))
        cached = self._tech_td_expanded_cache.get(key)
        if cached is not None:
            return cached
        tex = self._get_tech_temporal_exchange(int(year), int(act_idx), int(prod_idx))
        if tex is None:
            self._tech_td_expanded_cache[key] = (None, [])
            return (None, [])
        offsets_and_weights = self._get_td_offsets(tex=tex, debug=debug)
        self._tech_td_expanded_cache[key] = (tex, offsets_and_weights)
        return tex, offsets_and_weights

    def _get_td_offsets(
        self, *, tex: TemporalExchange, debug: bool
    ) -> list[tuple[int, float]]:
        """get td offsets.

        :param tex: Value for `tex`.
        :type tex: TemporalExchange
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: list[tuple[int, float]]"""
        key = (
            int(tex.distribution),
            float(tex.loc) if tex.loc is not None else None,
            float(tex.scale) if tex.scale is not None else None,
            int(tex.offset_min),
            int(tex.offset_max),
            getattr(tex, "amount_source", "port"),
            tuple(getattr(tex, "offsets", ()) or ()),
            tuple(getattr(tex, "weights", ()) or ()),
        )
        cached = self._td_offsets_cache.get(key)
        if cached is not None:
            return cached
        td = TemporalDistribution(tex)
        offsets_and_weights = list(td.iter_offsets_and_weights(debug=debug))
        self._td_offsets_cache[key] = offsets_and_weights
        return offsets_and_weights

    def _get_tech_expansion_template(
        self,
        *,
        year: int,
        t: int,
        act_idx: int,
        use_temporal_distributions: bool,
        debug: bool,
    ) -> tuple[_TechExpansionEntry, ...]:
        """Return cached child-expansion entries for one technosphere row."""
        if self.A is None:
            return ()

        template_year = int(self._map_year_to_template_year(int(year)))
        cache_key = (
            int(template_year),
            int(t),
            int(act_idx),
            bool(use_temporal_distributions),
        )
        cached_template = self._tech_expansion_template_cache.get(cache_key)
        if cached_template is not None:
            return cached_template

        row_key = (int(t), int(act_idx))
        cached_row = self._A_row_cache.get(row_key)
        if cached_row is not None:
            product_indices, values = cached_row
        else:
            A_row = self.A[int(t), int(act_idx), :]
            if A_row.nnz == 0:
                self._tech_expansion_template_cache[cache_key] = ()
                return ()
            product_indices = A_row.coords[0]
            values = A_row.data
            self._A_row_cache[row_key] = (product_indices, values)

        entries: list[_TechExpansionEntry] = []
        for product_index_raw, exchange_value_raw in zip(product_indices, values):
            product_index = int(product_index_raw)
            exchange_value = float(exchange_value_raw)
            if exchange_value == 0.0 or product_index == int(act_idx):
                continue

            tex = None
            offsets_and_weights: list[tuple[int, float]] = []
            if use_temporal_distributions:
                tex, offsets_and_weights = self._get_tech_td_expanded(
                    year=int(year),
                    act_idx=int(act_idx),
                    prod_idx=product_index,
                    debug=debug,
                )

            if tex is None or not use_temporal_distributions:
                production_amount = self._production_amount(int(t), product_index)
                if production_amount == 0.0:
                    continue
                scale = self._child_amount(1.0, exchange_value) / production_amount
                if scale != 0.0:
                    entries.append(
                        _TechExpansionEntry(
                            child_act_idx=product_index,
                            scale=float(scale),
                            amount_source="none",
                        )
                    )
                continue

            if not offsets_and_weights:
                if debug:
                    logger.warning(
                        "expand_temporal_exchanges: TD produced no "
                        "offsets/weights for (year=%d prod=%d) -> dropping "
                        "exchange",
                        int(year),
                        product_index,
                    )
                continue

            offsets = tuple(int(offset) for offset, _weight in offsets_and_weights)
            weights = tuple(float(weight) for _offset, weight in offsets_and_weights)
            amount_source = str(getattr(tex, "amount_source", "port"))
            if amount_source == "matrix":
                entries.append(
                    _TechExpansionEntry(
                        child_act_idx=product_index,
                        scale=0.0,
                        amount_source="matrix",
                        offsets=offsets,
                        weights=weights,
                    )
                )
                continue

            production_amount = self._production_amount(int(t), product_index)
            if production_amount == 0.0:
                continue
            scale = self._child_amount(1.0, exchange_value) / production_amount
            if scale != 0.0:
                entries.append(
                    _TechExpansionEntry(
                        child_act_idx=product_index,
                        scale=float(scale),
                        amount_source="port",
                        offsets=offsets,
                        weights=weights,
                    )
                )

        result = tuple(entries)
        self._tech_expansion_template_cache[cache_key] = result
        return result

    def _expand_temporal_child_demands_fast(
        self,
        *,
        year: int,
        act_idx: int,
        amount: float,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> dict[tuple[int, int], float]:
        """Expand child demands as flat ``(year, activity)`` amounts."""
        context = self._get_scenario_context(int(year))
        if context is None:
            return {}
        _scenario_year, _scenario_label, t = context
        entries = self._get_tech_expansion_template(
            year=int(year),
            t=int(t),
            act_idx=int(act_idx),
            use_temporal_distributions=bool(use_temporal_distributions),
            debug=debug,
        )
        if not entries:
            return {}

        out: dict[tuple[int, int], float] = {}

        year_int = int(year)
        act_int = int(act_idx)
        parent_amount = float(amount)
        out_get = out.get
        scenario_index_get = self.scenario_index.get
        for entry in entries:
            child_act = int(entry.child_act_idx)
            if entry.amount_source == "matrix":
                for offset, weight in zip(entry.offsets, entry.weights):
                    if weight == 0.0:
                        continue
                    raw_year = year_int + offset
                    y_eff = int(self._map_year_to_scenario_year(raw_year))
                    t_eff = scenario_index_get(str(y_eff))
                    if t_eff is None:
                        continue
                    exchange_value = float(self.A[int(t_eff), act_int, child_act])
                    if exchange_value == 0.0:
                        continue
                    product_amount = self._child_amount(1.0, exchange_value)
                    production_amount = self._production_amount(int(t_eff), child_act)
                    if production_amount == 0.0:
                        continue
                    child_amount = (
                        parent_amount
                        * float(product_amount)
                        / float(production_amount)
                        * float(weight)
                    )
                    if child_amount != 0.0:
                        key = (raw_year, child_act)
                        out[key] = out_get(key, 0.0) + child_amount
                continue

            child_base_amount = parent_amount * float(entry.scale)
            if child_base_amount == 0.0:
                continue
            if not entry.offsets:
                key = (year_int, child_act)
                out[key] = out_get(key, 0.0) + child_base_amount
                continue
            for offset, weight in zip(entry.offsets, entry.weights):
                if weight == 0.0:
                    continue
                child_amount = child_base_amount * weight
                if child_amount != 0.0:
                    raw_year = year_int + offset
                    key = (raw_year, child_act)
                    out[key] = out_get(key, 0.0) + child_amount

        return out

    def get_A_for_scenario(self, label: str) -> sparse.COO:
        """Get a for scenario.

        :param label: Value for `label`.
        :type label: str
        :returns: Return value.
        :rtype: sparse.COO"""
        t = self.scenario_index[label]
        return self.A[t, :, :]

    def get_B_for_scenario(self, label: str) -> sparse.COO:
        """Get b for scenario.

        :param label: Value for `label`.
        :type label: str
        :returns: Return value.
        :rtype: sparse.COO"""
        t = self.scenario_index[label]
        return self.B[t, :, :]

    def get_temporal_exchange(
        self, year: int, act_idx: int, prod_idx: int
    ) -> TemporalExchange | None:
        """Get temporal exchange.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param prod_idx: Value for `prod_idx`.
        :type prod_idx: int
        :returns: Return value.
        :rtype: TemporalExchange | None"""
        return self._interpolate_temporal_exchange(
            year,
            act_idx,
            prod_idx,
            self.temporal_technosphere_exchanges,
        )

    def get_temporal_distribution(
        self, year: int, act_idx: int, prod_idx: int
    ) -> TemporalDistribution | None:
        """Get temporal distribution.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param prod_idx: Value for `prod_idx`.
        :type prod_idx: int
        :returns: Return value.
        :rtype: TemporalDistribution | None"""
        tex = self._interpolate_temporal_exchange(
            year,
            act_idx,
            prod_idx,
            self.temporal_technosphere_exchanges,
        )
        if tex is None:
            return None
        return TemporalDistribution(tex)

    def lca(self, *args: Any, **kwargs: Any) -> Any:
        """Run temporal LCA using the stored routing graph.

        :param args: Variadic positional arguments.
        :type args: Any
        :param kwargs: Keyword arguments forwarded to ``trails.lca.lca``. If
            ``methods``, ``edges_methods``, or ``ei_version`` are omitted, the
            defaults configured on this ``Trails`` instance are used.
        :type kwargs: Any
        :returns: Return value.
        :rtype: Any"""
        from .lca import lca as lca_fn

        if "debug" in kwargs:
            self.debug = bool(kwargs.pop("debug"))

        return lca_fn(self, *args, **kwargs)

    def lci(self, *args: Any, **kwargs: Any) -> xr.DataArray:
        """Build the temporal life-cycle inventory from the routing graph.

        This phase performs the year-specific linear solves and always retains
        the finalized inventory. It does not load or apply LCIA methods.
        """
        from .lca import lci as lci_fn

        if "debug" in kwargs:
            self.debug = bool(kwargs.pop("debug"))
        return lci_fn(self, *args, **kwargs)

    def lcia(self, *args: Any, **kwargs: Any) -> xr.DataArray:
        """Characterize the finalized temporal inventory.

        Methods supplied here override optional constructor defaults. Repeated
        calls reuse the same inventory and never rerun the linear systems.
        """
        from .lcia import lcia as lcia_fn

        return lcia_fn(self, *args, **kwargs)

    def temporal_routing(
        self,
        *,
        start_year: int,
        start_act_idx: int,
        amount: float = 1.0,
        max_depth: int | None = None,
        min_amount: float = 1e-18,
        show_progress: bool = True,
        attribute_to_roots: bool = True,
        debug: bool = False,
        adaptive_methods: str | list[str] | tuple[str, ...] | None = None,
        adaptive_relative_score_cutoff: Any = (_DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF),
        adaptive_ei_version: str | None = None,
        adaptive_min_depth: int = 1,
        adaptive_use_cache: bool = True,
    ) -> None:
        """Build the temporal routing graph for a functional unit.

        Routing follows temporalized technosphere exchanges from the root
        demand and stores an explicit graph of process-year nodes. Branches
        stop at frontier nodes when they reach an optional ``max_depth``, fall
        below ``min_amount``, have no child demands, or meet the adaptive
        score-potential cutoff. Frontier demands are still solved by ``lci()``
        in the corresponding year-specific background matrices, so adaptive
        pruning changes how much of the graph is routed explicitly, not whether
        the remaining demand is included.

        Adaptive routing precomputes static activity scores for the selected
        ``adaptive_methods`` and estimates each child branch's impact potential
        as ``abs(reference product demand) * max(abs(static activity score))``.
        The relative cutoff is multiplied by the functional-unit static score
        potential to get the effective score-potential threshold. The
        default public routing mode is adaptive: ``max_depth=None`` and
        ``adaptive_relative_score_cutoff=1e-4``. Passing an integer
        ``max_depth`` without an adaptive cutoff selects fixed-depth routing.
        Passing both an integer ``max_depth`` and a relative adaptive cutoff
        combines adaptive pruning with a hard depth cap.

        :param start_year: Scenario/calendar year of the functional unit.
        :type start_year: int
        :param start_act_idx: Activity index of the functional unit provider.
        :type start_act_idx: int
        :param amount: Functional-unit amount, expressed in the reference
            product of ``start_act_idx``.
        :type amount: float
        :param max_depth: Maximum explicit routing depth. The default ``None``
            lets the adaptive score cutoff define the stopping depth. Passing
            an integer without an adaptive cutoff uses fixed-depth routing.
        :type max_depth: int | None
        :param min_amount: Absolute routed-demand cutoff. Branches with child
            amounts below this value become frontier nodes.
        :type min_amount: float
        :param show_progress: Show a progress bar while routing.
        :type show_progress: bool
        :param attribute_to_roots: Track first-level supplier attribution for
            frontier and direct biosphere amounts.
        :type attribute_to_roots: bool
        :param debug: Enable detailed routing logging.
        :type debug: bool
        :param adaptive_methods: Optional regular LCIA method or methods used
            to screen branch impact potential. This routing configuration is
            independent of methods later passed to ``lcia()``. A deprecated
            constructor-method fallback remains for compatibility.
        :type adaptive_methods: str | list[str] | tuple[str, ...] | None
        :param adaptive_relative_score_cutoff: Relative cutoff multiplied by the
            functional-unit static score potential. Defaults to ``1e-4`` when
            adaptive routing is selected by default. Set to ``None`` to disable
            adaptive routing when using a fixed ``max_depth``.
        :type adaptive_relative_score_cutoff: float | None
        :param adaptive_ei_version: Ecoinvent LCIA version used for adaptive
            screening factors. If omitted, ``Trails(..., ei_version=...)`` is
            used.
        :type adaptive_ei_version: str | None
        :param adaptive_min_depth: Minimum child depth before adaptive pruning can
            stop a branch.
        :type adaptive_min_depth: int
        :param adaptive_use_cache: Reuse and write internal static activity score
            cache files.
        :type adaptive_use_cache: bool
        :raises ValueError: If adaptive routing is requested without explicit or
            constructor-default regular LCIA methods, or if ``max_depth=None``
            is used while the adaptive relative cutoff is disabled.
        :raises RuntimeError: If required routing dependencies are missing."""
        try:
            import networkx as nx
        except Exception as exc:  # pragma: no cover - dependency guard
            raise RuntimeError(
                "networkx is required for temporal_routing(). "
                "Install it with `pip install networkx`."
            ) from exc

        if debug:
            self.debug = True

        self._invalidate_calculation_results(close_inventory=True)

        start_activity = int(start_act_idx)
        start_year_int = int(start_year)
        start_amount = float(amount)

        # Hot-path caches
        year_cache: dict[int, int] = {}
        public_node_key_cache: dict[tuple[int, int, int], tuple] = {}
        meta_cache: dict[tuple[int, int], dict] = {}
        bio_cache: dict[tuple[int, int], bool] = {}
        node_attrs: dict[tuple[int, int, int], dict[str, Any]] = {}
        edge_amounts: dict[tuple[tuple[int, int, int], tuple[int, int, int]], float] = (
            {}
        )
        frontier_amounts: dict[tuple[int, int, int], float] = {}
        frontier_reasons: dict[str, dict[tuple[int, int, int], float]] = {
            "max_depth": {},
            "leaf": {},
            "min_amount": {},
            "adaptive_relative_score_cutoff": {},
        }
        frontier_max_depth = frontier_reasons["max_depth"]
        frontier_leaf = frontier_reasons["leaf"]
        frontier_min_amount = frontier_reasons["min_amount"]
        frontier_adaptive_cutoff = frontier_reasons["adaptive_relative_score_cutoff"]
        frontier_roots: dict[tuple[tuple[int, int, int], int], float] = {}
        adaptive_cutoff_nodes: set[tuple[int, int, int]] = set()
        adaptive_cutoff_potentials: dict[tuple[int, int, int], float] = {}

        def map_year(y: int) -> int:
            """Map year.

            :param y: Value for `y`.
            :type y: int
            :returns: Return value.
            :rtype: int"""
            yi = int(y)
            if yi in year_cache:
                return year_cache[yi]
            mapped = int(self._map_year_to_scenario_year(yi))
            year_cache[yi] = mapped
            return mapped

        def _get_activity_meta(label: str, idx: int) -> dict:
            """get activity meta.

            :param label: Value for `label`.
            :type label: str
            :param idx: Value for `idx`.
            :type idx: int
            :returns: Return value.
            :rtype: dict"""
            key = (int(label), int(idx)) if label.isdigit() else None
            if key is not None and key in meta_cache:
                return meta_cache[key]
            mapping = self.activity_indices.get(label)
            if mapping and idx in mapping:
                meta = mapping.get(idx, {})
                if key is not None:
                    meta_cache[key] = meta
                return meta
            for _label, _mapping in self.activity_indices.items():
                if idx in _mapping:
                    meta = _mapping.get(idx, {})
                    if key is not None:
                        meta_cache[key] = meta
                    return meta
            return {}

        def _public_node_key(key: tuple[int, int, int]) -> tuple:
            """Return public graph key with activity metadata."""
            cached = public_node_key_cache.get(key)
            if cached is not None:
                return cached
            year, depth, act_idx = key
            scenario_year = map_year(year)
            label = str(scenario_year)
            meta = _get_activity_meta(label, int(act_idx))
            name = meta.get("name") or ""
            ref_prod = meta.get("reference product") or ""
            location = meta.get("location") or ""
            public_key = (
                int(year),
                int(depth),
                name,
                ref_prod,
                location,
                int(act_idx),
            )
            public_node_key_cache[key] = public_key
            return public_key

        def _new_node_attrs(year: int, depth: int, act_idx: int) -> dict[str, Any]:
            """Return initial graph attributes for one compact routing node."""
            return {
                "year": year,
                "depth": depth,
                "act_idx": act_idx,
                "name": "",
                "reference_product": "",
                "location": "",
                "amount": 0.0,
                "frontier_amount": 0.0,
                "direct_bio_amount": 0.0,
                "score_potential": 0.0,
                "score_potential_by_method": {},
                "adaptive_cutoff_reason": None,
                "adaptive_effective_score_cutoff": None,
                "adaptive_cutoff_potential": 0.0,
                "frontier_reasons": {},
                "frontier_roots": {},
                "direct_bio_roots": {},
            }

        def _add_root_amount(
            data: dict[int, float], root_act: int | None, amt: float, fallback: int
        ) -> None:
            """add root amount.

            :param data: Value for `data`.
            :type data: dict[int, float]
            :param root_act: Value for `root_act`.
            :type root_act: int | None
            :param amt: Value for `amt`.
            :type amt: float
            :param fallback: Value for `fallback`.
            :type fallback: int"""
            root = int(root_act) if root_act is not None else int(fallback)
            data[root] = float(data.get(root, 0.0)) + float(amt)

        def _add_frontier_amount(
            node_key: tuple[int, int, int],
            amt: float,
            root_act: int | None,
            fallback: int,
            reason_amounts: dict[tuple[int, int, int], float],
        ) -> None:
            frontier_amounts[node_key] = frontier_amounts.get(node_key, 0.0) + amt
            reason_amounts[node_key] = reason_amounts.get(node_key, 0.0) + amt
            if attribute_to_roots:
                root = int(root_act) if root_act is not None else int(fallback)
                root_key = (node_key, root)
                frontier_roots[root_key] = frontier_roots.get(root_key, 0.0) + amt

        def _add_adaptive_frontier_amount(
            node_key: tuple[int, int, int],
            amt: float,
            root_act: int | None,
            fallback: int,
            score_potential: float,
        ) -> None:
            _add_frontier_amount(
                node_key,
                amt,
                root_act,
                fallback,
                frontier_adaptive_cutoff,
            )
            adaptive_cutoff_nodes.add(node_key)
            adaptive_cutoff_potentials[node_key] = max(
                adaptive_cutoff_potentials.get(node_key, 0.0),
                score_potential,
            )

        def _flush_frontier_amounts() -> None:
            """Write deferred frontier bookkeeping onto graph node attributes."""
            for node_key, amount in frontier_amounts.items():
                node_attrs[node_key]["frontier_amount"] = float(amount)
            for reason, reason_amounts in frontier_reasons.items():
                for node_key, amount in reason_amounts.items():
                    reasons = node_attrs[node_key].setdefault("frontier_reasons", {})
                    reasons[reason] = float(amount)
            for node_key in adaptive_cutoff_nodes:
                node = node_attrs[node_key]
                node["adaptive_cutoff_reason"] = "adaptive_relative_score_cutoff"
                node["adaptive_effective_score_cutoff"] = float(adaptive_threshold)
                node["adaptive_cutoff_potential"] = float(
                    adaptive_cutoff_potentials.get(node_key, 0.0)
                )
            if attribute_to_roots:
                for (node_key, root), amount in frontier_roots.items():
                    roots = node_attrs[node_key].setdefault("frontier_roots", {})
                    roots[int(root)] = float(amount)

        score_unit_cache: dict[tuple[int, int], tuple[float, tuple[float, ...]]] = {}
        score_amounts: dict[tuple, float] = {}

        def _score_product_amount_abs(
            *,
            year: int,
            act_idx: int,
            amt: float,
        ) -> float:
            """Return absolute reference-product demand for a routed amount."""
            context = self._get_scenario_context(int(year))
            if context is None:
                return abs(float(amt))
            _scenario_year, _scenario_label, t = context
            production_amount = self._production_amount(int(t), int(act_idx))
            return abs(float(amt)) * abs(float(production_amount))

        def _cached_unit_score_potential(
            *,
            year: int,
            act_idx: int,
        ) -> tuple[float, tuple[float, ...]]:
            """Return cached absolute static score potential for one unit."""
            if adaptive_scores is None:
                return 0.0, ()
            cache_key = (int(year), int(act_idx))
            cached = score_unit_cache.get(cache_key)
            if cached is not None:
                return cached
            unit_max, unit_values = adaptive_scores.unit_score_potential(
                year=int(year),
                activity=int(act_idx),
            )
            cached = (
                float(unit_max),
                tuple(float(value) for value in unit_values),
            )
            score_unit_cache[cache_key] = cached
            return cached

        def _score_potential(
            *,
            year: int,
            act_idx: int,
            amt: float,
        ) -> float:
            """Return adaptive static score potential for one routed amount."""
            if adaptive_scores is None:
                return 0.0
            unit_max, _unit_values = _cached_unit_score_potential(
                year=int(year),
                act_idx=int(act_idx),
            )
            amount_abs = _score_product_amount_abs(
                year=int(year),
                act_idx=int(act_idx),
                amt=float(amt),
            )
            return amount_abs * float(unit_max)

        def _flush_score_potentials() -> None:
            """Write deferred score-potential diagnostics onto graph nodes."""
            if adaptive_scores is None:
                return
            methods = adaptive_scores.methods
            for node_key, amount_abs in score_amounts.items():
                node = node_attrs[node_key]
                unit_max, unit_values = _cached_unit_score_potential(
                    year=int(node["year"]),
                    act_idx=int(node["act_idx"]),
                )
                node["score_potential"] = float(amount_abs) * float(unit_max)
                if len(methods) == 1:
                    node["score_potential_by_method"] = {
                        methods[0]: float(amount_abs) * float(unit_values[0])
                    }
                    continue
                node["score_potential_by_method"] = {
                    method: float(amount_abs) * float(value)
                    for method, value in zip(methods, unit_values)
                }

        queue = deque()
        start_year_int = int(self._map_year_to_scenario_year(start_year_int))
        start_context = self._get_scenario_context(start_year_int)
        if start_context is None:
            return
        _, _, start_t = start_context
        start_activity_amount = self._activity_amount_from_product_demand(
            int(start_t),
            start_activity,
            start_amount,
        )

        default_relative_cutoff = (
            adaptive_relative_score_cutoff is _DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF
        )
        if default_relative_cutoff:
            if max_depth is None or adaptive_methods is not None:
                adaptive_relative_score_cutoff = DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF
            else:
                adaptive_relative_score_cutoff = None

        adaptive_requested = adaptive_relative_score_cutoff is not None
        if adaptive_methods is None and adaptive_requested:
            default_methods = getattr(self, "default_methods", None)
            default_backend = getattr(self, "default_method_backend", "auto")
            if (
                default_methods
                and default_backend != "edges"
                and all(isinstance(method, str) for method in default_methods)
            ):
                warnings.warn(
                    "Using Trails(..., methods=...) for adaptive routing is "
                    "deprecated; pass adaptive_methods=... explicitly to "
                    "temporal_routing().",
                    FutureWarning,
                    stacklevel=2,
                )
                adaptive_methods = list(default_methods)
            elif default_methods or getattr(self, "default_edges_methods", None):
                raise ValueError(
                    "Adaptive routing requires regular LCIA methods. EDGES "
                    "methods can be used for final lcia(), but provide "
                    "adaptive_methods=... for adaptive screening."
                )
            else:
                raise ValueError(
                    "adaptive_methods (preferred) or Trails(..., methods=...) "
                    "must be provided when an adaptive relative score cutoff "
                    "is set."
                )
        adaptive_enabled = adaptive_methods is not None
        if adaptive_enabled and not adaptive_requested:
            raise ValueError(
                "Set adaptive_relative_score_cutoff "
                "when adaptive_methods is provided "
                "for adaptive routing."
            )
        if max_depth is None and not adaptive_enabled:
            raise ValueError("max_depth=None is only supported in adaptive mode.")
        max_depth_int = None if max_depth is None else int(max_depth)
        min_amount_float = float(min_amount)
        adaptive_min_depth_int = int(adaptive_min_depth)
        adaptive_ei_version_effective = (
            str(adaptive_ei_version)
            if adaptive_ei_version is not None
            else str(getattr(self, "default_ei_version", "3.11"))
        )

        adaptive_scores = None
        adaptive_threshold = 0.0
        adaptive_root_potential = 0.0
        adaptive_method_names: list[str] = []
        if adaptive_enabled:
            if isinstance(adaptive_methods, str):
                adaptive_method_names = [adaptive_methods]
            else:
                adaptive_method_names = [str(method) for method in adaptive_methods]
            if not adaptive_method_names:
                raise ValueError("adaptive_methods cannot be empty.")
            adaptive_scores = _ensure_static_activity_scores(
                self,
                methods=adaptive_method_names,
                ei_version=adaptive_ei_version_effective,
                use_cache=bool(adaptive_use_cache),
            )
            adaptive_root_potential, _ = _activity_score_potential(
                adaptive_scores,
                year=start_year_int,
                activity=start_activity,
                amount=start_amount,
            )
            adaptive_threshold = abs(float(adaptive_relative_score_cutoff)) * float(
                adaptive_root_potential
            )

        track_score_amounts = adaptive_scores is not None
        debug_enabled = bool(self.debug)
        track_roots = bool(attribute_to_roots)

        queue.append((start_year_int, start_activity, start_activity_amount, 0, None))

        if show_progress:
            pbar = tqdm(
                total=None,
                desc="Temporal routing",
                unit="node",
                dynamic_ncols=True,
            )
        else:
            pbar = None

        nodes_processed = 0
        max_processed_depth = 0

        while queue:
            year, act, amt, depth, root_act = queue.popleft()
            year = map_year(year)

            if amt == 0.0:
                continue

            nodes_processed += 1
            if int(depth) > max_processed_depth:
                max_processed_depth = int(depth)
            if pbar is not None:
                pbar.update(1)

            if root_act is None and depth > 0:
                root_act = int(act)

            node_key = (year, depth, act)
            node = node_attrs.get(node_key)
            if node is None:
                node = _new_node_attrs(year, depth, act)
                node_attrs[node_key] = node
            node["amount"] = node["amount"] + amt
            if track_score_amounts and depth == 0:
                score_amounts[node_key] = score_amounts.get(
                    node_key, 0.0
                ) + _score_product_amount_abs(year=year, act_idx=act, amt=amt)

            scenario_year = year
            has_direct_bio = self._has_direct_biosphere(scenario_year, act, bio_cache)

            if max_depth_int is not None and depth >= max_depth_int:
                _add_frontier_amount(
                    node_key,
                    amt,
                    root_act,
                    act,
                    frontier_max_depth,
                )
                continue

            child_demands = self._expand_temporal_child_demands_fast(
                year=year,
                act_idx=act,
                amount=amt,
                use_temporal_distributions=True,
                debug=debug_enabled,
            )

            if not child_demands:
                _add_frontier_amount(
                    node_key,
                    amt,
                    root_act,
                    act,
                    frontier_leaf,
                )
                continue

            if has_direct_bio and depth > 0:
                node["direct_bio_amount"] = node["direct_bio_amount"] + amt
                if track_roots:
                    _add_root_amount(node["direct_bio_roots"], root_act, amt, act)

            for (child_year, child_act), child_amt in child_demands.items():
                child_amt = float(child_amt)
                if child_amt == 0.0:
                    continue

                child_year = map_year(child_year)
                child_act = int(child_act)
                child_depth = depth + 1
                child_node = (child_year, child_depth, child_act)
                child_attrs = node_attrs.get(child_node)
                if child_attrs is None:
                    child_attrs = _new_node_attrs(child_year, child_depth, child_act)
                    node_attrs[child_node] = child_attrs

                edge_key = (node_key, child_node)
                edge_amounts[edge_key] = edge_amounts.get(edge_key, 0.0) + child_amt
                if track_score_amounts:
                    score_amounts[child_node] = score_amounts.get(
                        child_node, 0.0
                    ) + _score_product_amount_abs(
                        year=child_year,
                        act_idx=child_act,
                        amt=child_amt,
                    )

                if depth == 0:
                    child_root = child_act
                else:
                    child_root = root_act

                if max_depth_int is not None and child_depth >= max_depth_int:
                    child_attrs["amount"] = child_attrs["amount"] + child_amt
                    frontier_amounts[child_node] = (
                        frontier_amounts.get(child_node, 0.0) + child_amt
                    )
                    frontier_max_depth[child_node] = (
                        frontier_max_depth.get(child_node, 0.0) + child_amt
                    )
                    if track_roots:
                        root = (
                            int(child_root)
                            if child_root is not None
                            else int(child_act)
                        )
                        root_key = (child_node, root)
                        frontier_roots[root_key] = (
                            frontier_roots.get(root_key, 0.0) + child_amt
                        )
                    continue

                if abs(child_amt) < min_amount_float:
                    _add_frontier_amount(
                        child_node,
                        child_amt,
                        child_root,
                        child_act,
                        frontier_min_amount,
                    )
                    continue

                if adaptive_enabled and child_depth >= adaptive_min_depth_int:
                    potential = _score_potential(
                        year=child_year,
                        act_idx=child_act,
                        amt=child_amt,
                    )
                    if potential <= adaptive_threshold:
                        child_attrs["amount"] = child_attrs["amount"] + child_amt
                        _add_adaptive_frontier_amount(
                            child_node,
                            child_amt,
                            child_root,
                            child_act,
                            potential,
                        )
                        continue

                queue.append(
                    (
                        child_year,
                        child_act,
                        child_amt,
                        child_depth,
                        child_root,
                    )
                )

        if pbar is not None:
            pbar.close()

        _flush_frontier_amounts()
        _flush_score_potentials()

        G = nx.DiGraph()
        public_node_keys: dict[tuple[int, int, int], tuple] = {}
        for key, attrs in node_attrs.items():
            public_key = _public_node_key(key)
            public_node_keys[key] = public_key
            graph_attrs = dict(attrs)
            _year, _depth, name, ref_prod, location, _act_idx = public_key
            graph_attrs["name"] = name
            graph_attrs["reference_product"] = ref_prod
            graph_attrs["location"] = location
            G.add_node(public_key, **graph_attrs)
        G.add_edges_from(
            (public_node_keys[parent], public_node_keys[child], {"amount": amount})
            for (parent, child), amount in edge_amounts.items()
        )

        self.graph = G
        self._routing_attribute_to_roots = bool(attribute_to_roots)
        self._routing_params = {
            "start_year": start_year_int,
            "start_act_idx": start_activity,
            "amount": start_amount,
            "max_depth": None if max_depth_int is None else int(max_depth_int),
            "min_amount": float(min_amount),
            "nodes_processed": int(nodes_processed),
            "max_processed_depth": int(max_processed_depth),
            "adaptive_enabled": bool(adaptive_enabled),
            "adaptive_methods": list(adaptive_method_names),
            "adaptive_relative_score_cutoff": (
                None
                if adaptive_relative_score_cutoff is None
                else float(adaptive_relative_score_cutoff)
            ),
            "adaptive_effective_score_cutoff": (
                float(adaptive_threshold) if adaptive_enabled else None
            ),
            "adaptive_root_score_potential": (
                float(adaptive_root_potential) if adaptive_enabled else None
            ),
            "adaptive_ei_version": adaptive_ei_version_effective,
            "adaptive_min_depth": int(adaptive_min_depth),
        }

        if debug:
            logger.info(
                "temporal_routing done: nodes=%d edges=%d",
                int(G.number_of_nodes()),
                int(G.number_of_edges()),
            )

    def plot_adaptive_sankey(self, **kwargs: Any) -> Any:
        """Plot the current adaptive routed graph as a depth/year Sankey.

        This is a convenience wrapper around
        :func:`trails.plotting.plot_adaptive_sankey`. Run
        ``temporal_routing(...)`` first so ``self.graph`` contains the routed
        graph and adaptive score-potential attributes.

        :param kwargs: Keyword arguments forwarded to
            ``plot_adaptive_sankey``.
        :type kwargs: Any
        :returns: Plotly figure.
        :rtype: Any"""
        from .plotting import plot_adaptive_sankey

        return plot_adaptive_sankey(self, **kwargs)

    def static_lca(
        self,
        year: int,
        act_idx: int,
        methods: list[str] | None = None,
        amount: float = 1.0,
        debug: bool = False,
        ei_version: str | None = None,
    ) -> None:
        """Static lca.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param methods: Regular LCIA methods. If omitted, uses
            ``Trails(..., methods=...)``.
        :type methods: list[str] | None
        :param amount: Value for `amount`.
        :type amount: float
        :param debug: Value for `debug`.
        :type debug: bool
        :param ei_version: LCIA data version. If omitted, uses
            ``Trails(..., ei_version=...)``.
        :type ei_version: str | None
        :returns: Return value.
        :rtype: None"""
        lca_static(
            trails=self,
            year=int(year),
            fu_act_idx=int(act_idx),
            methods=methods,
            amount=float(amount),
            debug=debug,
            ei_version=ei_version,
        )
        return None

    def expand_temporal_exchanges(
        self,
        year: int,
        act_idx: int,
        amount: float = 1.0,
        *,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> dict[int, dict[int, float]]:
        """Expand temporal exchanges.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param amount: Value for `amount`.
        :type amount: float
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: dict[int, dict[int, float]]"""
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

        cache_key = (int(t), int(act_idx))
        cached = self._A_row_cache.get(cache_key)
        if cached is not None:
            product_indices, values = cached
        else:
            A_row = self.A[t, act_idx, :]
            if A_row.nnz == 0:
                return demand
            product_indices = A_row.coords[0]
            values = A_row.data
            self._A_row_cache[cache_key] = (product_indices, values)

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

            # Fetch TD metadata + offsets (cached per exchange)
            tex, offsets_and_weights = self._get_tech_td_expanded(
                year=year, act_idx=act_idx, prod_idx=product_index, debug=debug
            )

            # ------------------------------------------------------------------
            # No temporal distribution (or disabled): status quo
            # ------------------------------------------------------------------
            if (tex is None) or (not use_temporal_distributions):
                n_no_td += 1
                child_amount = self._child_activity_amount(
                    t=int(t),
                    product_index=product_index,
                    parent_amount=amount,
                    exchange_value=exchange_value,
                )
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
                self._apply_temporal_distribution_matrix_sourced_to_demand_offsets(
                    year=year,
                    act_idx=act_idx,
                    product_index=product_index,
                    parent_amount=amount,
                    offsets_and_weights=offsets_and_weights,
                    demand=demand,
                    debug=debug,
                )
                continue

            # ------------------------------------------------------------------
            # TD + ported magnitude (default): distribute anchor-year child amount
            # ------------------------------------------------------------------
            n_td_ported += 1
            child_amount = self._child_activity_amount(
                t=int(t),
                product_index=product_index,
                parent_amount=amount,
                exchange_value=exchange_value,
            )
            if child_amount == 0.0:
                continue

            if not offsets_and_weights:
                if debug:
                    logger.warning(
                        "expand_temporal_exchanges: TD produced no offsets/weights for (year=%d prod=%d) -> dropping exchange",
                        year,
                        product_index,
                    )
                continue
            for offset, weight in offsets_and_weights:
                raw_year = year + offset
                self._add_demand_entry(
                    demand,
                    int(raw_year),
                    product_index,
                    child_amount * float(weight),
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
        """get tech temporal exchange.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param prod_idx: Value for `prod_idx`.
        :type prod_idx: int
        :returns: Return value.
        :rtype: Optional[TemporalExchange]"""
        if not self.temporal_technosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        key = (int(y_tpl), int(act_idx), int(prod_idx))
        if key in self._tech_td_cache:
            return self._tech_td_cache[key]
        tex = self.temporal_technosphere_exchanges.get((str(y_tpl), key[1], key[2]))
        self._tech_td_cache[key] = tex
        return tex

    def _get_bio_temporal_exchange(
        self, year: int, act_idx: int, flow_idx: int
    ) -> Optional[TemporalExchange]:
        """get bio temporal exchange.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param flow_idx: Value for `flow_idx`.
        :type flow_idx: int
        :returns: Return value.
        :rtype: Optional[TemporalExchange]"""
        if not self.temporal_biosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        return self.temporal_biosphere_exchanges.get(
            (str(y_tpl), int(act_idx), int(flow_idx))
        )

    def print_exchange_table(
        self,
        *,
        year: int,
        act_idx: int,
        max_rows: int | None = None,
        sort_by_amount: bool = True,
    ) -> None:
        """Print exchange table.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param max_rows: Value for `max_rows`.
        :type max_rows: int | None
        :param sort_by_amount: Value for `sort_by_amount`.
        :type sort_by_amount: bool"""
        context = self._get_scenario_context(int(year))
        if context is None:
            print(f"No scenario data available for year={year}")
            return
        scenario_year, scenario_label, t = context
        scenario_label = str(scenario_label)

        def _format_amount(value: float) -> str:
            """format amount.

            :param value: Value for `value`.
            :type value: float
            :returns: Return value.
            :rtype: str"""
            if value == 0.0:
                return "0"
            abs_v = abs(value)
            if abs_v >= 1e4 or abs_v < 1e-2:
                return f"{value:.2e}"
            return f"{value:.2f}".rstrip("0").rstrip(".")

        def _get_activity_meta(label: str, idx: int) -> dict:
            """get activity meta.

            :param label: Value for `label`.
            :type label: str
            :param idx: Value for `idx`.
            :type idx: int
            :returns: Return value.
            :rtype: dict"""
            mapping = self.activity_indices.get(label)
            if mapping and idx in mapping:
                return mapping.get(idx, {})
            # fallback to any available scenario mapping
            for _label, _mapping in self.activity_indices.items():
                if idx in _mapping:
                    return _mapping.get(idx, {})
            return {}

        def _get_bio_meta(label: str, idx: int) -> dict:
            """get bio meta.

            :param label: Value for `label`.
            :type label: str
            :param idx: Value for `idx`.
            :type idx: int
            :returns: Return value.
            :rtype: dict"""
            mapping = self.biosphere_indices.get(label)
            if mapping and idx in mapping:
                return mapping.get(idx, {})
            for _label, _mapping in self.biosphere_indices.items():
                if idx in _mapping:
                    return _mapping.get(idx, {})
            return {}

        prod_rows: list[dict[str, object]] = []
        tech_rows: list[dict[str, object]] = []
        bio_rows: list[dict[str, object]] = []

        # ---- Technosphere + production exchanges (A matrix row) ----
        if self.A is not None:
            A_row = self.A[t, int(act_idx), :]
            if A_row.nnz:
                prod_indices = A_row.coords[0]
                values = A_row.data
                for prod_idx, value in zip(prod_indices, values):
                    prod_idx = int(prod_idx)
                    amount = float(value)
                    if amount == 0.0:
                        continue
                    flow_type = "prod" if prod_idx == int(act_idx) else "tech"
                    meta = _get_activity_meta(scenario_label, prod_idx)
                    tex = self._get_tech_temporal_exchange(
                        int(year), int(act_idx), prod_idx
                    )
                    entry = {
                        "direction": "out" if amount > 0.0 else "in",
                        "flow type": flow_type,
                        "id": prod_idx,
                        "name": meta.get("name") or "",
                        "reference product": meta.get("reference product") or "",
                        "amount": _format_amount(amount),
                        "td type": getattr(tex, "distribution", None),
                        "loc": getattr(tex, "loc", None),
                        "scale": getattr(tex, "scale", None),
                        "min": getattr(tex, "offset_min", None),
                        "max": getattr(tex, "offset_max", None),
                        "offsets": getattr(tex, "offsets", None),
                        "weights": getattr(tex, "weights", None),
                        "td source": getattr(tex, "amount_source", None),
                    }
                    if flow_type == "prod":
                        prod_rows.append(entry)
                    else:
                        tech_rows.append(entry)

        # ---- Biosphere exchanges (B matrix row) ----
        if self.B is not None:
            B_row = self.B[t, int(act_idx), :]
            if B_row.nnz:
                flow_indices = B_row.coords[0]
                values = B_row.data
                for flow_idx, value in zip(flow_indices, values):
                    flow_idx = int(flow_idx)
                    amount = float(value)
                    if amount == 0.0:
                        continue
                    meta = _get_bio_meta(scenario_label, flow_idx)
                    tex = self._get_bio_temporal_exchange(
                        int(year), int(act_idx), flow_idx
                    )
                    bio_rows.append(
                        {
                            "direction": "emission" if amount > 0.0 else "uptake",
                            "flow type": "bio",
                            "id": flow_idx,
                            "name": meta.get("name") or "",
                            "reference product": "",
                            "amount": _format_amount(amount),
                            "td type": getattr(tex, "distribution", None),
                            "loc": getattr(tex, "loc", None),
                            "scale": getattr(tex, "scale", None),
                            "min": getattr(tex, "offset_min", None),
                            "max": getattr(tex, "offset_max", None),
                            "offsets": getattr(tex, "offsets", None),
                            "weights": getattr(tex, "weights", None),
                            "td source": getattr(tex, "amount_source", None),
                        }
                    )

        if not (prod_rows or tech_rows or bio_rows):
            print(f"No exchanges found for act={act_idx} in year={scenario_year}")
            return

        if sort_by_amount:

            def _sort_key(row: dict[str, object]) -> float:
                """sort key.

                :param row: Value for `row`.
                :type row: dict[str, object]
                :returns: Return value.
                :rtype: float"""
                raw = row.get("amount", "0")
                try:
                    return abs(float(raw))
                except ValueError:
                    try:
                        return abs(float(str(raw)))
                    except ValueError:
                        return 0.0

            prod_rows.sort(key=_sort_key, reverse=True)
            tech_rows.sort(key=_sort_key, reverse=True)
            bio_rows.sort(key=_sort_key, reverse=True)

        rows = prod_rows + tech_rows + bio_rows

        if max_rows is not None:
            rows = rows[: int(max_rows)]

        def _truncate(value: object, limit: int) -> str:
            """truncate.

            :param value: Value for `value`.
            :type value: object
            :param limit: Value for `limit`.
            :type limit: int
            :returns: Return value.
            :rtype: str"""
            text = "" if value is None else str(value)
            if len(text) <= limit:
                return text
            return text[:limit]

        col_limits = {
            "name": 32,
            "reference product": 32,
            "td type": 5,
            "loc": 5,
            "scale": 5,
            "min": 5,
            "max": 5,
            "offsets": 18,
            "weights": 18,
            "td source": 7,
        }

        headers = [
            "flow type",
            "id",
            "name",
            "reference product",
            "amount",
            "td type",
            "loc",
            "scale",
            "min",
            "max",
            "offsets",
            "weights",
            "td source",
        ]

        try:
            from prettytable import PrettyTable

            dist_table = PrettyTable()
            dist_table.field_names = ["code", "distribution"]
            dist_table.add_row([1, "discrete (all mass at loc)"])
            dist_table.add_row([2, "lognormal"])
            dist_table.add_row([3, "normal"])
            dist_table.add_row([4, "uniform"])
            dist_table.add_row([5, "triangular"])
            dist_table.add_row([6, "discrete empirical (explicit pulses)"])
            print("Temporal distribution codes:")
            print(dist_table)

            fields_table = PrettyTable()
            fields_table.field_names = ["field", "meaning"]
            fields_table.add_row(["temporal_distribution", "distribution code"])
            fields_table.add_row(
                ["temporal_loc", "location parameter (mean/median/mode)"]
            )
            fields_table.add_row(["temporal_scale", "scale parameter (stddev/sigma)"])
            fields_table.add_row(["temporal_min", "minimum integer offset (inclusive)"])
            fields_table.add_row(["temporal_max", "maximum integer offset (inclusive)"])
            fields_table.add_row(["temporal_offsets", "JSON list of pulse offsets"])
            fields_table.add_row(["temporal_weights", "JSON list of pulse weights"])
            fields_table.add_row(["temporal_amount_source", "ported value or matrix"])
            print("Temporal distribution fields:")
            print(fields_table)

            table = PrettyTable()
            table.field_names = headers
            for row in rows:
                for col, limit in col_limits.items():
                    row[col] = _truncate(row.get(col), limit)
                table.add_row([row.get(h) for h in headers])
            print(table)
        except Exception:
            # Fallback to a simple fixed-width printout
            widths = {h: max(len(h), 12) for h in headers}
            for row in rows:
                for col, limit in col_limits.items():
                    row[col] = _truncate(row.get(col), limit)
                for h in headers:
                    widths[h] = max(widths[h], len(str(row.get(h, ""))))
            line = " | ".join(h.ljust(widths[h]) for h in headers)
            sep = "-+-".join("-" * widths[h] for h in headers)
            print(line)
            print(sep)
            for row in rows:
                print(" | ".join(str(row.get(h, "")).ljust(widths[h]) for h in headers))

    def list_lcia_methods(self, ei_version: str = "3.11") -> list[str]:
        """List lcia methods.

        :param ei_version: Value for `ei_version`.
        :type ei_version: str
        :returns: Return value.
        :rtype: list[str]"""
        from .lcia import get_lcia_method_names

        return get_lcia_method_names(ei_version=ei_version)

    def _get_biosphere_slice(
        self, base_year: int, debug: bool
    ) -> tuple[int, int, sparse.COO, int] | None:
        """get biosphere slice.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: tuple[int, int, sparse.COO, int] | None"""
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
        """get b csr for t.

        :param t: Value for `t`.
        :type t: int
        :returns: Return value.
        :rtype: sparse.GCXS
        :raises RuntimeError: If an error occurs."""
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
        """get b row cache for t.

        :param t: Value for `t`.
        :type t: int
        :returns: Return value.
        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
        :raises RuntimeError: If an error occurs."""
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
        """get b row index map for t act.

        :param t: Value for `t`.
        :type t: int
        :param act: Value for `act`.
        :type act: int
        :returns: Return value.
        :rtype: np.ndarray
        :raises RuntimeError: If an error occurs."""
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
        """row values for flows sorted.

        :param row_flows_sorted: Value for `row_flows_sorted`.
        :type row_flows_sorted: np.ndarray
        :param row_vals_sorted: Value for `row_vals_sorted`.
        :type row_vals_sorted: np.ndarray
        :param query_flows_sorted: Value for `query_flows_sorted`.
        :type query_flows_sorted: np.ndarray
        :returns: Return value.
        :rtype: np.ndarray"""
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
        """build bio accumulation context.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: _BioAccumulationContext | None"""
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
        """filter idx with keep.

        :param idx_full: Value for `idx_full`.
        :type idx_full: np.ndarray
        :param keep_full: Value for `keep_full`.
        :type keep_full: np.ndarray
        :returns: Return value.
        :rtype: np.ndarray"""
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
        """get b cf activity vector.

        :param t: Value for `t`.
        :type t: int
        :param cf: Value for `cf`.
        :type cf: np.ndarray
        :returns: Return value.
        :rtype: np.ndarray"""
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
        method_idx: int | None = None,
    ) -> None:
        """Accumulate temporalized biosphere score matrix.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param supply_matrix: Value for `supply_matrix`.
        :type supply_matrix: np.ndarray
        :param root_activities: Value for `root_activities`.
        :type root_activities: np.ndarray
        :param cf: Value for `cf`.
        :type cf: np.ndarray
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :raises ValueError: If an error occurs."""
        if self.B is None:
            return
        if supply_matrix.size == 0:
            return

        base_year = int(base_year)

        # Ensure score builders exist
        if not hasattr(self, "_score_year_index") or not hasattr(
            self, "_score_bulk_value"
        ):
            self.reset_scores(
                attribute_to_roots=True,
                methods=getattr(self, "_score_methods", None),
            )

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
                method_idx=method_idx,
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
            """Td key.

            :param tex: Value for `tex`.
            :type tex: TemporalExchange
            :returns: Return value.
            :rtype: tuple"""
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
                tuple(getattr(tex, "offsets", ()) or ()),
                tuple(getattr(tex, "weights", ()) or ()),
            )

        def pulses_from_key(k: tuple) -> list[tuple[int, float]]:
            """Pulses from key.

            :param k: Value for `k`.
            :type k: tuple
            :returns: Return value.
            :rtype: list[tuple[int, float]]"""
            dist, loc, scale, off_min, off_max, amt_src, offsets, weights = k
            tex = TemporalExchange(
                distribution=dist,
                loc=loc,
                scale=scale,
                offset_min=off_min,
                offset_max=off_max,
                amount_source=amt_src,
                offsets=offsets,
                weights=weights,
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
                    """Map year cached.

                    :param raw_year: Value for `raw_year`.
                    :type raw_year: int
                    :returns: Return value.
                    :rtype: int"""
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

            self._append_scores_bulk(
                act,
                yr,
                val,
                root_activity=root,
                method_idx=method_idx,
            )

    def accumulate_temporalized_biosphere_score_matrix_multi(
        self,
        base_year: int,
        supply_matrix: np.ndarray,
        root_activities: np.ndarray,
        cf_matrix: np.ndarray,
        *,
        method_indices: np.ndarray | None = None,
        min_amount: float = 0.0,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """Accumulate temporalized scores for several LCIA methods in one pass."""
        if self.B is None:
            return
        if supply_matrix.size == 0:
            return

        base_year = int(base_year)

        if not hasattr(self, "_score_year_index") or not hasattr(
            self, "_score_bulk_value"
        ):
            self.reset_scores(
                attribute_to_roots=True,
                methods=getattr(self, "_score_methods", None),
            )

        biosphere_slice = self._get_biosphere_slice(base_year, debug)
        if biosphere_slice is None:
            return
        _scenario_year, t, _B_t, _ = biosphere_slice

        year_to_idx = self._score_year_index
        base_year_idx = year_to_idx.get(base_year)
        if base_year_idx is None:
            return

        C = np.asarray(cf_matrix, dtype=np.float64)
        if C.ndim == 1:
            C = C[None, :]
        if C.ndim != 2 or C.shape[1] != int(self.B.shape[2]):
            raise ValueError("cf_matrix must have shape (methods, B flows)")

        n_methods = int(C.shape[0])
        if method_indices is None:
            method_idx_values = np.arange(n_methods, dtype=np.int64)
        else:
            method_idx_values = np.asarray(method_indices, dtype=np.int64)
            if method_idx_values.shape != (n_methods,):
                raise ValueError("method_indices must match cf_matrix method axis")

        X = np.asarray(supply_matrix, dtype=np.float64)
        if X.ndim != 2:
            raise ValueError("supply_matrix must be 2D (n_acts, n_roots)")
        roots = np.asarray(root_activities, dtype=np.int64)
        if roots.ndim != 1 or roots.size != X.shape[1]:
            raise ValueError("root_activities must be 1D length n_roots")

        n_acts = int(self.B.shape[1])
        if X.shape[0] != n_acts:
            raise ValueError(f"supply_matrix has {X.shape[0]} acts; expected {n_acts}")

        out_act: list[np.ndarray] = []
        out_year: list[np.ndarray] = []
        out_root: list[np.ndarray] = []
        out_method: list[np.ndarray] = []
        out_val: list[np.ndarray] = []

        def append_coefficients(
            *,
            act: int,
            year_idx: int,
            coeffs: np.ndarray,
            x_row: np.ndarray,
        ) -> None:
            coeff_arr = np.asarray(coeffs, dtype=np.float64)
            if coeff_arr.shape != (n_methods,):
                raise ValueError("coefficient vector must match cf_matrix methods")
            if not np.any(coeff_arr) or not np.any(x_row):
                return
            values = coeff_arr[:, None] * x_row[None, :]
            method_pos, root_pos = np.nonzero(values != 0.0)
            if method_pos.size == 0:
                return
            out_act.append(np.full(method_pos.size, int(act), dtype=np.int64))
            out_year.append(np.full(method_pos.size, int(year_idx), dtype=np.int64))
            out_root.append(roots[root_pos].astype(np.int64, copy=False))
            out_method.append(method_idx_values[method_pos])
            out_val.append(values[method_pos, root_pos].astype(np.float64, copy=False))

        row_ptr, flow_sorted, data_sorted = self._get_B_row_cache_for_t(int(t))

        if not use_temporal_distributions or not self.temporal_biosphere_exchanges:
            active_acts = np.where(np.any(X != 0.0, axis=1))[0]
            for a in active_acts:
                start = int(row_ptr[a])
                end = int(row_ptr[a + 1])
                if start == end:
                    continue
                flows = flow_sorted[start:end].astype(np.intp, copy=False)
                vals = data_sorted[start:end].astype(np.float64, copy=False)
                coeffs = C[:, flows] @ vals
                append_coefficients(
                    act=int(a),
                    year_idx=int(base_year_idx),
                    coeffs=coeffs,
                    x_row=X[a, :],
                )

            if out_val:
                self._append_scores_bulk(
                    np.concatenate(out_act),
                    np.concatenate(out_year),
                    np.concatenate(out_val).astype(np.float64, copy=False),
                    root_activity=np.concatenate(out_root),
                    method_idx=np.concatenate(out_method),
                )
            return

        tpl_label = str(self._map_year_to_template_year(base_year))
        bio_td_get = self.temporal_biosphere_exchanges.get

        if not hasattr(self, "_td_pulse_cache"):
            self._td_pulse_cache = {}
        pulse_cache = self._td_pulse_cache

        if not hasattr(self, "_bio_score_row_char_matrix_cache"):
            self._bio_score_row_char_matrix_cache = {}
        row_char_cache = self._bio_score_row_char_matrix_cache

        def td_key(tex: TemporalExchange) -> tuple:
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
                tuple(getattr(tex, "offsets", ()) or ()),
                tuple(getattr(tex, "weights", ()) or ()),
            )

        def pulses_from_key(k: tuple) -> list[tuple[int, float]]:
            dist, loc, scale, off_min, off_max, amt_src, offsets, weights = k
            tex = TemporalExchange(
                distribution=dist,
                loc=loc,
                scale=scale,
                offset_min=off_min,
                offset_max=off_max,
                amount_source=amt_src,
                offsets=offsets,
                weights=weights,
            )
            return [
                (int(o), float(w))
                for o, w in TemporalDistribution(tex).iter_offsets_and_weights(
                    debug=False
                )
            ]

        year_map_cache: dict[int, int] = {}
        t_eff_cache: dict[int, int | None] = {}
        B_row_cache_local: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        scenario_index_get = self.scenario_index.get

        def map_year_cached(raw_year: int) -> int:
            y = year_map_cache.get(raw_year)
            if y is None:
                y = int(self._map_year_to_scenario_year(raw_year))
                year_map_cache[raw_year] = y
            return y

        active_acts = np.where(np.any(X != 0.0, axis=1))[0]
        cf_key = int(id(C))
        for a in active_acts:
            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            flows_full = flow_sorted[start:end].astype(np.intp, copy=False)
            vals_full = data_sorted[start:end].astype(np.float64, copy=False)

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
                        port_groups_pos.setdefault(td_key(tex), []).append(p)

                if no_td_pos:
                    pos = np.asarray(no_td_pos, dtype=np.intp)
                    no_td_coeff = C[:, flows_full[pos]] @ vals_full[pos]
                else:
                    no_td_coeff = np.zeros(n_methods, dtype=np.float64)

                ported_coeffs = {}
                for k, plist in port_groups_pos.items():
                    pos = np.asarray(plist, dtype=np.intp)
                    ported_coeffs[k] = (
                        C[:, flows_full[pos]] @ vals_full[pos]
                        if pos.size
                        else np.zeros(n_methods, dtype=np.float64)
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
            x_row = X[a, :]

            append_coefficients(
                act=int(a),
                year_idx=int(base_year_idx),
                coeffs=no_td_coeff,
                x_row=x_row,
            )

            for k, coeff_k in ported_coeffs.items():
                if not np.any(coeff_k):
                    continue
                pulses = pulse_cache.get(k)
                if pulses is None:
                    pulses = pulses_from_key(k)
                    pulse_cache[k] = pulses
                for offset, weight in pulses:
                    if weight == 0.0:
                        continue
                    raw = int(base_year + offset)
                    y_clamped = self._clamp_year_to_scores(raw)
                    yidx = year_to_idx[int(y_clamped)]
                    append_coefficients(
                        act=int(a),
                        year_idx=int(yidx),
                        coeffs=coeff_k * float(weight),
                        x_row=x_row,
                    )

            if not matrix_groups:
                continue

            for k, idx_full in matrix_groups.items():
                if idx_full.size == 0:
                    continue
                f_arr = flows_full[idx_full]
                if f_arr.size == 0:
                    continue
                ord_f = np.argsort(f_arr, kind="mergesort")
                f_sorted = f_arr[ord_f].astype(np.intp, copy=False)

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

                    row_index_map = self._get_B_row_index_map_for_t_act(t_eff_i, a)
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
                        f_use = f_sorted[valid]
                    else:
                        pos = row_index_map[f_sorted]
                        valid = pos >= 0
                        if not np.any(valid):
                            continue
                        vals_eff = row_vals_eff[pos[valid]]
                        f_use = f_sorted[valid]

                    score_per_supply = (C[:, f_use] @ vals_eff) * float(weight)
                    append_coefficients(
                        act=int(a),
                        year_idx=int(yidx),
                        coeffs=score_per_supply,
                        x_row=x_row,
                    )

        if out_val:
            self._append_scores_bulk(
                np.concatenate(out_act),
                np.concatenate(out_year),
                np.concatenate(out_val).astype(np.float64, copy=False),
                root_activity=np.concatenate(out_root),
                method_idx=np.concatenate(out_method),
            )

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
        method_idx: int | None = None,
    ) -> None:
        """Accumulate temporalized biosphere score.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param supply_by_activity: Value for `supply_by_activity`.
        :type supply_by_activity: Dict[int, float]
        :param cf: Value for `cf`.
        :type cf: np.ndarray
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param store_activity: Value for `store_activity`.
        :type store_activity: int | None
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :raises ValueError: If an error occurs."""
        # Ensure score builders exist
        if not hasattr(self, "_score_year_index") or not hasattr(
            self, "_score_chunk_value"
        ):
            self.reset_scores(
                attribute_to_roots=bool(getattr(self, "_scores_has_root", False)),
                methods=getattr(self, "_score_methods", None),
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
            """Map year cached.

            :param raw_year: Value for `raw_year`.
            :type raw_year: int
            :returns: Return value.
            :rtype: int"""
            y = year_map_cache.get(raw_year)
            if y is None:
                y = int(map_year_to_scenario(raw_year))
                year_map_cache[raw_year] = y
            return y

        def td_key(tex: TemporalExchange) -> tuple:
            """Td key.

            :param tex: Value for `tex`.
            :type tex: TemporalExchange
            :returns: Return value.
            :rtype: tuple"""
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
                tuple(getattr(tex, "offsets", ()) or ()),
                tuple(getattr(tex, "weights", ()) or ()),
            )

        def pulses_from_key(k: tuple) -> list[tuple[int, float]]:
            """Pulses from key.

            :param k: Value for `k`.
            :type k: tuple
            :returns: Return value.
            :rtype: list[tuple[int, float]]"""
            dist, loc, scale, off_min, off_max, amt_src, offsets, weights = k
            tex = TemporalExchange(
                distribution=dist,
                loc=loc,
                scale=scale,
                offset_min=off_min,
                offset_max=off_max,
                amount_source=amt_src,
                offsets=offsets,
                weights=weights,
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
            """safe filter positions.

            :param pos: Value for `pos`.
            :type pos: np.ndarray
            :param keep: Value for `keep`.
            :type keep: np.ndarray | None
            :param row_len: Value for `row_len`.
            :type row_len: int
            :returns: Return value.
            :rtype: np.ndarray"""
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
                        score_act,
                        {base_year_idx: score},
                        root_activity=root_activity,
                        method_idx=method_idx,
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
                    score_act,
                    acc_yearidx,
                    root_activity=root_activity,
                    method_idx=method_idx,
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
        include_no_td: bool = True,
        include_ported_td: bool = True,
        include_matrix_td: bool = True,
    ) -> None:
        """accumulate temporalized biosphere inventory core.

        :param ctx: Value for `ctx`.
        :type ctx: _BioAccumulationContext
        :param supply_by_activity: Value for `supply_by_activity`.
        :type supply_by_activity: Dict[int, float]
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param store_activity: Value for `store_activity`.
        :type store_activity: int | None
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :param workers: Optional worker count used by the no-TD batch fast path.
        :type workers: int | None"""
        dbg = self._get_debug_flow_filters(debug=debug)

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
        if not hasattr(self, "_td_pulse_array_cache"):
            self._td_pulse_array_cache = {}  # type: ignore[attr-defined]
        pulse_array_cache = self._td_pulse_array_cache  # type: ignore[attr-defined]

        year_map_cache = ctx.year_map_cache
        t_eff_cache = ctx.t_eff_cache
        B_row_cache_local = ctx.B_row_cache_local

        if use_temporal_distributions and ctx.bio_td_get is None:
            use_temporal_distributions = False

        def map_year_cached(raw_year: int) -> int:
            """Map year cached.

            :param raw_year: Value for `raw_year`.
            :type raw_year: int
            :returns: Return value.
            :rtype: int"""
            y = year_map_cache.get(raw_year)
            if y is None:
                y = int(map_year_to_scenario(raw_year))
                year_map_cache[raw_year] = y
            return y

        def td_key(tex: TemporalExchange) -> tuple:
            """Td key.

            :param tex: Value for `tex`.
            :type tex: TemporalExchange
            :returns: Return value.
            :rtype: tuple"""
            return (
                tex.distribution,
                tex.loc,
                tex.scale,
                tex.offset_min,
                tex.offset_max,
                getattr(tex, "amount_source", "port"),
                tuple(getattr(tex, "offsets", ()) or ()),
                tuple(getattr(tex, "weights", ()) or ()),
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

        append_bulk = self._append_inventory_entries_bulk
        has_root = bool(getattr(self, "_inventory_has_root", False))
        store_act = int(store_activity) if store_activity is not None else None
        base_t = int(ctx.t)

        dbg_flow_id = None if dbg is None else dbg.get("flow_id")
        dbg_year = None if dbg is None else dbg.get("year")
        dbg_act = None if dbg is None else dbg.get("act")
        dbg_max_pulses = 0 if dbg is None else int(dbg.get("max_pulses", 12))
        dbg_max_matches = 0 if dbg is None else int(dbg.get("max_matches", 50))
        ported_kernel = self._get_numba_ported_kernel()

        queued_act_values: list[int] = []
        queued_year_parts: list[np.ndarray | int] = []
        queued_flow_parts: list[np.ndarray] = []
        queued_value_parts: list[np.ndarray] = []
        queued_root_values: list[int] = []
        queued_nnz = 0
        flush_nnz = int(getattr(self, "_bio_inventory_flush_nnz", 2_000_000))
        if flush_nnz < 100_000:
            flush_nnz = 100_000
        matrix_year_groups_cache: dict[
            tuple, list[tuple[int, np.ndarray, np.ndarray]]
        ] = {}

        def flush_queued_inventory() -> None:
            nonlocal queued_nnz
            if not queued_flow_parts:
                return
            if len(queued_flow_parts) == 1:
                flow_all = queued_flow_parts[0]
                value_all = queued_value_parts[0]
                total_nnz = int(flow_all.size)
                year_part = queued_year_parts[0]
                if isinstance(year_part, np.ndarray):
                    year_all = year_part
                else:
                    year_all = np.empty(total_nnz, dtype=np.int64)
                    year_all.fill(int(year_part))
                act_all = np.empty(total_nnz, dtype=np.int64)
                act_all.fill(int(queued_act_values[0]))
                if has_root:
                    root_all = np.empty(total_nnz, dtype=np.int64)
                    root_all.fill(int(queued_root_values[0]))
                else:
                    root_all = None
            else:
                flow_all = np.concatenate(queued_flow_parts)
                value_all = np.concatenate(queued_value_parts)
                total_nnz = int(flow_all.size)

                year_all = np.empty(total_nnz, dtype=np.int64)
                act_all = np.empty(total_nnz, dtype=np.int64)
                root_all = np.empty(total_nnz, dtype=np.int64) if has_root else None
                cursor = 0
                for i, flow_part in enumerate(queued_flow_parts):
                    n = int(flow_part.size)
                    year_part = queued_year_parts[i]
                    if isinstance(year_part, np.ndarray):
                        year_all[cursor : cursor + n] = year_part
                    else:
                        year_all[cursor : cursor + n] = int(year_part)
                    act_all[cursor : cursor + n] = int(queued_act_values[i])
                    if has_root and root_all is not None:
                        root_all[cursor : cursor + n] = int(queued_root_values[i])
                    cursor += n

            append_bulk(
                act_all,
                year_all,
                flow_all,
                value_all,
                root_activity=root_all,
            )

            queued_act_values.clear()
            queued_year_parts.clear()
            queued_flow_parts.clear()
            queued_value_parts.clear()
            queued_root_values.clear()
            queued_nnz = 0

        for act_idx, supply_amt in supply_by_activity.items():
            supply_amt = float(supply_amt)
            if supply_amt == 0.0:
                continue

            a = int(act_idx)
            if has_root:
                inventory_act = a
                root_activity = store_act if store_act is not None else a
            else:
                inventory_act = store_act if store_act is not None else a
                root_activity = None
            if a < 0 or a + 1 >= len(row_ptr):
                continue

            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            flows_full = np.asarray(flow_sorted[start:end], dtype=np.intp)
            vals_full = data_sorted[start:end]

            scaled_full = supply_amt * vals_full.astype(np.float64, copy=False)
            flows_full_i64 = flows_full.astype(np.int64, copy=False)

            anchor_flow_parts: list[np.ndarray] = []
            anchor_value_parts: list[np.ndarray] = []
            row_flow_parts: list[np.ndarray] = []
            row_year_parts: list[np.ndarray] = []
            row_value_parts: list[np.ndarray] = []

            if dbg is not None and dbg_flow_id is not None:
                if dbg_act is None or int(dbg_act) == int(a):
                    if dbg_year is None or int(dbg_year) == int(base_year):
                        if dbg.get("matches", 0) < dbg_max_matches:
                            pos = np.where(flows_full == int(dbg_flow_id))[0]
                            if pos.size:
                                logger.debug(
                                    "bio_inv_row: base_year=%d act=%d flow=%d pos=%s raw_vals=%s scaled_vals=%s",
                                    int(base_year),
                                    int(a),
                                    int(dbg_flow_id),
                                    pos.tolist(),
                                    [float(vals_full[p]) for p in pos[:5]],
                                    [float(scaled_full[p]) for p in pos[:5]],
                                )
                                dbg["matches"] = int(dbg.get("matches", 0)) + 1

            thr = float(min_amount)
            if thr > 0.0:
                temporalize = np.abs(scaled_full) >= thr
            else:
                temporalize = None

            # TD metadata is stable for a given template/scenario-slice/activity row signature.
            cache_key = (tpl_label, base_t, int(a), int(start), int(end))

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

            if include_no_td and no_td_idx is not None:
                idx = no_td_idx
                if idx.size:
                    anchor_flow_parts.append(flows_full[idx])
                    anchor_value_parts.append(scaled_full[idx])
                    if dbg is not None and dbg_flow_id is not None:
                        if dbg_act is None or int(dbg_act) == int(a):
                            pos = idx[flows_full[idx] == int(dbg_flow_id)]
                            if pos.size and (
                                dbg_year is None or int(dbg_year) == int(base_year)
                            ):
                                logger.debug(
                                    "bio_inv_no_td: year=%d act=%d flow=%d contrib=%s",
                                    int(base_year),
                                    int(a),
                                    int(dbg_flow_id),
                                    float(np.sum(scaled_full[pos])),
                                )

            # -------------------------
            # PORTED TD groups (min_amount controls temporalization only)
            # -------------------------
            if include_ported_td and ported_groups:
                for k, idx_full in ported_groups.items():
                    if idx_full is None or idx_full.size == 0:
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
                        anchor_flow_parts.append(flows_full[idx_anchor])
                        anchor_value_parts.append(scaled_full[idx_anchor])

                    # 2) Temporalize above-threshold contributions
                    if idx_td.size == 0:
                        continue

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        # Reconstruct a TemporalExchange from td_key tuple (your existing convention)
                        (
                            dist,
                            loc,
                            scale,
                            off_min,
                            off_max,
                            amt_src,
                            offsets,
                            weights,
                        ) = k
                        tex0 = TemporalExchange(
                            distribution=dist,
                            loc=loc,
                            scale=scale,
                            offset_min=off_min,
                            offset_max=off_max,
                            amount_source=amt_src,
                            offsets=offsets,
                            weights=weights,
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
                        anchor_flow_parts.append(flows_full[idx_td])
                        anchor_value_parts.append(scaled_full[idx_td])
                        continue

                    if dbg is not None and dbg_flow_id is not None:
                        if dbg_act is None or int(dbg_act) == int(a):
                            flow_mask = flows_full[idx_td] == int(dbg_flow_id)
                            if np.any(flow_mask):
                                flow_vals = scaled_full[idx_td][flow_mask]
                                for offset, weight in pulses[:dbg_max_pulses]:
                                    raw_year = int(base_year + int(offset))
                                    if dbg_year is None or raw_year == int(dbg_year):
                                        contrib = float(flow_vals.sum()) * float(weight)
                                        logger.debug(
                                            "bio_inv_ported: base_year=%d act=%d flow=%d raw_year=%d weight=%.6g contrib=%.6g",
                                            int(base_year),
                                            int(a),
                                            int(dbg_flow_id),
                                            int(raw_year),
                                            float(weight),
                                            float(contrib),
                                        )

                    # Expand pulses in a vectorized way:
                    pulse_arrays = pulse_array_cache.get(k)
                    if pulse_arrays is None:
                        offsets_arr = np.fromiter(
                            (o for o, _ in pulses), dtype=np.int64, count=len(pulses)
                        )
                        weights_arr = np.fromiter(
                            (w for _, w in pulses), dtype=np.float64, count=len(pulses)
                        )
                        pulse_array_cache[k] = (offsets_arr, weights_arr)
                    else:
                        offsets_arr, weights_arr = pulse_arrays
                    if ported_kernel is not None:
                        flows_use, years_use, contrib = ported_kernel(
                            flows_full_i64,
                            scaled_full,
                            idx_td.astype(np.int64, copy=False),
                            offsets_arr,
                            weights_arr,
                            int(base_year),
                        )
                    else:
                        idx_rep = np.repeat(idx_td, len(pulses))
                        offsets_rep = np.tile(offsets_arr, idx_td.size)
                        weights_rep = np.tile(weights_arr, idx_td.size)

                        # Build contributions
                        flows_use = flows_full[idx_rep]
                        years_use = base_year + offsets_rep
                        contrib = scaled_full[idx_rep] * weights_rep

                    row_flow_parts.append(flows_use)
                    row_year_parts.append(years_use)
                    row_value_parts.append(contrib)

            if include_matrix_td and matrix_entries:
                matrix_kernel = self._get_numba_matrix_kernel()
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
                            anchor_flow_parts.append(flows_full[idx_anchor])
                            anchor_value_parts.append(scaled_full[idx_anchor])
                        # TD only for above-threshold
                        idx = idx[temporalize[idx]]
                        if idx.size == 0:
                            continue

                    f_arr = flows_full[idx]
                    if f_arr.size == 0:
                        continue
                    f_arr_i64 = f_arr.astype(np.int64, copy=False)

                    year_groups_cached = matrix_year_groups_cache.get(k)
                    if year_groups_cached is None:
                        year_groups: dict[int, list[tuple[int, float]]] = {}
                        for offset, weight in zip(offsets_arr, weights_arr):
                            if weight == 0.0:
                                continue
                            raw_year = base_year + int(offset)
                            y_eff = map_year_cached(raw_year)
                            year_groups.setdefault(int(y_eff), []).append(
                                (int(raw_year), float(weight))
                            )
                        year_groups_cached = []
                        for y_eff, year_weights in year_groups.items():
                            years_vec = np.fromiter(
                                (yw[0] for yw in year_weights),
                                dtype=np.int64,
                                count=len(year_weights),
                            )
                            weights_vec = np.fromiter(
                                (yw[1] for yw in year_weights),
                                dtype=np.float64,
                                count=len(year_weights),
                            )
                            year_groups_cached.append(
                                (int(y_eff), years_vec, weights_vec)
                            )
                        matrix_year_groups_cache[k] = year_groups_cached

                    for y_eff, years_vec, weights_vec in year_groups_cached:
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
                        row_ptr_eff, _flow_sorted_eff, data_sorted_eff = cached_eff

                        start_eff = int(row_ptr_eff[a])
                        end_eff = int(row_ptr_eff[a + 1])
                        if start_eff == end_eff:
                            continue

                        row_vals_eff = data_sorted_eff[start_eff:end_eff]
                        row_index_map = self._get_B_row_index_map_for_t_act(t_eff_i, a)
                        if matrix_kernel is not None:
                            vals_eff, valid = matrix_kernel(
                                f_arr_i64,
                                row_index_map,
                                row_vals_eff,
                            )
                            if not np.any(valid):
                                continue
                            vals_eff = vals_eff[valid]
                            f_use = f_arr[valid]
                        else:
                            pos = row_index_map[f_arr_i64]
                            valid = pos >= 0
                            if not np.any(valid):
                                continue
                            vals_eff = row_vals_eff[pos[valid]]
                            f_use = f_arr[valid]

                        if dbg is not None and dbg_flow_id is not None:
                            if dbg_act is None or int(dbg_act) == int(a):
                                flow_mask = f_use == int(dbg_flow_id)
                                if np.any(flow_mask):
                                    for raw_year, weight, v in zip(
                                        years_vec,
                                        weights_vec,
                                        vals_eff[flow_mask],
                                    ):
                                        if dbg_year is None or int(raw_year) == int(
                                            dbg_year
                                        ):
                                            contrib = (
                                                float(supply_amt)
                                                * float(v)
                                                * float(weight)
                                            )
                                            logger.debug(
                                                "bio_inv_matrix: base_year=%d act=%d flow=%d raw_year=%d weight=%.6g Bval=%.6g contrib=%.6g",
                                                int(base_year),
                                                int(a),
                                                int(dbg_flow_id),
                                                int(raw_year),
                                                float(weight),
                                                float(v),
                                                float(contrib),
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

                        row_flow_parts.append(flows_rep)
                        row_year_parts.append(years_rep)
                        row_value_parts.append(contrib)

            if anchor_flow_parts:
                if len(anchor_flow_parts) == 1:
                    flows_anchor = anchor_flow_parts[0]
                    values_anchor = anchor_value_parts[0]
                else:
                    flows_anchor = np.concatenate(anchor_flow_parts)
                    values_anchor = np.concatenate(anchor_value_parts)
                if flows_anchor.size:
                    nnz_anchor = int(flows_anchor.size)
                    queued_act_values.append(int(inventory_act))
                    queued_year_parts.append(int(base_year))
                    if (
                        flows_anchor.dtype == np.int64
                        and flows_anchor.flags.c_contiguous
                    ):
                        queued_flow_parts.append(flows_anchor)
                    else:
                        queued_flow_parts.append(
                            np.asarray(flows_anchor, dtype=np.int64, order="C")
                        )
                    if values_anchor.dtype == np.float64:
                        queued_value_parts.append(values_anchor)
                    else:
                        queued_value_parts.append(
                            np.asarray(values_anchor, dtype=np.float64)
                        )
                    if has_root:
                        queued_root_values.append(int(root_activity))
                    queued_nnz += nnz_anchor
                    if queued_nnz >= flush_nnz:
                        flush_queued_inventory()

            if row_flow_parts:
                if len(row_flow_parts) == 1:
                    flows_all = row_flow_parts[0]
                    years_all = row_year_parts[0]
                    values_all = row_value_parts[0]
                else:
                    flows_all = np.concatenate(row_flow_parts)
                    years_all = np.concatenate(row_year_parts)
                    values_all = np.concatenate(row_value_parts)

                if flows_all.size:
                    nnz_all = int(flows_all.size)
                    queued_act_values.append(int(inventory_act))
                    if years_all.dtype == np.int64 and years_all.flags.c_contiguous:
                        queued_year_parts.append(years_all)
                    else:
                        queued_year_parts.append(
                            np.asarray(years_all, dtype=np.int64, order="C")
                        )
                    if flows_all.dtype == np.int64 and flows_all.flags.c_contiguous:
                        queued_flow_parts.append(flows_all)
                    else:
                        queued_flow_parts.append(
                            np.asarray(flows_all, dtype=np.int64, order="C")
                        )
                    if values_all.dtype == np.float64:
                        queued_value_parts.append(values_all)
                    else:
                        queued_value_parts.append(
                            np.asarray(values_all, dtype=np.float64)
                        )
                    if has_root:
                        queued_root_values.append(int(root_activity))
                    queued_nnz += nnz_all
                    if queued_nnz >= flush_nnz:
                        flush_queued_inventory()

        flush_queued_inventory()

    def _get_numba_ported_kernel(self) -> Callable | None:
        """get numba ported kernel.

        :returns: Return value.
        :rtype: Callable | None"""
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
            flow_idx: np_local.ndarray,
            offsets: np_local.ndarray,
            weights: np_local.ndarray,
            base_year: int,
        ) -> tuple[np_local.ndarray, np_local.ndarray, np_local.ndarray]:
            """ported kernel.

            :param flows_full: Value for `flows_full`.
            :type flows_full: np_local.ndarray
            :param scaled_full: Value for `scaled_full`.
            :type scaled_full: np_local.ndarray
            :param flow_idx: Value for `flow_idx`.
            :type flow_idx: np_local.ndarray
            :param offsets: Value for `offsets`.
            :type offsets: np_local.ndarray
            :param weights: Value for `weights`.
            :type weights: np_local.ndarray
            :param base_year: Value for `base_year`.
            :type base_year: int
            :returns: Return value.
            :rtype: tuple[np_local.ndarray, np_local.ndarray, np_local.ndarray]"""
            n_flows = flow_idx.size
            n_pulses = offsets.size
            n = n_flows * n_pulses
            flows_out = np_local.empty(n, dtype=np_local.int64)
            years_out = np_local.empty(n, dtype=np_local.int64)
            contrib_out = np_local.empty(n, dtype=np_local.float64)
            k = 0
            for i in range(n_flows):
                idx = flow_idx[i]
                flow = flows_full[idx]
                scaled = scaled_full[idx]
                for j in range(n_pulses):
                    flows_out[k] = flow
                    years_out[k] = base_year + offsets[j]
                    contrib_out[k] = scaled * weights[j]
                    k += 1
            return flows_out, years_out, contrib_out

        setattr(self, "_numba_ported_kernel", _ported_kernel)
        return _ported_kernel

    def _get_numba_matrix_kernel(self) -> Callable | None:
        """get numba matrix kernel.

        :returns: Return value.
        :rtype: Callable | None"""
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
            """matrix kernel.

            :param flows: Value for `flows`.
            :type flows: np_local.ndarray
            :param row_index_map: Value for `row_index_map`.
            :type row_index_map: np_local.ndarray
            :param row_vals: Value for `row_vals`.
            :type row_vals: np_local.ndarray
            :returns: Return value.
            :rtype: tuple[np_local.ndarray, np_local.ndarray]"""
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

    def _accumulate_no_td_supply_matrix(
        self,
        ctx: _BioAccumulationContext,
        supply_matrix: np.ndarray,
        root_ids: np.ndarray,
        *,
        min_amount: float,
        excluded_activities: set[int] | None = None,
    ) -> None:
        """Accumulate non-temporal biosphere rows for all roots in bounded chunks."""
        matrix = np.asarray(supply_matrix, dtype=np.float64)
        roots = np.asarray(root_ids, dtype=np.int64)
        if matrix.ndim != 2:
            raise ValueError("supply_matrix must be two-dimensional")
        if matrix.shape[0] != int(ctx.n_acts):
            raise ValueError(
                f"supply_matrix has {matrix.shape[0]} activities; "
                f"expected {ctx.n_acts}"
            )
        if roots.ndim != 1 or roots.size != matrix.shape[1]:
            raise ValueError(
                "root_activities must be one-dimensional and aligned with "
                "supply_matrix columns"
            )
        n_roots = int(roots.size)
        nnz = int(ctx.data.shape[0])
        if n_roots == 0 or nnz == 0:
            return

        # Bound the dense B-row x root work array plus its sparse coordinate
        # extraction. Keeping roughly one million candidate contributions per
        # chunk caps temporary memory near 64 MiB with int64 coordinates.
        configured_budget = int(
            getattr(self, "_inventory_memory_budget", DEFAULT_INVENTORY_MEMORY_BUDGET)
        )
        temporary_budget = max(
            8 * 2**20,
            min(64 * 2**20, configured_budget // 2),
        )
        candidate_limit = max(n_roots, temporary_budget // 64)
        chunk_size = max(1, min(200_000, candidate_limit // n_roots))
        excluded = excluded_activities or set()
        excluded_array = np.fromiter(excluded, dtype=np.int64, count=len(excluded))

        for start in range(0, nnz, chunk_size):
            end = min(start + chunk_size, nnz)
            act_chunk = ctx.act_coords[start:end].astype(np.int64, copy=False)
            flow_chunk = ctx.flow_coords[start:end].astype(np.int64, copy=False)
            data_chunk = ctx.data[start:end].astype(np.float64, copy=False)

            if excluded:
                keep_rows = ~np.isin(act_chunk, excluded_array)
                if not np.any(keep_rows):
                    continue
                act_chunk = act_chunk[keep_rows]
                flow_chunk = flow_chunk[keep_rows]
                data_chunk = data_chunk[keep_rows]

            supply_chunk = matrix[act_chunk, :]
            active_supply = supply_chunk != 0.0
            if not np.any(active_supply):
                continue

            values = data_chunk[:, None] * supply_chunk
            values[~active_supply] = 0.0
            row_indices, root_columns = np.nonzero(values)
            if row_indices.size == 0:
                continue
            self._append_inventory_entries_bulk(
                act_chunk[row_indices],
                int(ctx.base_year),
                flow_chunk[row_indices],
                values[row_indices, root_columns],
                root_activity=roots[root_columns],
            )

    def _accumulate_no_td_batch(
        self,
        ctx: _BioAccumulationContext,
        supplies: List[tuple[Dict[int, float], int | None]],
        *,
        min_amount: float,
        debug: bool,
        workers: int | None = None,
    ) -> bool:
        """accumulate no td batch.

        :param ctx: Value for `ctx`.
        :type ctx: _BioAccumulationContext
        :param supplies: Value for `supplies`.
        :type supplies: List[tuple[Dict[int, float], int | None]]
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: bool"""
        if not supplies:
            return True

        has_root = bool(getattr(self, "_inventory_has_root", False))
        if not has_root:
            return False

        if any(store_activity is None for _, store_activity in supplies):
            return False

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
        self._accumulate_no_td_supply_matrix(
            ctx,
            supply_matrix,
            root_ids,
            min_amount=min_amount,
        )
        return True

    def accumulate_temporalized_biosphere_inventory_matrix(
        self,
        base_year: int,
        supply_matrix: np.ndarray,
        root_activities: np.ndarray,
        *,
        min_amount: float = 0.0,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> None:
        """Accumulate a root-attributed supply matrix without Python dictionaries.

        Biosphere rows without temporal distributions are expanded across all
        root columns in bounded NumPy chunks. Activities owning at least one
        temporal biosphere exchange retain the established row-wise path, so
        ported and matrix-sourced temporal semantics remain unchanged.
        """
        if not bool(getattr(self, "_inventory_has_root", False)):
            raise ValueError(
                "Matrix inventory accumulation requires root-attributed storage"
            )
        ctx = self._build_bio_accumulation_context(
            base_year,
            use_temporal_distributions=use_temporal_distributions,
            debug=debug,
        )
        if ctx is None:
            return

        matrix = np.asarray(supply_matrix, dtype=np.float64)
        roots = np.asarray(root_activities, dtype=np.int64)
        temporal_activities: set[int] = set()
        if use_temporal_distributions and ctx.tpl_label is not None:
            temporal_activities = {
                int(key[1])
                for key in self.temporal_biosphere_exchanges
                if len(key) >= 3 and str(key[0]) == str(ctx.tpl_label)
            }

        builder = self._inventory_builder
        if isinstance(builder, FactorizedInventoryBuilder):
            inventory_year = self._clamp_year_to_inventory(int(ctx.base_year))
            year_index = self._inventory_year_index.get(inventory_year)
            if year_index is None:
                raise RuntimeError(
                    "Inventory year index is missing for factorized accumulation"
                )
            n_flows = int(self.B.shape[2]) if self.B is not None else 0
            row_codes = ctx.act_coords.astype(
                np.int64, copy=False
            ) * n_flows + ctx.flow_coords.astype(np.int64, copy=False)
            temporal_by_code: dict[int, TemporalExchange] = {}
            if use_temporal_distributions and ctx.tpl_label is not None:
                temporal_by_code = {
                    int(key[1]) * n_flows + int(key[2]): tex
                    for key, tex in self.temporal_biosphere_exchanges.items()
                    if len(key) >= 3 and str(key[0]) == str(ctx.tpl_label)
                }

            port_groups: dict[tuple, list[int]] = {}
            matrix_positions: list[int] = []
            if temporal_by_code:
                temporal_codes = np.fromiter(
                    temporal_by_code, dtype=np.int64, count=len(temporal_by_code)
                )
                candidate_positions = np.flatnonzero(np.isin(row_codes, temporal_codes))

                def factor_td_key(tex: TemporalExchange) -> tuple:
                    return (
                        tex.distribution,
                        tex.loc,
                        tex.scale,
                        tex.offset_min,
                        tex.offset_max,
                        getattr(tex, "amount_source", "port"),
                        tuple(getattr(tex, "offsets", ()) or ()),
                        tuple(getattr(tex, "weights", ()) or ()),
                    )

                for position in candidate_positions:
                    tex = temporal_by_code[int(row_codes[int(position)])]
                    if getattr(tex, "amount_source", "port") == "matrix":
                        matrix_positions.append(int(position))
                    else:
                        port_groups.setdefault(factor_td_key(tex), []).append(
                            int(position)
                        )

            ordinary_rows = np.ones(ctx.data.size, dtype=bool)
            if matrix_positions:
                ordinary_rows[np.asarray(matrix_positions, dtype=np.int64)] = False

            ported_row_parts: list[np.ndarray] = []
            pulse_count_parts: list[np.ndarray] = []
            pulse_year_parts: list[np.ndarray] = []
            pulse_weight_parts: list[np.ndarray] = []
            if port_groups:
                years_axis = self._inventory_years
                if years_axis is None or not years_axis.size:
                    raise RuntimeError("Inventory years are not initialized")
                first_year = int(years_axis[0])
                last_year = int(years_axis[-1])
                for key, positions_list in port_groups.items():
                    pulses = ctx.pulse_cache.get(key)
                    if pulses is None:
                        (
                            dist,
                            loc,
                            scale,
                            off_min,
                            off_max,
                            amount_source,
                            offsets,
                            weights,
                        ) = key
                        tex = TemporalExchange(
                            distribution=dist,
                            loc=loc,
                            scale=scale,
                            offset_min=off_min,
                            offset_max=off_max,
                            amount_source=amount_source,
                            offsets=offsets,
                            weights=weights,
                        )
                        pulses = [
                            (int(offset), float(weight))
                            for offset, weight in TemporalDistribution(
                                tex
                            ).iter_offsets_and_weights(debug=False)
                        ]
                        ctx.pulse_cache[key] = pulses
                    positions = np.asarray(positions_list, dtype=np.int64)
                    if not pulses:
                        continue
                    offsets_array = np.fromiter(
                        (offset for offset, _ in pulses),
                        dtype=np.int64,
                        count=len(pulses),
                    )
                    weights_array = np.fromiter(
                        (weight for _, weight in pulses),
                        dtype=np.float64,
                        count=len(pulses),
                    )
                    pulse_years = (
                        np.clip(
                            int(ctx.base_year) + offsets_array,
                            first_year,
                            last_year,
                        )
                        - first_year
                    ).astype(np.int64, copy=False)
                    ordinary_rows[positions] = False
                    ported_row_parts.append(positions)
                    pulse_count_parts.append(
                        np.full(positions.size, len(pulses), dtype=np.int64)
                    )
                    pulse_year_parts.append(np.tile(pulse_years, positions.size))
                    pulse_weight_parts.append(np.tile(weights_array, positions.size))

            builder.append_factor(
                year_index=int(year_index),
                activities=ctx.act_coords,
                flows=ctx.flow_coords,
                biosphere_values=ctx.data,
                supply_matrix=matrix,
                roots=roots,
                included_rows=ordinary_rows,
            )
            if ported_row_parts:
                ported_rows = np.concatenate(ported_row_parts)
                pulse_counts = np.concatenate(pulse_count_parts)
                pulse_indptr = np.concatenate(
                    [
                        np.array([0], dtype=np.int64),
                        np.cumsum(pulse_counts, dtype=np.int64),
                    ]
                )
                builder.append_temporal_factor(
                    base_year_index=int(year_index),
                    activities=ctx.act_coords[ported_rows],
                    flows=ctx.flow_coords[ported_rows],
                    biosphere_values=ctx.data[ported_rows],
                    supply_matrix=matrix,
                    roots=roots,
                    pulse_indptr=pulse_indptr,
                    pulse_year_indices=np.concatenate(pulse_year_parts),
                    pulse_weights=np.concatenate(pulse_weight_parts),
                    min_amount=float(min_amount),
                )

            matrix_activities = {
                int(ctx.act_coords[position]) for position in matrix_positions
            }
            for activity in matrix_activities:
                if activity < 0 or activity >= matrix.shape[0]:
                    continue
                supplies = matrix[activity, :]
                root_columns = np.flatnonzero(supplies != 0.0)
                for column in root_columns:
                    self._accumulate_temporalized_biosphere_inventory_core(
                        ctx,
                        {int(activity): float(supplies[int(column)])},
                        min_amount=min_amount,
                        store_activity=int(roots[int(column)]),
                        use_temporal_distributions=use_temporal_distributions,
                        debug=debug,
                        include_no_td=False,
                        include_ported_td=False,
                        include_matrix_td=True,
                    )
            return
        else:
            self._accumulate_no_td_supply_matrix(
                ctx,
                matrix,
                roots,
                min_amount=min_amount,
                excluded_activities=temporal_activities,
            )

        for activity in temporal_activities:
            if activity < 0 or activity >= matrix.shape[0]:
                continue
            supplies = matrix[activity, :]
            root_columns = np.flatnonzero(supplies != 0.0)
            for column in root_columns:
                self._accumulate_temporalized_biosphere_inventory_core(
                    ctx,
                    {int(activity): float(supplies[int(column)])},
                    min_amount=min_amount,
                    store_activity=int(roots[int(column)]),
                    use_temporal_distributions=use_temporal_distributions,
                    debug=debug,
                )

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
        """Accumulate temporalized biosphere inventory.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param supply_by_activity: Value for `supply_by_activity`.
        :type supply_by_activity: Dict[int, float]
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param store_activity: Value for `store_activity`.
        :type store_activity: int | None
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool"""
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
        workers: int | None = None,
    ) -> None:
        """Accumulate temporalized biosphere inventory batch.

        :param base_year: Value for `base_year`.
        :type base_year: int
        :param supplies: Value for `supplies`.
        :type supplies: List[tuple[Dict[int, float], int | None]]
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool"""
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
                ctx,
                supplies,
                min_amount=min_amount,
                debug=debug,
                workers=workers,
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
        """map year to available.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: int"""
        return self._map_year_to_scenario_year(year)

    @staticmethod
    def _estimate_total_from_depth(max_depth: int) -> int | None:
        """estimate total from depth.

        :param max_depth: Value for `max_depth`.
        :type max_depth: int
        :returns: Return value.
        :rtype: int | None"""
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
        """record frontier.

        :param frontier_total: Value for `frontier_total`.
        :type frontier_total: dict[tuple[int, int], float]
        :param provenance_roots: Value for `provenance_roots`.
        :type provenance_roots: dict[tuple[int, int], dict[int, float]]
        :param y: Value for `y`.
        :type y: int
        :param a: Value for `a`.
        :type a: int
        :param x: Value for `x`.
        :type x: float
        :param r: Value for `r`.
        :type r: Optional[int]
        :param return_provenance: Value for `return_provenance`.
        :type return_provenance: bool"""
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
        """record direct bio.

        :param direct_bio_total: Value for `direct_bio_total`.
        :type direct_bio_total: dict[tuple[int, int], float]
        :param direct_bio_roots: Value for `direct_bio_roots`.
        :type direct_bio_roots: dict[tuple[int, int], dict[int, float]]
        :param y: Value for `y`.
        :type y: int
        :param a: Value for `a`.
        :type a: int
        :param x: Value for `x`.
        :type x: float
        :param r: Value for `r`.
        :type r: Optional[int]
        :param return_provenance: Value for `return_provenance`.
        :type return_provenance: bool"""
        direct_bio_total[(int(y), int(a))] += float(x)
        if return_provenance and (r is not None):
            direct_bio_roots[(int(y), int(a))][int(r)] += float(x)

    def _has_direct_biosphere(
        self, scenario_year: int, act: int, bio_cache: dict
    ) -> bool:
        """has direct biosphere.

        :param scenario_year: Value for `scenario_year`.
        :type scenario_year: int
        :param act: Value for `act`.
        :type act: int
        :param bio_cache: Value for `bio_cache`.
        :type bio_cache: dict
        :returns: Return value.
        :rtype: bool"""
        label = str(scenario_year)
        if label not in self.scenario_index or self.B is None:
            return False
        key = (scenario_year, act)
        if key in bio_cache:
            return bio_cache[key]
        t = self.scenario_index[label]
        # Build per-year cache once: boolean array of direct biosphere presence.
        if t not in self._direct_bio_cache_by_year:
            # Convert row indices to presence flags.
            rows = self.B[t, :, :].coords[0] if self.B is not None else np.array([])
            present = np.zeros(int(self.B.shape[1]), dtype=bool)
            if rows.size:
                present[rows] = True
            self._direct_bio_cache_by_year[t] = present
        present = self._direct_bio_cache_by_year[t]
        has_direct_bio = bool(present[int(act)])
        bio_cache[key] = has_direct_bio
        return has_direct_bio

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
        """Temporal traversal.

        :param start_year: Value for `start_year`.
        :type start_year: int
        :param start_act_idx: Value for `start_act_idx`.
        :type start_act_idx: int
        :param amount: Value for `amount`.
        :type amount: float
        :param max_depth: Value for `max_depth`.
        :type max_depth: int
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :param return_provenance: Value for `return_provenance`.
        :type return_provenance: bool
        :param show_progress: Value for `show_progress`.
        :type show_progress: bool
        :param use_temporal_distributions: Value for `use_temporal_distributions`.
        :type use_temporal_distributions: bool
        :param debug: Value for `debug`.
        :type debug: bool
        :returns: Return value.
        :rtype: tuple[dict, dict] | tuple[dict, dict, dict, dict]"""

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
            """Estimate total from branching.

            :param branching_samples: Value for `branching_samples`.
            :type branching_samples: list[int]
            :returns: Return value.
            :rtype: int"""
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
            """pbar step."""
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
            """pbar finalize."""
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

            # If we expand a node, we must still account for its direct biosphere flows.
            if has_direct_bio and depth > 0:
                direct_root = root_act if root_act is not None else int(act)
                self._record_direct_bio(
                    direct_bio_total,
                    direct_bio_roots,
                    year,
                    act,
                    amt,
                    direct_root,
                    return_provenance,
                )

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
        """Frontier to demand vectors.

        :param frontier: Value for `frontier`.
        :type frontier: dict
        :returns: Return value.
        :rtype: dict[int, np.ndarray]
        :raises ValueError: If an error occurs."""
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
        """Collect traversal edges.

        :param start_year: Value for `start_year`.
        :type start_year: int
        :param start_act_idx: Value for `start_act_idx`.
        :type start_act_idx: int
        :param amount: Value for `amount`.
        :type amount: float
        :param max_depth: Value for `max_depth`.
        :type max_depth: int
        :param min_amount: Value for `min_amount`.
        :type min_amount: float
        :returns: Return value.
        :rtype: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]]"""

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
        """apply temporal distribution matrix sourced to demand.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param product_index: Value for `product_index`.
        :type product_index: int
        :param parent_amount: Value for `parent_amount`.
        :type parent_amount: float
        :param tex: Value for `tex`.
        :type tex: TemporalExchange
        :param demand: Value for `demand`.
        :type demand: dict[int, dict[int, float]]
        :param debug: Value for `debug`.
        :type debug: bool"""
        if self.A is None:
            return

        if int(product_index) == int(act_idx):
            return

        offsets_and_weights = self._get_td_offsets(tex=tex, debug=debug)
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
            child_amount = self._child_activity_amount(
                t=int(t_eff),
                product_index=int(product_index),
                parent_amount=float(parent_amount),
                exchange_value=exchange_value,
            )
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

    def _apply_temporal_distribution_matrix_sourced_to_demand_offsets(
        self,
        *,
        year: int,
        act_idx: int,
        product_index: int,
        parent_amount: float,
        offsets_and_weights: list[tuple[int, float]],
        demand: dict[int, dict[int, float]],
        debug: bool,
    ) -> None:
        """apply temporal distribution matrix sourced to demand offsets.

        :param year: Value for `year`.
        :type year: int
        :param act_idx: Value for `act_idx`.
        :type act_idx: int
        :param product_index: Value for `product_index`.
        :type product_index: int
        :param parent_amount: Value for `parent_amount`.
        :type parent_amount: float
        :param offsets_and_weights: Value for `offsets_and_weights`.
        :type offsets_and_weights: list[tuple[int, float]]
        :param demand: Value for `demand`.
        :type demand: dict[int, dict[int, float]]
        :param debug: Value for `debug`.
        :type debug: bool"""
        if self.A is None:
            return

        if int(product_index) == int(act_idx):
            return

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

            child_amount = self._child_activity_amount(
                t=int(t_eff),
                product_index=int(product_index),
                parent_amount=float(parent_amount),
                exchange_value=exchange_value,
            )
            if child_amount == 0.0:
                continue

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
