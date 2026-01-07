# trails.py

from typing import Any, Dict, List, Optional
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

        self.scenario_labels: List[str] = []
        self.scenario_index: Dict[str, int] = {}

        self.A: Optional[sparse.COO] = None
        self.B: Optional[sparse.COO] = None
        self.inventory: Optional[xr.DataArray] = None
        self.characterized_inventory: Optional[xr.DataArray] = None
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

    def reset_inventory(self) -> None:
        """Initialize inventory builders for sparse 3D inventory storage."""
        years = np.array(self.years_int, dtype=int)
        self._inventory_years = years
        self._inventory_year_index = {int(y): int(i) for i, y in enumerate(years)}
        self._inventory_coords = [[], [], []]
        self._inventory_data = []
        self.inventory = None
        self.characterized_inventory = None
        self.provenance = None

    def _append_inventory_entries(
        self, act_idx: int, year: int, flows: np.ndarray, values: np.ndarray
    ) -> None:
        if self._inventory_coords is None or self._inventory_data is None:
            raise RuntimeError(
                "Inventory builders not initialized. Call reset_inventory()."
            )

        year_idx = self._inventory_year_index.get(int(year))
        if year_idx is None:
            return

        flows = np.asarray(flows, dtype=np.int64)
        values = np.asarray(values, dtype=np.float64)
        if flows.size == 0:
            return

        mask = values != 0.0
        if not np.any(mask):
            return

        flows = flows[mask]
        values = values[mask]

        self._inventory_coords[0].append(
            np.full(flows.size, int(act_idx), dtype=np.int64)
        )
        self._inventory_coords[1].append(flows.astype(np.int64, copy=False))
        self._inventory_coords[2].append(
            np.full(flows.size, int(year_idx), dtype=np.int64)
        )
        self._inventory_data.append(values.astype(self.value_dtype, copy=False))

    def finalize_inventory(self) -> xr.DataArray:
        """Finalize and store sparse inventory as a 3D xarray."""
        if self.A is None or self.B is None:
            raise ValueError("Cannot finalize inventory: A or B is None.")

        n_activities = int(self.A.shape[1])
        n_flows = int(self.B.shape[2])
        years = self._inventory_years
        if years is None:
            raise RuntimeError(
                "Inventory years not initialized. Call reset_inventory()."
            )

        if self._inventory_coords is None or self._inventory_data is None:
            raise RuntimeError(
                "Inventory builders not initialized. Call reset_inventory()."
            )

        if self._inventory_data:
            coords = np.vstack(
                [np.concatenate(part) for part in self._inventory_coords]
            )
            data = np.concatenate(self._inventory_data)
            inv = sparse.COO(coords, data, shape=(n_activities, n_flows, len(years)))
        else:
            inv = sparse.COO.zeros(
                (n_activities, n_flows, len(years)), dtype=self.value_dtype
            )

        self.inventory = xr.DataArray(
            inv,
            dims=("activity", "flow", "year"),
            coords={
                "activity": np.arange(n_activities, dtype=int),
                "flow": np.arange(n_flows, dtype=int),
                "year": years,
            },
        )
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
            y_eff = self._map_year_to_scenario_year(raw_year)

            self._add_demand_entry(
                demand,
                y_eff,
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
        """Run temporal LCA using this Trails instance."""
        from .lca import lca as lca_fn

        return lca_fn(self, *args, **kwargs)

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

        context = self._get_scenario_context(year)
        if context is None:
            return demand
        scenario_year, scenario_label, t = context

        if debug:
            logger.info(
                "expand_tech: year=%d scenario_year=%d t=%d act=%d amount=%g",
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
            product_index = int(product_index)
            exchange_value = float(exchange_value)

            if exchange_value == 0.0:
                continue

            # Skip canonical production exchange (A[act, act] = 1)
            if product_index == act_idx and abs(exchange_value) == 1.0:
                continue

            # Fetch TD metadata (template-year lookup; stable across interpolation)
            tex = self._get_tech_temporal_exchange(year, act_idx, product_index)

            # ------------------------------------------------------------------
            # No temporal distribution (or disabled): status quo
            # ------------------------------------------------------------------
            if (tex is None) or (not use_temporal_distributions):
                child_amount = self._child_amount(amount, exchange_value)
                if child_amount != 0.0:
                    self._add_demand_entry(
                        demand, int(scenario_year), product_index, float(child_amount)
                    )
                continue

            amount_source = getattr(tex, "amount_source", "port")

            # ------------------------------------------------------------------
            # TD + matrix-sourced magnitude: read A at each pulse year
            # ------------------------------------------------------------------
            if amount_source == "matrix":
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

    def accumulate_temporalized_biosphere_inventory(
        self,
        base_year: int,
        supply_by_activity: Dict[int, float],
        *,
        min_amount: float = 0.0,
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
        """
        # ---------------------------
        # Early exits / slice resolve
        # ---------------------------
        biosphere_slice = self._get_biosphere_slice(base_year, debug)
        if biosphere_slice is None:
            return
        scenario_year, t, B_t, n_flows = biosphere_slice

        if not supply_by_activity:
            return

        # ---------------------------
        # Localize hot attrs / methods
        # ---------------------------
        value_dtype = self.value_dtype
        scenario_index_get = self.scenario_index.get
        map_year_to_scenario = self._map_year_to_scenario_year

        base_scenario_year = int(scenario_year)
        min_amt = float(min_amount) if min_amount else 0.0

        # ---------------------------
        # Temporal metadata context
        # ---------------------------
        bio_td = (
            self.temporal_biosphere_exchanges if use_temporal_distributions else None
        )
        if bio_td:
            tpl_label = str(self._map_year_to_template_year(base_year))
            bio_td_get = bio_td.get
        else:
            tpl_label = None
            bio_td_get = None  # type: ignore[assignment]

        # ---------------------------
        # Per-call caches
        # ---------------------------
        year_map_cache: dict[int, int] = {}
        t_eff_cache: dict[int, int | None] = {}

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

        # ---------------------------
        # Global caches on self (persist across calls)
        # ---------------------------
        # Cache TD row classification per (tpl_label, t, act):
        #   (no_td_idx, port_groups_idx, matrix_entries)
        # where:
        #   no_td_idx: np.ndarray[intp] | None  positions in FULL row arrays
        #   port_groups_idx: dict[td_key -> np.ndarray[intp]] positions in FULL row arrays
        #   matrix_entries: list[(pos:int, tex:TemporalExchange)] positions in FULL row arrays
        if not hasattr(self, "_bio_td_row_cache"):
            self._bio_td_row_cache = {}  # type: ignore[attr-defined]
        row_td_cache = self._bio_td_row_cache  # type: ignore[attr-defined]

        # Cache pulses per td_key across calls
        if not hasattr(self, "_td_pulse_cache"):
            self._td_pulse_cache = {}  # type: ignore[attr-defined]
        pulse_cache = self._td_pulse_cache  # type: ignore[attr-defined]

        # ---------------------------
        # Fast row-structure cache (per t)
        # ---------------------------
        if not hasattr(self, "_B_row_cache"):
            self._B_row_cache = {}  # type: ignore[attr-defined]

        row_cache = self._B_row_cache  # type: ignore[attr-defined]
        cached = row_cache.get(int(t))

        if cached is None:
            act_coords = B_t.coords[0].astype(np.int32, copy=False)
            flow_coords = B_t.coords[1].astype(np.int32, copy=False)
            data = (
                B_t.data.astype(np.float32, copy=False)
                if value_dtype == np.float32
                else B_t.data
            )

            n_acts = int(B_t.shape[0])
            nnz = int(getattr(B_t, "nnz", 0))
            if nnz == 0:
                row_cache[int(t)] = (
                    np.zeros(n_acts + 1, dtype=np.int64),
                    flow_coords,
                    data,
                )
                return

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

        # ---------------------------
        # Main accumulation
        # ---------------------------
        for act_idx, supply_amt in supply_by_activity.items():
            supply_amt = float(supply_amt)
            if supply_amt == 0.0:
                continue

            a = int(act_idx)
            if a < 0 or a + 1 >= len(row_ptr):
                continue

            start = int(row_ptr[a])
            end = int(row_ptr[a + 1])
            if start == end:
                continue

            # FULL row arrays (positions refer to these)
            flows_full = flow_sorted[start:end].astype(
                np.intp, copy=False
            )  # for np.add.at
            vals_full = data_sorted[start:end]

            # Precompute scaled contributions for this activity row (anchor-year)
            scaled_full = supply_amt * vals_full.astype(np.float64, copy=False)

            # Optional min_amount filter mask on FULL arrays (so cached indices still apply)
            if min_amt:
                keep_full = np.abs(scaled_full) >= min_amt
                if not keep_full.any():
                    continue
            else:
                keep_full = None

            # ---------------------------
            # No TD at all => fast vectorized anchor add
            # ---------------------------
            if not bio_td:
                if keep_full is None:
                    self._append_inventory_entries(
                        a, base_scenario_year, flows_full, scaled_full
                    )
                else:
                    self._append_inventory_entries(
                        a,
                        base_scenario_year,
                        flows_full[keep_full],
                        scaled_full[keep_full],
                    )
                continue

            # ---------------------------
            # TD enabled: use cached TD classification per (tpl_label, t, act)
            # ---------------------------
            cache_key = (tpl_label, int(t), int(a))
            td_struct = row_td_cache.get(cache_key)

            if td_struct is None:
                no_td_pos: list[int] = []
                port_groups_pos: dict[tuple, list[int]] = {}
                matrix_entries: list[tuple[int, TemporalExchange]] = []

                # Classify on FULL row once
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
                port_groups_idx = {
                    k: np.array(v, dtype=np.intp) for k, v in port_groups_pos.items()
                }

                td_struct = (no_td_idx, port_groups_idx, matrix_entries)
                row_td_cache[cache_key] = td_struct

            no_td_idx, port_groups_idx, matrix_entries = td_struct

            # ---------------------------
            # 1) No TD: batch anchor add (using cached indices)
            # ---------------------------
            if no_td_idx is not None:
                if keep_full is None:
                    idx = no_td_idx
                else:
                    idx = no_td_idx[keep_full[no_td_idx]]
                if idx.size:
                    self._append_inventory_entries(
                        a, base_scenario_year, flows_full[idx], scaled_full[idx]
                    )

            # ---------------------------
            # 2) Ported TD: grouped vector math + add.at per pulse (using cached indices)
            # ---------------------------
            if port_groups_idx:
                for k, idx_full in port_groups_idx.items():
                    # Apply min_amount mask after selecting idx_full
                    if keep_full is None:
                        idx = idx_full
                    else:
                        idx = idx_full[keep_full[idx_full]]
                        if idx.size == 0:
                            continue

                    f_arr = flows_full[idx]
                    s_arr = scaled_full[idx]

                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        # Need a representative tex to compute pulses.
                        # We can look up tex for the first flow in this group.
                        f0 = int(f_arr[0])
                        tex0 = bio_td_get((tpl_label, a, f0))  # type: ignore[misc]
                        if tex0 is None:
                            continue
                        pulses = [
                            (int(o), float(w))
                            for o, w in TemporalDistribution(
                                tex0
                            ).iter_offsets_and_weights(debug=False)
                        ]
                        pulse_cache[k] = pulses

                    for offset, weight in pulses:
                        if weight == 0.0:
                            continue

                        y_eff = map_year_cached(base_scenario_year + int(offset))
                        contrib = s_arr * float(weight)

                        if min_amt:
                            m = np.abs(contrib) >= min_amt
                            if not m.any():
                                continue
                            f_use = f_arr[m]
                            c_use = contrib[m]
                        else:
                            f_use = f_arr
                            c_use = contrib

                        self._append_inventory_entries(a, y_eff, f_use, c_use)

            # ---------------------------
            # 3) Matrix-sourced TD: keep semantics (year-dependent values)
            # ---------------------------
            if matrix_entries:
                for p, tex in matrix_entries:
                    if keep_full is not None and not keep_full[p]:
                        continue

                    f = int(flows_full[p])

                    k = td_key(tex)
                    pulses = pulse_cache.get(k)
                    if pulses is None:
                        pulses = [
                            (int(o), float(w))
                            for o, w in TemporalDistribution(
                                tex
                            ).iter_offsets_and_weights(debug=False)
                        ]
                        pulse_cache[k] = pulses

                    for offset, weight in pulses:
                        if weight == 0.0:
                            continue

                        y_eff = map_year_cached(base_scenario_year + int(offset))

                        t_eff = t_eff_cache.get(y_eff)
                        if t_eff is None and y_eff not in t_eff_cache:
                            t_eff = scenario_index_get(str(y_eff))
                            t_eff_cache[y_eff] = t_eff
                        if t_eff is None:
                            continue

                        value_eff = float(self.B[int(t_eff), a, f])
                        if value_eff == 0.0:
                            continue

                        contrib = supply_amt * value_eff * float(weight)
                        if contrib == 0.0:
                            continue
                        if min_amt and abs(contrib) < min_amt:
                            continue

                        self._append_inventory_entries(
                            a,
                            y_eff,
                            np.array([f], dtype=np.int64),
                            np.array([contrib]),
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
            if total_est is None:
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

            if abs(amt) < min_amount:
                continue

            _pbar_step()

            # Map to scenario year for "has direct biosphere" test (fast cutoff logic)
            scenario_year = self._map_year_to_scenario_year(year)
            has_direct_bio = self._has_direct_biosphere(scenario_year, act, bio_cache)

            # Helper: record a node into frontier + provenance
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

            # Expand this node
            child_demands = self.expand_temporal_exchanges(
                year=year,
                act_idx=act,
                amount=amt,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
            )

            # --------------------------------------------------------------
            # Warm-up: if tqdm started indeterminate (total=None),
            # estimate a total after WARMUP_LIMIT processed nodes using
            # observed branching and then set pbar.total.
            # --------------------------------------------------------------
            if show_progress and pbar is not None and pbar.total is None:
                # Branching sample = number of children edges we would enqueue
                # for this node (after min_amount filtering, consistent with traversal).
                if nodes_processed <= WARMUP_LIMIT:
                    # Count children that would actually be enqueued
                    cnt = 0
                    if child_demands:
                        for _cy, _mapping in child_demands.items():
                            for _ca, _camt in _mapping.items():
                                if abs(float(_camt)) >= float(min_amount):
                                    cnt += 1
                    branching_samples.append(cnt)

                # Once warm-up complete, set a total estimate
                if nodes_processed == WARMUP_LIMIT:
                    est = estimate_total_from_branching(branching_samples)
                    # Keep a bit conservative so it doesn't finish early
                    est = int(max(est, pbar.n + 1) * EMPIRICAL_SAFETY_FACTOR)
                    pbar.total = est
                    pbar.refresh()

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

            # IMPORTANT BEHAVIOR:
            # If this node has direct biosphere flows, we record it as part of the frontier.
            # This matches your existing “solve nodes with direct biosphere” design.
            # (It is NOT a full “score technosphere exchange at its own year” algorithm.)
            if has_direct_bio and depth > 0:
                self._record_direct_bio(
                    direct_bio_total,
                    direct_bio_roots,
                    year,
                    act,
                    amt,
                    root_act,
                    return_provenance,
                )

            # Enqueue children
            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    child_amt = float(child_amt)
                    if abs(child_amt) < min_amount:
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
        :param min_amount: Minimum magnitude to include.
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

            if abs(amt) < min_amount:
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
                    if abs(child_amt) < min_amount:
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

            self._add_demand_entry(
                demand, y_eff, int(product_index), weighted_child_amount
            )
