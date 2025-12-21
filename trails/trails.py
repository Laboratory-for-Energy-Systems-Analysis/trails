# trails.py

from typing import Dict, List, Optional
from collections import defaultdict, deque

import numpy as np
import sparse

from tqdm import tqdm

from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package
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
        package,
        interpolate_annual: bool = True,
        value_dtype=np.float32,
        index_dtype=np.int32,
        debug: bool = False,
    ):
        self.package = package
        self.value_dtype = value_dtype
        self.index_dtype = index_dtype

        self.scenario_labels: List[str] = []
        self.scenario_index: Dict[str, int] = {}

        self.A: Optional[sparse.COO] = None
        self.B: Optional[sparse.COO] = None

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
        self.template_years_int = np.array([int(lbl) for lbl in self.template_labels], dtype=int)

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

            self.years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
            self.min_year = int(self.years_int.min())
            self.max_year = int(self.years_int.max())

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
                    or (tex0.scale_mode or "") != (tex1.scale_mode or "")
                ):
                    return tex0 if (year - y0) <= (y1 - year) else tex1

                w = (year - y0) / (y1 - y0)

                def interp_optional(v0, v1):
                    if v0 is None or v1 is None:
                        return v0 if (year - y0) <= (y1 - year) else v1
                    return float(v0) + (float(v1) - float(v0)) * w

                return TemporalExchange(
                    distribution=tex0.distribution,
                    loc=interp_optional(tex0.loc, tex1.loc),
                    scale=interp_optional(tex0.scale, tex1.scale),
                    offset_min=tex0.offset_min,
                    offset_max=tex0.offset_max,
                    scale_mode=tex0.scale_mode,
                    scale_base=interp_optional(tex0.scale_base, tex1.scale_base),
                    scale_rate=interp_optional(tex0.scale_rate, tex1.scale_rate),
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

    def get_A_for_scenario(self, label: str) -> sparse.COO:
        """Return the 2D A matrix (activity x product) for a given scenario label."""
        t = self.scenario_index[label]
        return self.A[t, :, :]

    def get_B_for_scenario(self, label: str) -> sparse.COO:
        """Return the 2D B matrix (activity x flow) for a given scenario label."""
        t = self.scenario_index[label]
        return self.B[t, :, :]

    def get_temporal_exchange(self, year: int, act_idx: int, prod_idx: int):
        """
        Return TemporalExchange for (year, act_idx, prod_idx), or None.
        Keys are based on the original scenario label strings.
        """
        return self._interpolate_temporal_exchange(
            year,
            act_idx,
            prod_idx,
            self.temporal_technosphere_exchanges,
        )

    def get_temporal_distribution(self, year: int, act_idx: int, prod_idx: int):
        """
        Return a TemporalDistribution object for (year, act_idx, prod_idx),
        or None if this exchange has no temporal metadata.
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

    def expand_temporal_exchanges(
            self,
            year: int,
            act_idx: int,
            amount: float = 1.0,
            *,
            use_temporal_distributions: bool = True,
            debug: bool = False,
    ):
        """
        Expand one activity-year demand into temporally distributed multi-year
        demands for its *direct* exchanges.

        If use_temporal_distributions=False, treat all exchanges as occurring in
        the (mapped) scenario year (i.e. no temporal shifting).
        """
        demand: dict[int, dict[int, float]] = {}

        scenario_year = self._map_year_to_scenario_year(year)
        scenario_label = str(scenario_year)
        t = self.scenario_index[scenario_label]

        if debug:
            logger.info(
                "expand_tech: year=%d scenario_year=%d t=%d act=%d amount=%g",
                int(year), int(scenario_year), int(t), int(act_idx), float(amount)
            )

            logger.debug(
                "expand_temporal_exchanges: year=%d act=%d amount=%g scenario_year=%d use_td=%s",
                year, act_idx, amount, scenario_year, use_temporal_distributions
            )

        A_row = self.A[t, act_idx, :]
        if A_row.nnz == 0:
            if debug:
                logger.debug(
                    "expand_temporal_exchanges: EMPTY A_row (year=%d mapped=%d t=%d act=%d) -> no children",
                    year, scenario_year, t, act_idx
                )
            return demand


        product_indices = A_row.coords[0]
        values = A_row.data

        if debug:
            logger.debug(
                "expand_temporal_exchanges: A_row.nnz=%d (before skipping production/diag outputs)",
                int(A_row.nnz),
            )

        def _add_to_demand(target_year: int, product_index: int, exchange_amount: float) -> None:
            demand.setdefault(target_year, {})
            demand[target_year][product_index] = demand[target_year].get(product_index, 0.0) + exchange_amount

        def _child_amount(parent_amount: float, exchange_value: float) -> float:
            if exchange_value < 0.0:
                return parent_amount * (-exchange_value)
            return parent_amount * exchange_value

        for product_index, exchange_value in zip(product_indices, values):
            exchange_value = float(exchange_value)
            product_index = int(product_index)

            # Always skip the canonical production exchange, if present in your data
            # (common convention: A[act, act] = 1)
            if product_index == act_idx and abs(exchange_value) == 1.0:
                continue

            # Now compute the propagated requirement amount:
            # - Inputs are negative -> convert to positive requirement magnitude
            # - Off-diagonal positive entries (e.g., waste treatment services) are kept as-is
            child_amount = _child_amount(amount, exchange_value)

            tex = self._get_tech_temporal_exchange(year, act_idx, product_index)

            if (tex is None) or (not use_temporal_distributions):
                y_eff = scenario_year
                _add_to_demand(y_eff, product_index, child_amount)
                continue

            if debug:
                logger.debug(
                    "expand_temporal_exchanges: applying TD for (year=%d act=%d prod=%d) child_amt=%g",
                    year, act_idx, product_index, child_amount
                )

            td = TemporalDistribution(tex)

            offsets_and_weights = list(td.iter_offsets_and_weights(debug=debug))
            if not offsets_and_weights:
                if debug:
                    logger.warning(
                        "expand_temporal_exchanges: TD produced no offsets/weights for (year=%d act=%d prod=%d) -> dropping exchange",
                        year, act_idx, product_index
                    )
            else:
                if debug:
                    logger.debug(
                        "expand_temporal_exchanges: TD offsets=%s (sum_w=%g)",
                        [p[0] for p in offsets_and_weights],
                        float(sum(p[1] for p in offsets_and_weights)),
                    )

            for offset, weight in offsets_and_weights:
                raw_year = year + offset
                y_eff = self._map_year_to_scenario_year(raw_year)

                factor = td.scale_factor(offset)
                _add_to_demand(
                    y_eff,
                    product_index,
                    child_amount * float(weight) * float(factor),
                )

        total_children = sum(len(m) for m in demand.values())
        if debug:
            logger.debug(
                "expand_temporal_exchanges: produced years=%d total_children=%d",
                len(demand), total_children
            )

        return demand

    def _get_tech_temporal_exchange(self, year: int, act_idx: int, prod_idx: int) -> Optional[TemporalExchange]:
        """
        For technosphere TD metadata, do NOT interpolate across years.
        Instead, map to the nearest template year and do a direct lookup.
        This prevents TD metadata from 'bleeding' into years where it wasn't specified
        and makes TD availability stable across interpolated scenario years.
        """
        if not self.temporal_technosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        return self.temporal_technosphere_exchanges.get((str(y_tpl), int(act_idx), int(prod_idx)))

    def _get_bio_temporal_exchange(self, year: int, act_idx: int, flow_idx: int) -> Optional[TemporalExchange]:
        """
        For biosphere TD metadata, do NOT interpolate across years.
        Instead, map to the nearest template year and do a direct lookup.
        This prevents TD metadata from 'bleeding' into years where it wasn't specified.
        """
        if not self.temporal_biosphere_exchanges:
            return None

        y_tpl = self._map_year_to_template_year(year)
        return self.temporal_biosphere_exchanges.get((str(y_tpl), int(act_idx), int(flow_idx)))

    def accumulate_temporalized_biosphere_inventory(
            self,
            base_year: int,
            supply_by_activity: Dict[int, float],
            inventory_by_year: Dict[int, np.ndarray],
            *,
            min_amount: float = 0.0,
            use_temporal_distributions: bool = True,
            debug: bool = False,
    ) -> None:
        """
        Accumulate temporally shifted biosphere emissions resulting from a
        solved technosphere supply vector for a given calendar year.

        supply_by_activity maps Trails activity indices -> supply.
        """
        if debug:
            logger.info(
                "accumulate_bio: base_year=%d acts_in=%d min_amount=%g use_td=%s",
                int(base_year), len(supply_by_activity), float(min_amount), bool(use_temporal_distributions)
            )

        if self.B is None:
            if debug:
                logger.warning("accumulate_bio: B is None -> nothing to accumulate")
            return

        scenario_year = self._map_year_to_scenario_year(base_year)
        scenario_label = str(scenario_year)
        if scenario_label not in self.scenario_index:
            if debug:
                logger.error(
                    "accumulate_bio: scenario_label=%s not in scenario_index (base_year=%d) -> abort",
                    scenario_label, int(base_year)
                )
            return
        t = self.scenario_index[scenario_label]

        if debug:
            logger.info(
                "accumulate_bio: base_year=%d scenario_year=%d t=%d",
                int(base_year), int(scenario_year), int(t)
            )

        B_t = self.B[t, :, :]
        B_t_nnz = int(getattr(B_t, "nnz", 0))
        if debug:
            logger.info(
                "accumulate_bio: B slice t=%d nnz=%d shape=%s",
                int(t), B_t_nnz, getattr(B_t, "shape", None)
            )
        if B_t_nnz == 0:
            if debug:
                logger.warning("accumulate_bio: B_t.nnz==0 -> nothing to accumulate")
            return

        n_flows = int(self.B.shape[2])

        act_indices = B_t.coords[0].astype(int)
        flow_indices = B_t.coords[1].astype(int)
        values = B_t.data

        # Use template year for temporal metadata (distributions defined on original years)
        template_year = self._map_year_to_template_year(base_year)

        for act_idx, flow_idx, value in zip(act_indices, flow_indices, values):
            act_idx = int(act_idx)
            flow_idx = int(flow_idx)

            supply_amt = float(supply_by_activity.get(act_idx, 0.0))

            if supply_amt == 0.0:
                continue

            # IMPORTANT: biosphere sign is not flipped; value is already the correct sign
            scaled = supply_amt * float(value)
            if scaled == 0.0:
                continue

            if min_amount and abs(scaled) < min_amount:
                continue

            # Fetch temporal metadata using template year (not base_year)
            tex = self._get_bio_temporal_exchange(
                base_year,
                act_idx,
                flow_idx,
            )

            if (tex is None) or (not use_temporal_distributions):
                # Anchor inventory to the same scenario year used for B slice
                y_eff = int(scenario_year)
                inventory_by_year.setdefault(
                    y_eff,
                    np.zeros(n_flows, dtype=self.value_dtype),
                )
                inventory_by_year[y_eff][flow_idx] += scaled
                continue

            td = TemporalDistribution(tex)

            # Anchor offsets to scenario_year (because the coefficient came from B[t])
            any_pair = False
            for offset, weight in td.iter_offsets_and_weights(debug=debug):
                any_pair = True

                raw_year = int(scenario_year + int(offset))
                y_eff = int(self._map_year_to_scenario_year(raw_year))

                factor = float(td.scale_factor(offset))

                inventory_by_year.setdefault(
                    y_eff,
                    np.zeros(n_flows, dtype=self.value_dtype),
                )
                inventory_by_year[y_eff][flow_idx] += scaled * float(weight) * factor

            if not any_pair:
                if debug:
                    logger.warning(
                        "accumulate_bio: TD produced no offsets/weights (template_year=%d scenario_year=%d act=%d flow=%d) -> dropped",
                        int(template_year), int(scenario_year), act_idx, flow_idx
                    )

    def _map_year_to_available(self, year: int) -> int:
        """
        Backwards-compatible alias for _map_year_to_scenario_year.
        """
        return self._map_year_to_scenario_year(year)


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
    ):
        """
        Traverse the temporal-technosphere graph starting from (start_year, start_act_idx).

        Progress:
          - Prefer an empirical total estimate based on max_depth (DEPTH_TOTALS).
          - If not available, fall back to a short warm-up branching estimate.
        """

        if debug:
            logger.info(
                "temporal_traversal start: start_year=%d start_act=%d amount=%g max_depth=%d min_amount=%g use_td=%s",
                start_year, start_act_idx, amount, max_depth, min_amount, use_temporal_distributions
            )

        # ------------------------------------------------------------------
        # 1) Empirical totals by depth (YOU populate these)
        #    Values should represent the *typical nodes_processed* for that depth.
        # ------------------------------------------------------------------
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

        # Conservative multiplier so the bar doesn't finish early
        EMPIRICAL_SAFETY_FACTOR = 1.05  # tune: 1.05–1.30

        # ------------------------------------------------------------------
        # 2) Warm-up fallback params (only used if DEPTH_TOTALS lacks max_depth)
        # ------------------------------------------------------------------
        WARMUP_LIMIT = 1000
        BRANCHING_PERCENTILE = 95.0
        BRANCHING_SAFETY_FACTOR = 1.2

        def estimate_total_from_branching(branching_samples):
            if not branching_samples:
                return 1
            s = sorted(branching_samples)
            k = int((BRANCHING_PERCENTILE / 100.0) * (len(s) - 1))
            b = max(1.0, float(s[k]) * BRANCHING_SAFETY_FACTOR)
            if abs(b - 1.0) < 1e-9:
                return max_depth + 1
            return int((b ** (max_depth + 1) - 1.0) / (b - 1.0))

        def estimate_total_from_depth():
            # Exact depth available
            if max_depth in DEPTH_TOTALS:
                return int(max(1, DEPTH_TOTALS[max_depth] * EMPIRICAL_SAFETY_FACTOR))

            # If you have surrounding depths, interpolate in log-space (often more stable)
            depths = sorted(DEPTH_TOTALS.keys())
            if len(depths) >= 2:
                lo = max([d for d in depths if d < max_depth], default=None)
                hi = min([d for d in depths if d > max_depth], default=None)
                if lo is not None and hi is not None and DEPTH_TOTALS[lo] > 0 and DEPTH_TOTALS[hi] > 0:
                    import math
                    y0 = math.log(DEPTH_TOTALS[lo])
                    y1 = math.log(DEPTH_TOTALS[hi])
                    t = (max_depth - lo) / (hi - lo)
                    est = math.exp(y0 + t * (y1 - y0))
                    return int(max(1, est * EMPIRICAL_SAFETY_FACTOR))

            # Otherwise, no empirical estimate
            return None

        # ------------------------------------------------------------------
        # Progress bar setup
        # ------------------------------------------------------------------
        pbar = None
        branching_samples = []
        nodes_processed = 0

        total_est = None
        if show_progress:
            total_est = estimate_total_from_depth()
            try:
                if total_est is None:
                    # Indeterminate until warm-up can estimate
                    pbar = tqdm(total=None, desc="Temporal traversal", unit="node", dynamic_ncols=True)
                else:
                    pbar = tqdm(total=total_est, desc="Temporal traversal", unit="node", dynamic_ncols=True)
            except Exception:
                pbar = None

        # ------------------------------------------------------------------
        # Traversal state
        # ------------------------------------------------------------------
        queue = deque()
        queue.append((int(start_year), int(start_act_idx), float(amount), 0, (), None))

        frontier_total = defaultdict(float)  # (year, act) -> amt
        provenance_roots = defaultdict(lambda: defaultdict(float))  # (year, act) -> {root_act: amt}
        direct_bio_total = defaultdict(
            float)  # (year, act) -> amt (nodes with direct biosphere we do NOT want to solve)
        direct_bio_roots = defaultdict(lambda: defaultdict(float))  # same but by root, if provenance requested

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

            if pbar is not None:
                pbar.update(1)

            # Map to scenario year for "has direct biosphere" test (fast cutoff logic)
            scenario_year = self._map_year_to_scenario_year(year)
            label = str(scenario_year)

            has_direct_bio = False
            if label in self.scenario_index and (self.B is not None):
                t = self.scenario_index[label]
                key = (scenario_year, act)
                if key in bio_cache:
                    has_direct_bio = bio_cache[key]
                else:
                    has_direct_bio = (self.B[t, act, :].nnz > 0)
                    bio_cache[key] = has_direct_bio

            # Helper: record a node into frontier + provenance
            def _record_frontier(y: int, a: int, x: float, r: Optional[int]):
                frontier_total[(int(y), int(a))] += float(x)
                if return_provenance and (r is not None):
                    provenance_roots[(int(y), int(a))][int(r)] += float(x)

            def _record_direct_bio(y: int, a: int, x: float, r: Optional[int]):
                direct_bio_total[(int(y), int(a))] += float(x)
                if return_provenance and (r is not None):
                    direct_bio_roots[(int(y), int(a))][int(r)] += float(x)

            # Stop expanding at max_depth
            if depth >= max_depth:
                _record_frontier(year, act, amt, root_act)
                continue

            # Expand this node
            child_demands = self.expand_temporal_exchanges(
                year=year,
                act_idx=act,
                amount=amt,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
            )

            # Leaf: record it
            if not child_demands:
                _record_frontier(year, act, amt, root_act)
                continue

            # IMPORTANT BEHAVIOR:
            # If this node has direct biosphere flows, we record it as part of the frontier.
            # This matches your existing “solve nodes with direct biosphere” design.
            # (It is NOT a full “score technosphere exchange at its own year” algorithm.)
            if has_direct_bio and depth > 0:
                _record_direct_bio(year, act, amt, root_act)

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

                    queue.append((child_year, child_act, child_amt, depth + 1, child_path, child_root))

        if pbar is not None:
            try:
                pbar.close()
            except Exception:
                pass

        # Normalize provenance to plain dicts
        if return_provenance:
            provenance_roots = {k: dict(v) for k, v in provenance_roots.items()}
            direct_bio_roots = {k: dict(v) for k, v in direct_bio_roots.items()}
            return dict(frontier_total), provenance_roots, dict(direct_bio_total), direct_bio_roots

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
                raise ValueError(f"Frontier key must be a tuple (year, act). Got {type(key)}: {key}")

            if len(key) != 2:
                raise ValueError(f"Frontier key must be (year, act). Got len={len(key)}: {key}")

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
        """
        Traverse the temporal-technosphere graph starting from
        (start_year, start_act_idx) and record edges by depth.

        Returns
        -------
        edges_by_depth : dict[int, dict[((int, int), (int, int)), float]]
            {depth: {((year_from, act_from), (year_to, act_to)): amount, ...}, ...}

        - 'depth' is the depth of the *parent* node.
        - 'amount' is the flow leaving (year_from, act_from) towards (year_to, act_to)
          at that depth, starting from a functional unit of `amount` at the root.
        """
        queue = deque()
        # root_act is the first-level child activity (not year), or None for the FU itself
        queue.append((start_year, start_act_idx, float(amount), 0, (), None))

        edges_by_depth: dict[int, dict[tuple[tuple[int, int], tuple[int, int]], float]] = (
            defaultdict(lambda: defaultdict(float))
        )

        while queue:
            year, act, amt, depth = queue.popleft()

            if abs(amt) < min_amount:
                continue

            if depth >= max_depth:
                # We still record the node as existing at this depth,
                # but we don't expand children anymore.
                continue

            # Expand this node
            child_demands = self.expand_temporal_exchanges(
                year=year, act_idx=act, amount=amt
            )

            if not child_demands:
                continue

            parent_node = (year, act)

            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue
                    child_node = (child_year, int(child_act))
                    # edge stored at depth of the parent
                    edges_by_depth[depth][(parent_node, child_node)] += float(child_amt)

                    # enqueue child at next depth
                    queue.append((child_year, child_act, child_amt, depth + 1))

        # make inner dicts plain dicts
        return {d: dict(edges) for d, edges in edges_by_depth.items()}
