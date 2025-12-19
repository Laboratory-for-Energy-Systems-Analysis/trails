# trails.py

from typing import Dict, List, Optional
from collections import defaultdict, deque

import numpy as np
import sparse
import pyprind

from tqdm import tqdm

from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package
)

from .temporal_distributions import TemporalDistribution

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
        )

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
            )

            self.years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
            self.min_year = int(self.years_int.min())
            self.max_year = int(self.years_int.max())

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------

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
        label = str(year)
        return self.temporal_technosphere_exchanges.get((label, act_idx, prod_idx))

    def get_temporal_distribution(self, year: int, act_idx: int, prod_idx: int):
        """
        Return a TemporalDistribution object for (year, act_idx, prod_idx),
        or None if this exchange has no temporal metadata.
        """
        label = str(year)
        tex = self.temporal_technosphere_exchanges.get((label, act_idx, prod_idx))
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

        logger.info(
            "expand_tech: year=%d scenario_year=%d t=%d act=%d amount=%g",
            int(year), int(scenario_year), int(t), int(act_idx), float(amount)
        )

        template_year = self._map_year_to_template_year(year)
        template_label = str(template_year)

        logger.debug(
            "expand_temporal_exchanges: year=%d act=%d amount=%g scenario_year=%d template_year=%d use_td=%s",
            year, act_idx, amount, scenario_year, template_year, use_temporal_distributions
        )

        A_row = self.A[t, act_idx, :]
        if A_row.nnz == 0:
            logger.debug(
                "expand_temporal_exchanges: EMPTY A_row (year=%d mapped=%d t=%d act=%d) -> no children",
                year, scenario_year, t, act_idx
            )
            return demand


        prod_indices = A_row.coords[0]
        values = A_row.data

        logger.debug(
            "expand_temporal_exchanges: A_row.nnz=%d (before skipping production/diag outputs)",
            int(A_row.nnz),
        )

        for prod_idx, value in zip(prod_indices, values):
            value = float(value)
            prod_idx = int(prod_idx)

            # Always skip the canonical production exchange, if present in your data
            # (common convention: A[act, act] = 1)
            if prod_idx == act_idx and abs(value) == 1.0:
                continue

            # Now compute the propagated requirement amount:
            # - Inputs are negative -> convert to positive requirement magnitude
            # - Off-diagonal positive entries (e.g., waste treatment services) are kept as-is
            if value < 0.0:
                child_amt = amount * (-value)  # positive requirement magnitude
            else:
                child_amt = amount * (value)  # keep positive sign (supply-driven service)

            tex = self.temporal_technosphere_exchanges.get((template_label, act_idx, prod_idx))

            if (tex is None) or (not use_temporal_distributions):
                y_eff = scenario_year
                demand.setdefault(y_eff, {})
                demand[y_eff][prod_idx] = demand[y_eff].get(prod_idx, 0.0) + child_amt
                continue

            logger.debug(
                "expand_temporal_exchanges: applying TD for (template=%s act=%d prod=%d) child_amt=%g",
                template_label, act_idx, prod_idx, child_amt
            )

            td = TemporalDistribution(tex)

            pairs = list(td.iter_offsets_and_weights())
            if not pairs:
                logger.warning(
                    "expand_temporal_exchanges: TD produced no offsets/weights for (template=%s act=%d prod=%d) -> dropping exchange",
                    template_label, act_idx, prod_idx
                )
            else:
                logger.debug(
                    "expand_temporal_exchanges: TD offsets=%s (sum_w=%g)",
                    [p[0] for p in pairs],
                    float(sum(p[1] for p in pairs)),
                )

            for offset, weight in td.iter_offsets_and_weights():
                raw_year = year + offset
                y_eff = self._map_year_to_scenario_year(raw_year)

                factor = td.scale_factor(offset)
                demand.setdefault(y_eff, {})
                demand[y_eff][prod_idx] = (
                        demand[y_eff].get(prod_idx, 0.0)
                        + child_amt * float(weight) * float(factor)
                )

        total_children = sum(len(m) for m in demand.values())
        logger.debug(
            "expand_temporal_exchanges: produced years=%d total_children=%d",
            len(demand), total_children
        )

        return demand

    def accumulate_temporalized_biosphere_inventory(
            self,
            base_year: int,
            supply_by_activity: Dict[int, float],
            inventory_by_year: Dict[int, np.ndarray],
            *,
            min_amount: float = 0.0,
            use_temporal_distributions: bool = True,
    ) -> None:
        """
        Accumulate temporally shifted biosphere emissions resulting from a
        solved technosphere supply vector for a given calendar year.

        supply_by_activity maps Trails activity indices -> supply.
        """
        # -----------------------------
        # (A) Entry diagnostics
        # -----------------------------
        logger.info(
            "accumulate_bio: base_year=%d acts_in=%d min_amount=%g use_td=%s",
            int(base_year), len(supply_by_activity), float(min_amount), bool(use_temporal_distributions)
        )

        if self.B is None:
            logger.warning("accumulate_bio: B is None -> nothing to accumulate")
            return

        # -----------------------------
        # (B) Scenario/template mapping diagnostics
        # -----------------------------
        scenario_year = self._map_year_to_scenario_year(base_year)
        scenario_label = str(scenario_year)
        if scenario_label not in self.scenario_index:
            logger.error(
                "accumulate_bio: scenario_label=%s not in scenario_index (base_year=%d) -> abort",
                scenario_label, int(base_year)
            )
            return
        t = self.scenario_index[scenario_label]

        template_year = self._map_year_to_template_year(base_year)
        template_label = str(template_year)

        logger.info(
            "accumulate_bio: base_year=%d scenario_year=%d t=%d template_year=%d",
            int(base_year), int(scenario_year), int(t), int(template_year)
        )

        # -----------------------------
        # (C) B slice diagnostics
        # -----------------------------
        B_t = self.B[t, :, :]
        B_t_nnz = int(getattr(B_t, "nnz", 0))
        logger.info(
            "accumulate_bio: B slice t=%d nnz=%d shape=%s",
            int(t), B_t_nnz, getattr(B_t, "shape", None)
        )
        if B_t_nnz == 0:
            logger.warning("accumulate_bio: B_t.nnz==0 -> nothing to accumulate")
            return

        n_flows = int(self.B.shape[2])

        act_indices = B_t.coords[0].astype(int)
        flow_indices = B_t.coords[1].astype(int)
        values = B_t.data

        # -----------------------------
        # (D) Main accumulation with “why zero?” instrumentation
        # -----------------------------
        n_triplets = 0
        n_supply_nonzero = 0
        n_scaled_nonzero = 0
        n_below_min = 0
        n_no_td = 0
        n_with_td = 0

        # for a few examples only (avoid log spam)
        logged_nonzero_rows = 0
        LOG_EXAMPLES = 10

        for act_idx, flow_idx, value in zip(act_indices, flow_indices, values):
            n_triplets += 1
            act_idx = int(act_idx)
            flow_idx = int(flow_idx)

            supply_amt = float(supply_by_activity.get(act_idx, 0.0))
            if supply_amt == 0.0:
                continue

            n_supply_nonzero += 1

            scaled = supply_amt * float(value)
            if scaled == 0.0:
                continue

            if min_amount and abs(scaled) < min_amount:
                n_below_min += 1
                continue

            n_scaled_nonzero += 1

            # Optional: log a handful of actual contributing rows
            if logged_nonzero_rows < LOG_EXAMPLES:
                logger.info(
                    "accumulate_bio: contribute base_year=%d act=%d flow=%d supply=%g coeff=%g scaled=%g",
                    int(base_year), act_idx, flow_idx, supply_amt, float(value), float(scaled)
                )
                logged_nonzero_rows += 1

            tex = self.temporal_biosphere_exchanges.get((template_label, act_idx, flow_idx))

            if (tex is None) or (not use_temporal_distributions):
                n_no_td += 1
                y_eff = int(base_year)
                inventory_by_year.setdefault(
                    y_eff,
                    np.zeros(n_flows, dtype=self.value_dtype),
                )
                inventory_by_year[y_eff][flow_idx] += scaled
            else:
                n_with_td += 1
                td = TemporalDistribution(tex)

                any_pair = False
                for offset, weight in td.iter_offsets_and_weights():
                    any_pair = True
                    y_eff = int(base_year + offset)
                    factor = td.scale_factor(offset)

                    inventory_by_year.setdefault(
                        y_eff,
                        np.zeros(n_flows, dtype=self.value_dtype),
                    )
                    inventory_by_year[y_eff][flow_idx] += scaled * float(weight) * float(factor)

                if not any_pair:
                    logger.warning(
                        "accumulate_bio: TD produced no offsets/weights (template=%s act=%d flow=%d) -> dropped",
                        template_label, act_idx, flow_idx
                    )



    def _map_year_to_available(self, year: int) -> int:
        """
        Map an arbitrary year to the nearest available scenario year,
        clipped to [min_year, max_year].

        If you later interpolate annually, and every year in [min,max]
        is present, this is effectively just clipping.
        """
        # 1) Clip to global range
        y = max(self.min_year, min(self.max_year, int(year)))

        # 2) If we have a full annual grid, just return y
        if len(self.years_int) == (self.max_year - self.min_year + 1):
            return y

        # 3) Otherwise: snap to nearest available scenario year
        idx = int(np.abs(self.years_int - y).argmin())
        return int(self.years_int[idx])


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
    ):
        """
        Traverse the temporal-technosphere graph starting from (start_year, start_act_idx).

        Progress:
          - Prefer an empirical total estimate based on max_depth (DEPTH_TOTALS).
          - If not available, fall back to a short warm-up branching estimate.
        """

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
        queue.append((start_year, start_act_idx, float(amount), 0, ()))

        demand = defaultdict(float)
        provenance = defaultdict(lambda: defaultdict(float))
        bio_cache: dict[tuple[int, int], bool] = {}

        # ------------------------------------------------------------------
        # Main traversal loop
        # ------------------------------------------------------------------
        while queue:
            year, act, amt, depth, path = queue.popleft()

            if abs(amt) < min_amount:
                continue

            nodes_processed += 1
            if pbar is not None:
                pbar.update(1)

            # If we didn't have an empirical total, switch to an estimated total after warm-up
            if (
                    show_progress
                    and pbar is not None
                    and pbar.total is None
                    and nodes_processed >= WARMUP_LIMIT
            ):
                total_est = estimate_total_from_branching(branching_samples)
                total_est = max(total_est, nodes_processed + 1)
                pbar.total = total_est
                pbar.refresh()

            # Map to scenario year for looking into B
            scenario_year = self._map_year_to_scenario_year(year)
            label = str(scenario_year)

            if label in self.scenario_index:
                t = self.scenario_index[label]
                key = (scenario_year, act)
                if key in bio_cache:
                    has_direct_bio = bio_cache[key]
                else:
                    has_direct_bio = self.B[t, act, :].nnz > 0
                    bio_cache[key] = has_direct_bio

                    logger.info(
                        "temporal_traversal: node year=%d (scenario=%d t=%d) act=%d amt=%g depth=%d has_direct_bio=%s",
                        int(year), int(scenario_year), int(t), int(act), float(amt), int(depth), bool(has_direct_bio)
                    )
            else:
                has_direct_bio = False

            # Cutoff at max_depth
            if depth >= max_depth:
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt
                continue

            # Expand this node one step
            child_demands = self.expand_temporal_exchanges(
                year=year,
                act_idx=act,
                amount=amt,
                use_temporal_distributions=use_temporal_distributions,
            )

            if not child_demands:
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt
                continue

            if has_direct_bio:
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt

            # Enqueue children (+ warm-up branching sample if needed)
            children_enqueued = 0
            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue

                    child_node = (child_year, child_act)
                    if depth == 0:
                        child_path = (child_node,)
                    else:
                        child_path = path + (child_node,)

                    queue.append((child_year, child_act, child_amt, depth + 1, child_path))
                    children_enqueued += 1

            if show_progress and (pbar is not None) and (pbar.total is None) and nodes_processed <= WARMUP_LIMIT:
                branching_samples.append(children_enqueued)

        # ------------------------------------------------------------------
        # Force-complete progress bar (so it always ends at 100%)
        # ------------------------------------------------------------------
        if pbar is not None:
            try:
                # Case 1: total was never set (indeterminate bar) -> set it to actual work
                if pbar.total is None:
                    pbar.total = nodes_processed

                # Case 2: total exists but we processed fewer nodes than estimated -> fill to 100%
                if pbar.n < pbar.total:
                    pbar.update(pbar.total - pbar.n)

                pbar.refresh()
            except Exception:
                pass

            pbar.close()

        if pbar is not None:
            pbar.close()

        provenance = {}
        if return_provenance:
            provenance = {k: dict(v) for k, v in provenance.items()}
        return dict(demand), provenance

    def frontier_to_demand_vectors(self, frontier: dict) -> dict[int, np.ndarray]:
        """
        Convert a (year, activity) -> amount frontier into per-year demand vectors.

        IMPORTANT:
        - Keep calendar years as-is (do NOT map to scenario years here), otherwise
          all temporal links collapse to the nearest scenario year (e.g., 2050).
        """
        n_activities = self.A.shape[1]
        dtype = self.value_dtype

        f_by_year: dict[int, np.ndarray] = {}

        for (year, act_idx), amt in frontier.items():
            y = int(year)

            if y not in f_by_year:
                f_by_year[y] = np.zeros(n_activities, dtype=dtype)

            f_by_year[y][act_idx] += dtype(amt)

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
        queue.append((start_year, start_act_idx, float(amount), 0))

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