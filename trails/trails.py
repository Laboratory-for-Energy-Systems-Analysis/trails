# trails.py

from typing import Dict, List, Optional
from collections import defaultdict, deque

import numpy as np
import sparse
import pyprind

from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package
)

from .temporal_distributions import TemporalDistribution


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

        self.temporal_exchanges: Dict = {}
        self.temporal_biosphere_exchanges: Dict = {}

        (
            self.A,
            self.B,
            self.scenario_labels,
            self.scenario_index,
            self.temporal_exchanges,
            self.temporal_biosphere_exchanges,
        ) = load_matrices_from_package(
            package=self.package,
            value_dtype=self.value_dtype,
            index_dtype=self.index_dtype,
        )


        # Backward-compatible support for tagged temporal exchange keys.
        # If a single dict contains ('tech', ...) and ('bio', ...) keys, split them.
        if self.temporal_biosphere_exchanges is None:
            self.temporal_biosphere_exchanges = {}
        if self.temporal_exchanges:
            sample_key = next(iter(self.temporal_exchanges.keys()))
            if isinstance(sample_key, tuple) and len(sample_key) >= 1 and sample_key[0] in ("tech", "bio"):
                tagged = self.temporal_exchanges
                self.temporal_exchanges = {}
                for k, v in tagged.items():
                    if not (isinstance(k, tuple) and len(k) == 4):
                        continue
                    kind, label, act, idx = k
                    if kind == "bio":
                        self.temporal_biosphere_exchanges[(label, act, idx)] = v
                    else:
                        self.temporal_exchanges[(label, act, idx)] = v



        self.template_labels = list(self.scenario_labels)
        self.template_years_int = np.array([int(lbl) for lbl in self.template_labels], dtype=int)

        self.years_int = np.array([int(lbl) for lbl in self.scenario_labels], dtype=int)
        self.min_year = int(self.years_int.min())
        self.max_year = int(self.years_int.max())

        # Load indices/metadata
        (
            self.activity_indices,
            self.biosphere_indices,
        ) = load_indices_from_package(self.package)

        # Optional temporal interpolation to annual resolution
        if interpolate_annual and self.scenario_labels:
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
        return self.temporal_exchanges.get((label, act_idx, prod_idx))

    def get_temporal_distribution(self, year: int, act_idx: int, prod_idx: int):
        """
        Return a TemporalDistribution object for (year, act_idx, prod_idx),
        or None if this exchange has no temporal metadata.
        """
        label = str(year)
        tex = self.temporal_exchanges.get((label, act_idx, prod_idx))
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

        template_year = self._map_year_to_template_year(year)
        template_label = str(template_year)

        A_row = self.A[t, act_idx, :]
        if A_row.nnz == 0:
            return demand

        prod_indices = A_row.coords[0]
        values = A_row.data

        from .temporal_distributions import TemporalDistribution

        for prod_idx, value in zip(prod_indices, values):
            prod_idx = int(prod_idx)

            if prod_idx == act_idx and value == -1.0:
                continue

            tex = self.temporal_exchanges.get((template_label, act_idx, prod_idx))

            # NEW: bypass temporal shifting
            if (tex is None) or (not use_temporal_distributions):
                y_eff = scenario_year
                demand.setdefault(y_eff, {})
                demand[y_eff][prod_idx] = demand[y_eff].get(prod_idx, 0.0) + amount * float(value)
                continue

            td = TemporalDistribution(tex)

            for offset, weight in td.iter_offsets_and_weights():
                raw_year = year + offset
                y_eff = self._map_year_to_scenario_year(raw_year)

                factor = td.scale_factor(offset)
                demand.setdefault(y_eff, {})
                demand[y_eff][prod_idx] = (
                        demand[y_eff].get(prod_idx, 0.0)
                        + amount * float(value) * float(weight) * float(factor)
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
        if self.B is None:
            return

        # Choose which B slice to read (scenario coefficients)
        scenario_year = self._map_year_to_scenario_year(base_year)
        scenario_label = str(scenario_year)
        if scenario_label not in self.scenario_index:
            return
        t = self.scenario_index[scenario_label]

        # Choose which temporal metadata year to use
        template_year = self._map_year_to_template_year(base_year)
        template_label = str(template_year)

        B_t = self.B[t, :, :]
        if getattr(B_t, "nnz", 0) == 0:
            return

        n_flows = int(self.B.shape[2])

        act_indices = B_t.coords[0].astype(int)
        flow_indices = B_t.coords[1].astype(int)
        values = B_t.data

        for act_idx, flow_idx, value in zip(act_indices, flow_indices, values):
            act_idx = int(act_idx)
            flow_idx = int(flow_idx)

            supply_amt = float(supply_by_activity.get(act_idx, 0.0))
            if supply_amt == 0.0:
                continue

            scaled = supply_amt * float(value)
            if min_amount and abs(scaled) < min_amount:
                continue

            tex = self.temporal_biosphere_exchanges.get((template_label, act_idx, flow_idx))

            if (tex is None) or (not use_temporal_distributions):
                # No temporal shift -> assign to the base calendar year
                y_eff = int(base_year)
                inventory_by_year.setdefault(
                    y_eff,
                    np.zeros(n_flows, dtype=self.value_dtype),
                )
                inventory_by_year[y_eff][flow_idx] += scaled

            else:
                td = TemporalDistribution(tex)

                for offset, weight in td.iter_offsets_and_weights():
                    y_eff = int(base_year + offset)
                    factor = td.scale_factor(offset)

                    # Ensure the vector exists for this year
                    inventory_by_year.setdefault(
                        y_eff,
                        np.zeros(n_flows, dtype=self.value_dtype),
                    )

                    inventory_by_year[y_eff][flow_idx] += scaled * float(weight) * float(factor)

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


    def _count_traversal_nodes(
        self,
        start_year: int,
        start_act_idx: int,
        amount: float,
        max_depth: int,
        min_amount: float,
    ) -> int:
        """
        Dry-run traversal just to count how many nodes would be processed.
        Used to size the progress bar.
        """
        queue = deque()
        queue.append((start_year, start_act_idx, float(amount), 0))

        count = 0

        while queue:
            year, act, amt, depth = queue.popleft()

            if abs(amt) < min_amount:
                continue

            count += 1

            if depth >= max_depth:
                continue

            child_demands = self.expand_temporal_exchanges(
                year=year, act_idx=act, amount=amt
            )

            if not child_demands:
                continue

            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue
                    queue.append((child_year, child_act, child_amt, depth + 1))

        return count

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
        Traverse the temporal-technosphere graph starting from
        (start_year, start_act_idx).

        Nodes with direct biosphere flows (non-zero B row) are always kept
        in the frontier, even if they have temporal technosphere children.

        Provenance tracks full *time-stamped* paths:

            provenance[(year, act)][path_tuple] = amount

        where path_tuple is a tuple of (year, act) pairs from the
        *first-level child* down to this node, e.g.:

            ((y1, act_child1), (y2, act_child2), ..., (yk, act_here))
        """

        from collections import deque, defaultdict

        # --------------------------------------------------------------
        # Optional: pre-count nodes to size the progress bar
        # --------------------------------------------------------------
        bar = None
        if show_progress:
            total_nodes = self._count_traversal_nodes(
                start_year=start_year,
                start_act_idx=start_act_idx,
                amount=amount,
                max_depth=max_depth,
                min_amount=min_amount,
            )
            if total_nodes <= 0:
                total_nodes = 1
            bar = pyprind.ProgBar(total_nodes, title='Temporal traversal')

        queue = deque()
        # path: tuple[(year, act), ...] from first-level child onward
        queue.append((start_year, start_act_idx, float(amount), 0, ()))

        # Final demand at cutoff (year, act)
        demand = defaultdict(float)

        # provenance[(year, act)][path_tuple] = amount
        provenance = defaultdict(lambda: defaultdict(float))

        # Cache: (scenario_year, act_idx) -> bool(has_direct_bio)
        bio_cache: dict[tuple[int, int], bool] = {}

        while queue:
            year, act, amt, depth, path = queue.popleft()

            if abs(amt) < min_amount:
                continue

            if bar is not None:
                bar.update()

            # Map to scenario year for looking into B
            scenario_year = self._map_year_to_scenario_year(year)
            label = str(scenario_year)

            if label in self.scenario_index:
                t = self.scenario_index[label]
                key = (scenario_year, act)
                if key in bio_cache:
                    has_direct_bio = bio_cache[key]
                else:
                    # B[t, act, :] is a 1D sparse row (flows)
                    B_row = self.B[t, act, :]
                    has_direct_bio = B_row.nnz > 0
                    bio_cache[key] = has_direct_bio
            else:
                has_direct_bio = False

            # If we've reached max_depth, stop expanding and store demand here
            if depth >= max_depth:
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt
                continue

            # Expand this node one step
            child_demands = self.expand_temporal_exchanges(
                year=year, act_idx=act, amount=amt, use_temporal_distributions=use_temporal_distributions,
            )

            if not child_demands:
                # No outgoing exchanges → final node
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt
                continue

            if has_direct_bio:
                demand[(year, act)] += amt
                if path:
                    provenance[(year, act)][path] += amt

            # Enqueue children
            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue

                    # Build child path:
                    # - at depth 0 (root → first level): start the path with the child node
                    # - otherwise: extend the existing path
                    child_node = (child_year, child_act)
                    if depth == 0:
                        child_path = (child_node,)
                    else:
                        child_path = path + (child_node,)

                    queue.append(
                        (child_year, child_act, child_amt, depth + 1, child_path)
                    )

        if return_provenance:
            provenance = {key: dict(inner) for key, inner in provenance.items()}
            return dict(demand), provenance

        return dict(demand)

    def frontier_to_demand_vectors(self, frontier: dict) -> dict[int, np.ndarray]:
        """
        Convert a (year, activity) -> amount frontier into per-year demand vectors.

        Parameters
        ----------
        frontier : dict[(int, int), float]
            Mapping from (year, act_idx) to demanded amount.

        Returns
        -------
        f_by_year : dict[int, np.ndarray]
            {year: demand_vector}, where demand_vector has length = n_activities
            and dtype = self.value_dtype.
        """
        n_activities = self.A.shape[1]
        dtype = self.value_dtype

        # Use defaultdict of arrays to accumulate
        f_by_year: dict[int, np.ndarray] = {}

        for (year, act_idx), amt in frontier.items():
            # Map to nearest scenario year just in case
            y_eff = self._map_year_to_scenario_year(year)
            label = str(y_eff)

            # Skip anything that cannot be mapped (shouldn't happen, but safe)
            if label not in self.scenario_index:
                continue

            if y_eff not in f_by_year:
                f_by_year[y_eff] = np.zeros(n_activities, dtype=dtype)

            f_by_year[y_eff][act_idx] += dtype(amt)

        return f_by_year

    def collect_traversal_nodes(
            self,
            start_year: int,
            start_act_idx: int,
            amount: float = 1.0,
            max_depth: int = 3,
            min_amount: float = 1e-12,
    ) -> dict[int, dict[tuple[int, int], float]]:
        """
        Traverse the temporal-technosphere graph starting from
        (start_year, start_act_idx), and record visited nodes by depth.

        Returns
        -------
        nodes_by_depth : dict[int, dict[(int, int), float]]
            {depth: {(year, act_idx): cumulative_amount, ...}, ...}

        Here, 'amount' is the total flow of the functional unit that
        reaches that node (year, act) at that depth.
        """
        queue = deque()
        queue.append((start_year, start_act_idx, float(amount), 0))

        # nodes_by_depth[depth][(year, act)] = amount
        nodes_by_depth: dict[int, dict[tuple[int, int], float]] = defaultdict(
            lambda: defaultdict(float)
        )

        while queue:
            year, act, amt, depth = queue.popleft()

            if abs(amt) < min_amount:
                continue

            # record this node at this depth
            nodes_by_depth[depth][(year, act)] += amt

            if depth >= max_depth:
                continue

            # Expand this node one step in time + technosphere
            child_demands = self.expand_temporal_exchanges(
                year=year, act_idx=act, amount=amt
            )

            if not child_demands:
                continue

            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue
                    queue.append((child_year, child_act, child_amt, depth + 1))

        return nodes_by_depth

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