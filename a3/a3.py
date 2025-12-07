# a3.py

from typing import Dict, List, Optional
from collections import defaultdict, deque

import numpy as np
import sparse
import pyprind

from .datapackage import (
    load_matrices_from_package,
    interpolate_to_annual,
    load_indices_from_package,
    load_temporal_distributions,
)

from .temporal_distributions import TemporalDistribution


class A3:
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

        # Global classification-based temporal distributions (old design, maybe unused later)
        self.temporal_distributions = load_temporal_distributions(
            package, resource_name="temporal_distributions.csv"
        )

        # --- NEW: this will hold per-exchange temporal info
        self.temporal_exchanges: Dict = {}

        # Load core matrices + scenarios + per-exchange temporal metadata
        (
            self.A,
            self.B,
            self.scenario_labels,
            self.scenario_index,
            self.temporal_exchanges,   # <<< NEW
        ) = load_matrices_from_package(
            package=self.package,
            value_dtype=self.value_dtype,
            index_dtype=self.index_dtype,
        )

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

    def expand_temporal_exchanges(self, year: int, act_idx: int, amount: float = 1.0):
        """
        Expand one activity-year demand into temporally distributed multi-year
        demands for its *direct* exchanges.

        - A/B rows are always taken from the nearest available scenario year.
        - Temporal metadata (if any) is taken from the nearest original "template" year.
        """
        demand: dict[int, dict[int, float]] = {}

        # Which A/B slice do we use for this node?
        scenario_year = self._map_year_to_scenario_year(year)
        scenario_label = str(scenario_year)
        t = self.scenario_index[scenario_label]

        # And which year's temporal pattern should we use?
        template_year = self._map_year_to_template_year(year)
        template_label = str(template_year)

        # A[slice, act_idx, :]
        A_row = self.A[t, act_idx, :]

        if A_row.nnz == 0:
            return demand

        prod_indices = A_row.coords[0]
        values = A_row.data

        from .temporal_distributions import TemporalDistribution

        for prod_idx, value in zip(prod_indices, values):
            prod_idx = int(prod_idx)

            # Optional: skip self-production diagonal -1.0
            if prod_idx == act_idx and value == -1.0:
                continue

            # Temporal metadata is defined per (template_label, act_idx, prod_idx)
            tex = self.temporal_exchanges.get((template_label, act_idx, prod_idx))

            if tex is None:
                # No temporal shift → happens at the scenario_year
                y_eff = scenario_year
                demand.setdefault(y_eff, {})
                demand[y_eff][prod_idx] = demand[y_eff].get(prod_idx, 0.0) + amount * float(value)
            else:
                td = TemporalDistribution(tex)

                for offset, weight in td.iter_offsets_and_weights():
                    raw_year = year + offset  # conceptual calendar year
                    y_eff = self._map_year_to_scenario_year(raw_year)  # nearest scenario year

                    demand.setdefault(y_eff, {})
                    demand[y_eff][prod_idx] = (
                            demand[y_eff].get(prod_idx, 0.0)
                            + amount * float(value) * weight
                    )

        return demand

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
    ):
        """
        Traverse the temporal-technosphere graph starting from
        (start_year, start_act_idx) following direct exchanges up to
        `max_depth` layers.

        Parameters
        ----------
        start_year : int
        start_act_idx : int
        amount : float
        max_depth : int
        min_amount : float
        return_provenance : bool
            If True, also return a provenance mapping that tells,
            for each frontier node, which first-level exchange
            (direct child of the root activity) it came from.

        Returns
        -------
        If return_provenance is False:
            demand : dict[(year, act_idx), amount]

        If return_provenance is True:
            demand : dict[(year, act_idx), amount]
            provenance : dict[(year, act_idx),
                              {root_child_act_idx: amount_contribution}]
        """

        # queue entries: (year, act, amount, depth, root_child)
        #   root_child is:
        #     - None for the root node itself
        #     - act_idx of the first-level supplier for all descendants
        from collections import deque, defaultdict

        queue = deque()
        queue.append((start_year, start_act_idx, float(amount), 0, None))

        # Final demand at cutoff (year, act)
        demand = defaultdict(float)

        # provenance[(year, act)][root_child_act] = amount
        provenance = defaultdict(lambda: defaultdict(float))

        while queue:
            year, act, amt, depth, root_child = queue.popleft()

            if abs(amt) < min_amount:
                continue

            # If we've reached max_depth, we stop expanding and store demand here
            if depth >= max_depth:
                demand[(year, act)] += amt
                if root_child is not None:
                    provenance[(year, act)][root_child] += amt
                continue

            # Expand this node one step
            child_demands = self.expand_temporal_exchanges(
                year=year, act_idx=act, amount=amt
            )

            if not child_demands:
                # No outgoing exchanges → treat as a final demand node
                demand[(year, act)] += amt
                if root_child is not None:
                    provenance[(year, act)][root_child] += amt
                continue

            # Enqueue children
            for child_year, mapping in child_demands.items():
                for child_act, child_amt in mapping.items():
                    if abs(child_amt) < min_amount:
                        continue

                    # At depth 0 we are expanding the root activity;
                    # the first-level children define the root_child.
                    if depth == 0:
                        child_root = child_act
                    else:
                        child_root = root_child

                    queue.append(
                        (child_year, child_act, child_amt, depth + 1, child_root)
                    )

        # Return provenance if requested
        if return_provenance:
            provenance = {
                key: dict(inner) for key, inner in provenance.items()
            }
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
