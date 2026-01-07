import warnings
from collections import defaultdict
from typing import Any, Dict, List

import bw2calc as bc
import numpy as np
from scikits.umfpack import UmfpackWarning
from tqdm.auto import tqdm

from .bw_interface import (
    _extract_supply_fast,
    _extract_supply_fast_cached,
    _get_datapackage,
    _reference_product_id_from_activity_id,
    build_datapackage_for_year_from_trails,
)
from .characterization import (
    _characterize_impact_years,
    _top_flow_contributions,
    top_activity_contributions_from_cfvec,
)
from .lcia import fill_characterization_factors_matrices
from .trails import Trails

warnings.filterwarnings("ignore", category=UmfpackWarning)


def lca_static_simple(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    debug: bool = False,
) -> Dict[str, Any]:
    """Run a static LCA for a single functional unit and year.

    :param trails: Trails instance providing matrices and metadata.
    :type trails: Trails
    :param year: Calendar year to solve.
    :type year: int
    :param fu_act_idx: Functional unit activity index.
    :type fu_act_idx: int
    :param methods: LCIA method identifiers.
    :type methods: list[str]
    :param amount: Functional unit amount.
    :type amount: float
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: Result dictionaries for solve-year and impact-year views.
    :rtype: dict[str, Any]
    """
    # Build datapackage for that year with biosphere enabled
    dp, _, bio_idx, _ = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        debug=debug,
    )

    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    # NEW: rebuild LCA using product-id demand (correct space)
    fu_prod_id = _reference_product_id_from_activity_id(lca_obj, int(fu_act_idx))

    # Recreate (or redo) with correct demand; recreating is simplest & robust
    lca_obj = bc.LCA(demand={int(fu_prod_id): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    # Build characterization matrix in BW order (static)
    bw_bio_map = lca_obj.dicts.biosphere  # flow_id -> row position
    biosphere_matrix_dict = {
        int(flow_id): int(pos) for flow_id, pos in bw_bio_map.items()
    }

    biosphere_dict_simple = {
        (name, comp, subcomp): int(flow_id)
        for (name, comp, subcomp, unit), flow_id in bio_idx.items()
    }

    C = fill_characterization_factors_matrices(
        methods=methods,
        biosphere_matrix_dict=biosphere_matrix_dict,
        biosphere_dict=biosphere_dict_simple,
        debug=debug,
    )

    inv_vec = np.asarray(lca_obj.inventory.sum(axis=1)).reshape((-1, 1))
    score = float(np.sum(C.dot(inv_vec)))

    # If C is a sparse/dense matrix (n_char, n_flows), reduce to 1D CF per flow row
    # Common case: C is (1, n_flows)
    cf_vec = np.asarray(C.sum(axis=0)).ravel()

    score = float(cf_vec @ inv_vec.ravel())

    top_flows = _top_flow_contributions(
        lca_obj, biosphere_dict_simple, cf_vec, top_n=25
    )

    top_acts = top_activity_contributions_from_cfvec(
        lca_obj=lca_obj,
        cf_vec_bw_bio_rows=cf_vec,
        trails=trails,
        year_label=str(year),  # or label_for_matrix; see below
        top_n=25,
        min_abs=1e-12,
    )

    return {
        "results_by_solve_year": {
            int(year): {
                "fu": {int(fu_act_idx): float(amount)},
                "lca": lca_obj,
                "scores": score,
                "scores_per_root": {int(fu_act_idx): score},  # trivial "root"
            }
        },
        "results_by_impact_year": {
            int(year): {
                "scores": score,
                "scores_per_root": {int(fu_act_idx): score},
            }
        },
        "fu_activity": {int(fu_act_idx): float(amount)},
        "fu_product": {int(fu_prod_id): float(amount)},
        "top_flow_contributions": top_flows,
        "top_activity_contributions": top_acts,
    }


def lca(
    trails: Trails,
    start_year: int,
    start_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    max_depth: int = 2,
    min_amount: float = 1e-18,
    show_progress: bool = True,
    debug: bool = False,
    return_provenance: bool = False,
    use_temporal_distributions: bool = True,
) -> Dict[str, Any]:
    """Run temporal LCA for a functional unit and year.

    Returns:
      - results_by_solve_year: diagnostic per-year solves
      - results_by_impact_year: time series of impacts booked in impact years
    """

    def _normalize_root(fu0: int, root: int, legacy_root: int) -> int:
        """Normalize legacy root identifiers to the functional unit index.

        :param fu0: Functional unit activity index.
        :type fu0: int
        :param root: Root activity index to normalize.
        :type root: int
        :param legacy_root: Legacy sentinel root to replace.
        :type legacy_root: int
        :returns: Normalized root index.
        :rtype: int
        """
        root = int(root)
        return fu0 if root == legacy_root else root

    def _run_temporal_traversal(
        trails: Trails,
        y0: int,
        fu0: int,
        amt0: float,
        max_depth: int,
        min_amount: float,
        show_progress: bool,
        debug: bool,
    ) -> tuple[dict, dict, dict, dict]:
        """Run temporal traversal and normalize optional return structures.

        :param trails: Trails instance to traverse.
        :type trails: Trails
        :param y0: Start year.
        :type y0: int
        :param fu0: Functional unit activity index.
        :type fu0: int
        :param amt0: Functional unit amount.
        :type amt0: float
        :param max_depth: Maximum traversal depth.
        :type max_depth: int
        :param min_amount: Minimum amount threshold for traversal.
        :type min_amount: float
        :param show_progress: Whether to show traversal progress.
        :type show_progress: bool
        :param debug: Whether to emit debug logging.
        :type debug: bool
        :returns: Traversal frontier, provenance, and injected supply mappings.
        :rtype: tuple[dict, dict, dict, dict]
        """
        (
            frontier,
            provenance,
            injected_supply_by_year_act,
            injected_supply_prov_by_year_act,
        ) = trails.temporal_traversal(
            start_year=y0,
            start_act_idx=fu0,
            amount=amt0,
            max_depth=int(max_depth),
            min_amount=float(min_amount),
            return_provenance=True,
            show_progress=bool(show_progress),
            use_temporal_distributions=True,
            debug=debug,
        )

        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        if injected_supply_prov_by_year_act is None:
            injected_supply_prov_by_year_act = {}

        return (
            frontier,
            provenance,
            injected_supply_by_year_act,
            injected_supply_prov_by_year_act,
        )

    def _apply_fu_direct_injection(
        injected_supply_by_year_act: dict[tuple[int, int], float],
        injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]],
        y0: int,
        fu0: int,
        amt0: float,
        legacy_root: int,
    ) -> None:
        """Record functional-unit direct biosphere injection entries.

        :param injected_supply_by_year_act: Mapping of injected supply by year/activity.
        :type injected_supply_by_year_act: dict[tuple[int, int], float]
        :param injected_supply_prov_by_year_act: Provenance mapping by year/activity.
        :type injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]]
        :param y0: Start year.
        :type y0: int
        :param fu0: Functional unit activity index.
        :type fu0: int
        :param amt0: Functional unit amount.
        :type amt0: float
        :param legacy_root: Legacy sentinel root to normalize.
        :type legacy_root: int
        """
        injected_supply_by_year_act[(y0, fu0)] = (
            float(injected_supply_by_year_act.get((y0, fu0), 0.0)) + amt0
        )

        injected_supply_prov_by_year_act.setdefault((y0, fu0), {})
        injected_supply_prov_by_year_act[(y0, fu0)][fu0] = (
            float(injected_supply_prov_by_year_act[(y0, fu0)].get(fu0, 0.0)) + amt0
        )

        if legacy_root in injected_supply_prov_by_year_act[(y0, fu0)]:
            injected_supply_prov_by_year_act[(y0, fu0)][fu0] = float(
                injected_supply_prov_by_year_act[(y0, fu0)].get(fu0, 0.0)
            ) + float(injected_supply_prov_by_year_act[(y0, fu0)].pop(legacy_root))

    def _build_rooted_frontier(
        frontier: dict[tuple[int, int], float],
        provenance: dict[tuple[int, int], dict],
        fu0: int,
        root_tol: float,
    ) -> tuple[dict[tuple[int, int, int], float], set[int]]:
        """Build a frontier keyed by activity roots using provenance paths.

        :param frontier: Mapping of ``(year, activity)`` to amount.
        :type frontier: dict[tuple[int, int], float]
        :param provenance: Mapping of provenance path to amount shares.
        :type provenance: dict[tuple[int, int], dict]
        :param fu0: Functional unit activity index.
        :type fu0: int
        :param root_tol: Tolerance for residual attribution.
        :type root_tol: float
        :returns: Rooted frontier and set of roots seen.
        :rtype: tuple[dict[tuple[int, int, int], float], set[int]]
        """
        rooted_frontier = defaultdict(float)
        roots_seen = set()

        def _iter_path_nodes(path: tuple[tuple[int, int], ...] | None) -> list[int]:
            """Yield activity indices from a provenance path.

            :param path: Provenance path tuple to iterate.
            :type path: tuple[tuple[int, int], ...] | None
            :returns: Activity indices extracted from the path.
            :rtype: list[int]
            """
            if path is None:
                return
            if isinstance(path, (int, np.integer)):
                yield int(path)
                return
            try:
                for node in path:
                    if isinstance(node, (int, np.integer)):
                        yield int(node)
                    elif isinstance(node, (tuple, list)) and len(node) >= 2:
                        yield int(node[1])
            except TypeError:
                return

        def _root_from_path(
            path: tuple[tuple[int, int], ...] | int | None,
            fu_act: int,
            fallback_root: int,
        ) -> int:
            """Infer the root activity from a provenance path.

            :param path: Provenance path representation.
            :param fu_act: Functional unit activity index.
            :param fallback_root: Fallback root when none is found.
            :type path: object
            :type fu_act: int
            :type fallback_root: int
            :returns: Root activity index.
            :rtype: int
            """
            if isinstance(path, (int, np.integer)):
                r = int(path)
                return r if r != fu_act else int(fallback_root)

            for act in _iter_path_nodes(path):
                if act != fu_act:
                    return int(act)

            return int(fallback_root)

        for (year, act), total_amt in frontier.items():
            year = int(year)
            act = int(act)
            total_amt = float(total_amt)

            prov = provenance.get((year, act), {})

            if not prov:
                rooted_frontier[(year, act, act)] += total_amt
                roots_seen.add(act)
                continue

            prov_sum = 0.0
            for path, amt_share in prov.items():
                amt_share = float(amt_share)
                prov_sum += amt_share

                root_act = _root_from_path(path, fu_act=fu0, fallback_root=act)
                rooted_frontier[(year, act, int(root_act))] += amt_share
                roots_seen.add(int(root_act))

            residual = total_amt - prov_sum
            if abs(residual) > root_tol:
                rooted_frontier[(year, act, act)] += residual
                roots_seen.add(act)

        return rooted_frontier, roots_seen

    def _extend_roots_with_injected_supply(
        roots_seen: set[int],
        injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]],
        normalize_root: Any,
    ) -> None:
        """Extend the root set with injected supply provenance entries.

        :param roots_seen: Set of root activity indices.
        :type roots_seen: set[int]
        :param injected_supply_prov_by_year_act: Provenance mapping for injected supply.
        :type injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]]
        :param normalize_root: Function to normalize root identifiers.
        :type normalize_root: callable
        """
        for (_y, _a), roots_map in injected_supply_prov_by_year_act.items():
            for r in roots_map.keys():
                roots_seen.add(normalize_root(int(r)))

    def _build_root_demand_vectors(
        trails: Trails,
        rooted_frontier: dict[tuple[int, int, int], float],
        roots_seen: set[int],
        normalize_root: Any,
    ) -> dict[int, dict[int, np.ndarray]]:
        """Build per-root demand vectors from a rooted frontier.

        :param trails: Trails instance with matrix shapes.
        :type trails: Trails
        :param rooted_frontier: Rooted frontier mapping.
        :type rooted_frontier: dict[tuple[int, int, int], float]
        :param roots_seen: Iterable of roots to include.
        :type roots_seen: iterable[int]
        :param normalize_root: Function to normalize root identifiers.
        :type normalize_root: callable
        :returns: Mapping of root to per-year demand vectors.
        :rtype: dict[int, dict[int, numpy.ndarray]]
        """
        n_activities = int(trails.A.shape[1])
        dtype = trails.value_dtype
        f_by_year_by_root = {r: {} for r in roots_seen}

        for (year, act, root_act), amt in rooted_frontier.items():
            y = int(year)
            a = int(act)
            r = normalize_root(int(root_act))

            if y not in f_by_year_by_root[r]:
                f_by_year_by_root[r][y] = np.zeros(n_activities, dtype=dtype)
            f_by_year_by_root[r][y][a] += dtype(amt)

        return f_by_year_by_root

    def _build_demand_by_first_level_child(
        f_by_year_by_root: dict[int, dict[int, np.ndarray]],
        solve_year: int,
        min_amount: float,
    ) -> Dict[int, Dict[int, float]]:
        """Build demand mappings by first-level child for a solve year.

        :param f_by_year_by_root: Demand vectors by root and year.
        :type f_by_year_by_root: dict[int, dict[int, numpy.ndarray]]
        :param solve_year: Year to extract demand for.
        :type solve_year: int
        :param min_amount: Minimum amount threshold to include.
        :type min_amount: float
        :returns: Demand mapping by root and child index.
        :rtype: dict[int, dict[int, float]]
        """
        demand_by_first_level_child: Dict[int, Dict[int, float]] = {}
        for root_idx, by_year in f_by_year_by_root.items():
            vec = by_year.get(solve_year)
            if vec is None:
                continue
            arr_r = np.asarray(vec)
            nz_r = np.where(np.abs(arr_r) > float(min_amount))[0]
            if nz_r.size == 0:
                continue
            demand_by_first_level_child[int(root_idx)] = {
                int(i): float(arr_r[i]) for i in nz_r
            }
        return demand_by_first_level_child

    def _assert_rooted_closure(
        arr: np.ndarray,
        demand_by_first_level_child: dict[int, dict[int, float]],
        solve_year: int,
    ) -> None:
        """Assert that rooted demands sum to total demand for a solve year.

        :param arr: Total demand array for the solve year.
        :type arr: numpy.ndarray
        :param demand_by_first_level_child: Rooted demand mapping.
        :type demand_by_first_level_child: dict[int, dict[int, float]]
        :param solve_year: Solve year for error reporting.
        :type solve_year: int
        :raises ValueError: If rooted demand does not sum to total demand.
        """
        summed = np.zeros_like(arr, dtype=float)
        for dct in demand_by_first_level_child.values():
            for i, v in dct.items():
                summed[i] += v
        abs_err = float(np.max(np.abs(summed - arr))) if arr.size else 0.0
        scale = float(np.max(np.abs(arr))) if arr.size else 0.0
        tol = max(1e-9, 1e-6 * max(1.0, scale))
        if abs_err > tol:
            raise ValueError(
                f"Rooted demand does not sum to total demand in year {solve_year} (abs_err={abs_err:g}, tol={tol:g})."
            )

    def _build_injected_supply(
        injected_supply_by_year_act: dict[tuple[int, int], float],
        solve_year: int,
        min_amount: float,
    ) -> Dict[int, float]:
        """Build injected supply mapping for a solve year.

        :param injected_supply_by_year_act: Injected supply mapping by year/activity.
        :type injected_supply_by_year_act: dict[tuple[int, int], float]
        :param solve_year: Solve year to filter by.
        :type solve_year: int
        :param min_amount: Minimum magnitude to include.
        :type min_amount: float
        :returns: Mapping of activity index to injected supply amount.
        :rtype: dict[int, float]
        """
        injected_supply: Dict[int, float] = {}
        for (y, a), v in injected_supply_by_year_act.items():
            if int(y) != solve_year:
                continue
            v = float(v)
            if abs(v) <= float(min_amount):
                continue
            injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v
        return injected_supply

    def _build_injected_supply_by_root(
        injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]],
        solve_year: int,
        min_amount: float,
        normalize_root: Any,
    ) -> Dict[int, Dict[int, float]]:
        """Build injected supply mapping by root for a solve year.

        :param injected_supply_prov_by_year_act: Provenance mapping for injected supply.
        :type injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]]
        :param solve_year: Solve year to filter by.
        :type solve_year: int
        :param min_amount: Minimum magnitude to include.
        :type min_amount: float
        :param normalize_root: Function to normalize root identifiers.
        :type normalize_root: callable
        :returns: Mapping of root to injected supply entries.
        :rtype: dict[int, dict[int, float]]
        """
        injected_supply_by_first_level_child: Dict[int, Dict[int, float]] = {}
        for (y, act), roots_map in injected_supply_prov_by_year_act.items():
            if int(y) != solve_year:
                continue
            for r, share in roots_map.items():
                r = normalize_root(int(r))
                share = float(share)
                if abs(share) <= float(min_amount):
                    continue
                injected_supply_by_first_level_child.setdefault(r, {})
                injected_supply_by_first_level_child[r][int(act)] = (
                    injected_supply_by_first_level_child[r].get(int(act), 0.0) + share
                )
        return injected_supply_by_first_level_child

    # Legacy sentinel (may still appear in some provenance artifacts)
    LEGACY_FU_DIRECT_ROOT = -1
    ROOT_CLOSURE_TOL = 1e-12

    fu0 = int(start_act_idx)
    y0 = int(start_year)
    amt0 = float(amount)

    # All FU-direct biosphere will be booked under this root
    FU_DIRECT_ROOT = fu0

    def _normalize_root_local(r: int) -> int:
        return _normalize_root(fu0, int(r), LEGACY_FU_DIRECT_ROOT)

    # ---------------------------------------------------------
    # STATIC SHORTCUT
    # ---------------------------------------------------------
    if not use_temporal_distributions:
        return lca_static_simple(
            trails=trails,
            year=y0,
            fu_act_idx=fu0,
            methods=methods,
            amount=amt0,
            debug=debug,
        )

    # -----------------------------
    # 1) Temporal traversal (+ provenance for rooting)
    # -----------------------------
    need_roots = True

    if return_provenance or need_roots:
        (
            frontier,
            provenance,
            injected_supply_by_year_act,
            injected_supply_prov_by_year_act,
        ) = _run_temporal_traversal(
            trails=trails,
            y0=y0,
            fu0=fu0,
            amt0=amt0,
            max_depth=max_depth,
            min_amount=min_amount,
            show_progress=show_progress,
            debug=debug,
        )
    else:
        frontier, injected_supply_by_year_act = trails.temporal_traversal(
            start_year=y0,
            start_act_idx=fu0,
            amount=amt0,
            max_depth=int(max_depth),
            min_amount=float(min_amount),
            return_provenance=False,
            show_progress=bool(show_progress),
            use_temporal_distributions=True,
            debug=debug,
        )
        provenance = {}
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        injected_supply_prov_by_year_act = {}

    # -----------------------------
    # 1a) Inject FU self-supply for FU-direct biosphere only
    # -----------------------------
    _apply_fu_direct_injection(
        injected_supply_by_year_act=injected_supply_by_year_act,
        injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
        y0=y0,
        fu0=fu0,
        amt0=amt0,
        legacy_root=LEGACY_FU_DIRECT_ROOT,
    )

    # -----------------------------
    # 1b) Frontier -> per-year demand vectors
    # -----------------------------
    f_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(f_by_year.keys())

    # -----------------------------
    # 1c) Rooted frontier using provenance paths
    # -----------------------------
    rooted_frontier, roots_seen = _build_rooted_frontier(
        frontier=frontier,
        provenance=provenance,
        fu0=fu0,
        root_tol=ROOT_CLOSURE_TOL,
    )
    _extend_roots_with_injected_supply(
        roots_seen=roots_seen,
        injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
        normalize_root=_normalize_root_local,
    )
    roots_seen = sorted(roots_seen)

    # Convert rooted frontier to per-year demand vectors by root (technosphere only)
    f_by_year_by_root = _build_root_demand_vectors(
        trails=trails,
        rooted_frontier=rooted_frontier,
        roots_seen=roots_seen,
        normalize_root=_normalize_root_local,
    )

    # -----------------------------
    # Results + caches
    # -----------------------------
    results_by_solve_year: Dict[int, Dict[str, Any]] = {}
    dp_cache: Dict[tuple, Any] = {}
    char_cache: Dict[tuple, Any] = {}

    inventory_total_by_impact_year: Dict[int, np.ndarray] = {}
    inventory_by_root_by_impact_year: Dict[int, Dict[int, np.ndarray]] = defaultdict(
        dict
    )

    # -----------------------------
    # 2) Solve-year loop
    # -----------------------------
    solve_iter = candidate_years
    if show_progress:
        solve_iter = tqdm(
            candidate_years, desc="Temporal LCA: solve years", unit="year"
        )

    for solve_year in solve_iter:
        solve_year = int(solve_year)
        arr = np.asarray(f_by_year[solve_year])

        nz_idx = np.where(np.abs(arr) > float(min_amount))[0]
        if nz_idx.size == 0:
            continue

        fu_demand = {int(i): float(arr[i]) for i in nz_idx}

        # ------------------------------------------------------------------
        # Build rooted demand mapping for this solve year (MUST exist before per-root loop)
        # ------------------------------------------------------------------
        demand_by_first_level_child = _build_demand_by_first_level_child(
            f_by_year_by_root=f_by_year_by_root,
            solve_year=solve_year,
            min_amount=min_amount,
        )
        _assert_rooted_closure(arr, demand_by_first_level_child, solve_year)

        # ------------------------------------------------------------------
        # Build datapackage (zero biosphere) + one LCA object per solve_year
        # ------------------------------------------------------------------
        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=solve_year,
            zero_bio=True,
            debug=debug,
        )

        lca_obj = bc.LCA(demand=fu_demand, data_objs=[dp])
        # factorize=True helps when reusing via redo_lci
        lca_obj.lci(factorize=True)

        # supply extraction (total)
        # ---------------------------------------------------------
        # Cache activity->position mapping once per solve_year
        # ---------------------------------------------------------
        act_map = getattr(lca_obj.dicts, "activity", None)
        if not act_map:
            act_ids = None
            positions = None
        else:
            act_ids = np.fromiter(act_map.keys(), dtype=np.int64, count=len(act_map))
            positions = np.fromiter(
                act_map.values(), dtype=np.int64, count=len(act_map)
            )

        # supply extraction (total) using cached mapping when available
        if act_ids is None or positions is None:
            supply_total = _extract_supply_fast(lca_obj, min_amount)
        else:
            supply_total = _extract_supply_fast_cached(
                lca_obj.supply_array, act_ids, positions, min_amount
            )

        # (1) solved-supply contribution
        trails.accumulate_temporalized_biosphere_inventory(
            base_year=solve_year,
            supply_by_activity=supply_total,
            inventory_by_year=inventory_total_by_impact_year,
            min_amount=float(min_amount),
            use_temporal_distributions=True,
            debug=debug,
        )

        # (2) injected supply for this solve_year (no solve)
        injected_supply = _build_injected_supply(
            injected_supply_by_year_act=injected_supply_by_year_act,
            solve_year=solve_year,
            min_amount=min_amount,
        )
        if injected_supply:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=injected_supply,
                inventory_by_year=inventory_total_by_impact_year,
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

        # Build injected supply by root, normalizing legacy sentinel -> fu0
        injected_supply_by_first_level_child = _build_injected_supply_by_root(
            injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
            solve_year=solve_year,
            min_amount=min_amount,
            normalize_root=_normalize_root_local,
        )

        # Book FU-direct inventory under FU activity root (fu0)
        fu_direct_injected = injected_supply_by_first_level_child.get(
            FU_DIRECT_ROOT, {}
        )
        if fu_direct_injected:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=fu_direct_injected,
                inventory_by_year=inventory_by_root_by_impact_year[FU_DIRECT_ROOT],
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

        # ------------------------------------------------------------------
        # Per-root inventories: reuse LCA object (redo_lci), fall back if missing
        # ------------------------------------------------------------------
        has_redo = hasattr(lca_obj, "redo_lci")

        for root_idx, root_demand in demand_by_first_level_child.items():
            if has_redo:
                lca_obj.redo_lci(demand=root_demand)
            else:
                # conservative fallback
                lca_obj.demand = root_demand
                lca_obj.lci()

            if act_ids is None or positions is None:
                supply_root = _extract_supply_fast(lca_obj, min_amount)
            else:
                supply_root = _extract_supply_fast_cached(
                    lca_obj.supply_array, act_ids, positions, min_amount
                )

            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=supply_root,
                inventory_by_year=inventory_by_root_by_impact_year[int(root_idx)],
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

            root_injected = injected_supply_by_first_level_child.get(int(root_idx), {})
            if root_injected:
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=solve_year,
                    supply_by_activity=root_injected,
                    inventory_by_year=inventory_by_root_by_impact_year[int(root_idx)],
                    min_amount=float(min_amount),
                    use_temporal_distributions=True,
                    debug=debug,
                )

        # Diagnostics roots present either via rooted demand or injected supply
        all_diag_roots = set(int(k) for k in demand_by_first_level_child.keys()) | set(
            int(k) for k in injected_supply_by_first_level_child.keys()
        )

        results_by_solve_year[solve_year] = {
            "fu_demand": fu_demand,
            "n_nonzero_demand": int(len(fu_demand)),
            "sum_abs_demand": float(np.sum(np.abs(arr[nz_idx]))),
            "max_abs_demand": float(np.max(np.abs(arr[nz_idx]))),
            "n_injected_supply": int(len(injected_supply)),
            "injected_supply": injected_supply,
            "roots": sorted(int(k) for k in all_diag_roots),
            "demand_by_first_level_child": demand_by_first_level_child,
            "injected_supply_by_first_level_child": injected_supply_by_first_level_child,
            "n_supply_total": int(len(supply_total)),
            "supply_total": supply_total if debug else None,
            "lca": lca_obj,  # last solved state is last root; keep for introspection only
        }

    # -----------------------------
    # 3) Impact-year characterization
    # -----------------------------
    results_by_impact_year = _characterize_impact_years(
        trails=trails,
        inventory_total_by_impact_year=inventory_total_by_impact_year,
        inventory_by_root_by_impact_year=inventory_by_root_by_impact_year,
        dp_cache=dp_cache,
        char_cache=char_cache,
        methods=methods,
        min_amount=min_amount,
        normalize_root=_normalize_root_local,
        debug=debug,
    )

    out = {
        "results_by_solve_year": results_by_solve_year,
        "results_by_impact_year": results_by_impact_year,
    }
    return (out, provenance) if return_provenance else out
