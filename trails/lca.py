from __future__ import annotations
import warnings
from typing import Any, Dict, List, TYPE_CHECKING

import bw2calc as bc
import numpy as np
from scikits.umfpack import UmfpackWarning
from tqdm import tqdm

from scikits.umfpack import UmfpackContext, UMFPACK_A
from scipy import sparse as sp

from .bw_interface import (
    _extract_supply_fast,
    _extract_supply_fast_cached,
    _get_datapackage,
    _reference_product_id_from_activity_id,
    build_datapackage_for_year_from_trails,
)
from .characterization import build_characterized_inventory
from .characterization import get_cf_vector

if TYPE_CHECKING:
    from .trails import Trails

warnings.filterwarnings("ignore", category=UmfpackWarning)
warnings.filterwarnings("ignore", module="scikits")

_CHAR_CACHE: dict = {}


def _get_mapping_arrays(mapping) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Return (ids, positions) arrays for a bw2calc dict mapping-like object."""
    if mapping:
        try:
            ids = np.fromiter(mapping.keys(), dtype=np.int64, count=len(mapping))
            pos = np.fromiter(mapping.values(), dtype=np.int64, count=len(mapping))
            return ids, pos
        except Exception:
            return None, None
    return None, None


def _map_activity_demands_to_products(
    lca_obj,
    activity_demands: dict[int, float],
) -> dict[int, float]:
    """Map activity-indexed demands to product ids using LCA metadata."""
    if not activity_demands:
        return {}

    act_map = getattr(lca_obj.dicts, "activity", None)
    prod_map = getattr(lca_obj.dicts, "product", None)
    if not act_map or not prod_map:
        raise ValueError("LCA object missing dicts.activity or dicts.product")

    mapped: dict[int, float] = {}
    prod_map_keys = prod_map.keys()
    cache: dict[int, int] = {}
    for act_id, amount in activity_demands.items():
        amt = float(amount)
        act_id = int(act_id)
        prod_id = cache.get(act_id)
        if prod_id is None:
            if act_id in prod_map_keys:
                prod_id = act_id
            else:
                prod_id = _reference_product_id_from_activity_id(lca_obj, act_id)
            cache[act_id] = int(prod_id)
        mapped[int(prod_id)] = mapped.get(int(prod_id), 0.0) + amt

    return mapped


def _build_rhs_matrix_from_root_demands(
    *,
    per_root_demands: dict[int, dict[int, float]],
    product_dict: dict,
    n: int,
    min_amount: float,
) -> tuple[list[int], np.ndarray]:
    """
    Build dense RHS matrix B (n, k) for all roots in `per_root_demands`.
    Columns correspond to roots in returned `roots` list.
    """
    roots = [int(r) for r in per_root_demands.keys()]
    k = len(roots)
    if k == 0:
        return [], np.zeros((n, 0), dtype=np.float64)

    B = np.zeros((n, k), dtype=np.float64)

    # Fill columns by mapping product IDs -> row indices (product_dict maps id -> row index)
    for j, root in enumerate(roots):
        demand = per_root_demands[root]
        for prod_id, v in demand.items():
            v = float(v)
            if v == 0.0:
                continue
            # Assume prod_id keys are valid; skip check_demand to avoid overhead
            try:
                i = product_dict[int(prod_id)]
            except KeyError:
                # Keep behavior strict if desired; otherwise `continue`
                raise
            B[int(i), j] += v

    return roots, B


def solve_many_rhs_umfpack_factorized(
    A_csc: sp.csc_matrix, B: np.ndarray
) -> np.ndarray:
    """
    Solve A X = B using a single UMFPACK factorization.

    Notes on scikits.umfpack behavior:
      - UmfpackContext.symbolic() and numeric() often return None on success.
      - ctx.solve(...) may return None and write into `x`, depending on version.
      - Therefore, rely on exceptions instead of status codes.
    """
    if not sp.isspmatrix_csc(A_csc):
        A_csc = A_csc.tocsc()

    # UMFPACK is float64-centric; keep this strict to avoid silent slow paths or failures.
    if A_csc.dtype != np.float64:
        A_csc = A_csc.astype(np.float64)

    # UMFPACK can be sensitive to unsorted indices; enforce canonical CSC.
    A_csc.sort_indices()

    B = np.asarray(B)
    if B.ndim != 2:
        raise ValueError("B must be 2D (n, k)")
    if B.dtype != np.float64:
        B = B.astype(np.float64, copy=False)

    n, k = B.shape
    if A_csc.shape != (n, n):
        raise ValueError(f"Shape mismatch: A {A_csc.shape}, B {B.shape}")

    ctx = UmfpackContext()

    # These calls frequently return None on success; failures raise exceptions.
    ctx.symbolic(A_csc)
    ctx.numeric(A_csc)

    X = np.empty((n, k), dtype=np.float64)

    # Reuse a preallocated x to reduce allocations
    for j in range(k):
        # In this scikits.umfpack version: solve(system, A, b, autoTranspose=True)
        X[:, j] = ctx.solve(UMFPACK_A, A_csc, B[:, j])

    return X


def lca_static_simple(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    debug: bool = False,
) -> None:
    """Run a static LCA for a single functional unit and year."""
    prev_inventory = trails.inventory
    prev_characterized = trails.characterized_inventory

    trails.reset_inventory(reset_scores=False)
    trails.static_score = None

    dp, _, _, _ = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        debug=debug,
    )

    # Build supply
    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    fu_prod_id = _reference_product_id_from_activity_id(lca_obj, int(fu_act_idx))

    lca_obj = bc.LCA(demand={int(fu_prod_id): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    supply_total = _extract_supply_fast(lca_obj, min_amount=0.0)

    # Accumulate inventory (no TD in static run)
    trails.accumulate_temporalized_biosphere_inventory(
        base_year=int(year),
        supply_by_activity=supply_total,
        min_amount=0.0,
        use_temporal_distributions=False,
        debug=debug,
    )

    trails.finalize_inventory()

    # Reuse the shared characterization cache to avoid rebuilding CF vectors.
    characterized = build_characterized_inventory(
        trails=trails,
        methods=methods,
        char_cache=_CHAR_CACHE,
    )
    trails.static_score = float(characterized.data.sum())

    trails.inventory = prev_inventory
    trails.characterized_inventory = prev_characterized


def lca(
    trails: Trails,
    methods: List[str],
    show_progress: bool = True,
    attribute_to_roots: bool = True,
    *,
    store_inventory: bool = False,
    compute_score: bool = True,
    ei_version: str = "3.11",
) -> None:
    """Run temporal LCA for a functional unit and year.

    When ``attribute_to_roots`` is enabled, biosphere impacts are accumulated under
    the first-level root activities while stored in the Trails inventory arrays with
    an added "root activity" dimension.
    """
    debug = bool(getattr(trails, "debug", False))

    trails.reset_inventory(attribute_to_roots=attribute_to_roots)

    cf = None
    if compute_score:
        cf = get_cf_vector(
            trails=trails,
            methods=methods,
            char_cache=_CHAR_CACHE,
            debug=debug,
            ei_version=ei_version,
        )

    # Ensure inventory builders are ready when we intend to store inventory data.
    if store_inventory:
        required = (
            hasattr(trails, "_inventory_years")
            and trails._inventory_years is not None
            and hasattr(trails, "_inv_chunk_flows")
            and hasattr(trails, "_inv_chunk_values")
            and hasattr(trails, "_inv_chunk_len")
        )
        if not required:
            raise RuntimeError(
                "BUG: store_inventory=True but reset_inventory() did not initialize "
                "chunk-based inventory builders."
            )

    # Routing must be run explicitly before LCA.
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing not initialized; run trails.temporal_routing(...) "
            "before trails.lca()."
        )

    routing_params = getattr(trails, "_routing_params", {}) or {}
    routing_attr_to_roots = getattr(trails, "_routing_attribute_to_roots", None)
    if routing_attr_to_roots is not None and routing_attr_to_roots != bool(
        attribute_to_roots
    ):
        raise RuntimeError(
            "Temporal routing was computed with attribute_to_roots="
            f"{routing_attr_to_roots}; rerun temporal_routing with "
            f"attribute_to_roots={attribute_to_roots}."
        )

    if not routing_params:
        raise RuntimeError(
            "Temporal routing parameters missing; rerun trails.temporal_routing(...)."
        )

    start_year_int = int(routing_params["start_year"])
    start_activity = int(routing_params["start_act_idx"])
    start_amount = float(routing_params["amount"])

    frontier: dict[tuple[int, int], float] = {}
    provenance: dict[tuple[int, int], dict[int, float]] = {}
    injected_supply_by_year_act: dict[tuple[int, int], float] = {}
    injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]] = {}

    for node, data in graph.nodes(data=True):
        year = int(data.get("year"))
        act = int(data.get("act_idx"))
        frontier_amt = float(data.get("frontier_amount") or 0.0)
        direct_bio_amt = float(data.get("direct_bio_amount") or 0.0)

        if frontier_amt:
            key = (year, act)
            frontier[key] = float(frontier.get(key, 0.0)) + frontier_amt
            if attribute_to_roots:
                roots = data.get("frontier_roots") or {}
                bucket = provenance.setdefault(key, {})
                for root_act, amt in roots.items():
                    bucket[int(root_act)] = float(bucket.get(int(root_act), 0.0)) + float(
                        amt
                    )

        if direct_bio_amt:
            key = (year, act)
            injected_supply_by_year_act[key] = float(
                injected_supply_by_year_act.get(key, 0.0)
            ) + direct_bio_amt
            if attribute_to_roots:
                roots = data.get("direct_bio_roots") or {}
                bucket = injected_supply_prov_by_year_act.setdefault(key, {})
                for root_act, amt in roots.items():
                    bucket[int(root_act)] = float(bucket.get(int(root_act), 0.0)) + float(
                        amt
                    )

    # Inject FU directly only when it is not already in the frontier (e.g., max_depth=0).
    if (start_year_int, start_activity) not in frontier:
        injected_supply_by_year_act[(start_year_int, start_activity)] = (
            float(injected_supply_by_year_act.get((start_year_int, start_activity), 0.0))
            + start_amount
        )
        injected_supply_prov_by_year_act.setdefault((start_year_int, start_activity), {})
        injected_supply_prov_by_year_act[(start_year_int, start_activity)][
            start_activity
        ] = (
            float(
                injected_supply_prov_by_year_act[(start_year_int, start_activity)].get(
                    start_activity, 0.0
                )
            )
            + start_amount
        )

    # Frontier -> demand vectors (calendar years preserved)
    frontier_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(frontier_by_year.keys())

    datapackage_cache: Dict[tuple, Any] = {}

    root_demands_by_year: dict[int, dict[int, dict[int, float]]] = {}
    root_injected_by_year: dict[int, dict[int, dict[int, float]]] = {}

    if attribute_to_roots:
        for (y, a), total_amt in frontier.items():
            root_map = provenance.get((y, a))
            if not root_map:
                root_map = {int(a): float(total_amt)}
            year_bucket = root_demands_by_year.setdefault(int(y), {})
            for root, amt in root_map.items():
                root_bucket = year_bucket.setdefault(int(root), {})
                root_bucket[int(a)] = root_bucket.get(int(a), 0.0) + float(amt)

        for (y, a), total_amt in injected_supply_by_year_act.items():
            root_map = injected_supply_prov_by_year_act.get((y, a))
            if not root_map:
                root_map = {int(a): float(total_amt)}
            year_bucket = root_injected_by_year.setdefault(int(y), {})
            for root, amt in root_map.items():
                root_bucket = year_bucket.setdefault(int(root), {})
                root_bucket[int(a)] = root_bucket.get(int(a), 0.0) + float(amt)

    pbar = None
    if show_progress:
        pbar = tqdm(
            total=len(candidate_years),
            desc="Temporal LCA: solve years",
            unit="year",
            leave=True,
        )

    for solve_year in candidate_years:
        root_ids = None
        root_supply_matrix = None

        solve_year = int(solve_year)
        demand_vector = np.asarray(frontier_by_year[solve_year])
        nonzero_indices = np.where(demand_vector != 0.0)[0]
        if nonzero_indices.size == 0:
            if pbar is not None:
                pbar.update(1)
            continue

        dp, _, _, _ = _get_datapackage(
            dp_cache=datapackage_cache,
            trails=trails,
            year=solve_year,
            zero_bio=True,
            debug=debug,
        )

        activity_demand = {int(i): float(demand_vector[i]) for i in nonzero_indices}

        lca_obj = bc.LCA(demand=activity_demand, data_objs=[dp])
        lca_obj.load_lci_data()  # build matrices + dicts

        # Cache mappings once per year
        act_map = getattr(lca_obj.dicts, "activity", None)
        prod_map = getattr(lca_obj.dicts, "product", None)

        act_ids, positions = _get_mapping_arrays(act_map)
        # For RHS building we need the product-id -> row-index dict.
        product_dict = prod_map  # mapping-like

        functional_unit_demand = _map_activity_demands_to_products(
            lca_obj,
            activity_demand,
        )
        if not functional_unit_demand:
            if pbar is not None:
                pbar.update(1)
            continue

        # Prepare the technosphere matrix once per year (UMFPACK expects CSC).
        A_csc = lca_obj.technosphere_matrix
        if not sp.isspmatrix_csc(A_csc):
            A_csc = A_csc.tocsc()

        supplies: list[tuple[Dict[int, float], int | None]] = []

        if attribute_to_roots:
            # Build one dense RHS matrix for all roots in this year
            per_root_demands_raw = root_demands_by_year.get(solve_year, {})
            per_root_demands: dict[int, dict[int, float]] = {}
            for root_act, demand in per_root_demands_raw.items():
                mapped = _map_activity_demands_to_products(
                    lca_obj,
                    demand,
                )
                if mapped:
                    per_root_demands[int(root_act)] = mapped
            roots, rhs_matrix = _build_rhs_matrix_from_root_demands(
                per_root_demands=per_root_demands,
                product_dict=product_dict,
                n=A_csc.shape[0],
                min_amount=float(min_amount),
            )

            if roots:
                if bc.PYPARDISO:
                    # Keep the PARDISO path for environments that rely on it.
                    for root_act in roots:
                        root_demand = per_root_demands[root_act]
                        lca_obj.build_demand_array(root_demand)
                        lca_obj.supply_array = lca_obj.solve_linear_system()
                        if act_ids is None or positions is None:
                            supply_total = _extract_supply_fast(lca_obj, min_amount)
                        else:
                            supply_total = _extract_supply_fast_cached(
                                lca_obj.supply_array, act_ids, positions, min_amount
                            )
                        if supply_total:
                            supplies.append((supply_total, int(root_act)))
                else:
                    # UMFPACK: factorize once and solve all RHS vectors.
                    root_supply_matrix = solve_many_rhs_umfpack_factorized(
                        A_csc, rhs_matrix
                    )
                    root_ids = np.asarray(
                        roots, dtype=np.int64
                    )  # Column order of RHS and solution.

                    for j, root_act in enumerate(roots):
                        supply_vec = root_supply_matrix[:, j]
                        if act_ids is None or positions is None:
                            # Fallback: emulate lca_obj.supply_array for extractor
                            lca_obj.supply_array = supply_vec
                            supply_total = _extract_supply_fast(lca_obj, min_amount)
                        else:
                            supply_total = _extract_supply_fast_cached(
                                supply_vec, act_ids, positions, min_amount
                            )
                        if supply_total:
                            supplies.append((supply_total, int(root_act)))

        else:
            # Original single-demand path (no per-root attribution)
            lca_obj.build_demand_array(functional_unit_demand)
            if not bc.PYPARDISO:
                lca_obj.decompose_technosphere()
            lca_obj.supply_array = lca_obj.solve_linear_system()

            if act_ids is None or positions is None:
                supply_total = _extract_supply_fast(lca_obj, min_amount)
            else:
                supply_total = _extract_supply_fast_cached(
                    lca_obj.supply_array, act_ids, positions, min_amount
                )
            if supply_total:
                supplies.append((supply_total, None))

        # Injected supply (direct additions outside the solved system).
        if attribute_to_roots:
            per_root_injected = root_injected_by_year.get(solve_year, {})
            for root_act, injected_supply in per_root_injected.items():
                if not injected_supply:
                    continue
                supplies.append((injected_supply, int(root_act)))
        else:
            injected_supply: Dict[int, float] = {}
            for (y, a), v in injected_supply_by_year_act.items():
                if int(y) != solve_year:
                    continue
                v = float(v)
                injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v

            if injected_supply:
                supplies.append((injected_supply, None))

        # ---- Inventory ----
        if supplies and store_inventory:
            for supply_dict, root_act in supplies:
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=solve_year,
                    supply_by_activity=supply_dict,
                    min_amount=float(min_amount),
                    store_activity=root_act,
                    debug=debug,
                )

        # ---- Scores ----
        if compute_score:
            # Build dense supply matrix (n_acts, n_roots_this_year).
            # Injected supplies remain dictionary-based because they are usually small.
            if attribute_to_roots:
                # Only call the matrix scorer when we solved a multi-RHS system
                # and the root ordering matches the solution columns.
                if root_ids is not None and root_supply_matrix is not None:
                    trails.accumulate_temporalized_biosphere_score_matrix(
                        base_year=solve_year,
                        supply_matrix=root_supply_matrix,
                        root_activities=root_ids,
                        cf=cf,
                        min_amount=float(min_amount),
                        use_temporal_distributions=True,
                        debug=debug,
                    )

                # Injected supplies: keep cheap dict scoring (usually small)
                per_root_injected = root_injected_by_year.get(solve_year, {})
                for root_act, injected_supply in per_root_injected.items():
                    if not injected_supply:
                        continue
                    trails.accumulate_temporalized_biosphere_score(
                        base_year=solve_year,
                        supply_by_activity=injected_supply,
                        cf=cf,
                        min_amount=float(min_amount),
                        store_activity=int(root_act),
                        debug=debug,
                    )
            else:
                # Non-root mode: single column supply array already exists
                # Use the dict path to keep the scorer simple.
                supply_total = supplies[0][0] if supplies else {}
                if supply_total:
                    trails.accumulate_temporalized_biosphere_score(
                        base_year=solve_year,
                        supply_by_activity=supply_total,
                        cf=cf,
                        min_amount=float(min_amount),
                        store_activity=None,
                        debug=debug,
                    )

        if pbar is not None:
            pbar.update(1)
    if pbar is not None:
        pbar.close()

    if store_inventory:
        trails.finalize_inventory()

    if compute_score:
        trails.finalize_scores()

        if trails.scores is None:
            raise RuntimeError(
                "compute_score=True but trails.scores is still None. "
                "This indicates lca() did not finalize scores correctly."
            )

    # Characterized inventory is optional.
    if store_inventory and (not compute_score):
        build_characterized_inventory(
            trails=trails, methods=methods, char_cache=_CHAR_CACHE
        )
