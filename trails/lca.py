from __future__ import annotations
import warnings
from typing import Any, Dict, List, TYPE_CHECKING

import bw2calc as bc
import numpy as np
from scikits.umfpack import UmfpackWarning
from tqdm import tqdm

from scikits.umfpack import UmfpackContext, UMFPACK_A
from scipy import sparse as sp
import sparse
import xarray as xr

from .bw_interface import (
    _extract_supply_fast,
    _extract_supply_fast_cached,
    _get_datapackage,
    _reference_product_id_from_activity_id,
    _reference_product_from_activity_id,
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
    cache: dict[int, tuple[int, float]] = {}
    for act_id, amount in activity_demands.items():
        amt = float(amount)
        act_id = int(act_id)
        cached = cache.get(act_id)
        if cached is None:
            prod_id, prod_value = _reference_product_from_activity_id(lca_obj, act_id)
            cache[act_id] = (int(prod_id), float(prod_value))
        else:
            prod_id, prod_value = cached
        sign = -1.0 if prod_value < 0.0 else 1.0
        mapped[int(prod_id)] = mapped.get(int(prod_id), 0.0) + amt * sign

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
    A_csc: sp.csc_matrix,
    B: np.ndarray,
    *,
    cache: dict | None = None,
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

    ctx = None
    symbolic = None
    pattern_sig = None
    if cache is not None:
        ctx = cache.get("ctx")
        symbolic = cache.get("symbolic")
        pattern_sig = cache.get("pattern_sig")
    if ctx is None:
        ctx = UmfpackContext()
        if cache is not None:
            cache["ctx"] = ctx

    # Reuse symbolic factorization when sparsity pattern is unchanged.
    sig = (
        A_csc.indptr.size,
        A_csc.indices.size,
        int(np.bitwise_xor.reduce(A_csc.indptr)),
        int(np.bitwise_xor.reduce(A_csc.indices)),
    )
    if cache is None or sig != pattern_sig or symbolic is None:
        symbolic = ctx.symbolic(A_csc)
        if cache is not None:
            cache["symbolic"] = symbolic
            cache["pattern_sig"] = sig

    # These calls frequently return None on success; failures raise exceptions.
    try:
        ctx.numeric(A_csc, symbolic)
    except TypeError:
        ctx.numeric(A_csc)

    X = np.empty((n, k), dtype=np.float64)

    # Reuse a preallocated x to reduce allocations
    for j in range(k):
        # In this scikits.umfpack version: solve(system, A, b, autoTranspose=True)
        X[:, j] = ctx.solve(UMFPACK_A, A_csc, B[:, j])

    return X


def lca_static(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    debug: bool = False,
    ei_version: str = "3.11",
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

    # Build supply using the same activity->product mapping as temporal LCA
    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    inv = lca_obj.inventory  # SciPy sparse (flow_pos x act_pos)

    inv_coo = sparse.COO.from_scipy_sparse(inv)
    flow_pos = inv_coo.coords[0].astype(np.int64, copy=False)
    act_pos = inv_coo.coords[1].astype(np.int64, copy=False)

    act_map = getattr(lca_obj.dicts, "activity", None) or {}
    bio_map = getattr(lca_obj.dicts, "biosphere", None) or {}
    pos_to_act = {int(pos): int(act_id) for act_id, pos in act_map.items()}
    pos_to_flow = {int(pos): int(flow_id) for flow_id, pos in bio_map.items()}

    act_ids = np.array([pos_to_act.get(int(p), -1) for p in act_pos], dtype=np.int64)
    flow_ids = np.array([pos_to_flow.get(int(p), -1) for p in flow_pos], dtype=np.int64)

    valid = (act_ids >= 0) & (flow_ids >= 0)
    act_ids = act_ids[valid]
    flow_ids = flow_ids[valid]
    data = inv_coo.data[valid]

    n_acts = int(trails.A.shape[1]) if trails.A is not None else int(act_ids.max() + 1)
    n_flows = (
        int(trails.B.shape[2]) if trails.B is not None else int(flow_ids.max() + 1)
    )

    inv_coo = sparse.COO(
        coords=np.vstack([act_ids, flow_ids]),
        data=data,
        shape=(n_acts, n_flows),
    )

    trails.inventory = xr.DataArray(
        inv_coo[:, :, None],  # add year axis
        dims=("activity", "flow", "year"),
        coords={
            "activity": np.arange(n_acts, dtype=int),
            "flow": np.arange(n_flows, dtype=int),
            "year": np.array([int(year)], dtype=int),
        },
    )
    trails.demand = lca_obj.demand

    characterized_inventory = build_characterized_inventory(
        trails=trails, methods=methods, char_cache=_CHAR_CACHE
    )

    if "method" in characterized_inventory.dims:
        inv = trails.inventory
        if inv is None:
            raise RuntimeError("Static inventory missing while scoring methods.")
        inv_data = inv.data
        if not isinstance(inv_data, sparse.COO):
            inv_data = sparse.COO.from_numpy(np.asarray(inv_data))
        flow_coords = inv_data.coords[1]
        vals = inv_data.data
        scores: dict[str, float] = {}
        for m in methods:
            cf = get_cf_vector(
                trails=trails,
                methods=[m],
                char_cache=_CHAR_CACHE,
                debug=debug,
                ei_version=ei_version,
            )
            score = float(np.dot(vals, cf[flow_coords]))
            scores[str(m)] = score
        trails.static_score = scores
    else:
        trails.static_score = float(characterized_inventory.data.sum())

    trails.inventory = prev_inventory
    trails.characterized_inventory = prev_characterized


def lca(
    trails: Trails,
    methods: List[str] | None = None,
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

    umfpack_cache: dict | None = (
        {} if (attribute_to_roots and not bc.PYPARDISO) else None
    )

    cf = None
    if compute_score:
        if not methods:
            raise ValueError("methods must be provided when compute_score=True.")
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
    min_amount = float(routing_params.get("min_amount", 1e-18))

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
                    bucket[int(root_act)] = float(
                        bucket.get(int(root_act), 0.0)
                    ) + float(amt)
        # Inject direct biosphere only for expanded nodes (not frontier nodes).
        if direct_bio_amt and not frontier_amt:
            key = (year, act)
            injected_supply_by_year_act[key] = (
                float(injected_supply_by_year_act.get(key, 0.0)) + direct_bio_amt
            )
            if attribute_to_roots:
                roots = data.get("direct_bio_roots") or {}
                bucket = injected_supply_prov_by_year_act.setdefault(key, {})
                for root_act, amt in roots.items():
                    bucket[int(root_act)] = float(
                        bucket.get(int(root_act), 0.0)
                    ) + float(amt)

    # Inject FU directly only when it is not already in the frontier (e.g., max_depth=0).
    if (start_year_int, start_activity) not in frontier:
        injected_supply_by_year_act[(start_year_int, start_activity)] = (
            float(
                injected_supply_by_year_act.get((start_year_int, start_activity), 0.0)
            )
            + start_amount
        )
        injected_supply_prov_by_year_act.setdefault(
            (start_year_int, start_activity), {}
        )
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
                        A_csc, rhs_matrix, cache=umfpack_cache
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

    # Characterized inventory is optional and only meaningful when scoring.
    if store_inventory and compute_score:
        build_characterized_inventory(
            trails=trails, methods=methods, char_cache=_CHAR_CACHE
        )


def build_temporal_sankey_tree(
    trails: Trails,
    *,
    root_year: int | None = None,
    root_act_idx: int | None = None,
    max_depth: int | None = None,
    min_amount: float = 0.0,
    sort_children: bool = True,
) -> dict[str, Any]:
    """Build a nested dict from the temporal routing graph for Sankey-style plots.

    The returned structure is a recursive tree:
        {
          "node": {...},
          "children": [
              {"edge_amount": float, "node": {...}, "children": [...]},
              ...
          ],
        }

    Nodes are taken from ``trails.graph`` built by ``trails.temporal_routing()``.

    :param trails: Trails instance with a populated temporal routing graph.
    :type trails: Trails
    :param root_year: Optional root year override (defaults to routing params).
    :type root_year: int | None
    :param root_act_idx: Optional root activity override (defaults to routing params).
    :type root_act_idx: int | None
    :param max_depth: Optional depth cutoff (inclusive of root at depth 0).
    :type max_depth: int | None
    :param min_amount: Filter edges with abs(amount) below this threshold.
    :type min_amount: float
    :param sort_children: Sort children by abs(edge_amount) descending.
    :type sort_children: bool
    :returns: Nested Sankey-ready tree dict.
    :rtype: dict
    """
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing graph missing; run trails.temporal_routing(...) first."
        )

    routing_params = getattr(trails, "_routing_params", {}) or {}
    if root_year is None:
        if "start_year" not in routing_params:
            raise RuntimeError(
                "Root year not provided and routing params missing; "
                "rerun temporal_routing or pass root_year."
            )
        root_year = int(routing_params["start_year"])
    if root_act_idx is None:
        if "start_act_idx" not in routing_params:
            raise RuntimeError(
                "Root activity not provided and routing params missing; "
                "rerun temporal_routing or pass root_act_idx."
            )
        root_act_idx = int(routing_params["start_act_idx"])

    min_amount = float(min_amount)

    root_key = None
    for node, data in graph.nodes(data=True):
        if int(data.get("depth", -1)) != 0:
            continue
        if int(data.get("year", -9999)) != int(root_year):
            continue
        if int(data.get("act_idx", -1)) != int(root_act_idx):
            continue
        root_key = node
        break

    if root_key is None:
        raise ValueError(
            "Root node not found in temporal routing graph for "
            f"(year={int(root_year)}, act_idx={int(root_act_idx)})."
        )

    def _node_payload(node_key: tuple) -> dict[str, Any]:
        data = graph.nodes[node_key]
        return {
            "key": node_key,
            "year": int(data.get("year")),
            "depth": int(data.get("depth")),
            "act_idx": int(data.get("act_idx")),
            "name": data.get("name") or "",
            "reference_product": data.get("reference_product") or "",
            "location": data.get("location") or "",
            "amount": float(data.get("amount") or 0.0),
            "frontier_amount": float(data.get("frontier_amount") or 0.0),
            "direct_bio_amount": float(data.get("direct_bio_amount") or 0.0),
        }

    def _build_tree(node_key: tuple) -> dict[str, Any]:
        node_data = graph.nodes[node_key]
        depth = int(node_data.get("depth"))
        if max_depth is not None and depth >= int(max_depth):
            return {"node": _node_payload(node_key), "children": []}

        children = []
        for child in graph.successors(node_key):
            edge_data = graph.edges[node_key, child]
            edge_amt = float(edge_data.get("amount") or 0.0)
            if abs(edge_amt) < min_amount:
                continue
            child_tree = _build_tree(child)
            children.append(
                {
                    "edge_amount": edge_amt,
                    "node": child_tree["node"],
                    "children": child_tree["children"],
                }
            )

        if sort_children:
            children.sort(key=lambda c: abs(c["edge_amount"]), reverse=True)

        return {"node": _node_payload(node_key), "children": children}

    return _build_tree(root_key)


def score_temporal_graph_nodes(
    trails: Trails,
    methods: List[str],
    *,
    min_amount: float = 0.0,
    show_progress: bool = True,
    ei_version: str = "3.11",
) -> dict[tuple, float]:
    """Score nodes in the temporal routing graph for Sankey-style weighting.

    Rules:
      - Non-frontier nodes: score direct biosphere only.
      - Frontier nodes: run a full LCA solve for that node-year demand.

    This function does not require calling ``trails.lca()`` beforehand.
    It only requires a populated routing graph and valid matrices in ``trails``.

    Node scores are keyed by the graph node key (tuple used by networkx).
    """
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing graph missing; run trails.temporal_routing(...) first."
        )

    if trails.A is None or trails.B is None:
        raise RuntimeError("Trails matrices missing; run setup before scoring.")

    cf = get_cf_vector(
        trails=trails,
        methods=methods,
        char_cache=_CHAR_CACHE,
        debug=bool(getattr(trails, "debug", False)),
        ei_version=ei_version,
    )

    min_amount = float(min_amount)

    # Caches
    dp_cache: Dict[tuple, Any] = {}
    char_row_cache: dict[int, np.ndarray] = {}  # t -> per-activity coeffs
    has_bio_cache: dict[tuple[int, int], bool] = {}
    supply_cache: dict[tuple[int, int], dict[int, float]] = {}

    def _char_row_for_year(year: int) -> np.ndarray:
        context = trails._get_scenario_context(int(year))
        if context is None:
            return np.zeros(int(trails.A.shape[1]), dtype=np.float64)
        _scenario_year, _label, t = context
        if t in char_row_cache:
            return char_row_cache[t]
        B_t = trails.B[t, :, :]
        # B_t is (activity, flow); multiply by CF (flow) -> per-activity coefficients
        coeff = B_t @ cf  # type: ignore[operator]
        coeff = np.asarray(coeff, dtype=np.float64)
        char_row_cache[t] = coeff
        return coeff

    def _score_direct_td(year: int, act: int, amount: float) -> float:
        if amount == 0.0:
            return 0.0

        # Intercept score appends to avoid mutating trails.scores
        total = 0.0
        original_append = trails._append_score_entry
        original_append_bulk = getattr(trails, "_append_scores_bulk", None)

        def _capture_append(
            act_idx: int,
            year_val: int,
            value: float,
            *,
            root_activity: int | None = None,
        ) -> None:
            nonlocal total
            total += float(value)

        def _capture_append_bulk(
            act_idx: int,
            year_idx: int,
            value: float,
            *,
            root_activity: int | None = None,
        ) -> None:
            nonlocal total
            total += float(value)

        trails._append_score_entry = _capture_append  # type: ignore[assignment]
        if original_append_bulk is not None:
            trails._append_scores_bulk = _capture_append_bulk  # type: ignore[assignment]

        try:
            trails.accumulate_temporalized_biosphere_score(
                base_year=int(year),
                supply_by_activity={int(act): float(amount)},
                cf=cf,
                min_amount=min_amount,
                store_activity=None,
                use_temporal_distributions=True,
                debug=bool(getattr(trails, "debug", False)),
            )
        finally:
            trails._append_score_entry = original_append  # type: ignore[assignment]
            if original_append_bulk is not None:
                trails._append_scores_bulk = original_append_bulk  # type: ignore[assignment]

        return float(total)

    def _solve_supply(year: int, act: int, amount: float) -> dict[int, float]:
        key = (int(year), int(act))
        if key in supply_cache:
            if amount == 1.0:
                return supply_cache[key]
            return {a: float(v) * float(amount) for a, v in supply_cache[key].items()}

        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=int(year),
            zero_bio=False,
            debug=bool(getattr(trails, "debug", False)),
        )

        lca_obj = bc.LCA(demand={int(act): 1.0}, data_objs=[dp])
        lca_obj.load_lci_data()

        activity_demand = {int(act): 1.0}
        functional_unit_demand = _map_activity_demands_to_products(
            lca_obj, activity_demand
        )
        if not functional_unit_demand:
            supply_cache[key] = {}
            return {}

        act_map = getattr(lca_obj.dicts, "activity", None)
        act_ids, positions = _get_mapping_arrays(act_map)

        lca_obj.build_demand_array(functional_unit_demand)
        if not bc.PYPARDISO:
            lca_obj.decompose_technosphere()
        lca_obj.supply_array = lca_obj.solve_linear_system()

        if act_ids is None or positions is None:
            supply_total = _extract_supply_fast(lca_obj, 0.0)
        else:
            supply_total = _extract_supply_fast_cached(
                lca_obj.supply_array, act_ids, positions, 0.0
            )

        supply_cache[key] = supply_total
        if amount == 1.0:
            return supply_total
        return {a: float(v) * float(amount) for a, v in supply_total.items()}

    def _score_frontier(year: int, act: int, amount: float) -> float:
        if amount == 0.0:
            return 0.0
        supply = _solve_supply(year, act, amount)
        if not supply:
            return 0.0
        coeff = _char_row_for_year(year)
        total = 0.0
        for a, v in supply.items():
            a = int(a)
            if a < 0 or a >= coeff.size:
                continue
            total += float(v) * float(coeff[a])
        return float(total)

    def _is_frontier(node_data: dict) -> bool:
        frontier_amt = float(node_data.get("frontier_amount") or 0.0)
        if frontier_amt != 0.0:
            return True
        return False

    def _direct_amount(node_data: dict) -> float:
        direct_amt = float(node_data.get("direct_bio_amount") or 0.0)
        if direct_amt != 0.0:
            return direct_amt
        depth = int(node_data.get("depth", 0))
        if depth == 0:
            year = int(node_data.get("year"))
            act = int(node_data.get("act_idx"))
            if trails._has_direct_biosphere(
                int(trails._map_year_to_scenario_year(year)), act, has_bio_cache
            ):
                return float(node_data.get("amount") or 0.0)
        return 0.0

    nodes = list(graph.nodes(data=True))
    pbar = None
    if show_progress:
        pbar = tqdm(
            total=len(nodes),
            desc="Score temporal nodes",
            unit="node",
            leave=True,
        )

    node_scores: dict[tuple, float] = {}

    for node_key, data in nodes:
        year = int(data.get("year"))
        act = int(data.get("act_idx"))
        if _is_frontier(data):
            amount = float(data.get("frontier_amount") or 0.0)
            if abs(amount) >= min_amount:
                node_scores[node_key] = _score_frontier(year, act, amount)
            else:
                node_scores[node_key] = 0.0
        else:
            amount = _direct_amount(data)
            if abs(amount) >= min_amount:
                node_scores[node_key] = _score_direct_td(year, act, amount)
            else:
                node_scores[node_key] = 0.0

        if pbar is not None:
            pbar.update(1)

    if pbar is not None:
        pbar.close()

    return node_scores
