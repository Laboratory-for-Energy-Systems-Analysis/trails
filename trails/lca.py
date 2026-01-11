from __future__ import annotations
import sys

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

    # IMPORTANT: reuse cache instead of passing {}
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
    start_year: int,
    start_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    max_depth: int = 2,
    min_amount: float = 1e-18,
    show_progress: bool = True,
    attribute_to_roots: bool = False,
    debug: bool = False,
    *,
    store_inventory: bool = True,
    compute_score: bool = True,
    ei_version: str = "3.11",
) -> None:
    """Run temporal LCA for a functional unit and year.

    When ``attribute_to_roots`` is enabled, biosphere impacts are accumulated under
    the first-level root activities while stored in the Trails inventory arrays with
    an added "root activity" dimension.
    """
    min_amount = 0.0
    if store_inventory:
        trails.reset_inventory(attribute_to_roots=attribute_to_roots)
    else:
        trails.reset_scores(attribute_to_roots=attribute_to_roots, methods=methods)

    cf = None
    if compute_score:
        cf = get_cf_vector(
            trails=trails,
            methods=methods,
            char_cache=_CHAR_CACHE,
            debug=debug,
            ei_version=ei_version,
        )

    # Fail fast with a precise message if inventory builders are not live
    # Fail fast: chunk-based inventory builder must be initialized
    # Only require inventory builders if we will store inventory
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

    fu0 = int(start_act_idx)
    y0 = int(start_year)
    amt0 = float(amount)

    # Traversal
    if attribute_to_roots:
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
            debug=debug,
        )
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        if injected_supply_prov_by_year_act is None:
            injected_supply_prov_by_year_act = {}
    else:
        frontier, injected_supply_by_year_act = trails.temporal_traversal(
            start_year=y0,
            start_act_idx=fu0,
            amount=amt0,
            max_depth=int(max_depth),
            min_amount=float(min_amount),
            return_provenance=False,
            show_progress=bool(show_progress),
            debug=debug,
        )
        provenance = {}
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        injected_supply_prov_by_year_act = {}

    # Always inject FU directly
    injected_supply_by_year_act[(y0, fu0)] = (
        float(injected_supply_by_year_act.get((y0, fu0), 0.0)) + amt0
    )
    injected_supply_prov_by_year_act.setdefault((y0, fu0), {})
    injected_supply_prov_by_year_act[(y0, fu0)][fu0] = (
        float(injected_supply_prov_by_year_act[(y0, fu0)].get(fu0, 0.0)) + amt0
    )

    # Frontier -> demand vectors (calendar years preserved)
    f_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(f_by_year.keys())

    dp_cache: Dict[tuple, Any] = {}

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
        roots_arr = None
        X_roots = None

        solve_year = int(solve_year)
        arr = np.asarray(f_by_year[solve_year])
        nz_idx = np.where(arr != 0.0)[0]
        if nz_idx.size == 0:
            if pbar is not None:
                pbar.update(1)
            continue

        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=solve_year,
            zero_bio=True,
            debug=debug,
        )

        activity_demand = {int(i): float(arr[i]) for i in nz_idx}

        lca_obj = bc.LCA(demand=activity_demand, data_objs=[dp])
        lca_obj.load_lci_data()  # build matrices + dicts

        # Cache mappings once per year
        act_map = getattr(lca_obj.dicts, "activity", None)
        prod_map = getattr(lca_obj.dicts, "product", None)

        act_ids, positions = _get_mapping_arrays(act_map)
        # For RHS building we need the product id -> row index dict
        product_dict = prod_map  # mapping-like

        fu_demand = _map_activity_demands_to_products(
            lca_obj,
            activity_demand,
        )
        if not fu_demand:
            if pbar is not None:
                pbar.update(1)
            continue

        # Prepare A (CSC) once per year; UMFPACK expects CSC
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
            roots, B = _build_rhs_matrix_from_root_demands(
                per_root_demands=per_root_demands,
                product_dict=product_dict,
                n=A_csc.shape[0],
                min_amount=float(min_amount),
            )

            if roots:
                if bc.PYPARDISO:
                    # Keep old path for PARDISO; your Mac path is UMFPACK so this usually won't run
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
                    # UMFPACK: factorize once and solve all RHS vectors
                    X = solve_many_rhs_umfpack_factorized(A_csc, B)
                    roots_arr = np.asarray(
                        roots, dtype=np.int64
                    )  # column order of B and X
                    X_roots = X

                    for j, root_act in enumerate(roots):
                        supply_vec = X[:, j]
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
            lca_obj.build_demand_array(fu_demand)
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

        # Injected supply
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

        # ---- Inventory (keep your existing path for now) ----
        if supplies and store_inventory:
            for supply_dict, root_act in supplies:
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=solve_year,
                    supply_by_activity=supply_dict,
                    min_amount=float(min_amount),
                    store_activity=root_act,
                    debug=debug,
                )

        # ---- Scores: NEW dense multi-root path ----
        if compute_score:
            # Build dense supply matrix (n_acts, n_roots_this_year)
            # We only include "root demands" supplies here. Injected supplies are already dicts;
            # keep them as-is (small) or merge them into the dense matrix if you want.
            if attribute_to_roots:
                # Only call the matrix scorer if we actually solved a multi-RHS system this year
                # AND the roots ordering matches X columns.
                if roots_arr is not None and X_roots is not None:
                    trails.accumulate_temporalized_biosphere_score_matrix(
                        base_year=solve_year,
                        supply_matrix=X_roots,
                        root_activities=roots_arr,
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
                # You can keep the dict path or extend matrix scorer for no-root.
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

    # Characterized inventory is optional now
    if store_inventory and (not compute_score):
        build_characterized_inventory(
            trails=trails, methods=methods, char_cache=_CHAR_CACHE
        )
