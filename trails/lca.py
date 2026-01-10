from __future__ import annotations
import sys

import warnings
from typing import Any, Dict, List, TYPE_CHECKING

import bw2calc as bc
import numpy as np
from scikits.umfpack import UmfpackWarning
from tqdm import tqdm

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
                if abs(float(amt)) <= float(min_amount):
                    continue
                root_bucket = year_bucket.setdefault(int(root), {})
                root_bucket[int(a)] = root_bucket.get(int(a), 0.0) + float(amt)

        for (y, a), total_amt in injected_supply_by_year_act.items():
            root_map = injected_supply_prov_by_year_act.get((y, a))
            if not root_map:
                root_map = {int(a): float(total_amt)}
            year_bucket = root_injected_by_year.setdefault(int(y), {})
            for root, amt in root_map.items():
                if abs(float(amt)) <= float(min_amount):
                    continue
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
        solve_year = int(solve_year)
        arr = np.asarray(f_by_year[solve_year])
        nz_idx = np.where(np.abs(arr) > float(min_amount))[0]
        if nz_idx.size == 0:
            if pbar is not None:
                pbar.update(1)
            continue

        fu_demand = {int(i): float(arr[i]) for i in nz_idx}

        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=solve_year,
            zero_bio=True,
            debug=debug,
        )

        lca_obj = bc.LCA(demand=fu_demand, data_objs=[dp])
        lca_obj.lci(factorize=True)

        act_map = getattr(lca_obj.dicts, "activity", None)
        if act_map:
            act_ids = np.fromiter(
                act_map.keys(), dtype=np.int64, count=len(act_map)
            )
            positions = np.fromiter(
                act_map.values(), dtype=np.int64, count=len(act_map)
            )
        else:
            act_ids = None
            positions = None

        if attribute_to_roots:
            supplies: list[tuple[Dict[int, float], int | None]] = []
            per_root_demands = root_demands_by_year.get(solve_year, {})
            for root_act, root_demand in per_root_demands.items():
                lca_obj.redo_lci(demand=root_demand)
                if act_ids is None or positions is None:
                    supply_total = _extract_supply_fast(lca_obj, min_amount)
                else:
                    supply_total = _extract_supply_fast_cached(
                        lca_obj.supply_array, act_ids, positions, min_amount
                    )
                if supply_total:
                    supplies.append((supply_total, int(root_act)))
        else:
            supplies = []
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
                if abs(v) <= float(min_amount):
                    continue
                injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v

            if injected_supply:
                supplies.append((injected_supply, None))

        if supplies:
            for supply_dict, root_act in supplies:
                if store_inventory:
                    trails.accumulate_temporalized_biosphere_inventory(
                        base_year=solve_year,
                        supply_by_activity=supply_dict,
                        min_amount=float(min_amount),
                        store_activity=root_act,
                        debug=debug,
                    )

                if compute_score:
                    trails.accumulate_temporalized_biosphere_score(
                        base_year=solve_year,
                        supply_by_activity=supply_dict,
                        cf=cf,
                        min_amount=float(min_amount),
                        store_activity=root_act,
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
        build_characterized_inventory(trails=trails, methods=methods, char_cache=_CHAR_CACHE)

