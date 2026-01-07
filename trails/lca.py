import warnings
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
from .characterization import build_characterized_inventory
from .trails import Trails

warnings.filterwarnings("ignore", category=UmfpackWarning)


def lca_static_simple(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    debug: bool = False,
) -> None:
    """Run a static LCA for a single functional unit and year."""
    trails.reset_inventory()

    dp, _, _, _ = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        debug=debug,
    )

    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    fu_prod_id = _reference_product_id_from_activity_id(lca_obj, int(fu_act_idx))

    lca_obj = bc.LCA(demand={int(fu_prod_id): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    supply_total = _extract_supply_fast(lca_obj, min_amount=0.0)

    trails.accumulate_temporalized_biosphere_inventory(
        base_year=int(year),
        supply_by_activity=supply_total,
        min_amount=0.0,
        use_temporal_distributions=True,
        debug=debug,
    )

    trails.finalize_inventory()
    build_characterized_inventory(trails=trails, methods=methods, char_cache={})


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
) -> None:
    """Run temporal LCA for a functional unit and year."""

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

    def _build_injected_supply(
        injected_supply_by_year_act: dict[tuple[int, int], float],
        solve_year: int,
        min_amount: float,
    ) -> Dict[int, float]:
        injected_supply: Dict[int, float] = {}
        for (y, a), v in injected_supply_by_year_act.items():
            if int(y) != solve_year:
                continue
            v = float(v)
            if abs(v) <= float(min_amount):
                continue
            injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v
        return injected_supply

    if not use_temporal_distributions:
        return lca_static_simple(
            trails=trails,
            year=int(start_year),
            fu_act_idx=int(start_act_idx),
            methods=methods,
            amount=float(amount),
            debug=debug,
        )

    trails.reset_inventory()

    LEGACY_FU_DIRECT_ROOT = -1

    fu0 = int(start_act_idx)
    y0 = int(start_year)
    amt0 = float(amount)

    if return_provenance:
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

    _apply_fu_direct_injection(
        injected_supply_by_year_act=injected_supply_by_year_act,
        injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
        y0=y0,
        fu0=fu0,
        amt0=amt0,
        legacy_root=LEGACY_FU_DIRECT_ROOT,
    )

    f_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(f_by_year.keys())

    dp_cache: Dict[tuple, Any] = {}

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
        if not act_map:
            act_ids = None
            positions = None
        else:
            act_ids = np.fromiter(act_map.keys(), dtype=np.int64, count=len(act_map))
            positions = np.fromiter(
                act_map.values(), dtype=np.int64, count=len(act_map)
            )

        if act_ids is None or positions is None:
            supply_total = _extract_supply_fast(lca_obj, min_amount)
        else:
            supply_total = _extract_supply_fast_cached(
                lca_obj.supply_array, act_ids, positions, min_amount
            )

        trails.accumulate_temporalized_biosphere_inventory(
            base_year=solve_year,
            supply_by_activity=supply_total,
            min_amount=float(min_amount),
            use_temporal_distributions=True,
            debug=debug,
        )

        injected_supply = _build_injected_supply(
            injected_supply_by_year_act=injected_supply_by_year_act,
            solve_year=solve_year,
            min_amount=min_amount,
        )
        if injected_supply:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=injected_supply,
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

    trails.finalize_inventory()
    build_characterized_inventory(trails=trails, methods=methods, char_cache={})

    if return_provenance:
        trails.provenance = provenance
