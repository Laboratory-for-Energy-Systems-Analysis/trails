from __future__ import annotations

import importlib.util
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import sparse


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

DIAG_PATH = REPO_ROOT / "dev" / "diagnose_polyol_propylene_gwp.py"
OUTPUT_DIR = REPO_ROOT / "dev" / "notebook_runs" / "polyol_propylene_score_jump"

ROOT_IDX = 13395
EMISSION_YEARS = {2024, 2025}
TARGET_ACTIVITIES = [
    31291,  # supply of forest residue | forest residue | EUR
    28125,  # farming and supply of rapeseed | rapeseed cultivation | RER
    29431,  # softwood forestry, spruce | EUR
    29563,  # softwood forestry, pine | EUR
    31298,  # supply of forest residue | forest residue | REF
]


def _load_diag() -> Any:
    spec = importlib.util.spec_from_file_location(
        "diagnose_polyol_propylene_gwp", DIAG_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load diagnostic helpers: {DIAG_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


def _timed(label: str, func: Any, *args: Any, **kwargs: Any) -> Any:
    print(label, flush=True)
    start = time.perf_counter()
    result = func(*args, **kwargs)
    print(f"  done in {time.perf_counter() - start:.1f}s", flush=True)
    return result


def _flow_metadata(trails: Any) -> dict[int, dict[str, Any]]:
    if not trails.biosphere_indices:
        return {}
    first_label = next(iter(trails.biosphere_indices))
    return {
        int(key): value
        for key, value in trails.biosphere_indices[first_label].items()
        if isinstance(value, dict)
    }


def _decode_method_score_dims(trails: Any) -> tuple[np.ndarray, dict[int, int]]:
    years = np.asarray(trails._score_years, dtype=int)
    return years, {int(i): int(year) for i, year in enumerate(years)}


def _install_score_capture(
    trails: Any,
    *,
    root_idx: int,
    target_activities: list[int],
    emission_years: set[int],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], Any]:
    captured_scores: list[dict[str, Any]] = []
    captured_supplies: list[dict[str, Any]] = []
    current_base_year: list[int | None] = [None]

    original_append = trails._append_scores_bulk
    original_matrix_multi = trails.accumulate_temporalized_biosphere_score_matrix_multi
    original_single = trails.accumulate_temporalized_biosphere_score

    target_arr = np.asarray(target_activities, dtype=np.int64)

    def append_wrapper(
        act_idx: np.ndarray,
        year_idx: np.ndarray,
        values: np.ndarray,
        *,
        root_activity: np.ndarray | None = None,
        method_idx: int | np.ndarray | None = None,
    ) -> None:
        if root_activity is not None:
            years, idx_to_year = _decode_method_score_dims(trails)
            acts = np.asarray(act_idx, dtype=np.int64)
            yidx = np.asarray(year_idx, dtype=np.int64)
            vals = np.asarray(values, dtype=np.float64)
            roots = np.asarray(root_activity, dtype=np.int64)
            mask = roots == int(root_idx)
            if mask.any():
                emission = np.array([idx_to_year[int(i)] for i in yidx], dtype=np.int64)
                mask &= np.isin(emission, np.fromiter(emission_years, dtype=np.int64))
                if mask.any():
                    for a, y, v in zip(acts[mask], emission[mask], vals[mask]):
                        captured_scores.append(
                            {
                                "base_year": current_base_year[0],
                                "emission_year": int(y),
                                "root": int(root_idx),
                                "activity": int(a),
                                "score": float(v),
                            }
                        )

        return original_append(
            act_idx,
            year_idx,
            values,
            root_activity=root_activity,
            method_idx=method_idx,
        )

    def matrix_multi_wrapper(
        *,
        base_year: int,
        supply_matrix: np.ndarray,
        root_activities: np.ndarray,
        cf_matrix: np.ndarray,
        method_indices: np.ndarray | None = None,
        min_amount: float = 0.0,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> Any:
        roots = np.asarray(root_activities, dtype=np.int64)
        root_pos = np.where(roots == int(root_idx))[0]
        if root_pos.size:
            column = int(root_pos[0])
            supplies = np.asarray(supply_matrix, dtype=np.float64)
            for activity in target_arr:
                amount = float(supplies[int(activity), column])
                if amount:
                    captured_supplies.append(
                        {
                            "base_year": int(base_year),
                            "root": int(root_idx),
                            "activity": int(activity),
                            "supply_amount": amount,
                        }
                    )
        previous = current_base_year[0]
        current_base_year[0] = int(base_year)
        try:
            return original_matrix_multi(
                base_year=base_year,
                supply_matrix=supply_matrix,
                root_activities=root_activities,
                cf_matrix=cf_matrix,
                method_indices=method_indices,
                min_amount=min_amount,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
            )
        finally:
            current_base_year[0] = previous

    def single_wrapper(
        base_year: int,
        supply_by_activity: dict[int, float],
        cf: np.ndarray,
        *,
        min_amount: float = 0.0,
        store_activity: int | None = None,
        use_temporal_distributions: bool = True,
        debug: bool = False,
        method_idx: int | None = None,
    ) -> Any:
        previous = current_base_year[0]
        current_base_year[0] = int(base_year)
        try:
            return original_single(
                base_year,
                supply_by_activity,
                cf,
                min_amount=min_amount,
                store_activity=store_activity,
                use_temporal_distributions=use_temporal_distributions,
                debug=debug,
                method_idx=method_idx,
            )
        finally:
            current_base_year[0] = previous

    trails._append_scores_bulk = append_wrapper
    trails.accumulate_temporalized_biosphere_score_matrix_multi = matrix_multi_wrapper
    trails.accumulate_temporalized_biosphere_score = single_wrapper

    def restore() -> None:
        trails._append_scores_bulk = original_append
        trails.accumulate_temporalized_biosphere_score_matrix_multi = (
            original_matrix_multi
        )
        trails.accumulate_temporalized_biosphere_score = original_single

    return captured_scores, captured_supplies, restore


def _td_offsets(trails: Any, tex: Any) -> list[tuple[int, float]]:
    if tex is None:
        return [(0, 1.0)]
    return [
        (int(offset), float(weight))
        for offset, weight in trails._get_td_offsets(tex=tex, debug=False)
    ]


def _decompose_target_flows(
    trails: Any,
    *,
    supplies: pd.DataFrame,
    emission_years: set[int],
    target_activities: list[int],
    method: str,
) -> pd.DataFrame:
    from trails.characterization import get_cf_vector

    cf = get_cf_vector(
        trails=trails,
        methods=[method],
        char_cache={},
        ei_version="3.12",
    )
    flows_meta = _flow_metadata(trails)

    rows: list[dict[str, Any]] = []
    supplies = supplies[supplies["activity"].isin(target_activities)].copy()

    for row in supplies.itertuples(index=False):
        base_year = int(row.base_year)
        activity = int(row.activity)
        supply_amount = float(row.supply_amount)
        context = trails._get_scenario_context(base_year)
        if context is None:
            continue
        _scenario_year, scenario_label, t = context
        b_row = trails.B[int(t), activity, :]
        if getattr(b_row, "nnz", 0) == 0:
            continue
        for flow_idx, value in zip(b_row.coords[0], b_row.data):
            flow_idx = int(flow_idx)
            value = float(value)
            factor = float(cf[flow_idx])
            if factor == 0.0 or value == 0.0:
                continue
            tex = trails.temporal_biosphere_exchanges.get(
                (str(trails._map_year_to_template_year(base_year)), activity, flow_idx)
            )
            amount_source = getattr(tex, "amount_source", "port") if tex else "none"
            for offset, weight in _td_offsets(trails, tex):
                raw_year = int(base_year + offset)
                emission_year = int(trails._clamp_year_to_scores(raw_year))
                if emission_year not in emission_years:
                    continue

                # Matrix-sourced temporal biosphere exchanges are uncommon in these
                # target activities; read the effective-year B value if needed.
                effective_value = value
                if tex is not None and amount_source == "matrix":
                    effective_context = trails._get_scenario_context(
                        int(trails._map_year_to_scenario_year(raw_year))
                    )
                    if effective_context is None:
                        continue
                    _, _, t_eff = effective_context
                    b_eff = trails.B[int(t_eff), activity, :]
                    if getattr(b_eff, "nnz", 0) == 0:
                        continue
                    matches = np.where(b_eff.coords[0].astype(int) == flow_idx)[0]
                    if not matches.size:
                        continue
                    effective_value = float(b_eff.data[int(matches[0])])

                contribution = supply_amount * effective_value * factor * float(weight)
                if contribution == 0.0:
                    continue
                meta = flows_meta.get(flow_idx, {})
                rows.append(
                    {
                        "base_year": base_year,
                        "emission_year": emission_year,
                        "activity": activity,
                        "supply_amount": supply_amount,
                        "flow": flow_idx,
                        "flow_name": meta.get("name", ""),
                        "flow_compartment": meta.get("compartment", ""),
                        "flow_subcompartment": meta.get("subcompartment", ""),
                        "flow_amount_per_supply": effective_value,
                        "cf": factor,
                        "td_distribution": getattr(tex, "distribution", None),
                        "td_min": getattr(tex, "offset_min", None),
                        "td_max": getattr(tex, "offset_max", None),
                        "offset": int(offset),
                        "weight": float(weight),
                        "score": float(contribution),
                    }
                )

    return pd.DataFrame(rows)


def main() -> None:
    diag = _load_diag()
    if diag.LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(diag.LCIA_JSON)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    runner = diag._load_runner()
    runner._validate_paths(diag.DATAPACKAGE, runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS)

    trails = _timed(
        "Loading Trails and importing OneDrive foreground inventories",
        runner._load_trails,
        datapackage=diag.DATAPACKAGE,
        interpolation_cache_dir=None,
        inventory_paths=[
            path.resolve() for path in runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS
        ],
        import_before_interpolation=False,
        remove_base_temporal_distributions=False,
        no_cache_interpolation=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )
    by_idx = diag._metadata_by_idx(trails)
    activity_index = diag._find_activity(
        by_idx,
        name=diag.ACTIVITY_NAME,
        reference_product=diag.REFERENCE_PRODUCT,
        location=diag.LOCATION,
    )

    _timed(
        "Temporal routing",
        trails.temporal_routing,
        start_year=diag.REFERENCE_YEAR,
        start_act_idx=activity_index,
        amount=diag.FUNCTIONAL_UNIT_AMOUNT,
        max_depth=diag.DEPTH,
        min_amount=diag.MIN_AMOUNT,
        show_progress=False,
        attribute_to_roots=True,
    )

    captured_scores, captured_supplies, restore = _install_score_capture(
        trails,
        root_idx=ROOT_IDX,
        target_activities=TARGET_ACTIVITIES,
        emission_years=EMISSION_YEARS,
    )
    try:
        _timed(
            "Temporal LCA",
            trails.lca,
            methods=[diag.METHOD],
            show_progress=False,
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode="iterative",
            iterative_rtol=1e-3,
            iterative_atol=0.0,
            iterative_restart=100,
            iterative_maxiter=1000,
            iterative_use_guess=True,
            iterative_preconditioner="jacobi",
            iterative_ilu_drop_tol=1e-4,
            iterative_ilu_fill_factor=10.0,
            ei_version="3.12",
        )
    finally:
        restore()

    scores = pd.DataFrame(captured_scores)
    supplies = pd.DataFrame(captured_supplies)
    scores.to_csv(OUTPUT_DIR / "captured_scores_root_13395_2024_2025.csv", index=False)
    supplies.to_csv(OUTPUT_DIR / "captured_supplies_target_activities.csv", index=False)

    if scores.empty:
        raise RuntimeError("No scores captured.")

    activity_meta_rows = []
    for idx in TARGET_ACTIVITIES:
        metadata = by_idx.get(idx, {})
        activity_meta_rows.append(
            {
                "activity": idx,
                "label": diag._activity_label(by_idx, idx),
                "name": metadata.get("name", ""),
                "reference_product": metadata.get("reference product", ""),
                "location": metadata.get("location", ""),
            }
        )
    pd.DataFrame(activity_meta_rows).to_csv(
        OUTPUT_DIR / "target_activity_metadata.csv", index=False
    )

    activity_year = (
        scores.groupby(["emission_year", "activity"], as_index=False)["score"]
        .sum()
        .sort_values(["emission_year", "score"])
    )
    activity_year.to_csv(OUTPUT_DIR / "activity_scores_2024_2025.csv", index=False)

    print("\nCaptured root=propylene activity scores in 2024/2025:", flush=True)
    for year in sorted(EMISSION_YEARS):
        subset = activity_year[activity_year["emission_year"] == year]
        print(f"\nYear {year}, top negatives:", flush=True)
        for row in subset.head(12).itertuples(index=False):
            print(
                f"  {row.score: .6g}  {diag._activity_label(by_idx, int(row.activity))}",
                flush=True,
            )
        print(f"Year {year}, top positives:", flush=True)
        for row in subset.sort_values("score", ascending=False).head(8).itertuples(
            index=False
        ):
            print(
                f"  {row.score: .6g}  {diag._activity_label(by_idx, int(row.activity))}",
                flush=True,
            )

    if not supplies.empty:
        supply_summary = (
            supplies.groupby(["base_year", "activity"], as_index=False)[
                "supply_amount"
            ]
            .sum()
            .sort_values(["activity", "base_year"])
        )
        supply_summary.to_csv(OUTPUT_DIR / "supply_summary_by_base_year.csv", index=False)
        print("\nTarget activity supply amounts under propylene root:", flush=True)
        for activity in TARGET_ACTIVITIES:
            subset = supply_summary[supply_summary["activity"] == activity]
            if subset.empty:
                continue
            print(f"\n{diag._activity_label(by_idx, activity)}", flush=True)
            print(
                "  first/last base year: "
                f"{int(subset.base_year.min())}/{int(subset.base_year.max())}; "
                f"sum={subset.supply_amount.sum():.12g}; "
                f"min={subset.supply_amount.min():.12g}; "
                f"max={subset.supply_amount.max():.12g}",
                flush=True,
            )
            around = subset[
                (subset["base_year"] >= 2024) & (subset["base_year"] <= 2036)
            ]
            if not around.empty:
                print(around.to_string(index=False), flush=True)

    flow_decomp = _decompose_target_flows(
        trails,
        supplies=supplies,
        emission_years=EMISSION_YEARS,
        target_activities=TARGET_ACTIVITIES,
        method=diag.METHOD,
    )
    flow_decomp.to_csv(OUTPUT_DIR / "flow_decomposition_2024_2025.csv", index=False)
    if flow_decomp.empty:
        print("\nNo flow-level decomposition rows produced.", flush=True)
        return

    flow_summary = (
        flow_decomp.groupby(
            [
                "emission_year",
                "activity",
                "flow",
                "flow_name",
                "flow_compartment",
                "flow_subcompartment",
                "cf",
                "td_distribution",
                "td_min",
                "td_max",
            ],
            dropna=False,
            as_index=False,
        )["score"]
        .sum()
        .sort_values(["emission_year", "score"])
    )
    flow_summary.to_csv(OUTPUT_DIR / "flow_summary_2024_2025.csv", index=False)

    print("\nFlow-level explanation for target activities:", flush=True)
    for year in sorted(EMISSION_YEARS):
        print(f"\nYear {year}, top negative flows:", flush=True)
        subset = flow_summary[flow_summary["emission_year"] == year].head(12)
        for row in subset.itertuples(index=False):
            print(
                f"  {row.score: .6g}  act={row.activity} "
                f"{row.flow_name} [{row.flow_compartment}/{row.flow_subcompartment}] "
                f"cf={row.cf:g} TD={row.td_distribution} min={row.td_min} max={row.td_max}",
                flush=True,
            )
        print(f"Year {year}, top positive flows:", flush=True)
        subset = flow_summary[flow_summary["emission_year"] == year].sort_values(
            "score", ascending=False
        ).head(8)
        for row in subset.itertuples(index=False):
            print(
                f"  {row.score: .6g}  act={row.activity} "
                f"{row.flow_name} [{row.flow_compartment}/{row.flow_subcompartment}] "
                f"cf={row.cf:g} TD={row.td_distribution} min={row.td_min} max={row.td_max}",
                flush=True,
            )

    print(f"\nWrote diagnostics to {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
