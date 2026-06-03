"""Validate Trails against limiting cases for the method manuscript.

The script uses the simple car example datapackage and compares Trails results
against independent dense-matrix calculations. It writes a detailed JSON report
and a compact CSV summary under ``dev/publication`` by default.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import sparse
from datapackage import Package

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from trails import Trails  # noqa: E402
from trails.temporal_distributions import (  # noqa: E402
    TemporalDistribution,
    TemporalExchange,
)

PACKAGE_PATH = REPO_ROOT / "examples" / "example data package" / "datapackage.json"
ROOT_ACTIVITY = 13  # transport, passenger car, ICEV
TOY_CHILD_ACTIVITY = 11  # gasoline production
START_YEAR = 2050
FUNCTIONAL_UNIT_AMOUNT = 1.0
TOY_YEAR = 2030
RTOL = 1e-9
ATOL = 1e-9


def load_example_trails(*, interpolate_annual: bool) -> Trails:
    """Load the example datapackage with deterministic settings."""
    package = Package(str(PACKAGE_PATH))
    return Trails(
        package,
        interpolate_annual=interpolate_annual,
        cache_interpolation=False,
        value_dtype=np.float64,
    )


def clear_internal_caches(trails: Trails) -> None:
    """Clear caches after mutating matrices or temporal exchange metadata."""
    for name in (
        "_tech_td_cache",
        "_tech_td_expanded_cache",
        "_td_offsets_cache",
        "_A_row_cache",
        "_production_amount_cache",
        "_direct_bio_cache_by_year",
    ):
        cache = getattr(trails, name, None)
        if cache is not None:
            cache.clear()
    cache = getattr(trails, "_bio_td_expanded_cache", None)
    if cache is not None:
        cache.clear()


def zero_temporal_exchange_like(tex: TemporalExchange) -> TemporalExchange:
    """Return an explicit zero-offset temporal exchange."""
    return TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        amount_source=getattr(tex, "amount_source", "port"),
        offsets=[0],
        weights=[1.0],
    )


def set_all_temporal_offsets_to_zero(trails: Trails) -> None:
    """Replace every temporal distribution with a single zero-year pulse."""
    trails.temporal_technosphere_exchanges = {
        key: zero_temporal_exchange_like(tex)
        for key, tex in trails.temporal_technosphere_exchanges.items()
    }
    trails.temporal_biosphere_exchanges = {
        key: zero_temporal_exchange_like(tex)
        for key, tex in trails.temporal_biosphere_exchanges.items()
    }
    clear_internal_caches(trails)


def keep_only_foreground_temporal_exchanges(
    trails: Trails,
    *,
    root_activity: int,
) -> None:
    """Keep temporal distributions only on the functional-unit activity."""
    root = int(root_activity)
    trails.temporal_technosphere_exchanges = {
        key: tex if int(key[1]) == root else zero_temporal_exchange_like(tex)
        for key, tex in trails.temporal_technosphere_exchanges.items()
    }
    trails.temporal_biosphere_exchanges = {
        key: tex if int(key[1]) == root else zero_temporal_exchange_like(tex)
        for key, tex in trails.temporal_biosphere_exchanges.items()
    }
    clear_internal_caches(trails)


def set_all_temporal_amount_sources(trails: Trails, source: str) -> None:
    """Set all temporal exchanges to use either anchor or target coefficients."""
    if source not in {"port", "matrix"}:
        raise ValueError("source must be 'port' or 'matrix'")

    def replace(tex: TemporalExchange) -> TemporalExchange:
        return TemporalExchange(
            distribution=tex.distribution,
            loc=tex.loc,
            scale=tex.scale,
            offset_min=tex.offset_min,
            offset_max=tex.offset_max,
            amount_source=source,
            offsets=getattr(tex, "offsets", None),
            weights=getattr(tex, "weights", None),
        )

    trails.temporal_technosphere_exchanges = {
        key: replace(tex) for key, tex in trails.temporal_technosphere_exchanges.items()
    }
    trails.temporal_biosphere_exchanges = {
        key: replace(tex) for key, tex in trails.temporal_biosphere_exchanges.items()
    }
    clear_internal_caches(trails)


def freeze_all_matrix_years_at(trails: Trails, year: int) -> None:
    """Make every scenario-year matrix slice identical to one base year."""
    context = trails._get_scenario_context(int(year))
    if context is None:
        raise ValueError(f"No scenario context for year {year}")
    _, _, t = context

    base_a = np.asarray(trails.A[int(t), :, :].todense(), dtype=np.float64)
    base_b = np.asarray(trails.B[int(t), :, :].todense(), dtype=np.float64)
    n_years = len(trails.scenario_labels)

    trails.A = sparse.COO.from_numpy(np.repeat(base_a[None, :, :], n_years, axis=0))
    trails.B = sparse.COO.from_numpy(np.repeat(base_b[None, :, :], n_years, axis=0))
    clear_internal_caches(trails)


def dense_a_activity_product(trails: Trails, year: int) -> np.ndarray:
    """Return A as activity-by-product for one scenario year."""
    context = trails._get_scenario_context(int(year))
    if context is None:
        raise ValueError(f"No scenario context for year {year}")
    _, _, t = context
    return np.asarray(trails.A[int(t), :, :].todense(), dtype=np.float64)


def dense_b_activity_flow(trails: Trails, year: int) -> np.ndarray:
    """Return B as activity-by-flow for one scenario year."""
    context = trails._get_scenario_context(int(year))
    if context is None:
        raise ValueError(f"No scenario context for year {year}")
    _, _, t = context
    return np.asarray(trails.B[int(t), :, :].todense(), dtype=np.float64)


def static_inventory_by_flow(
    trails: Trails,
    *,
    year: int,
    activity_demands: dict[int, float],
) -> np.ndarray:
    """Independent static matrix inventory by biosphere flow."""
    a_product_activity = dense_a_activity_product(trails, year).T
    demand = np.zeros(a_product_activity.shape[0], dtype=np.float64)
    for act, amount in activity_demands.items():
        demand[int(act)] += float(amount)

    supply = np.linalg.solve(a_product_activity, demand)
    b_activity_flow = dense_b_activity_flow(trails, year)
    return b_activity_flow.T @ supply


def temporal_pulses(tex: TemporalExchange | None) -> list[tuple[int, float]]:
    """Return normalized temporal pulses, or a single zero-offset pulse."""
    if tex is None:
        return [(0, 1.0)]
    return list(TemporalDistribution(tex).iter_offsets_and_weights(debug=False))


def run_trails_inventory(
    trails: Trails,
    *,
    start_year: int,
    root_activity: int,
    amount: float,
    max_depth: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Run Trails and return flow-by-year inventory totals."""
    trails.temporal_routing(
        start_year=int(start_year),
        start_act_idx=int(root_activity),
        amount=float(amount),
        max_depth=int(max_depth),
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
    )
    trails.lca(
        show_progress=False,
        compute_score=False,
        store_inventory=True,
        solver_mode="direct",
    )

    reduced = trails.inventory.sum(dim="activity")
    data = reduced.data
    values = (
        np.asarray(data.todense(), dtype=np.float64)
        if hasattr(data, "todense")
        else np.asarray(data, dtype=np.float64)
    )
    years = np.asarray(reduced.coords["year"].values, dtype=int)
    return years, values


def align_to_years(
    years: np.ndarray,
    values: np.ndarray,
    target_years: list[int],
) -> np.ndarray:
    """Align a flow-by-year matrix to a common year axis."""
    out = np.zeros((values.shape[0], len(target_years)), dtype=np.float64)
    year_pos = {int(year): pos for pos, year in enumerate(years)}
    for target_pos, year in enumerate(target_years):
        source_pos = year_pos.get(int(year))
        if source_pos is not None:
            out[:, target_pos] = values[:, source_pos]
    return out


def dict_inventory_to_arrays(
    inventory: dict[int, np.ndarray],
    *,
    n_flows: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a year -> flow-vector inventory dictionary to arrays."""
    years = np.asarray(sorted(inventory), dtype=int)
    values = np.zeros((n_flows, len(years)), dtype=np.float64)
    for pos, year in enumerate(years):
        values[:, pos] = inventory[int(year)]
    return years, values


def max_abs_diff_for_aligned(
    left_years: np.ndarray,
    left_values: np.ndarray,
    right_years: np.ndarray,
    right_values: np.ndarray,
) -> float:
    """Return maximum absolute difference after aligning year axes."""
    years = sorted(set(map(int, left_years)) | set(map(int, right_years)))
    left = align_to_years(left_years, left_values, years)
    right = align_to_years(right_years, right_values, years)
    return float(np.max(np.abs(left - right))) if years else 0.0


def manual_biosphere_for_activity(
    trails: Trails,
    *,
    base_year: int,
    activity: int,
    supply_amount: float,
    out: dict[int, np.ndarray],
) -> None:
    """Accumulate one activity's biosphere exchanges with temporal logic."""
    b_anchor = dense_b_activity_flow(trails, base_year)
    n_flows = b_anchor.shape[1]

    for flow in range(n_flows):
        anchor_value = float(b_anchor[int(activity), int(flow)])
        if anchor_value == 0.0:
            continue
        tex = trails._get_bio_temporal_exchange(base_year, activity, flow)
        for offset, weight in temporal_pulses(tex):
            raw_year = int(base_year) + int(offset)
            if tex is not None and getattr(tex, "amount_source", "port") == "matrix":
                coefficient = float(
                    dense_b_activity_flow(trails, raw_year)[int(activity), int(flow)]
                )
            else:
                coefficient = anchor_value
            out.setdefault(raw_year, np.zeros(n_flows, dtype=np.float64))
            out[raw_year][flow] += float(supply_amount) * coefficient * float(weight)


def manual_child_demands(
    trails: Trails,
    *,
    base_year: int,
    activity: int,
    supply_amount: float,
) -> dict[int, dict[int, float]]:
    """Expand one foreground technosphere row independently."""
    a_anchor = dense_a_activity_product(trails, base_year)
    out: dict[int, dict[int, float]] = defaultdict(lambda: defaultdict(float))
    _, _, anchor_t = trails._get_scenario_context(base_year)

    for product, exchange_value in enumerate(a_anchor[int(activity), :]):
        exchange_value = float(exchange_value)
        product = int(product)
        if exchange_value == 0.0 or product == int(activity):
            continue

        tex = trails._get_tech_temporal_exchange(base_year, activity, product)
        if tex is not None and getattr(tex, "amount_source", "port") == "matrix":
            for offset, weight in temporal_pulses(tex):
                raw_year = int(base_year) + int(offset)
                context = trails._get_scenario_context(raw_year)
                if context is None:
                    continue
                _, _, target_t = context
                target_value = float(
                    dense_a_activity_product(trails, raw_year)[activity, product]
                )
                production = trails._production_amount(target_t, product)
                child_amount = (
                    float(supply_amount)
                    * abs(target_value)
                    / production
                    * float(weight)
                )
                mapped_year = trails._map_year_to_scenario_year(raw_year)
                out[int(mapped_year)][product] += child_amount
        else:
            production = trails._production_amount(anchor_t, product)
            child_amount_anchor = (
                float(supply_amount) * abs(exchange_value) / production
            )
            for offset, weight in temporal_pulses(tex):
                raw_year = int(base_year) + int(offset)
                mapped_year = trails._map_year_to_scenario_year(raw_year)
                out[int(mapped_year)][product] += child_amount_anchor * float(weight)

    return {year: dict(mapping) for year, mapping in out.items()}


def manual_foreground_only_inventory(
    trails: Trails,
    *,
    start_year: int,
    root_activity: int,
    amount: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Manual reference for foreground-only temporal routing."""
    n_flows = int(trails.B.shape[2])
    out: dict[int, np.ndarray] = {}
    _, _, start_t = trails._get_scenario_context(start_year)
    root_supply = float(amount) / trails._production_amount(start_t, root_activity)

    manual_biosphere_for_activity(
        trails,
        base_year=start_year,
        activity=root_activity,
        supply_amount=root_supply,
        out=out,
    )

    frontier_demands = manual_child_demands(
        trails,
        base_year=start_year,
        activity=root_activity,
        supply_amount=root_supply,
    )
    for year, demands in frontier_demands.items():
        inv = static_inventory_by_flow(
            trails,
            year=int(year),
            activity_demands={int(k): float(v) for k, v in demands.items()},
        )
        out.setdefault(int(year), np.zeros(n_flows, dtype=np.float64))
        out[int(year)] += inv

    return dict_inventory_to_arrays(out, n_flows=n_flows)


def patch_toy_anchor_target_system(trails: Trails, *, source: str) -> None:
    """Patch the example package in memory into the manuscript toy system."""
    if source not in {"port", "matrix"}:
        raise ValueError("source must be 'port' or 'matrix'")

    a_data = np.zeros_like(np.asarray(trails.A.todense(), dtype=np.float64))
    n_activities = a_data.shape[1]
    for t in range(a_data.shape[0]):
        np.fill_diagonal(a_data[t], 1.0)

    target_values = {2030: 10.0, 2031: 8.0, 2033: 6.0}
    for label in trails.scenario_labels:
        year = int(label)
        value = target_values.get(year, 10.0)
        a_data[
            trails.scenario_index[label],
            ROOT_ACTIVITY,
            TOY_CHILD_ACTIVITY,
        ] = value

    b_data = np.zeros_like(np.asarray(trails.B.todense(), dtype=np.float64))
    b_data[:, TOY_CHILD_ACTIVITY, 0] = 1.0

    trails.A = sparse.COO.from_numpy(a_data)
    trails.B = sparse.COO.from_numpy(b_data)

    tex = TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=3,
        amount_source=source,
        offsets=[0, 1, 3],
        weights=[0.5, 0.3, 0.2],
    )
    trails.temporal_technosphere_exchanges = {
        (label, ROOT_ACTIVITY, TOY_CHILD_ACTIVITY): tex
        for label in trails.template_labels
    }
    trails.temporal_biosphere_exchanges = {}

    if n_activities <= max(ROOT_ACTIVITY, TOY_CHILD_ACTIVITY):
        raise RuntimeError("Toy activities are outside the example matrix shape.")
    clear_internal_caches(trails)


def nonzero_years(years: np.ndarray, values: np.ndarray) -> list[int]:
    """Return years with any non-zero flow contribution."""
    mask = np.abs(values).sum(axis=0) > 1e-12
    return [int(year) for year in years[mask]]


def case_static_equivalence() -> dict[str, Any]:
    """Case 1: zero temporal offsets and a single matrix year reproduce static LCA."""
    trails = load_example_trails(interpolate_annual=False)
    set_all_temporal_offsets_to_zero(trails)

    observed_years, observed_values = run_trails_inventory(
        trails,
        start_year=START_YEAR,
        root_activity=ROOT_ACTIVITY,
        amount=FUNCTIONAL_UNIT_AMOUNT,
        max_depth=3,
    )
    static = static_inventory_by_flow(
        trails,
        year=START_YEAR,
        activity_demands={ROOT_ACTIVITY: FUNCTIONAL_UNIT_AMOUNT},
    )

    expected = np.zeros_like(observed_values)
    year_pos = {int(year): pos for pos, year in enumerate(observed_years)}
    expected[:, year_pos[START_YEAR]] = static

    max_diff = float(np.max(np.abs(observed_values - expected)))
    return {
        "case": "zero_offsets_single_year_static_equivalence",
        "description": (
            "All temporal distributions are replaced by a single offset-zero "
            "pulse and the non-interpolated 2050 matrix is used."
        ),
        "max_abs_diff": max_diff,
        "tolerance": ATOL,
        "status": "pass" if max_diff <= ATOL else "fail",
        "observed_flow_totals": observed_values.sum(axis=1).tolist(),
        "static_reference_flow_totals": static.tolist(),
        "nonzero_years": nonzero_years(observed_years, observed_values),
    }


def case_foreground_only() -> dict[str, Any]:
    """Case 2: first-tier temporal routing matches a manual reference."""
    trails = load_example_trails(interpolate_annual=True)
    keep_only_foreground_temporal_exchanges(
        trails,
        root_activity=ROOT_ACTIVITY,
    )

    expected_years, expected_values = manual_foreground_only_inventory(
        trails,
        start_year=START_YEAR,
        root_activity=ROOT_ACTIVITY,
        amount=FUNCTIONAL_UNIT_AMOUNT,
    )
    observed_years, observed_values = run_trails_inventory(
        trails,
        start_year=START_YEAR,
        root_activity=ROOT_ACTIVITY,
        amount=FUNCTIONAL_UNIT_AMOUNT,
        max_depth=1,
    )

    max_diff = max_abs_diff_for_aligned(
        observed_years,
        observed_values,
        expected_years,
        expected_values,
    )
    return {
        "case": "foreground_only_temporalized",
        "description": (
            "Only exchanges attached to the functional-unit activity are "
            "temporalized; first-tier child demands are solved statically in "
            "their target years."
        ),
        "max_abs_diff": max_diff,
        "tolerance": ATOL,
        "status": "pass" if max_diff <= ATOL else "fail",
        "observed_flow_totals": observed_values.sum(axis=1).tolist(),
        "manual_reference_flow_totals": expected_values.sum(axis=1).tolist(),
        "nonzero_years": nonzero_years(observed_years, observed_values),
    }


def case_identical_matrices_timing_only() -> dict[str, Any]:
    """Case 3: with identical matrices, anchor and target modes coincide."""
    anchor = load_example_trails(interpolate_annual=True)
    freeze_all_matrix_years_at(anchor, START_YEAR)
    set_all_temporal_amount_sources(anchor, "port")
    anchor_years, anchor_values = run_trails_inventory(
        anchor,
        start_year=START_YEAR,
        root_activity=ROOT_ACTIVITY,
        amount=FUNCTIONAL_UNIT_AMOUNT,
        max_depth=3,
    )

    target = load_example_trails(interpolate_annual=True)
    freeze_all_matrix_years_at(target, START_YEAR)
    set_all_temporal_amount_sources(target, "matrix")
    target_years, target_values = run_trails_inventory(
        target,
        start_year=START_YEAR,
        root_activity=ROOT_ACTIVITY,
        amount=FUNCTIONAL_UNIT_AMOUNT,
        max_depth=3,
    )

    max_diff = max_abs_diff_for_aligned(
        anchor_years,
        anchor_values,
        target_years,
        target_values,
    )
    return {
        "case": "identical_matrices_timing_only",
        "description": (
            "All scenario-year A and B slices are replaced by the 2050 slice. "
            "Anchor-year and target-year coefficient modes should then produce "
            "identical time-resolved inventories."
        ),
        "max_abs_diff": max_diff,
        "tolerance": ATOL,
        "status": "pass" if max_diff <= ATOL else "fail",
        "anchor_flow_totals": anchor_values.sum(axis=1).tolist(),
        "target_flow_totals": target_values.sum(axis=1).tolist(),
        "nonzero_years": nonzero_years(anchor_years, anchor_values),
    }


def case_anchor_target_toy() -> dict[str, Any]:
    """Case 4: patched toy system matches hand-calculated pulse values."""
    expected_by_source = {
        "port": {2030: 5.0, 2031: 3.0, 2033: 2.0},
        "matrix": {2030: 5.0, 2031: 2.4, 2033: 1.2},
    }
    observed_by_source: dict[str, dict[int, float]] = {}
    max_diff = 0.0

    for source, expected in expected_by_source.items():
        trails = load_example_trails(interpolate_annual=True)
        patch_toy_anchor_target_system(trails, source=source)
        years, values = run_trails_inventory(
            trails,
            start_year=TOY_YEAR,
            root_activity=ROOT_ACTIVITY,
            amount=1.0,
            max_depth=1,
        )
        observed = {
            int(year): float(values[0, pos])
            for pos, year in enumerate(years)
            if abs(values[0, pos]) > 1e-12
        }
        observed_by_source[source] = observed
        all_years = sorted(set(expected) | set(observed))
        max_diff = max(
            max_diff,
            max(
                abs(float(observed.get(year, 0.0)) - float(expected.get(year, 0.0)))
                for year in all_years
            ),
        )

    return {
        "case": "anchor_vs_target_toy_hand_calculation",
        "description": (
            "The example package is patched in memory to contain one foreground "
            "exchange with coefficient values 10, 8, and 6 in years 2030, "
            "2031, and 2033 and temporal weights 0.5, 0.3, and 0.2."
        ),
        "max_abs_diff": float(max_diff),
        "tolerance": ATOL,
        "status": "pass" if max_diff <= ATOL else "fail",
        "expected_by_source": expected_by_source,
        "observed_by_source": {
            source: {str(year): value for year, value in observed.items()}
            for source, observed in observed_by_source.items()
        },
    }


def write_reports(results: list[dict[str, Any]], output_dir: Path) -> dict[str, Path]:
    """Write JSON and CSV validation summaries."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "validation_limiting_cases_results.json"
    csv_path = output_dir / "validation_limiting_cases_summary.csv"

    payload = {
        "package_path": str(PACKAGE_PATH),
        "root_activity_index": ROOT_ACTIVITY,
        "start_year": START_YEAR,
        "functional_unit_amount": FUNCTIONAL_UNIT_AMOUNT,
        "results": results,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("case", "status", "max_abs_diff", "tolerance", "description"),
        )
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
                    "case": result["case"],
                    "status": result["status"],
                    "max_abs_diff": f"{float(result['max_abs_diff']):.17g}",
                    "tolerance": f"{float(result['tolerance']):.17g}",
                    "description": result["description"],
                }
            )

    return {"json": json_path, "csv": csv_path}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate Trails limiting cases with the simple car datapackage."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for JSON and CSV reports.",
    )
    return parser.parse_args()


def main() -> int:
    """Run all validation cases."""
    args = parse_args()

    results = [
        case_static_equivalence(),
        case_foreground_only(),
        case_identical_matrices_timing_only(),
        case_anchor_target_toy(),
    ]
    paths = write_reports(results, args.output_dir)

    print("Validation results")
    for result in results:
        print(
            f"- {result['case']}: {result['status']} "
            f"(max_abs_diff={float(result['max_abs_diff']):.3e})"
        )
    print(f"JSON: {paths['json']}")
    print(f"CSV:  {paths['csv']}")

    failed = [result for result in results if result["status"] != "pass"]
    if failed:
        print("Failed cases: " + ", ".join(result["case"] for result in failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
