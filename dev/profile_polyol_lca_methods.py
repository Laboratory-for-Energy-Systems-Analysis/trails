from __future__ import annotations

import argparse
import importlib.util
import importlib
import os
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DATAPACKAGE = REPO_ROOT / "dev" / "trails_2026-05-18.zip"
LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")

METHODS = [
    (
        "IPCC 2021 (incl. biogenic CO2) - climate change: total "
        "(incl. biogenic CO2) - global warming potential (GWP100)"
    ),
    "EF v3.1 - ecotoxicity: freshwater - comparative toxic unit for ecosystems (CTUe)",
    "EF v3.1 - human toxicity: carcinogenic - comparative toxic unit for human (CTUh)",
    "EF v3.1 - human toxicity: non-carcinogenic - comparative toxic unit for human (CTUh)",
    (
        "EF v3.1 - material resources: metals/minerals - abiotic depletion "
        "potential (ADP): elements (ultimate reserves)"
    ),
    "EF v3.1 - particulate matter formation - impact on human health",
    (
        "EF v3.1 - water use - user deprivation potential "
        "(deprivation-weighted water consumption)"
    ),
]


def _load_runner() -> Any:
    spec = importlib.util.spec_from_file_location(
        "plot_terminal_lci_td_comparison", RUNNER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load runner script: {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


def _parse_counts(value: str) -> tuple[int, ...]:
    counts = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not counts:
        raise argparse.ArgumentTypeError("At least one method count is required.")
    if any(count < 1 or count > len(METHODS) for count in counts):
        raise argparse.ArgumentTypeError(
            f"Method counts must be between 1 and {len(METHODS)}."
        )
    return counts


def _timed(label: str, func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    print(label, flush=True)
    start = time.perf_counter()
    result = func(*args, **kwargs)
    print(f"  done in {time.perf_counter() - start:.1f}s", flush=True)
    return result


def _wrap_timed_method(obj: Any, name: str, stats: dict[str, dict[str, float]]) -> None:
    original = getattr(obj, name)

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        start = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            elapsed = time.perf_counter() - start
            bucket = stats.setdefault(name, {"calls": 0.0, "seconds": 0.0})
            bucket["calls"] += 1.0
            bucket["seconds"] += elapsed

    setattr(obj, name, wrapped)


def _wrap_timed_function(
    module: Any,
    name: str,
    stats: dict[str, dict[str, float]],
) -> None:
    original = getattr(module, name)

    def wrapped(*args: Any, **kwargs: Any) -> Any:
        start = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            elapsed = time.perf_counter() - start
            bucket = stats.setdefault(name, {"calls": 0.0, "seconds": 0.0})
            bucket["calls"] += 1.0
            bucket["seconds"] += elapsed

    setattr(module, name, wrapped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Profile temporal LCA method scaling for the polyol depth-1 case."
        )
    )
    parser.add_argument("--method-counts", type=_parse_counts, default=(1,))
    parser.add_argument("--depth", type=int, default=1)
    parser.add_argument("--amount", type=float, default=50_000_000_000.0)
    parser.add_argument("--min-amount", type=float, default=1e-3)
    parser.add_argument("--solver-mode", choices=("direct", "iterative"), default="direct")
    args = parser.parse_args()

    if LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(LCIA_JSON)

    runner = _load_runner()
    helpers = runner.helpers
    polyol = helpers.ActivityDef(
        "polyol precursor production from captured CO2",
        "polyol precursor",
        "RER",
    )

    trails = _timed(
        "Loading Trails and importing foreground inventories",
        runner._load_trails,
        datapackage=DATAPACKAGE,
        interpolation_cache_dir=None,
        inventory_paths=[path.resolve() for path in runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS],
        import_before_interpolation=False,
        remove_base_temporal_distributions=False,
        no_cache_interpolation=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )

    activity_index = helpers._match_activity_indices(trails, [polyol])[polyol]
    print(f"Polyol activity index: {activity_index}", flush=True)

    for count in args.method_counts:
        selected = METHODS[: int(count)]
        print(
            f"\nMethod count={len(selected)}, depth={args.depth}, "
            f"amount={args.amount:g}, min_amount={args.min_amount:g}",
            flush=True,
        )
        _timed(
            "Routing",
            trails.temporal_routing,
            start_year=2025,
            start_act_idx=int(activity_index),
            amount=float(args.amount),
            max_depth=int(args.depth),
            min_amount=float(args.min_amount),
            show_progress=False,
            attribute_to_roots=True,
        )
        print(f"  graph: {helpers._routing_graph_summary(trails)}", flush=True)

        stats: dict[str, dict[str, float]] = {}
        _wrap_timed_method(
            trails,
            "accumulate_temporalized_biosphere_score_matrix",
            stats,
        )
        _wrap_timed_method(
            trails,
            "accumulate_temporalized_biosphere_score_matrix_multi",
            stats,
        )
        _wrap_timed_method(
            trails,
            "accumulate_temporalized_biosphere_score",
            stats,
        )
        _wrap_timed_method(trails, "finalize_scores", stats)
        lca_module = importlib.import_module("trails.lca")
        _wrap_timed_function(lca_module, "_build_direct_technosphere_for_year", stats)
        _wrap_timed_function(lca_module, "solve_many_rhs_umfpack_factorized", stats)
        _wrap_timed_function(lca_module, "solve_many_rhs_jacobi_gmres", stats)

        _timed(
            "Temporal LCA",
            trails.lca,
            methods=selected,
            show_progress=False,
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode=str(args.solver_mode),
            ei_version="3.12",
        )
        for key, value in sorted(stats.items()):
            print(
                f"  {key}: calls={int(value['calls'])}, "
                f"seconds={value['seconds']:.1f}",
                flush=True,
            )


if __name__ == "__main__":
    main()
