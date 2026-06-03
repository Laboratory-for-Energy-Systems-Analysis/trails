from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DEFAULT_DATAPACKAGE = Path(
    "/Users/romain/GitHub/trails/dev/trails_remind_SSP2-PkBudg1000.zip"
)
DEFAULT_LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")
DEFAULT_OUTPUT_CSV = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_depth_sweep.csv"
)
DEFAULT_METHOD = "EF v3.1 - particulate matter formation - impact on human health"
DEFAULT_REFERENCE_YEAR = 2025
DEFAULT_AMOUNT = 20_000_000_000.0
DEFAULT_DEPTHS = [1, 2, 3, 4, 5, 6]
DEFAULT_ROUTING_MIN_AMOUNT = 1e-3


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "plot_terminal_lci_td_comparison", RUNNER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load runner script: {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


runner = _load_runner()
helpers = runner.helpers
DEFAULT_ACTIVITY = runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS["daccs"]
LOCAL_INVENTORY_PATHS = [
    REPO_ROOT / "dev" / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    REPO_ROOT / "dev" / "lci-case-study-daccs_storage_risk.xlsx",
    REPO_ROOT / "dev" / "lci-case-study-marine_fuel_switch.xlsx",
    REPO_ROOT / "dev" / "lci-pass_cars.xlsx",
]
DEFAULT_INVENTORY_PATHS = (
    LOCAL_INVENTORY_PATHS
    if all(path.exists() for path in LOCAL_INVENTORY_PATHS)
    else list(runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS)
)


def _validate_paths(datapackage: Path, inventory_paths: list[Path]) -> None:
    datapackage = Path(datapackage).expanduser()
    if not datapackage.exists():
        raise FileNotFoundError(
            "Datapackage not found: "
            f"{datapackage}\nPass --datapackage PATH to use another package."
        )

    missing = [path for path in inventory_paths if not Path(path).expanduser().exists()]
    if missing:
        raise FileNotFoundError(
            "Missing inventory file(s):\n- "
            + "\n- ".join(str(path) for path in missing)
        )


def _method_unit(method: str, ei_version: str) -> str:
    from trails.lcia import _get_lcia_methods_filepath

    path = _get_lcia_methods_filepath(str(ei_version))
    with path.open("r", encoding="utf-8") as handle:
        import json

        data = json.load(handle)
    return {
        " - ".join(item["name"]): str(item.get("unit") or "")
        for item in data
    }.get(method, "")


def _graph_stats(trails: Any) -> dict[str, int]:
    graph = getattr(trails, "graph", None)
    if graph is None:
        return {
            "graph_nodes": 0,
            "graph_edges": 0,
            "graph_max_depth": 0,
            "frontier_nodes": 0,
            "direct_bio_nodes": 0,
        }

    max_depth = 0
    frontier_nodes = 0
    direct_bio_nodes = 0
    nodes_by_depth: dict[int, int] = {}
    for _, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        max_depth = max(max_depth, depth)
        nodes_by_depth[depth] = nodes_by_depth.get(depth, 0) + 1
        if data.get("frontier_amount"):
            frontier_nodes += 1
        if data.get("direct_bio_amount"):
            direct_bio_nodes += 1

    out = {
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "graph_max_depth": int(max_depth),
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
    }
    for depth, count in sorted(nodes_by_depth.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    return out


def _select_method(data: Any, method: str) -> Any:
    if "method" not in data.dims:
        return data
    methods = [str(value) for value in data.coords["method"].values.tolist()]
    try:
        position = methods.index(str(method))
    except ValueError as exc:
        raise ValueError(f"Method {method!r} not present in score tensor.") from exc
    return data.isel(method=position, drop=True)


def _temporal_total_score(trails: Any, method: str) -> float:
    if getattr(trails, "scores", None) is None:
        raise RuntimeError("trails.scores is missing after temporal LCA.")

    data = _select_method(trails.scores, method)
    reduce_dims = list(data.dims)
    if reduce_dims:
        data = data.sum(dim=reduce_dims)
    values = data.data
    if hasattr(values, "todense"):
        return float(np.asarray(values.todense(), dtype=float).sum())
    return float(np.asarray(data.values, dtype=float).sum())


def _score_to_float(score: object) -> float:
    values = np.asarray(score, dtype=float).ravel()
    if values.size != 1:
        raise ValueError(f"Expected one static score, got {values.size}.")
    return float(values[0])


def _run_static_lca(
    trails: Any,
    *,
    activity_index: int,
    method: str,
    reference_year: int,
    amount: float,
    ei_version: str,
) -> tuple[float, float]:
    start = time.perf_counter()
    trails.static_lca(
        year=int(reference_year),
        act_idx=int(activity_index),
        amount=float(amount),
        methods=[method],
        ei_version=str(ei_version),
    )
    seconds = time.perf_counter() - start
    return _score_to_float(trails.static_score), seconds


def _run_temporal_lca(
    trails: Any,
    *,
    method: str,
    solver_mode: str,
    fallback_solver_mode: str,
    iterative_rtol: float,
    iterative_atol: float,
    iterative_restart: int | None,
    iterative_maxiter: int | None,
    iterative_use_guess: bool,
    iterative_preconditioner: str,
    iterative_ilu_drop_tol: float,
    iterative_ilu_fill_factor: float,
    show_progress: bool,
    ei_version: str,
) -> tuple[float, str]:
    def call(mode: str) -> None:
        trails.lca(
            methods=[method],
            show_progress=bool(show_progress),
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode=str(mode),
            iterative_rtol=float(iterative_rtol),
            iterative_atol=float(iterative_atol),
            iterative_restart=iterative_restart,
            iterative_maxiter=iterative_maxiter,
            iterative_use_guess=bool(iterative_use_guess),
            iterative_preconditioner=str(iterative_preconditioner),
            iterative_ilu_drop_tol=float(iterative_ilu_drop_tol),
            iterative_ilu_fill_factor=float(iterative_ilu_fill_factor),
            ei_version=str(ei_version),
        )

    start = time.perf_counter()
    try:
        call(str(solver_mode))
        return time.perf_counter() - start, str(solver_mode)
    except RuntimeError as exc:
        fallback = str(fallback_solver_mode).strip().lower()
        if fallback in {"", "none"} or fallback == str(solver_mode).strip().lower():
            raise
        print(
            f"    solver_mode={solver_mode!r} failed: {exc}; "
            f"retrying with solver_mode={fallback_solver_mode!r}",
            flush=True,
        )
        call(str(fallback_solver_mode))
        return time.perf_counter() - start, str(fallback_solver_mode)


def _deviation(score: float, static_score: float) -> tuple[float, float]:
    absolute = float(score) - float(static_score)
    if static_score == 0:
        return absolute, float("nan")
    return absolute, absolute / float(static_score)


def _write_csv(rows: list[dict[str, Any]], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float_or_nan(value: Any) -> float:
    try:
        text = str(value).strip()
        if text == "":
            return float("nan")
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def _int_or_none(value: Any) -> int | None:
    try:
        text = str(value).strip()
        if text == "":
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _existing_static_row(
    rows: list[dict[str, Any]],
    *,
    method: str,
    activity_index: int,
    reference_year: int,
    amount: float,
) -> dict[str, Any] | None:
    for row in rows:
        if str(row.get("mode", "")).strip() != "static":
            continue
        if str(row.get("method", "")).strip() != str(method):
            continue
        if _int_or_none(row.get("activity_index")) != int(activity_index):
            continue
        if _int_or_none(row.get("reference_year")) != int(reference_year):
            continue
        if not np.isclose(_float_or_nan(row.get("amount")), float(amount)):
            continue
        return row
    return None


def _existing_depths(
    rows: list[dict[str, Any]],
    *,
    method: str,
    activity_index: int,
    reference_year: int,
    amount: float,
) -> set[int]:
    depths: set[int] = set()
    for row in rows:
        if str(row.get("mode", "")).strip() != "temporal":
            continue
        if str(row.get("method", "")).strip() != str(method):
            continue
        if _int_or_none(row.get("activity_index")) != int(activity_index):
            continue
        if _int_or_none(row.get("reference_year")) != int(reference_year):
            continue
        if not np.isclose(_float_or_nan(row.get("amount")), float(amount)):
            continue
        depth = _int_or_none(row.get("depth"))
        if depth is not None:
            depths.add(depth)
    return depths


def run(args: argparse.Namespace) -> int:
    datapackage = Path(args.datapackage).expanduser().resolve()
    inventory_paths = [Path(path).expanduser().resolve() for path in args.inventories]
    _validate_paths(datapackage, inventory_paths)
    output_csv = Path(args.output_csv).expanduser().resolve()

    lcia_json = None if args.lcia_json is None else Path(args.lcia_json).expanduser()
    if lcia_json is not None:
        if not lcia_json.exists():
            raise FileNotFoundError(f"LCIA JSON not found: {lcia_json}")
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json.resolve())

    available = runner.get_lcia_method_names(ei_version=str(args.ei_version))
    if args.method not in available:
        raise ValueError(
            f"LCIA method not found for ecoinvent {args.ei_version}: {args.method}"
        )

    print(f"Datapackage: {datapackage}", flush=True)
    print("Inventory files:", flush=True)
    for path in inventory_paths:
        print(f" - {path}", flush=True)
    print(f"Method: {args.method}", flush=True)

    load_start = time.perf_counter()
    trails = runner._load_trails(
        datapackage=datapackage,
        interpolation_cache_dir=None,
        inventory_paths=inventory_paths,
        import_before_interpolation=bool(args.import_before_interpolation),
        remove_base_temporal_distributions=False,
        no_cache_interpolation=bool(args.no_cache_interpolation),
        interpolation_start_year_offset=int(args.interpolation_start_year_offset),
        interpolation_end_year_offset=int(args.interpolation_end_year_offset),
    )
    load_seconds = time.perf_counter() - load_start

    activity_maps = helpers._match_activity_indices(trails, [DEFAULT_ACTIVITY])
    if DEFAULT_ACTIVITY not in activity_maps:
        raise RuntimeError(
            "Could not match the notebook DAC activity after inventory import: "
            f"{DEFAULT_ACTIVITY}"
        )
    activity_index = int(activity_maps[DEFAULT_ACTIVITY])
    activity_label = helpers._activity_label(
        trails, activity_index, int(args.reference_year)
    )
    print(f"Activity: {activity_label} (idx={activity_index})", flush=True)

    unit = _method_unit(args.method, str(args.ei_version))
    rows: list[dict[str, Any]] = _read_csv_rows(output_csv) if args.append else []

    existing_static = (
        _existing_static_row(
            rows,
            method=args.method,
            activity_index=activity_index,
            reference_year=int(args.reference_year),
            amount=float(args.amount),
        )
        if args.append
        else None
    )
    if existing_static is not None:
        static_score = _float_or_nan(
            existing_static.get("static_score") or existing_static.get("score")
        )
        static_seconds = _float_or_nan(existing_static.get("static_lca_seconds"))
        print(
            f"Reusing existing static score={static_score:.12g} "
            f"from {output_csv}",
            flush=True,
        )
    else:
        print("Running static LCA", flush=True)
        static_score, static_seconds = _run_static_lca(
            trails,
            activity_index=activity_index,
            method=args.method,
            reference_year=int(args.reference_year),
            amount=float(args.amount),
            ei_version=str(args.ei_version),
        )
        rows.append(
            {
                "mode": "static",
                "depth": "",
                "activity_index": activity_index,
                "activity": activity_label,
                "reference_year": int(args.reference_year),
                "amount": float(args.amount),
                "method": args.method,
                "unit": unit,
                "score": static_score,
                "static_score": static_score,
                "score_deviation_from_static": 0.0,
                "relative_deviation_from_static": 0.0,
                "load_seconds": load_seconds,
                "static_lca_seconds": static_seconds,
                "routing_seconds": "",
                "temporal_lca_seconds": "",
                "total_depth_step_seconds": static_seconds,
                "solver_mode": "",
                "routing_min_amount": "",
                "graph_nodes": "",
                "graph_edges": "",
                "graph_max_depth": "",
                "frontier_nodes": "",
                "direct_bio_nodes": "",
            }
        )
        print(
            f"  static score={static_score:.12g} in {static_seconds:.1f}s",
            flush=True,
        )

    existing_depths = _existing_depths(
        rows,
        method=args.method,
        activity_index=activity_index,
        reference_year=int(args.reference_year),
        amount=float(args.amount),
    )

    for depth in args.depths:
        depth = int(depth)
        if args.append and not args.force and depth in existing_depths:
            print(f"Skipping existing temporal depth {depth}", flush=True)
            continue

        print(f"Running temporal depth {depth}", flush=True)
        step_start = time.perf_counter()

        routing_start = time.perf_counter()
        trails.temporal_routing(
            start_year=int(args.reference_year),
            start_act_idx=activity_index,
            amount=float(args.amount),
            max_depth=depth,
            min_amount=float(args.routing_min_amount),
            show_progress=bool(args.show_progress),
            attribute_to_roots=True,
        )
        routing_seconds = time.perf_counter() - routing_start
        graph_stats = _graph_stats(trails)
        print(
            "  routing "
            f"{routing_seconds:.1f}s, nodes={graph_stats['graph_nodes']:,}, "
            f"edges={graph_stats['graph_edges']:,}",
            flush=True,
        )

        lca_seconds, actual_solver = _run_temporal_lca(
            trails,
            method=args.method,
            solver_mode=str(args.solver_mode),
            fallback_solver_mode=str(args.fallback_solver_mode),
            iterative_rtol=float(args.iterative_rtol),
            iterative_atol=float(args.iterative_atol),
            iterative_restart=args.iterative_restart,
            iterative_maxiter=args.iterative_maxiter,
            iterative_use_guess=bool(args.iterative_use_guess),
            iterative_preconditioner=str(args.iterative_preconditioner),
            iterative_ilu_drop_tol=float(args.iterative_ilu_drop_tol),
            iterative_ilu_fill_factor=float(args.iterative_ilu_fill_factor),
            show_progress=bool(args.show_progress),
            ei_version=str(args.ei_version),
        )
        score = _temporal_total_score(trails, args.method)
        absolute_deviation, relative_deviation = _deviation(score, static_score)
        step_seconds = time.perf_counter() - step_start

        row = {
            "mode": "temporal",
            "depth": depth,
            "activity_index": activity_index,
            "activity": activity_label,
            "reference_year": int(args.reference_year),
            "amount": float(args.amount),
            "method": args.method,
            "unit": unit,
            "score": score,
            "static_score": static_score,
            "score_deviation_from_static": absolute_deviation,
            "relative_deviation_from_static": relative_deviation,
            "load_seconds": load_seconds,
            "static_lca_seconds": static_seconds,
            "routing_seconds": routing_seconds,
            "temporal_lca_seconds": lca_seconds,
            "total_depth_step_seconds": step_seconds,
            "solver_mode": actual_solver,
            "routing_min_amount": float(args.routing_min_amount),
            **graph_stats,
        }
        rows.append(row)
        _write_csv(rows, output_csv)
        print(
            f"  score={score:.12g}, deviation={absolute_deviation:.12g}, "
            f"relative={relative_deviation:.6g}, lca={lca_seconds:.1f}s",
            flush=True,
        )

        trails.inventory = None
        trails.characterized_inventory = None
        trails.scores = None
        trails.graph = None
        gc.collect()

    _write_csv(rows, output_csv)
    print(f"Wrote CSV: {output_csv}", flush=True)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Assess the notebook DAC activity with the EF v3.1 particulate "
            "matter indicator for static LCA and temporal routing depths."
        )
    )
    parser.add_argument("--datapackage", type=Path, default=DEFAULT_DATAPACKAGE)
    parser.add_argument(
        "--inventories",
        type=Path,
        nargs="+",
        default=DEFAULT_INVENTORY_PATHS,
        help="Excel inventories to import, defaulting to the depth-sweep notebook set.",
    )
    parser.add_argument("--lcia-json", type=Path, default=DEFAULT_LCIA_JSON)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    parser.add_argument("--method", default=DEFAULT_METHOD)
    parser.add_argument("--ei-version", default="3.12")
    parser.add_argument("--reference-year", type=int, default=DEFAULT_REFERENCE_YEAR)
    parser.add_argument("--amount", type=float, default=DEFAULT_AMOUNT)
    parser.add_argument("--depths", type=int, nargs="+", default=DEFAULT_DEPTHS)
    parser.add_argument(
        "--routing-min-amount",
        type=float,
        default=DEFAULT_ROUTING_MIN_AMOUNT,
    )
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--no-cache-interpolation", action="store_true")
    parser.add_argument("--import-before-interpolation", action="store_true")
    parser.add_argument("--interpolation-start-year-offset", type=int, default=-20)
    parser.add_argument("--interpolation-end-year-offset", type=int, default=20)
    parser.add_argument("--solver-mode", default="iterative")
    parser.add_argument("--fallback-solver-mode", default="direct")
    parser.add_argument("--iterative-rtol", type=float, default=1e-3)
    parser.add_argument("--iterative-atol", type=float, default=0.0)
    parser.add_argument("--iterative-restart", type=int, default=100)
    parser.add_argument("--iterative-maxiter", type=int, default=1000)
    parser.add_argument(
        "--iterative-use-guess",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--iterative-preconditioner", default="jacobi")
    parser.add_argument("--iterative-ilu-drop-tol", type=float, default=1e-4)
    parser.add_argument("--iterative-ilu-fill-factor", type=float, default=10.0)
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append new temporal depths to an existing CSV and reuse its static row.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="With --append, rerun requested depths even if they are already present.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
