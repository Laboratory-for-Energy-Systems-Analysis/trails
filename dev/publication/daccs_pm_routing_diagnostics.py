from __future__ import annotations


import argparse
from dataclasses import dataclass
import csv
import json
import gc
import os
import re
import textwrap
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from datapackage import Package
from trails import Trails, get_lcia_method_names
from trails.datapackage import interpolate_to_annual

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

PUBLICATION_DIR = Path(__file__).resolve().parent
DEFAULT_DATAPACKAGE = PUBLICATION_DIR / "trails_remind_SSP2-PkBudg1000.zip"
DEFAULT_LCIA_JSON: Path | None = None


@dataclass(frozen=True)
class ActivityDef:
    name: str
    reference_product: str
    location: str


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _metadata_by_idx(trails: Trails) -> dict[int, dict]:
    if not getattr(trails, "activity_indices", None):
        return {}
    first_label = next(iter(trails.activity_indices))
    return {
        int(key): value for key, value in trails.activity_indices[first_label].items()
    }


def _match_activity_indices(
    trails: Trails,
    targets: list[ActivityDef],
) -> dict[ActivityDef, int]:
    by_idx = _metadata_by_idx(trails)
    out: dict[ActivityDef, int] = {}

    for target in targets:
        exact: list[int] = []
        for index, metadata in by_idx.items():
            if not isinstance(metadata, dict):
                continue
            if _clean(metadata.get("name")) != target.name:
                continue
            if _clean(metadata.get("reference product")) != target.reference_product:
                continue
            if _clean(metadata.get("location")) != target.location:
                continue
            exact.append(index)
        if len(exact) == 1:
            out[target] = exact[0]
            continue

        by_name = [
            index
            for index, metadata in by_idx.items()
            if isinstance(metadata, dict)
            and _clean(metadata.get("name")) == target.name
        ]
        if len(by_name) == 1:
            out[target] = by_name[0]

    return out


def _activity_metadata(
    trails: Trails, activity_index: int, reference_year: int
) -> dict[str, Any]:
    labels: list[str] = []
    try:
        labels.append(str(trails._map_year_to_scenario_year(int(reference_year))))
    except Exception:
        pass
    labels.extend(str(label) for label in getattr(trails, "scenario_labels", []))

    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        mapping = trails.activity_indices.get(label, {})
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)

    for mapping in trails.activity_indices.values():
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)

    return {}


def _activity_label(
    trails: Trails,
    activity_index: int,
    reference_year: int = 2025,
    *,
    max_length: int | None = None,
) -> str:
    meta = _activity_metadata(trails, activity_index, reference_year)
    name = _clean(meta.get("name")) or f"Activity {activity_index}"
    reference_product = _clean(meta.get("reference product"))
    location = _clean(meta.get("location"))
    parts = [name]
    if reference_product:
        parts.append(reference_product)
    if location:
        parts.append(location)
    label = " | ".join(parts)
    if max_length is not None:
        return textwrap.shorten(label, width=max_length, placeholder="...")
    return label


def _clear_base_temporal_distributions(trails: Trails) -> tuple[int, int]:
    tech_count = len(trails.temporal_technosphere_exchanges or {})
    bio_count = len(trails.temporal_biosphere_exchanges or {})
    trails.temporal_technosphere_exchanges = {}
    trails.temporal_biosphere_exchanges = {}
    for cache_name in (
        "_td_offsets_cache",
        "_tech_td_cache",
        "_tech_td_expanded_cache",
        "_direct_bio_cache_by_year",
        "_bio_td_expanded_cache",
        "_bio_score_row_char_cache",
        "_bio_score_row_char_matrix_cache",
    ):
        cache = getattr(trails, cache_name, None)
        if hasattr(cache, "clear"):
            cache.clear()
    return tech_count, bio_count


def _interpolate_trails_after_import(
    trails: Trails,
    *,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> None:
    print("Interpolating foreground-augmented matrices to annual resolution", flush=True)
    trails.A, trails.B, trails.scenario_labels, trails.scenario_index = (
        interpolate_to_annual(
            trails.A,
            trails.B,
            trails.scenario_labels,
            value_dtype=trails.value_dtype,
            start_year_offset=int(interpolation_start_year_offset),
            end_year_offset=int(interpolation_end_year_offset),
        )
    )
    trails.years_int = np.array(
        [int(label) for label in trails.scenario_labels],
        dtype=int,
    )
    trails.min_year = int(trails.years_int.min())
    trails.max_year = int(trails.years_int.max())
    for cache_name in (
        "_td_offsets_cache",
        "_tech_td_cache",
        "_tech_td_expanded_cache",
        "_direct_bio_cache_by_year",
    ):
        cache = getattr(trails, cache_name, None)
        if hasattr(cache, "clear"):
            cache.clear()


def _load_trails(
    *,
    datapackage: Path,
    interpolation_cache_dir: Path | None,
    inventory_paths: list[Path],
    import_before_interpolation: bool,
    remove_base_temporal_distributions: bool,
    no_cache_interpolation: bool,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> Trails:
    if interpolation_cache_dir is not None:
        raise ValueError("Explicit interpolation cache loading is not supported here.")

    trails = Trails(
        Package(str(datapackage)),
        interpolate_annual=not bool(import_before_interpolation),
        cache_interpolation=not bool(no_cache_interpolation),
        interpolation_start_year_offset=int(interpolation_start_year_offset),
        interpolation_end_year_offset=int(interpolation_end_year_offset),
    )
    if remove_base_temporal_distributions:
        tech_count, bio_count = _clear_base_temporal_distributions(trails)
        print(
            "Cleared base temporal distributions before LCI import: "
            f"technosphere={tech_count:,}, biosphere={bio_count:,}",
            flush=True,
        )

    print(
        "Importing foreground inventories together: "
        + ", ".join(path.name for path in inventory_paths),
        flush=True,
    )
    trails.import_excel_inventory([str(path) for path in inventory_paths])
    if import_before_interpolation:
        _interpolate_trails_after_import(
            trails,
            interpolation_start_year_offset=int(interpolation_start_year_offset),
            interpolation_end_year_offset=int(interpolation_end_year_offset),
        )
    return trails


DEFAULT_CASE_STUDY_ACTIVITY_KEYS = {
    "bev": ActivityDef(
        "transport, passenger, car, battery electric",
        "transport, passenger, car",
        "RER",
    ),
    "polyol": ActivityDef(
        "polyol precursor from captured CO2",
        "polyol precursor",
        "RER",
    ),
    "marine": ActivityDef(
        "marine freight service, temporal fuel transition",
        "transport service",
        "RER",
    ),
    "daccs": ActivityDef(
        "carbon dioxide, captured, with a solvent-based direct air capture system, 1MtCO2",
        "carbon dioxide, captured",
        "Europe",
    ),
}

DEFAULT_ACTIVITY = DEFAULT_CASE_STUDY_ACTIVITY_KEYS["daccs"]
LOCAL_INVENTORY_PATHS = [
    PUBLICATION_DIR / "LCIs" / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-case-study-daccs_storage_risk.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-case-study-marine_fuel_switch.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-pass_cars.xlsx",
]
DEFAULT_INVENTORY_PATHS = LOCAL_INVENTORY_PATHS
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
    return {" - ".join(item["name"]): str(item.get("unit") or "") for item in data}.get(
        method, ""
    )


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


def run_depth_sweep(args: argparse.Namespace) -> int:
    datapackage = Path(args.datapackage).expanduser().resolve()
    inventory_paths = [Path(path).expanduser().resolve() for path in args.inventories]
    _validate_paths(datapackage, inventory_paths)
    output_csv = Path(args.output_csv).expanduser().resolve()

    lcia_json = None if args.lcia_json is None else Path(args.lcia_json).expanduser()
    if lcia_json is not None:
        if not lcia_json.exists():
            raise FileNotFoundError(f"LCIA JSON not found: {lcia_json}")
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json.resolve())

    available = get_lcia_method_names(ei_version=str(args.ei_version))
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
    trails = _load_trails(
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

    activity_maps = _match_activity_indices(trails, [DEFAULT_ACTIVITY])
    if DEFAULT_ACTIVITY not in activity_maps:
        raise RuntimeError(
            "Could not match the notebook DAC activity after inventory import: "
            f"{DEFAULT_ACTIVITY}"
        )
    activity_index = int(activity_maps[DEFAULT_ACTIVITY])
    activity_label = _activity_label(
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
            f"Reusing existing static score={static_score:.12g} " f"from {output_csv}",
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


def parse_depth_sweep_args() -> argparse.Namespace:
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


# Adaptive-routing diagnostics share the depth-sweep configuration in this file.
# The adaptive code expects a depth_runner-like module; in this merged file,
# the current module provides those defaults and helper functions.
depth_runner = sys.modules[__name__]
DEFAULT_DEPTH_SWEEP_CSV = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_depth_sweep.csv"
)
DEFAULT_ADAPTIVE_OUTPUT_CSV = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_adaptive_routing.csv"
)
DEFAULT_COMPARISON_CSV = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_adaptive_vs_depth.csv"
)
DEFAULT_SANKEY_DIR = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "sankey"
    / "adaptive_relative_cutoffs"
)
DEFAULT_ACTIVITY_SANKEY_DIR = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "sankey"
    / "adaptive_activity_year_score"
)
DEFAULT_SANKEY_SUMMARY_CSV = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "daccs_pm_adaptive_sankey_depth_summary.csv"
)
DEFAULT_SCORE_PLOT_DIR = (
    REPO_ROOT
    / "dev"
    / "notebook_runs"
    / "daccs_pm_depth_sweep"
    / "temporal_score_plots"
    / "adaptive_relative_cutoffs"
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _float(value: Any) -> float:
    try:
        text = str(value).strip()
        if not text:
            return float("nan")
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def _int(value: Any) -> int | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _depth_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in rows:
        if str(row.get("mode", "")).strip() != "temporal":
            continue
        depth = _int(row.get("depth"))
        if depth is None:
            continue
        parsed = dict(row)
        parsed["_depth"] = depth
        parsed["_score"] = _float(row.get("score"))
        parsed["_routing_seconds"] = _float(row.get("routing_seconds"))
        parsed["_temporal_lca_seconds"] = _float(row.get("temporal_lca_seconds"))
        parsed["_graph_nodes"] = _int(row.get("graph_nodes")) or 0
        parsed["_graph_edges"] = _int(row.get("graph_edges")) or 0
        parsed["_relative_deviation"] = _float(
            row.get("relative_deviation_from_static")
        )
        out.append(parsed)
    return sorted(out, key=lambda row: row["_depth"])


def _static_score(rows: list[dict[str, str]]) -> float | None:
    for row in rows:
        if str(row.get("mode", "")).strip() == "static":
            value = _float(row.get("static_score") or row.get("score"))
            if np.isfinite(value):
                return value
    return None


def _method_unit(method: str, ei_version: str) -> str:
    from trails.lcia import _get_lcia_methods_filepath

    path = _get_lcia_methods_filepath(str(ei_version))
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return {" - ".join(item["name"]): str(item.get("unit") or "") for item in data}.get(
        method, ""
    )


def _graph_stats(trails: Any) -> dict[str, Any]:
    graph = getattr(trails, "graph", None)
    if graph is None:
        return {}

    nodes_by_depth: dict[int, int] = {}
    frontier_reasons: dict[str, float] = {}
    frontier_nodes = 0
    direct_bio_nodes = 0
    adaptive_pruned_nodes = 0
    adaptive_pruned_amount_abs = 0.0
    adaptive_cutoff_potentials: list[float] = []
    frontier_score_potential_sum = 0.0

    for _node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        nodes_by_depth[depth] = nodes_by_depth.get(depth, 0) + 1
        frontier_amount = float(data.get("frontier_amount") or 0.0)
        if frontier_amount:
            frontier_nodes += 1
            frontier_score_potential_sum += float(data.get("score_potential") or 0.0)
            reasons = data.get("frontier_reasons") or {}
            for reason, amount in reasons.items():
                key = f"frontier_reason_{reason}_amount"
                frontier_reasons[key] = frontier_reasons.get(key, 0.0) + float(amount)
        if data.get("direct_bio_amount"):
            direct_bio_nodes += 1
        if data.get("adaptive_cutoff_reason") == "adaptive_relative_score_cutoff":
            adaptive_pruned_nodes += 1
            adaptive_pruned_amount_abs += abs(frontier_amount)
            adaptive_cutoff_potentials.append(
                float(data.get("adaptive_cutoff_potential") or 0.0)
            )

    out: dict[str, Any] = {
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "graph_max_depth": int(max(nodes_by_depth, default=0)),
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
        "adaptive_pruned_nodes": int(adaptive_pruned_nodes),
        "adaptive_pruned_amount_abs": float(adaptive_pruned_amount_abs),
        "frontier_score_potential_sum": float(frontier_score_potential_sum),
        "adaptive_cutoff_potential_sum": float(sum(adaptive_cutoff_potentials)),
        "adaptive_cutoff_potential_max": (
            float(max(adaptive_cutoff_potentials))
            if adaptive_cutoff_potentials
            else 0.0
        ),
    }
    for depth, count in sorted(nodes_by_depth.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    out.update(frontier_reasons)
    return out


def _safe_slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return slug.replace(".", "p").replace("-", "m")


def _reason_label(reason: str) -> str:
    labels = {
        "adaptive_relative_score_cutoff": "Adaptive cutoff",
        "max_depth": "Max depth",
        "min_amount": "Minimum amount",
        "leaf": "Leaf",
    }
    return labels.get(str(reason), str(reason).replace("_", " ").title())


def _short_reason_label(reason: str) -> str:
    labels = {
        "adaptive_relative_score_cutoff": "Cutoff",
        "max_depth": "Max depth",
        "min_amount": "Min. amount",
        "leaf": "Leaf",
    }
    return labels.get(str(reason), str(reason).replace("_", " ").title())


def _reason_color(reason: str, alpha: float = 0.55) -> str:
    colors = {
        "adaptive_relative_score_cutoff": (228, 87, 86),
        "max_depth": (90, 90, 90),
        "min_amount": (139, 102, 190),
        "leaf": (89, 161, 79),
    }
    r, g, b = colors.get(str(reason), (242, 142, 43))
    return f"rgba({r},{g},{b},{alpha})"


def _routing_depth_summary_rows(
    trails: Any,
    *,
    case: str,
    relative_cutoff: float,
    effective_cutoff: float | None,
) -> list[dict[str, Any]]:
    graph = getattr(trails, "graph", None)
    if graph is None:
        return []

    nodes_by_depth: dict[int, int] = {}
    score_potential_by_depth: dict[int, float] = {}
    outgoing_edges_by_depth: dict[int, int] = {}
    outgoing_abs_amount_by_depth: dict[int, float] = {}
    frontier_counts: dict[tuple[int, str], int] = {}
    frontier_abs_amounts: dict[tuple[int, str], float] = {}
    frontier_potentials: dict[tuple[int, str], float] = {}

    for _node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        nodes_by_depth[depth] = nodes_by_depth.get(depth, 0) + 1
        score_potential_by_depth[depth] = score_potential_by_depth.get(
            depth, 0.0
        ) + float(data.get("score_potential") or 0.0)

        reasons = data.get("frontier_reasons") or {}
        for reason, amount in reasons.items():
            key = (depth, str(reason))
            frontier_counts[key] = frontier_counts.get(key, 0) + 1
            frontier_abs_amounts[key] = frontier_abs_amounts.get(key, 0.0) + abs(
                float(amount)
            )
            frontier_potentials[key] = frontier_potentials.get(key, 0.0) + float(
                data.get("score_potential") or 0.0
            )

    for source, _target, data in graph.edges(data=True):
        source_depth = int(graph.nodes[source].get("depth", 0))
        outgoing_edges_by_depth[source_depth] = (
            outgoing_edges_by_depth.get(source_depth, 0) + 1
        )
        outgoing_abs_amount_by_depth[source_depth] = outgoing_abs_amount_by_depth.get(
            source_depth, 0.0
        ) + abs(float(data.get("amount") or 0.0))

    rows: list[dict[str, Any]] = []
    reasons = sorted({reason for _depth, reason in frontier_counts})
    for depth in sorted(nodes_by_depth):
        row: dict[str, Any] = {
            "case": case,
            "adaptive_relative_score_cutoff": float(relative_cutoff),
            "adaptive_effective_score_cutoff": effective_cutoff,
            "depth": int(depth),
            "nodes": int(nodes_by_depth.get(depth, 0)),
            "score_potential_sum": float(score_potential_by_depth.get(depth, 0.0)),
            "outgoing_edge_count": int(outgoing_edges_by_depth.get(depth, 0)),
            "outgoing_abs_amount_sum": float(
                outgoing_abs_amount_by_depth.get(depth, 0.0)
            ),
            "frontier_nodes_total": int(
                sum(
                    count
                    for (frontier_depth, _reason), count in frontier_counts.items()
                    if frontier_depth == depth
                )
            ),
        }
        for reason in reasons:
            key = (depth, reason)
            safe_reason = str(reason).replace(" ", "_")
            row[f"frontier_nodes_{safe_reason}"] = int(frontier_counts.get(key, 0))
            row[f"frontier_abs_amount_{safe_reason}"] = float(
                frontier_abs_amounts.get(key, 0.0)
            )
            row[f"frontier_score_potential_{safe_reason}"] = float(
                frontier_potentials.get(key, 0.0)
            )
        rows.append(row)
    return rows


def _write_routing_depth_sankey(
    trails: Any,
    *,
    case: str,
    relative_cutoff: float,
    static_score: float,
    score: float,
    relative_deviation: float,
    graph_stats: dict[str, Any],
    routing_params: dict[str, Any],
    output_dir: Path,
    write_png: bool,
) -> tuple[str, str]:
    graph = getattr(trails, "graph", None)
    if graph is None:
        return "", ""

    import plotly.graph_objects as go

    output_dir.mkdir(parents=True, exist_ok=True)

    nodes_by_depth: dict[int, int] = {}
    edge_counts: dict[tuple[int, int], int] = {}
    edge_amounts: dict[tuple[int, int], float] = {}
    frontier_counts: dict[tuple[int, str], int] = {}
    frontier_amounts: dict[tuple[int, str], float] = {}
    frontier_potentials: dict[tuple[int, str], float] = {}

    for _node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        nodes_by_depth[depth] = nodes_by_depth.get(depth, 0) + 1
        for reason, amount in (data.get("frontier_reasons") or {}).items():
            key = (depth, str(reason))
            frontier_counts[key] = frontier_counts.get(key, 0) + 1
            frontier_amounts[key] = frontier_amounts.get(key, 0.0) + abs(float(amount))
            frontier_potentials[key] = frontier_potentials.get(key, 0.0) + float(
                data.get("score_potential") or 0.0
            )

    for source, target, data in graph.edges(data=True):
        source_depth = int(graph.nodes[source].get("depth", 0))
        target_depth = int(graph.nodes[target].get("depth", source_depth + 1))
        key = (source_depth, target_depth)
        edge_counts[key] = edge_counts.get(key, 0) + 1
        edge_amounts[key] = edge_amounts.get(key, 0.0) + abs(
            float(data.get("amount") or 0.0)
        )

    labels: list[str] = []
    colors: list[str] = []
    node_index: dict[tuple[Any, ...], int] = {}

    def add_node(key: tuple[Any, ...], label: str, color: str) -> int:
        if key not in node_index:
            node_index[key] = len(labels)
            labels.append(label)
            colors.append(color)
        return node_index[key]

    for depth in sorted(nodes_by_depth):
        add_node(
            ("depth", depth),
            f"D{depth}<br>{nodes_by_depth[depth]:,} nodes",
            "rgba(76,120,168,0.82)",
        )

    sources: list[int] = []
    targets: list[int] = []
    values: list[float] = []
    link_colors: list[str] = []
    customdata: list[str] = []

    for (source_depth, target_depth), count in sorted(edge_counts.items()):
        sources.append(add_node(("depth", source_depth), f"D{source_depth}", ""))
        targets.append(add_node(("depth", target_depth), f"D{target_depth}", ""))
        values.append(float(count))
        link_colors.append("rgba(76,120,168,0.32)")
        customdata.append(
            "Continued routing"
            f"<br>From depth {source_depth} to depth {target_depth}"
            f"<br>Graph edges: {count:,}"
            f"<br>Sum abs(edge amounts): {edge_amounts[(source_depth, target_depth)]:.4g}"
        )

    for (depth, reason), count in sorted(frontier_counts.items()):
        stop_label = f"{_short_reason_label(reason)} d{depth}<br>{count:,} nodes"
        stop_node = add_node(
            ("stop", depth, reason),
            stop_label,
            _reason_color(reason, alpha=0.82),
        )
        sources.append(add_node(("depth", depth), f"D{depth}", ""))
        targets.append(stop_node)
        values.append(float(count))
        link_colors.append(_reason_color(reason, alpha=0.42))
        customdata.append(
            f"Stopped: {_reason_label(reason)}"
            f"<br>Depth: {depth}"
            f"<br>Frontier nodes: {count:,}"
            f"<br>Sum abs(frontier amounts): {frontier_amounts[(depth, reason)]:.4g}"
            f"<br>Sum score potential: {frontier_potentials[(depth, reason)]:.4g}"
        )

    effective_cutoff = routing_params.get("adaptive_effective_score_cutoff")
    root_potential = routing_params.get("adaptive_root_score_potential")
    score_text = "not recalculated"
    if np.isfinite(score):
        score_text = (
            f"{score:.5g}; deviation {relative_deviation:.3%} vs static "
            f"{static_score:.5g}"
        )

    title = (
        "Adaptive routing depth Sankey"
        f"<br><sup>{case}; relative cutoff={relative_cutoff:g}; "
        f"effective cutoff={float(effective_cutoff):.4g}; "
        f"FU potential={float(root_potential):.4g}; "
        f"score={score_text}; "
        f"nodes={int(graph_stats.get('graph_nodes', 0)):,}; "
        f"edges={int(graph_stats.get('graph_edges', 0)):,}. "
        "Link widths are counts: continued graph edges and stopped frontier nodes."
        "</sup>"
    )

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=18,
                    thickness=18,
                    line=dict(color="rgba(40,40,40,0.35)", width=0.5),
                    label=labels,
                    color=colors,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                    customdata=customdata,
                    hovertemplate="%{customdata}<extra></extra>",
                ),
            )
        ]
    )
    fig.update_layout(
        title=title,
        font=dict(size=12),
        width=1600,
        height=1100,
        margin=dict(l=30, r=30, t=110, b=30),
    )

    filename = f"{_safe_slug(case)}_routing_depth_sankey"
    html_path = output_dir / f"{filename}.html"
    fig.write_html(str(html_path), include_plotlyjs="cdn", full_html=True)

    png_path = ""
    if write_png:
        candidate = output_dir / f"{filename}.png"
        try:
            png_fig = go.Figure(fig)
            png_fig.update_layout(
                title=f"Adaptive routing depth Sankey: {case}",
                margin=dict(l=30, r=30, t=70, b=30),
            )
            png_fig.write_image(str(candidate), width=1600, height=1100, scale=2)
            png_path = str(candidate)
        except Exception as exc:
            print(f"  could not write Sankey PNG for {case}: {exc}", flush=True)

    return str(html_path), png_path


def _wrap_plotly_label(
    label: str,
    *,
    line_chars: int = 18,
    max_lines: int = 3,
) -> str:
    """Wrap a Plotly HTML label without forcing a wide Sankey layout."""
    lines = textwrap.wrap(
        re.sub(r"\s+", " ", label).strip(),
        width=max(8, int(line_chars)),
        break_long_words=False,
        break_on_hyphens=False,
    )
    if not lines:
        return label
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(". ") + "..."
    return "<br>".join(lines)


def _short_activity_label(
    data: dict[str, Any],
    max_chars: int = 48,
    *,
    wrap: bool = False,
    line_chars: int = 18,
    max_lines: int = 3,
) -> str:
    name = str(data.get("name") or f"activity {data.get('act_idx', '')}").strip()
    ref = str(data.get("reference_product") or "").strip()
    location = str(data.get("location") or "").strip()
    label = name
    if ref:
        label = f"{label} | {ref}"
    if location:
        label = f"{label} | {location}"
    if len(label) > max_chars:
        label = label[: max_chars - 1].rstrip() + "..."
    if wrap:
        return _wrap_plotly_label(
            label,
            line_chars=line_chars,
            max_lines=max_lines,
        )
    return label


def _year_color(year: int, min_year: int, max_year: int, alpha: float = 0.82) -> str:
    if max_year <= min_year:
        t = 0.5
    else:
        t = (int(year) - int(min_year)) / float(int(max_year) - int(min_year))
    # Early years are blue; later years move toward red.
    r = int(round(76 + t * (228 - 76)))
    g = int(round(120 + t * (87 - 120)))
    b = int(round(168 + t * (86 - 168)))
    return f"rgba({r},{g},{b},{alpha})"


def _display_year_y(year: int, min_year: int, max_year: int) -> float:
    if max_year <= min_year:
        return 0.5
    return (int(year) - int(min_year)) / float(int(max_year) - int(min_year))


def _write_activity_year_score_sankey(
    trails: Any,
    *,
    case: str,
    relative_cutoff: float,
    static_score: float,
    score: float,
    relative_deviation: float,
    graph_stats: dict[str, Any],
    routing_params: dict[str, Any],
    output_dir: Path,
    write_png: bool,
    max_nodes: int,
    max_links: int,
    min_link_relative_score: float,
    min_label_relative_score: float,
) -> tuple[str, str]:
    """Write a score-potential-weighted activity/year Sankey for one run.

    The full routing graph can contain hundreds of thousands of nodes, so the
    figure groups low-score activities by year/depth and filters low-score
    links. To avoid mistaking filtering for a routing cap, one weighted path to
    the deepest routed depth is kept and labeled explicitly.
    """
    graph = getattr(trails, "graph", None)
    if graph is None:
        return "", ""

    import plotly.graph_objects as go

    output_dir.mkdir(parents=True, exist_ok=True)

    node_weight = {
        node: abs(float(data.get("score_potential") or 0.0))
        for node, data in graph.nodes(data=True)
    }
    if not node_weight:
        return "", ""

    years = [int(data.get("year", 0)) for _node, data in graph.nodes(data=True)]
    depths = [int(data.get("depth", 0)) for _node, data in graph.nodes(data=True)]
    min_year = min(years)
    max_year = max(years)
    max_depth_seen = max(depths) if depths else 0
    score_reference = abs(float(score)) if np.isfinite(score) and score else 0.0
    if score_reference <= 0:
        score_reference = abs(
            float(routing_params.get("adaptive_root_score_potential") or 0.0)
        )
    if score_reference <= 0:
        score_reference = max(node_weight.values(), default=1.0)

    root_nodes = {
        node for node, data in graph.nodes(data=True) if int(data.get("depth", 0)) == 0
    }

    node_depth = {
        node: int(data.get("depth", 0)) for node, data in graph.nodes(data=True)
    }

    def deepest_weighted_path() -> list[Any]:
        """Return one explicit root-to-deepest path for visual orientation."""
        if max_depth_seen <= 0:
            return []

        best_score: dict[Any, float] = {}
        best_predecessor: dict[Any, Any | None] = {}
        for node in sorted(
            graph.nodes,
            key=lambda graph_node: (
                node_depth.get(graph_node, 0),
                -node_weight.get(graph_node, 0.0),
            ),
        ):
            depth = node_depth.get(node, 0)
            predecessors = [
                predecessor
                for predecessor in graph.predecessors(node)
                if node_depth.get(predecessor, -1) < depth and predecessor in best_score
            ]
            if predecessors:
                predecessor = max(
                    predecessors,
                    key=lambda candidate: (
                        best_score[candidate],
                        node_weight.get(candidate, 0.0),
                    ),
                )
                best_score[node] = best_score[predecessor] + node_weight.get(node, 0.0)
                best_predecessor[node] = predecessor
            elif depth == 0:
                best_score[node] = node_weight.get(node, 0.0)
                best_predecessor[node] = None

        deepest_reachable_depth = max(
            (node_depth[node] for node in best_score),
            default=0,
        )
        deepest_nodes = [
            node
            for node, depth in node_depth.items()
            if depth == deepest_reachable_depth and node in best_score
        ]
        if not deepest_nodes:
            return []

        end_node = max(
            deepest_nodes,
            key=lambda node: (best_score[node], node_weight.get(node, 0.0)),
        )
        path: list[Any] = []
        seen: set[Any] = set()
        while end_node is not None and end_node not in seen:
            path.append(end_node)
            seen.add(end_node)
            end_node = best_predecessor.get(end_node)
        path.reverse()
        return path

    depth_spine_nodes = deepest_weighted_path()
    ranked_nodes = sorted(node_weight, key=lambda node: node_weight[node], reverse=True)
    selected_nodes = set(ranked_nodes[: max(1, int(max_nodes))])
    selected_nodes.update(root_nodes)
    selected_nodes.update(depth_spine_nodes)

    incoming_abs_amount: dict[Any, float] = {node: 0.0 for node in graph.nodes}
    for _source, target, data in graph.edges(data=True):
        incoming_abs_amount[target] = incoming_abs_amount.get(target, 0.0) + abs(
            float(data.get("amount") or 0.0)
        )

    display_meta: dict[tuple[Any, ...], dict[str, Any]] = {}

    def depth_x(depth: int, *, stop: bool = False) -> float:
        if max_depth_seen <= 0:
            return 0.05 if not stop else 0.95
        offset = 0.6 if stop else 0.0
        return min(0.98, max(0.02, (float(depth) + offset) / float(max_depth_seen + 1)))

    def display_node(graph_node: Any) -> tuple[Any, ...]:
        data = graph.nodes[graph_node]
        depth = int(data.get("depth", 0))
        year = int(data.get("year", 0))
        if graph_node in selected_nodes:
            key = ("activity", graph_node)
            if key not in display_meta:
                potential = node_weight.get(graph_node, 0.0)
                label = _short_activity_label(
                    data,
                    max_chars=44,
                    wrap=True,
                    line_chars=18,
                    max_lines=3,
                )
                display_meta[key] = {
                    "label": (
                        f"{label}" f"<br>{year}, d{depth}" f"<br>{potential:.3g}"
                    ),
                    "hover": (
                        f"<b>{_short_activity_label(data, max_chars=90)}</b>"
                        f"<br>Year: {year}"
                        f"<br>Depth: {depth}"
                        f"<br>Activity index: {data.get('act_idx')}"
                        f"<br>Score potential: {potential:.6g}"
                        f"<br>Share of temporal score: "
                        f"{potential / score_reference:.3%}"
                    ),
                    "year": year,
                    "depth": depth,
                    "x": depth_x(depth),
                    "color": _year_color(year, min_year, max_year),
                }
            return key

        key = ("other", depth, year)
        if key not in display_meta:
            display_meta[key] = {
                "label": f"Other<br>{year}, d{depth}",
                "hover": (
                    f"<b>Other routed activities</b>"
                    f"<br>Year: {year}"
                    f"<br>Depth: {depth}"
                    "<br>Low-score nodes grouped for readability"
                ),
                "year": year,
                "depth": depth,
                "x": depth_x(depth),
                "color": _year_color(year, min_year, max_year, alpha=0.35),
            }
        return key

    def stop_node(depth: int, year: int, reason: str) -> tuple[Any, ...]:
        key = ("stop", str(reason), int(depth), int(year))
        if key not in display_meta:
            display_meta[key] = {
                "label": f"{_short_reason_label(reason)}<br>{year}, d{depth}",
                "hover": (
                    f"<b>{_reason_label(reason)}</b>"
                    f"<br>Year: {year}"
                    f"<br>Depth: {depth}"
                    "<br>Frontier demand solved by matrix/background"
                ),
                "year": int(year),
                "depth": int(depth),
                "x": depth_x(int(depth), stop=True),
                "color": _reason_color(str(reason), alpha=0.82),
            }
        return key

    link_data: dict[tuple[tuple[Any, ...], tuple[Any, ...]], dict[str, Any]] = {}

    def add_link(
        source_key: tuple[Any, ...],
        target_key: tuple[Any, ...],
        *,
        value: float,
        color: str,
        kind: str,
        amount: float = 0.0,
    ) -> None:
        if source_key == target_key or value <= 0:
            return
        key = (source_key, target_key)
        row = link_data.setdefault(
            key,
            {
                "value": 0.0,
                "amount_abs": 0.0,
                "count": 0,
                "kind": kind,
                "color": color,
            },
        )
        row["value"] = float(row["value"]) + float(value)
        row["amount_abs"] = float(row["amount_abs"]) + abs(float(amount))
        row["count"] = int(row["count"]) + 1

    for source, target, data in graph.edges(data=True):
        target_weight = node_weight.get(target, 0.0)
        amount_abs = abs(float(data.get("amount") or 0.0))
        incoming = incoming_abs_amount.get(target, 0.0)
        if incoming > 0:
            value = target_weight * amount_abs / incoming
        else:
            value = target_weight
        if value <= 0:
            continue
        source_key = display_node(source)
        target_key = display_node(target)
        target_year = int(graph.nodes[target].get("year", 0))
        add_link(
            source_key,
            target_key,
            value=value,
            color=_year_color(target_year, min_year, max_year, alpha=0.34),
            kind="routed demand",
            amount=amount_abs,
        )

    for node, data in graph.nodes(data=True):
        reasons = data.get("frontier_reasons") or {}
        if not reasons:
            continue
        frontier_total = sum(abs(float(amount)) for amount in reasons.values())
        if frontier_total <= 0:
            frontier_total = 1.0
        source_key = display_node(node)
        depth = int(data.get("depth", 0))
        year = int(data.get("year", 0))
        potential = node_weight.get(node, 0.0)
        for reason, amount in reasons.items():
            amount_abs = abs(float(amount))
            value = potential * amount_abs / frontier_total
            add_link(
                source_key,
                stop_node(depth, year, str(reason)),
                value=value,
                color=_reason_color(str(reason), alpha=0.42),
                kind=f"frontier: {_reason_label(str(reason))}",
                amount=amount_abs,
            )

    min_link_value = float(min_link_relative_score) * float(score_reference)
    all_links = [
        (source_key, target_key, row)
        for (source_key, target_key), row in link_data.items()
    ]
    depth_spine_link_keys = set()
    for source, target in zip(depth_spine_nodes, depth_spine_nodes[1:]):
        source_key = display_node(source)
        target_key = display_node(target)
        if (source_key, target_key) in link_data:
            depth_spine_link_keys.add((source_key, target_key))

    filtered_links = [
        (source_key, target_key, row)
        for source_key, target_key, row in all_links
        if float(row["value"]) >= min_link_value
        or (source_key, target_key) in depth_spine_link_keys
    ]
    if len(filtered_links) > int(max_links):
        spine_links = [
            item
            for item in filtered_links
            if (item[0], item[1]) in depth_spine_link_keys
        ]
        non_spine_links = [
            item
            for item in filtered_links
            if (item[0], item[1]) not in depth_spine_link_keys
        ]
        non_spine_links = sorted(
            non_spine_links,
            key=lambda item: float(item[2]["value"]),
            reverse=True,
        )
        filtered_links = (
            spine_links + non_spine_links[: max(0, int(max_links) - len(spine_links))]
        )

    filtered_keys = {
        (source_key, target_key) for source_key, target_key, _row in filtered_links
    }
    for depth in range(int(max_depth_seen) + 1):
        if any(
            int(display_meta[source_key]["depth"]) == depth
            or int(display_meta[target_key]["depth"]) == depth
            for source_key, target_key, _row in filtered_links
        ):
            continue
        candidates = [
            (source_key, target_key, row)
            for source_key, target_key, row in all_links
            if int(display_meta[source_key]["depth"]) == depth
            or int(display_meta[target_key]["depth"]) == depth
        ]
        if not candidates:
            continue
        source_key, target_key, row = max(
            candidates,
            key=lambda item: float(item[2]["value"]),
        )
        if (source_key, target_key) not in filtered_keys:
            filtered_links.append((source_key, target_key, row))
            filtered_keys.add((source_key, target_key))

    used_keys = {
        key
        for source_key, target_key, _row in filtered_links
        for key in (source_key, target_key)
    }
    ordered_keys = sorted(
        used_keys,
        key=lambda key: (
            int(display_meta[key]["depth"]),
            int(display_meta[key]["year"]),
            str(display_meta[key]["label"]),
        ),
    )
    key_to_index = {key: pos for pos, key in enumerate(ordered_keys)}
    if ordered_keys:
        plot_years = [int(display_meta[key]["year"]) for key in ordered_keys]
        plot_min_year = min(plot_years)
        plot_max_year = max(plot_years)
    else:
        plot_min_year = min_year
        plot_max_year = max_year

    groups: dict[tuple[int, int], list[tuple[Any, ...]]] = {}
    for key in ordered_keys:
        depth = int(display_meta[key]["depth"])
        year = int(display_meta[key]["year"])
        groups.setdefault((depth, year), []).append(key)

    node_values = {key: 0.0 for key in ordered_keys}
    for source_key, target_key, row in filtered_links:
        value = float(row["value"])
        node_values[source_key] = node_values.get(source_key, 0.0) + value
        node_values[target_key] = node_values.get(target_key, 0.0) + value

    node_x: list[float] = []
    node_y_by_key: dict[tuple[Any, ...], float] = {}
    year_count = max(1, plot_max_year - plot_min_year + 1)
    year_band = min(0.06, 0.85 / float(year_count))
    for (_depth, year), keys in groups.items():
        keys_sorted = sorted(
            keys, key=lambda key: node_values.get(key, 0.0), reverse=True
        )
        n_keys = len(keys_sorted)
        for pos, key in enumerate(keys_sorted):
            base_y = _display_year_y(year, plot_min_year, plot_max_year)
            if n_keys <= 1:
                offset = 0.0
            else:
                offset = ((pos / float(n_keys - 1)) - 0.5) * year_band
            node_y_by_key[key] = min(0.98, max(0.02, base_y + offset))

    min_label_value = float(min_label_relative_score) * float(score_reference)
    depth_spine_keys = {display_node(node) for node in depth_spine_nodes}
    labels = []
    for key in ordered_keys:
        value = float(node_values.get(key, 0.0))
        meta = display_meta[key]
        if (
            value >= min_label_value
            or key[0] in {"stop", "other"}
            or key in depth_spine_keys
        ):
            labels.append(str(meta["label"]))
        elif key[0] == "activity":
            labels.append("")
        else:
            labels.append("")
    colors = [str(display_meta[key]["color"]) for key in ordered_keys]
    node_x = [float(display_meta[key]["x"]) for key in ordered_keys]
    node_y = [float(node_y_by_key.get(key, 0.5)) for key in ordered_keys]
    custom_nodes = [
        str(display_meta[key]["hover"])
        + f"<br>Displayed Sankey value: {node_values.get(key, 0.0):.6g}"
        + (
            "<br>Retained as deepest-depth display path"
            if key in depth_spine_keys
            else ""
        )
        for key in ordered_keys
    ]

    sources: list[int] = []
    targets: list[int] = []
    values: list[float] = []
    link_colors: list[str] = []
    custom_links: list[str] = []
    for source_key, target_key, row in sorted(
        filtered_links,
        key=lambda item: float(item[2]["value"]),
        reverse=True,
    ):
        value = float(row["value"])
        sources.append(key_to_index[source_key])
        targets.append(key_to_index[target_key])
        values.append(value)
        link_colors.append(str(row["color"]))
        custom_links.append(
            f"<b>{row['kind']}</b>"
            f"<br>Score potential: {value:.6g}"
            f"<br>Share of temporal score: {value / score_reference:.3%}"
            f"<br>Aggregated graph links/nodes: {int(row['count']):,}"
            f"<br>Sum abs(amount): {float(row['amount_abs']):.6g}"
        )

    effective_cutoff = routing_params.get("adaptive_effective_score_cutoff")
    root_potential = routing_params.get("adaptive_root_score_potential")
    depth_spine_depth = max(
        (node_depth.get(node, 0) for node in depth_spine_nodes),
        default=0,
    )
    title = (
        "Adaptive activity-year Sankey"
        f"<br><span style='font-size:11px'>case={case}</span>"
        f"<br><span style='font-size:11px'>relative cutoff={relative_cutoff:g}; "
        f"effective cutoff={float(effective_cutoff):.4g}; "
        f"FU potential={float(root_potential):.4g}</span>"
        f"<br><span style='font-size:11px'>temporal score={score:.5g}; "
        f"static score={static_score:.5g}; "
        f"deviation={relative_deviation:.3%}</span>"
        f"<br><span style='font-size:11px'>shown: {len(ordered_keys):,}/"
        f"{graph.number_of_nodes():,} nodes and {len(filtered_links):,}/"
        f"{len(link_data):,} links; graph max depth={max_depth_seen}; "
        f"display path reaches d{depth_spine_depth}</span>"
        "<br><span style='font-size:11px'>Node and link sizes use static LCIA "
        "score potential. Vertical position is year.</span>"
    )

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="fixed",
                node=dict(
                    pad=14,
                    thickness=16,
                    line=dict(color="rgba(45,45,45,0.35)", width=0.4),
                    label=labels,
                    color=colors,
                    x=node_x,
                    y=node_y,
                    customdata=custom_nodes,
                    hovertemplate="%{customdata}<extra></extra>",
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                    customdata=custom_links,
                    hovertemplate="%{customdata}<extra></extra>",
                ),
            )
        ]
    )
    fig.update_layout(
        title=dict(text=title, x=0.01, xanchor="left"),
        font=dict(size=10),
        width=1900,
        height=1300,
        margin=dict(l=35, r=35, t=150, b=35),
    )

    filename = f"{_safe_slug(case)}_activity_year_score_sankey"
    html_path = output_dir / f"{filename}.html"
    fig.write_html(str(html_path), include_plotlyjs="cdn", full_html=True)

    png_path = ""
    if write_png:
        candidate = output_dir / f"{filename}.png"
        try:
            png_fig = go.Figure(fig)
            png_fig.update_layout(
                title=f"Adaptive activity-year score Sankey: {case}",
                margin=dict(l=35, r=35, t=75, b=35),
            )
            png_fig.write_image(str(candidate), width=1900, height=1300, scale=2)
            png_path = str(candidate)
        except Exception as exc:
            print(
                f"  could not write activity Sankey PNG for {case}: {exc}", flush=True
            )

    return str(html_path), png_path


def _write_temporal_score_plot(
    trails: Any,
    *,
    case: str,
    method: str,
    unit: str,
    static_score: float,
    output_dir: Path,
    write_png: bool,
    reference_year: int,
) -> tuple[str, str]:
    if getattr(trails, "scores", None) is None:
        return "", ""

    from trails.plotting import plot_temporal_scores

    output_dir.mkdir(parents=True, exist_ok=True)
    label = unit or "Impact score"
    fig = plot_temporal_scores(
        trails,
        title=f"Temporal score by year: {case}",
        method_label=label,
        method=method,
        cumulative=False,
        stacked=True,
        legend_top_n=8,
        show_flow_contributions=False,
        width=1500,
        height=850,
        year_tick=5,
        reference_year=int(reference_year),
        show_cumulative_axis=True,
        cumulative_axis_label=f"Cumulative {label}",
        show_cumulative_in_legend=True,
        static_score=float(static_score),
        static_score_label="Static score",
    )

    filename = f"{_safe_slug(case)}_temporal_score"
    html_path = output_dir / f"{filename}.html"
    fig.write_html(str(html_path), include_plotlyjs="cdn", full_html=True)

    png_path = ""
    if write_png:
        candidate = output_dir / f"{filename}.png"
        try:
            fig.write_image(str(candidate), width=1500, height=850, scale=2)
            png_path = str(candidate)
        except Exception as exc:
            print(f"  could not write temporal score PNG for {case}: {exc}", flush=True)
    return str(html_path), png_path


def _temporal_total_score(trails: Any, method: str) -> float:
    if getattr(trails, "scores", None) is None:
        raise RuntimeError("trails.scores is missing after temporal LCA.")
    data = trails.scores
    if "method" in data.dims:
        methods = [str(value) for value in data.coords["method"].values.tolist()]
        data = data.isel(method=methods.index(str(method)), drop=True)
    if data.dims:
        data = data.sum(dim=list(data.dims))
    values = data.data
    if hasattr(values, "todense"):
        return float(np.asarray(values.todense(), dtype=float).sum())
    return float(np.asarray(data.values, dtype=float).sum())


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


def _comparison_rows(
    adaptive_rows: list[dict[str, Any]],
    fixed_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    if not fixed_rows:
        return []
    depth7 = next((row for row in fixed_rows if row["_depth"] == 7), fixed_rows[-1])
    out = []
    for row in adaptive_rows:
        adaptive_nodes = int(row["graph_nodes"])
        nearest = min(
            fixed_rows,
            key=lambda fixed: abs(int(fixed["_graph_nodes"]) - adaptive_nodes),
        )
        score = float(row["score"])
        out.append(
            {
                "adaptive_case": row["case"],
                "adaptive_relative_cutoff": row["adaptive_relative_score_cutoff"],
                "adaptive_effective_cutoff": row["adaptive_effective_score_cutoff"],
                "adaptive_score": score,
                "adaptive_relative_deviation_from_static": row[
                    "relative_deviation_from_static"
                ],
                "adaptive_graph_nodes": adaptive_nodes,
                "adaptive_graph_edges": row["graph_edges"],
                "adaptive_routing_seconds": row["routing_seconds"],
                "adaptive_lca_seconds": row["temporal_lca_seconds"],
                "nearest_depth_by_nodes": nearest["_depth"],
                "nearest_depth_score": nearest["_score"],
                "nearest_depth_graph_nodes": nearest["_graph_nodes"],
                "nearest_depth_routing_seconds": nearest["_routing_seconds"],
                "score_minus_nearest_depth": score - float(nearest["_score"]),
                "routing_speedup_vs_nearest_depth": (
                    float(nearest["_routing_seconds"]) / float(row["routing_seconds"])
                    if float(row["routing_seconds"]) > 0
                    else float("nan")
                ),
                "score_minus_depth7": score - float(depth7["_score"]),
                "routing_speedup_vs_depth7": (
                    float(depth7["_routing_seconds"]) / float(row["routing_seconds"])
                    if float(row["routing_seconds"]) > 0
                    else float("nan")
                ),
                "node_ratio_vs_depth7": adaptive_nodes / float(depth7["_graph_nodes"]),
            }
        )
    return out


def run_adaptive(args: argparse.Namespace) -> int:
    datapackage = Path(args.datapackage).expanduser().resolve()
    inventory_paths = [Path(path).expanduser().resolve() for path in args.inventories]
    depth_runner._validate_paths(datapackage, inventory_paths)

    lcia_json = None if args.lcia_json is None else Path(args.lcia_json).expanduser()
    if lcia_json is not None:
        if not lcia_json.exists():
            raise FileNotFoundError(f"LCIA JSON not found: {lcia_json}")
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json.resolve())

    depth_rows_raw = _read_csv(Path(args.depth_sweep_csv).expanduser().resolve())
    fixed_rows = _depth_rows(depth_rows_raw)
    static_score = _static_score(depth_rows_raw)
    if static_score is None:
        raise RuntimeError(
            "No static row found in depth-sweep CSV; run this script with the depth-sweep subcommand first."
        )

    output_csv = Path(args.output_csv).expanduser().resolve()
    comparison_csv = Path(args.comparison_csv).expanduser().resolve()
    sankey_dir = Path(args.sankey_dir).expanduser().resolve()
    depth_sankey_dir = Path(args.depth_sankey_dir).expanduser().resolve()
    score_plot_dir = Path(args.score_plot_dir).expanduser().resolve()
    sankey_summary_csv = Path(args.sankey_summary_csv).expanduser().resolve()
    existing_rows = _read_csv(output_csv) if args.append else []
    existing_by_case = {str(row.get("case", "")): dict(row) for row in existing_rows}
    existing_cases = {str(row.get("case", "")) for row in existing_rows}
    adaptive_rows: list[dict[str, Any]] = [dict(row) for row in existing_rows]
    sankey_summary_rows: list[dict[str, Any]] = (
        [dict(row) for row in _read_csv(sankey_summary_csv)]
        if args.append and sankey_summary_csv.exists()
        else []
    )

    print(f"Datapackage: {datapackage}", flush=True)
    print(f"Depth sweep CSV: {args.depth_sweep_csv}", flush=True)
    print(f"Method: {args.method}", flush=True)
    print(f"Static score from depth sweep: {static_score:.12g}", flush=True)

    load_start = time.perf_counter()
    trails = _load_trails(
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

    activity_maps = _match_activity_indices(
        trails, [depth_runner.DEFAULT_ACTIVITY]
    )
    if depth_runner.DEFAULT_ACTIVITY not in activity_maps:
        raise RuntimeError(
            f"Could not match DACCS activity: {depth_runner.DEFAULT_ACTIVITY}"
        )
    activity_index = int(activity_maps[depth_runner.DEFAULT_ACTIVITY])
    activity_label = _activity_label(
        trails,
        activity_index,
        int(args.reference_year),
    )
    unit = _method_unit(args.method, str(args.ei_version))
    print(f"Activity: {activity_label} (idx={activity_index})", flush=True)

    for relative_cutoff in args.relative_cutoffs:
        relative_cutoff = float(relative_cutoff)
        routing_max_depth = None if bool(args.no_max_depth) else int(args.max_depth)
        depth_label = (
            "nomaxdepth"
            if routing_max_depth is None
            else f"maxdepth_{routing_max_depth}"
        )
        case = f"adaptive_rel_{relative_cutoff:g}_{depth_label}"
        if args.append and not args.force and case in existing_cases:
            print(f"Skipping existing case {case}", flush=True)
            continue
        if args.append and args.force and case in existing_cases:
            adaptive_rows = [
                row for row in adaptive_rows if str(row.get("case", "")) != case
            ]
            sankey_summary_rows = [
                row for row in sankey_summary_rows if str(row.get("case", "")) != case
            ]
            existing_cases.discard(case)

        print(f"Running {case}", flush=True)
        step_start = time.perf_counter()
        routing_start = time.perf_counter()
        trails.temporal_routing(
            start_year=int(args.reference_year),
            start_act_idx=activity_index,
            amount=float(args.amount),
            max_depth=routing_max_depth,
            min_amount=float(args.routing_min_amount),
            show_progress=bool(args.show_progress),
            attribute_to_roots=True,
            adaptive_methods=[args.method],
            adaptive_relative_score_cutoff=relative_cutoff,
            adaptive_ei_version=str(args.ei_version),
            adaptive_min_depth=int(args.adaptive_min_depth),
            adaptive_use_cache=not bool(args.no_static_score_cache),
        )
        routing_seconds = time.perf_counter() - routing_start
        graph_stats = _graph_stats(trails)
        routing_params = getattr(trails, "_routing_params", {}) or {}
        print(
            "  routing "
            f"{routing_seconds:.1f}s, nodes={graph_stats['graph_nodes']:,}, "
            f"edges={graph_stats['graph_edges']:,}, "
            f"pruned={graph_stats['adaptive_pruned_nodes']:,}",
            flush=True,
        )

        previous_row = existing_by_case.get(case, {})
        if args.skip_lca:
            score = _float(previous_row.get("score"))
            deviation = _float(previous_row.get("score_deviation_from_static"))
            relative_deviation = _float(
                previous_row.get("relative_deviation_from_static")
            )
            if np.isfinite(score) and not np.isfinite(deviation):
                deviation = float(score) - float(static_score)
            if np.isfinite(deviation) and not np.isfinite(relative_deviation):
                relative_deviation = deviation / float(static_score)
            lca_seconds = _float(previous_row.get("temporal_lca_seconds"))
            actual_solver = str(previous_row.get("solver_mode") or "not_run")
            print("  skipped temporal LCA; reusing existing score metadata", flush=True)
        else:
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
            deviation = float(score) - float(static_score)
            relative_deviation = deviation / float(static_score)

        sankey_html = str(previous_row.get("sankey_html") or "")
        sankey_png = str(previous_row.get("sankey_png") or "")
        depth_sankey_html = str(previous_row.get("depth_sankey_html") or "")
        depth_sankey_png = str(previous_row.get("depth_sankey_png") or "")
        score_plot_html = str(previous_row.get("score_plot_html") or "")
        score_plot_png = str(previous_row.get("score_plot_png") or "")
        if args.write_sankey:
            sankey_html, sankey_png = _write_activity_year_score_sankey(
                trails,
                case=case,
                relative_cutoff=relative_cutoff,
                static_score=float(static_score),
                score=float(score),
                relative_deviation=float(relative_deviation),
                graph_stats=graph_stats,
                routing_params=routing_params,
                output_dir=sankey_dir,
                write_png=bool(args.write_sankey_png),
                max_nodes=int(args.sankey_max_nodes),
                max_links=int(args.sankey_max_links),
                min_link_relative_score=float(args.sankey_min_link_relative_score),
                min_label_relative_score=float(args.sankey_min_label_relative_score),
            )
            print(f"  activity sankey={sankey_html}", flush=True)

        if args.write_depth_sankey:
            depth_sankey_html, depth_sankey_png = _write_routing_depth_sankey(
                trails,
                case=case,
                relative_cutoff=relative_cutoff,
                static_score=float(static_score),
                score=float(score),
                relative_deviation=float(relative_deviation),
                graph_stats=graph_stats,
                routing_params=routing_params,
                output_dir=depth_sankey_dir,
                write_png=bool(args.write_sankey_png),
            )
            print(f"  depth sankey={depth_sankey_html}", flush=True)

        if args.write_sankey or args.write_depth_sankey:
            sankey_summary_rows.extend(
                _routing_depth_summary_rows(
                    trails,
                    case=case,
                    relative_cutoff=relative_cutoff,
                    effective_cutoff=routing_params.get(
                        "adaptive_effective_score_cutoff"
                    ),
                )
            )
            _write_csv(sankey_summary_rows, sankey_summary_csv)

        if args.write_score_plots:
            if args.skip_lca:
                print(
                    "  skipped temporal score plot because --skip-lca is set",
                    flush=True,
                )
            else:
                score_plot_html, score_plot_png = _write_temporal_score_plot(
                    trails,
                    case=case,
                    method=args.method,
                    unit=unit,
                    static_score=float(static_score),
                    output_dir=score_plot_dir,
                    write_png=bool(args.write_score_plot_png),
                    reference_year=int(args.reference_year),
                )
                print(f"  temporal score plot={score_plot_html}", flush=True)

        step_seconds = time.perf_counter() - step_start

        row = {
            "case": case,
            "mode": "adaptive",
            "activity_index": activity_index,
            "activity": activity_label,
            "reference_year": int(args.reference_year),
            "amount": float(args.amount),
            "method": args.method,
            "unit": unit,
            "max_depth": "" if routing_max_depth is None else int(routing_max_depth),
            "max_depth_disabled": bool(routing_max_depth is None),
            "routing_min_amount": float(args.routing_min_amount),
            "adaptive_min_depth": int(args.adaptive_min_depth),
            "adaptive_relative_score_cutoff": relative_cutoff,
            "adaptive_effective_score_cutoff": routing_params.get(
                "adaptive_effective_score_cutoff"
            ),
            "adaptive_root_score_potential": routing_params.get(
                "adaptive_root_score_potential"
            ),
            "score": score,
            "static_score": static_score,
            "score_deviation_from_static": deviation,
            "relative_deviation_from_static": relative_deviation,
            "load_seconds": load_seconds,
            "routing_seconds": routing_seconds,
            "temporal_lca_seconds": lca_seconds,
            "total_step_seconds": step_seconds,
            "solver_mode": actual_solver,
            "sankey_html": sankey_html,
            "sankey_png": sankey_png,
            "depth_sankey_html": depth_sankey_html,
            "depth_sankey_png": depth_sankey_png,
            "score_plot_html": score_plot_html,
            "score_plot_png": score_plot_png,
            **graph_stats,
        }
        adaptive_rows.append(row)
        _write_csv(adaptive_rows, output_csv)
        comparison = _comparison_rows(adaptive_rows, fixed_rows)
        _write_csv(comparison, comparison_csv)

        print(
            f"  score={score:.12g}, rel_dev={relative_deviation:.6g}, "
            f"lca={lca_seconds:.1f}s",
            flush=True,
        )

        trails.inventory = None
        trails.characterized_inventory = None
        trails.scores = None
        trails.graph = None
        gc.collect()

    comparison = _comparison_rows(adaptive_rows, fixed_rows)
    _write_csv(adaptive_rows, output_csv)
    _write_csv(comparison, comparison_csv)
    if args.write_sankey or args.write_depth_sankey:
        _write_csv(sankey_summary_rows, sankey_summary_csv)
    print(f"Wrote adaptive CSV: {output_csv}", flush=True)
    print(f"Wrote comparison CSV: {comparison_csv}", flush=True)
    if args.write_sankey:
        print(f"Wrote activity Sankey directory: {sankey_dir}", flush=True)
    if args.write_depth_sankey:
        print(f"Wrote depth Sankey directory: {depth_sankey_dir}", flush=True)
    if args.write_sankey or args.write_depth_sankey:
        print(f"Wrote Sankey depth summary CSV: {sankey_summary_csv}", flush=True)
    if args.write_score_plots:
        print(f"Wrote temporal score plot directory: {score_plot_dir}", flush=True)
    return 0


def parse_adaptive_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run adaptive routing for the DACCS particulate-matter case."
    )
    parser.add_argument(
        "--datapackage", type=Path, default=depth_runner.DEFAULT_DATAPACKAGE
    )
    parser.add_argument(
        "--inventories",
        type=Path,
        nargs="+",
        default=depth_runner.DEFAULT_INVENTORY_PATHS,
    )
    parser.add_argument(
        "--lcia-json", type=Path, default=depth_runner.DEFAULT_LCIA_JSON
    )
    parser.add_argument("--depth-sweep-csv", type=Path, default=DEFAULT_DEPTH_SWEEP_CSV)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_ADAPTIVE_OUTPUT_CSV)
    parser.add_argument("--comparison-csv", type=Path, default=DEFAULT_COMPARISON_CSV)
    parser.add_argument("--sankey-dir", type=Path, default=DEFAULT_ACTIVITY_SANKEY_DIR)
    parser.add_argument("--depth-sankey-dir", type=Path, default=DEFAULT_SANKEY_DIR)
    parser.add_argument("--score-plot-dir", type=Path, default=DEFAULT_SCORE_PLOT_DIR)
    parser.add_argument(
        "--sankey-summary-csv",
        type=Path,
        default=DEFAULT_SANKEY_SUMMARY_CSV,
    )
    parser.add_argument("--method", default=depth_runner.DEFAULT_METHOD)
    parser.add_argument("--ei-version", default="3.12")
    parser.add_argument(
        "--reference-year", type=int, default=depth_runner.DEFAULT_REFERENCE_YEAR
    )
    parser.add_argument("--amount", type=float, default=depth_runner.DEFAULT_AMOUNT)
    parser.add_argument("--max-depth", type=int, default=7)
    parser.add_argument(
        "--no-max-depth",
        action="store_true",
        help=(
            "Disable the routing depth cap in adaptive mode. Branches then stop "
            "only by adaptive score cutoff, min_amount, or leaf behavior."
        ),
    )
    parser.add_argument(
        "--relative-cutoffs",
        type=float,
        nargs="+",
        default=[1e-2, 1e-3, 1e-4],
    )
    parser.add_argument("--adaptive-min-depth", type=int, default=1)
    parser.add_argument(
        "--routing-min-amount",
        type=float,
        default=depth_runner.DEFAULT_ROUTING_MIN_AMOUNT,
    )
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--no-cache-interpolation", action="store_true")
    parser.add_argument("--no-static-score-cache", action="store_true")
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
        "--write-sankey",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Write one activity-year Sankey HTML file per adaptive cutoff. "
            "Node and link sizes use static LCIA score potential."
        ),
    )
    parser.add_argument(
        "--write-depth-sankey",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also write the previous depth-count Sankey diagnostic.",
    )
    parser.add_argument(
        "--write-sankey-png",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also try to export PNG copies of the Sankey diagrams with kaleido.",
    )
    parser.add_argument(
        "--sankey-max-nodes",
        type=int,
        default=70,
        help="Maximum number of explicit graph nodes in each activity Sankey.",
    )
    parser.add_argument(
        "--sankey-max-links",
        type=int,
        default=220,
        help="Maximum number of aggregated links shown in each activity Sankey.",
    )
    parser.add_argument(
        "--sankey-min-link-relative-score",
        type=float,
        default=2e-4,
        help=(
            "Minimum activity Sankey link value as a fraction of the temporal "
            "score; lower values show more links."
        ),
    )
    parser.add_argument(
        "--sankey-min-label-relative-score",
        type=float,
        default=1e-2,
        help=(
            "Minimum displayed node value as a fraction of the temporal score "
            "before full activity labels are shown."
        ),
    )
    parser.add_argument(
        "--write-score-plots",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Write one plot_temporal_scores() HTML file per adaptive cutoff.",
    )
    parser.add_argument(
        "--write-score-plot-png",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Also try to export PNG copies of the temporal score plots.",
    )
    parser.add_argument(
        "--skip-lca",
        action="store_true",
        help=(
            "Only rerun routing and diagram generation. Existing score columns are "
            "reused when present."
        ),
    )
    parser.add_argument("--append", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()




def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in {"-h", "--help"}:
        print(
            "usage: daccs_pm_routing_diagnostics.py {depth-sweep,adaptive} [options]\n\n"
            "subcommands:\n"
            "  depth-sweep  reproduce the fixed-depth DACCS PM sensitivity table\n"
            "  adaptive     run adaptive routing diagnostics and optional Sankey plots\n\n"
            "Run a subcommand with --help for its options."
        )
        return 0

    command = argv.pop(0) if argv else "depth-sweep"
    if command == "depth-sweep":
        original_argv = sys.argv
        sys.argv = [original_argv[0], *argv]
        try:
            return run_depth_sweep(parse_depth_sweep_args())
        finally:
            sys.argv = original_argv
    if command == "adaptive":
        original_argv = sys.argv
        sys.argv = [original_argv[0], *argv]
        try:
            return run_adaptive(parse_adaptive_args())
        finally:
            sys.argv = original_argv

    raise SystemExit(
        f"Unknown subcommand {command!r}. Use 'depth-sweep' or 'adaptive'."
    )


if __name__ == "__main__":
    raise SystemExit(main())
