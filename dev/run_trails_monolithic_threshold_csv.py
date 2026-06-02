#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import gc
import os
import re
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import psutil
from datapackage import Package

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from trails import Trails

import benchmark_trails_matrix_solves as bench


DEFAULT_PACKAGE = Path(
    "/Users/romain/GitHub/trails/dev/trails_2026-03-18/datapackage.json"
)
DEFAULT_OUTPUT = Path("/Users/romain/GitHub/trails/dev/monolithic_vs_trails.csv")
SIZE_UNITS = {
    "B": 1,
    "KiB": 1024,
    "MiB": 1024**2,
    "GiB": 1024**3,
    "TiB": 1024**4,
}
UMFPACK_ESTIMATE_RE = re.compile(
    r"(numeric-size|peak-memory) estimate=([0-9.]+) (B|KiB|MiB|GiB|TiB)"
)
FIELDNAMES = [
    "run_id",
    "phase",
    "status",
    "horizon_start",
    "horizon_end",
    "year_blocks",
    "matrix_shape_rows",
    "matrix_shape_cols",
    "matrix_nnz",
    "temporal_exchanges_selected",
    "temporal_offdiag_entries_added",
    "temporal_offsets_dropped",
    "csc_storage_bytes",
    "matrix_build_seconds",
    "symbolic_seconds",
    "numeric_seconds",
    "solve_seconds",
    "total_factorize_solve_seconds",
    "routing_seconds",
    "lca_seconds",
    "routing_max_depth",
    "graph_nodes",
    "graph_edges",
    "umfpack_symbolic_size_bytes",
    "umfpack_symbolic_peak_bytes",
    "umfpack_numeric_size_bytes",
    "umfpack_peak_memory_bytes",
    "umfpack_variable_peak_bytes",
    "rss_before_bytes",
    "rss_after_bytes",
    "rss_delta_bytes",
    "rss_peak_bytes",
    "rss_peak_delta_bytes",
    "failure_type",
    "failure_message",
]


def _rss(process: psutil.Process) -> int:
    return int(process.memory_info().rss)


def _blank_row(run_id: str, phase: str, status: str) -> dict[str, Any]:
    return {key: "" for key in FIELDNAMES} | {
        "run_id": run_id,
        "phase": phase,
        "status": status,
    }


def _memory_fields(memory: bench.MemoryWindow) -> dict[str, Any]:
    return {
        "rss_before_bytes": memory.before,
        "rss_after_bytes": memory.after,
        "rss_delta_bytes": memory.delta,
        "rss_peak_bytes": memory.peak,
        "rss_peak_delta_bytes": memory.peak_delta,
    }


def _umfpack_failure_estimate_fields(message: str) -> dict[str, int]:
    fields: dict[str, int] = {}
    for estimate_name, value, unit in UMFPACK_ESTIMATE_RE.findall(message):
        bytes_value = int(round(float(value) * SIZE_UNITS[unit]))
        if estimate_name == "numeric-size":
            fields["umfpack_numeric_size_bytes"] = bytes_value
        elif estimate_name == "peak-memory":
            fields["umfpack_peak_memory_bytes"] = bytes_value
    return fields


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def _load_trails(package_path: Path) -> Trails:
    return Trails(Package(str(package_path)), interpolate_annual=True, debug=False)


def _cases_from_routing_with_stats(
    *,
    trails: Trails,
    start_year: int,
    activity_index: int,
    amount: float,
    max_depth: int,
    selected_years: set[int],
) -> tuple[list[bench.YearCase], float, int, int]:
    start = time.perf_counter()
    trails.temporal_routing(
        start_year=int(start_year),
        start_act_idx=int(activity_index),
        amount=float(amount),
        max_depth=int(max_depth),
        show_progress=False,
        attribute_to_roots=False,
    )
    routing_seconds = time.perf_counter() - start

    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError("Temporal routing did not initialize trails.graph")

    frontier: dict[tuple[int, int], float] = {}
    for _node, data in graph.nodes(data=True):
        frontier_amount = float(data.get("frontier_amount") or 0.0)
        if not frontier_amount:
            continue
        year = int(data.get("year"))
        if year not in selected_years:
            continue
        activity = int(data.get("act_idx"))
        key = (year, activity)
        frontier[key] = float(frontier.get(key, 0.0)) + frontier_amount

    frontier_by_year = trails.frontier_to_demand_vectors(frontier)
    cases: list[bench.YearCase] = []
    for raw_year in sorted(frontier_by_year):
        demand = np.asarray(frontier_by_year[raw_year], dtype=np.float64)
        if np.count_nonzero(demand) == 0:
            continue
        cases.append(
            bench.YearCase(
                raw_year=int(raw_year),
                time_index=bench._scenario_time_index(trails, int(raw_year)),
                activity_demand=demand,
            )
        )

    if not cases:
        raise RuntimeError("Temporal routing produced no selected frontier demands.")

    return (
        cases,
        routing_seconds,
        int(graph.number_of_nodes()),
        int(graph.number_of_edges()),
    )


def _run_paired_horizon(
    *,
    trails: Trails,
    process: psutil.Process,
    horizon_start: int,
    horizon_end: int,
    activity_index: int,
    amount: float,
    max_depth: int,
    rows: list[dict[str, Any]],
) -> bool:
    run_id = f"{horizon_start}-{horizon_end}"
    years = list(range(int(horizon_start), int(horizon_end) + 1))
    time_indices = bench._selected_time_indices(trails, years)

    print(f"\nPaired horizon {run_id} ({len(years)} blocks)", flush=True)

    before = _rss(process)
    start = time.perf_counter()
    with bench.PeakRSSMonitor(process) as monitor:
        matrix, expanded_stats = bench._build_time_expanded_csc(trails, time_indices)
        rhs = bench._build_single_start_rhs(
            matrix,
            trails,
            time_indices,
            start_year=int(horizon_end),
            activity_index=int(activity_index),
            amount=float(amount),
        )
    matrix_build_seconds = time.perf_counter() - start
    after = _rss(process)
    build_memory = bench.MemoryWindow(before=before, after=after, peak=monitor.peak)
    csc_storage = bench._storage_total(bench._csc_storage_bytes(matrix))

    row = _blank_row(run_id, "monolithic_matrix_build", "ok")
    row.update(
        {
            "horizon_start": horizon_start,
            "horizon_end": horizon_end,
            "year_blocks": len(years),
            "matrix_shape_rows": matrix.shape[0],
            "matrix_shape_cols": matrix.shape[1],
            "matrix_nnz": int(matrix.nnz),
            "temporal_exchanges_selected": expanded_stats.td_exchanges_selected,
            "temporal_offdiag_entries_added": expanded_stats.td_offsets_added,
            "temporal_offsets_dropped": expanded_stats.td_offsets_dropped,
            "csc_storage_bytes": csc_storage,
            "matrix_build_seconds": matrix_build_seconds,
        }
    )
    row.update(_memory_fields(build_memory))
    rows.append(row)
    print(
        "  built matrix "
        f"shape={matrix.shape}, nnz={matrix.nnz:,}, "
        f"storage={bench._format_bytes(csc_storage)} in {matrix_build_seconds:.3f}s",
        flush=True,
    )

    monolithic_ok = True
    monolithic_solution = None
    try:
        monolithic_solution, solve_stats = bench._solve_umfpack(
            matrix,
            rhs,
            repeat=1,
            process=process,
        )
    except (MemoryError, RuntimeError) as exc:
        monolithic_ok = False
        row = _blank_row(run_id, "monolithic_lu_solve", "failed")
        row.update(
            {
                "horizon_start": horizon_start,
                "horizon_end": horizon_end,
                "year_blocks": len(years),
                "matrix_shape_rows": matrix.shape[0],
                "matrix_shape_cols": matrix.shape[1],
                "matrix_nnz": int(matrix.nnz),
                "csc_storage_bytes": csc_storage,
                "failure_type": type(exc).__name__,
                "failure_message": str(exc),
                "rss_after_bytes": _rss(process),
            }
        )
        row.update(_umfpack_failure_estimate_fields(str(exc)))
        rows.append(row)
        print(f"  monolithic LU failed: {exc}", flush=True)
    else:
        row = _blank_row(run_id, "monolithic_lu_solve", "ok")
        row.update(
            {
                "horizon_start": horizon_start,
                "horizon_end": horizon_end,
                "year_blocks": len(years),
                "matrix_shape_rows": matrix.shape[0],
                "matrix_shape_cols": matrix.shape[1],
                "matrix_nnz": int(matrix.nnz),
                "csc_storage_bytes": csc_storage,
                "symbolic_seconds": solve_stats.symbolic_seconds,
                "numeric_seconds": solve_stats.numeric_seconds,
                "solve_seconds": sum(solve_stats.solve_seconds),
                "total_factorize_solve_seconds": (
                    solve_stats.first_factorize_solve_seconds
                ),
                "umfpack_symbolic_size_bytes": solve_stats.symbolic_size,
                "umfpack_symbolic_peak_bytes": solve_stats.symbolic_peak,
                "umfpack_numeric_size_bytes": solve_stats.numeric_size,
                "umfpack_peak_memory_bytes": solve_stats.peak_memory,
                "umfpack_variable_peak_bytes": solve_stats.variable_peak,
            }
        )
        row.update(_memory_fields(solve_stats.memory))
        rows.append(row)
        print(
            "  monolithic LU solved in "
            f"{solve_stats.first_factorize_solve_seconds:.3f}s",
            flush=True,
        )

    del matrix, rhs
    if monolithic_solution is not None:
        del monolithic_solution
    gc.collect()

    cases, routing_seconds, graph_nodes, graph_edges = (
        _cases_from_routing_with_stats(
            trails=trails,
            start_year=int(horizon_end),
            activity_index=int(activity_index),
            amount=float(amount),
            max_depth=int(max_depth),
            selected_years=set(years),
        )
    )
    print(
        "  TRAILS routing finished in "
        f"{routing_seconds:.3f}s over {graph_nodes:,} nodes",
        flush=True,
    )

    sequential_solution, sequential_stats = bench._solve_sequential_umfpack(
        trails,
        cases,
        process=process,
        reuse_symbolic=True,
    )
    del sequential_solution
    row = _blank_row(run_id, "trails_yearwise_selected_solve", "ok")
    row.update(
        {
            "horizon_start": horizon_start,
            "horizon_end": horizon_end,
            "year_blocks": len(years),
            "matrix_shape_rows": int(trails.A.shape[2]),
            "matrix_shape_cols": int(trails.A.shape[1]),
            "csc_storage_bytes": sequential_stats.max_csc_storage,
            "matrix_build_seconds": sequential_stats.matrix_build_seconds,
            "routing_seconds": routing_seconds,
            "routing_max_depth": int(max_depth),
            "graph_nodes": graph_nodes,
            "graph_edges": graph_edges,
            "symbolic_seconds": sequential_stats.symbolic_seconds,
            "numeric_seconds": sequential_stats.numeric_seconds,
            "solve_seconds": sequential_stats.solve_seconds,
            "total_factorize_solve_seconds": (
                sequential_stats.symbolic_seconds
                + sequential_stats.numeric_seconds
                + sequential_stats.solve_seconds
            ),
            "umfpack_symbolic_size_bytes": sequential_stats.max_symbolic_size,
            "umfpack_symbolic_peak_bytes": sequential_stats.max_symbolic_peak,
            "umfpack_numeric_size_bytes": sequential_stats.max_numeric_size,
            "umfpack_peak_memory_bytes": sequential_stats.max_factorization_peak,
            "umfpack_variable_peak_bytes": sequential_stats.max_variable_peak,
        }
    )
    row.update(_memory_fields(sequential_stats.memory))
    rows.append(row)
    print(
        "  TRAILS-style selected year-wise solve finished in "
        f"{row['total_factorize_solve_seconds']:.3f}s",
        flush=True,
    )

    gc.collect()
    return monolithic_ok


def _run_full_trails(
    *,
    trails: Trails,
    process: psutil.Process,
    start_year: int,
    activity_index: int,
    amount: float,
    max_depth: int,
    rows: list[dict[str, Any]],
) -> None:
    print("\nFinal full TRAILS run", flush=True)
    full_time_indices = list(range(len(trails.scenario_labels)))
    temporal_stats = bench._time_expanded_temporal_stats(trails, full_time_indices)

    before = _rss(process)
    routing_seconds: float | str = ""
    lca_seconds: float | str = ""
    status = "ok"
    failure_type = ""
    failure_message = ""
    with bench.PeakRSSMonitor(process) as monitor:
        try:
            start = time.perf_counter()
            trails.temporal_routing(
                start_year=int(start_year),
                start_act_idx=int(activity_index),
                amount=float(amount),
                max_depth=int(max_depth),
                show_progress=False,
                attribute_to_roots=True,
            )
            routing_seconds = time.perf_counter() - start

            start = time.perf_counter()
            trails.lca(
                methods=[],
                show_progress=False,
                compute_score=False,
                store_inventory=False,
                attribute_to_roots=True,
                solver_mode="direct",
            )
            lca_seconds = time.perf_counter() - start
        except Exception as exc:
            status = "failed"
            failure_type = type(exc).__name__
            failure_message = str(exc)

    after = _rss(process)
    memory = bench.MemoryWindow(before=before, after=after, peak=monitor.peak)
    graph = getattr(trails, "graph", None)
    graph_message = (
        f"graph_nodes={graph.number_of_nodes()} "
        f"graph_edges={graph.number_of_edges()}"
        if graph is not None
        else ""
    )
    if failure_message and graph_message:
        failure_message = f"{failure_message}; {graph_message}"
    elif graph_message:
        failure_message = graph_message

    rows.append(
        _blank_row("full", "trails_full_routing_yearwise_solve", status)
        | {
            "horizon_start": min(int(label) for label in trails.scenario_labels),
            "horizon_end": max(int(label) for label in trails.scenario_labels),
            "year_blocks": len(trails.scenario_labels),
            "matrix_shape_rows": int(trails.A.shape[2]),
            "matrix_shape_cols": int(trails.A.shape[1]),
            "matrix_nnz": int(trails.A.nnz),
            "temporal_exchanges_selected": (
                temporal_stats.td_exchanges_selected
            ),
            "temporal_offdiag_entries_added": temporal_stats.td_offsets_added,
            "temporal_offsets_dropped": temporal_stats.td_offsets_dropped,
            "csc_storage_bytes": bench._storage_total(
                bench._coo_storage_bytes(trails.A)
            ),
            "routing_seconds": routing_seconds,
            "lca_seconds": lca_seconds,
            "routing_max_depth": int(max_depth),
            "graph_nodes": graph.number_of_nodes() if graph is not None else "",
            "graph_edges": graph.number_of_edges() if graph is not None else "",
            "failure_type": failure_type,
            "failure_message": failure_message,
        }
        | _memory_fields(memory)
    )
    if status == "ok":
        print(
            f"  full TRAILS routing={routing_seconds:.3f}s, "
            f"lca={lca_seconds:.3f}s",
            flush=True,
        )
    else:
        print(f"  full TRAILS failed: {failure_type}: {failure_message}", flush=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run paired monolithic time-expanded LU and TRAILS year-wise solves "
            "until monolithic LU fails, then run a final full TRAILS case."
        )
    )
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--activity-index", type=int, default=31387)
    parser.add_argument("--start-year", type=int, default=2025)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument(
        "--skip-paired",
        action="store_true",
        help="Skip the monolithic threshold pairs and run only the full TRAILS case.",
    )
    parser.add_argument(
        "--horizon-sizes",
        nargs="+",
        type=int,
        default=[1, 2, 3],
        help="Contiguous horizon sizes ending at --start-year.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    process = psutil.Process(os.getpid())
    rows: list[dict[str, Any]] = []
    package_path = Path(args.package).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    print(f"PID: {process.pid}", flush=True)
    print(f"Package: {package_path}", flush=True)
    print(f"Output CSV: {output_path}", flush=True)

    load_before = _rss(process)
    start = time.perf_counter()
    with bench.PeakRSSMonitor(process) as monitor:
        trails = _load_trails(package_path)
    load_after = _rss(process)
    load_memory = bench.MemoryWindow(
        before=load_before,
        after=load_after,
        peak=monitor.peak,
    )
    row = _blank_row("load", "load_trails", "ok")
    row.update(
        {
            "year_blocks": len(trails.scenario_labels),
            "matrix_shape_rows": int(trails.A.shape[2]),
            "matrix_shape_cols": int(trails.A.shape[1]),
            "matrix_nnz": int(trails.A.nnz),
            "matrix_build_seconds": time.perf_counter() - start,
            "csc_storage_bytes": bench._storage_total(
                bench._coo_storage_bytes(trails.A)
            ),
        }
    )
    row.update(_memory_fields(load_memory))
    rows.append(row)

    if args.skip_paired:
        print("Skipping paired monolithic threshold runs.", flush=True)
    else:
        first_failure_seen = False
        for size in args.horizon_sizes:
            horizon_end = int(args.start_year)
            horizon_start = horizon_end - int(size) + 1
            ok = _run_paired_horizon(
                trails=trails,
                process=process,
                horizon_start=horizon_start,
                horizon_end=horizon_end,
            activity_index=int(args.activity_index),
            amount=float(args.amount),
            max_depth=int(args.max_depth),
            rows=rows,
        )
            _write_rows(output_path, rows)
            if not ok:
                first_failure_seen = True
                break

        if not first_failure_seen:
            print(
                "No monolithic failure observed in requested horizon sizes.",
                flush=True,
            )

    _run_full_trails(
        trails=trails,
        process=process,
        start_year=int(args.start_year),
        activity_index=int(args.activity_index),
        amount=float(args.amount),
        max_depth=int(args.max_depth),
        rows=rows,
    )
    _write_rows(output_path, rows)
    print(f"\nWrote CSV: {output_path}", flush=True)


if __name__ == "__main__":
    main()
