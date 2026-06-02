#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import psutil
from datapackage import Package
from scipy import sparse as sp

try:
    from scikits.umfpack import UMFPACK_A, UmfpackContext
    import scikits.umfpack as umfpack
except ImportError as exc:  # pragma: no cover - environment dependent
    raise RuntimeError(
        "scikit-umfpack is required for this benchmark. Run it in the `trails` "
        "conda environment or install scikit-umfpack first."
    ) from exc

from trails import Trails

DEFAULT_PACKAGE = Path("/Users/romain/GitHub/trails/dev/trails_2026-03-18.zip")
DEFAULT_METHOD = (
    "IPCC 2021 (incl. biogenic CO2) - climate change: total "
    "(incl. biogenic CO2) - global warming potential (GWP100)"
)


@dataclass
class MemoryWindow:
    before: int
    after: int
    peak: int

    @property
    def delta(self) -> int:
        return self.after - self.before

    @property
    def peak_delta(self) -> int:
        return self.peak - self.before


@dataclass
class UmfpackStats:
    symbolic_seconds: float
    numeric_seconds: float
    solve_seconds: list[float]
    symbolic_size: int | None
    symbolic_peak: int | None
    numeric_size: int | None
    peak_memory: int | None
    variable_peak: int | None
    memory: MemoryWindow

    @property
    def first_factorize_solve_seconds(self) -> float:
        return self.symbolic_seconds + self.numeric_seconds + self.solve_seconds[0]


@dataclass
class YearCase:
    raw_year: int
    time_index: int
    activity_demand: np.ndarray


@dataclass
class TimeExpandedBuildStats:
    td_exchanges_selected: int
    td_offsets_added: int
    td_offsets_dropped: int
    td_same_year_cancellations: int


@dataclass
class SequentialStats:
    years: list[int]
    matrix_build_seconds: float
    symbolic_seconds: float
    numeric_seconds: float
    solve_seconds: float
    max_csc_storage: int
    sum_csc_storage: int
    max_symbolic_size: int | None
    max_symbolic_peak: int | None
    max_numeric_size: int | None
    sum_numeric_size: int
    max_factorization_peak: int | None
    max_variable_peak: int | None
    memory: MemoryWindow


class PeakRSSMonitor:
    def __init__(self, process: psutil.Process, interval_s: float = 0.01) -> None:
        self.process = process
        self.interval_s = float(interval_s)
        self.peak = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _loop(self) -> None:
        while not self._stop.is_set():
            try:
                rss = int(self.process.memory_info().rss)
            except psutil.Error:
                break
            self.peak = max(self.peak, rss)
            self._stop.wait(self.interval_s)

    def __enter__(self) -> "PeakRSSMonitor":
        self.peak = int(self.process.memory_info().rss)
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *args: object) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)


def _rss(process: psutil.Process) -> int:
    return int(process.memory_info().rss)


def _format_bytes(value: int | float | None) -> str:
    if value is None:
        return "n/a"
    value = float(value)
    sign = "-" if value < 0 else ""
    value = abs(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{sign}{value:.3f} {unit}"
        value /= 1024.0
    return f"{sign}{value:.3f} TiB"


def _csc_storage_bytes(matrix: sp.csc_matrix) -> dict[str, int]:
    return {
        "data": int(matrix.data.nbytes),
        "indices": int(matrix.indices.nbytes),
        "indptr": int(matrix.indptr.nbytes),
    }


def _coo_storage_bytes(matrix: Any) -> dict[str, int]:
    return {
        "data": int(getattr(matrix.data, "nbytes", 0)),
        "coords": int(getattr(matrix.coords, "nbytes", 0)),
    }


def _storage_total(parts: dict[str, int]) -> int:
    return int(sum(parts.values()))


def _umfpack_info_bytes(ctx: UmfpackContext, name: str) -> int | None:
    index = getattr(umfpack, name, None)
    unit_index = getattr(umfpack, "UMFPACK_SIZE_OF_UNIT", None)
    if index is None or unit_index is None:
        return None
    value = float(ctx.info[index])
    unit_size = float(ctx.info[unit_index])
    if value < 0 or unit_size <= 0:
        return None
    return int(round(value * unit_size))


def _matrix_signature(matrix: sp.csc_matrix) -> tuple[int, int, int, int]:
    return (
        int(matrix.indptr.size),
        int(matrix.indices.size),
        int(np.bitwise_xor.reduce(matrix.indptr)),
        int(np.bitwise_xor.reduce(matrix.indices)),
    )


def _selected_time_indices(trails: Trails, years: list[int] | None) -> list[int]:
    labels = [str(label) for label in trails.scenario_labels]
    if years is None:
        return list(range(len(labels)))

    label_to_t = {label: idx for idx, label in enumerate(labels)}
    selected: list[int] = []
    missing: list[str] = []
    for year in years:
        label = str(int(year))
        if label in label_to_t:
            selected.append(int(label_to_t[label]))
        else:
            missing.append(label)
    if missing:
        raise ValueError(f"Year labels not present in trails.A: {', '.join(missing)}")
    return selected


def _time_labels(trails: Trails, time_indices: list[int]) -> list[int]:
    return [int(trails.scenario_labels[t]) for t in time_indices]


def _scenario_time_index(trails: Trails, raw_year: int) -> int:
    context = trails._get_scenario_context(int(raw_year))
    if context is None:
        raise RuntimeError(f"No scenario context available for year={raw_year}")
    _scenario_year, _label, time_index = context
    return int(time_index)


def _cases_from_routing(
    trails: Trails,
    *,
    start_year: int,
    activity_index: int,
    amount: float,
    max_depth: int,
    selected_years: set[int] | None,
    show_progress: bool,
) -> list[YearCase]:
    start = time.perf_counter()
    trails.temporal_routing(
        start_year=int(start_year),
        start_act_idx=int(activity_index),
        amount=float(amount),
        max_depth=int(max_depth),
        show_progress=bool(show_progress),
        attribute_to_roots=False,
    )
    print(f"Temporal routing completed in {time.perf_counter() - start:.6f} s")

    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError("Temporal routing did not initialize trails.graph")

    frontier: dict[tuple[int, int], float] = {}
    for _node, data in graph.nodes(data=True):
        frontier_amount = float(data.get("frontier_amount") or 0.0)
        if not frontier_amount:
            continue
        year = int(data.get("year"))
        activity = int(data.get("act_idx"))
        key = (year, activity)
        frontier[key] = float(frontier.get(key, 0.0)) + frontier_amount

    frontier_by_year = trails.frontier_to_demand_vectors(frontier)
    cases: list[YearCase] = []
    for raw_year in sorted(frontier_by_year):
        if selected_years is not None and int(raw_year) not in selected_years:
            continue
        demand = np.asarray(frontier_by_year[raw_year], dtype=np.float64)
        if np.count_nonzero(demand) == 0:
            continue
        cases.append(
            YearCase(
                raw_year=int(raw_year),
                time_index=_scenario_time_index(trails, int(raw_year)),
                activity_demand=demand,
            )
        )

    if not cases:
        raise RuntimeError("Temporal routing produced no nonzero frontier demands.")
    return cases


def _cases_from_activity_each_year(
    trails: Trails,
    time_indices: list[int],
    *,
    activity_index: int,
    amount: float,
) -> list[YearCase]:
    if trails.A is None:
        raise RuntimeError("trails.A is None")
    n_activities = int(trails.A.shape[1])
    cases: list[YearCase] = []
    for time_index in time_indices:
        demand = np.zeros(n_activities, dtype=np.float64)
        demand[int(activity_index)] = float(amount)
        cases.append(
            YearCase(
                raw_year=int(trails.scenario_labels[int(time_index)]),
                time_index=int(time_index),
                activity_demand=demand,
            )
        )
    return cases


def _reference_product_from_slice(
    trails: Trails,
    time_index: int,
    activity_index: int,
) -> tuple[int, float]:
    if trails.A is None:
        raise RuntimeError("trails.A is None")

    A_t = trails.A[int(time_index), :, :]
    coords = A_t.coords
    if coords.shape[0] == 3:
        act_idx = np.asarray(coords[1], dtype=np.int64)
        prod_idx = np.asarray(coords[2], dtype=np.int64)
    elif coords.shape[0] == 2:
        act_idx = np.asarray(coords[0], dtype=np.int64)
        prod_idx = np.asarray(coords[1], dtype=np.int64)
    else:
        raise ValueError(f"Unsupported A_t coords ndim={coords.shape[0]}")

    mask = act_idx == int(activity_index)
    if not np.any(mask):
        raise ValueError(
            f"Activity index {activity_index} has no technosphere column entries "
            f"in time index {time_index}."
        )

    rows = prod_idx[mask]
    vals = np.asarray(A_t.data[mask], dtype=np.float64)

    diagonal = rows == int(activity_index)
    if np.any(diagonal):
        pos = int(np.flatnonzero(diagonal)[0])
    else:
        pos = int(np.argmin(np.abs(np.abs(vals) - 1.0)))

    return int(rows[pos]), float(vals[pos])


def _reference_product_from_csc(
    matrix: sp.csc_matrix,
    activity_col: int,
    *,
    preferred_row: int,
) -> tuple[int, float]:
    col = matrix.getcol(int(activity_col))
    rows = col.indices
    vals = np.asarray(col.data, dtype=np.float64)
    if rows.size == 0:
        raise ValueError(f"Activity column {activity_col} has no entries")

    preferred = rows == int(preferred_row)
    if np.any(preferred):
        pos = int(np.flatnonzero(preferred)[0])
    else:
        pos = int(np.argmin(np.abs(np.abs(vals) - 1.0)))
    return int(rows[pos]), float(vals[pos])


def _activity_demand_to_product_rhs(
    matrix: sp.csc_matrix,
    activity_demand: np.ndarray,
    *,
    activity_offset: int = 0,
    product_offset: int = 0,
) -> dict[int, float]:
    rhs: dict[int, float] = {}
    nonzero = np.flatnonzero(activity_demand != 0.0)
    for activity in nonzero:
        global_activity = int(activity_offset + int(activity))
        preferred_product = int(product_offset + int(activity))
        product, product_value = _reference_product_from_csc(
            matrix,
            global_activity,
            preferred_row=preferred_product,
        )
        sign = -1.0 if product_value < 0.0 else 1.0
        rhs[int(product)] = rhs.get(int(product), 0.0) + float(
            activity_demand[int(activity)]
        ) * sign
    return rhs


def _build_block_rhs(
    matrix: sp.csc_matrix,
    cases: list[YearCase],
    *,
    n_activities: int,
    n_products: int,
) -> np.ndarray:
    rhs = np.zeros(matrix.shape[0], dtype=np.float64)
    for local_t, case in enumerate(cases):
        mapped = _activity_demand_to_product_rhs(
            matrix,
            case.activity_demand,
            activity_offset=local_t * n_activities,
            product_offset=local_t * n_products,
        )
        for product, amount in mapped.items():
            rhs[int(product)] += float(amount)
    return rhs


def _build_single_start_rhs(
    matrix: sp.csc_matrix,
    trails: Trails,
    time_indices: list[int],
    *,
    start_year: int,
    activity_index: int,
    amount: float,
) -> np.ndarray:
    if trails.A is None:
        raise RuntimeError("trails.A is None")
    n_activities = int(trails.A.shape[1])
    n_products = int(trails.A.shape[2])
    start_t = _scenario_time_index(trails, int(start_year))
    try:
        local_t = time_indices.index(start_t)
    except ValueError as exc:
        raise ValueError(
            f"Mapped start year {start_year} to time index {start_t}, "
            "which is not in the selected matrix horizon."
        ) from exc

    activity_demand = np.zeros(n_activities, dtype=np.float64)
    activity_demand[int(activity_index)] = float(amount)
    mapped = _activity_demand_to_product_rhs(
        matrix,
        activity_demand,
        activity_offset=local_t * n_activities,
        product_offset=local_t * n_products,
    )
    rhs = np.zeros(matrix.shape[0], dtype=np.float64)
    for product, value in mapped.items():
        rhs[int(product)] += float(value)
    return rhs


def _build_year_csc(trails: Trails, time_index: int) -> sp.csc_matrix:
    if trails.A is None:
        raise RuntimeError("trails.A is None")

    A_t = trails.A[int(time_index), :, :]
    coords = A_t.coords
    if coords.shape[0] == 3:
        act_idx = np.asarray(coords[1], dtype=np.int64)
        prod_idx = np.asarray(coords[2], dtype=np.int64)
    elif coords.shape[0] == 2:
        act_idx = np.asarray(coords[0], dtype=np.int64)
        prod_idx = np.asarray(coords[1], dtype=np.int64)
    else:
        raise ValueError(f"Unsupported A_t coords ndim={coords.shape[0]}")

    data = np.asarray(A_t.data, dtype=np.float64)
    n_activities = int(trails.A.shape[1])
    n_products = int(trails.A.shape[2])
    matrix = sp.coo_matrix(
        (data, (prod_idx, act_idx)), shape=(n_products, n_activities)
    ).tocsc()
    matrix.sort_indices()
    return matrix


def _build_block_csc(trails: Trails, cases: list[YearCase]) -> sp.csc_matrix:
    if trails.A is None:
        raise RuntimeError("trails.A is None")
    if trails.A.shape[1] != trails.A.shape[2]:
        raise ValueError(f"trails.A slices are not square: {trails.A.shape}")

    A = trails.A
    n_activities = int(A.shape[1])
    n_products = int(A.shape[2])
    time_indices = [int(case.time_index) for case in cases]

    if len(set(time_indices)) == len(time_indices):
        selected = np.asarray(time_indices, dtype=np.int64)
        position_by_t = np.full(int(A.shape[0]), -1, dtype=np.int64)
        position_by_t[selected] = np.arange(selected.size, dtype=np.int64)

        all_t = np.asarray(A.coords[0], dtype=np.int64)
        local_t = position_by_t[all_t]
        mask = local_t >= 0

        local_t = local_t[mask]
        act_idx = np.asarray(A.coords[1, mask], dtype=np.int64)
        prod_idx = np.asarray(A.coords[2, mask], dtype=np.int64)
        data = np.asarray(A.data[mask], dtype=np.float64)
    else:
        row_parts: list[np.ndarray] = []
        col_parts: list[np.ndarray] = []
        data_parts: list[np.ndarray] = []
        all_t = np.asarray(A.coords[0], dtype=np.int64)
        for local_pos, time_index in enumerate(time_indices):
            mask = all_t == int(time_index)
            row_parts.append(
                local_pos * n_products
                + np.asarray(A.coords[2, mask], dtype=np.int64)
            )
            col_parts.append(
                local_pos * n_activities
                + np.asarray(A.coords[1, mask], dtype=np.int64)
            )
            data_parts.append(np.asarray(A.data[mask], dtype=np.float64))

        row = np.concatenate(row_parts)
        col = np.concatenate(col_parts)
        data = np.concatenate(data_parts)
        shape = (len(cases) * n_products, len(cases) * n_activities)
        matrix = sp.coo_matrix((data, (row, col)), shape=shape).tocsc()
        matrix.sort_indices()
        return matrix

    row = local_t * n_products + prod_idx
    col = local_t * n_activities + act_idx

    shape = (len(cases) * n_products, len(cases) * n_activities)
    matrix = sp.coo_matrix((data, (row, col)), shape=shape).tocsc()
    matrix.sort_indices()
    return matrix


def _build_time_expanded_csc(
    trails: Trails,
    time_indices: list[int],
) -> tuple[sp.csc_matrix, TimeExpandedBuildStats]:
    """Build a bw_timex-like all-year matrix with temporal off-diagonal links."""
    if trails.A is None:
        raise RuntimeError("trails.A is None")
    if trails.A.shape[1] != trails.A.shape[2]:
        raise ValueError(f"trails.A slices are not square: {trails.A.shape}")

    A = trails.A
    n_activities = int(A.shape[1])
    n_products = int(A.shape[2])
    selected = np.asarray(time_indices, dtype=np.int64)
    position_by_t = np.full(int(A.shape[0]), -1, dtype=np.int64)
    position_by_t[selected] = np.arange(selected.size, dtype=np.int64)

    all_t = np.asarray(A.coords[0], dtype=np.int64)
    local_t = position_by_t[all_t]
    mask = local_t >= 0

    base_local_t = local_t[mask]
    base_act = np.asarray(A.coords[1, mask], dtype=np.int64)
    base_prod = np.asarray(A.coords[2, mask], dtype=np.int64)
    base_data = np.asarray(A.data[mask], dtype=np.float64)

    row_parts = [base_local_t * n_products + base_prod]
    col_parts = [base_local_t * n_activities + base_act]
    data_parts = [base_data]

    label_to_t = {str(label): int(t) for t, label in enumerate(trails.scenario_labels)}
    t_to_local = {int(t): i for i, t in enumerate(time_indices)}

    cancel_rows: list[int] = []
    cancel_cols: list[int] = []
    cancel_data: list[float] = []
    td_rows: list[int] = []
    td_cols: list[int] = []
    td_data: list[float] = []

    td_selected = 0
    td_offsets_added = 0
    td_offsets_dropped = 0

    temporal_exchanges = trails.temporal_technosphere_exchanges.items()
    for (label, act_idx, prod_idx), tex in temporal_exchanges:
        source_t = label_to_t.get(str(label))
        if source_t is None or source_t not in t_to_local:
            continue

        source_local = int(t_to_local[source_t])
        act_idx = int(act_idx)
        prod_idx = int(prod_idx)
        anchor_value = float(A[int(source_t), act_idx, prod_idx])
        if anchor_value == 0.0:
            continue

        td_selected += 1
        source_year = int(trails.scenario_labels[int(source_t)])
        col = source_local * n_activities + act_idx
        same_year_row = source_local * n_products + prod_idx

        cancel_rows.append(same_year_row)
        cancel_cols.append(col)
        cancel_data.append(-anchor_value)

        offsets = trails._get_td_offsets(tex=tex, debug=False)
        amount_source = getattr(tex, "amount_source", "port")
        for offset, weight in offsets:
            raw_year = int(source_year + int(offset))
            target_year = int(trails._map_year_to_scenario_year(raw_year))
            target_t = label_to_t.get(str(target_year))
            if target_t is None or target_t not in t_to_local:
                td_offsets_dropped += 1
                continue

            if amount_source == "matrix":
                exchange_value = float(A[int(target_t), act_idx, prod_idx])
                if exchange_value == 0.0:
                    td_offsets_dropped += 1
                    continue
            else:
                exchange_value = anchor_value

            target_local = int(t_to_local[target_t])
            td_rows.append(target_local * n_products + prod_idx)
            td_cols.append(col)
            td_data.append(float(exchange_value) * float(weight))
            td_offsets_added += 1

    if cancel_rows:
        row_parts.append(np.asarray(cancel_rows, dtype=np.int64))
        col_parts.append(np.asarray(cancel_cols, dtype=np.int64))
        data_parts.append(np.asarray(cancel_data, dtype=np.float64))
    if td_rows:
        row_parts.append(np.asarray(td_rows, dtype=np.int64))
        col_parts.append(np.asarray(td_cols, dtype=np.int64))
        data_parts.append(np.asarray(td_data, dtype=np.float64))

    rows = np.concatenate(row_parts)
    cols = np.concatenate(col_parts)
    data = np.concatenate(data_parts)
    shape = (len(time_indices) * n_products, len(time_indices) * n_activities)
    matrix = sp.coo_matrix((data, (rows, cols)), shape=shape).tocsc()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    stats = TimeExpandedBuildStats(
        td_exchanges_selected=td_selected,
        td_offsets_added=td_offsets_added,
        td_offsets_dropped=td_offsets_dropped,
        td_same_year_cancellations=len(cancel_rows),
    )
    return matrix, stats


def _time_expanded_temporal_stats(
    trails: Trails,
    time_indices: list[int],
) -> TimeExpandedBuildStats:
    """Count temporal links without materializing the time-expanded matrix."""
    if trails.A is None:
        raise RuntimeError("trails.A is None")

    A = trails.A
    label_to_t = {str(label): int(t) for t, label in enumerate(trails.scenario_labels)}
    t_to_local = {int(t): i for i, t in enumerate(time_indices)}

    td_selected = 0
    td_offsets_added = 0
    td_offsets_dropped = 0

    temporal_exchanges = trails.temporal_technosphere_exchanges.items()
    for (label, act_idx, prod_idx), tex in temporal_exchanges:
        source_t = label_to_t.get(str(label))
        if source_t is None or source_t not in t_to_local:
            continue

        act_idx = int(act_idx)
        prod_idx = int(prod_idx)
        anchor_value = float(A[int(source_t), act_idx, prod_idx])
        if anchor_value == 0.0:
            continue

        td_selected += 1
        source_year = int(trails.scenario_labels[int(source_t)])
        offsets = trails._get_td_offsets(tex=tex, debug=False)
        amount_source = getattr(tex, "amount_source", "port")
        for offset, _weight in offsets:
            raw_year = int(source_year + int(offset))
            target_year = int(trails._map_year_to_scenario_year(raw_year))
            target_t = label_to_t.get(str(target_year))
            if target_t is None or target_t not in t_to_local:
                td_offsets_dropped += 1
                continue

            if amount_source == "matrix":
                exchange_value = float(A[int(target_t), act_idx, prod_idx])
                if exchange_value == 0.0:
                    td_offsets_dropped += 1
                    continue

            td_offsets_added += 1

    return TimeExpandedBuildStats(
        td_exchanges_selected=td_selected,
        td_offsets_added=td_offsets_added,
        td_offsets_dropped=td_offsets_dropped,
        td_same_year_cancellations=td_selected,
    )


def _solve_umfpack(
    matrix: sp.csc_matrix,
    rhs: np.ndarray,
    *,
    repeat: int,
    process: psutil.Process,
) -> tuple[np.ndarray, UmfpackStats]:
    ctx = UmfpackContext()
    before = _rss(process)
    with PeakRSSMonitor(process) as monitor:
        start = time.perf_counter()
        symbolic = ctx.symbolic(matrix)
        symbolic_seconds = time.perf_counter() - start

        start = time.perf_counter()
        try:
            ctx.numeric(matrix, symbolic)
        except TypeError:
            try:
                ctx.numeric(matrix)
            except RuntimeError as exc:
                peak_estimate = _umfpack_info_bytes(
                    ctx, "UMFPACK_PEAK_MEMORY_ESTIMATE"
                )
                numeric_estimate = _umfpack_info_bytes(
                    ctx, "UMFPACK_NUMERIC_SIZE_ESTIMATE"
                )
                raise RuntimeError(
                    f"{exc}; UMFPACK numeric-size estimate="
                    f"{_format_bytes(numeric_estimate)}, peak-memory estimate="
                    f"{_format_bytes(peak_estimate)}"
                ) from exc
        except RuntimeError as exc:
            peak_estimate = _umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY_ESTIMATE")
            numeric_estimate = _umfpack_info_bytes(ctx, "UMFPACK_NUMERIC_SIZE_ESTIMATE")
            raise RuntimeError(
                f"{exc}; UMFPACK numeric-size estimate="
                f"{_format_bytes(numeric_estimate)}, peak-memory estimate="
                f"{_format_bytes(peak_estimate)}"
            ) from exc
        numeric_seconds = time.perf_counter() - start

        solve_seconds: list[float] = []
        solution = np.empty_like(rhs)
        for _ in range(max(1, int(repeat))):
            start = time.perf_counter()
            solution = ctx.solve(UMFPACK_A, matrix, rhs)
            solve_seconds.append(time.perf_counter() - start)

    after = _rss(process)
    stats = UmfpackStats(
        symbolic_seconds=symbolic_seconds,
        numeric_seconds=numeric_seconds,
        solve_seconds=solve_seconds,
        symbolic_size=_umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_SIZE"),
        symbolic_peak=_umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_PEAK_MEMORY"),
        numeric_size=_umfpack_info_bytes(ctx, "UMFPACK_NUMERIC_SIZE"),
        peak_memory=_umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY"),
        variable_peak=_umfpack_info_bytes(ctx, "UMFPACK_VARIABLE_PEAK"),
        memory=MemoryWindow(before=before, after=after, peak=monitor.peak),
    )
    return solution, stats


def _solve_sequential_umfpack(
    trails: Trails,
    cases: list[YearCase],
    *,
    process: psutil.Process,
    reuse_symbolic: bool,
) -> tuple[np.ndarray, SequentialStats]:
    years = [int(case.raw_year) for case in cases]
    if trails.A is None:
        raise RuntimeError("trails.A is None")
    n_products = int(trails.A.shape[2])
    n_activities = int(trails.A.shape[1])
    solution_parts: list[np.ndarray] = []

    ctx = UmfpackContext()
    symbolic = None
    pattern_sig: tuple[int, int, int, int] | None = None
    before = _rss(process)

    matrix_build_seconds = 0.0
    symbolic_seconds = 0.0
    numeric_seconds = 0.0
    solve_seconds = 0.0
    max_csc_storage = 0
    sum_csc_storage = 0
    max_symbolic_size: int | None = None
    max_symbolic_peak: int | None = None
    max_numeric_size: int | None = None
    sum_numeric_size = 0
    max_factorization_peak: int | None = None
    max_variable_peak: int | None = None

    with PeakRSSMonitor(process) as monitor:
        for local_t, case in enumerate(cases):
            start = time.perf_counter()
            matrix = _build_year_csc(trails, int(case.time_index))
            matrix_build_seconds += time.perf_counter() - start

            csc_storage = _storage_total(_csc_storage_bytes(matrix))
            max_csc_storage = max(max_csc_storage, csc_storage)
            sum_csc_storage += csc_storage

            rhs = np.zeros(n_products, dtype=np.float64)
            mapped = _activity_demand_to_product_rhs(matrix, case.activity_demand)
            for product, amount in mapped.items():
                rhs[int(product)] += float(amount)

            current_sig = _matrix_signature(matrix)
            if (not reuse_symbolic) or symbolic is None or current_sig != pattern_sig:
                start = time.perf_counter()
                symbolic = ctx.symbolic(matrix)
                symbolic_seconds += time.perf_counter() - start
                pattern_sig = current_sig

            start = time.perf_counter()
            try:
                ctx.numeric(matrix, symbolic)
            except TypeError:
                ctx.numeric(matrix)
            numeric_seconds += time.perf_counter() - start

            start = time.perf_counter()
            solution_parts.append(ctx.solve(UMFPACK_A, matrix, rhs))
            solve_seconds += time.perf_counter() - start

            sym_size = _umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_SIZE")
            sym_peak = _umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_PEAK_MEMORY")
            num_size = _umfpack_info_bytes(ctx, "UMFPACK_NUMERIC_SIZE")
            fact_peak = _umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY")
            variable_peak = _umfpack_info_bytes(ctx, "UMFPACK_VARIABLE_PEAK")

            max_symbolic_size = max(
                max_symbolic_size or 0, sym_size or 0
            ) or None
            max_symbolic_peak = max(
                max_symbolic_peak or 0, sym_peak or 0
            ) or None
            max_numeric_size = max(max_numeric_size or 0, num_size or 0) or None
            sum_numeric_size += int(num_size or 0)
            max_factorization_peak = max(
                max_factorization_peak or 0, fact_peak or 0
            ) or None
            max_variable_peak = max(
                max_variable_peak or 0, variable_peak or 0
            ) or None

            try:
                ctx.free_numeric()
            except Exception:
                pass

            if local_t % 10 == 0:
                print(
                    f"  sequential year {years[local_t]} "
                    f"({local_t + 1}/{len(cases)})",
                    flush=True,
                )

    after = _rss(process)
    try:
        ctx.free_symbolic()
    except Exception:
        pass

    solution = (
        np.concatenate(solution_parts)
        if solution_parts
        else np.zeros(0, dtype=np.float64)
    )
    if solution.size != len(cases) * n_activities:
        raise RuntimeError(
            f"Unexpected sequential solution size {solution.size}; "
            f"expected {len(cases) * n_activities}."
        )

    return solution, SequentialStats(
        years=years,
        matrix_build_seconds=matrix_build_seconds,
        symbolic_seconds=symbolic_seconds,
        numeric_seconds=numeric_seconds,
        solve_seconds=solve_seconds,
        max_csc_storage=max_csc_storage,
        sum_csc_storage=sum_csc_storage,
        max_symbolic_size=max_symbolic_size,
        max_symbolic_peak=max_symbolic_peak,
        max_numeric_size=max_numeric_size,
        sum_numeric_size=sum_numeric_size,
        max_factorization_peak=max_factorization_peak,
        max_variable_peak=max_variable_peak,
        memory=MemoryWindow(before=before, after=after, peak=monitor.peak),
    )


def _run_full_trails_lca(trails: Trails, args: argparse.Namespace) -> None:
    process = psutil.Process(os.getpid())
    before = _rss(process)
    with PeakRSSMonitor(process) as monitor:
        start = time.perf_counter()
        trails.temporal_routing(
            start_year=int(args.start_year),
            start_act_idx=int(args.activity_index),
            amount=float(args.amount),
            max_depth=int(args.max_depth),
            show_progress=bool(args.show_progress),
            attribute_to_roots=bool(args.attribute_to_roots),
        )
        routing_seconds = time.perf_counter() - start

        start = time.perf_counter()
        trails.lca(
            methods=args.method or [DEFAULT_METHOD],
            show_progress=bool(args.show_progress),
            compute_score=not bool(args.no_compute_score),
            store_inventory=not bool(args.no_store_inventory),
            attribute_to_roots=bool(args.attribute_to_roots),
            solver_mode=str(args.solver_mode),
            iterative_rtol=float(args.iterative_rtol),
        )
        lca_seconds = time.perf_counter() - start

    after = _rss(process)
    memory = MemoryWindow(before=before, after=after, peak=monitor.peak)
    print("\nFull Trails temporal_routing + lca")
    print(f"  routing time: {routing_seconds:.6f} s")
    print(f"  lca time: {lca_seconds:.6f} s")
    print(
        "  process RSS: "
        f"after={_format_bytes(memory.after)}, "
        f"delta={_format_bytes(memory.delta)}, "
        f"peak_delta={_format_bytes(memory.peak_delta)}"
    )


def _print_umfpack_stats(prefix: str, stats: UmfpackStats) -> None:
    print(prefix)
    print(f"  UMFPACK symbolic time: {stats.symbolic_seconds:.6f} s")
    print(f"  UMFPACK numeric time: {stats.numeric_seconds:.6f} s")
    print(
        "  UMFPACK solve times: "
        + ", ".join(f"{value:.6f} s" for value in stats.solve_seconds)
    )
    print(
        "  UMFPACK total factorize+first-solve time: "
        f"{stats.first_factorize_solve_seconds:.6f} s"
    )
    print(f"  UMFPACK symbolic object size: {_format_bytes(stats.symbolic_size)}")
    print(f"  UMFPACK symbolic peak memory: {_format_bytes(stats.symbolic_peak)}")
    print(f"  UMFPACK LU numeric object size: {_format_bytes(stats.numeric_size)}")
    print(f"  UMFPACK factorization peak memory: {_format_bytes(stats.peak_memory)}")
    print(f"  UMFPACK variable peak memory: {_format_bytes(stats.variable_peak)}")
    print(
        "  process RSS: "
        f"after={_format_bytes(stats.memory.after)}, "
        f"delta={_format_bytes(stats.memory.delta)}, "
        f"peak_delta={_format_bytes(stats.memory.peak_delta)}"
    )


def _print_sequential_stats(stats: SequentialStats) -> None:
    print("\nSequential year-by-year UMFPACK solves")
    print(f"  years solved: {len(stats.years)}")
    print(f"  matrix build time: {stats.matrix_build_seconds:.6f} s")
    print(f"  UMFPACK symbolic time: {stats.symbolic_seconds:.6f} s")
    print(f"  UMFPACK numeric time: {stats.numeric_seconds:.6f} s")
    print(f"  UMFPACK solve time: {stats.solve_seconds:.6f} s")
    print(
        "  total factorize+solve time: "
        f"{stats.symbolic_seconds + stats.numeric_seconds + stats.solve_seconds:.6f} s"
    )
    print(f"  max per-year CSC storage: {_format_bytes(stats.max_csc_storage)}")
    print(f"  sum per-year CSC storage: {_format_bytes(stats.sum_csc_storage)}")
    print(f"  max symbolic object size: {_format_bytes(stats.max_symbolic_size)}")
    print(f"  max symbolic peak memory: {_format_bytes(stats.max_symbolic_peak)}")
    print(f"  max LU numeric object size: {_format_bytes(stats.max_numeric_size)}")
    print(
        "  sum LU numeric sizes if retained: "
        f"{_format_bytes(stats.sum_numeric_size)}"
    )
    print(
        "  max UMFPACK factorization peak memory: "
        f"{_format_bytes(stats.max_factorization_peak)}"
    )
    print(
        "  max UMFPACK variable peak memory: "
        f"{_format_bytes(stats.max_variable_peak)}"
    )
    print(
        "  process RSS: "
        f"after={_format_bytes(stats.memory.after)}, "
        f"delta={_format_bytes(stats.memory.delta)}, "
        f"peak_delta={_format_bytes(stats.memory.peak_delta)}"
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a block-diagonal sparse matrix from trails.A and compare it "
            "with year-by-year TRAILS-style technosphere solves."
        )
    )
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument(
        "--interpolate-annual",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--years",
        nargs="+",
        type=int,
        default=None,
        help=(
            "Years to include. In routing mode this filters routed frontier "
            "years; in activity-each-year mode these are scenario labels."
        ),
    )
    parser.add_argument(
        "--horizon-start",
        type=int,
        default=None,
        help=(
            "First scenario year to include in the matrix horizon. Ignored when "
            "--years is given."
        ),
    )
    parser.add_argument(
        "--horizon-end",
        type=int,
        default=None,
        help=(
            "Last scenario year to include in the matrix horizon. Ignored when "
            "--years is given."
        ),
    )
    parser.add_argument(
        "--max-years",
        type=int,
        default=None,
        help="Use only the first N selected years.",
    )
    parser.add_argument("--activity-index", type=int, default=31387)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument(
        "--demand-mode",
        choices=("routing", "activity-each-year"),
        default="routing",
        help=(
            "routing uses temporal_routing frontier demands for an exact "
            "monolithic-vs-sequential comparison. activity-each-year applies "
            "the activity demand independently in each selected year."
        ),
    )
    parser.add_argument(
        "--matrix-mode",
        choices=("block-diagonal", "time-expanded"),
        default="time-expanded",
        help=(
            "time-expanded replaces temporally distributed same-year exchanges "
            "with off-diagonal activity-year links. block-diagonal keeps years "
            "independent and is useful for exact routed-frontier comparisons."
        ),
    )
    parser.add_argument(
        "--run-monolithic-factorization",
        action="store_true",
        help=(
            "Factorize and solve the all-year block matrix. This can require "
            "many GiB for all annual years."
        ),
    )
    parser.add_argument(
        "--skip-block-matrix",
        action="store_true",
        help="Skip materializing the block matrix; still run sequential solves.",
    )
    parser.add_argument(
        "--no-sequential-solve",
        action="store_true",
        help="Only build/report the block matrix.",
    )
    parser.add_argument(
        "--no-reuse-symbolic",
        action="store_true",
        help="Do not reuse UMFPACK symbolic factorization across year slices.",
    )
    parser.add_argument("--run-trails-lca", action="store_true")
    parser.add_argument("--start-year", type=int, default=2025)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument(
        "--method",
        action="append",
        default=None,
        help="LCIA method for --run-trails-lca. Can be repeated.",
    )
    parser.add_argument(
        "--solver-mode",
        choices=("bw2calc", "direct", "iterative"),
        default="direct",
    )
    parser.add_argument("--iterative-rtol", type=float, default=1e-3)
    parser.add_argument(
        "--attribute-to-roots",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--no-store-inventory", action="store_true")
    parser.add_argument("--no-compute-score", action="store_true")
    parser.add_argument("--show-progress", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    process = psutil.Process(os.getpid())
    package_path = Path(args.package).expanduser().resolve()

    print(f"PID: {process.pid}")
    print(f"Package: {package_path}")
    print(f"interpolate_annual: {bool(args.interpolate_annual)}")

    start = time.perf_counter()
    before = _rss(process)
    with PeakRSSMonitor(process) as monitor:
        trails = Trails(
            Package(str(package_path)),
            interpolate_annual=bool(args.interpolate_annual),
            debug=False,
        )
    after = _rss(process)
    load_memory = MemoryWindow(before=before, after=after, peak=monitor.peak)
    print(f"Loaded Trails in {time.perf_counter() - start:.6f} s")
    print(
        "Load RSS: "
        f"after={_format_bytes(load_memory.after)}, "
        f"delta={_format_bytes(load_memory.delta)}, "
        f"peak_delta={_format_bytes(load_memory.peak_delta)}"
    )

    if trails.A is None:
        raise RuntimeError("trails.A is None after loading package")
    if trails.A.shape[1] != trails.A.shape[2]:
        raise RuntimeError(f"trails.A slices are not square: {trails.A.shape}")

    A_storage = _coo_storage_bytes(trails.A)
    print("\ntrails.A")
    print(f"  shape: {trails.A.shape}")
    print(f"  nnz: {int(trails.A.nnz):,}")
    print(f"  dtype: {trails.A.dtype}")
    print(
        "  sparse.COO storage: "
        f"total={_format_bytes(_storage_total(A_storage))} "
        f"(data={_format_bytes(A_storage['data'])}, "
        f"coords={_format_bytes(A_storage['coords'])})"
    )

    requested_years = args.years
    if requested_years is None and (
        args.horizon_start is not None or args.horizon_end is not None
    ):
        all_years = [int(label) for label in trails.scenario_labels]
        start_year = (
            min(all_years) if args.horizon_start is None else int(args.horizon_start)
        )
        end_year = max(all_years) if args.horizon_end is None else int(args.horizon_end)
        requested_years = [
            year for year in all_years if start_year <= int(year) <= end_year
        ]
        if not requested_years:
            raise ValueError(
                f"No scenario years found in requested horizon {start_year}..{end_year}"
            )

    selected_years = (
        {int(year) for year in requested_years} if requested_years else None
    )
    selected_time_indices = _selected_time_indices(trails, requested_years)
    if args.max_years is not None:
        selected_time_indices = selected_time_indices[: int(args.max_years)]
    if args.matrix_mode == "time-expanded":
        start_t = _scenario_time_index(trails, int(args.start_year))
        if start_t not in selected_time_indices:
            selected_time_indices = sorted([*selected_time_indices, start_t])

    cases: list[YearCase] = []
    if args.demand_mode == "routing":
        cases = _cases_from_routing(
            trails,
            start_year=int(args.start_year),
            activity_index=int(args.activity_index),
            amount=float(args.amount),
            max_depth=int(args.max_depth),
            selected_years=selected_years,
            show_progress=bool(args.show_progress),
        )
    else:
        cases = _cases_from_activity_each_year(
            trails,
            selected_time_indices,
            activity_index=int(args.activity_index),
            amount=float(args.amount),
        )

    if args.max_years is not None and args.demand_mode == "routing":
        cases = cases[: int(args.max_years)]
    if not cases:
        raise RuntimeError("No year cases selected for benchmark.")

    raw_years = [case.raw_year for case in cases]
    scenario_years = [int(trails.scenario_labels[case.time_index]) for case in cases]
    print("\nBenchmark demand cases")
    print(f"  demand mode: {args.demand_mode}")
    print(f"  count: {len(cases)}")
    print(f"  raw first/last: {raw_years[0]} / {raw_years[-1]}")
    print(f"  scenario first/last: {scenario_years[0]} / {scenario_years[-1]}")
    print(
        "  total activity-demand nonzeros: "
        f"{sum(int(np.count_nonzero(case.activity_demand)) for case in cases):,}"
    )

    matrix_time_indices = (
        selected_time_indices
        if args.matrix_mode == "time-expanded"
        else [case.time_index for case in cases]
    )
    matrix_years = _time_labels(trails, matrix_time_indices)
    print("\nMatrix horizon")
    print(f"  matrix mode: {args.matrix_mode}")
    print(f"  year blocks: {len(matrix_time_indices)}")
    print(f"  first/last: {matrix_years[0]} / {matrix_years[-1]}")

    block_matrix: sp.csc_matrix | None = None
    mono_solution: np.ndarray | None = None
    rhs: np.ndarray | None = None
    if not args.skip_block_matrix:
        start = time.perf_counter()
        before = _rss(process)
        with PeakRSSMonitor(process) as monitor:
            if args.matrix_mode == "time-expanded":
                block_matrix, expanded_stats = _build_time_expanded_csc(
                    trails,
                    matrix_time_indices,
                )
                rhs = _build_single_start_rhs(
                    block_matrix,
                    trails,
                    matrix_time_indices,
                    start_year=int(args.start_year),
                    activity_index=int(args.activity_index),
                    amount=float(args.amount),
                )
            else:
                expanded_stats = None
                block_matrix = _build_block_csc(trails, cases)
                rhs = _build_block_rhs(
                    block_matrix,
                    cases,
                    n_activities=int(trails.A.shape[1]),
                    n_products=int(trails.A.shape[2]),
                )
        after = _rss(process)
        build_memory = MemoryWindow(before=before, after=after, peak=monitor.peak)
        csc_parts = _csc_storage_bytes(block_matrix)
        title = (
            "Time-expanded technosphere matrix"
            if args.matrix_mode == "time-expanded"
            else "Block-diagonal routed-frontier technosphere matrix"
        )
        print(f"\n{title}")
        print(f"  build time: {time.perf_counter() - start:.6f} s")
        print(f"  shape: {block_matrix.shape}")
        print(f"  nnz: {block_matrix.nnz:,}")
        if expanded_stats is not None:
            print(
                "  temporal exchanges selected: "
                f"{expanded_stats.td_exchanges_selected:,}"
            )
            print(
                "  temporal off-diagonal entries added: "
                f"{expanded_stats.td_offsets_added:,}"
            )
            print(
                "  temporal offsets dropped outside horizon: "
                f"{expanded_stats.td_offsets_dropped:,}"
            )
        print(
            "  CSC storage: "
            f"total={_format_bytes(_storage_total(csc_parts))} "
            f"(data={_format_bytes(csc_parts['data'])}, "
            f"indices={_format_bytes(csc_parts['indices'])}, "
            f"indptr={_format_bytes(csc_parts['indptr'])})"
        )
        print(f"  RHS nonzero entries: {int(np.count_nonzero(rhs)):,}")
        print(
            "  process RSS: "
            f"after={_format_bytes(build_memory.after)}, "
            f"delta={_format_bytes(build_memory.delta)}, "
            f"peak_delta={_format_bytes(build_memory.peak_delta)}"
        )

        if args.run_monolithic_factorization:
            try:
                mono_solution, mono_stats = _solve_umfpack(
                    block_matrix,
                    rhs,
                    repeat=max(1, int(args.repeat)),
                    process=process,
                )
            except (MemoryError, RuntimeError) as exc:
                print("\nMonolithic block-matrix UMFPACK solve failed")
                print(f"  failure type: {type(exc).__name__}")
                print(f"  failure message: {exc}")
                print(f"  process RSS at failure: {_format_bytes(_rss(process))}")
            else:
                _print_umfpack_stats(
                    "\nMonolithic block-matrix UMFPACK solve",
                    mono_stats,
                )
                residual = block_matrix @ mono_solution - rhs
                print(
                    "  residual infinity norm: "
                    f"{float(np.linalg.norm(residual, ord=np.inf)):.3e}"
                )
        else:
            print(
                "  monolithic factorization skipped; pass "
                "--run-monolithic-factorization to attempt it."
            )

    if not args.no_sequential_solve:
        if args.matrix_mode == "time-expanded":
            print(
                "\nSequential year-by-year solve skipped for time-expanded "
                "matrix mode. Use --run-trails-lca to time TRAILS' graph + "
                "year-wise workflow, or --matrix-mode block-diagonal for an "
                "exact same-RHS solution comparison."
            )
        else:
            sequential_solution, sequential_stats = _solve_sequential_umfpack(
                trails,
                cases,
                process=process,
                reuse_symbolic=not bool(args.no_reuse_symbolic),
            )
            _print_sequential_stats(sequential_stats)
            if mono_solution is not None:
                diff = sequential_solution - mono_solution
                print("\nMonolithic vs sequential solution comparison")
                print(f"  solution size: {sequential_solution.size:,}")
                print(f"  max abs difference: {float(np.max(np.abs(diff))):.3e}")
                print(f"  L2 difference: {float(np.linalg.norm(diff)):.3e}")
                print(
                    "  sequential solution L1 norm: "
                    f"{float(np.linalg.norm(sequential_solution, ord=1)):.12g}"
                )

    if args.run_trails_lca:
        _run_full_trails_lca(trails, args)


if __name__ == "__main__":
    main()
