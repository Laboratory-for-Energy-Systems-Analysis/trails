#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import gc
import os
import re
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import psutil
from scipy import sparse as sp

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

if TYPE_CHECKING:
    from trails import Trails

DEFAULT_PACKAGE = SCRIPT_DIR / "trails_remind_SSP2-PkBudg1000.zip"
DEFAULT_OUTPUT = SCRIPT_DIR / "monolithic_vs_trails_remind_SSP2-PkBudg1000.csv"
DEFAULT_METHOD = (
    "IPCC 2021 (incl. biogenic CO2) - climate change: total "
    "(incl. biogenic CO2) - global warming potential (GWP100)"
)
UMFPACK_A: Any | None = None
UmfpackContext: Any | None = None
umfpack: Any | None = None
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
    "routing_nodes_processed",
    "routing_max_processed_depth",
    "adaptive_relative_score_cutoff",
    "adaptive_methods",
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


def _require_umfpack() -> tuple[Any, Any]:
    global UMFPACK_A, UmfpackContext, umfpack
    if UmfpackContext is None:
        try:
            from scikits.umfpack import UMFPACK_A as loaded_umfpack_a
            from scikits.umfpack import UmfpackContext as loaded_context
            import scikits.umfpack as loaded_umfpack
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "scikit-umfpack is required for this benchmark. Run it in the "
                "`trails` conda environment or install scikit-umfpack first."
            ) from exc

        UMFPACK_A = loaded_umfpack_a
        UmfpackContext = loaded_context
        umfpack = loaded_umfpack
    return UMFPACK_A, UmfpackContext


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


def _umfpack_info_bytes(ctx: Any, name: str) -> int | None:
    if umfpack is None:
        _require_umfpack()
    if umfpack is None:
        return None
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


def _scenario_time_index(trails: Trails, raw_year: int) -> int:
    context = trails._get_scenario_context(int(raw_year))
    if context is None:
        raise RuntimeError(f"No scenario context available for year={raw_year}")
    _scenario_year, _label, time_index = context
    return int(time_index)


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
        rhs[int(product)] = (
            rhs.get(int(product), 0.0) + float(activity_demand[int(activity)]) * sign
        )
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
    umfpack_a, context_cls = _require_umfpack()
    ctx = context_cls()
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
                peak_estimate = _umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY_ESTIMATE")
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
            solution = ctx.solve(umfpack_a, matrix, rhs)
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

    umfpack_a, context_cls = _require_umfpack()
    ctx = context_cls()
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
            solution_parts.append(ctx.solve(umfpack_a, matrix, rhs))
            solve_seconds += time.perf_counter() - start

            sym_size = _umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_SIZE")
            sym_peak = _umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_PEAK_MEMORY")
            num_size = _umfpack_info_bytes(ctx, "UMFPACK_NUMERIC_SIZE")
            fact_peak = _umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY")
            variable_peak = _umfpack_info_bytes(ctx, "UMFPACK_VARIABLE_PEAK")

            max_symbolic_size = max(max_symbolic_size or 0, sym_size or 0) or None
            max_symbolic_peak = max(max_symbolic_peak or 0, sym_peak or 0) or None
            max_numeric_size = max(max_numeric_size or 0, num_size or 0) or None
            sum_numeric_size += int(num_size or 0)
            max_factorization_peak = (
                max(max_factorization_peak or 0, fact_peak or 0) or None
            )
            max_variable_peak = max(max_variable_peak or 0, variable_peak or 0) or None

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


def _blank_row(run_id: str, phase: str, status: str) -> dict[str, Any]:
    return {key: "" for key in FIELDNAMES} | {
        "run_id": run_id,
        "phase": phase,
        "status": status,
    }


def _memory_fields(memory: MemoryWindow) -> dict[str, Any]:
    return {
        "rss_before_bytes": memory.before,
        "rss_after_bytes": memory.after,
        "rss_delta_bytes": memory.delta,
        "rss_peak_bytes": memory.peak,
        "rss_peak_delta_bytes": memory.peak_delta,
    }


def _routing_depth_field(max_depth: int | None) -> str | int:
    return "" if max_depth is None else int(max_depth)


def _adaptive_methods_field(methods: list[str] | None) -> str:
    return "|".join(methods or [])


def _adaptive_cutoff_field(cutoff: float | None) -> str | float:
    return "" if cutoff is None else float(cutoff)


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
    from datapackage import Package
    from trails import Trails

    return Trails(Package(str(package_path)), interpolate_annual=True, debug=False)


def _cases_from_routing_with_stats(
    *,
    trails: Trails,
    start_year: int,
    activity_index: int,
    amount: float,
    max_depth: int | None,
    adaptive_relative_score_cutoff: float | None,
    adaptive_methods: list[str] | None,
    selected_years: set[int],
) -> tuple[list[YearCase], float, int, int, int | str, int | str]:
    start = time.perf_counter()
    trails.temporal_routing(
        start_year=int(start_year),
        start_act_idx=int(activity_index),
        amount=float(amount),
        max_depth=max_depth,
        adaptive_relative_score_cutoff=adaptive_relative_score_cutoff,
        adaptive_methods=adaptive_methods,
        show_progress=False,
        attribute_to_roots=False,
    )
    routing_seconds = time.perf_counter() - start
    routing_params = getattr(trails, "_routing_params", {}) or {}

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
    cases: list[YearCase] = []
    for raw_year in sorted(frontier_by_year):
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
        raise RuntimeError("Temporal routing produced no selected frontier demands.")

    return (
        cases,
        routing_seconds,
        int(graph.number_of_nodes()),
        int(graph.number_of_edges()),
        routing_params.get("nodes_processed", ""),
        routing_params.get("max_processed_depth", ""),
    )


def _run_paired_horizon(
    *,
    trails: Trails,
    process: psutil.Process,
    horizon_start: int,
    horizon_end: int,
    activity_index: int,
    amount: float,
    max_depth: int | None,
    adaptive_relative_score_cutoff: float | None,
    adaptive_methods: list[str] | None,
    rows: list[dict[str, Any]],
) -> bool:
    run_id = f"{horizon_start}-{horizon_end}"
    years = list(range(int(horizon_start), int(horizon_end) + 1))
    time_indices = _selected_time_indices(trails, years)

    print(f"\nPaired horizon {run_id} ({len(years)} blocks)", flush=True)

    before = _rss(process)
    start = time.perf_counter()
    with PeakRSSMonitor(process) as monitor:
        matrix, expanded_stats = _build_time_expanded_csc(trails, time_indices)
        rhs = _build_single_start_rhs(
            matrix,
            trails,
            time_indices,
            start_year=int(horizon_end),
            activity_index=int(activity_index),
            amount=float(amount),
        )
    matrix_build_seconds = time.perf_counter() - start
    after = _rss(process)
    build_memory = MemoryWindow(before=before, after=after, peak=monitor.peak)
    csc_storage = _storage_total(_csc_storage_bytes(matrix))

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
        f"storage={_format_bytes(csc_storage)} in {matrix_build_seconds:.3f}s",
        flush=True,
    )

    monolithic_ok = True
    monolithic_solution = None
    try:
        monolithic_solution, solve_stats = _solve_umfpack(
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

    (
        cases,
        routing_seconds,
        graph_nodes,
        graph_edges,
        routing_nodes_processed,
        routing_max_processed_depth,
    ) = _cases_from_routing_with_stats(
        trails=trails,
        start_year=int(horizon_end),
        activity_index=int(activity_index),
        amount=float(amount),
        max_depth=max_depth,
        adaptive_relative_score_cutoff=adaptive_relative_score_cutoff,
        adaptive_methods=adaptive_methods,
        selected_years=set(years),
    )
    print(
        "  TRAILS routing finished in "
        f"{routing_seconds:.3f}s over {graph_nodes:,} nodes",
        flush=True,
    )

    sequential_solution, sequential_stats = _solve_sequential_umfpack(
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
            "routing_max_depth": _routing_depth_field(max_depth),
            "routing_nodes_processed": routing_nodes_processed,
            "routing_max_processed_depth": routing_max_processed_depth,
            "adaptive_relative_score_cutoff": _adaptive_cutoff_field(
                adaptive_relative_score_cutoff
            ),
            "adaptive_methods": _adaptive_methods_field(adaptive_methods),
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
    max_depth: int | None,
    adaptive_relative_score_cutoff: float | None,
    adaptive_methods: list[str] | None,
    rows: list[dict[str, Any]],
) -> None:
    print("\nFinal full TRAILS run", flush=True)
    full_time_indices = list(range(len(trails.scenario_labels)))
    temporal_stats = _time_expanded_temporal_stats(trails, full_time_indices)

    before = _rss(process)
    routing_seconds: float | str = ""
    lca_seconds: float | str = ""
    status = "ok"
    failure_type = ""
    failure_message = ""
    with PeakRSSMonitor(process) as monitor:
        try:
            start = time.perf_counter()
            trails.temporal_routing(
                start_year=int(start_year),
                start_act_idx=int(activity_index),
                amount=float(amount),
                max_depth=max_depth,
                adaptive_relative_score_cutoff=adaptive_relative_score_cutoff,
                adaptive_methods=adaptive_methods,
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
    memory = MemoryWindow(before=before, after=after, peak=monitor.peak)
    graph = getattr(trails, "graph", None)
    routing_params = getattr(trails, "_routing_params", {}) or {}
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
            "temporal_exchanges_selected": (temporal_stats.td_exchanges_selected),
            "temporal_offdiag_entries_added": temporal_stats.td_offsets_added,
            "temporal_offsets_dropped": temporal_stats.td_offsets_dropped,
            "csc_storage_bytes": _storage_total(_coo_storage_bytes(trails.A)),
            "routing_seconds": routing_seconds,
            "lca_seconds": lca_seconds,
            "routing_max_depth": _routing_depth_field(max_depth),
            "routing_nodes_processed": routing_params.get("nodes_processed", ""),
            "routing_max_processed_depth": routing_params.get(
                "max_processed_depth", ""
            ),
            "adaptive_relative_score_cutoff": _adaptive_cutoff_field(
                adaptive_relative_score_cutoff
            ),
            "adaptive_methods": _adaptive_methods_field(adaptive_methods),
            "graph_nodes": graph.number_of_nodes() if graph is not None else "",
            "graph_edges": graph.number_of_edges() if graph is not None else "",
            "failure_type": failure_type,
            "failure_message": failure_message,
        }
        | _memory_fields(memory)
    )
    if status == "ok":
        print(
            f"  full TRAILS routing={routing_seconds:.3f}s, " f"lca={lca_seconds:.3f}s",
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
        "--no-max-depth",
        action="store_true",
        help="Run adaptive routing without a hard depth cap.",
    )
    parser.add_argument(
        "--adaptive-relative-score-cutoff",
        type=float,
        default=None,
        help=(
            "Enable adaptive score-potential routing with this relative cutoff. "
            "If no --adaptive-method is given, the benchmark's default GWP100 "
            "method is used."
        ),
    )
    parser.add_argument(
        "--adaptive-method",
        action="append",
        default=None,
        help="Regular LCIA method used for adaptive routing. Can be repeated.",
    )
    parser.add_argument(
        "--skip-paired",
        action="store_true",
        help="Skip the monolithic threshold pairs and run only the full TRAILS case.",
    )
    parser.add_argument(
        "--horizon-sizes",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4],
        help="Contiguous horizon sizes ending at --start-year.",
    )
    return parser.parse_args()


def run_monolithic_vs_trails_benchmark(
    *,
    package: Path | str = DEFAULT_PACKAGE,
    output: Path | str = DEFAULT_OUTPUT,
    activity_index: int = 31387,
    start_year: int = 2025,
    amount: float = 1.0,
    max_depth: int | None = 4,
    adaptive_relative_score_cutoff: float | None = None,
    adaptive_methods: list[str] | None = None,
    skip_paired: bool = False,
    horizon_sizes: tuple[int, ...] | list[int] = (1, 2, 3, 4),
) -> list[dict[str, Any]]:
    """Run the monolithic-vs-TRAILS benchmark and write the result CSV."""
    process = psutil.Process(os.getpid())
    rows: list[dict[str, Any]] = []
    package_path = Path(package).expanduser().resolve()
    output_path = Path(output).expanduser().resolve()
    if adaptive_methods:
        resolved_adaptive_methods = list(adaptive_methods)
    elif adaptive_relative_score_cutoff is not None:
        resolved_adaptive_methods = [DEFAULT_METHOD]
    else:
        resolved_adaptive_methods = None

    print(f"PID: {process.pid}", flush=True)
    print(f"Package: {package_path}", flush=True)
    print(f"Output CSV: {output_path}", flush=True)

    load_before = _rss(process)
    start = time.perf_counter()
    with PeakRSSMonitor(process) as monitor:
        trails = _load_trails(package_path)
    load_after = _rss(process)
    load_memory = MemoryWindow(
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
            "csc_storage_bytes": _storage_total(_coo_storage_bytes(trails.A)),
        }
    )
    row.update(_memory_fields(load_memory))
    rows.append(row)

    if skip_paired:
        print("Skipping paired monolithic threshold runs.", flush=True)
    else:
        first_failure_seen = False
        for size in horizon_sizes:
            horizon_end = int(start_year)
            horizon_start = horizon_end - int(size) + 1
            ok = _run_paired_horizon(
                trails=trails,
                process=process,
                horizon_start=horizon_start,
                horizon_end=horizon_end,
                activity_index=int(activity_index),
                amount=float(amount),
                max_depth=max_depth,
                adaptive_relative_score_cutoff=adaptive_relative_score_cutoff,
                adaptive_methods=resolved_adaptive_methods,
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
        start_year=int(start_year),
        activity_index=int(activity_index),
        amount=float(amount),
        max_depth=max_depth,
        adaptive_relative_score_cutoff=adaptive_relative_score_cutoff,
        adaptive_methods=resolved_adaptive_methods,
        rows=rows,
    )
    _write_rows(output_path, rows)
    print(f"\nWrote CSV: {output_path}", flush=True)
    return rows


def main() -> None:
    args = _parse_args()
    max_depth = None if args.no_max_depth else int(args.max_depth)
    run_monolithic_vs_trails_benchmark(
        package=args.package,
        output=args.output,
        activity_index=int(args.activity_index),
        start_year=int(args.start_year),
        amount=float(args.amount),
        max_depth=max_depth,
        adaptive_relative_score_cutoff=args.adaptive_relative_score_cutoff,
        adaptive_methods=args.adaptive_method,
        skip_paired=bool(args.skip_paired),
        horizon_sizes=tuple(int(size) for size in args.horizon_sizes),
    )


if __name__ == "__main__":
    main()
