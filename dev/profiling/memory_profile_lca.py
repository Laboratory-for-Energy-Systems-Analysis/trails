#!/usr/bin/env python
"""Memory diagnostics for trails.lca() with focus on finalize_inventory()."""

from __future__ import annotations

import argparse
import gc
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import psutil
from datapackage import Package

from trails import Trails

DEFAULT_METHOD = (
    "IPCC 2021 - climate change: total (excl. biogenic CO2) - "
    "global warming potential (GWP100)"
)


def _fmt_bytes(n: int | float) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(n)
    for unit in units:
        if abs(value) < 1024.0 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{value:.2f} TiB"


def _log(msg: str) -> None:
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)


@dataclass
class MemorySample:
    rss: int
    vms: int


class PeakRSSMonitor:
    """Background RSS sampler to capture short-lived peaks."""

    def __init__(self, pid: int, interval_s: float = 0.02) -> None:
        self._process = psutil.Process(pid)
        self._interval_s = float(interval_s)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.peak_rss = 0

    def _loop(self) -> None:
        while not self._stop.is_set():
            try:
                rss = int(self._process.memory_info().rss)
            except psutil.Error:
                break
            if rss > self.peak_rss:
                self.peak_rss = rss
            self._stop.wait(self._interval_s)

    def start(self) -> None:
        if self._thread is not None:
            return
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop.set()
        self._thread.join(timeout=2.0)


def _mem_sample(process: psutil.Process) -> MemorySample:
    info = process.memory_info()
    return MemorySample(rss=int(info.rss), vms=int(info.vms))


def _sum_nbytes(parts: Iterable[Any]) -> int:
    total = 0
    for part in parts:
        total += int(getattr(part, "nbytes", 0))
    return total


def _sum_entries(parts: Iterable[Any]) -> int:
    total = 0
    for part in parts:
        size = getattr(part, "size", None)
        if size is None:
            continue
        total += int(size)
    return total


def _builder_stats(trails: Trails) -> dict[str, int]:
    key_parts = list(getattr(trails, "_inv_key_parts", []) or [])
    value_parts = list(getattr(trails, "_inv_value_parts", []) or [])
    return {
        "key_parts": len(key_parts),
        "value_parts": len(value_parts),
        "key_entries": _sum_entries(key_parts),
        "value_entries": _sum_entries(value_parts),
        "key_nbytes": _sum_nbytes(key_parts),
        "value_nbytes": _sum_nbytes(value_parts),
    }


def _log_builder_stats(trails: Trails, *, label: str) -> None:
    stats = _builder_stats(trails)
    _log(
        (
            f"{label}: _inv_key_parts={stats['key_parts']} "
            f"({_fmt_bytes(stats['key_nbytes'])}, entries={stats['key_entries']}), "
            f"_inv_value_parts={stats['value_parts']} "
            f"({_fmt_bytes(stats['value_nbytes'])}, entries={stats['value_entries']})"
        )
    )


def _log_inventory_stats(trails: Trails) -> None:
    inv = getattr(trails, "inventory", None)
    if inv is None:
        _log("inventory: None")
        return

    arr = inv.data
    shape = getattr(arr, "shape", None)
    nnz = int(getattr(arr, "nnz", 0))
    data_nbytes = int(getattr(getattr(arr, "data", None), "nbytes", 0))
    coords_nbytes = int(getattr(getattr(arr, "coords", None), "nbytes", 0))
    _log(
        (
            f"inventory sparse array: shape={shape}, nnz={nnz}, "
            f"data={_fmt_bytes(data_nbytes)}, coords={_fmt_bytes(coords_nbytes)}"
        )
    )


def _instrument_finalize_inventory(trails: Trails, process: psutil.Process) -> None:
    original = trails.finalize_inventory

    def wrapped_finalize_inventory(*args: Any, **kwargs: Any) -> Any:
        gc.collect()
        before = _mem_sample(process)
        _log(
            f"finalize_inventory:start rss={_fmt_bytes(before.rss)} "
            f"vms={_fmt_bytes(before.vms)}"
        )
        _log_builder_stats(trails, label="builders before finalize")
        monitor = PeakRSSMonitor(pid=process.pid, interval_s=0.01)
        monitor.start()
        t0 = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            dt = time.perf_counter() - t0
            monitor.stop()
            gc.collect()
            after = _mem_sample(process)
            _log(
                (
                    "finalize_inventory:end "
                    f"elapsed={dt:.3f}s "
                    f"rss={_fmt_bytes(after.rss)} "
                    f"(delta={_fmt_bytes(after.rss - before.rss)}), "
                    f"peak_during_finalize={_fmt_bytes(monitor.peak_rss)} "
                    f"(delta_peak={_fmt_bytes(monitor.peak_rss - before.rss)})"
                )
            )
            _log_builder_stats(trails, label="builders after finalize")
            _log_inventory_stats(trails)

    trails.finalize_inventory = wrapped_finalize_inventory  # type: ignore[method-assign]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Memory diagnostics for trails.lca() around finalize_inventory()."
    )
    parser.add_argument(
        "--datapackage",
        default="/Users/romain/GitHub/premise/dev/trails_2026-02-22.zip",
        help="Path to Frictionless datapackage zip.",
    )
    parser.add_argument(
        "--inventory-xlsx",
        default="/Users/romain/GitHub/trails/dev/lci-case-study-daccs_storage_risk.xlsx",
        help="Path to Excel inventory.",
    )
    parser.add_argument("--ref-year", type=int, default=2025, help="Routing start year.")
    parser.add_argument(
        "--act-idx", type=int, default=41792, help="Routing start activity index."
    )
    parser.add_argument("--amount", type=float, default=1.0, help="Functional unit amount.")
    parser.add_argument("--max-depth", type=int, default=4, help="Routing max depth.")
    parser.add_argument(
        "--method",
        action="append",
        default=None,
        help=(
            "LCIA method name; may be passed multiple times. "
            f"Default: {DEFAULT_METHOD}"
        ),
    )
    parser.add_argument(
        "--solver-mode",
        default="iterative",
        choices=("bw2calc", "direct", "iterative"),
        help="Solver mode passed to trails.lca().",
    )
    parser.add_argument(
        "--iterative-rtol",
        type=float,
        default=1e-4,
        help="iterative_rtol passed to trails.lca().",
    )
    parser.add_argument(
        "--attribute-to-roots",
        action="store_true",
        help="Enable root attribution in temporal routing and lca.",
    )
    parser.add_argument(
        "--show-progress",
        action="store_true",
        help="Show progress bars during routing and lca.",
    )
    parser.add_argument(
        "--no-store-inventory",
        action="store_true",
        help="Disable store_inventory to compare behavior.",
    )
    parser.add_argument(
        "--no-compute-score",
        action="store_true",
        help="Disable score computation (inventory only).",
    )
    parser.add_argument(
        "--repeat-lca",
        type=int,
        default=1,
        help="Number of repeated lca() runs in the same process.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    methods = args.method or [DEFAULT_METHOD]
    store_inventory = not bool(args.no_store_inventory)
    compute_score = not bool(args.no_compute_score)

    process = psutil.Process(os.getpid())
    _log(f"PID={process.pid}")
    _log(
        f"input datapackage={args.datapackage}, inventory_xlsx={args.inventory_xlsx}, "
        f"start_year={args.ref_year}, act_idx={args.act_idx}, "
        f"max_depth={args.max_depth}, amount={args.amount}"
    )

    t0 = time.perf_counter()
    dp = Package(str(Path(args.datapackage).expanduser().resolve()))
    trails = Trails(
        dp,
        interpolate_annual=True,
        debug=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )
    _log(f"initialized Trails in {time.perf_counter() - t0:.3f}s")

    _instrument_finalize_inventory(trails, process)
    _log("instrumented Trails.finalize_inventory()")

    t1 = time.perf_counter()
    trails.import_excel_inventory(str(Path(args.inventory_xlsx).expanduser().resolve()))
    _log(f"import_excel_inventory done in {time.perf_counter() - t1:.3f}s")

    t2 = time.perf_counter()
    trails.temporal_routing(
        start_year=int(args.ref_year),
        start_act_idx=int(args.act_idx),
        amount=float(args.amount),
        max_depth=int(args.max_depth),
        show_progress=bool(args.show_progress),
        attribute_to_roots=bool(args.attribute_to_roots),
    )
    _log(f"temporal_routing done in {time.perf_counter() - t2:.3f}s")
    _log(f"rss before lca={_fmt_bytes(_mem_sample(process).rss)}")

    repeats = max(1, int(args.repeat_lca))
    for run_idx in range(1, repeats + 1):
        _log(f"lca run {run_idx}/{repeats}: start")
        _log(f"rss pre-run={_fmt_bytes(_mem_sample(process).rss)}")
        t3 = time.perf_counter()
        trails.lca(
            methods=methods,
            show_progress=bool(args.show_progress),
            compute_score=bool(compute_score),
            store_inventory=bool(store_inventory),
            attribute_to_roots=bool(args.attribute_to_roots),
            solver_mode=str(args.solver_mode),
            iterative_rtol=float(args.iterative_rtol),
        )
        lca_dt = time.perf_counter() - t3
        _log(f"lca run {run_idx}/{repeats}: done in {lca_dt:.3f}s")
        _log(f"rss post-run={_fmt_bytes(_mem_sample(process).rss)}")

        if compute_score and getattr(trails, "scores", None) is not None:
            scores = trails.scores.data
            nnz = int(getattr(scores, "nnz", 0))
            data_nbytes = int(getattr(getattr(scores, "data", None), "nbytes", 0))
            coords_nbytes = int(getattr(getattr(scores, "coords", None), "nbytes", 0))
            _log(
                (
                    f"scores sparse array: shape={getattr(scores, 'shape', None)}, "
                    f"nnz={nnz}, data={_fmt_bytes(data_nbytes)}, "
                    f"coords={_fmt_bytes(coords_nbytes)}"
                )
            )
        gc.collect()
        _log(f"rss post-run gc={_fmt_bytes(_mem_sample(process).rss)}")

    gc.collect()
    _log(f"rss after final gc={_fmt_bytes(_mem_sample(process).rss)}")


if __name__ == "__main__":
    main()
