#!/usr/bin/env python
"""Benchmark finalize_inventory variants for speed/memory/score parity."""

from __future__ import annotations

import argparse
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import psutil
import sparse
import xarray as xr
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
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


@dataclass
class VariantResult:
    variant: str
    setup_seconds: float
    lca_seconds: float
    finalize_seconds: float
    rss_before_lca: int
    rss_after_lca: int
    peak_rss_lca: int
    peak_rss_finalize: int
    score_sum: float
    inventory_nnz: int
    score_nnz: int
    builders_after: tuple[int, int]


class PeakRSSMonitor:
    """Background RSS sampler."""

    def __init__(self, pid: int, interval_s: float = 0.02) -> None:
        self._process = psutil.Process(pid)
        self._interval_s = float(interval_s)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.peak_rss = 0

    def _run(self) -> None:
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
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop.set()
        self._thread.join(timeout=2.0)


def _build_dataarray_from_inv(
    self: Trails,
    inv: sparse.COO,
    years: np.ndarray,
    n_activities: int,
    n_flows: int,
    has_root: bool,
) -> xr.DataArray:
    dims = ("activity", "flow", "year")
    coords_xr: dict[str, np.ndarray] = {
        "activity": np.arange(n_activities, dtype=int),
        "flow": np.arange(n_flows, dtype=int),
        "year": years,
    }
    if has_root:
        dims = ("activity", "flow", "year", "root activity")
        coords_xr["root activity"] = np.arange(n_activities, dtype=int)

    self.inventory = xr.DataArray(inv, dims=dims, coords=coords_xr)
    self._inv_key_parts = []
    self._inv_value_parts = []
    return self.inventory


def _finalize_legacy_vstack(self: Trails) -> xr.DataArray:
    if self.A is None or self.B is None:
        raise ValueError("Cannot finalize inventory: A or B is None.")

    years = self._inventory_years
    if years is None:
        raise RuntimeError("Inventory years not initialized. Call reset_inventory().")

    if not hasattr(self, "_inv_key_parts"):
        raise RuntimeError(
            "Inventory builders not initialized. Call reset_inventory()."
        )

    n_activities = int(self.A.shape[1])
    n_flows = int(self.B.shape[2])
    has_root = bool(getattr(self, "_inventory_has_root", False))

    if self._inv_key_parts:
        keys = np.concatenate(self._inv_key_parts).astype(np.int64, copy=False)
        data = np.concatenate(self._inv_value_parts).astype(
            self.value_dtype, copy=False
        )

        if keys.size:
            order = np.argsort(keys, kind="quicksort")
            keys_sorted = keys[order]
            data_sorted = data[order]

            first = np.empty(keys_sorted.size, dtype=bool)
            first[0] = True
            first[1:] = keys_sorted[1:] != keys_sorted[:-1]
            group_starts = np.flatnonzero(first)
            data_agg = np.add.reduceat(data_sorted, group_starts).astype(
                self.value_dtype, copy=False
            )
            keys_agg = keys_sorted[group_starts]

            keep = data_agg != 0.0
            keys_agg = keys_agg[keep]
            data_agg = data_agg[keep]
        else:
            keys_agg = np.empty(0, dtype=np.int64)
            data_agg = np.empty(0, dtype=self.value_dtype)

        n_years = int(len(years))
        if has_root:
            shape = (n_activities, n_flows, n_years, n_activities)
        else:
            shape = (n_activities, n_flows, n_years)

        coord_dtype = np.int64
        if keys_agg.size:
            if has_root:
                q, root = np.divmod(keys_agg, n_activities)
            else:
                root = None
                q = keys_agg

            q, year_idx = np.divmod(q, n_years)
            act_idx, flow_idx = np.divmod(q, n_flows)

            if has_root and root is not None:
                coords = np.vstack([act_idx, flow_idx, year_idx, root]).astype(
                    coord_dtype, copy=False
                )
            else:
                coords = np.vstack([act_idx, flow_idx, year_idx]).astype(
                    coord_dtype, copy=False
                )
            inv = sparse.COO(
                coords,
                data_agg,
                shape=shape,
                has_duplicates=False,
                sorted=True,
                idx_dtype=coord_dtype,
            )
        else:
            inv = sparse.zeros(shape, dtype=self.value_dtype)
    else:
        if has_root:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years), n_activities),
                dtype=self.value_dtype,
            )
        else:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years)), dtype=self.value_dtype
            )

    return _build_dataarray_from_inv(self, inv, years, n_activities, n_flows, has_root)


def _finalize_legacy_inplace_decode(self: Trails) -> xr.DataArray:
    if self.A is None or self.B is None:
        raise ValueError("Cannot finalize inventory: A or B is None.")

    years = self._inventory_years
    if years is None:
        raise RuntimeError("Inventory years not initialized. Call reset_inventory().")

    if not hasattr(self, "_inv_key_parts"):
        raise RuntimeError(
            "Inventory builders not initialized. Call reset_inventory()."
        )

    n_activities = int(self.A.shape[1])
    n_flows = int(self.B.shape[2])
    has_root = bool(getattr(self, "_inventory_has_root", False))

    if self._inv_key_parts:
        keys = np.concatenate(self._inv_key_parts).astype(np.int64, copy=False)
        data = np.concatenate(self._inv_value_parts).astype(
            self.value_dtype, copy=False
        )

        if keys.size:
            order = np.argsort(keys, kind="quicksort")
            keys_sorted = keys[order]
            data_sorted = data[order]

            first = np.empty(keys_sorted.size, dtype=bool)
            first[0] = True
            first[1:] = keys_sorted[1:] != keys_sorted[:-1]
            group_starts = np.flatnonzero(first)
            data_agg = np.add.reduceat(data_sorted, group_starts).astype(
                self.value_dtype, copy=False
            )
            keys_agg = keys_sorted[group_starts].copy()

            keep = data_agg != 0.0
            keys_agg = keys_agg[keep]
            data_agg = data_agg[keep]
        else:
            keys_agg = np.empty(0, dtype=np.int64)
            data_agg = np.empty(0, dtype=self.value_dtype)

        n_years = int(len(years))
        if has_root:
            shape = (n_activities, n_flows, n_years, n_activities)
        else:
            shape = (n_activities, n_flows, n_years)

        coord_dtype = np.int64
        if keys_agg.size:
            if has_root:
                coords = np.empty((4, keys_agg.size), dtype=coord_dtype)
                q = keys_agg
                np.remainder(q, n_activities, out=coords[3])
                np.floor_divide(q, n_activities, out=q)
                np.remainder(q, n_years, out=coords[2])
                np.floor_divide(q, n_years, out=q)
                np.remainder(q, n_flows, out=coords[1])
                np.floor_divide(q, n_flows, out=coords[0])
            else:
                coords = np.empty((3, keys_agg.size), dtype=coord_dtype)
                q = keys_agg
                np.remainder(q, n_years, out=coords[2])
                np.floor_divide(q, n_years, out=q)
                np.remainder(q, n_flows, out=coords[1])
                np.floor_divide(q, n_flows, out=coords[0])
            inv = sparse.COO(
                coords,
                data_agg,
                shape=shape,
                has_duplicates=False,
                sorted=True,
                idx_dtype=coord_dtype,
            )
        else:
            inv = sparse.zeros(shape, dtype=self.value_dtype)
    else:
        if has_root:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years), n_activities),
                dtype=self.value_dtype,
            )
        else:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years)), dtype=self.value_dtype
            )

    return _build_dataarray_from_inv(self, inv, years, n_activities, n_flows, has_root)


def _finalize_legacy_reassign(self: Trails) -> xr.DataArray:
    """Fast path using argsort + reassignment to drop intermediates sooner."""
    if self.A is None or self.B is None:
        raise ValueError("Cannot finalize inventory: A or B is None.")

    years = self._inventory_years
    if years is None:
        raise RuntimeError("Inventory years not initialized. Call reset_inventory().")

    if not hasattr(self, "_inv_key_parts"):
        raise RuntimeError(
            "Inventory builders not initialized. Call reset_inventory()."
        )

    n_activities = int(self.A.shape[1])
    n_flows = int(self.B.shape[2])
    has_root = bool(getattr(self, "_inventory_has_root", False))

    if self._inv_key_parts:
        keys = np.concatenate(self._inv_key_parts).astype(np.int64, copy=False)
        data = np.concatenate(self._inv_value_parts).astype(
            self.value_dtype, copy=False
        )

        if keys.size:
            order = np.argsort(keys, kind="quicksort")
            keys = keys[order]
            data = data[order]
            del order

            first = np.empty(keys.size, dtype=bool)
            first[0] = True
            first[1:] = keys[1:] != keys[:-1]
            group_starts = np.flatnonzero(first)
            data_agg = np.add.reduceat(data, group_starts).astype(
                self.value_dtype, copy=False
            )
            keys_agg = keys[group_starts].copy()

            keep = data_agg != 0.0
            keys_agg = keys_agg[keep]
            data_agg = data_agg[keep]
        else:
            keys_agg = np.empty(0, dtype=np.int64)
            data_agg = np.empty(0, dtype=self.value_dtype)

        n_years = int(len(years))
        if has_root:
            shape = (n_activities, n_flows, n_years, n_activities)
        else:
            shape = (n_activities, n_flows, n_years)

        coord_dtype = np.int64
        if keys_agg.size:
            if has_root:
                coords = np.empty((4, keys_agg.size), dtype=coord_dtype)
                q = keys_agg
                np.remainder(q, n_activities, out=coords[3])
                np.floor_divide(q, n_activities, out=q)
                np.remainder(q, n_years, out=coords[2])
                np.floor_divide(q, n_years, out=q)
                np.remainder(q, n_flows, out=coords[1])
                np.floor_divide(q, n_flows, out=coords[0])
            else:
                coords = np.empty((3, keys_agg.size), dtype=coord_dtype)
                q = keys_agg
                np.remainder(q, n_years, out=coords[2])
                np.floor_divide(q, n_years, out=q)
                np.remainder(q, n_flows, out=coords[1])
                np.floor_divide(q, n_flows, out=coords[0])
            inv = sparse.COO(
                coords,
                data_agg,
                shape=shape,
                has_duplicates=False,
                sorted=True,
                idx_dtype=coord_dtype,
            )
        else:
            inv = sparse.zeros(shape, dtype=self.value_dtype)
    else:
        if has_root:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years), n_activities),
                dtype=self.value_dtype,
            )
        else:
            inv = sparse.zeros(
                (n_activities, n_flows, len(years)), dtype=self.value_dtype
            )

    return _build_dataarray_from_inv(self, inv, years, n_activities, n_flows, has_root)


def _get_variant_finalize(
    trails_obj: Trails, name: str, current_finalize: Callable[[], xr.DataArray]
) -> Callable[[], xr.DataArray]:
    if name == "current_compact":
        return current_finalize
    if name == "legacy_vstack":
        return lambda: _finalize_legacy_vstack(trails_obj)
    if name == "legacy_inplace":
        return lambda: _finalize_legacy_inplace_decode(trails_obj)
    if name == "legacy_reassign":
        return lambda: _finalize_legacy_reassign(trails_obj)
    raise ValueError(f"Unknown variant: {name}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark finalize_inventory implementations."
    )
    parser.add_argument(
        "--datapackage",
        default="/Users/romain/GitHub/premise/dev/trails_2026-02-22.zip",
    )
    parser.add_argument(
        "--inventory-xlsx",
        default="/Users/romain/GitHub/trails/dev/lci-case-study-daccs_storage_risk.xlsx",
    )
    parser.add_argument("--ref-year", type=int, default=2025)
    parser.add_argument("--act-idx", type=int, default=41792)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument(
        "--variant",
        action="append",
        default=None,
        help=(
            "Variant to benchmark. Can be passed multiple times. "
            "Choices: current_compact, legacy_vstack, legacy_inplace, legacy_reassign"
        ),
    )
    parser.add_argument("--method", action="append", default=None)
    parser.add_argument(
        "--solver-mode",
        default="iterative",
        choices=("bw2calc", "direct", "iterative"),
    )
    parser.add_argument("--iterative-rtol", type=float, default=1e-4)
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--attribute-to-roots", action="store_true")
    return parser.parse_args()


def run_variant(
    dp: Package,
    process: psutil.Process,
    variant: str,
    *,
    methods: list[str],
    inventory_xlsx: str,
    ref_year: int,
    act_idx: int,
    amount: float,
    max_depth: int,
    show_progress: bool,
    attribute_to_roots: bool,
    solver_mode: str,
    iterative_rtol: float,
) -> VariantResult:
    _log(f"running variant={variant}")
    setup_t0 = time.perf_counter()
    trails_obj = Trails(
        dp,
        interpolate_annual=True,
        debug=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )
    trails_obj.import_excel_inventory(inventory_xlsx)
    trails_obj.temporal_routing(
        start_year=int(ref_year),
        start_act_idx=int(act_idx),
        amount=float(amount),
        max_depth=int(max_depth),
        show_progress=bool(show_progress),
        attribute_to_roots=bool(attribute_to_roots),
    )
    setup_seconds = time.perf_counter() - setup_t0

    rss_before_lca = int(process.memory_info().rss)
    lca_monitor = PeakRSSMonitor(process.pid, interval_s=0.01)
    lca_monitor.start()

    finalize_elapsed = 0.0
    finalize_peak = rss_before_lca
    finalize_base = rss_before_lca
    current_finalize = trails_obj.finalize_inventory
    variant_finalize = _get_variant_finalize(trails_obj, variant, current_finalize)

    def wrapped_finalize() -> xr.DataArray:
        nonlocal finalize_elapsed, finalize_peak, finalize_base
        finalize_base = int(process.memory_info().rss)
        fin_mon = PeakRSSMonitor(process.pid, interval_s=0.005)
        fin_mon.start()
        t0 = time.perf_counter()
        try:
            return variant_finalize()
        finally:
            finalize_elapsed = time.perf_counter() - t0
            fin_mon.stop()
            finalize_peak = fin_mon.peak_rss

    trails_obj.finalize_inventory = wrapped_finalize  # type: ignore[method-assign]

    lca_t0 = time.perf_counter()
    trails_obj.lca(
        methods=methods,
        show_progress=bool(show_progress),
        compute_score=True,
        store_inventory=True,
        attribute_to_roots=bool(attribute_to_roots),
        solver_mode=solver_mode,
        iterative_rtol=float(iterative_rtol),
    )
    lca_seconds = time.perf_counter() - lca_t0

    lca_monitor.stop()
    rss_after_lca = int(process.memory_info().rss)

    if trails_obj.scores is None:
        raise RuntimeError("scores is None after lca run")
    score_sum = float(trails_obj.scores.sum().item())

    inv_nnz = int(getattr(trails_obj.inventory.data, "nnz", 0))
    score_nnz = int(getattr(trails_obj.scores.data, "nnz", 0))
    builders_after = (
        len(getattr(trails_obj, "_inv_key_parts", []) or []),
        len(getattr(trails_obj, "_inv_value_parts", []) or []),
    )

    _log(
        (
            f"variant={variant} done: setup={setup_seconds:.3f}s lca={lca_seconds:.3f}s "
            f"finalize={finalize_elapsed:.3f}s peak_lca={_fmt_bytes(lca_monitor.peak_rss)} "
            f"peak_finalize={_fmt_bytes(finalize_peak)} score={score_sum:.12f} "
            f"inv_nnz={inv_nnz} score_nnz={score_nnz}"
        )
    )

    return VariantResult(
        variant=variant,
        setup_seconds=setup_seconds,
        lca_seconds=lca_seconds,
        finalize_seconds=finalize_elapsed,
        rss_before_lca=rss_before_lca,
        rss_after_lca=rss_after_lca,
        peak_rss_lca=lca_monitor.peak_rss,
        peak_rss_finalize=finalize_peak,
        score_sum=score_sum,
        inventory_nnz=inv_nnz,
        score_nnz=score_nnz,
        builders_after=builders_after,
    )


def main() -> None:
    args = parse_args()
    variants = args.variant or [
        "current_compact",
        "legacy_vstack",
        "legacy_inplace",
        "legacy_reassign",
    ]
    methods = args.method or [DEFAULT_METHOD]
    process = psutil.Process(os.getpid())

    _log(
        (
            f"benchmark start variants={variants} depth={args.max_depth} "
            f"solver={args.solver_mode} attr_roots={bool(args.attribute_to_roots)}"
        )
    )
    dp = Package(str(Path(args.datapackage).expanduser().resolve()))
    inv_path = str(Path(args.inventory_xlsx).expanduser().resolve())

    results: list[VariantResult] = []
    for variant in variants:
        result = run_variant(
            dp,
            process,
            variant,
            methods=methods,
            inventory_xlsx=inv_path,
            ref_year=int(args.ref_year),
            act_idx=int(args.act_idx),
            amount=float(args.amount),
            max_depth=int(args.max_depth),
            show_progress=bool(args.show_progress),
            attribute_to_roots=bool(args.attribute_to_roots),
            solver_mode=str(args.solver_mode),
            iterative_rtol=float(args.iterative_rtol),
        )
        results.append(result)

    ref = results[0]
    print("\n=== Finalize Variant Benchmark ===")
    print(
        "variant | lca_s | finalize_s | peak_lca | peak_finalize | "
        "score_sum | d_score_abs | inv_nnz | score_nnz | builders_after"
    )
    for row in results:
        d_score_abs = abs(row.score_sum - ref.score_sum)
        print(
            f"{row.variant} | {row.lca_seconds:.3f} | {row.finalize_seconds:.3f} | "
            f"{_fmt_bytes(row.peak_rss_lca)} | {_fmt_bytes(row.peak_rss_finalize)} | "
            f"{row.score_sum:.12f} | {d_score_abs:.6e} | {row.inventory_nnz} | "
            f"{row.score_nnz} | {row.builders_after}"
        )


if __name__ == "__main__":
    main()
