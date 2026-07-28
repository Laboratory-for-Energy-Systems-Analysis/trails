"""Bounded-memory builders and helpers for temporal sparse inventories."""

from __future__ import annotations

from dataclasses import dataclass
import heapq
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Any, Iterable, Literal
import uuid

import dask
import dask.array as da
import numpy as np
import sparse
from tqdm import tqdm

try:  # sparse currently depends on numba, but keep a functional fallback.
    import numba
except ImportError:  # pragma: no cover
    numba = None


InventoryBackend = Literal["auto", "coo", "chunked"]

DEFAULT_INVENTORY_MEMORY_BUDGET = 256 * 2**20
DEFAULT_ACTIVITY_BLOCK_SIZE = 4_096
DEFAULT_ROOT_BLOCK_SIZE = 1_024
DEFAULT_RUN_TARGET_ENTRIES = 2_000_000
MIN_RUN_TARGET_ENTRIES = 10_000
MIN_FREE_DISK_RESERVE = 512 * 2**20


def _dtype_bytes(dtype: np.dtype | type) -> int:
    return int(np.dtype(dtype).itemsize)


def estimate_flush_peak_bytes(
    n_entries: int,
    *,
    value_dtype: np.dtype | type,
    safety_factor: float = 1.15,
) -> int:
    """Conservatively estimate peak bytes for sorting one buffered run."""
    n = max(0, int(n_entries))
    value_bytes = _dtype_bytes(value_dtype)
    # Existing key/value parts, contiguous inputs, sort order, sorted copies,
    # duplicate mask/group indices, aggregate output, and zero-filter mask.
    per_entry = (
        (8 + value_bytes)
        + (8 + value_bytes)
        + 8
        + (8 + value_bytes)
        + 1
        + 8
        + (8 + value_bytes)
        + 1
    )
    return int(math.ceil(n * per_entry * float(safety_factor)))


def estimate_merge_peak_bytes(
    left_entries: int,
    right_entries: int,
    *,
    value_dtype: np.dtype | type,
    window_entries: int = 65_536,
    safety_factor: float = 1.15,
) -> int:
    """Estimate resident memory for a memory-mapped streaming merge."""
    value_bytes = _dtype_bytes(value_dtype)
    resident_entries = min(
        int(left_entries) + int(right_entries),
        max(1, int(window_entries)) * 3,
    )
    # Two input windows and one output window; file capacity is disk-backed.
    return int(
        math.ceil(
            resident_entries * (8 + value_bytes) * float(safety_factor)
        )
    )


def estimate_decode_peak_bytes(
    n_entries: int,
    *,
    ndim: int,
    value_dtype: np.dtype | type,
    coordinate_dtype: np.dtype | type = np.int32,
    safety_factor: float = 1.15,
) -> int:
    """Estimate memory needed to decode one final sparse block."""
    n = max(0, int(n_entries))
    per_entry = (
        8
        + _dtype_bytes(value_dtype)
        + int(ndim) * _dtype_bytes(coordinate_dtype)
    )
    return int(math.ceil(n * per_entry * float(safety_factor)))


def estimate_materialization_peak_bytes(
    n_entries: int,
    *,
    ndim: int,
    value_dtype: np.dtype | type,
    coordinate_dtype: np.dtype | type = np.int64,
    safety_factor: float = 1.5,
) -> int:
    """Estimate a global COO materialization, including concatenation."""
    n = max(0, int(n_entries))
    sparse_bytes = n * (
        _dtype_bytes(value_dtype)
        + int(ndim) * _dtype_bytes(coordinate_dtype)
    )
    return int(math.ceil(sparse_bytes * float(safety_factor)))


def safe_run_target_entries(
    memory_budget: int,
    *,
    value_dtype: np.dtype | type,
) -> int:
    """Choose a run size whose predicted flush fits the memory budget."""
    budget = int(memory_budget)
    if budget <= 0:
        raise ValueError("inventory_memory_budget must be a positive integer")
    minimum_budget = estimate_flush_peak_bytes(
        MIN_RUN_TARGET_ENTRIES, value_dtype=value_dtype
    )
    if budget < minimum_budget:
        raise ValueError(
            "inventory_memory_budget is too small for the minimum bounded "
            f"run ({minimum_budget:,} bytes required)"
        )

    lo = 1
    hi = DEFAULT_RUN_TARGET_ENTRIES
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if estimate_flush_peak_bytes(mid, value_dtype=value_dtype) <= budget:
            lo = mid
        else:
            hi = mid - 1
    return max(MIN_RUN_TARGET_ENTRIES, min(DEFAULT_RUN_TARGET_ENTRIES, lo))


@dataclass(frozen=True)
class SparseRun:
    """One sorted unique key/value run stored in two memory-mapped files."""

    key_path: Path
    value_path: Path
    count: int
    capacity: int
    level: int


@dataclass(frozen=True)
class InventoryBlock:
    """Manifest entry for one finalized sparse block."""

    activity_block: int
    year_index: int
    root_block: int | None
    activity_start: int
    activity_stop: int
    root_start: int | None
    root_stop: int | None
    run: SparseRun


def _merge_sorted_unique_python(
    left_keys: np.ndarray,
    left_values: np.ndarray,
    right_keys: np.ndarray,
    right_values: np.ndarray,
    out_keys: np.ndarray,
    out_values: np.ndarray,
) -> int:
    i = 0
    j = 0
    k = 0
    n_left = int(left_keys.size)
    n_right = int(right_keys.size)
    while i < n_left and j < n_right:
        left_key = int(left_keys[i])
        right_key = int(right_keys[j])
        if left_key < right_key:
            out_keys[k] = left_key
            out_values[k] = left_values[i]
            i += 1
        elif right_key < left_key:
            out_keys[k] = right_key
            out_values[k] = right_values[j]
            j += 1
        else:
            value = left_values[i] + right_values[j]
            if value != 0:
                out_keys[k] = left_key
                out_values[k] = value
                k += 1
            i += 1
            j += 1
            continue
        k += 1
    while i < n_left:
        out_keys[k] = left_keys[i]
        out_values[k] = left_values[i]
        i += 1
        k += 1
    while j < n_right:
        out_keys[k] = right_keys[j]
        out_values[k] = right_values[j]
        j += 1
        k += 1
    return k


if numba is not None:
    _merge_sorted_unique = numba.njit(cache=False)(_merge_sorted_unique_python)
else:  # pragma: no cover
    _merge_sorted_unique = _merge_sorted_unique_python


def _chunk_lengths(size: int, block_size: int) -> tuple[int, ...]:
    size_i = int(size)
    block_i = max(1, int(block_size))
    return tuple(
        min(block_i, size_i - start) for start in range(0, size_i, block_i)
    )


def _axis_chunks_with_sparse_singletons(
    size: int,
    nonempty_indices: Iterable[int],
) -> tuple[int, ...]:
    """Make singleton chunks at used indices and coalesce unused gaps."""
    size_i = int(size)
    used = sorted({int(i) for i in nonempty_indices if 0 <= int(i) < size_i})
    if not used:
        return (size_i,)
    chunks: list[int] = []
    cursor = 0
    for idx in used:
        if idx > cursor:
            chunks.append(idx - cursor)
        chunks.append(1)
        cursor = idx + 1
    if cursor < size_i:
        chunks.append(size_i - cursor)
    return tuple(chunks)


def _chunk_starts(chunks: tuple[int, ...]) -> tuple[int, ...]:
    starts: list[int] = []
    cursor = 0
    for length in chunks:
        starts.append(cursor)
        cursor += int(length)
    return tuple(starts)


def _load_sparse_block(
    key_path: str,
    value_path: str,
    count: int,
    shape: tuple[int, ...],
    *,
    has_root: bool,
    n_flows: int,
    root_width: int,
    value_dtype: str,
) -> sparse.COO:
    keys_map = np.load(key_path, mmap_mode="r")
    values_map = np.load(value_path, mmap_mode="r")
    keys = np.asarray(keys_map[: int(count)], dtype=np.int64)
    values = np.asarray(values_map[: int(count)], dtype=np.dtype(value_dtype))
    if not keys.size:
        return sparse.zeros(shape, dtype=np.dtype(value_dtype))

    q = keys.copy()
    if has_root:
        coords = np.empty((4, keys.size), dtype=np.int32)
        np.remainder(q, int(root_width), out=coords[3])
        np.floor_divide(q, int(root_width), out=q)
        np.remainder(q, int(n_flows), out=coords[1])
        np.floor_divide(q, int(n_flows), out=q)
        coords[0] = q.astype(np.int32, copy=False)
        coords[2].fill(0)
    else:
        coords = np.empty((3, keys.size), dtype=np.int32)
        np.remainder(q, int(n_flows), out=coords[1])
        np.floor_divide(q, int(n_flows), out=q)
        coords[0] = q.astype(np.int32, copy=False)
        coords[2].fill(0)

    return sparse.COO(
        coords,
        values,
        shape=shape,
        has_duplicates=False,
        sorted=True,
        idx_dtype=np.int32,
    )


def _empty_sparse_block(shape: tuple[int, ...], dtype: str) -> sparse.COO:
    return sparse.zeros(shape, dtype=np.dtype(dtype))


class ChunkedInventoryBuilder:
    """Accumulate and canonicalize a large sparse inventory within a budget."""

    def __init__(
        self,
        *,
        n_activities: int,
        n_flows: int,
        n_years: int,
        has_root: bool,
        value_dtype: np.dtype | type,
        memory_budget: int = DEFAULT_INVENTORY_MEMORY_BUDGET,
        store: str | Path | None = None,
        activity_block_size: int = DEFAULT_ACTIVITY_BLOCK_SIZE,
        root_block_size: int = DEFAULT_ROOT_BLOCK_SIZE,
    ) -> None:
        self.n_activities = int(n_activities)
        self.n_flows = int(n_flows)
        self.n_years = int(n_years)
        self.has_root = bool(has_root)
        self.value_dtype = np.dtype(value_dtype)
        self.memory_budget = int(memory_budget)
        if self.memory_budget <= 0:
            raise ValueError("inventory_memory_budget must be a positive integer")
        self.activity_block_size = max(1, int(activity_block_size))
        self.root_block_size = max(1, int(root_block_size))
        self.run_target = safe_run_target_entries(
            self.memory_budget, value_dtype=self.value_dtype
        )

        self.owned_store = store is None
        if store is None:
            self.store_path = Path(tempfile.mkdtemp(prefix="trails-inventory-"))
        else:
            root = Path(store).expanduser().resolve()
            root.mkdir(parents=True, exist_ok=True)
            self.store_path = root / f"run-{uuid.uuid4().hex}"
            self.store_path.mkdir(parents=True, exist_ok=False)

        self._buffers: dict[tuple[int, ...], list[tuple[np.ndarray, np.ndarray]]] = {}
        self._buffer_entries: dict[tuple[int, ...], int] = {}
        self._buffered_total_entries = 0
        self._levels: dict[tuple[int, ...], dict[int, SparseRun]] = {}
        self._final_runs: dict[tuple[int, ...], SparseRun] = {}
        self._run_counter = 0
        self._closed = False
        self._finalized = False
        self.raw_entries = 0
        self.canonical_entries = 0
        self.bytes_written = 0
        self.bytes_merged = 0
        self.peak_buffer_bytes = 0
        self.append_seconds = 0.0
        self.online_flush_seconds = 0.0
        self.finalize_flush_seconds = 0.0
        self.online_merge_seconds = 0.0
        self.finalize_merge_seconds = 0.0
        self._finishing = False

    def _ensure_disk_capacity(self, additional_bytes: int, operation: str) -> None:
        """Fail before creating a run when the backing store is too full."""
        required = max(0, int(additional_bytes))
        free = int(shutil.disk_usage(self.store_path).free)
        reserve = max(MIN_FREE_DISK_RESERVE, 2 * self.memory_budget)
        if free < required + reserve:
            raise OSError(
                f"Not enough free space to {operation}: "
                f"{required / 2**30:.2f} GiB additional space is required, "
                f"with a {reserve / 2**30:.2f} GiB safety reserve, but only "
                f"{free / 2**30:.2f} GiB is free in {self.store_path}."
            )

    @property
    def buffered_entries(self) -> int:
        return int(self._buffered_total_entries)

    @property
    def buffered_bytes(self) -> int:
        return self._buffered_total_entries * (8 + self.value_dtype.itemsize)

    def _spill_largest_buffers(self, *, target_fraction: float = 0.75) -> None:
        """Flush a batch of large buffers, amortizing the partition scan."""
        if not self._buffer_entries:
            raise RuntimeError("No buffered inventory partition is available")
        bytes_per_entry = 8 + self.value_dtype.itemsize
        target_entries = int(
            self.memory_budget * float(target_fraction) / bytes_per_entry
        )
        heap = [
            (-int(count), partition)
            for partition, count in self._buffer_entries.items()
        ]
        heapq.heapify(heap)
        while self._buffered_total_entries > target_entries and heap:
            _, partition = heapq.heappop(heap)
            self._flush_partition(partition)
        if self._buffered_total_entries > target_entries:
            raise RuntimeError("Unable to spill inventory buffers below target")

    @property
    def nnz(self) -> int:
        if not self._finalized:
            return int(self.raw_entries)
        return int(sum(run.count for run in self._final_runs.values()))

    def _partition_and_local_keys(
        self,
        activities: np.ndarray,
        flows: np.ndarray,
        year_indices: np.ndarray,
        roots: np.ndarray | None,
    ) -> tuple[np.ndarray, np.ndarray]:
        activity_blocks = activities // self.activity_block_size
        activity_local = activities % self.activity_block_size
        if self.has_root:
            if roots is None:
                raise ValueError("root coordinates required for root-attributed inventory")
            root_blocks = roots // self.root_block_size
            root_local = roots % self.root_block_size
            partition_ids = (
                (activity_blocks * self.n_years + year_indices)
                * math.ceil(self.n_activities / self.root_block_size)
                + root_blocks
            )
            local_keys = (
                (activity_local * self.n_flows + flows) * self.root_block_size
                + root_local
            )
        else:
            partition_ids = activity_blocks * self.n_years + year_indices
            local_keys = activity_local * self.n_flows + flows
        return (
            partition_ids.astype(np.int64, copy=False),
            local_keys.astype(np.int64, copy=False),
        )

    def _partition_tuple(self, partition_id: int) -> tuple[int, ...]:
        pid = int(partition_id)
        if self.has_root:
            n_root_blocks = math.ceil(self.n_activities / self.root_block_size)
            root_block = pid % n_root_blocks
            q = pid // n_root_blocks
            year_index = q % self.n_years
            activity_block = q // self.n_years
            return (activity_block, year_index, root_block)
        year_index = pid % self.n_years
        activity_block = pid // self.n_years
        return (activity_block, year_index)

    def append(
        self,
        activities: np.ndarray,
        flows: np.ndarray,
        year_indices: np.ndarray,
        values: np.ndarray,
        *,
        roots: np.ndarray | None = None,
    ) -> None:
        append_started = time.perf_counter()
        if self._finalized or self._closed:
            raise RuntimeError("Cannot append to a finalized inventory builder")
        acts = np.asarray(activities, dtype=np.int64)
        flows_arr = np.asarray(flows, dtype=np.int64)
        years = np.asarray(year_indices, dtype=np.int64)
        vals = np.asarray(values, dtype=self.value_dtype)
        roots_arr = None if roots is None else np.asarray(roots, dtype=np.int64)
        if not (acts.size == flows_arr.size == years.size == vals.size):
            raise ValueError("Inventory coordinate and value arrays must align")
        if roots_arr is not None and roots_arr.size != acts.size:
            raise ValueError("Root coordinate array must align with inventory entries")
        if not acts.size:
            self.append_seconds += time.perf_counter() - append_started
            return

        keep = vals != 0
        if not np.any(keep):
            self.append_seconds += time.perf_counter() - append_started
            return
        acts = acts[keep]
        flows_arr = flows_arr[keep]
        years = years[keep]
        vals = vals[keep]
        if roots_arr is not None:
            roots_arr = roots_arr[keep]

        partition_ids, local_keys = self._partition_and_local_keys(
            acts, flows_arr, years, roots_arr
        )
        order = np.argsort(partition_ids, kind="stable")
        partition_ids = partition_ids[order]
        local_keys = local_keys[order]
        vals = vals[order]
        starts = np.r_[0, np.flatnonzero(partition_ids[1:] != partition_ids[:-1]) + 1]
        stops = np.r_[starts[1:], partition_ids.size]
        self.raw_entries += int(partition_ids.size)

        for start, stop in zip(starts, stops):
            partition = self._partition_tuple(int(partition_ids[int(start)]))
            cursor = int(start)
            stop_i = int(stop)
            while cursor < stop_i:
                if self.buffered_bytes >= self.memory_budget:
                    self._spill_largest_buffers()
                current = int(self._buffer_entries.get(partition, 0))
                room = max(1, self.run_target - current)
                bytes_per_entry = 8 + self.value_dtype.itemsize
                budget_room = max(
                    1,
                    (self.memory_budget - self.buffered_bytes)
                    // bytes_per_entry,
                )
                take = min(room, int(budget_room), stop_i - cursor)
                key_part = local_keys[cursor : cursor + take].copy()
                value_part = vals[cursor : cursor + take].copy()
                self._buffers.setdefault(partition, []).append((key_part, value_part))
                new_size = current + take
                self._buffer_entries[partition] = new_size
                self._buffered_total_entries += int(take)
                cursor += take
                self.peak_buffer_bytes = max(
                    self.peak_buffer_bytes, self.buffered_bytes
                )
                if self._buffer_entries[partition] >= self.run_target:
                    self._flush_partition(partition)
        self.append_seconds += time.perf_counter() - append_started

    def append_linear_global(self, keys: np.ndarray, values: np.ndarray) -> None:
        """Decode legacy global keys and append them to bounded partitions."""
        keys_arr = np.asarray(keys, dtype=np.int64)
        vals = np.asarray(values, dtype=self.value_dtype)
        if not keys_arr.size:
            return
        q = keys_arr.copy()
        if self.has_root:
            roots = np.empty(q.size, dtype=np.int64)
            np.remainder(q, self.n_activities, out=roots)
            np.floor_divide(q, self.n_activities, out=q)
        else:
            roots = None
        years = np.empty(q.size, dtype=np.int64)
        flows = np.empty(q.size, dtype=np.int64)
        np.remainder(q, self.n_years, out=years)
        np.floor_divide(q, self.n_years, out=q)
        np.remainder(q, self.n_flows, out=flows)
        np.floor_divide(q, self.n_flows, out=q)
        self.append(q, flows, years, vals, roots=roots)

    def _new_run_paths(self, partition: tuple[int, ...], level: int) -> tuple[Path, Path]:
        self._run_counter += 1
        label = "-".join(str(x) for x in partition)
        stem = f"run-{label}-l{level}-{self._run_counter}"
        return self.store_path / f"{stem}-keys.npy", self.store_path / f"{stem}-values.npy"

    def _write_run(
        self,
        partition: tuple[int, ...],
        keys: np.ndarray,
        values: np.ndarray,
        *,
        level: int,
    ) -> SparseRun:
        key_path, value_path = self._new_run_paths(partition, level)
        count = int(keys.size)
        self._ensure_disk_capacity(
            count * (8 + self.value_dtype.itemsize) + 512,
            "write an inventory run",
        )
        np.save(key_path, np.asarray(keys, dtype=np.int64), allow_pickle=False)
        np.save(
            value_path,
            np.asarray(values, dtype=self.value_dtype),
            allow_pickle=False,
        )
        self.bytes_written += count * (8 + self.value_dtype.itemsize)
        return SparseRun(key_path, value_path, count, count, int(level))

    def _flush_partition(self, partition: tuple[int, ...]) -> None:
        flush_started = time.perf_counter()
        parts = self._buffers.pop(partition, [])
        removed_entries = int(self._buffer_entries.pop(partition, 0))
        self._buffered_total_entries -= removed_entries
        if not parts:
            return
        keys = np.concatenate([part[0] for part in parts]).astype(np.int64, copy=False)
        values = np.concatenate([part[1] for part in parts]).astype(
            self.value_dtype, copy=False
        )
        order = np.argsort(keys, kind="quicksort")
        keys = keys[order]
        values = values[order]
        first = np.empty(keys.size, dtype=bool)
        first[0] = True
        first[1:] = keys[1:] != keys[:-1]
        starts = np.flatnonzero(first)
        values_agg = np.add.reduceat(values, starts).astype(
            self.value_dtype, copy=False
        )
        keys_agg = keys[starts]
        keep = values_agg != 0
        run = self._write_run(
            partition,
            keys_agg[keep],
            values_agg[keep],
            level=0,
        )
        self._add_run(partition, run)
        seconds = time.perf_counter() - flush_started
        if self._finishing:
            self.finalize_flush_seconds += seconds
        else:
            self.online_flush_seconds += seconds

    def _add_run(self, partition: tuple[int, ...], run: SparseRun) -> None:
        levels = self._levels.setdefault(partition, {})
        current = run
        while current.level in levels:
            previous = levels.pop(current.level)
            current = self._merge_runs(partition, previous, current)
        levels[current.level] = current

    def _merge_runs(
        self,
        partition: tuple[int, ...],
        left: SparseRun,
        right: SparseRun,
    ) -> SparseRun:
        merge_started = time.perf_counter()
        level = max(left.level, right.level) + 1
        key_path, value_path = self._new_run_paths(partition, level)
        capacity = int(left.count + right.count)
        merge_peak = estimate_merge_peak_bytes(
            left.count,
            right.count,
            value_dtype=self.value_dtype,
        )
        if merge_peak > self.memory_budget:
            raise MemoryError(
                "The estimated streaming-merge peak exceeds the inventory "
                f"memory budget ({merge_peak:,} > {self.memory_budget:,} bytes)."
            )
        self._ensure_disk_capacity(
            capacity * (8 + self.value_dtype.itemsize) + 512,
            "merge inventory runs",
        )
        out_keys = np.lib.format.open_memmap(
            key_path, mode="w+", dtype=np.int64, shape=(capacity,)
        )
        out_values = np.lib.format.open_memmap(
            value_path,
            mode="w+",
            dtype=self.value_dtype,
            shape=(capacity,),
        )
        left_keys_map = np.load(left.key_path, mmap_mode="r")
        left_values_map = np.load(left.value_path, mmap_mode="r")
        right_keys_map = np.load(right.key_path, mmap_mode="r")
        right_values_map = np.load(right.value_path, mmap_mode="r")
        count = int(
            _merge_sorted_unique(
                left_keys_map[: left.count],
                left_values_map[: left.count],
                right_keys_map[: right.count],
                right_values_map[: right.count],
                out_keys,
                out_values,
            )
        )
        out_keys.flush()
        out_values.flush()
        del out_keys, out_values
        del left_keys_map, left_values_map, right_keys_map, right_values_map
        self.bytes_merged += capacity * (8 + self.value_dtype.itemsize)
        for path in (
            left.key_path,
            left.value_path,
            right.key_path,
            right.value_path,
        ):
            path.unlink(missing_ok=True)
        result = SparseRun(key_path, value_path, count, capacity, level)
        seconds = time.perf_counter() - merge_started
        if self._finishing:
            self.finalize_merge_seconds += seconds
        else:
            self.online_merge_seconds += seconds
        return result

    def _finish_runs(self, *, show_progress: bool = False) -> None:
        self._finishing = True
        pending_buffers = list(self._buffers)
        buffer_iter = tqdm(
            pending_buffers,
            desc="TRAILS LCA [3/5] Flush inventory buffers",
            unit="partition",
            leave=True,
            disable=not show_progress,
        )
        n_buffers = len(pending_buffers)
        for index, partition in enumerate(buffer_iter, start=1):
            self._flush_partition(partition)
            if index == 1 or index % 128 == 0 or index == n_buffers:
                buffer_iter.set_postfix(
                    raw=f"{self.raw_entries:,}",
                    buffered=f"{self.buffered_entries:,}",
                )
        level_items = list(self._levels.items())
        level_iter = tqdm(
            level_items,
            desc="TRAILS LCA [3/5] Merge inventory runs",
            unit="partition",
            leave=True,
            disable=not show_progress,
        )
        n_level_items = len(level_items)
        for index, (partition, levels) in enumerate(level_iter, start=1):
            runs = [levels[level] for level in sorted(levels)]
            if not runs:
                continue
            current = runs[0]
            for run in runs[1:]:
                current = self._merge_runs(partition, current, run)
            self._final_runs[partition] = current
            if index == 1 or index % 128 == 0 or index == n_level_items:
                level_iter.set_postfix(
                    written=f"{self.bytes_written / 2**20:.0f} MiB",
                    merged=f"{self.bytes_merged / 2**20:.0f} MiB",
                )
        self._levels.clear()
        self.canonical_entries = int(
            sum(run.count for run in self._final_runs.values())
        )
        self._finishing = False

    def _block_manifest(self) -> list[InventoryBlock]:
        blocks: list[InventoryBlock] = []
        for partition, run in sorted(self._final_runs.items()):
            activity_block = int(partition[0])
            year_index = int(partition[1])
            root_block = int(partition[2]) if self.has_root else None
            activity_start = activity_block * self.activity_block_size
            activity_stop = min(
                activity_start + self.activity_block_size, self.n_activities
            )
            root_start = (
                root_block * self.root_block_size if root_block is not None else None
            )
            root_stop = (
                min(root_start + self.root_block_size, self.n_activities)
                if root_start is not None
                else None
            )
            blocks.append(
                InventoryBlock(
                    activity_block=activity_block,
                    year_index=year_index,
                    root_block=root_block,
                    activity_start=activity_start,
                    activity_stop=activity_stop,
                    root_start=root_start,
                    root_stop=root_stop,
                    run=run,
                )
            )
        return blocks

    def _write_manifest(self, blocks: list[InventoryBlock]) -> None:
        payload = {
            "version": 1,
            "shape": [
                self.n_activities,
                self.n_flows,
                self.n_years,
                *([self.n_activities] if self.has_root else []),
            ],
            "value_dtype": self.value_dtype.str,
            "has_root": self.has_root,
            "raw_entries": self.raw_entries,
            "canonical_entries": self.canonical_entries,
            "peak_buffer_bytes": self.peak_buffer_bytes,
            "bytes_written": self.bytes_written,
            "bytes_merged": self.bytes_merged,
            "blocks": [
                {
                    "activity_block": block.activity_block,
                    "year_index": block.year_index,
                    "root_block": block.root_block,
                    "activity_start": block.activity_start,
                    "activity_stop": block.activity_stop,
                    "root_start": block.root_start,
                    "root_stop": block.root_stop,
                    "key_path": block.run.key_path.name,
                    "value_path": block.run.value_path.name,
                    "count": block.run.count,
                    "capacity": block.run.capacity,
                }
                for block in blocks
            ],
        }
        tmp_path = self.store_path / "manifest.json.tmp"
        final_path = self.store_path / "manifest.json"
        tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        os.replace(tmp_path, final_path)

    def _delayed_block(
        self,
        block: InventoryBlock | None,
        *,
        shape: tuple[int, ...],
        root_width: int,
    ) -> da.Array:
        if block is None:
            delayed_value = dask.delayed(_empty_sparse_block)(
                shape, self.value_dtype.str
            )
        else:
            decode_peak = estimate_decode_peak_bytes(
                block.run.count,
                ndim=len(shape),
                value_dtype=self.value_dtype,
            )
            if decode_peak > self.memory_budget:
                raise MemoryError(
                    "A finalized inventory block is too large to decode within "
                    f"the configured budget ({decode_peak:,} > "
                    f"{self.memory_budget:,} bytes). Reduce activity_block_size "
                    "or root_block_size when constructing the builder."
                )
            delayed_value = dask.delayed(_load_sparse_block)(
                str(block.run.key_path),
                str(block.run.value_path),
                int(block.run.count),
                shape,
                has_root=self.has_root,
                n_flows=self.n_flows,
                root_width=int(root_width),
                value_dtype=self.value_dtype.str,
            )
        meta_shape = tuple(0 for _ in shape)
        meta = sparse.zeros(meta_shape, dtype=self.value_dtype)
        return da.from_delayed(
            delayed_value,
            shape=shape,
            dtype=self.value_dtype,
            meta=meta,
        )

    def _assemble_dask_array(
        self, blocks: list[InventoryBlock], *, show_progress: bool = False
    ) -> da.Array:
        block_map = {
            (
                block.activity_block,
                block.year_index,
                block.root_block if self.has_root else None,
            ): block
            for block in blocks
        }
        activity_chunks = _chunk_lengths(
            self.n_activities, self.activity_block_size
        )
        activity_starts = _chunk_starts(activity_chunks)
        used_years = {block.year_index for block in blocks}
        year_chunks = _axis_chunks_with_sparse_singletons(self.n_years, used_years)
        year_starts = _chunk_starts(year_chunks)

        if self.has_root:
            root_chunks = _chunk_lengths(self.n_activities, self.root_block_size)
            root_starts = _chunk_starts(root_chunks)
            activity_arrays: list[da.Array] = []
            activity_iter = tqdm(
                list(zip(activity_starts, activity_chunks)),
                desc="TRAILS LCA [4/5] Assemble inventory blocks",
                unit="activity-block",
                leave=True,
                disable=not show_progress,
            )
            for activity_start, activity_width in activity_iter:
                activity_block = activity_start // self.activity_block_size
                year_arrays: list[da.Array] = []
                for year_start, year_width in zip(year_starts, year_chunks):
                    root_arrays: list[da.Array] = []
                    for root_start, root_width in zip(root_starts, root_chunks):
                        root_block = root_start // self.root_block_size
                        block = None
                        if year_width == 1:
                            block = block_map.get(
                                (activity_block, year_start, root_block)
                            )
                        shape = (
                            activity_width,
                            self.n_flows,
                            year_width,
                            root_width,
                        )
                        root_arrays.append(
                            self._delayed_block(
                                block,
                                shape=shape,
                                root_width=self.root_block_size,
                            )
                        )
                    year_arrays.append(da.concatenate(root_arrays, axis=3))
                activity_arrays.append(da.concatenate(year_arrays, axis=2))
            return da.concatenate(activity_arrays, axis=0)

        activity_arrays = []
        activity_iter = tqdm(
            list(zip(activity_starts, activity_chunks)),
            desc="TRAILS LCA [4/5] Assemble inventory blocks",
            unit="activity-block",
            leave=True,
            disable=not show_progress,
        )
        for activity_start, activity_width in activity_iter:
            activity_block = activity_start // self.activity_block_size
            year_arrays = []
            for year_start, year_width in zip(year_starts, year_chunks):
                block = (
                    block_map.get((activity_block, year_start, None))
                    if year_width == 1
                    else None
                )
                shape = (activity_width, self.n_flows, year_width)
                year_arrays.append(
                    self._delayed_block(block, shape=shape, root_width=1)
                )
            activity_arrays.append(da.concatenate(year_arrays, axis=2))
        return da.concatenate(activity_arrays, axis=0)

    def finalize(self, *, show_progress: bool = False) -> da.Array:
        if self._closed:
            raise RuntimeError("Inventory builder is closed")
        if self._finalized:
            raise RuntimeError("Inventory builder was already finalized")
        self._finish_runs(show_progress=show_progress)
        blocks = self._block_manifest()
        self._write_manifest(blocks)
        result = self._assemble_dask_array(blocks, show_progress=show_progress)
        self._finalized = True
        return result

    def diagnostics(self) -> dict[str, int | float | str | bool]:
        return {
            "backend": "chunked",
            "store": str(self.store_path),
            "owned_store": self.owned_store,
            "raw_entries": int(self.raw_entries),
            "canonical_entries": int(self.canonical_entries),
            "peak_buffer_bytes": int(self.peak_buffer_bytes),
            "bytes_written": int(self.bytes_written),
            "bytes_merged": int(self.bytes_merged),
            "run_target_entries": int(self.run_target),
            "memory_budget": int(self.memory_budget),
            "free_disk_bytes": int(shutil.disk_usage(self.store_path).free),
            "append_seconds": float(self.append_seconds),
            "online_flush_seconds": float(self.online_flush_seconds),
            "finalize_flush_seconds": float(self.finalize_flush_seconds),
            "online_merge_seconds": float(self.online_merge_seconds),
            "finalize_merge_seconds": float(self.finalize_merge_seconds),
        }

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._buffers.clear()
        self._buffer_entries.clear()
        self._buffered_total_entries = 0
        if self.owned_store:
            shutil.rmtree(self.store_path, ignore_errors=True)

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass


def is_chunked_sparse(data: Any) -> bool:
    return isinstance(data, da.Array)


def iter_sparse_blocks(
    array: da.Array,
    *,
    primary_axis: int | None = None,
) -> Iterable[tuple[tuple[slice, ...], sparse.COO]]:
    """Yield non-empty Dask chunks sequentially with global axis slices."""
    delayed_blocks = array.to_delayed(optimize_graph=True)
    axis_starts = [
        _chunk_starts(tuple(int(length) for length in axis_chunks))
        for axis_chunks in array.chunks
    ]
    block_indices = list(np.ndindex(delayed_blocks.shape))
    if primary_axis is not None:
        axis = int(primary_axis)
        block_indices.sort(
            key=lambda index: (index[axis],) + tuple(
                value for i, value in enumerate(index) if i != axis
            )
        )
    for block_index in block_indices:
        block = delayed_blocks[block_index].compute(scheduler="synchronous")
        if not isinstance(block, sparse.COO):
            block = sparse.COO.from_numpy(np.asarray(block))
        if block.nnz:
            slices = tuple(
                slice(
                    int(axis_starts[axis][chunk_index]),
                    int(axis_starts[axis][chunk_index] + block.shape[axis]),
                )
                for axis, chunk_index in enumerate(block_index)
            )
            yield slices, block
