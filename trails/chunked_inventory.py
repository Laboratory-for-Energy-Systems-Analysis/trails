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


InventoryBackend = Literal["auto", "coo", "chunked", "factorized"]

DEFAULT_INVENTORY_MEMORY_BUDGET = 256 * 2**20
DEFAULT_ACTIVITY_BLOCK_SIZE = 4_096
DEFAULT_ROOT_BLOCK_SIZE = 4_096
DEFAULT_YEAR_BLOCK_SIZE = 8
DEFAULT_RUN_TARGET_ENTRIES = 2_000_000
MIN_RUN_TARGET_ENTRIES = 10_000
MIN_FREE_DISK_RESERVE = 512 * 2**20
DEFAULT_SHARD_BUCKET_COUNT = 16
DEFAULT_BUCKET_COMPACTION_BYTES = 64 * 2**20


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
    return int(math.ceil(resident_entries * (8 + value_bytes) * float(safety_factor)))


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
        8 + _dtype_bytes(value_dtype) + int(ndim) * _dtype_bytes(coordinate_dtype)
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
        _dtype_bytes(value_dtype) + int(ndim) * _dtype_bytes(coordinate_dtype)
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
    offset: int = 0
    exclusive: bool = False


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
    return tuple(min(block_i, size_i - start) for start in range(0, size_i, block_i))


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


def _load_sparse_block_group(
    runs: tuple[tuple[str, str, int, int, int], ...],
    shape: tuple[int, ...],
    *,
    has_root: bool,
    n_flows: int,
    root_width: int,
    year_stride: int,
    value_dtype: str,
) -> sparse.COO:
    total = sum(int(item[3]) for item in runs)
    if total == 0:
        return sparse.zeros(shape, dtype=np.dtype(value_dtype))

    ndim = 4 if has_root else 3
    coords = np.empty((ndim, total), dtype=np.int64)
    values = np.empty(total, dtype=np.dtype(value_dtype))
    cursor = 0
    for key_path, value_path, offset, count, year_offset in runs:
        count_i = int(count)
        stop = cursor + count_i
        keys_map = np.memmap(
            key_path,
            mode="r",
            dtype=np.int64,
            offset=int(offset) * np.dtype(np.int64).itemsize,
            shape=(count_i,),
        )
        values_map = np.memmap(
            value_path,
            mode="r",
            dtype=np.dtype(value_dtype),
            offset=int(offset) * np.dtype(value_dtype).itemsize,
            shape=(count_i,),
        )
        try:
            q = np.asarray(keys_map, dtype=np.int64).copy()
            values[cursor:stop] = values_map
        finally:
            keys_map._mmap.close()
            values_map._mmap.close()
        if has_root:
            np.remainder(q, int(root_width), out=coords[3, cursor:stop])
            np.floor_divide(q, int(root_width), out=q)
            np.remainder(q, int(year_stride), out=coords[2, cursor:stop])
            coords[2, cursor:stop] += int(year_offset)
            np.floor_divide(q, int(year_stride), out=q)
            np.remainder(q, int(n_flows), out=coords[1, cursor:stop])
            np.floor_divide(q, int(n_flows), out=q)
        else:
            np.remainder(q, int(year_stride), out=coords[2, cursor:stop])
            coords[2, cursor:stop] += int(year_offset)
            np.floor_divide(q, int(year_stride), out=q)
            np.remainder(q, int(n_flows), out=coords[1, cursor:stop])
            np.floor_divide(q, int(n_flows), out=q)
        coords[0, cursor:stop] = q
        cursor = stop

    return sparse.COO(
        coords,
        values,
        shape=shape,
        has_duplicates=False,
        sorted=len(runs) <= 1,
        idx_dtype=np.int64,
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
        year_block_size: int = DEFAULT_YEAR_BLOCK_SIZE,
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
        self.year_block_size = max(1, int(year_block_size))
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
        self._pending_runs: dict[tuple[int, ...], list[SparseRun]] = {}
        self._final_runs: dict[tuple[int, ...], SparseRun] = {}
        self._input_shard_paths: set[Path] = set()
        self._bucket_paths: dict[int, tuple[Path, Path]] = {}
        self._bucket_entries: dict[int, int] = {}
        self._bucket_appended_bytes: dict[int, int] = {}
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
        self.online_compaction_count = 0
        self.finalize_compaction_count = 0
        self._finishing = False
        self.dask_block_count = 0
        self._activity_reduction_cache: dict[tuple[int, ...], sparse.COO] = {}

    @staticmethod
    def _partition_bucket(partition: tuple[int, ...]) -> int:
        """Assign related runs to one of a small number of stable disk buckets."""
        value = 0
        # Keep all year blocks for one activity/root pair together. This gives
        # stable range-like buckets without allowing the long year axis to
        # multiply the number of files.
        bucket_items = (
            (partition[0], partition[-1]) if len(partition) == 3 else (partition[0],)
        )
        for item in bucket_items:
            value = (value * 1_000_003 + int(item)) & 0x7FFF_FFFF_FFFF_FFFF
        return int(value % DEFAULT_SHARD_BUCKET_COUNT)

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
        selected: list[tuple[int, ...]] = []
        selected_entries = 0
        while self._buffered_total_entries > target_entries and heap:
            negative_count, partition = heapq.heappop(heap)
            selected.append(partition)
            count = -int(negative_count)
            selected_entries += count
            self._buffered_total_entries -= count
        if self._buffered_total_entries > target_entries or not selected:
            raise RuntimeError("Unable to spill inventory buffers below target")
        # Restore the counter until the selected buffers are removed atomically.
        self._buffered_total_entries += selected_entries
        self._flush_partitions_to_shard(selected)

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
        n_year_blocks = math.ceil(self.n_years / self.year_block_size)
        year_blocks = year_indices // self.year_block_size
        year_local = year_indices % self.year_block_size
        if self.has_root:
            if roots is None:
                raise ValueError(
                    "root coordinates required for root-attributed inventory"
                )
            root_blocks = roots // self.root_block_size
            root_local = roots % self.root_block_size
            partition_ids = (activity_blocks * n_year_blocks + year_blocks) * math.ceil(
                self.n_activities / self.root_block_size
            ) + root_blocks
            local_keys = (
                (activity_local * self.n_flows + flows) * self.year_block_size
                + year_local
            ) * self.root_block_size + root_local
        else:
            partition_ids = activity_blocks * n_year_blocks + year_blocks
            local_keys = (
                activity_local * self.n_flows + flows
            ) * self.year_block_size + year_local
        return (
            partition_ids.astype(np.int64, copy=False),
            local_keys.astype(np.int64, copy=False),
        )

    def _partition_tuple(self, partition_id: int) -> tuple[int, ...]:
        pid = int(partition_id)
        n_year_blocks = math.ceil(self.n_years / self.year_block_size)
        if self.has_root:
            n_root_blocks = math.ceil(self.n_activities / self.root_block_size)
            root_block = pid % n_root_blocks
            q = pid // n_root_blocks
            year_block = q % n_year_blocks
            activity_block = q // n_year_blocks
            year_index = year_block * self.year_block_size
            return (activity_block, year_index, root_block)
        year_block = pid % n_year_blocks
        activity_block = pid // n_year_blocks
        year_index = year_block * self.year_block_size
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
                    (self.memory_budget - self.buffered_bytes) // bytes_per_entry,
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
                    self._flush_partitions_to_shard([partition])
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

    def _new_run_paths(self, label: str) -> tuple[Path, Path]:
        self._run_counter += 1
        stem = f"{label}-{self._run_counter}"
        return (
            self.store_path / f"{stem}-keys.bin",
            self.store_path / f"{stem}-values.bin",
        )

    def _canonicalize_parts(
        self, parts: list[tuple[np.ndarray, np.ndarray]]
    ) -> tuple[np.ndarray, np.ndarray]:
        keys = np.concatenate([part[0] for part in parts]).astype(np.int64, copy=False)
        values = np.concatenate([part[1] for part in parts]).astype(
            self.value_dtype, copy=False
        )
        return self._canonicalize_arrays(keys, values)

    def _flush_partitions_to_shard(self, partitions: Iterable[tuple[int, ...]]) -> None:
        """Append canonical runs to a small fixed set of bucket shard files."""
        flush_started = time.perf_counter()
        selected: list[
            tuple[tuple[int, ...], list[tuple[np.ndarray, np.ndarray]], int]
        ] = []
        raw_count = 0
        for partition in partitions:
            parts = self._buffers.pop(partition, [])
            count = int(self._buffer_entries.pop(partition, 0))
            if not parts:
                continue
            selected.append((partition, parts, count))
            raw_count += count
        self._buffered_total_entries -= raw_count
        if not selected:
            return

        self._ensure_disk_capacity(
            raw_count * (8 + self.value_dtype.itemsize),
            "write an inventory shard",
        )
        by_bucket: dict[
            int,
            list[tuple[tuple[int, ...], list[tuple[np.ndarray, np.ndarray]], int]],
        ] = {}
        for item in selected:
            by_bucket.setdefault(self._partition_bucket(item[0]), []).append(item)

        compact_buckets: list[int] = []
        written_entries = 0
        for bucket, bucket_items in by_bucket.items():
            paths = self._bucket_paths.get(bucket)
            if paths is None:
                paths = self._new_run_paths(f"bucket-{bucket}")
                self._bucket_paths[bucket] = paths
                self._bucket_entries[bucket] = 0
                self._bucket_appended_bytes[bucket] = 0
                self._input_shard_paths.update(paths)
            key_path, value_path = paths
            offset = int(self._bucket_entries[bucket])
            bucket_written = 0
            with (
                key_path.open("ab") as key_stream,
                value_path.open("ab") as value_stream,
            ):
                for partition, parts, _ in bucket_items:
                    keys, values = self._canonicalize_parts(parts)
                    count = int(keys.size)
                    if not count:
                        continue
                    keys.tofile(key_stream)
                    values.tofile(value_stream)
                    run = SparseRun(
                        key_path,
                        value_path,
                        count,
                        count,
                        0,
                        offset=offset,
                        exclusive=False,
                    )
                    self._pending_runs.setdefault(partition, []).append(run)
                    offset += count
                    bucket_written += count
            self._bucket_entries[bucket] = offset
            written_entries += bucket_written
            appended_bytes = bucket_written * (8 + self.value_dtype.itemsize)
            self._bucket_appended_bytes[bucket] += appended_bytes
            bytes_per_entry = 8 + self.value_dtype.itemsize
            existing_compacted_bytes = max(
                0,
                int(self._bucket_entries[bucket]) * bytes_per_entry
                - int(self._bucket_appended_bytes[bucket]),
            )
            compaction_threshold = max(
                int(DEFAULT_BUCKET_COMPACTION_BYTES),
                existing_compacted_bytes,
            )
            if (
                not self._finishing
                and self._bucket_appended_bytes[bucket]
                >= compaction_threshold
            ):
                compact_buckets.append(bucket)

        self.bytes_written += written_entries * (8 + self.value_dtype.itemsize)

        seconds = time.perf_counter() - flush_started
        if self._finishing:
            self.finalize_flush_seconds += seconds
        else:
            self.online_flush_seconds += seconds

        for bucket in compact_buckets:
            self._compact_bucket(bucket)

    def _run_memmaps(self, run: SparseRun) -> tuple[np.memmap, np.memmap]:
        keys = np.memmap(
            run.key_path,
            mode="r",
            dtype=np.int64,
            offset=int(run.offset) * np.dtype(np.int64).itemsize,
            shape=(int(run.count),),
        )
        values = np.memmap(
            run.value_path,
            mode="r",
            dtype=self.value_dtype,
            offset=int(run.offset) * self.value_dtype.itemsize,
            shape=(int(run.count),),
        )
        return keys, values

    def _canonicalize_runs(
        self, runs: list[SparseRun]
    ) -> tuple[np.ndarray, np.ndarray]:
        """Load runs into owned arrays and close every mapping immediately.

        Keeping ``np.memmap`` instances alive while canonicalizing thousands of
        partitions can leave gigabytes of file-backed pages resident in a
        long-lived Jupyter kernel. Copy each run into one bounded pair of owned
        arrays and explicitly close the mapping before sorting or reducing.
        """
        total = sum(int(run.count) for run in runs)
        keys = np.empty(total, dtype=np.int64)
        values = np.empty(total, dtype=self.value_dtype)
        cursor = 0
        for run in runs:
            count = int(run.count)
            stop = cursor + count
            keys_map, values_map = self._run_memmaps(run)
            try:
                keys[cursor:stop] = keys_map
                values[cursor:stop] = values_map
            finally:
                keys_map._mmap.close()
                values_map._mmap.close()
            cursor = stop

        return self._canonicalize_arrays(keys, values)

    def _canonicalize_arrays(
        self, keys: np.ndarray, values: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Sort and aggregate one owned key/value pair."""
        if not keys.size:
            return keys, values
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
        return keys_agg[keep], values_agg[keep]

    def _copy_run_to_streams(
        self,
        run: SparseRun,
        key_stream: Any,
        value_stream: Any,
        *,
        window_entries: int = 1_000_000,
        source_streams: dict[Path, Any] | None = None,
    ) -> None:
        """Copy one disk-backed run without materializing it in memory."""
        streams = source_streams if source_streams is not None else {}
        window = max(1, int(window_entries))
        for path, itemsize, destination in (
            (run.key_path, np.dtype(np.int64).itemsize, key_stream),
            (run.value_path, self.value_dtype.itemsize, value_stream),
        ):
            source = streams.get(path)
            close_source = source is None and source_streams is None
            if source is None:
                source = path.open("rb")
                streams[path] = source
            source.seek(int(run.offset) * itemsize)
            remaining = int(run.count) * itemsize
            window_bytes = window * itemsize
            while remaining:
                payload = source.read(min(remaining, window_bytes))
                if not payload:
                    raise EOFError(f"Unexpected end of inventory shard: {path}")
                destination.write(payload)
                remaining -= len(payload)
            if close_source:
                source.close()

    def _compact_bucket(self, bucket: int) -> None:
        """Canonicalize one shard bucket and atomically replace its old files."""
        bucket_i = int(bucket)
        old_paths = self._bucket_paths.get(bucket_i)
        if old_paths is None:
            return
        items = [
            (partition, runs)
            for partition, runs in self._pending_runs.items()
            if self._partition_bucket(partition) == bucket_i and runs
        ]
        if not items:
            return

        compaction_started = time.perf_counter()
        total_entries = sum(run.count for _, runs in items for run in runs)
        self._ensure_disk_capacity(
            total_entries * (8 + self.value_dtype.itemsize),
            f"compact inventory bucket {bucket_i}",
        )
        new_key_path, new_value_path = self._new_run_paths(f"bucket-{bucket_i}-compact")
        replacements: dict[tuple[int, ...], SparseRun] = {}
        offset = 0
        try:
            with (
                new_key_path.open("wb") as key_stream,
                new_value_path.open("wb") as value_stream,
            ):
                for partition, runs in items:
                    if len(runs) == 1:
                        count = int(runs[0].count)
                        self._copy_run_to_streams(runs[0], key_stream, value_stream)
                    elif (
                        estimate_flush_peak_bytes(
                            sum(run.count for run in runs),
                            value_dtype=self.value_dtype,
                        )
                        <= self.memory_budget
                    ):
                        keys, values = self._canonicalize_runs(runs)
                        count = int(keys.size)
                        if count:
                            keys.tofile(key_stream)
                            values.tofile(value_stream)
                        del keys, values
                    else:
                        heap: list[tuple[int, int, SparseRun]] = [
                            (int(run.count), serial, run)
                            for serial, run in enumerate(runs)
                        ]
                        heapq.heapify(heap)
                        serial = len(heap)
                        while len(heap) > 1:
                            _, _, left = heapq.heappop(heap)
                            _, _, right = heapq.heappop(heap)
                            merged = self._merge_runs(partition, left, right)
                            heapq.heappush(heap, (int(merged.count), serial, merged))
                            serial += 1
                        merged = heap[0][2]
                        count = int(merged.count)
                        if count:
                            self._copy_run_to_streams(merged, key_stream, value_stream)
                        if merged.exclusive:
                            merged.key_path.unlink(missing_ok=True)
                            merged.value_path.unlink(missing_ok=True)

                    if count:
                        replacements[partition] = SparseRun(
                            new_key_path,
                            new_value_path,
                            count,
                            count,
                            0,
                            offset=offset,
                            exclusive=False,
                        )
                        offset += count
        except Exception:
            new_key_path.unlink(missing_ok=True)
            new_value_path.unlink(missing_ok=True)
            raise

        for partition, _ in items:
            replacement = replacements.get(partition)
            if replacement is None:
                self._pending_runs.pop(partition, None)
            else:
                self._pending_runs[partition] = [replacement]

        for path in old_paths:
            path.unlink(missing_ok=True)
            self._input_shard_paths.discard(path)
        if offset:
            self._bucket_paths[bucket_i] = (new_key_path, new_value_path)
            self._bucket_entries[bucket_i] = offset
            self._bucket_appended_bytes[bucket_i] = 0
            self._input_shard_paths.update((new_key_path, new_value_path))
            self.bytes_written += offset * (8 + self.value_dtype.itemsize)
            self.bytes_merged += total_entries * (8 + self.value_dtype.itemsize)
        else:
            new_key_path.unlink(missing_ok=True)
            new_value_path.unlink(missing_ok=True)
            self._bucket_paths.pop(bucket_i, None)
            self._bucket_entries.pop(bucket_i, None)
            self._bucket_appended_bytes.pop(bucket_i, None)

        seconds = time.perf_counter() - compaction_started
        if self._finishing:
            self.finalize_merge_seconds += seconds
            self.finalize_compaction_count += 1
        else:
            self.online_merge_seconds += seconds
            self.online_compaction_count += 1

    def _merge_runs(
        self,
        partition: tuple[int, ...],
        left: SparseRun,
        right: SparseRun,
    ) -> SparseRun:
        merge_started = time.perf_counter()
        level = max(left.level, right.level) + 1
        label = "merge-" + "-".join(str(item) for item in partition)
        key_path, value_path = self._new_run_paths(label)
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
        out_keys = np.memmap(key_path, mode="w+", dtype=np.int64, shape=(capacity,))
        out_values = np.memmap(
            value_path, mode="w+", dtype=self.value_dtype, shape=(capacity,)
        )
        left_keys_map, left_values_map = self._run_memmaps(left)
        right_keys_map, right_values_map = self._run_memmaps(right)
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
        for mapping in (
            out_keys,
            out_values,
            left_keys_map,
            left_values_map,
            right_keys_map,
            right_values_map,
        ):
            mapping._mmap.close()
        self.bytes_merged += capacity * (8 + self.value_dtype.itemsize)
        with key_path.open("r+b") as stream:
            stream.truncate(count * np.dtype(np.int64).itemsize)
        with value_path.open("r+b") as stream:
            stream.truncate(count * self.value_dtype.itemsize)
        for run in (left, right):
            if run.exclusive:
                run.key_path.unlink(missing_ok=True)
                run.value_path.unlink(missing_ok=True)
        result = SparseRun(
            key_path,
            value_path,
            count,
            count,
            level,
            offset=0,
            exclusive=True,
        )
        seconds = time.perf_counter() - merge_started
        if self._finishing:
            self.finalize_merge_seconds += seconds
        else:
            self.online_merge_seconds += seconds
        return result

    def _finish_runs(self, *, show_progress: bool = False) -> None:
        self._finishing = True
        if self._buffers:
            self._flush_partitions_to_shard(list(self._buffers))

        buckets = sorted(self._bucket_paths)
        bucket_iter = tqdm(
            buckets,
            desc="TRAILS LCI [3/4] Compact inventory buckets",
            unit="bucket",
            leave=True,
            disable=not show_progress,
        )
        for bucket in bucket_iter:
            needs_compaction = any(
                len(runs) > 1
                for partition, runs in self._pending_runs.items()
                if self._partition_bucket(partition) == bucket
            )
            if needs_compaction:
                self._compact_bucket(bucket)

        self._final_runs = {
            partition: runs[0] for partition, runs in self._pending_runs.items() if runs
        }
        self._pending_runs.clear()
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
            "version": 2,
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
                    "offset": block.run.offset,
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
        blocks: list[InventoryBlock],
        *,
        shape: tuple[int, ...],
        root_width: int,
        year_start: int,
    ) -> da.Array:
        if not blocks:
            delayed_value = dask.delayed(_empty_sparse_block)(
                shape, self.value_dtype.str
            )
        else:
            total_entries = sum(int(block.run.count) for block in blocks)
            decode_peak = estimate_decode_peak_bytes(
                total_entries,
                ndim=len(shape),
                value_dtype=self.value_dtype,
                coordinate_dtype=np.int64,
            )
            if decode_peak > self.memory_budget:
                raise MemoryError(
                    "A finalized inventory block is too large to decode within "
                    f"the configured budget ({decode_peak:,} > "
                    f"{self.memory_budget:,} bytes). Reduce activity_block_size "
                    "root_block_size, or year_block_size when constructing the builder."
                )
            run_specs = tuple(
                (
                    str(block.run.key_path),
                    str(block.run.value_path),
                    int(block.run.offset),
                    int(block.run.count),
                    int(block.year_index) - int(year_start),
                )
                for block in sorted(blocks, key=lambda item: item.year_index)
            )
            delayed_value = dask.delayed(_load_sparse_block_group)(
                run_specs,
                shape,
                has_root=self.has_root,
                n_flows=self.n_flows,
                root_width=int(root_width),
                year_stride=self.year_block_size,
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
        activity_chunks = _chunk_lengths(self.n_activities, self.activity_block_size)
        activity_starts = _chunk_starts(activity_chunks)
        year_chunks = _chunk_lengths(self.n_years, self.year_block_size)
        year_starts = _chunk_starts(year_chunks)

        if self.has_root:
            root_chunks = _chunk_lengths(self.n_activities, self.root_block_size)
            root_starts = _chunk_starts(root_chunks)
            activity_arrays: list[da.Array] = []
            activity_iter = tqdm(
                list(zip(activity_starts, activity_chunks)),
                desc="TRAILS LCI [4/4] Assemble inventory blocks",
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
                        grouped_blocks = [
                            block_map[(activity_block, year_index, root_block)]
                            for year_index in range(year_start, year_start + year_width)
                            if (activity_block, year_index, root_block) in block_map
                        ]
                        shape = (
                            activity_width,
                            self.n_flows,
                            year_width,
                            root_width,
                        )
                        root_arrays.append(
                            self._delayed_block(
                                grouped_blocks,
                                shape=shape,
                                root_width=self.root_block_size,
                                year_start=year_start,
                            )
                        )
                    year_arrays.append(da.concatenate(root_arrays, axis=3))
                activity_arrays.append(da.concatenate(year_arrays, axis=2))
            return da.concatenate(activity_arrays, axis=0)

        activity_arrays = []
        activity_iter = tqdm(
            list(zip(activity_starts, activity_chunks)),
            desc="TRAILS LCI [4/4] Assemble inventory blocks",
            unit="activity-block",
            leave=True,
            disable=not show_progress,
        )
        for activity_start, activity_width in activity_iter:
            activity_block = activity_start // self.activity_block_size
            year_arrays = []
            for year_start, year_width in zip(year_starts, year_chunks):
                grouped_blocks = [
                    block_map[(activity_block, year_index, None)]
                    for year_index in range(year_start, year_start + year_width)
                    if (activity_block, year_index, None) in block_map
                ]
                shape = (activity_width, self.n_flows, year_width)
                year_arrays.append(
                    self._delayed_block(
                        grouped_blocks,
                        shape=shape,
                        root_width=1,
                        year_start=year_start,
                    )
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
        self.dask_block_count = int(np.prod(result.numblocks, dtype=np.int64))
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
            "online_compaction_count": int(self.online_compaction_count),
            "finalize_compaction_count": int(self.finalize_compaction_count),
            "dask_block_count": int(self.dask_block_count),
            "logical_partition_count": int(len(self._final_runs)),
            "storage_file_count": int(
                sum(1 for path in self.store_path.iterdir() if path.is_file())
            ),
        }

    def reduce_activity_for_flows(
        self,
        flow_indices: Iterable[int],
        *,
        show_progress: bool = False,
        window_entries: int = 1_000_000,
    ) -> sparse.COO:
        """Stream selected flows into a compact flow/year/root inventory.

        This avoids constructing every activity-resolved Dask block when a
        downstream model, such as FaIR, needs only a small subset of elementary
        flows reduced over activities.
        """
        if not self._finalized:
            raise RuntimeError("Inventory must be finalized before reduction")
        if not self.has_root:
            raise ValueError("Flow reduction requires a root-attributed inventory")

        selected = tuple(
            int(value)
            for value in np.unique(np.asarray(list(flow_indices), dtype=np.int64))
            if 0 <= int(value) < self.n_flows
        )
        shape = (self.n_flows, self.n_years, self.n_activities)
        if not selected:
            return sparse.zeros(shape, dtype=np.float64)
        cached = self._activity_reduction_cache.get(selected)
        if cached is not None:
            return cached

        selected_mask = np.zeros(self.n_flows, dtype=bool)
        selected_mask[np.asarray(selected, dtype=np.int64)] = True
        totals: dict[int, float] = {}
        items = sorted(
            self._final_runs.items(),
            key=lambda item: (
                str(item[1].key_path),
                int(item[1].offset),
            ),
        )
        total_entries = sum(int(run.count) for _, run in items)
        progress = tqdm(
            total=total_entries,
            desc="Prepare FaIR inventory",
            unit="entry",
            unit_scale=True,
            leave=True,
            disable=not show_progress,
        )
        window = max(1, int(window_entries))
        try:
            for partition, run in items:
                keys_map, values_map = self._run_memmaps(run)
                try:
                    for start in range(0, int(run.count), window):
                        stop = min(start + window, int(run.count))
                        q = np.asarray(keys_map[start:stop], dtype=np.int64).copy()
                        values = np.asarray(
                            values_map[start:stop], dtype=np.float64
                        ).copy()

                        root_local = np.remainder(q, self.root_block_size)
                        np.floor_divide(q, self.root_block_size, out=q)
                        year_local = np.remainder(q, self.year_block_size)
                        np.floor_divide(q, self.year_block_size, out=q)
                        flows = np.remainder(q, self.n_flows)
                        keep = selected_mask[flows]
                        if np.any(keep):
                            flows_kept = flows[keep]
                            years_kept = int(partition[1]) + year_local[keep]
                            roots_kept = (
                                int(partition[2]) * self.root_block_size
                                + root_local[keep]
                            )
                            values_kept = values[keep]
                            linear = (
                                flows_kept * self.n_years + years_kept
                            ) * self.n_activities + roots_kept
                            order = np.argsort(linear, kind="quicksort")
                            linear = linear[order]
                            values_kept = values_kept[order]
                            first = np.empty(linear.size, dtype=bool)
                            first[0] = True
                            first[1:] = linear[1:] != linear[:-1]
                            starts = np.flatnonzero(first)
                            reduced_values = np.add.reduceat(values_kept, starts)
                            for key, value in zip(linear[starts], reduced_values):
                                key_i = int(key)
                                totals[key_i] = totals.get(key_i, 0.0) + float(value)
                        progress.update(stop - start)
                finally:
                    keys_map._mmap.close()
                    values_map._mmap.close()
        finally:
            progress.close()

        if not totals:
            result = sparse.zeros(shape, dtype=np.float64)
        else:
            linear = np.fromiter(totals.keys(), dtype=np.int64)
            values = np.fromiter(totals.values(), dtype=np.float64)
            keep = values != 0.0
            linear = linear[keep]
            values = values[keep]
            q = linear.copy()
            coords = np.empty((3, linear.size), dtype=np.int64)
            np.remainder(q, self.n_activities, out=coords[2])
            np.floor_divide(q, self.n_activities, out=q)
            np.remainder(q, self.n_years, out=coords[1])
            np.floor_divide(q, self.n_years, out=coords[0])
            result = sparse.COO(coords, values, shape=shape)

        self._activity_reduction_cache[selected] = result
        return result

    def iter_entries_for_flows(
        self,
        flow_indices: Iterable[int],
        *,
        show_progress: bool = False,
        progress_desc: str = "Scan inventory",
        progress: Any | None = None,
        window_entries: int = 1_000_000,
    ) -> Iterable[
        tuple[
            int,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
        ]
    ]:
        """Yield selected finalized entries in bounded, year-ordered windows.

        Each item contains ``(year_block_start, activities, flows, years,
        roots, values)`` with global integer coordinates. Only the requested
        flows are copied out of the disk-backed runs.
        """
        if not self._finalized:
            raise RuntimeError("Inventory must be finalized before streaming")
        if not self.has_root:
            raise ValueError("Selected entry streaming requires root attribution")

        selected = np.unique(np.asarray(list(flow_indices), dtype=np.int64))
        selected = selected[(selected >= 0) & (selected < self.n_flows)]
        if not selected.size:
            return
        selected_mask = np.zeros(self.n_flows, dtype=bool)
        selected_mask[selected] = True
        items = sorted(
            self._final_runs.items(),
            key=lambda item: (
                int(item[0][1]),
                str(item[1].key_path),
                int(item[1].offset),
            ),
        )
        owns_progress = progress is None
        progress_bar = progress
        if progress_bar is None:
            progress_bar = tqdm(
                total=sum(int(run.count) for _, run in items),
                desc=str(progress_desc),
                unit="entry",
                unit_scale=True,
                leave=True,
                disable=not show_progress,
            )
        window = max(1, int(window_entries))
        try:
            for partition, run in items:
                keys_map, values_map = self._run_memmaps(run)
                try:
                    for start in range(0, int(run.count), window):
                        stop = min(start + window, int(run.count))
                        q = np.asarray(keys_map[start:stop], dtype=np.int64).copy()
                        roots = np.remainder(q, self.root_block_size)
                        np.floor_divide(q, self.root_block_size, out=q)
                        years = np.remainder(q, self.year_block_size)
                        np.floor_divide(q, self.year_block_size, out=q)
                        flows = np.remainder(q, self.n_flows)
                        np.floor_divide(q, self.n_flows, out=q)
                        keep = selected_mask[flows]
                        progress_bar.update(stop - start)
                        if not np.any(keep):
                            continue
                        activities = (
                            int(partition[0]) * self.activity_block_size + q[keep]
                        )
                        years_global = int(partition[1]) + years[keep]
                        roots_global = (
                            int(partition[2]) * self.root_block_size + roots[keep]
                        )
                        values = np.asarray(values_map[start:stop])[keep].copy()
                        yield (
                            int(partition[1]),
                            activities.astype(np.int64, copy=False),
                            flows[keep].astype(np.int64, copy=False),
                            years_global.astype(np.int64, copy=False),
                            roots_global.astype(np.int64, copy=False),
                            values.astype(self.value_dtype, copy=False),
                        )
                finally:
                    keys_map._mmap.close()
                    values_map._mmap.close()
        finally:
            if owns_progress:
                progress_bar.close()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._buffers.clear()
        self._buffer_entries.clear()
        self._buffered_total_entries = 0
        self._pending_runs.clear()
        self._final_runs.clear()
        self._input_shard_paths.clear()
        self._bucket_paths.clear()
        self._bucket_entries.clear()
        self._bucket_appended_bytes.clear()
        self._activity_reduction_cache.clear()
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
            key=lambda index: (index[axis],)
            + tuple(value for i, value in enumerate(index) if i != axis)
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
