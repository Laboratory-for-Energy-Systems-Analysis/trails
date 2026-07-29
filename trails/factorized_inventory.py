"""Lazy factorized storage for root-attributed temporal inventories."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Any, Iterable
import uuid

import dask
import dask.array as da
import numpy as np
import sparse
from tqdm import tqdm

from .chunked_inventory import (
    ChunkedInventoryBuilder,
    DEFAULT_ACTIVITY_BLOCK_SIZE,
    DEFAULT_INVENTORY_MEMORY_BUDGET,
    DEFAULT_ROOT_BLOCK_SIZE,
    DEFAULT_YEAR_BLOCK_SIZE,
    MIN_FREE_DISK_RESERVE,
)


@dataclass(frozen=True)
class FactorizedInventoryRecord:
    """Disk-backed annual ``B * supply`` factor metadata."""

    year_index: int
    activities_path: Path
    flows_path: Path
    biosphere_values_path: Path
    supply_path: Path
    roots_path: Path
    biosphere_entries: int
    supply_shape: tuple[int, int]
    candidate_entries: int

    def delayed_spec(self) -> tuple[int, str, str, str, str, str]:
        """Return a scheduler-safe representation containing plain paths."""
        return (
            int(self.year_index),
            str(self.activities_path),
            str(self.flows_path),
            str(self.biosphere_values_path),
            str(self.supply_path),
            str(self.roots_path),
        )


def _chunk_lengths(size: int, width: int) -> tuple[int, ...]:
    """Return Dask chunk lengths for one dimension."""
    full, remainder = divmod(int(size), int(width))
    chunks = [int(width)] * full
    if remainder:
        chunks.append(int(remainder))
    return tuple(chunks) or (0,)


def _chunk_starts(chunks: Iterable[int]) -> tuple[int, ...]:
    """Return cumulative starts for a sequence of chunk lengths."""
    starts: list[int] = []
    cursor = 0
    for length in chunks:
        starts.append(cursor)
        cursor += int(length)
    return tuple(starts)


def _empty_sparse_block(shape: tuple[int, ...], dtype: str) -> sparse.COO:
    """Construct one empty sparse block for a delayed task."""
    return sparse.zeros(shape, dtype=np.dtype(dtype))


def _load_factorized_block(
    specs: tuple[tuple[int, str, str, str, str, str], ...],
    shape: tuple[int, int, int, int],
    *,
    activity_start: int,
    year_start: int,
    root_start: int,
    value_dtype: str,
    candidate_limit: int,
) -> sparse.COO:
    """Materialize one bounded inventory block from annual factors."""
    coords_parts: list[np.ndarray] = []
    value_parts: list[np.ndarray] = []
    activity_stop = int(activity_start) + int(shape[0])
    year_stop = int(year_start) + int(shape[2])
    root_stop = int(root_start) + int(shape[3])

    for (
        year_index,
        activities_path,
        flows_path,
        biosphere_values_path,
        supply_path,
        roots_path,
    ) in specs:
        if not year_start <= int(year_index) < year_stop:
            continue
        activities = np.load(activities_path, mmap_mode="r")
        in_activity_block = (activities >= activity_start) & (
            activities < activity_stop
        )
        if not np.any(in_activity_block):
            continue

        roots = np.load(roots_path, mmap_mode="r")
        root_columns = np.flatnonzero((roots >= root_start) & (roots < root_stop))
        if not root_columns.size:
            continue

        selected_rows = np.flatnonzero(in_activity_block)
        flows = np.load(flows_path, mmap_mode="r")
        biosphere_values = np.load(biosphere_values_path, mmap_mode="r")
        supply = np.load(supply_path, mmap_mode="r")
        rows_per_chunk = max(
            1,
            int(candidate_limit) // max(1, int(root_columns.size)),
        )

        for cursor in range(0, int(selected_rows.size), rows_per_chunk):
            row_indices = selected_rows[cursor : cursor + rows_per_chunk]
            activity_indices = np.asarray(activities[row_indices], dtype=np.int64)
            supply_values = np.asarray(
                supply[np.ix_(activity_indices, root_columns)],
                dtype=np.float64,
            )
            values = (
                np.asarray(biosphere_values[row_indices], dtype=np.float64)[:, None]
                * supply_values
            )
            bio_rows, root_offsets = np.nonzero(values)
            if not bio_rows.size:
                continue

            count = int(bio_rows.size)
            coords = np.empty((4, count), dtype=np.int64)
            coords[0] = activity_indices[bio_rows] - int(activity_start)
            coords[1] = np.asarray(flows[row_indices], dtype=np.int64)[bio_rows]
            coords[2].fill(int(year_index) - int(year_start))
            coords[3] = np.asarray(roots[root_columns], dtype=np.int64)[
                root_offsets
            ] - int(root_start)
            coords_parts.append(coords)
            value_parts.append(
                values[bio_rows, root_offsets].astype(np.dtype(value_dtype), copy=False)
            )

    if not value_parts:
        return sparse.zeros(shape, dtype=np.dtype(value_dtype))
    return sparse.COO(
        np.concatenate(coords_parts, axis=1),
        np.concatenate(value_parts),
        shape=shape,
        has_duplicates=True,
        sorted=False,
        idx_dtype=np.int64,
    )


class FactorizedInventoryBuilder:
    """Store annual supply factors and only explicit temporal corrections.

    The ordinary inventory identity ``inventory = B * supply`` remains lazy.
    Activities with temporal biosphere exchanges, plus direct biosphere inputs,
    are accumulated by a bounded :class:`ChunkedInventoryBuilder` and added as
    a sparse correction array during finalization.
    """

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
        if not has_root:
            raise ValueError(
                "The factorized inventory backend requires root attribution"
            )
        self.n_activities = int(n_activities)
        self.n_flows = int(n_flows)
        self.n_years = int(n_years)
        self.has_root = True
        self.value_dtype = np.dtype(value_dtype)
        self.memory_budget = int(memory_budget)
        if self.memory_budget <= 0:
            raise ValueError("inventory_memory_budget must be a positive integer")
        self.activity_block_size = max(1, int(activity_block_size))
        self.root_block_size = max(1, int(root_block_size))
        self.year_block_size = max(1, int(year_block_size))

        self.owned_store = store is None
        if store is None:
            self.store_path = Path(
                tempfile.mkdtemp(prefix="trails-factorized-inventory-")
            )
        else:
            root = Path(store).expanduser().resolve()
            root.mkdir(parents=True, exist_ok=True)
            self.store_path = root / f"factorized-run-{uuid.uuid4().hex}"
            self.store_path.mkdir(parents=True, exist_ok=False)

        self._corrections = ChunkedInventoryBuilder(
            n_activities=self.n_activities,
            n_flows=self.n_flows,
            n_years=self.n_years,
            has_root=True,
            value_dtype=self.value_dtype,
            memory_budget=self.memory_budget,
            store=self.store_path / "corrections",
            activity_block_size=self.activity_block_size,
            root_block_size=self.root_block_size,
            year_block_size=self.year_block_size,
        )
        self._records: list[FactorizedInventoryRecord] = []
        self._closed = False
        self._finalized = False
        self._final_array: da.Array | None = None
        self.factor_candidate_entries = 0
        self.factor_bytes_written = 0
        self.factor_write_seconds = 0.0
        self.dask_block_count = 0
        self._activity_reduction_cache: dict[tuple[int, ...], sparse.COO] = {}

    def _ensure_disk_capacity(self, additional_bytes: int) -> None:
        free = int(shutil.disk_usage(self.store_path).free)
        reserve = max(MIN_FREE_DISK_RESERVE, 2 * self.memory_budget)
        required = max(0, int(additional_bytes))
        if free < required + reserve:
            raise OSError(
                "Not enough free space to write a factorized inventory: "
                f"{required / 2**30:.2f} GiB additional space is required, "
                f"with a {reserve / 2**30:.2f} GiB safety reserve, but only "
                f"{free / 2**30:.2f} GiB is free in {self.store_path}."
            )

    def append(
        self,
        activities: np.ndarray,
        flows: np.ndarray,
        year_indices: np.ndarray,
        values: np.ndarray,
        *,
        roots: np.ndarray | None = None,
    ) -> None:
        """Append an explicit temporal or direct-biosphere correction."""
        self._corrections.append(
            activities,
            flows,
            year_indices,
            values,
            roots=roots,
        )

    def append_factor(
        self,
        *,
        year_index: int,
        activities: np.ndarray,
        flows: np.ndarray,
        biosphere_values: np.ndarray,
        supply_matrix: np.ndarray,
        roots: np.ndarray,
        excluded_activities: set[int] | None = None,
    ) -> None:
        """Persist one annual non-temporal ``B * supply`` factor."""
        if self._closed or self._finalized:
            raise RuntimeError("Cannot append to a finalized inventory builder")
        started = time.perf_counter()
        acts = np.asarray(activities, dtype=np.int64)
        flow_arr = np.asarray(flows, dtype=np.int64)
        bio_values = np.asarray(biosphere_values)
        supply = np.asarray(supply_matrix, dtype=np.float64)
        root_arr = np.asarray(roots, dtype=np.int64)
        if not (acts.size == flow_arr.size == bio_values.size):
            raise ValueError("Biosphere factor arrays must align")
        if supply.ndim != 2 or supply.shape[0] != self.n_activities:
            raise ValueError("supply_matrix must have shape (n_activities, n_roots)")
        if root_arr.ndim != 1 or root_arr.size != supply.shape[1]:
            raise ValueError("roots must align with supply_matrix columns")

        keep = bio_values != 0
        excluded = excluded_activities or set()
        if excluded:
            excluded_arr = np.fromiter(excluded, dtype=np.int64, count=len(excluded))
            keep &= ~np.isin(acts, excluded_arr)
        acts = acts[keep]
        flow_arr = flow_arr[keep]
        bio_values = bio_values[keep]
        if not acts.size or not root_arr.size:
            return

        nonzero_supply = np.count_nonzero(supply, axis=1).astype(np.int64, copy=False)
        biosphere_rows = np.bincount(acts, minlength=self.n_activities).astype(
            np.int64, copy=False
        )
        candidate_entries = int(np.dot(nonzero_supply, biosphere_rows))
        if candidate_entries == 0:
            return

        serial = len(self._records)
        stem = self.store_path / f"factor-{serial:04d}"
        paths = {
            "activities": stem.with_name(stem.name + "-activities.npy"),
            "flows": stem.with_name(stem.name + "-flows.npy"),
            "biosphere_values": stem.with_name(stem.name + "-biosphere-values.npy"),
            "supply": stem.with_name(stem.name + "-supply.npy"),
            "roots": stem.with_name(stem.name + "-roots.npy"),
        }
        arrays = {
            "activities": acts,
            "flows": flow_arr,
            "biosphere_values": bio_values,
            "supply": supply,
            "roots": root_arr,
        }
        estimated_bytes = sum(int(array.nbytes) + 256 for array in arrays.values())
        self._ensure_disk_capacity(estimated_bytes)
        written_paths: list[Path] = []
        try:
            for name, array in arrays.items():
                path = paths[name]
                np.save(path, array, allow_pickle=False)
                written_paths.append(path)
        except Exception:
            for path in written_paths:
                path.unlink(missing_ok=True)
            raise

        bytes_written = sum(path.stat().st_size for path in written_paths)
        self.factor_bytes_written += int(bytes_written)
        self.factor_candidate_entries += candidate_entries
        self._records.append(
            FactorizedInventoryRecord(
                year_index=int(year_index),
                activities_path=paths["activities"],
                flows_path=paths["flows"],
                biosphere_values_path=paths["biosphere_values"],
                supply_path=paths["supply"],
                roots_path=paths["roots"],
                biosphere_entries=int(acts.size),
                supply_shape=(int(supply.shape[0]), int(supply.shape[1])),
                candidate_entries=candidate_entries,
            )
        )
        self.factor_write_seconds += time.perf_counter() - started

    def _factor_array(self, *, show_progress: bool) -> da.Array:
        specs = tuple(record.delayed_spec() for record in self._records)
        activity_chunks = _chunk_lengths(self.n_activities, self.activity_block_size)
        year_chunks = _chunk_lengths(self.n_years, self.year_block_size)
        root_chunks = _chunk_lengths(self.n_activities, self.root_block_size)
        activity_starts = _chunk_starts(activity_chunks)
        year_starts = _chunk_starts(year_chunks)
        root_starts = _chunk_starts(root_chunks)
        candidate_limit = max(
            100_000,
            min(2_000_000, self.memory_budget // 64),
        )
        meta = sparse.zeros((0, 0, 0, 0), dtype=self.value_dtype)
        activity_arrays: list[da.Array] = []
        iterator = tqdm(
            list(zip(activity_starts, activity_chunks)),
            desc="TRAILS LCI [4/4] Assemble factorized inventory",
            unit="activity-block",
            leave=True,
            disable=not show_progress,
        )
        for activity_start, activity_width in iterator:
            year_arrays: list[da.Array] = []
            for year_start, year_width in zip(year_starts, year_chunks):
                relevant_specs = tuple(
                    spec
                    for spec in specs
                    if year_start <= int(spec[0]) < year_start + year_width
                )
                root_arrays: list[da.Array] = []
                for root_start, root_width in zip(root_starts, root_chunks):
                    shape = (
                        int(activity_width),
                        self.n_flows,
                        int(year_width),
                        int(root_width),
                    )
                    if relevant_specs:
                        delayed_value = dask.delayed(_load_factorized_block)(
                            relevant_specs,
                            shape,
                            activity_start=int(activity_start),
                            year_start=int(year_start),
                            root_start=int(root_start),
                            value_dtype=self.value_dtype.str,
                            candidate_limit=int(candidate_limit),
                        )
                    else:
                        delayed_value = dask.delayed(_empty_sparse_block)(
                            shape, self.value_dtype.str
                        )
                    root_arrays.append(
                        da.from_delayed(
                            delayed_value,
                            shape=shape,
                            dtype=self.value_dtype,
                            meta=meta,
                        )
                    )
                year_arrays.append(da.concatenate(root_arrays, axis=3))
            activity_arrays.append(da.concatenate(year_arrays, axis=2))
        return da.concatenate(activity_arrays, axis=0)

    def finalize(self, *, show_progress: bool = False) -> da.Array:
        """Build a lazy Dask inventory from factors plus sparse corrections."""
        if self._closed:
            raise RuntimeError("Inventory builder is closed")
        if self._finalized:
            raise RuntimeError("Inventory builder was already finalized")
        factors = self._factor_array(show_progress=show_progress)
        if self._corrections.raw_entries:
            corrections = self._corrections.finalize(show_progress=show_progress)
            result = factors + corrections
        else:
            result = factors
        self._final_array = result
        self.dask_block_count = int(np.prod(result.numblocks, dtype=np.int64))
        self._write_manifest()
        self._finalized = True
        return result

    def _write_manifest(self) -> None:
        payload = {
            "version": 1,
            "backend": "factorized",
            "shape": [
                self.n_activities,
                self.n_flows,
                self.n_years,
                self.n_activities,
            ],
            "value_dtype": self.value_dtype.str,
            "factor_candidate_entries": self.factor_candidate_entries,
            "records": [
                {
                    "year_index": record.year_index,
                    "activities": record.activities_path.name,
                    "flows": record.flows_path.name,
                    "biosphere_values": record.biosphere_values_path.name,
                    "supply": record.supply_path.name,
                    "roots": record.roots_path.name,
                    "biosphere_entries": record.biosphere_entries,
                    "supply_shape": list(record.supply_shape),
                    "candidate_entries": record.candidate_entries,
                }
                for record in self._records
            ],
        }
        temporary = self.store_path / "manifest.json.tmp"
        final = self.store_path / "manifest.json"
        temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        os.replace(temporary, final)

    @property
    def nnz(self) -> int:
        """Return a conservative upper bound used by materialization guards."""
        correction_nnz = int(self._corrections.nnz)
        return int(self.factor_candidate_entries + correction_nnz)

    def diagnostics(self) -> dict[str, int | float | str | bool | dict]:
        correction = self._corrections.diagnostics()
        return {
            "backend": "factorized",
            "store": str(self.store_path),
            "owned_store": self.owned_store,
            "factor_count": int(len(self._records)),
            "factor_candidate_entries": int(self.factor_candidate_entries),
            "factor_bytes_written": int(self.factor_bytes_written),
            "factor_write_seconds": float(self.factor_write_seconds),
            "explicit_correction_entries": int(
                correction.get("canonical_entries", correction.get("raw_entries", 0))
            ),
            "dask_block_count": int(self.dask_block_count),
            "corrections": correction,
        }

    def reduce_activity_for_flows(
        self,
        flow_indices: Iterable[int],
        *,
        show_progress: bool = False,
        window_entries: int = 1_000_000,
    ) -> sparse.COO:
        """Reduce selected flows without expanding the activity dimension."""
        if not self._finalized:
            raise RuntimeError("Inventory must be finalized before reduction")
        selected = tuple(
            int(value)
            for value in np.unique(np.asarray(list(flow_indices), dtype=np.int64))
            if 0 <= int(value) < self.n_flows
        )
        cached = self._activity_reduction_cache.get(selected)
        if cached is not None:
            return cached
        shape = (self.n_flows, self.n_years, self.n_activities)
        if not selected:
            result = sparse.zeros(shape, dtype=self.value_dtype)
            self._activity_reduction_cache[selected] = result
            return result

        selected_mask = np.zeros(self.n_flows, dtype=bool)
        selected_mask[np.asarray(selected, dtype=np.int64)] = True
        coords_parts: list[np.ndarray] = []
        value_parts: list[np.ndarray] = []
        iterator = tqdm(
            self._records,
            desc="Factorized inventory flow reduction",
            unit="year",
            leave=True,
            disable=not show_progress,
        )
        for record in iterator:
            flows = np.load(record.flows_path, mmap_mode="r")
            selected_rows = np.flatnonzero(selected_mask[flows])
            if not selected_rows.size:
                continue
            activities = np.load(record.activities_path, mmap_mode="r")
            biosphere_values = np.load(record.biosphere_values_path, mmap_mode="r")
            supply = np.load(record.supply_path, mmap_mode="r")
            roots = np.load(record.roots_path, mmap_mode="r")
            rows_per_chunk = max(
                1,
                int(window_entries) // max(1, int(roots.size)),
            )
            root_columns = np.arange(roots.size, dtype=np.int64)
            for cursor in range(0, int(selected_rows.size), rows_per_chunk):
                row_indices = selected_rows[cursor : cursor + rows_per_chunk]
                activity_indices = np.asarray(activities[row_indices], dtype=np.int64)
                values = np.asarray(biosphere_values[row_indices], dtype=np.float64)[
                    :, None
                ] * np.asarray(
                    supply[np.ix_(activity_indices, root_columns)],
                    dtype=np.float64,
                )
                bio_rows, root_offsets = np.nonzero(values)
                if not bio_rows.size:
                    continue
                count = int(bio_rows.size)
                coords = np.empty((3, count), dtype=np.int64)
                coords[0] = np.asarray(flows[row_indices], dtype=np.int64)[bio_rows]
                coords[1].fill(int(record.year_index))
                coords[2] = np.asarray(roots, dtype=np.int64)[root_offsets]
                coords_parts.append(coords)
                value_parts.append(
                    values[bio_rows, root_offsets].astype(self.value_dtype, copy=False)
                )

        if value_parts:
            result = sparse.COO(
                np.concatenate(coords_parts, axis=1),
                np.concatenate(value_parts),
                shape=shape,
                has_duplicates=True,
                sorted=False,
                idx_dtype=np.int64,
            )
        else:
            result = sparse.zeros(shape, dtype=self.value_dtype)
        if self._corrections.raw_entries:
            correction = self._corrections.reduce_activity_for_flows(
                selected,
                show_progress=show_progress,
                window_entries=window_entries,
            )
            result = result + correction
        self._activity_reduction_cache[selected] = result
        return result

    def iter_entries_for_flows(
        self,
        flow_indices: Iterable[int],
        *,
        show_progress: bool = False,
        progress_desc: str = "Factorized inventory",
        window_entries: int = 1_000_000,
        progress_bar: Any | None = None,
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
        """Yield selected explicit inventory entries directly from factors."""
        if not self._finalized:
            raise RuntimeError("Inventory must be finalized before iteration")
        selected = np.asarray(
            [
                int(value)
                for value in np.unique(np.asarray(list(flow_indices), dtype=np.int64))
                if 0 <= int(value) < self.n_flows
            ],
            dtype=np.int64,
        )
        if not selected.size:
            return
        selected_mask = np.zeros(self.n_flows, dtype=bool)
        selected_mask[selected] = True
        owns_progress = progress_bar is None
        bar = progress_bar
        if bar is None:
            bar = tqdm(
                total=len(self._records),
                desc=progress_desc,
                unit="year",
                leave=True,
                disable=not show_progress,
            )
        try:
            for record in self._records:
                flows = np.load(record.flows_path, mmap_mode="r")
                selected_rows = np.flatnonzero(selected_mask[flows])
                if selected_rows.size:
                    activities = np.load(record.activities_path, mmap_mode="r")
                    biosphere_values = np.load(
                        record.biosphere_values_path, mmap_mode="r"
                    )
                    supply = np.load(record.supply_path, mmap_mode="r")
                    roots = np.load(record.roots_path, mmap_mode="r")
                    root_columns = np.arange(roots.size, dtype=np.int64)
                    rows_per_chunk = max(
                        1,
                        int(window_entries) // max(1, int(roots.size)),
                    )
                    for cursor in range(0, int(selected_rows.size), rows_per_chunk):
                        row_indices = selected_rows[cursor : cursor + rows_per_chunk]
                        activity_indices = np.asarray(
                            activities[row_indices], dtype=np.int64
                        )
                        values = np.asarray(
                            biosphere_values[row_indices],
                            dtype=np.float64,
                        )[:, None] * np.asarray(
                            supply[np.ix_(activity_indices, root_columns)],
                            dtype=np.float64,
                        )
                        bio_rows, root_offsets = np.nonzero(values)
                        if not bio_rows.size:
                            continue
                        count = int(bio_rows.size)
                        years = np.full(count, int(record.year_index), dtype=np.int64)
                        yield (
                            int(record.year_index),
                            activity_indices[bio_rows],
                            np.asarray(flows[row_indices], dtype=np.int64)[bio_rows],
                            years,
                            np.asarray(roots, dtype=np.int64)[root_offsets],
                            values[bio_rows, root_offsets].astype(
                                self.value_dtype, copy=False
                            ),
                        )
                bar.update(1)
        finally:
            if owns_progress:
                bar.close()

        if self._corrections.raw_entries:
            yield from self._corrections.iter_entries_for_flows(
                selected,
                show_progress=show_progress,
                progress_desc=f"{progress_desc} corrections",
                window_entries=window_entries,
            )

    def close(self) -> None:
        """Close correction storage and remove a managed factor store."""
        if self._closed:
            return
        self._closed = True
        self._corrections.close()
        self._records.clear()
        self._final_array = None
        self._activity_reduction_cache.clear()
        if self.owned_store:
            shutil.rmtree(self.store_path, ignore_errors=True)

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass
