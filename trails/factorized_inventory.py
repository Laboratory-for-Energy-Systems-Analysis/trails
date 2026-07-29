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


@dataclass(frozen=True)
class TemporalFactorizedInventoryRecord:
    """Disk-backed ``B * supply * temporal kernel`` metadata."""

    base_year_index: int
    activities_path: Path
    flows_path: Path
    biosphere_values_path: Path
    supply_path: Path
    roots_path: Path
    pulse_indptr_path: Path
    pulse_year_indices_path: Path
    pulse_weights_path: Path
    biosphere_entries: int
    supply_shape: tuple[int, int]
    candidate_entries: int
    pulse_entries: int
    max_pulses: int
    min_amount: float
    relevant_year_indices: tuple[int, ...]

    def delayed_spec(
        self,
    ) -> tuple[int, str, str, str, str, str, str, str, str, float, int]:
        """Return a scheduler-safe representation containing plain paths."""
        return (
            int(self.base_year_index),
            str(self.activities_path),
            str(self.flows_path),
            str(self.biosphere_values_path),
            str(self.supply_path),
            str(self.roots_path),
            str(self.pulse_indptr_path),
            str(self.pulse_year_indices_path),
            str(self.pulse_weights_path),
            float(self.min_amount),
            int(self.max_pulses),
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


def _expand_temporal_candidates(
    *,
    row_indices: np.ndarray,
    activity_indices: np.ndarray,
    flow_indices: np.ndarray,
    values: np.ndarray,
    roots: np.ndarray,
    pulse_indptr: np.ndarray,
    pulse_year_indices: np.ndarray,
    pulse_weights: np.ndarray,
    base_year_index: int,
    min_amount: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Expand bounded nonzero ``B * supply`` candidates through TD kernels."""
    bio_rows, root_offsets = np.nonzero(values)
    if not bio_rows.size:
        empty_i = np.empty(0, dtype=np.int64)
        empty_v = np.empty(0, dtype=values.dtype)
        return empty_i, empty_i, empty_i, empty_i, empty_v

    candidate_values = values[bio_rows, root_offsets]
    output_activities: list[np.ndarray] = []
    output_flows: list[np.ndarray] = []
    output_years: list[np.ndarray] = []
    output_roots: list[np.ndarray] = []
    output_values: list[np.ndarray] = []

    if min_amount > 0.0:
        anchor = np.abs(candidate_values) < float(min_amount)
        if np.any(anchor):
            output_activities.append(activity_indices[bio_rows[anchor]])
            output_flows.append(flow_indices[bio_rows[anchor]])
            output_years.append(
                np.full(int(np.count_nonzero(anchor)), base_year_index, dtype=np.int64)
            )
            output_roots.append(roots[root_offsets[anchor]])
            output_values.append(candidate_values[anchor])
        temporal = ~anchor
    else:
        temporal = np.ones(candidate_values.size, dtype=bool)

    if np.any(temporal):
        temporal_bio_rows = bio_rows[temporal]
        starts = pulse_indptr[row_indices[temporal_bio_rows]]
        stops = pulse_indptr[row_indices[temporal_bio_rows] + 1]
        counts = stops - starts
        valid = counts > 0
        if np.any(valid):
            temporal_bio_rows = temporal_bio_rows[valid]
            temporal_root_offsets = root_offsets[temporal][valid]
            temporal_values = candidate_values[temporal][valid]
            starts = starts[valid]
            counts = counts[valid]
            expanded_candidates = np.repeat(
                np.arange(counts.size, dtype=np.int64), counts
            )
            group_starts = np.repeat(np.cumsum(counts) - counts, counts)
            pulse_positions = (
                np.arange(int(counts.sum()), dtype=np.int64)
                - group_starts
                + np.repeat(starts, counts)
            )
            output_activities.append(
                activity_indices[temporal_bio_rows][expanded_candidates]
            )
            output_flows.append(flow_indices[temporal_bio_rows][expanded_candidates])
            output_years.append(pulse_year_indices[pulse_positions])
            output_roots.append(roots[temporal_root_offsets][expanded_candidates])
            output_values.append(
                temporal_values[expanded_candidates] * pulse_weights[pulse_positions]
            )

    if not output_values:
        empty_i = np.empty(0, dtype=np.int64)
        empty_v = np.empty(0, dtype=values.dtype)
        return empty_i, empty_i, empty_i, empty_i, empty_v
    return (
        np.concatenate(output_activities),
        np.concatenate(output_flows),
        np.concatenate(output_years),
        np.concatenate(output_roots),
        np.concatenate(output_values),
    )


def _load_temporal_factorized_block(
    specs: tuple[tuple[int, str, str, str, str, str, str, str, str, float, int], ...],
    shape: tuple[int, int, int, int],
    *,
    activity_start: int,
    year_start: int,
    root_start: int,
    value_dtype: str,
    candidate_limit: int,
) -> sparse.COO:
    """Materialize one bounded block from compact ported temporal factors."""
    coords_parts: list[np.ndarray] = []
    value_parts: list[np.ndarray] = []
    activity_stop = int(activity_start) + int(shape[0])
    year_stop = int(year_start) + int(shape[2])
    root_stop = int(root_start) + int(shape[3])

    for (
        base_year_index,
        activities_path,
        flows_path,
        biosphere_values_path,
        supply_path,
        roots_path,
        pulse_indptr_path,
        pulse_year_indices_path,
        pulse_weights_path,
        min_amount,
        max_pulses,
    ) in specs:
        activities = np.load(activities_path, mmap_mode="r")
        selected_rows = np.flatnonzero(
            (activities >= activity_start) & (activities < activity_stop)
        )
        if not selected_rows.size:
            continue
        roots = np.load(roots_path, mmap_mode="r")
        root_columns = np.flatnonzero((roots >= root_start) & (roots < root_stop))
        if not root_columns.size:
            continue

        flows = np.load(flows_path, mmap_mode="r")
        biosphere_values = np.load(biosphere_values_path, mmap_mode="r")
        supply = np.load(supply_path, mmap_mode="r")
        pulse_indptr = np.load(pulse_indptr_path, mmap_mode="r")
        pulse_year_indices = np.load(pulse_year_indices_path, mmap_mode="r")
        pulse_weights = np.load(pulse_weights_path, mmap_mode="r")
        rows_per_chunk = max(
            1,
            int(candidate_limit)
            // max(1, int(root_columns.size) * max(1, int(max_pulses))),
        )
        for cursor in range(0, int(selected_rows.size), rows_per_chunk):
            row_indices = selected_rows[cursor : cursor + rows_per_chunk]
            activity_indices = np.asarray(activities[row_indices], dtype=np.int64)
            flow_indices = np.asarray(flows[row_indices], dtype=np.int64)
            values = np.asarray(biosphere_values[row_indices], dtype=np.float64)[
                :, None
            ] * np.asarray(
                supply[np.ix_(activity_indices, root_columns)], dtype=np.float64
            )
            out_acts, out_flows, out_years, out_roots, out_values = (
                _expand_temporal_candidates(
                    row_indices=row_indices,
                    activity_indices=activity_indices,
                    flow_indices=flow_indices,
                    values=values,
                    roots=np.asarray(roots[root_columns], dtype=np.int64),
                    pulse_indptr=pulse_indptr,
                    pulse_year_indices=pulse_year_indices,
                    pulse_weights=pulse_weights,
                    base_year_index=int(base_year_index),
                    min_amount=float(min_amount),
                )
            )
            keep = (out_years >= year_start) & (out_years < year_stop)
            if not np.any(keep):
                continue
            count = int(np.count_nonzero(keep))
            coords = np.empty((4, count), dtype=np.int64)
            coords[0] = out_acts[keep] - int(activity_start)
            coords[1] = out_flows[keep]
            coords[2] = out_years[keep] - int(year_start)
            coords[3] = out_roots[keep] - int(root_start)
            coords_parts.append(coords)
            value_parts.append(
                out_values[keep].astype(np.dtype(value_dtype), copy=False)
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


def _iter_temporal_record_entries(
    record: TemporalFactorizedInventoryRecord,
    selected_mask: np.ndarray,
    *,
    window_entries: int,
) -> Iterable[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Yield bounded selected-flow entries from one temporal factor record."""
    flows = np.load(record.flows_path, mmap_mode="r")
    selected_rows = np.flatnonzero(selected_mask[flows])
    if not selected_rows.size:
        return
    activities = np.load(record.activities_path, mmap_mode="r")
    biosphere_values = np.load(record.biosphere_values_path, mmap_mode="r")
    supply = np.load(record.supply_path, mmap_mode="r")
    roots = np.load(record.roots_path, mmap_mode="r")
    pulse_indptr = np.load(record.pulse_indptr_path, mmap_mode="r")
    pulse_year_indices = np.load(record.pulse_year_indices_path, mmap_mode="r")
    pulse_weights = np.load(record.pulse_weights_path, mmap_mode="r")
    root_columns = np.arange(roots.size, dtype=np.int64)
    rows_per_chunk = max(
        1,
        int(window_entries) // max(1, int(roots.size) * max(1, int(record.max_pulses))),
    )
    for cursor in range(0, int(selected_rows.size), rows_per_chunk):
        row_indices = selected_rows[cursor : cursor + rows_per_chunk]
        activity_indices = np.asarray(activities[row_indices], dtype=np.int64)
        flow_indices = np.asarray(flows[row_indices], dtype=np.int64)
        values = np.asarray(biosphere_values[row_indices], dtype=np.float64)[
            :, None
        ] * np.asarray(supply[np.ix_(activity_indices, root_columns)], dtype=np.float64)
        expanded = _expand_temporal_candidates(
            row_indices=row_indices,
            activity_indices=activity_indices,
            flow_indices=flow_indices,
            values=values,
            roots=np.asarray(roots, dtype=np.int64),
            pulse_indptr=pulse_indptr,
            pulse_year_indices=pulse_year_indices,
            pulse_weights=pulse_weights,
            base_year_index=int(record.base_year_index),
            min_amount=float(record.min_amount),
        )
        if expanded[-1].size:
            yield expanded


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
        root_block_size: int | None = None,
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
        # Factor records contain only the active root columns. A full-width
        # sparse root chunk therefore avoids thousands of empty Dask tasks
        # without creating a dense n_activity x n_root temporary.
        self.root_block_size = (
            self.n_activities
            if root_block_size is None
            else max(1, int(root_block_size))
        )
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
        self._temporal_records: list[TemporalFactorizedInventoryRecord] = []
        self._closed = False
        self._finalized = False
        self._final_array: da.Array | None = None
        self.factor_candidate_entries = 0
        self.factor_bytes_written = 0
        self.factor_write_seconds = 0.0
        self.temporal_factor_candidate_entries = 0
        self.temporal_factor_pulse_entries = 0
        self.temporal_factor_bytes_written = 0
        self.temporal_factor_write_seconds = 0.0
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
        included_rows: np.ndarray | None = None,
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
        if included_rows is not None:
            row_mask = np.asarray(included_rows, dtype=bool)
            if row_mask.shape != keep.shape:
                raise ValueError(
                    "included_rows must align with biosphere factor arrays"
                )
            keep &= row_mask
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

    def append_temporal_factor(
        self,
        *,
        base_year_index: int,
        activities: np.ndarray,
        flows: np.ndarray,
        biosphere_values: np.ndarray,
        supply_matrix: np.ndarray,
        roots: np.ndarray,
        pulse_indptr: np.ndarray,
        pulse_year_indices: np.ndarray,
        pulse_weights: np.ndarray,
        min_amount: float = 0.0,
    ) -> None:
        """Persist one annual ported temporal ``B * supply * kernel`` factor."""
        if self._closed or self._finalized:
            raise RuntimeError("Cannot append to a finalized inventory builder")
        started = time.perf_counter()
        acts = np.asarray(activities, dtype=np.int64)
        flow_arr = np.asarray(flows, dtype=np.int64)
        bio_values = np.asarray(biosphere_values)
        supply = np.asarray(supply_matrix, dtype=np.float64)
        root_arr = np.asarray(roots, dtype=np.int64)
        indptr = np.asarray(pulse_indptr, dtype=np.int64)
        year_indices = np.asarray(pulse_year_indices, dtype=np.int64)
        weights = np.asarray(pulse_weights, dtype=np.float64)
        if not (acts.size == flow_arr.size == bio_values.size):
            raise ValueError("Temporal biosphere factor arrays must align")
        if indptr.ndim != 1 or indptr.size != acts.size + 1:
            raise ValueError("pulse_indptr must contain one segment per exchange")
        if not year_indices.size == weights.size == int(indptr[-1]):
            raise ValueError("Temporal pulse arrays do not match pulse_indptr")
        if supply.ndim != 2 or supply.shape[0] != self.n_activities:
            raise ValueError("supply_matrix must have shape (n_activities, n_roots)")
        if root_arr.ndim != 1 or root_arr.size != supply.shape[1]:
            raise ValueError("roots must align with supply_matrix columns")
        if np.any(year_indices < 0) or np.any(year_indices >= self.n_years):
            raise ValueError("pulse_year_indices fall outside the inventory")

        keep = (bio_values != 0) & (np.diff(indptr) > 0)
        if not np.all(keep):
            kept_rows = np.flatnonzero(keep)
            if not kept_rows.size:
                return
            new_counts = np.diff(indptr)[kept_rows]
            pulse_positions = np.concatenate(
                [
                    np.arange(indptr[row], indptr[row + 1], dtype=np.int64)
                    for row in kept_rows
                ]
            )
            acts = acts[kept_rows]
            flow_arr = flow_arr[kept_rows]
            bio_values = bio_values[kept_rows]
            indptr = np.concatenate(
                [np.array([0], dtype=np.int64), np.cumsum(new_counts, dtype=np.int64)]
            )
            year_indices = year_indices[pulse_positions]
            weights = weights[pulse_positions]
        if not acts.size or not root_arr.size:
            return

        nonzero_supply = np.count_nonzero(supply, axis=1).astype(np.int64, copy=False)
        row_supply_counts = nonzero_supply[acts]
        pulse_counts = np.diff(indptr)
        candidate_entries = int(np.sum(row_supply_counts, dtype=np.int64))
        pulse_entries = int(np.sum(row_supply_counts * pulse_counts, dtype=np.int64))
        if candidate_entries == 0:
            return

        serial = len(self._temporal_records)
        stem = self.store_path / f"temporal-factor-{serial:04d}"
        paths = {
            "activities": stem.with_name(stem.name + "-activities.npy"),
            "flows": stem.with_name(stem.name + "-flows.npy"),
            "biosphere_values": stem.with_name(stem.name + "-biosphere-values.npy"),
            "supply": stem.with_name(stem.name + "-supply.npy"),
            "roots": stem.with_name(stem.name + "-roots.npy"),
            "pulse_indptr": stem.with_name(stem.name + "-pulse-indptr.npy"),
            "pulse_year_indices": stem.with_name(stem.name + "-pulse-year-indices.npy"),
            "pulse_weights": stem.with_name(stem.name + "-pulse-weights.npy"),
        }
        arrays = {
            "activities": acts,
            "flows": flow_arr,
            "biosphere_values": bio_values,
            "supply": supply,
            "roots": root_arr,
            "pulse_indptr": indptr,
            "pulse_year_indices": year_indices,
            "pulse_weights": weights,
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
        counts = np.diff(indptr)
        relevant_years = np.unique(
            np.concatenate(
                [year_indices, np.array([int(base_year_index)], dtype=np.int64)]
            )
        )
        self.temporal_factor_candidate_entries += candidate_entries
        self.temporal_factor_pulse_entries += pulse_entries
        self.temporal_factor_bytes_written += int(bytes_written)
        self._temporal_records.append(
            TemporalFactorizedInventoryRecord(
                base_year_index=int(base_year_index),
                activities_path=paths["activities"],
                flows_path=paths["flows"],
                biosphere_values_path=paths["biosphere_values"],
                supply_path=paths["supply"],
                roots_path=paths["roots"],
                pulse_indptr_path=paths["pulse_indptr"],
                pulse_year_indices_path=paths["pulse_year_indices"],
                pulse_weights_path=paths["pulse_weights"],
                biosphere_entries=int(acts.size),
                supply_shape=(int(supply.shape[0]), int(supply.shape[1])),
                candidate_entries=candidate_entries,
                pulse_entries=pulse_entries,
                max_pulses=int(counts.max(initial=0)),
                min_amount=float(min_amount),
                relevant_year_indices=tuple(int(value) for value in relevant_years),
            )
        )
        self.temporal_factor_write_seconds += time.perf_counter() - started

    def _factor_array(self, *, show_progress: bool) -> da.Array:
        specs = tuple(record.delayed_spec() for record in self._records)
        temporal_specs = tuple(
            (record.relevant_year_indices, record.delayed_spec())
            for record in self._temporal_records
        )
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
                relevant_temporal_specs = tuple(
                    spec
                    for years, spec in temporal_specs
                    if any(
                        year_start <= int(year) < year_start + year_width
                        for year in years
                    )
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
                        delayed_base = dask.delayed(_load_factorized_block)(
                            relevant_specs,
                            shape,
                            activity_start=int(activity_start),
                            year_start=int(year_start),
                            root_start=int(root_start),
                            value_dtype=self.value_dtype.str,
                            candidate_limit=int(candidate_limit),
                        )
                    else:
                        delayed_base = dask.delayed(_empty_sparse_block)(
                            shape, self.value_dtype.str
                        )
                    block = da.from_delayed(
                        delayed_base,
                        shape=shape,
                        dtype=self.value_dtype,
                        meta=meta,
                    )
                    if relevant_temporal_specs:
                        delayed_temporal = dask.delayed(
                            _load_temporal_factorized_block
                        )(
                            relevant_temporal_specs,
                            shape,
                            activity_start=int(activity_start),
                            year_start=int(year_start),
                            root_start=int(root_start),
                            value_dtype=self.value_dtype.str,
                            candidate_limit=int(candidate_limit),
                        )
                        block = block + da.from_delayed(
                            delayed_temporal,
                            shape=shape,
                            dtype=self.value_dtype,
                            meta=meta,
                        )
                    root_arrays.append(block)
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
            "temporal_factor_candidate_entries": (
                self.temporal_factor_candidate_entries
            ),
            "temporal_factor_pulse_entries": self.temporal_factor_pulse_entries,
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
            "temporal_records": [
                {
                    "base_year_index": record.base_year_index,
                    "activities": record.activities_path.name,
                    "flows": record.flows_path.name,
                    "biosphere_values": record.biosphere_values_path.name,
                    "supply": record.supply_path.name,
                    "roots": record.roots_path.name,
                    "pulse_indptr": record.pulse_indptr_path.name,
                    "pulse_year_indices": record.pulse_year_indices_path.name,
                    "pulse_weights": record.pulse_weights_path.name,
                    "biosphere_entries": record.biosphere_entries,
                    "supply_shape": list(record.supply_shape),
                    "candidate_entries": record.candidate_entries,
                    "pulse_entries": record.pulse_entries,
                    "max_pulses": record.max_pulses,
                    "min_amount": record.min_amount,
                    "relevant_year_indices": list(record.relevant_year_indices),
                }
                for record in self._temporal_records
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
        return int(
            self.factor_candidate_entries
            + self.temporal_factor_pulse_entries
            + correction_nnz
        )

    def diagnostics(self) -> dict[str, int | float | str | bool | dict]:
        correction = self._corrections.diagnostics()
        return {
            "backend": "factorized",
            "store": str(self.store_path),
            "owned_store": self.owned_store,
            "factor_count": int(len(self._records)),
            "temporal_factor_count": int(len(self._temporal_records)),
            "factor_candidate_entries": int(self.factor_candidate_entries),
            "temporal_factor_candidate_entries": int(
                self.temporal_factor_candidate_entries
            ),
            "temporal_factor_pulse_entries": int(self.temporal_factor_pulse_entries),
            "factor_bytes_written": int(self.factor_bytes_written),
            "temporal_factor_bytes_written": int(self.temporal_factor_bytes_written),
            "factor_write_seconds": float(self.factor_write_seconds),
            "temporal_factor_write_seconds": float(self.temporal_factor_write_seconds),
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

        temporal_iterator = tqdm(
            self._temporal_records,
            desc="Factorized temporal inventory flow reduction",
            unit="year",
            leave=True,
            disable=not show_progress,
        )
        for record in temporal_iterator:
            for (
                activities,
                flows,
                years,
                roots,
                values,
            ) in _iter_temporal_record_entries(
                record,
                selected_mask,
                window_entries=window_entries,
            ):
                coords = np.vstack([flows, years, roots]).astype(np.int64, copy=False)
                coords_parts.append(coords)
                value_parts.append(values.astype(self.value_dtype, copy=False))

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
        progress: Any | None = None,
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
        if progress_bar is not None and progress is not None:
            raise ValueError("Pass only one of progress_bar or progress")
        external_progress = progress_bar if progress_bar is not None else progress
        owns_progress = external_progress is None
        bar = external_progress
        if bar is None:
            bar = tqdm(
                total=len(self._records) + len(self._temporal_records),
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

            for record in self._temporal_records:
                for (
                    activities,
                    flows,
                    years,
                    roots,
                    values,
                ) in _iter_temporal_record_entries(
                    record,
                    selected_mask,
                    window_entries=window_entries,
                ):
                    yield (
                        int(record.base_year_index),
                        activities,
                        flows,
                        years,
                        roots,
                        values.astype(self.value_dtype, copy=False),
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
        self._temporal_records.clear()
        self._final_array = None
        self._activity_reduction_cache.clear()
        if self.owned_store:
            shutil.rmtree(self.store_path, ignore_errors=True)

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass
