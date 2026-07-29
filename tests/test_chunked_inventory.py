from __future__ import annotations

import importlib
from pathlib import Path

import dask.array as da
from datapackage import Package
import numpy as np
import pytest
from scipy import sparse as sp
import sparse
import xarray as xr

from trails import Trails
import trails.chunked_inventory as chunked_inventory
from trails.chunked_inventory import (
    ChunkedInventoryBuilder,
    estimate_decode_peak_bytes,
    estimate_flush_peak_bytes,
    estimate_materialization_peak_bytes,
    estimate_merge_peak_bytes,
)
from trails.factorized_inventory import FactorizedInventoryBuilder
from trails.lcia import get_lcia_method_names
from trails.edges_matrix import score_inventory_with_edges
from trails.fair_rf import _inventory_flow_year_root
from trails.plotting import _characterized_inventory_to_root_results

lca_module = importlib.import_module("trails.lca")


def test_memory_estimators_scale_with_entries_and_dimensions() -> None:
    small = estimate_flush_peak_bytes(1_000, value_dtype=np.float32)
    large = estimate_flush_peak_bytes(2_000, value_dtype=np.float32)
    assert large == 2 * small
    assert estimate_flush_peak_bytes(1_000, value_dtype=np.float64) > small
    assert estimate_decode_peak_bytes(
        1_000, ndim=4, value_dtype=np.float32
    ) > estimate_decode_peak_bytes(1_000, ndim=3, value_dtype=np.float32)
    assert estimate_materialization_peak_bytes(
        1_000, ndim=4, value_dtype=np.float32
    ) > estimate_materialization_peak_bytes(1_000, ndim=3, value_dtype=np.float32)
    assert estimate_merge_peak_bytes(1_000, 1_000, value_dtype=np.float32) > 0


@pytest.mark.parametrize("has_root", [False, True])
def test_chunked_builder_merges_duplicates_and_cancellations(
    has_root: bool,
) -> None:
    builder = ChunkedInventoryBuilder(
        n_activities=5,
        n_flows=3,
        n_years=2,
        has_root=has_root,
        value_dtype=np.float32,
        memory_budget=2**20,
    )
    kwargs = {"roots": np.array([2, 2, 4, 2])} if has_root else {}
    builder.append(
        np.array([0, 0, 4, 0]),
        np.array([1, 1, 2, 1]),
        np.array([0, 0, 1, 0]),
        np.array([1.0, 2.0, 3.0, -3.0], dtype=np.float32),
        **kwargs,
    )
    result = builder.finalize().compute(scheduler="synchronous")
    assert isinstance(result, sparse.COO)
    assert result.nnz == 1
    assert float(result.sum()) == pytest.approx(3.0)
    diagnostics = builder.diagnostics()
    assert diagnostics["raw_entries"] == 4
    assert diagnostics["canonical_entries"] == 1
    store = builder.store_path
    assert store.exists()
    builder.close()
    assert not store.exists()


def test_factorized_builder_keeps_base_inventory_lazy_and_adds_corrections() -> None:
    builder = FactorizedInventoryBuilder(
        n_activities=3,
        n_flows=2,
        n_years=2,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
        activity_block_size=2,
        root_block_size=2,
        year_block_size=1,
    )
    builder.append_factor(
        year_index=0,
        activities=np.array([0, 1, 1]),
        flows=np.array([0, 0, 1]),
        biosphere_values=np.array([2.0, 3.0, 4.0]),
        supply_matrix=np.array(
            [
                [1.0, 2.0],
                [0.5, 0.0],
                [0.0, 1.0],
            ]
        ),
        roots=np.array([0, 2]),
    )
    builder.append(
        np.array([1]),
        np.array([1]),
        np.array([1]),
        np.array([7.0], dtype=np.float32),
        roots=np.array([2]),
    )
    result = builder.finalize()
    store = builder.store_path
    try:
        assert isinstance(result, da.Array)
        actual = result.compute(scheduler="synchronous")
        expected = np.zeros((3, 2, 2, 3), dtype=np.float32)
        expected[0, 0, 0, 0] = 2.0
        expected[0, 0, 0, 2] = 4.0
        expected[1, 0, 0, 0] = 1.5
        expected[1, 1, 0, 0] = 2.0
        expected[1, 1, 1, 2] = 7.0
        assert np.allclose(actual.todense(), expected)
        reduced = builder.reduce_activity_for_flows([1])
        expected_reduced = expected.sum(axis=0)
        expected_reduced[0, :, :] = 0.0
        assert np.allclose(reduced.todense(), expected_reduced)
        diagnostics = builder.diagnostics()
        assert diagnostics["backend"] == "factorized"
        assert diagnostics["factor_count"] == 1
        assert diagnostics["explicit_correction_entries"] == 1
    finally:
        builder.close()
    assert not store.exists()


def test_factorized_builder_keeps_ported_temporal_kernels_lazy() -> None:
    builder = FactorizedInventoryBuilder(
        n_activities=3,
        n_flows=2,
        n_years=4,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
        activity_block_size=2,
        root_block_size=2,
        year_block_size=2,
    )
    builder.append_temporal_factor(
        base_year_index=1,
        activities=np.array([0, 1]),
        flows=np.array([0, 1]),
        biosphere_values=np.array([2.0, 4.0]),
        supply_matrix=np.array(
            [
                [1.0, 2.0],
                [0.1, 0.0],
                [0.0, 0.0],
            ]
        ),
        roots=np.array([0, 2]),
        pulse_indptr=np.array([0, 2, 4]),
        pulse_year_indices=np.array([1, 3, 0, 2]),
        pulse_weights=np.array([0.25, 0.75, 0.5, 0.5]),
        min_amount=1.0,
    )
    result = builder.finalize()
    try:
        actual = result.compute(scheduler="synchronous")
        expected = np.zeros((3, 2, 4, 3), dtype=np.float32)
        expected[0, 0, 1, 0] = 0.5
        expected[0, 0, 3, 0] = 1.5
        expected[0, 0, 1, 2] = 1.0
        expected[0, 0, 3, 2] = 3.0
        # Below-threshold values remain anchored at the base year.
        expected[1, 1, 1, 0] = 0.4
        assert np.allclose(actual.todense(), expected)

        reduced = builder.reduce_activity_for_flows([0, 1])
        assert np.allclose(reduced.todense(), expected.sum(axis=0))

        streamed = list(builder.iter_entries_for_flows([0, 1]))
        coords = np.concatenate(
            [np.vstack([part[1], part[2], part[3], part[4]]) for part in streamed],
            axis=1,
        )
        values = np.concatenate([part[5] for part in streamed])
        streamed_inventory = sparse.COO(
            coords,
            values,
            shape=expected.shape,
            has_duplicates=True,
        )
        assert np.allclose(streamed_inventory.todense(), expected)

        diagnostics = builder.diagnostics()
        assert diagnostics["temporal_factor_count"] == 1
        assert diagnostics["temporal_factor_candidate_entries"] == 3
        assert diagnostics["temporal_factor_pulse_entries"] == 6
        assert diagnostics["explicit_correction_entries"] == 0
    finally:
        builder.close()


def test_explicit_builder_store_is_not_deleted(tmp_path: Path) -> None:
    builder = ChunkedInventoryBuilder(
        n_activities=1,
        n_flows=1,
        n_years=1,
        has_root=False,
        value_dtype=np.float32,
        memory_budget=2**20,
        store=tmp_path,
    )
    builder.append(
        np.array([0]),
        np.array([0]),
        np.array([0]),
        np.array([1.0], dtype=np.float32),
    )
    builder.finalize()
    store = builder.store_path
    builder.close()
    assert store.exists()


def test_many_runs_compact_to_one_bounded_result() -> None:
    builder = ChunkedInventoryBuilder(
        n_activities=2,
        n_flows=2,
        n_years=1,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
    )
    batch_size = 10_000
    batches = 25
    for _ in range(batches):
        builder.append(
            np.zeros(batch_size, dtype=np.int64),
            np.ones(batch_size, dtype=np.int64),
            np.zeros(batch_size, dtype=np.int64),
            np.ones(batch_size, dtype=np.float32),
            roots=np.ones(batch_size, dtype=np.int64),
        )
    result = builder.finalize().compute(scheduler="synchronous")
    try:
        assert result.nnz == 1
        assert float(result[0, 1, 0, 1]) == pytest.approx(batch_size * batches)
        diagnostics = builder.diagnostics()
        assert diagnostics["raw_entries"] == batch_size * batches
        assert diagnostics["canonical_entries"] == 1
        assert diagnostics["peak_buffer_bytes"] <= 2**20
        assert diagnostics["bytes_merged"] > 0
    finally:
        builder.close()


def test_spills_to_shards_and_coalesces_year_chunks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        chunked_inventory, "DEFAULT_BUCKET_COMPACTION_BYTES", 256 * 2**10
    )
    builder = ChunkedInventoryBuilder(
        n_activities=16,
        n_flows=4,
        n_years=17,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
    )
    count = 200_000
    rng = np.random.default_rng(42)
    builder.append(
        rng.integers(0, 16, count),
        rng.integers(0, 4, count),
        rng.integers(0, 17, count),
        np.ones(count, dtype=np.float32),
        roots=rng.integers(0, 16, count),
    )
    result = builder.finalize()
    try:
        assert result.chunks[2] == (8, 8, 1)
        assert float(result.sum().compute(scheduler="synchronous")) == pytest.approx(
            count
        )
        diagnostics = builder.diagnostics()
        assert diagnostics["online_merge_seconds"] > 0.0
        assert diagnostics["bytes_merged"] > 0
        assert diagnostics["storage_file_count"] == 3
        assert diagnostics["dask_block_count"] == 3
        assert diagnostics["logical_partition_count"] == 3
    finally:
        builder.close()


def test_geometric_bucket_compaction_bounds_write_amplification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        chunked_inventory,
        "DEFAULT_BUCKET_COMPACTION_BYTES",
        8 * 2**10,
    )
    builder = ChunkedInventoryBuilder(
        n_activities=1_024,
        n_flows=256,
        n_years=1,
        has_root=False,
        value_dtype=np.float32,
        memory_budget=2**21,
    )
    batch_size = 7_500
    batches = 20
    try:
        for batch in range(batches):
            linear = np.arange(
                batch * batch_size,
                (batch + 1) * batch_size,
                dtype=np.int64,
            )
            builder.append(
                linear // 256,
                linear % 256,
                np.zeros(batch_size, dtype=np.int64),
                np.ones(batch_size, dtype=np.float32),
            )

        builder._finish_runs()
        diagnostics = builder.diagnostics()
        assert diagnostics["canonical_entries"] == batch_size * batches
        final_payload = int(diagnostics["canonical_entries"]) * (
            np.dtype(np.int64).itemsize + np.dtype(np.float32).itemsize
        )
        assert diagnostics["online_compaction_count"] <= 5
        assert diagnostics["bytes_written"] <= 4 * final_payload
    finally:
        builder.close()


def test_oversized_single_run_is_copied_to_stable_final_shard() -> None:
    builder = ChunkedInventoryBuilder(
        n_activities=4_096,
        n_flows=10,
        n_years=1,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**21,
    )
    count = 30_000
    linear = np.arange(count, dtype=np.int64)
    builder.append(
        np.zeros(count, dtype=np.int64),
        linear // 4_096,
        np.zeros(count, dtype=np.int64),
        np.ones(count, dtype=np.float32),
        roots=linear % 4_096,
    )
    result = builder.finalize()
    try:
        assert float(result.sum().compute(scheduler="synchronous")) == pytest.approx(
            count
        )
        assert builder.diagnostics()["storage_file_count"] == 3
    finally:
        builder.close()


def test_selective_activity_reduction_streams_and_caches_selected_flows() -> None:
    builder = ChunkedInventoryBuilder(
        n_activities=10,
        n_flows=5,
        n_years=9,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
        activity_block_size=4,
        root_block_size=4,
        year_block_size=4,
    )
    rng = np.random.default_rng(7)
    count = 25_000
    builder.append(
        rng.integers(0, 10, count),
        rng.integers(0, 5, count),
        rng.integers(0, 9, count),
        rng.normal(size=count).astype(np.float32),
        roots=rng.integers(0, 10, count),
    )
    inventory = builder.finalize()
    try:
        reduced = builder.reduce_activity_for_flows([1, 3])
        expected = inventory.compute(scheduler="synchronous").sum(axis=0)
        expected_dense = expected.todense()
        expected_dense[[0, 2, 4], :, :] = 0.0

        assert np.allclose(reduced.todense(), expected_dense, atol=1e-5)
        assert builder.reduce_activity_for_flows([3, 1]) is reduced
    finally:
        builder.close()


def _run_example_lca(
    package_path: Path,
    *,
    backend: str,
    memory_budget: int = 2**20,
) -> Trails:
    trails = Trails(Package(str(package_path)), interpolate_annual=False)
    activity_indices = next(iter(trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))
    method = get_lcia_method_names(ei_version="3.11")[:1]
    trails.temporal_routing(
        start_year=2005,
        start_act_idx=start_act_idx,
        max_depth=1,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=True,
    )
    trails.lca(
        methods=method,
        show_progress=False,
        compute_score=True,
        store_inventory=True,
        attribute_to_roots=True,
        solver_mode="direct",
        inventory_backend=backend,
        inventory_memory_budget=memory_budget,
    )
    return trails


def test_chunked_lca_preserves_inventory_and_characterization() -> None:
    package_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "example data package"
        / "datapackage.json"
    )
    eager = _run_example_lca(package_path, backend="coo")
    chunked = _run_example_lca(package_path, backend="chunked")
    try:
        assert isinstance(eager.inventory.data, sparse.COO)
        assert isinstance(chunked.inventory.data, da.Array)
        assert isinstance(chunked.characterized_inventory.data, da.Array)
        chunked_inventory = chunked.inventory.data.compute(scheduler="synchronous")
        chunked_characterized = chunked.characterized_inventory.data.compute(
            scheduler="synchronous"
        )
        assert np.allclose(chunked_inventory.todense(), eager.inventory.data.todense())
        assert np.allclose(
            chunked_characterized.todense(),
            eager.characterized_inventory.data.todense(),
        )
        assert float(chunked.scores.data.sum()) == pytest.approx(
            float(eager.scores.data.sum())
        )
        score_data = chunked.scores.data
        reduced_data = chunked.characterized_inventory.sum(dim="flow").data.compute(
            scheduler="synchronous"
        )
        assert np.allclose(
            score_data.todense(),
            reduced_data.todense(),
        )
        with pytest.raises(MemoryError):
            chunked.materialize_inventory(memory_budget=1)
    finally:
        eager.close()
        chunked.close()


def test_factorized_lca_preserves_inventory_and_characterization() -> None:
    package_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "example data package"
        / "datapackage.json"
    )
    eager = _run_example_lca(package_path, backend="coo")
    factorized = _run_example_lca(package_path, backend="factorized")
    try:
        assert isinstance(factorized.inventory.data, da.Array)
        assert isinstance(factorized.characterized_inventory.data, da.Array)
        actual = factorized.inventory.data.compute(scheduler="synchronous")
        characterized = factorized.characterized_inventory.data.compute(
            scheduler="synchronous"
        )
        assert np.allclose(actual.todense(), eager.inventory.data.todense())
        assert np.allclose(
            characterized.todense(),
            eager.characterized_inventory.data.todense(),
        )
        assert float(factorized.scores.data.sum()) == pytest.approx(
            float(eager.scores.data.sum())
        )
        assert factorized.inventory_diagnostics["backend"] == "factorized"
    finally:
        eager.close()
        factorized.close()


def test_auto_lca_selects_factorized_when_estimate_exceeds_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    package_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "example data package"
        / "datapackage.json"
    )
    monkeypatch.setattr(
        lca_module,
        "estimate_materialization_peak_bytes",
        lambda *args, **kwargs: 2**21,
    )
    automatic = _run_example_lca(package_path, backend="auto")
    try:
        diagnostics = automatic.inventory_diagnostics
        assert diagnostics["backend"] == "factorized"
        selection = diagnostics["backend_selection"]
        assert selection["requested"] == "auto"
        assert selection["selected"] == "factorized"
        assert selection["estimated_peak_bytes"] > selection["memory_budget"]
    finally:
        automatic.close()


def test_coo_backend_raises_before_unsafe_finalize(example_trails: Trails) -> None:
    example_trails.configure_inventory_storage(backend="coo", memory_budget=1)
    example_trails.reset_inventory()
    example_trails.accumulate_temporalized_biosphere_inventory(
        base_year=2005,
        supply_by_activity={0: 2.0},
        use_temporal_distributions=False,
    )
    with pytest.raises(MemoryError, match="Eager inventory finalization"):
        example_trails.finalize_inventory()


def test_auto_backend_promotes_before_unsafe_finalize(example_trails: Trails) -> None:
    example_trails.configure_inventory_storage(backend="auto", memory_budget=2**20)
    example_trails.reset_inventory(attribute_to_roots=True)
    assert example_trails._inventory_years is not None
    count = 20_000
    example_trails._append_inventory_entries_bulk(
        np.zeros(count, dtype=np.int64),
        np.full(count, int(example_trails._inventory_years[0]), dtype=np.int64),
        np.zeros(count, dtype=np.int64),
        np.ones(count, dtype=np.float32),
        root_activity=np.zeros(count, dtype=np.int64),
    )
    inventory = example_trails.finalize_inventory()
    assert isinstance(inventory.data, da.Array)
    assert float(
        inventory.data.sum().compute(scheduler="synchronous")
    ) == pytest.approx(float(count))


def _small_lazy_inventory() -> tuple[ChunkedInventoryBuilder, xr.DataArray]:
    builder = ChunkedInventoryBuilder(
        n_activities=2,
        n_flows=2,
        n_years=2,
        has_root=True,
        value_dtype=np.float32,
        memory_budget=2**20,
    )
    builder.append(
        np.array([0, 1, 1]),
        np.array([0, 0, 1]),
        np.array([0, 0, 1]),
        np.array([2.0, 3.0, 4.0], dtype=np.float32),
        roots=np.array([0, 1, 1]),
    )
    inventory = xr.DataArray(
        builder.finalize(),
        dims=("activity", "flow", "year", "root activity"),
        coords={
            "activity": [0, 1],
            "flow": [0, 1],
            "year": [2000, 2001],
            "root activity": [0, 1],
        },
    )
    return builder, inventory


def test_selected_flow_stream_can_reuse_one_external_progress_bar(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    builder, _inventory = _small_lazy_inventory()

    class ExternalProgress:
        def __init__(self) -> None:
            self.n = 0
            self.closed = False

        def update(self, amount: int) -> None:
            self.n += int(amount)

        def close(self) -> None:
            self.closed = True

    progress = ExternalProgress()
    monkeypatch.setattr(
        "trails.chunked_inventory.tqdm",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("an external progress bar should prevent nested tqdm")
        ),
    )
    try:
        rows = list(
            builder.iter_entries_for_flows(
                [0, 1],
                show_progress=True,
                progress=progress,
            )
        )
        assert rows
        assert progress.n == builder.nnz
        assert progress.closed is False
    finally:
        builder.close()


def test_fair_and_plotting_reduce_lazy_blocks_sequentially() -> None:
    builder, inventory = _small_lazy_inventory()
    try:
        reduced = _inventory_flow_year_root(inventory)
        assert reduced.shape == (2, 2, 2)
        assert float(reduced.sum()) == pytest.approx(9.0)

        characterized = inventory.expand_dims(method=["test"])
        results = _characterized_inventory_to_root_results(characterized)
        assert results[2000]["scores"] == pytest.approx(5.0)
        assert results[2001]["scores"] == pytest.approx(4.0)
    finally:
        builder.close()


def test_edges_scores_lazy_inventory_without_materializing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    builder, inventory = _small_lazy_inventory()

    class DummyTrails:
        value_dtype = np.float32

        def __init__(self, value: xr.DataArray) -> None:
            self.inventory = value
            self.scores = None
            self._inventory_builder = builder

    matrix_calls: dict[int, int] = {}

    def fake_matrices(*args, year: int, progress=None, **kwargs):
        matrix_calls[int(year)] = matrix_calls.get(int(year), 0) + 1
        if progress is not None:
            progress.update(1)
        multiplier = 1.0 if int(year) == 2000 else 2.0
        return [sp.csr_matrix(np.full((2, 2), multiplier))]

    monkeypatch.setattr(
        "trails.edges_matrix._build_edges_characterization_matrices_for_year",
        fake_matrices,
    )
    monkeypatch.setattr(
        "trails.edges_matrix._eligible_biosphere_flow_ids",
        lambda *args, **kwargs: np.array([0, 1], dtype=np.int64),
    )
    monkeypatch.setattr(
        "trails.edges_matrix.iter_sparse_blocks",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("direct finalized-run path should avoid Dask blocks")
        ),
    )
    trails = DummyTrails(inventory)
    try:
        scores = score_inventory_with_edges(
            trails,
            [{"name": "test"}],
            show_progress=False,
        )
        assert float(scores.data.sum()) == pytest.approx(13.0)
        assert float(scores.sel(year=2000).data.sum()) == pytest.approx(5.0)
        assert float(scores.sel(year=2001).data.sum()) == pytest.approx(8.0)
        assert matrix_calls == {2000: 1, 2001: 1}
    finally:
        builder.close()
