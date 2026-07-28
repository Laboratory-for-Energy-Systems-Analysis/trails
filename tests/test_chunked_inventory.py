from __future__ import annotations

from pathlib import Path

import dask.array as da
from datapackage import Package
import numpy as np
import pytest
from scipy import sparse as sp
import sparse
import xarray as xr

from trails import Trails
from trails.chunked_inventory import (
    ChunkedInventoryBuilder,
    estimate_decode_peak_bytes,
    estimate_flush_peak_bytes,
    estimate_materialization_peak_bytes,
    estimate_merge_peak_bytes,
)
from trails.lcia import get_lcia_method_names
from trails.edges_matrix import score_inventory_with_edges
from trails.fair_rf import _inventory_flow_year_root
from trails.plotting import _characterized_inventory_to_root_results


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
    ) > estimate_materialization_peak_bytes(
        1_000, ndim=3, value_dtype=np.float32
    )
    assert estimate_merge_peak_bytes(
        1_000, 1_000, value_dtype=np.float32
    ) > 0


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


def _run_example_lca(package_path: Path, *, backend: str) -> Trails:
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
        inventory_memory_budget=2**20,
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
        assert np.allclose(
            chunked_inventory.todense(), eager.inventory.data.todense()
        )
        assert np.allclose(
            chunked_characterized.todense(),
            eager.characterized_inventory.data.todense(),
        )
        assert float(chunked.scores.data.sum()) == pytest.approx(
            float(eager.scores.data.sum())
        )
        with pytest.raises(MemoryError):
            chunked.materialize_inventory(memory_budget=1)
    finally:
        eager.close()
        chunked.close()


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
    example_trails.configure_inventory_storage(
        backend="auto", memory_budget=2**20
    )
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
    assert float(inventory.data.sum().compute(scheduler="synchronous")) == pytest.approx(
        float(count)
    )


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

    def fake_matrices(*args, year: int, progress=None, **kwargs):
        if progress is not None:
            progress.update(1)
        multiplier = 1.0 if int(year) == 2000 else 2.0
        return [sp.csr_matrix(np.full((2, 2), multiplier))]

    monkeypatch.setattr(
        "trails.edges_matrix._build_edges_characterization_matrices_for_year",
        fake_matrices,
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
    finally:
        builder.close()
