from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import sparse
import xarray as xr

from trails.lcia import lcia, resolve_lcia_methods
from trails.lcia import get_lcia_method_names
from trails.trails import Trails


def test_resolve_lcia_methods_requires_configuration() -> None:
    trails = SimpleNamespace(
        default_methods=None,
        default_edges_methods=None,
        default_method_backend="auto",
    )
    with pytest.raises(ValueError, match="No LCIA methods configured"):
        resolve_lcia_methods(trails, None, None)


def test_resolve_lcia_methods_infers_regular_and_edges() -> None:
    trails = SimpleNamespace(
        default_methods=None,
        default_edges_methods=None,
        default_method_backend="auto",
    )
    methods, backend = resolve_lcia_methods(trails, ["regular"], None)
    assert methods == ["regular"]
    assert backend == "regular"

    edge_method = {"name": "edge method", "exchanges": []}
    methods, backend = resolve_lcia_methods(trails, [edge_method], None)
    assert methods == [edge_method]
    assert backend == "edges"


def test_lcia_requires_lci_first() -> None:
    trails = SimpleNamespace(inventory=None)
    with pytest.raises(RuntimeError, match=r"run trails\.lci\(\)"):
        lcia(trails, methods=["regular"])


def test_lci_then_regular_lcia_reuses_inventory(example_trails: Trails) -> None:
    activity_indices = next(iter(example_trails.activity_indices.values()))
    start_act_idx = int(next(iter(activity_indices)))
    method = get_lcia_method_names(ei_version="3.11")[0]

    example_trails.temporal_routing(
        start_year=2005,
        start_act_idx=start_act_idx,
        max_depth=1,
        adaptive_relative_score_cutoff=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=True,
    )
    inventory = example_trails.lci(show_progress=False, solver_mode="direct")

    assert inventory is example_trails.inventory
    assert example_trails.scores is None
    assert example_trails.characterized_inventory is None

    inventory_id = id(example_trails.inventory)
    scores = example_trails.lcia(methods=[method], show_progress=False)

    assert id(example_trails.inventory) == inventory_id
    assert scores is example_trails.scores
    assert example_trails.characterized_inventory is not None
    assert example_trails.lcia_diagnostics["backend"] == "regular"
    assert example_trails.current_lcia_result in example_trails.lcia_results
    assert float(scores.data.sum()) == pytest.approx(
        float(example_trails.characterized_inventory.data.sum()), rel=1e-7
    )


def test_lcia_dispatches_edges_and_records_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array([[0], [0], [0], [0]]),
            data=np.array([2.0]),
            shape=(1, 1, 1, 1),
        ),
        dims=("activity", "flow", "year", "root activity"),
        coords={
            "activity": [0],
            "flow": [0],
            "year": [2025],
            "root activity": [0],
        },
    )
    expected = xr.DataArray(
        sparse.COO(
            coords=np.array([[0], [0], [0]]),
            data=np.array([3.0]),
            shape=(1, 1, 1),
        ),
        dims=("activity", "year", "root activity"),
    )
    calls: list[dict[str, object]] = []

    def fake_edges(_trails, methods, **kwargs):
        calls.append({"methods": methods, **kwargs})
        return expected

    monkeypatch.setattr("trails.edges_matrix.score_inventory_with_edges", fake_edges)
    trails = SimpleNamespace(
        inventory=inventory,
        characterized_inventory=object(),
        scores=None,
        default_methods=None,
        default_edges_methods=None,
        default_method_backend="auto",
        default_ei_version="3.11",
        lcia_results={},
        current_lcia_result=None,
        lcia_diagnostics={},
    )
    method = {"name": "AWARE", "exchanges": []}

    result = lcia(
        trails,
        methods=[method],
        reuse_mappings=False,
        show_progress=False,
    )

    assert result is expected
    assert trails.characterized_inventory is None
    assert calls[0]["reuse_cached_cfs"] is False
    assert trails.lcia_diagnostics["backend"] == "edges"
    assert trails.current_lcia_result == "edges:AWARE"
