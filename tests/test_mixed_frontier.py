"""Keep expanded biosphere supplies when routing merges frontier paths."""

import importlib

import networkx as nx
import numpy as np
import pytest

from trails.trails import Trails

lca_module = importlib.import_module("trails.lca")


@pytest.mark.parametrize("attribute_to_roots", [False, True])
@pytest.mark.parametrize("expanded", [2.0, -2.0, 0.0])
@pytest.mark.parametrize("frontier", [0.25, -0.25])
def test_mixed_node_preserves_disjoint_supplies(
    example_package, monkeypatch, attribute_to_roots, expanded, frontier
):
    """Inject only the expanded portion and solve only the frontier portion."""
    with Trails(example_package, interpolate_annual=False) as trails:
        indices = next(iter(trails.activity_indices.values()))
        fu, shared = [int(idx) for idx in list(indices)[:2]]
        graph = nx.DiGraph()
        graph.add_node(
            "fu",
            year=2005,
            depth=0,
            act_idx=fu,
            amount=1.0,
            frontier_amount=0.0,
            direct_bio_amount=0.0,
            frontier_roots={},
            direct_bio_roots={},
        )
        graph.add_node(
            "mixed",
            year=2005,
            depth=2,
            act_idx=shared,
            amount=expanded + frontier,
            frontier_amount=frontier,
            direct_bio_amount=expanded,
            frontier_roots={fu: frontier},
            direct_bio_roots={shared: expanded} if expanded else {},
        )
        trails.graph = graph
        trails._routing_attribute_to_roots = attribute_to_roots
        trails._routing_params = dict(
            start_year=2005, start_act_idx=fu, amount=1.0, min_amount=0.0
        )
        monkeypatch.setattr(
            trails, "_activity_amount_from_product_demand", lambda *args: 1.0
        )
        captured_frontier = {}
        captured_direct = {}

        def capture_frontier(demands):
            captured_frontier.update(demands)
            # Isolate the graph-to-supply handoff from matrix solving.
            return {}

        def capture_supply(*, base_year, supply_by_activity, **kwargs):
            for activity, amount in supply_by_activity.items():
                key = (base_year, activity, kwargs.get("root"))
                captured_direct[key] = captured_direct.get(key, 0.0) + amount

        def capture_batch(*, base_year, supplies, **kwargs):
            for supply, root in supplies:
                capture_supply(
                    base_year=base_year, supply_by_activity=supply, root=root
                )

        monkeypatch.setattr(trails, "frontier_to_demand_vectors", capture_frontier)
        monkeypatch.setattr(
            trails, "accumulate_temporalized_biosphere_inventory", capture_supply
        )
        monkeypatch.setattr(
            trails, "accumulate_temporalized_biosphere_inventory_batch", capture_batch
        )
        monkeypatch.setattr(trails, "finalize_inventory", lambda: None)
        lca_module.lca(
            trails=trails,
            show_progress=False,
            compute_score=False,
            store_inventory=True,
            attribute_to_roots=attribute_to_roots,
        )

        assert captured_frontier == {(2005, shared): frontier}
        expected = {(2005, fu, fu if attribute_to_roots else None): 1.0}
        if expanded:
            expected[(2005, shared, shared if attribute_to_roots else None)] = expanded
        assert captured_direct == expected


@pytest.mark.parametrize("attribute_to_roots", [False, True])
@pytest.mark.parametrize("solver_mode", ["direct", "iterative"])
def test_merging_frontier_and_expanded_nodes_preserves_inventory(
    example_package, attribute_to_roots, solver_mode
):
    """Real LCI must be invariant to merging disjoint paths into one node."""

    def calculate(merged):
        with Trails(example_package, interpolate_annual=False) as trails:
            indices = next(iter(trails.activity_indices.values()))
            shared = int(trails.B.coords[1, 0])
            fu = next(int(idx) for idx in indices if int(idx) != shared)
            graph = nx.DiGraph()
            graph.add_node(
                "fu",
                year=2005,
                depth=0,
                act_idx=fu,
                amount=1.0,
                frontier_amount=0.0,
                direct_bio_amount=0.0,
                frontier_roots={},
                direct_bio_roots={},
            )
            graph.add_node(
                "expanded",
                year=2005,
                depth=2,
                act_idx=shared,
                amount=2.25 if merged else 2.0,
                frontier_amount=0.25 if merged else 0.0,
                direct_bio_amount=2.0,
                frontier_roots={fu: 0.25} if merged else {},
                direct_bio_roots={shared: 2.0},
            )
            if not merged:
                graph.add_node(
                    "frontier",
                    year=2005,
                    depth=2,
                    act_idx=shared,
                    amount=0.25,
                    frontier_amount=0.25,
                    direct_bio_amount=0.0,
                    frontier_roots={fu: 0.25},
                    direct_bio_roots={},
                )
            trails.graph = graph
            trails._routing_attribute_to_roots = attribute_to_roots
            trails._routing_params = dict(
                start_year=2005, start_act_idx=fu, amount=1.0, min_amount=0.0
            )
            trails.lci(
                show_progress=False,
                attribute_to_roots=attribute_to_roots,
                solver_mode=solver_mode,
                inventory_backend="chunked",
            )
            inventory = trails.inventory
            assert inventory is not None
            values = inventory.data
            # Materialize this tiny fixture before closing its disk-backed store.
            if hasattr(values, "compute"):
                values = values.compute()
            return values.todense() if hasattr(values, "todense") else values.copy()

    split = calculate(False)
    merged = calculate(True)
    assert np.any(split)
    np.testing.assert_allclose(merged, split, rtol=1e-10, atol=1e-12)
