import numpy as np
import pytest

from trails.temporal_distributions import TemporalDistribution


def test_map_year_helpers(example_trails):
    assert example_trails._map_year_to_scenario_year(2010) == 2005
    assert example_trails._map_year_to_template_year(2010) == 2005
    assert example_trails._map_year_to_available(2010) == 2005


def test_get_matrices_for_scenario(example_trails):
    A_2005 = example_trails.get_A_for_scenario("2005")
    B_2005 = example_trails.get_B_for_scenario("2005")
    assert A_2005.shape[0] == example_trails.A.shape[1]
    assert B_2005.shape[0] == example_trails.B.shape[1]


def test_get_temporal_exchange(example_trails):
    tex = example_trails.get_temporal_exchange(2005, 2, 0)
    assert tex is not None


def test_get_temporal_distribution(example_trails):
    td = example_trails.get_temporal_distribution(2005, 2, 0)
    assert isinstance(td, TemporalDistribution)


def test_expand_temporal_exchanges_without_temporal_distributions(example_trails):
    demand = example_trails.expand_temporal_exchanges(
        year=2005, act_idx=1, amount=2.0, use_temporal_distributions=False
    )

    assert set(demand.keys()) == {2005}
    assert set(demand[2005].keys()) == {5, 7}

    assert demand[2005][5] == pytest.approx(1.4)
    assert demand[2005][7] == pytest.approx(0.6)



def test_expand_temporal_exchanges_with_temporal_distributions(example_trails):
    demand = example_trails.expand_temporal_exchanges(year=2005, act_idx=2, amount=1.0)
    assert any(0 in mapping for mapping in demand.values())


def test_accumulate_temporalized_biosphere_inventory(example_trails):
    inventory_by_year = {}
    example_trails.accumulate_temporalized_biosphere_inventory(
        base_year=2005,
        supply_by_activity={0: 2.0},
        inventory_by_year=inventory_by_year,
        use_temporal_distributions=False,
    )
    assert np.isclose(inventory_by_year[2005][0], 30000.0)


def test_temporal_traversal_basic(example_trails):
    out = example_trails.temporal_traversal(
        start_year=2005,
        start_act_idx=1,
        amount=1.0,
        max_depth=1,
        min_amount=0.0,
        return_provenance=True,
        show_progress=False,
        use_temporal_distributions=False,
    )

    # Backwards-compatible: allow (frontier, provenance, ...) tuples
    demand, provenance = out[0], out[1]

    assert (2005, 5) in demand
    assert (2005, 7) in demand
    assert demand[(2005, 5)] == pytest.approx(0.7)
    assert demand[(2005, 7)] == pytest.approx(0.3)



def test_frontier_to_demand_vectors(example_trails):
    frontier = {(2005, 1): 1.0, (2006, 2): 2.0}
    f_by_year = example_trails.frontier_to_demand_vectors(frontier)
    assert f_by_year[2005][1] == 1.0
    assert f_by_year[2006][2] == 2.0


def test_collect_traversal_edges(example_trails):
    edges = example_trails.collect_traversal_edges(
        start_year=2005, start_act_idx=1, amount=1.0, max_depth=1, min_amount=0.0
    )
    assert 0 in edges
    assert any(edge[0] == (2005, 1) for edge in edges[0])
