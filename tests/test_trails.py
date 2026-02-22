import numpy as np
import pytest

from trails.trails import Trails

from trails.temporal_distributions import TemporalDistribution, TemporalExchange


def test_map_year_helpers(example_trails: Trails) -> None:
    """Verify year mapping helpers for scenario and template years.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    assert example_trails._map_year_to_scenario_year(2010) == 2005
    assert example_trails._map_year_to_template_year(2010) == 2005
    assert example_trails._map_year_to_available(2010) == 2005


def test_get_matrices_for_scenario(example_trails: Trails) -> None:
    """Verify scenario matrix accessors.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    A_2005 = example_trails.get_A_for_scenario("2005")
    B_2005 = example_trails.get_B_for_scenario("2005")
    assert A_2005.shape[0] == example_trails.A.shape[1]
    assert B_2005.shape[0] == example_trails.B.shape[1]


def test_get_temporal_exchange(example_trails: Trails) -> None:
    """Verify temporal exchange lookup behavior.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    tex = example_trails.get_temporal_exchange(2005, 2, 0)
    assert tex is not None


def test_get_temporal_distribution(example_trails: Trails) -> None:
    """Verify temporal distribution construction.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    td = example_trails.get_temporal_distribution(2005, 2, 0)
    assert isinstance(td, TemporalDistribution)


def test_interpolate_temporal_exchange_keeps_explicit_pulses(
    example_trails: Trails,
) -> None:
    """Verify interpolation preserves explicit pulse definitions.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    exchanges = {
        ("2005", 1, 2): TemporalExchange(
            distribution=6,
            loc=None,
            scale=None,
            offset_min=0,
            offset_max=0,
            amount_source="port",
            offsets=[-1, 9],
            weights=[0.5, 0.5],
        ),
        ("2020", 1, 2): TemporalExchange(
            distribution=6,
            loc=None,
            scale=None,
            offset_min=0,
            offset_max=0,
            amount_source="port",
            offsets=[-1, 9],
            weights=[0.5, 0.5],
        ),
    }
    tex = example_trails._interpolate_temporal_exchange(2010, 1, 2, exchanges)
    assert tex is not None
    assert tex.offsets == [-1, 9]
    assert tex.weights == [0.5, 0.5]
    assert list(TemporalDistribution(tex).iter_offsets_and_weights()) == [
        (-1, pytest.approx(0.5)),
        (9, pytest.approx(0.5)),
    ]


def test_td_offsets_cache_distinguishes_explicit_pulse_sets(
    example_trails: Trails,
) -> None:
    """Verify TD cache key includes explicit pulse vectors.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    tex_a = TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        amount_source="port",
        offsets=[-1, 9],
        weights=[0.5, 0.5],
    )
    tex_b = TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        amount_source="port",
        offsets=[0, 1],
        weights=[0.9, 0.1],
    )
    out_a = example_trails._get_td_offsets(tex=tex_a, debug=False)
    out_b = example_trails._get_td_offsets(tex=tex_b, debug=False)
    assert out_a == [(-1, pytest.approx(0.5)), (9, pytest.approx(0.5))]
    assert out_b == [(0, pytest.approx(0.9)), (1, pytest.approx(0.1))]


def test_expand_temporal_exchanges_without_temporal_distributions(
    example_trails: Trails,
) -> None:
    """Verify exchange expansion without temporal distributions.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    demand = example_trails.expand_temporal_exchanges(
        year=2005, act_idx=1, amount=2.0, use_temporal_distributions=False
    )

    assert set(demand.keys()) == {2005}
    assert set(demand[2005].keys()) == {5, 7}

    assert demand[2005][5] == pytest.approx(1.4)
    assert demand[2005][7] == pytest.approx(0.6)


def test_expand_temporal_exchanges_with_temporal_distributions(
    example_trails: Trails,
) -> None:
    """Verify exchange expansion with temporal distributions.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    demand = example_trails.expand_temporal_exchanges(year=2005, act_idx=2, amount=1.0)
    assert any(0 in mapping for mapping in demand.values())


def test_accumulate_temporalized_biosphere_inventory(example_trails: Trails) -> None:
    """Verify temporalized biosphere inventory accumulation.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    example_trails.reset_inventory()
    example_trails.accumulate_temporalized_biosphere_inventory(
        base_year=2005,
        supply_by_activity={0: 2.0},
        use_temporal_distributions=False,
    )
    inv = example_trails.finalize_inventory()
    year_idx = int(np.where(inv.coords["year"].values == 2005)[0][0])
    value = inv.data[0, 0, year_idx]
    assert np.isclose(float(value), 30000.0)


def test_temporal_traversal_basic(example_trails: Trails) -> None:
    """Verify basic temporal traversal behavior.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
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


def test_frontier_to_demand_vectors(example_trails: Trails) -> None:
    """Verify conversion from frontier to demand vectors.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    frontier = {(2005, 1): 1.0, (2006, 2): 2.0}
    f_by_year = example_trails.frontier_to_demand_vectors(frontier)
    assert f_by_year[2005][1] == 1.0
    assert f_by_year[2006][2] == 2.0


def test_collect_traversal_edges(example_trails: Trails) -> None:
    """Verify traversal edge collection.

    :param example_trails: Trails fixture under test.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    edges = example_trails.collect_traversal_edges(
        start_year=2005, start_act_idx=1, amount=1.0, max_depth=1, min_amount=0.0
    )
    assert 0 in edges
    assert any(edge[0] == (2005, 1) for edge in edges[0])
