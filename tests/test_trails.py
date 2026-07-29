import numpy as np
import pytest
import sparse

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


def test_matrix_inventory_matches_temporal_dictionary_batches(
    example_package,
) -> None:
    """Matrix-native root accumulation must preserve temporal inventory values."""
    legacy = Trails(example_package, interpolate_annual=False)
    matrix_native = Trails(example_package, interpolate_annual=False)
    try:
        n_activities = int(legacy.A.shape[1])
        root_ids = np.array([0, 1], dtype=np.int64)
        supply_matrix = np.zeros((n_activities, root_ids.size), dtype=np.float64)
        supply_matrix[0, 0] = 2.0
        supply_matrix[2, 0] = 0.5
        supply_matrix[1, 1] = 1.25
        supply_matrix[2, 1] = 0.75
        supplies = [
            ({0: 2.0, 2: 0.5}, 0),
            ({1: 1.25, 2: 0.75}, 1),
        ]

        legacy.reset_inventory(attribute_to_roots=True)
        legacy.accumulate_temporalized_biosphere_inventory_batch(
            base_year=2005,
            supplies=supplies,
            min_amount=1e6,
            use_temporal_distributions=True,
        )
        expected = legacy.finalize_inventory().data

        matrix_native.reset_inventory(attribute_to_roots=True)
        matrix_native.accumulate_temporalized_biosphere_inventory_matrix(
            base_year=2005,
            supply_matrix=supply_matrix,
            root_activities=root_ids,
            min_amount=1e6,
            use_temporal_distributions=True,
        )
        actual = matrix_native.finalize_inventory().data

        assert np.allclose(actual.todense(), expected.todense())
    finally:
        legacy.close()
        matrix_native.close()


def test_reset_inventory_uses_inventory_offset_bounds(example_trails: Trails) -> None:
    """Verify inventory year axis uses both technosphere and biosphere offsets."""
    min_inv, max_inv = example_trails._inventory_offset_bounds()
    min_bio, _ = example_trails._biosphere_offset_bounds()

    example_trails.reset_inventory()

    years = example_trails._inventory_years
    assert years is not None and years.size > 0
    assert int(years[0]) == int(example_trails.min_year) + int(min_inv)
    assert int(years[-1]) == int(example_trails.max_year) + int(max_inv) + 500
    assert int(years[0]) <= int(example_trails.min_year) + int(min_bio)


def test_inventory_root_reduce_uses_int64_coords(example_trails: Trails) -> None:
    """Inventory reductions should work with large shapes and root attribution."""
    example_trails.A = sparse.COO(
        np.array([[0], [0], [0]], dtype=np.int64),
        np.array([1.0], dtype=np.float32),
        shape=(1, 42000, 1),
    )
    example_trails.B = sparse.COO(
        np.array([[0], [0], [0]], dtype=np.int64),
        np.array([1.0], dtype=np.float32),
        shape=(1, 42000, 50000),
    )

    example_trails.reset_inventory(attribute_to_roots=True)
    assert example_trails._inventory_years is not None
    y0 = int(example_trails._inventory_years[0])

    example_trails._append_inventory_entries_bulk(
        np.array([0, 1, 2], dtype=np.int64),
        np.array([y0, y0, y0], dtype=np.int64),
        np.array([0, 1, 2], dtype=np.int64),
        np.array([1.0, 2.0, 3.0], dtype=np.float64),
        root_activity=np.array([0, 1, 2], dtype=np.int64),
    )

    inv = example_trails.finalize_inventory()
    assert isinstance(inv.data, sparse.COO)
    assert inv.data.coords.dtype == np.int64

    reduced = inv.sum(dim=["activity"]).transpose("flow", "year", "root activity")
    assert isinstance(reduced.data, sparse.COO)
    assert reduced.data.coords.dtype == np.int64


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
