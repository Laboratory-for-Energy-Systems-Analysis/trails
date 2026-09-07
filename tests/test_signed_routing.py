"""Compare signed temporal routing with the signed technosphere equations."""

import numpy as np
import pytest
import sparse

from trails.temporal_distributions import TemporalExchange
from trails.trails import Trails


@pytest.fixture
def signed_model(example_package):
    with Trails(example_package, interpolate_annual=False) as model:
        n = model.A.shape[1]
        assert n >= 3
        a = np.eye(n)
        # Non-unit production, a normal input, an avoided supplier, and a cycle.
        a[:3, :3] = [[2, -3, 4], [0, 5, -2], [0, -0.2, 4]]
        model.A = sparse.COO.from_numpy(
            np.repeat(a[None, :, :], model.A.shape[0], axis=0)
        )
        b = np.zeros(model.B.shape)
        b[:, :3, 0] = [1, 10, 20]
        model.B = sparse.COO.from_numpy(b)
        model.temporal_technosphere_exchanges = {}
        model.temporal_biosphere_exchanges = {}
        yield model, a


@pytest.mark.parametrize("parent", [2.0, -2.0])
@pytest.mark.parametrize("coefficient", [-3.0, 0.0, 3.0])
def test_child_amount_preserves_signed_matrix_convention(parent, coefficient):
    assert Trails._child_amount(parent, coefficient) == -parent * coefficient


@pytest.mark.parametrize("attribute_to_roots", [False, True])
@pytest.mark.parametrize("solver_mode", ["direct", "iterative"])
@pytest.mark.parametrize("depth", [0, 1, 3, 6, "adaptive", "pruned"])
@pytest.mark.parametrize("demand", [2.0, -2.0])
def test_signed_routing_matches_matrix_inventory(
    signed_model, attribute_to_roots, solver_mode, depth, demand
):
    model, a = signed_model
    routing = {"max_depth": depth}
    if isinstance(depth, str):
        routing = dict(
            max_depth=None,
            adaptive_methods=[
                "CML v4.8 2016 no LT - climate change no LT - "
                "global warming potential (GWP100) no LT"
            ],
            adaptive_relative_score_cutoff=1e-4 if depth == "adaptive" else 10.0,
        )
    model.temporal_routing(
        start_year=2005,
        start_act_idx=0,
        amount=demand,
        **routing,
        attribute_to_roots=attribute_to_roots,
        show_progress=False,
    )
    model.lci(
        solver_mode=solver_mode,
        attribute_to_roots=attribute_to_roots,
        iterative_rtol=1e-10,
        show_progress=False,
    )
    result = model.inventory.sum(
        dim=[d for d in model.inventory.dims if d != "activity"]
    )
    values = result.data
    if hasattr(values, "compute"):
        values = values.compute()
    if hasattr(values, "todense"):
        values = values.todense()
    fu = np.zeros(a.shape[0])
    fu[0] = demand
    expected = np.linalg.solve(a.T, fu) * np.asarray(model.B[0].todense()).sum(axis=1)
    np.testing.assert_allclose(values, expected, rtol=1e-6, atol=1e-8)


@pytest.mark.parametrize("source", ["port", "matrix"])
@pytest.mark.parametrize("parent", [1.0, -1.0])
def test_signed_temporal_pulses_match_fast_and_regular_expansion(
    signed_model, source, parent
):
    model, _ = signed_model
    a = model.A.todense()
    # The supplier becomes an input rather than a credit in the later matrix.
    a[model.scenario_index["2020"], 0, 2] = -8
    model.A = sparse.COO.from_numpy(a)
    model.temporal_technosphere_exchanges = {
        ("2005", 0, 2): TemporalExchange(
            distribution=6,
            loc=None,
            scale=None,
            offset_min=0,
            offset_max=15,
            amount_source=source,
            offsets=[0, 15],
            weights=[0.25, 0.75],
        )
    }
    expected = {
        (2005, 1): parent * 3 / 5,
        (2005, 2): -parent * 0.25,
        (2020, 2): parent * (1.5 if source == "matrix" else -0.75),
    }
    fast = model._expand_temporal_child_demands_fast(
        year=2005,
        act_idx=0,
        amount=parent,
        use_temporal_distributions=True,
    )
    regular = model.expand_temporal_exchanges(2005, 0, parent)
    flat = {
        (year, act): value
        for year, acts in regular.items()
        for act, value in acts.items()
    }
    assert fast == pytest.approx(expected)
    assert flat == pytest.approx(expected)
    # Reuse the cached fast template with the opposite signed parent demand.
    again = model._expand_temporal_child_demands_fast(
        year=2005,
        act_idx=0,
        amount=-parent,
        use_temporal_distributions=True,
    )
    assert again == pytest.approx({key: -value for key, value in expected.items()})
