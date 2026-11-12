"""Empirical profiles must retain their dates beyond the default time margin."""

from types import SimpleNamespace

import numpy as np
import pytest
import sparse

from trails import Trails
from trails.temporal_distributions import TemporalExchange


@pytest.mark.parametrize(
    "offsets,declared,expected",
    [
        ([975], (None, None), (0, 975)),
        ([-750], (None, None), (-750, 0)),
        ([-750.5, 975.5], (0, 1), (-751, 976)),
        ([], (-10, 20), (-10, 20)),
        ([np.nan, np.inf, 15], (None, None), (0, 15)),
        (None, (None, None), (0, 0)),
    ],
)
def test_explicit_offsets_extend_declared_support(offsets, declared, expected):
    tex = SimpleNamespace(
        offsets=offsets, offset_min=declared[0], offset_max=declared[1]
    )
    assert Trails._exchange_offset_bounds({"profile": tex}) == expected


@pytest.mark.parametrize("roots", [False, True])
@pytest.mark.parametrize("store_inventory", [False, True])
@pytest.mark.parametrize("declared", [(None, None), (0, 1)])
def test_long_empirical_emissions_keep_dates_and_amounts(
    example_package, roots, store_inventory, declared
):
    model = Trails(example_package, interpolate_annual=True, cache_interpolation=False)
    try:
        a = np.repeat(np.eye(model.A.shape[1])[None, :, :], model.A.shape[0], axis=0)
        b = np.zeros(model.B.shape)
        b[:, 0, 0] = 1
        model.A = sparse.COO.from_numpy(a)
        model.B = sparse.COO.from_numpy(b)
        model.temporal_technosphere_exchanges = {}
        model.temporal_biosphere_exchanges = {
            (year, 0, 0): TemporalExchange(
                distribution=6,
                loc=None,
                scale=None,
                offset_min=declared[0],
                offset_max=declared[1],
                amount_source="port",
                offsets=[-750, 975],
                weights=[0.25, 0.75],
            )
            for year in model.scenario_index
        }
        model.temporal_routing(
            start_year=2025,
            start_act_idx=0,
            amount=1,
            max_depth=0,
            attribute_to_roots=roots,
            show_progress=False,
        )
        model.lca(
            methods=[
                "CML v4.8 2016 no LT - climate change no LT - global warming potential (GWP100) no LT"
            ],
            store_inventory=store_inventory,
            compute_score=True,
            solver_mode="direct",
            attribute_to_roots=roots,
            show_progress=False,
        )
        arrays = [model.scores]
        if store_inventory:
            arrays.append(model.inventory)
        for array in arrays:
            by_year = array.sum(dim=[d for d in array.dims if d != "year"])
            values = by_year.data
            if hasattr(values, "compute"):
                values = values.compute()
            if hasattr(values, "todense"):
                values = values.todense()
            nonzero = {
                int(y): float(v) for y, v in zip(by_year.year.values, values) if v != 0
            }
            assert nonzero == pytest.approx({1275: 0.25, 3000: 0.75})
    finally:
        if hasattr(model, "close"):
            model.close()
