"""Independent analytical protection for the real-pilot equivalence verifier."""

import numpy as np
import pytest
from scipy.sparse import csc_matrix

from dev.stock_vintage.check_water_pilot import (
    prove_scoped_rewrite,
    stream_inventory_totals,
)
from trails.factorized_inventory import FactorizedInventoryBuilder
from trails.chunked_inventory import ChunkedInventoryBuilder


def example():
    names = [
        "service",
        "capital",
        "disposal",
        "service [stock context water: service]",
        "capital [stock context water: capital]",
        "disposal [stock context water: lifecycle-test]",
    ]
    indices = {
        i: {
            "name": n,
            "reference product": names[i if i < 3 else i - 3],
            "unit": "unit",
            "location": "GLO",
        }
        for i, n in enumerate(names)
    }
    a = np.eye(6)
    a[0, 0] = a[3, 3] = 0.8
    a[1, 0] = a[4, 3] = -0.4
    a[2, 1] = 2
    a[5, 3] = 0.8
    a[2, 5] = -1
    b = np.array([[0, 3, 0, 0, 3, 0], [0, 0, 5, 0, 0, 0]])
    audit = {"context_id": "water", "caller": indices[3], "original_caller": indices[0]}
    return a, b, indices, audit


def test_scoped_operator_proof_checks_own_product_losses_and_lifted_sign():
    a, b, indices, audit = example()
    result = prove_scoped_rewrite(csc_matrix(a), csc_matrix(b), indices, audit)
    assert result["static_equivalence_passed"]
    assert result["net_caller_production"] == 0.8
    for index in (0, 3):
        demand = np.zeros(6)
        demand[index] = 1
        assert b @ np.linalg.solve(a, demand) == pytest.approx([1.5, -5])
    a[5, 3] *= -1
    with pytest.raises(AssertionError, match="changed the source inventory"):
        prove_scoped_rewrite(csc_matrix(a), csc_matrix(b), indices, audit)


def test_scoped_operator_proof_rejects_lost_distribution_loss_and_biosphere():
    a, b, indices, audit = example()
    a[3, 3] = 1
    with pytest.raises(AssertionError):
        prove_scoped_rewrite(csc_matrix(a), csc_matrix(b), indices, audit)
    a, b, indices, audit = example()
    b[0, 4] = 0
    with pytest.raises(AssertionError):
        prove_scoped_rewrite(csc_matrix(a), csc_matrix(b), indices, audit)


@pytest.mark.parametrize("backend", ["factorized", "chunked"])
@pytest.mark.parametrize("sign", [1, -2])
def test_streamed_full_flow_totals_and_dated_roles_match_analytical_inventory(
    backend, sign
):
    cls = (
        FactorizedInventoryBuilder
        if backend == "factorized"
        else ChunkedInventoryBuilder
    )
    builder = cls(
        n_activities=3,
        n_flows=3,
        n_years=2,
        has_root=True,
        value_dtype=np.float64,
        memory_budget=2**20,
    )
    try:
        if backend == "factorized":
            builder.append_factor(
                year_index=0,
                activities=np.array([0, 1, 1]),
                flows=np.array([0, 0, 1]),
                biosphere_values=np.array([2, 3, 4]) * sign,
                supply_matrix=np.array([[1, 2], [0.5, 0], [0, 1]]),
                roots=np.array([0, 2]),
            )
        else:
            builder.append(
                np.array([0, 0, 1, 1]),
                np.array([0, 0, 0, 1]),
                np.array([0, 0, 0, 0]),
                np.array([2, 4, 1.5, 2]) * sign,
                roots=np.array([0, 2, 0, 0]),
            )
        builder.append(
            np.array([1]),
            np.array([1]),
            np.array([1]),
            np.array([-7.0]) * sign,
            roots=np.array([2]),
        )
        builder.finalize()
        totals, roles, entries = stream_inventory_totals(
            builder, np.array([0, 1, 2]), 1, {0: "construction", 2: "disposal"}, 2
        )
        assert totals == pytest.approx(np.array([7.5, -5, 0]) * sign)
        assert roles["construction"] == pytest.approx(np.array([2, 0]) * sign)
        assert roles["disposal"] == pytest.approx(np.array([0, -7]) * sign)
        assert entries == 5
        assert builder._activity_reduction_cache == {}
    finally:
        builder.close()
