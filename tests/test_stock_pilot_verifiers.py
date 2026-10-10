"""Independent analytical protection for the real-pilot equivalence verifier."""

import numpy as np
import pytest
from types import SimpleNamespace
from scipy.sparse import csc_matrix

from dev.stock_vintage.check_water_pilot import (
    check_direct_biosphere_dates,
    prove_scoped_rewrite,
    stream_inventory_totals,
)
from dev.stock_vintage.check_ccgt_pilot import identity
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


@pytest.mark.parametrize("amount", [1, -2])
@pytest.mark.parametrize("event_index", [0, 1])
def test_land_occupation_is_checked_at_service_year_with_own_product_loss(
    amount, event_index
):
    supplier = {
        "name": "occupation",
        "reference product": "occupation",
        "unit": "square meter-year",
        "location": "GLO",
    }
    builder = ChunkedInventoryBuilder(
        n_activities=3,
        n_flows=2,
        n_years=2,
        has_root=True,
        value_dtype=np.float64,
        memory_budget=2**20,
    )
    try:
        # Unrelated background occupation is deliberately at another date/root.
        builder.append(
            np.array([0, 2]),
            np.array([1, 1]),
            np.array([0, event_index]),
            np.array([0.7, 0.3 * amount]),
            roots=np.array([0, 2]),
        )
        builder.finalize()
        model = SimpleNamespace(_inventory_builder=builder)
        inventory = SimpleNamespace(
            coords={"year": SimpleNamespace(values=np.array([2000, 2023]))}
        )
        b = csc_matrix(np.array([[0, 0, 0], [1, 0, 1]]))
        audit = {
            "lifted_biosphere": [
                {"supplier": supplier, "amount_per_calling_dataset": 0.24}
            ]
        }
        args = (model, inventory, b, {identity(supplier): 2}, audit, 2023, amount, 0.8)
        if event_index == 0:
            with pytest.raises(AssertionError, match="dated outside service"):
                check_direct_biosphere_dates(*args)
        else:
            result = check_direct_biosphere_dates(*args)
            assert result[0]["events"] == pytest.approx({2023: 0.3 * amount})
    finally:
        builder.close()


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
