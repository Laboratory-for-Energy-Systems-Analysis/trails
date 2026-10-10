"""Independent analytical protection for the real-pilot equivalence verifier."""

import numpy as np
import pytest
from scipy.sparse import csc_matrix

from dev.stock_vintage.check_water_pilot import prove_scoped_rewrite


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
