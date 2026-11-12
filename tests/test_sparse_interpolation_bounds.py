"""Sparse anchors must interpolate as zero outside their nonzero support."""

import numpy as np
import pytest
import sparse

from trails.datapackage import _interp_slice_union_vectorized


@pytest.mark.parametrize("value_dtype", [np.float32, np.float64])
@pytest.mark.parametrize("index_dtype", [np.int32, np.int64])
@pytest.mark.parametrize("weight", [0.0, 0.25, 0.5, 1.0])
@pytest.mark.parametrize(
    "left,right",
    [
        ([[1, 0, 0]], [[0, 0, 2]]),
        ([[0, 0, 2]], [[1, 0, 0]]),
        ([[1, 0, 3], [0, -4, 0]], [[0, 2, -3], [5, 0, 0]]),
        ([[0, 0, 0]], [[1, 0, 2]]),
        ([[1, 0, 2]], [[0, 0, 0]]),
        ([[0, 0, 0]], [[0, 0, 0]]),
        ([[1, -2, 3]], [[-1, 2, -3]]),
    ],
    ids=[
        "appearing_last_coordinate",
        "disappearing_last_coordinate",
        "different_supports",
        "empty_left",
        "empty_right",
        "both_empty",
        "signed_cancellation",
    ],
)
def test_sparse_interpolation_matches_dense_reference(
    left, right, weight, index_dtype, value_dtype
):
    """Check missing entries, endpoints, signs, precision and sparse pruning."""
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    actual = _interp_slice_union_vectorized(
        sparse.COO.from_numpy(left),
        sparse.COO.from_numpy(right),
        weight,
        np.dtype(index_dtype),
        np.dtype(value_dtype),
    )
    expected = ((1.0 - weight) * left + weight * right).astype(value_dtype)
    np.testing.assert_array_equal(actual.todense(), expected)
    assert actual.dtype == value_dtype
    if actual.nnz:
        assert actual.coords.dtype == index_dtype
    assert actual.nnz == np.count_nonzero(expected)
