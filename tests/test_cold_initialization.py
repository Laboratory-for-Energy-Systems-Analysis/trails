"""Regression coverage for cold-loading and annual interpolation optimizations."""

import numpy as np
import pandas as pd
import pytest
import sparse

from trails import datapackage as dp


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("weight", [0.0, 0.25, 0.5, 1.0])
def test_union_interpolation_different_sparsity(dtype, weight):
    # Different supports include a key beyond either anchor's last key;
    # overlapping opposite values cancel exactly at the midpoint.
    left = np.array([[0, 2, 0], [0, -4, 0]], dtype=dtype)
    right = np.array([[3, -2, 0], [0, 0, 8]], dtype=dtype)
    actual = dp._interp_slice_union_vectorized(
        sparse.COO.from_numpy(left),
        sparse.COO.from_numpy(right),
        weight,
        np.int32,
        dtype,
    )
    expected = (
        (1 - weight) * left.astype(np.float64) + weight * right.astype(np.float64)
    ).astype(dtype)
    np.testing.assert_array_equal(actual.todense(), expected)
    assert actual.dtype == dtype
    assert actual.coords.dtype == np.int32
    assert actual.nnz == np.count_nonzero(expected)


@pytest.mark.parametrize("empty_left", [True, False])
def test_union_interpolation_empty_anchor(empty_left):
    zero = sparse.COO.from_numpy(np.zeros((2, 3)))
    full = sparse.COO.from_numpy(np.array([[1.0, 0, 0], [0, 0, 2]]))
    left, right = (zero, full) if empty_left else (full, zero)
    actual = dp._interp_slice_union_vectorized(left, right, 0.5, np.int32, np.float64)
    np.testing.assert_array_equal(actual.todense(), full.todense() * 0.5)


def test_inventory_reader_preserves_consumed_fields(tmp_path):
    rows = []
    for index, value, dist in [
        (" 0 ", " 1e-9 ", ""),
        ("1.0", "invalid", "6"),
        ("2", "", "0"),
        ("\u00a03\u00a0", "\u20031.25\u2003", "\u00a04\u00a0"),
    ]:
        row = dict.fromkeys(dp.BASE_COLS_B + dp.TEMPORAL_COLS, "")
        row.update(
            {
                "index of activity": index,
                "index of biosphere flow": " 2.0 ",
                "value": value,
                "temporal_distribution": dist,
                "temporal_offsets": "[0, 2]",
                "temporal_weights": "[0.25, 0.75]",
                "temporal_amount_source": " matrix ",
            }
        )
        rows.append(row)
    path = tmp_path / "B.csv"
    pd.DataFrame(rows).to_csv(path, sep=";", index=False)
    original = dp._read_matrix_csv_fast(path, kind="B")
    selected = dp._read_matrix_csv_fast(path, kind="B", inventory_only=True)
    pd.testing.assert_frame_equal(selected, original[selected.columns])
    assert selected.loc[3, "index of activity"] == 3
    assert selected.loc[3, "value"] == 1.25
    assert selected.loc[3, "temporal_distribution"] == 4
    # Header validation remains strict even for columns not used in LCI.
    pd.DataFrame(rows).drop(columns="flip").to_csv(path, sep=";", index=False)
    with pytest.raises(ValueError, match="Missing required columns"):
        dp._read_matrix_csv_fast(path, kind="B", inventory_only=True)


def test_annual_interpolation_unsorted_anchors_and_padding():
    # Anchor years arrive out of order; supports change and include an empty year.
    values = np.array([[[0.0, 4.0]], [[2.0, 0.0]], [[0.0, 0.0]]], dtype=np.float64)
    matrix = sparse.COO.from_numpy(values)
    a, b, labels, index = dp.interpolate_to_annual(
        matrix,
        matrix,
        ["2004", "2000", "2002"],
        value_dtype=np.float64,
        start_year_offset=-2,
        end_year_offset=2,
    )
    expected = np.array(
        [
            [[2, 0]],
            [[2, 0]],
            [[2, 0]],
            [[1, 0]],
            [[0, 0]],
            [[0, 2]],
            [[0, 4]],
            [[0, 4]],
            [[0, 4]],
        ]
    )
    assert labels == [str(y) for y in range(1998, 2007)]
    assert index["2002"] == 4
    np.testing.assert_array_equal(a.todense(), expected)
    np.testing.assert_array_equal(b.todense(), expected)


def test_cold_cache_roundtrip_preserves_matrices_and_timing(
    example_package, tmp_path, monkeypatch
):
    import importlib
    from trails import Trails

    cache_module = importlib.import_module("trails.cache_interpolation")
    monkeypatch.setattr(
        cache_module, "cache_dir_for_package", lambda *a, **kw: tmp_path / "cache"
    )
    first = Trails(example_package)
    second = Trails(example_package)
    for name in ["A", "B"]:
        left, right = getattr(first, name), getattr(second, name)
        np.testing.assert_array_equal(left.coords, right.coords)
        np.testing.assert_array_equal(left.data, right.data)
    assert first.scenario_labels == second.scenario_labels
    assert first.template_labels == second.template_labels
    assert first.activity_indices == second.activity_indices
    assert first.biosphere_indices == second.biosphere_indices
    assert (
        first.temporal_technosphere_exchanges == second.temporal_technosphere_exchanges
    )
    assert first.temporal_biosphere_exchanges == second.temporal_biosphere_exchanges
    first.close()
    second.close()


def test_repeated_pulses_keep_independent_lists(tmp_path, monkeypatch):
    from types import SimpleNamespace

    for kind, columns in [("A", dp.BASE_COLS_A), ("B", dp.BASE_COLS_B)]:
        rows = []
        for activity in range(2):
            row = dict.fromkeys(columns, 0)
            row.update(
                {
                    "index of activity": activity,
                    "value": 1,
                    "temporal_distribution": 6,
                    "temporal_offsets": "[0, 2]",
                    "temporal_weights": "[0.25, 0.75]",
                }
            )
            rows.append(row)
        pd.DataFrame(rows).to_csv(tmp_path / f"{kind}_matrix.csv", sep=";", index=False)
    monkeypatch.setattr(
        dp,
        "_iter_inventory_resources",
        lambda package, filename: iter(
            [("2030", SimpleNamespace(source=tmp_path / filename))]
        ),
    )
    _, _, _, _, tech, bio = dp.load_matrices_from_package(None)
    for mapping in [tech, bio]:
        first, second = mapping[("2030", 0, 0)], mapping[("2030", 1, 0)]
        assert first.offsets == second.offsets == [0, 2]
        assert first.weights == second.weights == [0.25, 0.75]
        first.offsets.append(99)
        first.weights[0] = 1
        assert second.offsets == [0, 2]
        assert second.weights == [0.25, 0.75]
