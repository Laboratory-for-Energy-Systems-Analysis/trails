from __future__ import annotations

from pathlib import Path

import numpy as np
import sparse

from trails.importer import import_excel_inventory

from pathlib import Path as _Path
import importlib.util as _importlib_util

_helpers_path = _Path(__file__).resolve().parent / "helpers.py"
_spec = _importlib_util.spec_from_file_location("tests_helpers", _helpers_path)
_mod = _importlib_util.module_from_spec(_spec)  # type: ignore[arg-type]
assert _spec and _spec.loader
_spec.loader.exec_module(_mod)

DummyTrails = _mod.DummyTrails
_install_fake_bw2io = _mod._install_fake_bw2io
_sample_importer_data = _mod._sample_importer_data


def test_import_excel_inventory_all_template_years(tmp_path: Path) -> None:
    trails = DummyTrails()
    data = _sample_importer_data()
    _install_fake_bw2io(data)

    inv_path = tmp_path / "inv.xlsx"
    inv_path.write_text("stub")

    result = import_excel_inventory(trails, inv_path)

    assert result["new_activities"] == 2
    assert result["new_flows"] == 1

    dense_a = trails.A.todense()
    dense_b = trails.B.todense()

    # activity indices: A->0, B->1. flow index: CO2->0
    assert dense_a.shape == (3, 2, 2)
    assert dense_b.shape == (3, 2, 1)

    # Production exchange (A->A) and technosphere (A->B) should be interpolated
    np.testing.assert_allclose(dense_a[:, 0, 0], np.array([1.0, 1.0, 1.0]))
    np.testing.assert_allclose(dense_a[:, 0, 1], np.array([-2.0, -2.0, -2.0]))

    # Biosphere exchange should be interpolated across all years
    np.testing.assert_allclose(dense_b[:, 0, 0], np.array([3.0, 3.0, 3.0]))

    # Temporal exchange stored on template labels
    assert ("2000", 0, 0) in trails.temporal_technosphere_exchanges
    assert ("2010", 0, 0) in trails.temporal_technosphere_exchanges

    # Indices are consistent across template labels
    assert trails.activity_indices["2000"] == trails.activity_indices["2010"]
    assert trails.biosphere_indices["2000"] == trails.biosphere_indices["2010"]


def test_import_excel_inventory_single_year(tmp_path: Path) -> None:
    trails = DummyTrails()
    data = _sample_importer_data()
    _install_fake_bw2io(data)

    inv_path = tmp_path / "inv.xlsx"
    inv_path.write_text("stub")

    result = import_excel_inventory(trails, inv_path, year=2005)

    assert result["new_activities"] == 2
    assert result["new_flows"] == 1

    dense_a = trails.A.todense()
    dense_b = trails.B.todense()

    # Only year 2005 should be populated (scenario index 1)
    assert dense_a.shape == (3, 2, 2)
    assert dense_b.shape == (3, 2, 1)

    np.testing.assert_allclose(dense_a[1, 0, 0], 1.0)
    np.testing.assert_allclose(dense_a[1, 0, 1], -2.0)
    np.testing.assert_allclose(dense_b[1, 0, 0], 3.0)

    assert np.all(dense_a[0] == 0)
    assert np.all(dense_a[2] == 0)
    assert np.all(dense_b[0] == 0)
    assert np.all(dense_b[2] == 0)
