from __future__ import annotations

from pathlib import Path
import sys
import types

import numpy as np
import sparse

from trails.importer import import_excel_inventory


class DummyTrails:
    def __init__(self) -> None:
        self.value_dtype = np.float32
        self.index_dtype = np.int32

        self.scenario_labels = ["2000", "2005", "2010"]
        self.scenario_index = {label: i for i, label in enumerate(self.scenario_labels)}
        self.template_labels = ["2000", "2010"]
        self.years_int = np.array([2000, 2005, 2010], dtype=int)

        self.activity_indices: dict[str, dict[int, dict]] = {}
        self.biosphere_indices: dict[str, dict[int, dict]] = {}

        self.temporal_technosphere_exchanges: dict = {}
        self.temporal_biosphere_exchanges: dict = {}

        self.A = sparse.COO(
            coords=[np.array([], dtype=self.index_dtype)] * 3,
            data=np.array([], dtype=self.value_dtype),
            shape=(len(self.scenario_labels), 0, 0),
        )
        self.B = sparse.COO(
            coords=[np.array([], dtype=self.index_dtype)] * 3,
            data=np.array([], dtype=self.value_dtype),
            shape=(len(self.scenario_labels), 0, 0),
        )

        self._A_row_cache: dict = {}
        self._direct_bio_cache_by_year: dict = {}
        self._tech_td_cache: dict = {}
        self._tech_td_expanded_cache: dict = {}
        self._td_offsets_cache: dict = {}

        self.min_year = 2000
        self.max_year = 2010

    def _map_year_to_template_year(self, year: int) -> int:
        return 2000 if year < 2005 else 2010

    def _get_scenario_context(self, year: int):
        year = int(year)
        label = str(year)
        if label not in self.scenario_index:
            return None
        return year, label, self.scenario_index[label]


def _install_fake_bw2io(
    data: list[dict] | None = None,
    *,
    data_by_path: dict[str, list[dict]] | None = None,
) -> None:
    bw2io = types.ModuleType("bw2io")
    importers = types.ModuleType("bw2io.importers")
    excel = types.ModuleType("bw2io.importers.excel")

    class ExcelImporter:
        _data: list[dict] = []
        _data_by_path: dict[str, list[dict]] = {}

        def __init__(self, path: str) -> None:
            self.path = Path(path)
            self.data = list(self._data_by_path.get(str(self.path), self._data))

        def apply_strategies(self) -> None:
            return None

    ExcelImporter._data = data or []
    ExcelImporter._data_by_path = data_by_path or {}

    excel.ExcelImporter = ExcelImporter
    importers.excel = excel
    bw2io.importers = importers

    sys.modules["bw2io"] = bw2io
    sys.modules["bw2io.importers"] = importers
    sys.modules["bw2io.importers.excel"] = excel


def _sample_importer_data() -> list[dict]:
    return [
        {
            "name": "A",
            "reference product": "A",
            "location": "GLO",
            "unit": "kg",
            "database": "db",
            "code": "A",
            "exchanges": [
                {
                    "type": "production",
                    "name": "A",
                    "reference product": "A",
                    "location": "GLO",
                    "amount": 1.0,
                    "temporal_distribution": 3,
                    "temporal_loc": 0,
                    "temporal_scale": 1,
                    "temporal_min": -1,
                    "temporal_max": 1,
                    "temporal_amount_source": "port",
                },
                {
                    "type": "technosphere",
                    "name": "B",
                    "reference product": "B",
                    "location": "GLO",
                    "amount": 2.0,
                },
                {
                    "type": "biosphere",
                    "name": "CO2",
                    "categories": ("air", "urban"),
                    "amount": 3.0,
                },
            ],
        },
        {
            "name": "B",
            "reference product": "B",
            "location": "GLO",
            "unit": "kg",
            "database": "db",
            "code": "B",
            "exchanges": [
                {
                    "type": "production",
                    "name": "B",
                    "reference product": "B",
                    "location": "GLO",
                    "amount": 1.0,
                }
            ],
        },
    ]


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


def test_import_excel_inventory_year_specific_amounts(tmp_path: Path) -> None:
    trails = DummyTrails()
    data = [
        {
            "name": "A",
            "reference product": "A",
            "location": "GLO",
            "unit": "kg",
            "database": "db",
            "code": "A",
            "exchanges": [
                {
                    "type": "production",
                    "name": "A",
                    "reference product": "A",
                    "location": "GLO",
                    "amount": 0.0,
                    "2000": 1.0,
                    "2010": 3.0,
                },
                {
                    "type": "technosphere",
                    "name": "B",
                    "reference product": "B",
                    "location": "GLO",
                    "amount": 0.0,
                    "2000": 2.0,
                    "2010": 4.0,
                },
            ],
        },
        {
            "name": "B",
            "reference product": "B",
            "location": "GLO",
            "unit": "kg",
            "database": "db",
            "code": "B",
            "exchanges": [
                {
                    "type": "production",
                    "name": "B",
                    "reference product": "B",
                    "location": "GLO",
                    "amount": 1.0,
                }
            ],
        },
    ]
    _install_fake_bw2io(data)

    inv_path = tmp_path / "inv.xlsx"
    inv_path.write_text("stub")

    import_excel_inventory(trails, inv_path)

    dense_a = trails.A.todense()

    # Production exchange interpolated: 1.0 (2000) -> 3.0 (2010) => 2.0 at 2005
    np.testing.assert_allclose(dense_a[:, 0, 0], np.array([1.0, 2.0, 3.0]))
    # Technosphere exchange interpolated and sign-flipped
    np.testing.assert_allclose(dense_a[:, 0, 1], np.array([-2.0, -3.0, -4.0]))


def test_import_excel_inventory_multiple_files_concatenated(tmp_path: Path) -> None:
    trails = DummyTrails()

    inv1_path = tmp_path / "inv1.xlsx"
    inv2_path = tmp_path / "inv2.xlsx"
    inv1_path.write_text("stub")
    inv2_path.write_text("stub")

    data_by_path = {
        str(inv1_path): [
            {
                "name": "X",
                "reference product": "X",
                "location": "GLO",
                "unit": "kg",
                "database": "db",
                "code": "X",
                "exchanges": [
                    {
                        "type": "production",
                        "name": "X",
                        "reference product": "X",
                        "location": "GLO",
                        "amount": 1.0,
                    },
                    {
                        "type": "technosphere",
                        "name": "Y",
                        "reference product": "Y",
                        "location": "GLO",
                        "amount": 2.0,
                    },
                ],
            }
        ],
        str(inv2_path): [
            {
                "name": "Y",
                "reference product": "Y",
                "location": "GLO",
                "unit": "kg",
                "database": "db",
                "code": "Y",
                "exchanges": [
                    {
                        "type": "production",
                        "name": "Y",
                        "reference product": "Y",
                        "location": "GLO",
                        "amount": 1.0,
                    }
                ],
            }
        ],
    }
    _install_fake_bw2io(data_by_path=data_by_path)

    result = import_excel_inventory(trails, [str(inv1_path), str(inv2_path)], year=2005)

    assert result["unlinked"] == 0
    assert result["new_activities"] == 2
    assert result["production"] == 2
    assert result["technosphere"] == 1

    dense_a = trails.A.todense()
    np.testing.assert_allclose(dense_a[1, 0, 0], 1.0)  # X production
    np.testing.assert_allclose(dense_a[1, 0, 1], -2.0)  # X -> Y technosphere
    np.testing.assert_allclose(dense_a[1, 1, 1], 1.0)  # Y production
