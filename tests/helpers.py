from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import sparse


class DummyTrails:
    """Minimal Trails stub for importer tests."""

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
        """Map to template year."""
        return 2000 if year < 2005 else 2010

    def _get_scenario_context(self, year: int):
        """Return scenario context for a year."""
        year = int(year)
        label = str(year)
        if label not in self.scenario_index:
            return None
        return year, label, self.scenario_index[label]


def _install_fake_bw2io(data: list[dict]) -> None:
    """Install a fake bw2io ExcelImporter in sys.modules."""
    bw2io = types.ModuleType("bw2io")
    importers = types.ModuleType("bw2io.importers")
    excel = types.ModuleType("bw2io.importers.excel")

    class ExcelImporter:
        _data: list[dict] = []

        def __init__(self, path: str) -> None:
            self.path = Path(path)
            self.data = list(self._data)

        def apply_strategies(self) -> None:
            return None

    ExcelImporter._data = data

    excel.ExcelImporter = ExcelImporter
    importers.excel = excel
    bw2io.importers = importers

    sys.modules["bw2io"] = bw2io
    sys.modules["bw2io.importers"] = importers
    sys.modules["bw2io.importers.excel"] = excel


def _sample_importer_data() -> list[dict]:
    """Sample importer data with a linked technosphere exchange."""
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
