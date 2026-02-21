from __future__ import annotations

from pathlib import Path
import sys
import types
import warnings
from typing import Any

import numpy as np
import pandas as pd
import pytest
import sparse
import xarray as xr

from trails.cache_interpolation import cache_dir_for_package
from trails.fair_rf import _sanitize_emissions_year_values, run_fair_delta_rf
from trails.lca import lca_static
from trails.plotting import plot_temporal_sankey_graphlike


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
        return 2000 if year < 2005 else 2010

    def _get_scenario_context(self, year: int):
        year = int(year)
        label = str(year)
        if label not in self.scenario_index:
            return None
        return year, label, self.scenario_index[label]


def _install_fake_bw2io(data: list[dict]) -> None:
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


class DummyTrailsStatic:
    """Minimal Trails stub for static LCA tests."""

    def __init__(self) -> None:
        self.inventory = None
        self.characterized_inventory = None
        self.static_score = None
        self.A = sparse.COO(
            coords=[np.array([], dtype=int)] * 3, data=np.array([]), shape=(1, 2, 2)
        )
        self.B = sparse.COO(
            coords=[np.array([], dtype=int)] * 3, data=np.array([]), shape=(1, 2, 3)
        )

    def reset_inventory(self, reset_scores: bool = False) -> None:
        """Reset inventory for testing."""
        self.inventory = None
        if reset_scores:
            self.characterized_inventory = None


class DummyTrailsSankey:
    """Minimal Trails stub for Sankey plotting."""

    def __init__(self) -> None:
        import networkx as nx

        self.graph = nx.DiGraph()
        # Nodes: (year, act_idx) with depth
        self.graph.add_node((2050, 0), year=2050, depth=0, act_idx=0)
        self.graph.add_node((2051, 1), year=2051, depth=1, act_idx=1)
        self.graph.add_edge((2050, 0), (2051, 1), amount=1.0)

        self.activity_indices = {
            "2050": {
                0: {
                    "name": "Root",
                    "reference product": "Root product",
                    "location": "GLO",
                },
                1: {
                    "name": "Child",
                    "reference product": "Child product",
                    "location": "GLO",
                },
            }
        }
        self.characterized_inventory = None


class DummyFairRun:
    """Minimal FaIR run stub with forcing and temperature arrays."""

    def __init__(self, scenario: str, config_names: list[str], values: float) -> None:
        timebounds = np.array([2000.5, 2001.5], dtype=float)
        forcing = xr.DataArray(
            np.full((1, len(config_names), len(timebounds), 1), values, dtype=float),
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [scenario],
                "config": config_names,
                "timebounds": timebounds,
                "specie": ["CO2"],
            },
        )
        temperature = xr.DataArray(
            np.full((1, len(config_names), len(timebounds)), values, dtype=float),
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [scenario],
                "config": config_names,
                "timebounds": timebounds,
            },
        )
        self.forcing = forcing
        self.temperature = temperature


class DummyTrailsFair:
    """Minimal Trails stub for FaIR quantile tests."""

    def __init__(self) -> None:
        self.debug = False
        # Inventory with dims (activity, flow, year, root activity)
        data = np.array([1.0, 2.0], dtype=float)
        coords = np.array([[0, 0], [0, 0], [0, 1], [0, 0]], dtype=int)
        inv = sparse.COO(coords=coords, data=data, shape=(1, 1, 2, 1))
        self.inventory = xr.DataArray(
            inv,
            dims=("activity", "flow", "year", "root activity"),
            coords={
                "activity": [0],
                "flow": [0],
                "year": [2000, 2001],
                "root activity": [0],
            },
        )
        self.biosphere_indices = {
            "2000": {
                0: {
                    "name": "CO2",
                    "compartment": "air",
                    "subcompartment": "",
                }
            }
        }


def test_static_score_preserves_method_order(monkeypatch: pytest.MonkeyPatch) -> None:
    trails = DummyTrailsStatic()

    class DummyDicts:
        def __init__(self) -> None:
            self.activity = {0: 0, 1: 1}
            self.biosphere = {0: 0, 1: 1, 2: 2}

    class DummyLCA:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self.dicts = DummyDicts()
            self.demand = {0: 1.0}
            # inventory: flow x act
            from scipy import sparse as sp

            rows = np.array([0, 1, 2], dtype=int)
            cols = np.array([0, 0, 1], dtype=int)
            data = np.array([1.0, 2.0, 3.0], dtype=float)
            self.inventory = sp.coo_matrix((data, (rows, cols)), shape=(3, 2)).tocsr()

        def lci(self) -> None:
            return None

    def fake_build_dp(*args: Any, **kwargs: Any) -> tuple[Any, Any, Any, Any]:
        return None, None, None, None

    def fake_get_cf_vector(*args: Any, **kwargs: Any) -> np.ndarray:
        methods = kwargs.get("methods") or []
        method = methods[0] if methods else ""
        if method == "m1":
            return np.array([1.0, 0.0, 0.0], dtype=float)
        return np.array([0.0, 1.0, 0.0], dtype=float)

    import importlib

    lca_module = importlib.import_module("trails.lca")
    monkeypatch.setattr(
        lca_module, "build_datapackage_for_year_from_trails", fake_build_dp
    )
    monkeypatch.setattr(lca_module.bc, "LCA", DummyLCA)
    monkeypatch.setattr(lca_module, "get_cf_vector", fake_get_cf_vector)
    monkeypatch.setattr(
        lca_module,
        "build_characterized_inventory",
        lambda *args, **kwargs: xr.DataArray(
            np.zeros((2, 3, 1), dtype=float),
            dims=("method", "flow", "year"),
            coords={"method": ["m1", "m2"], "flow": [0, 1, 2], "year": [2000]},
        ),
    )

    lca_static(trails, year=2000, fu_act_idx=0, methods=["m1", "m2"], amount=1.0)

    assert isinstance(trails.static_score, list)
    assert trails.static_score == pytest.approx([1.0, 2.0])


def test_importer_unlinked_exchange_raises(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from trails.importer import import_excel_inventory

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
                    "amount": 1.0,
                },
                {
                    "type": "technosphere",
                    "name": "MISSING",
                    "reference product": "MISSING",
                    "location": "GLO",
                    "amount": 1.0,
                },
            ],
        }
    ]
    _install_fake_bw2io(data)
    trails = DummyTrails()

    inv_path = tmp_path / "inv.xlsx"
    inv_path.write_text("stub")

    with pytest.raises(ValueError, match="Unlinked exchanges detected"):
        import_excel_inventory(trails, inv_path)

    out = capsys.readouterr().out
    assert "Unlinked exchanges" in out
    assert "MISSING" in out


def test_fair_quantile_dimensions(monkeypatch: pytest.MonkeyPatch) -> None:
    trails = DummyTrailsFair()

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["CO2"],
                "unit": ["Gt CO2/yr"],
                "2000.5": [0.0],
                "2001.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {("CO2", "air", ""): "CO2"}, {}

    call_state = {"count": 0}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> DummyFairRun:
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1", "c2"]
        values = 0.0 if call_state["count"] == 0 else 1.0
        call_state["count"] += 1
        return DummyFairRun(str(scenario), list(config_names), values)

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    rf = run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1", "c2"],
        per_species_runs=True,
        validate_emissions_delta=False,
        scale_factor=1.0,
    )

    assert rf.dims == ("quantile", "year", "flow", "root activity")
    assert trails.delta_temperature.dims == (
        "quantile",
        "year",
        "flow",
        "root activity",
    )
    assert rf.shape[0] == 5
    assert int(getattr(rf.data, "nnz", 0)) > 0


def test_sanitize_emissions_year_values_fills_missing() -> None:
    df = pd.DataFrame(
        {
            "scenario": ["s"],
            "region": ["World"],
            "variable": ["CH4"],
            "unit": ["Mt CH4/yr"],
            "2000.5": [1.0],
            "2001.5": [np.nan],
            "2002.5": ["3.5"],
        }
    )

    out = _sanitize_emissions_year_values(df)

    assert float(out.loc[0, "2000.5"]) == pytest.approx(1.0)
    assert float(out.loc[0, "2001.5"]) == pytest.approx(0.0)
    assert float(out.loc[0, "2002.5"]) == pytest.approx(3.5)


def test_fair_quantiles_suppress_all_nan_warning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    trails = DummyTrailsFair()

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["CO2"],
                "unit": ["Gt CO2/yr"],
                "2000.5": [0.0],
                "2001.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {("CO2", "air", ""): "CO2"}, {}

    class DummyFairRunWithNaN:
        def __init__(
            self, scenario: str, config_names: list[str], baseline: bool
        ) -> None:
            timebounds = np.array([2000.5, 2001.5], dtype=float)
            if baseline:
                vals = np.array([[0.0, np.nan], [0.0, np.nan]], dtype=float)
            else:
                vals = np.array([[1.0, np.nan], [2.0, np.nan]], dtype=float)
            forcing = xr.DataArray(
                vals[None, :, :, None],
                dims=("scenario", "config", "timebounds", "specie"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": timebounds,
                    "specie": ["CO2"],
                },
            )
            temperature = xr.DataArray(
                vals[None, :, :],
                dims=("scenario", "config", "timebounds"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": timebounds,
                },
            )
            self.forcing = forcing
            self.temperature = temperature

    call_state = {"count": 0}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> DummyFairRunWithNaN:
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1", "c2"]
        baseline = call_state["count"] == 0
        call_state["count"] += 1
        return DummyFairRunWithNaN(str(scenario), list(config_names), baseline=baseline)

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        run_fair_delta_rf(
            trails,
            scenario="s",
            config_names=["c1", "c2"],
            per_species_runs=True,
            validate_emissions_delta=False,
            scale_factor=1.0,
        )

    assert not any("All-NaN slice encountered" in str(w.message) for w in caught)


def test_sankey_graphlike_writes_html(tmp_path: Path) -> None:
    trails = DummyTrailsSankey()
    out = tmp_path / "sankey.html"

    fig = plot_temporal_sankey_graphlike(
        trails,
        edge_weight="amount",
        orientation="depth_x_year_y",
        y_padding=0.03,
        filename=str(out),
    )

    assert fig is not None
    assert out.exists()
    html = out.read_text()
    assert "nouislider" in html.lower()


def test_cache_dir_uses_short_hash(example_package: Any) -> None:
    cache_dir = cache_dir_for_package(
        example_package, value_dtype="float32", index_dtype="int32"
    )
    name = cache_dir.name
    assert name.startswith("interp_")
    suffix = name.split("interp_", 1)[1]
    assert len(suffix) == 12
