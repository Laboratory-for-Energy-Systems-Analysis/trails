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

import trails.fair_rf as fair_rf_module
from trails.fair_io import DEFAULT_MAPPING_YAML, load_species_mapping
from trails.cache_interpolation import cache_dir_for_package
from trails.fair_rf import _sanitize_emissions_year_values, run_fair_delta_rf
from trails.lca import lca_static
from trails.plotting import plot_temporal_sankey_graphlike, plot_temporal_scores


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


class DummyTrailsScores:
    """Minimal Trails stub for temporal score plotting."""

    def __init__(self) -> None:
        years = np.array([2000, 2001, 2002, 2003, 2004, 2005, 2006], dtype=int)
        vals = np.array([[0.0, 0.0, 0.0, 2.0, 3.0, 0.0, 0.0]], dtype=float)
        self.scores = xr.DataArray(
            vals,
            dims=("activity", "year"),
            coords={"activity": [0], "year": years},
        )
        self.characterized_inventory = None
        self.activity_indices = {
            "2005": {
                0: {
                    "name": "A0",
                    "reference product": "P0",
                    "location": "GLO",
                }
            }
        }


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


def test_fair_fast_mode_temperature_not_duplicated_across_species(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsTwoFlows:
        def __init__(self) -> None:
            self.debug = False
            coords = np.array(
                [
                    [0, 0, 0, 0],  # activity
                    [0, 0, 1, 1],  # flow
                    [0, 1, 0, 1],  # year
                    [0, 0, 0, 0],  # root activity
                ],
                dtype=int,
            )
            data = np.array([1.0, 1.0, 1.0, 1.0], dtype=float)
            inv = sparse.COO(coords=coords, data=data, shape=(1, 2, 2, 1))
            self.inventory = xr.DataArray(
                inv,
                dims=("activity", "flow", "year", "root activity"),
                coords={
                    "activity": [0],
                    "flow": [0, 1],
                    "year": [2000, 2001],
                    "root activity": [0],
                },
            )
            self.biosphere_indices = {
                "2000": {
                    0: {"name": "CO2", "compartment": "air", "subcompartment": ""},
                    1: {"name": "CH4", "compartment": "air", "subcompartment": ""},
                }
            }

    trails = DummyTrailsTwoFlows()
    state: dict[str, np.ndarray] = {}

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s", "s"],
                "region": ["World", "World"],
                "variable": ["CO2", "CH4"],
                "unit": ["Gt CO2/yr", "Mt CH4/yr"],
                "2000.5": [0.0, 0.0],
                "2001.5": [0.0, 0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {
            ("CO2", "air", ""): "CO2",
            ("CH4", "air", ""): "CH4",
        }, {}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1"]

        local = fair_rf_module._normalize_emissions_columns(emissions_df)
        local = local[
            (local["scenario"] == scenario) & (local["region"].str.lower() == "world")
        ].copy()
        years = np.array([2000.5, 2001.5], dtype=float)
        species = ["CO2", "CH4"]
        forcing_vals = np.zeros((1, len(config_names), len(years), len(species)))

        for i, specie in enumerate(species):
            rows = local[local["variable"] == specie]
            if rows.empty:
                continue
            row = rows.iloc[0]
            forcing_vals[:, :, :, i] = np.array(
                [
                    float(pd.to_numeric(row.get("2000.5", 0.0), errors="coerce") or 0),
                    float(pd.to_numeric(row.get("2001.5", 0.0), errors="coerce") or 0),
                ],
                dtype=float,
            )[None, None, :]

        temperature_vals = np.sum(forcing_vals, axis=3)
        if "baseline_temp" not in state:
            state["baseline_temp"] = temperature_vals[0, 0, :].copy()
        else:
            state["pert_temp"] = temperature_vals[0, 0, :].copy()

        forcing = xr.DataArray(
            forcing_vals,
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
                "specie": species,
            },
        )
        temperature = xr.DataArray(
            temperature_vals,
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
            },
        )
        return types.SimpleNamespace(
            forcing=forcing,
            temperature=temperature,
            configs=list(config_names),
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1"],
        per_species_runs=False,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    allocated = np.asarray(trails.delta_temperature.data.todense(), dtype=float)
    allocated_total = allocated[0].sum(axis=(1, 2))
    expected_total = state["pert_temp"] - state["baseline_temp"]

    np.testing.assert_allclose(
        allocated_total,
        expected_total,
        rtol=1e-12,
        atol=1e-15,
    )


def test_fair_fast_mode_keeps_tiny_temperature_signal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsTwoFlowsTiny:
        def __init__(self) -> None:
            self.debug = False
            coords = np.array(
                [
                    [0, 0, 0, 0],  # activity
                    [0, 0, 1, 1],  # flow
                    [0, 1, 0, 1],  # year
                    [0, 0, 0, 0],  # root activity
                ],
                dtype=int,
            )
            data = np.array([1.0, 1.0, 1.0, 1.0], dtype=float)
            inv = sparse.COO(coords=coords, data=data, shape=(1, 2, 2, 1))
            self.inventory = xr.DataArray(
                inv,
                dims=("activity", "flow", "year", "root activity"),
                coords={
                    "activity": [0],
                    "flow": [0, 1],
                    "year": [2000, 2001],
                    "root activity": [0],
                },
            )
            self.biosphere_indices = {
                "2000": {
                    0: {"name": "CO2", "compartment": "air", "subcompartment": ""},
                    1: {"name": "CH4", "compartment": "air", "subcompartment": ""},
                }
            }

    trails = DummyTrailsTwoFlowsTiny()
    call_state = {"count": 0}

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s", "s"],
                "region": ["World", "World"],
                "variable": ["CO2", "CH4"],
                "unit": ["Gt CO2/yr", "Mt CH4/yr"],
                "2000.5": [0.0, 0.0],
                "2001.5": [0.0, 0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {
            ("CO2", "air", ""): "CO2",
            ("CH4", "air", ""): "CH4",
        }, {}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1"]
        years = np.array([2000.5, 2001.5], dtype=float)
        species = ["CO2", "CH4"]
        forcing_vals = np.zeros((1, len(config_names), len(years), len(species)))
        if call_state["count"] > 0:
            forcing_vals[0, :, :, 0] = np.array([1.0e-15, 2.0e-15])[None, :]
            forcing_vals[0, :, :, 1] = np.array([2.0e-15, 1.0e-15])[None, :]
        call_state["count"] += 1

        forcing = xr.DataArray(
            forcing_vals,
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
                "specie": species,
            },
        )
        temperature = xr.DataArray(
            np.sum(forcing_vals, axis=3),
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
            },
        )
        return types.SimpleNamespace(
            forcing=forcing,
            temperature=temperature,
            configs=list(config_names),
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1"],
        per_species_runs=False,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    assert int(getattr(trails.delta_temperature.data, "nnz", 0)) > 0
    assert float(trails.delta_temperature.sum().item()) != 0.0


def test_fair_fast_mode_allocates_negative_emissions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsNegativeOnly:
        def __init__(self) -> None:
            self.debug = False
            coords = np.array(
                [
                    [0, 0],  # activity
                    [0, 0],  # flow
                    [0, 1],  # year
                    [0, 0],  # root activity
                ],
                dtype=int,
            )
            data = np.array([1.0, 1.0], dtype=float)
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
                        "name": "Carbon dioxide, in air",
                        "compartment": "natural resource",
                        "subcompartment": "in air",
                    }
                }
            }

    trails = DummyTrailsNegativeOnly()

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["CO2 AFOLU"],
                "unit": ["Gt CO2/yr"],
                "2000.5": [0.0],
                "2001.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        key = ("Carbon dioxide, in air", "natural resource", "in air")
        return {key: "CO2 AFOLU"}, {key: -1.0}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1"]
        local = fair_rf_module._normalize_emissions_columns(emissions_df)
        local = local[
            (local["scenario"] == scenario) & (local["region"].str.lower() == "world")
        ].copy()
        years = np.array([2000.5, 2001.5], dtype=float)
        rows = local[local["variable"] == "CO2 AFOLU"]
        vals = np.zeros(2, dtype=float)
        if not rows.empty:
            row = rows.iloc[0]
            vals = np.array(
                [
                    float(pd.to_numeric(row.get("2000.5", 0.0), errors="coerce") or 0),
                    float(pd.to_numeric(row.get("2001.5", 0.0), errors="coerce") or 0),
                ],
                dtype=float,
            )

        # Mimic FaIR returning aggregate CO2 forcing channel, not CO2 AFOLU.
        forcing = xr.DataArray(
            vals[None, None, :, None],
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
                "specie": ["CO2"],
            },
        )
        temperature = xr.DataArray(
            vals[None, None, :],
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
            },
        )
        return types.SimpleNamespace(
            forcing=forcing,
            temperature=temperature,
            configs=list(config_names),
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1"],
        per_species_runs=False,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    rf = np.asarray(trails.instant_radiative_forcing.data.todense(), dtype=float)
    temp = np.asarray(trails.delta_temperature.data.todense(), dtype=float)

    assert np.any(rf < 0.0)
    assert np.any(temp < 0.0)


def test_fair_fast_mode_splits_aggregate_co2_between_emission_and_uptake(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsCO2Split:
        def __init__(self) -> None:
            self.debug = False
            coords = np.array(
                [
                    [0, 0, 0, 0],  # activity
                    [0, 0, 1, 1],  # flow
                    [0, 1, 0, 1],  # year
                    [0, 0, 0, 0],  # root activity
                ],
                dtype=int,
            )
            # flow 0 (fossil) has +1, flow 1 (uptake) has +2 then sign override -1
            data = np.array([1.0, 1.0, 2.0, 2.0], dtype=float)
            inv = sparse.COO(coords=coords, data=data, shape=(1, 2, 2, 1))
            self.inventory = xr.DataArray(
                inv,
                dims=("activity", "flow", "year", "root activity"),
                coords={
                    "activity": [0],
                    "flow": [0, 1],
                    "year": [2000, 2001],
                    "root activity": [0],
                },
            )
            self.biosphere_indices = {
                "2000": {
                    0: {
                        "name": "Carbon dioxide, fossil",
                        "compartment": "air",
                        "subcompartment": "non-urban air or from high stacks",
                    },
                    1: {
                        "name": "Carbon dioxide, in air",
                        "compartment": "natural resource",
                        "subcompartment": "in air",
                    },
                }
            }

    trails = DummyTrailsCO2Split()
    call_state = {"count": 0}

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s", "s"],
                "region": ["World", "World"],
                "variable": ["CO2 FFI", "CO2 AFOLU"],
                "unit": ["Gt CO2/yr", "Gt CO2/yr"],
                "2000.5": [0.0, 0.0],
                "2001.5": [0.0, 0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        fossil_key = (
            "Carbon dioxide, fossil",
            "air",
            "non-urban air or from high stacks",
        )
        uptake_key = ("Carbon dioxide, in air", "natural resource", "in air")
        return {fossil_key: "CO2 FFI", uptake_key: "CO2 AFOLU"}, {uptake_key: -1.0}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1"]
        local = fair_rf_module._normalize_emissions_columns(emissions_df)
        local = local[
            (local["scenario"] == scenario) & (local["region"].str.lower() == "world")
        ].copy()

        years = np.array([2000.5, 2001.5], dtype=float)
        if call_state["count"] == 0:
            net = np.zeros(2, dtype=float)
        else:
            ffi = local[local["variable"] == "CO2 FFI"]
            afolu = local[local["variable"] == "CO2 AFOLU"]
            ffi_vals = np.zeros(2, dtype=float)
            afolu_vals = np.zeros(2, dtype=float)
            if not ffi.empty:
                row = ffi.iloc[0]
                ffi_vals = np.array(
                    [
                        float(
                            pd.to_numeric(row.get("2000.5", 0.0), errors="coerce") or 0
                        ),
                        float(
                            pd.to_numeric(row.get("2001.5", 0.0), errors="coerce") or 0
                        ),
                    ],
                    dtype=float,
                )
            if not afolu.empty:
                row = afolu.iloc[0]
                afolu_vals = np.array(
                    [
                        float(
                            pd.to_numeric(row.get("2000.5", 0.0), errors="coerce") or 0
                        ),
                        float(
                            pd.to_numeric(row.get("2001.5", 0.0), errors="coerce") or 0
                        ),
                    ],
                    dtype=float,
                )
            # FaIR may expose CO2 as one aggregate channel.
            net = ffi_vals + afolu_vals
        call_state["count"] += 1

        forcing = xr.DataArray(
            net[None, None, :, None],
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
                "specie": ["CO2"],
            },
        )
        temperature = xr.DataArray(
            net[None, None, :],
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
            },
        )
        return types.SimpleNamespace(
            forcing=forcing,
            temperature=temperature,
            configs=list(config_names),
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1"],
        per_species_runs=False,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    rf = np.asarray(trails.instant_radiative_forcing.data.todense(), dtype=float)
    temp = np.asarray(trails.delta_temperature.data.todense(), dtype=float)
    rf_by_flow = rf[0].sum(axis=(0, 2))
    temp_by_flow = temp[0].sum(axis=(0, 2))

    assert rf_by_flow[0] > 0.0
    assert temp_by_flow[0] > 0.0
    assert rf_by_flow[1] < 0.0
    assert temp_by_flow[1] < 0.0


def test_fair_flow_mapping_filters_non_atmospheric_categories(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsSulfurFilter:
        def __init__(self) -> None:
            self.debug = False
            coords = np.array(
                [
                    [0, 0],  # activity
                    [0, 1],  # flow
                    [0, 0],  # year
                    [0, 0],  # root activity
                ],
                dtype=int,
            )
            data = np.array([1.0, 100.0], dtype=float)
            inv = sparse.COO(coords=coords, data=data, shape=(1, 2, 1, 1))
            self.inventory = xr.DataArray(
                inv,
                dims=("activity", "flow", "year", "root activity"),
                coords={
                    "activity": [0],
                    "flow": [0, 1],
                    "year": [2000],
                    "root activity": [0],
                },
            )
            self.biosphere_indices = {
                "2000": {
                    0: {
                        "name": "Sulfur",
                        "compartment": "air",
                        "subcompartment": "non-urban air or from high stacks",
                    },
                    1: {
                        "name": "Sulfur",
                        "compartment": "natural resource",
                        "subcompartment": "in ground",
                    },
                }
            }

    trails = DummyTrailsSulfurFilter()
    call_state = {"count": 0}

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["Sulfur"],
                "unit": ["Mt Sulfur/yr"],
                "2000.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {"Sulfur": "Sulfur"}, {}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1"]
        local = fair_rf_module._normalize_emissions_columns(emissions_df)
        local = local[
            (local["scenario"] == scenario) & (local["region"].str.lower() == "world")
        ].copy()
        years = np.array([2000.5], dtype=float)
        vals = np.zeros(1, dtype=float)
        if call_state["count"] > 0:
            rows = local[local["variable"] == "Sulfur"]
            if not rows.empty:
                row = rows.iloc[0]
                vals = np.array(
                    [float(pd.to_numeric(row.get("2000.5", 0.0), errors="coerce") or 0)],
                    dtype=float,
                )
        call_state["count"] += 1

        forcing = xr.DataArray(
            vals[None, None, :, None],
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
                "specie": ["Sulfur"],
            },
        )
        temperature = xr.DataArray(
            vals[None, None, :],
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": years,
            },
        )
        return types.SimpleNamespace(
            forcing=forcing,
            temperature=temperature,
            configs=list(config_names),
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    run_fair_delta_rf(
        trails,
        scenario="s",
        config_names=["c1"],
        per_species_runs=False,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    rf = np.asarray(trails.instant_radiative_forcing.data.todense(), dtype=float)
    rf_by_flow = rf[0].sum(axis=(0, 2))

    assert rf_by_flow[0] > 0.0
    assert rf_by_flow[1] == 0.0


def test_inventory_emissions_prefers_name_mapping_over_category_mapping() -> None:
    inv_years = [2000, 2001]
    # dims: (flow, year, root activity)
    inv_data = sparse.COO(
        coords=np.array(
            [
                [0, 0, 1, 1],  # flow
                [0, 1, 0, 1],  # year
                [0, 0, 0, 0],  # root activity
            ],
            dtype=int,
        ),
        data=np.array([1.0, 2.0, 3.0, 4.0], dtype=float),
        shape=(2, 2, 1),
    )
    flow_pos_to_key = {
        0: ("Methane, fossil", "air", "urban air close to ground"),
        1: ("Methane, fossil", "air", "non-urban air or from high stacks"),
    }
    species_map: dict[object, str] = {
        ("Methane, fossil", "air", "urban air close to ground"): "SF6",
        ("Methane, fossil", "air", "non-urban air or from high stacks"): "N2O",
        "Methane, fossil": "CH4",
    }
    signs: dict[object, float] = {}

    out = fair_rf_module._inventory_emissions_by_fair_species(
        inv_data=inv_data,
        inv_years=inv_years,
        n_flow=2,
        flow_pos_to_key=flow_pos_to_key,
        species_map=species_map,
        signs=signs,
    )

    assert out.columns.tolist() == ["CH4"]
    np.testing.assert_allclose(
        out["CH4"].to_numpy(dtype=float),
        np.array([4.0, 6.0], dtype=float),
    )


def test_fair_precursor_response_species_are_captured(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class DummyTrailsPrecursor:
        def __init__(self) -> None:
            self.debug = False
            data = np.array([1.0e9], dtype=float)
            coords = np.array([[0], [0], [1], [0]], dtype=int)
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
                        "name": "Nitrogen oxides",
                        "compartment": "air",
                        "subcompartment": "",
                    }
                }
            }

    class DummyFairRunPrecursor:
        def __init__(
            self,
            scenario: str,
            config_names: list[str],
            species: list[str],
            ozone_response: float,
        ) -> None:
            timebounds = np.array([2000.5, 2001.5], dtype=float)
            forcing_vals = np.zeros(
                (1, len(config_names), len(timebounds), len(species)),
                dtype=float,
            )
            if "Ozone" in species:
                forcing_vals[:, :, 1, species.index("Ozone")] = float(ozone_response)
            forcing = xr.DataArray(
                forcing_vals,
                dims=("scenario", "config", "timebounds", "specie"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": timebounds,
                    "specie": species,
                },
            )
            temperature_vals = np.zeros((1, len(config_names), len(timebounds)))
            temperature_vals[:, :, 1] = float(ozone_response)
            temperature = xr.DataArray(
                temperature_vals,
                dims=("scenario", "config", "timebounds"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": timebounds,
                },
            )
            self.forcing = forcing
            self.temperature = temperature

    trails = DummyTrailsPrecursor()

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["NOx"],
                "unit": ["Mt NO2/yr"],
                "2000.5": [0.0],
                "2001.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {("Nitrogen oxides", "air", ""): "NOx"}, {}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> DummyFairRunPrecursor:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1", "c2"]
        local = emissions_df[
            (emissions_df["scenario"] == scenario)
            & (emissions_df["region"].str.lower() == "world")
        ].copy()
        species = sorted(local["variable"].astype(str).unique().tolist())
        nox_rows = local[local["variable"] == "NOx"]
        if nox_rows.empty:
            ozone_response = 0.0
        else:
            ozone_response = float(
                pd.to_numeric(nox_rows["2001.5"], errors="coerce").fillna(0.0).iloc[0]
            )
        return DummyFairRunPrecursor(
            str(scenario),
            list(config_names),
            species,
            ozone_response,
        )

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
        quantiles=[50.0],
    )

    arr = np.asarray(rf.data.todense(), dtype=float)
    assert float(np.max(np.abs(arr))) > 0.0


def test_fair_all_mapped_species_return_non_null_rf(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    species_map, signs = load_species_mapping(DEFAULT_MAPPING_YAML)
    mapping_items = list(species_map.items())
    assert mapping_items

    n_flow = len(mapping_items)
    years = np.array([2000, 2001], dtype=int)
    flow_ids = np.arange(n_flow, dtype=int)
    data = np.full(n_flow, 1.0e9, dtype=float)
    year_idx = np.ones(n_flow, dtype=int)

    biosphere_meta: dict[int, dict[str, str]] = {}
    for i, (flow_key, _specie) in enumerate(mapping_items):
        if isinstance(flow_key, tuple):
            name = str(flow_key[0]) if len(flow_key) >= 1 else ""
            compartment = str(flow_key[1]) if len(flow_key) >= 2 else "air"
            subcompartment = str(flow_key[2]) if len(flow_key) >= 3 else ""
            name_key: object = name
        else:
            name = str(flow_key)
            compartment = "air"
            subcompartment = ""
            name_key = name

        sign_value = float(signs.get(flow_key, signs.get(name_key, 1.0)))
        year_idx[i] = 0 if sign_value < 0 else 1
        biosphere_meta[int(i)] = {
            "name": name,
            "compartment": compartment,
            "subcompartment": subcompartment,
        }

    coords = np.vstack(
        [
            np.zeros(n_flow, dtype=int),  # activity
            flow_ids,  # flow
            year_idx,  # year
            np.zeros(n_flow, dtype=int),  # root activity
        ]
    )
    inv = sparse.COO(coords=coords, data=data, shape=(1, n_flow, len(years), 1))

    class DummyTrailsAllMapped:
        def __init__(self) -> None:
            self.debug = False
            self.inventory = xr.DataArray(
                inv,
                dims=("activity", "flow", "year", "root activity"),
                coords={
                    "activity": [0],
                    "flow": flow_ids,
                    "year": years,
                    "root activity": [0],
                },
            )
            self.biosphere_indices = {"2000": biosphere_meta}

    trails = DummyTrailsAllMapped()

    mapped_species = sorted({str(v) for v in species_map.values()})
    base_rows = [
        {
            "scenario": "s",
            "region": "World",
            "variable": specie,
            "unit": f"Mt {specie}/yr",
            "2000": 0.0,
            "2001": 0.0,
        }
        for specie in mapped_species
    ]
    base_df = pd.DataFrame(base_rows)
    base_lookup = {
        str(row["variable"]): np.array([float(row["2000"]), float(row["2001"])])
        for _, row in base_df.iterrows()
    }

    response_species = {
        rsp
        for values in fair_rf_module._PRECURSOR_RESPONSE_SPECIES.values()
        for rsp in values
    }
    forcing_species = sorted(
        {
            ("CO2" if specie in {"CO2 FFI", "CO2 AFOLU"} else str(specie))
            for specie in mapped_species
        }
        | response_species
    )
    forcing_index = {specie: i for i, specie in enumerate(forcing_species)}

    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return base_df.copy()

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = kwargs.get("scenario") or args[1]
        config_names = kwargs.get("config_names") or ["c1", "c2"]
        local = fair_rf_module._normalize_emissions_columns(emissions_df)
        local = local[
            (local["scenario"] == scenario) & (local["region"].str.lower() == "world")
        ].copy()

        forcing_vals = np.zeros(
            (1, len(config_names), len(years), len(forcing_species)),
            dtype=float,
        )
        for _, row in local.iterrows():
            variable = str(row["variable"])
            if variable not in base_lookup:
                continue
            row_vals = np.array(
                [
                    float(pd.to_numeric(row.get("2000", 0.0), errors="coerce") or 0.0),
                    float(pd.to_numeric(row.get("2001", 0.0), errors="coerce") or 0.0),
                ],
                dtype=float,
            )
            delta = row_vals - base_lookup[variable]
            if not np.any(delta):
                continue

            primary = "CO2" if variable in {"CO2 FFI", "CO2 AFOLU"} else variable
            forcing_vals[:, :, :, forcing_index[primary]] += delta[None, None, :]
            for response in fair_rf_module._PRECURSOR_RESPONSE_SPECIES.get(primary, ()):
                forcing_vals[:, :, :, forcing_index[response]] += (
                    0.5 * delta[None, None, :]
                )

        forcing = xr.DataArray(
            forcing_vals,
            dims=("scenario", "config", "timebounds", "specie"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": np.array([2000.0, 2001.0], dtype=float),
                "specie": forcing_species,
            },
        )
        temperature = xr.DataArray(
            np.sum(forcing_vals, axis=3),
            dims=("scenario", "config", "timebounds"),
            coords={
                "scenario": [str(scenario)],
                "config": list(config_names),
                "timebounds": np.array([2000.0, 2001.0], dtype=float),
            },
        )
        return types.SimpleNamespace(forcing=forcing, temperature=temperature)

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    rf = run_fair_delta_rf(
        trails,
        scenario="s",
        mapping_yaml=DEFAULT_MAPPING_YAML,
        config_names=["c1", "c2"],
        per_species_runs=True,
        per_species_workers=1,
        validate_emissions_delta=False,
        scale_factor=1.0,
        quantiles=[50.0],
    )

    arr = np.asarray(rf.data.todense(), dtype=float)
    assert arr.shape == (1, len(years), n_flow, 1)

    for flow_pos, (flow_key, _specie) in enumerate(mapping_items):
        series = arr[0, :, flow_pos, 0]
        assert np.all(
            np.isfinite(series)
        ), f"non-finite RF for mapped flow {flow_key!r}"
        assert np.any(np.abs(series) > 0), f"null RF for mapped flow {flow_key!r}"


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


def test_sanitize_emissions_year_values_handles_string_dtype_column() -> None:
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
    df["2002.5"] = df["2002.5"].astype("string")

    out = _sanitize_emissions_year_values(df)

    assert float(out.loc[0, "2000.5"]) == pytest.approx(1.0)
    assert float(out.loc[0, "2001.5"]) == pytest.approx(0.0)
    assert float(out.loc[0, "2002.5"]) == pytest.approx(3.5)
    assert pd.api.types.is_numeric_dtype(out["2002.5"])


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


def test_plot_temporal_scores_auto_trims_year_window() -> None:
    trails = DummyTrailsScores()
    fig = plot_temporal_scores(
        trails=trails,
        year_range=None,
        show_flow_contributions=False,
        show_cumulative_axis=False,
        legend_top_n=1,
    )

    assert fig is not None
    assert len(fig.data) >= 1
    x = list(fig.data[0].x)
    assert x[0] == 2003
    assert x[-1] == 2004


def test_plot_temporal_scores_stacked_separates_sign_groups() -> None:
    class DummyTrailsMixedSigns:
        def __init__(self) -> None:
            years = np.array([2000, 2001, 2002], dtype=int)
            vals = np.array(
                [
                    [1.0, -1.0, 1.0],
                    [-2.0, 2.0, -2.0],
                ],
                dtype=float,
            )
            self.scores = xr.DataArray(
                vals,
                dims=("activity", "year"),
                coords={"activity": [0, 1], "year": years},
            )
            self.characterized_inventory = None
            self.activity_indices = {
                "2000": {
                    0: {
                        "name": "A0",
                        "reference product": "P0",
                        "location": "GLO",
                    },
                    1: {
                        "name": "A1",
                        "reference product": "P1",
                        "location": "GLO",
                    },
                }
            }

    trails = DummyTrailsMixedSigns()
    fig = plot_temporal_scores(
        trails=trails,
        stacked=True,
        year_range=(2000, 2002),
        show_flow_contributions=False,
        show_cumulative_axis=False,
        legend_top_n=2,
    )

    assert fig is not None
    stackgroups = {getattr(trace, "stackgroup", None) for trace in fig.data}
    assert "positive" in stackgroups
    assert "negative" in stackgroups
    assert "one" not in stackgroups

    for trace in fig.data:
        y = np.asarray(trace.y, dtype=float)
        stackgroup = getattr(trace, "stackgroup", None)
        if stackgroup == "positive":
            assert np.all(y >= 0.0)
        if stackgroup == "negative":
            assert np.all(y <= 0.0)


def test_cache_dir_uses_short_hash(example_package: Any) -> None:
    cache_dir = cache_dir_for_package(
        example_package, value_dtype="float32", index_dtype="int32"
    )
    name = cache_dir.name
    assert name.startswith("interp_")
    suffix = name.split("interp_", 1)[1]
    assert len(suffix) == 12
