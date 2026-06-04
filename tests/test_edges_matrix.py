import importlib

import numpy as np
import pytest
from scipy import sparse as sp
import sparse
import xarray as xr

from trails.edges_matrix import (
    _build_edges_characterization_matrices_for_year,
    score_inventory_with_edges,
)

lca_module = importlib.import_module("trails.lca")


class DummyTrails:
    def __init__(self, inventory: xr.DataArray) -> None:
        self.inventory = inventory
        self.value_dtype = np.float64
        self.scores = None

    def _map_year_to_scenario_year(self, year: int) -> int:
        return 2000


class DummyMatrixTrails:
    def __init__(self, n_activities: int = 1, n_flows: int = 1) -> None:
        self.A = np.zeros((1, n_activities, n_activities))
        self.B = np.zeros((1, n_activities, n_flows))
        self.B[0, 0, 0] = 1.0
        self.activity_indices = {}
        self.biosphere_indices = {}
        self.package = None

    def _get_scenario_context(self, year: int):
        return year, str(year), 0


class DummyScoringTrails(DummyTrails):
    def __init__(self, inventory: xr.DataArray) -> None:
        super().__init__(inventory)
        self.A = np.zeros((2, 2, 2))
        self.B = np.zeros((2, 2, 2))
        self.package = None
        self.activity_indices = {
            str(year): {
                0: {"name": "act0", "reference product": "prod0", "location": "CH"},
                1: {"name": "act1", "reference product": "prod1", "location": "FR"},
            }
            for year in (2000, 2001)
        }
        self.biosphere_indices = {
            str(year): {
                0: {"name": "flow0", "unit": "kg", "location": "CH"},
                1: {"name": "flow1", "unit": "kg", "location": "FR"},
            }
            for year in (2000, 2001)
        }

    def _map_year_to_scenario_year(self, year: int) -> int:
        return int(year)

    def _get_scenario_context(self, year: int):
        return year, str(year), int(year) - 2000


def test_build_edges_matrices_suppresses_edges_stdout(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    class FakeEdgeLCIA:
        def __init__(self, *args, **kwargs):
            self.raw_cfs_data = [{"supplier": {"matrix": "biosphere"}}]
            self.characterization_matrix = None

        def _preprocess_lookups(self, **kwargs):
            print("lookup output")

        def apply_strategies(self, strategies):
            print("strategy output")
            self.cfs_mapping.append(
                {
                    "supplier": {"matrix": "biosphere"},
                    "consumer": {"matrix": "technosphere"},
                    "positions": ((0, 0),),
                    "direction": "biosphere-technosphere",
                    "value": 4.0,
                }
            )

    monkeypatch.setattr(
        "trails.edges_matrix._ensure_edges_available",
        lambda: FakeEdgeLCIA,
    )

    matrices = _build_edges_characterization_matrices_for_year(
        DummyMatrixTrails(),
        ["edge-method"],
        year=2000,
    )

    assert np.asarray(matrices[0].todense()).tolist() == [[4.0]]
    assert capsys.readouterr().out == ""


def test_build_edges_matrices_restricts_edges_and_metadata(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen = {}

    class FakeEdgeLCIA:
        def __init__(self, *args, **kwargs):
            self.lca = kwargs["lca"]
            self.raw_cfs_data = [{"supplier": {"matrix": "biosphere"}}]
            self.characterization_matrix = None

        def _preprocess_lookups(self, **kwargs):
            seen["edges"] = set(self.biosphere_edges)
            seen["biosphere_positions"] = [
                flow["position"] for flow in self.biosphere_flows
            ]
            seen["activity_positions"] = [
                flow["position"] for flow in self.technosphere_flows
            ]
            seen["inventory_shape"] = self.lca.inventory.shape

        def apply_strategies(self, strategies):
            return None

        def evaluate_cfs(self, scenario_idx):
            self.characterization_matrix = sp.csr_matrix(self.lca.inventory.shape)

    monkeypatch.setattr(
        "trails.edges_matrix._ensure_edges_available",
        lambda: FakeEdgeLCIA,
    )

    _build_edges_characterization_matrices_for_year(
        DummyMatrixTrails(n_activities=3, n_flows=4),
        ["edge-method"],
        year=2000,
        biosphere_edges={(1, 2), (3, 0)},
    )

    assert seen["edges"] == {(1, 2), (3, 0)}
    assert seen["biosphere_positions"] == [1, 3]
    assert seen["activity_positions"] == [0, 2]
    assert seen["inventory_shape"] == (4, 3)


def test_score_inventory_with_edges_reuses_cached_mappings_for_next_year(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array(
                [
                    [0, 0, 1],
                    [0, 0, 1],
                    [0, 1, 1],
                ],
                dtype=np.int64,
            ),
            data=np.array([5.0, 6.0, 7.0]),
            shape=(2, 2, 2),
        ),
        dims=("activity", "flow", "year"),
        coords={"activity": [0, 1], "flow": [0, 1], "year": [2000, 2001]},
    )
    trails = DummyScoringTrails(inventory)
    mapped_edge_sets = []

    class FakeEdgeLCIA:
        def __init__(self, *args, **kwargs):
            self.lca = kwargs["lca"]
            self.raw_cfs_data = [{"supplier": {"matrix": "biosphere"}}]
            self.cfs_mapping = []
            self.characterization_matrix = None

        def _preprocess_lookups(self, **kwargs):
            return None

        def apply_strategies(self, strategies):
            mapped_edge_sets.append(set(self.biosphere_edges))
            for edge in sorted(self.biosphere_edges):
                self.cfs_mapping.append(
                    {
                        "supplier": {"matrix": "biosphere"},
                        "consumer": {"matrix": "technosphere"},
                        "positions": (edge,),
                        "direction": "biosphere-technosphere",
                        "value": 1.0,
                    }
                )

        def evaluate_cfs(self, scenario_idx):
            matrix = sp.lil_matrix(self.lca.inventory.shape)
            value = float(int(scenario_idx) - 1999)
            for cf in self.cfs_mapping:
                for supplier, consumer in cf["positions"]:
                    matrix[int(supplier), int(consumer)] = value
            self.characterization_matrix = matrix.tocsr()

        def _evaluate_cf_numeric_value(
            self,
            raw_value,
            *,
            resolved_params,
            scenario_idx,
        ):
            return float(raw_value) * float(int(scenario_idx) - 1999)

    monkeypatch.setattr(
        "trails.edges_matrix._ensure_edges_available",
        lambda: FakeEdgeLCIA,
    )

    scores = score_inventory_with_edges(trails, ["edge-method"], show_progress=False)

    assert mapped_edge_sets == [{(0, 0)}, {(1, 1)}]
    assert np.asarray(scores.data.todense()).tolist() == [
        [5.0, 12.0],
        [0.0, 14.0],
    ]


def test_score_inventory_with_edges_can_disable_cross_year_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array(
                [
                    [0, 0, 1],
                    [0, 0, 1],
                    [0, 1, 1],
                ],
                dtype=np.int64,
            ),
            data=np.array([5.0, 6.0, 7.0]),
            shape=(2, 2, 2),
        ),
        dims=("activity", "flow", "year"),
        coords={"activity": [0, 1], "flow": [0, 1], "year": [2000, 2001]},
    )
    trails = DummyScoringTrails(inventory)
    mapped_edge_sets = []

    class FakeEdgeLCIA:
        def __init__(self, *args, **kwargs):
            self.lca = kwargs["lca"]
            self.raw_cfs_data = [{"supplier": {"matrix": "biosphere"}}]
            self.cfs_mapping = []

        def _preprocess_lookups(self, **kwargs):
            return None

        def apply_strategies(self, strategies):
            mapped_edge_sets.append(set(self.biosphere_edges))
            for edge in sorted(self.biosphere_edges):
                self.cfs_mapping.append(
                    {
                        "supplier": {"matrix": "biosphere"},
                        "consumer": {"matrix": "technosphere"},
                        "positions": (edge,),
                        "direction": "biosphere-technosphere",
                        "value": 1.0,
                    }
                )

        def _evaluate_cf_numeric_value(
            self,
            raw_value,
            *,
            resolved_params,
            scenario_idx,
        ):
            return float(raw_value) * float(int(scenario_idx) - 1999)

    monkeypatch.setattr(
        "trails.edges_matrix._ensure_edges_available",
        lambda: FakeEdgeLCIA,
    )

    scores = score_inventory_with_edges(
        trails,
        ["edge-method"],
        reuse_cached_cfs=False,
        show_progress=False,
    )

    assert mapped_edge_sets == [{(0, 0)}, {(0, 0), (1, 1)}]
    assert np.asarray(scores.data.todense()).tolist() == [
        [5.0, 12.0],
        [0.0, 14.0],
    ]


def test_score_inventory_with_edges_uses_activity_specific_cfs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The same biosphere flow can receive different CFs by consuming activity."""
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array(
                [
                    [0, 1, 1],
                    [1, 1, 0],
                    [0, 0, 1],
                ],
                dtype=np.int64,
            ),
            data=np.array([10.0, 10.0, 2.0]),
            shape=(2, 2, 2),
        ),
        dims=("activity", "flow", "year"),
        coords={"activity": [0, 1], "flow": [0, 1], "year": [2000, 2001]},
    )
    trails = DummyTrails(inventory)
    calls = []

    def fake_builder(trails_arg, methods, *, year, **kwargs):
        calls.append((trails_arg, methods, year, kwargs))
        return [
            sp.csr_matrix(
                np.array(
                    [
                        [0.0, 11.0],
                        [3.0, 7.0],
                    ]
                )
            )
        ]

    monkeypatch.setattr(
        "trails.edges_matrix._build_edges_characterization_matrices_for_year",
        fake_builder,
    )

    scores = score_inventory_with_edges(trails, ["edge-method"], show_progress=False)

    assert scores.dims == ("activity", "year")
    assert np.asarray(scores.data.todense()).tolist() == [
        [30.0, 0.0],
        [70.0, 22.0],
    ]
    assert trails.scores is scores
    assert len(calls) == 1
    assert calls[0][1] == ["edge-method"]
    assert calls[0][2] == 2000
    assert calls[0][3]["biosphere_edges"] == {(1, 0), (1, 1), (0, 1)}


def test_score_inventory_with_edges_keeps_method_and_root_dimensions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array(
                [
                    [0, 0, 1],
                    [1, 1, 0],
                    [0, 0, 0],
                    [0, 1, 1],
                ],
                dtype=np.int64,
            ),
            data=np.array([5.0, 2.0, 3.0]),
            shape=(2, 2, 1, 2),
        ),
        dims=("activity", "flow", "year", "root activity"),
        coords={
            "activity": [0, 1],
            "flow": [0, 1],
            "year": [2000],
            "root activity": [0, 1],
        },
    )
    trails = DummyTrails(inventory)

    def fake_builder(*args, **kwargs):
        return [
            sp.csr_matrix(np.array([[0.0, 3.0], [2.0, 0.0]])),
            sp.csr_matrix(np.array([[0.0, 0.5], [10.0, 0.0]])),
        ]

    monkeypatch.setattr(
        "trails.edges_matrix._build_edges_characterization_matrices_for_year",
        fake_builder,
    )

    scores = score_inventory_with_edges(trails, ["m1", "m2"], show_progress=False)
    dense = np.asarray(scores.data.todense())

    assert scores.dims == ("method", "activity", "year", "root activity")
    assert scores.coords["method"].values.tolist() == ["m1", "m2"]
    assert dense[0, 0, 0, 0] == 10.0
    assert dense[0, 0, 0, 1] == 4.0
    assert dense[0, 1, 0, 1] == 9.0
    assert dense[1, 0, 0, 0] == 50.0
    assert dense[1, 0, 0, 1] == 20.0
    assert dense[1, 1, 0, 1] == 1.5


def test_score_inventory_with_edges_reports_matrix_build_progress(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inventory = xr.DataArray(
        sparse.COO(
            coords=np.array([[0], [0], [0]], dtype=np.int64),
            data=np.array([5.0]),
            shape=(1, 1, 1),
        ),
        dims=("activity", "flow", "year"),
        coords={"activity": [0], "flow": [0], "year": [2000]},
    )
    trails = DummyTrails(inventory)
    progress_instances = []

    class DummyProgress:
        def __init__(self, *, total, desc, unit):
            self.total = total
            self.desc = desc
            self.unit = unit
            self.n = 0
            progress_instances.append(self)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return None

        def set_postfix_str(self, value):
            self.postfix = value

        def update(self, amount):
            self.n += amount

    def fake_builder(trails_arg, methods, *, progress, **kwargs):
        progress.set_postfix_str("year=2000, method=1/2")
        progress.update(len(methods))
        return [
            sp.csr_matrix(np.array([[2.0]])),
            sp.csr_matrix(np.array([[3.0]])),
        ]

    monkeypatch.setattr("trails.edges_matrix.tqdm", DummyProgress)
    monkeypatch.setattr(
        "trails.edges_matrix._build_edges_characterization_matrices_for_year",
        fake_builder,
    )

    score_inventory_with_edges(trails, ["m1", "m2"], show_progress=True)

    assert len(progress_instances) == 1
    assert progress_instances[0].total == 2
    assert progress_instances[0].desc == "EDGES LCIA"
    assert progress_instances[0].unit == "method-year"
    assert progress_instances[0].n == 2


def test_lca_rejects_regular_and_edges_methods() -> None:
    with pytest.raises(ValueError, match="mutually exclusive"):
        lca_module.lca(
            object(),
            methods=["regular-method"],
            edges_methods=["edge-method"],
        )


def test_lca_defaults_prefer_edges_methods_for_final_scoring() -> None:
    class DummyTrails:
        default_methods = ["regular-method"]
        default_edges_methods = ["edge-method"]
        default_ei_version = "3.11"

    methods, edges_methods, ei_version = lca_module._resolve_lca_method_defaults(
        DummyTrails(),
        methods=None,
        edges_methods=None,
        ei_version=None,
    )
    assert methods is None
    assert edges_methods == ["edge-method"]
    assert ei_version == "3.11"

    methods, edges_methods, _ei_version = lca_module._resolve_lca_method_defaults(
        DummyTrails(),
        methods=["explicit-regular"],
        edges_methods=None,
        ei_version=None,
    )
    assert methods == ["explicit-regular"]
    assert edges_methods is None
