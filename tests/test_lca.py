import importlib
import numpy as np
import sparse
import pytest

from trails.trails import Trails

lca_module = importlib.import_module("trails.lca")


class DummyLCA:
    def __init__(self, demand: dict[int, float], data_objs: list[object]) -> None:
        """Initialize a minimal LCA stub for tests.

        :param demand: Demand mapping used by the dummy LCA.
        :type demand: dict[int, float]
        :param data_objs: Data objects passed to the LCA.
        :type data_objs: list[object]
        """
        self.demand = demand
        self.data_objs = data_objs
        self.dicts = type(
            "Dicts",
            (),
            {"product": {0: 0, 1: 1}, "biosphere": {0: 0, 1: 1}},
        )()
        self.inventory = np.array([[1.0], [2.0]])
        self.supply_array = np.array(
            [float(demand.get(0, 0.0)), float(demand.get(1, 0.0))]
        )

    def lci(self) -> None:
        """No-op LCI stub for tests.

        :returns: None.
        :rtype: None
        """
        return None


class DummyTrails:
    def __init__(self) -> None:
        """Initialize a minimal Trails stub for tests.

        :returns: None.
        :rtype: None
        """
        self.scenario_labels = ["2005"]
        self.scenario_index = {"2005": 0}
        self.value_dtype = np.float64
        self.A = sparse.COO(coords=[[0], [0], [0]], data=[1.0], shape=(1, 2, 2))
        self.B = sparse.COO(
            coords=[
                np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=int),
            ],
            data=np.array([], dtype=float),
            shape=(1, 2, 2),
        )
        self.activity_indices = {
            "2005": {
                0: {
                    "name": "act",
                    "reference product": "prod",
                    "unit": "u",
                    "location": "GLO",
                }
            }
        }
        self.biosphere_indices = {
            "2005": {
                0: {
                    "name": "flow",
                    "compartment": "air",
                    "subcompartment": "low",
                    "unit": "kg",
                }
            }
        }

    def temporal_traversal(self, **kwargs) -> tuple[dict[tuple[int, int], float], dict]:
        """Return a fixed frontier for traversal tests.

        :returns: Frontier mapping and empty direct biosphere mapping.
        :rtype: tuple[dict[tuple[int, int], float], dict]
        """
        return {(2005, 0): 1.0}, {}

    def frontier_to_demand_vectors(
        self, frontier: dict[tuple[int, int], float]
    ) -> dict[int, np.ndarray]:
        """Build a simple demand vector from a frontier mapping.

        :param frontier: Frontier mapping to convert.
        :type frontier: dict[tuple[int, int], float]
        :returns: Mapping of year to demand vectors.
        :rtype: dict[int, numpy.ndarray]
        """
        vec = np.zeros(2, dtype=float)
        for (year, act), amt in frontier.items():
            vec[act] += amt
        return {2005: vec}

    def expand_temporal_exchanges(self, **kwargs) -> dict:
        """Return empty temporal exchanges for tests.

        :returns: Empty demand mapping.
        :rtype: dict
        """
        return {}

    def _map_year_to_scenario_year(self, year: int) -> int:
        """Map any year to the stub scenario year.

        :param year: Calendar year to map.
        :type year: int
        :returns: Scenario year.
        :rtype: int
        """
        return 2005


def test_nearest_metadata_label_for_year() -> None:
    """Verify nearest metadata label selection.

    :returns: None.
    :rtype: None
    """
    trails = DummyTrails()
    trails.activity_indices = {"2000": {}, "2010": {}}
    assert lca_module._nearest_metadata_label_for_year(trails, 2003) == "2000"
    assert lca_module._nearest_metadata_label_for_year(trails, 2009) == "2010"


def test_build_datapackage_for_year_from_trails(example_trails: Trails) -> None:
    """Verify datapackage construction for a year.

    :param example_trails: Trails fixture for datapackage construction.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    dp, tech_idx, bio_idx, uncertain = (
        lca_module.build_datapackage_for_year_from_trails(example_trails, year=2005)
    )
    assert dp is not None
    assert any("battery electric vehicle, production" in key for key in tech_idx)
    assert any("Carbon dioxide, fossil" in key for key in bio_idx)
    assert uncertain == []


def test_lca_static_mode(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify static LCA behavior without temporal distributions.

    :param monkeypatch: Pytest monkeypatch fixture.
    :type monkeypatch: pytest.MonkeyPatch
    :returns: None.
    :rtype: None
    """
    trails = DummyTrails()

    def fake_build_dp(*args, **kwargs) -> tuple[object, dict, dict, list]:
        """Return a minimal datapackage tuple for tests.

        :returns: Datapackage tuple with empty metadata.
        :rtype: tuple[object, dict, dict, list]
        """
        return object(), {}, {}, []

    def fake_fill_characterization_factors_matrices(*args, **kwargs) -> np.ndarray:
        """Return a dummy characterization matrix for tests.

        :returns: Dummy characterization matrix.
        :rtype: numpy.ndarray
        """
        return np.ones((1, 2))

    monkeypatch.setattr(
        lca_module, "build_datapackage_for_year_from_trails", fake_build_dp
    )
    monkeypatch.setattr(lca_module.bc, "LCA", DummyLCA)
    monkeypatch.setattr(
        lca_module,
        "fill_characterization_factors_matrices",
        fake_fill_characterization_factors_matrices,
    )

    result = lca_module.lca(
        trails=trails,
        start_year=2005,
        start_act_idx=0,
        methods=["dummy"],
        amount=1.0,
        max_depth=1,
        min_amount=0.0,
        show_progress=False,
        return_provenance=False,
        use_temporal_distributions=False,
    )
    impact = result["results_by_impact_year"]
    assert 2005 in impact
    assert impact[2005]["scores"] == 3.0
