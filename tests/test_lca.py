import importlib
import numpy as np
import sparse

lca_module = importlib.import_module("trails.lca")


class DummyLCA:
    def __init__(self, demand, data_objs):
        """Initialize a minimal LCA stub for tests."""
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

    def lci(self):
        """No-op LCI stub for tests."""
        return None


class DummyTrails:
    def __init__(self):
        """Initialize a minimal Trails stub for tests."""
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

    def temporal_traversal(self, **kwargs):
        """Return a fixed frontier for traversal tests."""
        return {(2005, 0): 1.0}, {}

    def frontier_to_demand_vectors(self, frontier):
        """Build a simple demand vector from a frontier mapping."""
        vec = np.zeros(2, dtype=float)
        for (year, act), amt in frontier.items():
            vec[act] += amt
        return {2005: vec}

    def expand_temporal_exchanges(self, **kwargs):
        """Return empty temporal exchanges for tests."""
        return {}

    def _map_year_to_scenario_year(self, year):
        """Map any year to the stub scenario year."""
        return 2005


def test_nearest_metadata_label_for_year():
    """Verify nearest metadata label selection."""
    trails = DummyTrails()
    trails.activity_indices = {"2000": {}, "2010": {}}
    assert lca_module._nearest_metadata_label_for_year(trails, 2003) == "2000"
    assert lca_module._nearest_metadata_label_for_year(trails, 2009) == "2010"


def test_build_datapackage_for_year_from_trails(example_trails):
    """Verify datapackage construction for a year."""
    dp, tech_idx, bio_idx, uncertain = (
        lca_module.build_datapackage_for_year_from_trails(example_trails, year=2005)
    )
    assert dp is not None
    assert any("battery electric vehicle, production" in key for key in tech_idx)
    assert any("Carbon dioxide, fossil" in key for key in bio_idx)
    assert uncertain == []


def test_lca_static_mode(monkeypatch):
    """Verify static LCA behavior without temporal distributions."""
    trails = DummyTrails()

    def fake_build_dp(*args, **kwargs):
        """Return a minimal datapackage tuple for tests."""
        return object(), {}, {}, []

    def fake_fill_characterization_factors_matrices(*args, **kwargs):
        """Return a dummy characterization matrix for tests."""
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
