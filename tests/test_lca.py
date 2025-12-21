import importlib
import numpy as np
import sparse

lca_module = importlib.import_module("trails.lca")


class DummyLCA:
    def __init__(self, demand, data_objs):
        self.demand = demand
        self.data_objs = data_objs
        self.dicts = type(
            "Dicts",
            (),
            {"product": {0: 0, 1: 1}, "biosphere": {0: 0, 1: 1}},
        )()
        self.inventory = np.array([[1.0], [2.0]])
        self.supply_array = np.array([float(demand.get(0, 0.0)), float(demand.get(1, 0.0))])

    def lci(self):
        return None


class DummyTrails:
    def __init__(self):
        self.scenario_labels = ["2005"]
        self.scenario_index = {"2005": 0}
        self.value_dtype = np.float64
        self.A = sparse.COO(coords=[[0], [0], [0]], data=[1.0], shape=(1, 2, 2))
        self.B = sparse.COO(
            coords=[np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)],
            data=np.array([], dtype=float),
            shape=(1, 2, 2),
        )
        self.activity_indices = {"2005": {0: {"name": "act", "reference product": "prod", "unit": "u", "location": "GLO"}}}
        self.biosphere_indices = {"2005": {0: {"name": "flow", "compartment": "air", "subcompartment": "low", "unit": "kg"}}}

    def temporal_traversal(self, **kwargs):
        return {(2005, 0): 1.0}, {}

    def frontier_to_demand_vectors(self, frontier):
        vec = np.zeros(2, dtype=float)
        for (year, act), amt in frontier.items():
            vec[act] += amt
        return {2005: vec}

    def expand_temporal_exchanges(self, **kwargs):
        return {}

    def _map_year_to_scenario_year(self, year):
        return 2005


def test_nearest_metadata_label_for_year():
    trails = DummyTrails()
    trails.activity_indices = {"2000": {}, "2010": {}}
    assert lca_module._nearest_metadata_label_for_year(trails, 2003) == "2000"
    assert lca_module._nearest_metadata_label_for_year(trails, 2009) == "2010"


def test_build_datapackage_for_year_from_trails(example_trails):
    dp, tech_idx, bio_idx, uncertain = lca_module.build_datapackage_for_year_from_trails(
        example_trails, year=2005
    )
    assert dp is not None
    assert any("battery electric vehicle, production" in key for key in tech_idx)
    assert any("Carbon dioxide, fossil" in key for key in bio_idx)
    assert uncertain == []


def test_lca_static_mode(monkeypatch):
    trails = DummyTrails()

    def fake_build_dp(*args, **kwargs):
        return object(), {}, {}, []

    def fake_fill_characterization_factors_matrices(*args, **kwargs):
        return np.ones((1, 2))

    monkeypatch.setattr(lca_module, "build_datapackage_for_year_from_trails", fake_build_dp)
    monkeypatch.setattr(lca_module.bc, "LCA", DummyLCA)
    monkeypatch.setattr(
        lca_module, "fill_characterization_factors_matrices", fake_fill_characterization_factors_matrices
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


def test_compute_node_impact_intensities(monkeypatch):
    trails = DummyTrails()

    def fake_build_dp(*args, **kwargs):
        return object(), {}, {}, []

    def fake_fill_characterization_factors_matrices(*args, **kwargs):
        return np.ones((1, 2))

    monkeypatch.setattr(lca_module, "build_datapackage_for_year_from_trails", fake_build_dp)
    monkeypatch.setattr(lca_module.bc, "LCA", DummyLCA)
    monkeypatch.setattr(
        lca_module, "fill_characterization_factors_matrices", fake_fill_characterization_factors_matrices
    )
    lca_module.use_temporal_distributions = False

    result = lca_module.compute_node_impact_intensities(
        trails=trails, nodes=[(2005, 0)], methods=["dummy"]
    )
    assert result[(2005, 0)] == 3.0


def test_lca_example_package_gasoline_car(example_trails):
    method = (
        "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - "
        "global warming potential (GWP100)"
    )
    result = lca_module.lca(
        trails=example_trails,
        start_year=2005,
        start_act_idx=13,
        methods=[method],
        amount=1.0,
        max_depth=2,
        min_amount=0.0,
        show_progress=False,
        return_provenance=False,
        use_temporal_distributions=True,
    )

    impact = result["results_by_impact_year"]
    assert 2005 in impact
    assert 2020 in impact

    assert np.isclose(impact[2005]["scores"], 0.2491568943951279, rtol=1e-6)
    assert np.isclose(impact[2020]["scores"], 0.012761818594299257, rtol=1e-6)
    assert np.isclose(impact[2005]["scores_by_first_level_child"][13], 0.13472940237261355, rtol=1e-6)
    assert np.isclose(impact[2020]["scores_by_first_level_child"][13], 0.008235294255428016, rtol=1e-6)
