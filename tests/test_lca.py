import importlib
import numpy as np
import sparse

from trails.bw_interface import _nearest_metadata_label_for_year
from trails.trails import Trails

lca_module = importlib.import_module("trails.lca")


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
    assert _nearest_metadata_label_for_year(trails, 2003) == "2000"
    assert _nearest_metadata_label_for_year(trails, 2009) == "2010"


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


def test_lca_static_mode(example_trails: Trails) -> None:
    """Verify static LCA behavior without temporal distributions.

    :param example_trails: Trails fixture initialized from example data package datapackage.
    :type example_trails: trails.trails.Trails
    :returns: None.
    :rtype: None
    """
    activity_indices = next(iter(example_trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))

    example_trails.temporal_routing(
        start_year=2005,
        start_act_idx=start_act_idx,
        max_depth=1,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        debug=False,
    )

    lca_module.lca(
        trails=example_trails,
        show_progress=False,
        compute_score=False,
        store_inventory=False,
        attribute_to_roots=False,
    )
    assert example_trails.scores is None
