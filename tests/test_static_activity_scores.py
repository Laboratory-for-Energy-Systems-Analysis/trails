from __future__ import annotations

from pathlib import Path

from datapackage import Package
import numpy as np
import pytest

from trails.lca import lca
from trails.lcia import get_lcia_method_names
import trails.trails as trails_module
from trails.static_activity_scores import (
    StaticActivityScores,
    _activity_score_potential,
    _activity_scores_from_product_intensities,
    _ensure_static_activity_scores,
)
import trails.static_activity_scores as static_activity_scores
from trails.trails import Trails

EXAMPLE_PACKAGE = Path("examples/example data package/datapackage.json")
GWP_METHOD = (
    "CML v4.8 2016 no LT - climate change no LT - global warming potential "
    "(GWP100) no LT"
)
SYNTHETIC_METHOD = "synthetic adaptive method"
ADAPTIVE_CUTOFF_REASON = "adaptive_relative_score_cutoff"


def _load_example_trails() -> Trails:
    return Trails(Package(str(EXAMPLE_PACKAGE)), interpolate_annual=False)


class _SyntheticRoutingTrails(Trails):
    """Minimal Trails object for deterministic temporal-routing tests."""

    def __init__(
        self,
        expansions: dict[tuple[int, int], dict[tuple[int, int], float]],
        production_amounts: dict[int, float] | None = None,
    ) -> None:
        self.debug = False
        self.default_methods = None
        self.default_edges_methods = None
        self.default_ei_version = "3.11"
        self.activity_indices = {}
        self.scenario_index = {"2000": 0}
        self.expansions = {
            (int(year), int(act)): {
                (int(child_year), int(child_act)): float(amount)
                for (child_year, child_act), amount in children.items()
            }
            for (year, act), children in expansions.items()
        }
        self.production_amounts = {
            int(activity): float(amount)
            for activity, amount in (production_amounts or {}).items()
        }
        self.expand_calls: list[tuple[int, int, float]] = []

    def _map_year_to_scenario_year(self, year: int) -> int:
        return int(year)

    def _get_scenario_context(self, year: int) -> tuple[int, str, int]:
        return int(year), str(int(year)), 0

    def _activity_amount_from_product_demand(
        self,
        t: int,
        act_idx: int,
        amount: float,
    ) -> float:
        return float(amount) / self._production_amount(int(t), int(act_idx))

    def _production_amount(self, t: int, act_idx: int) -> float:
        return float(self.production_amounts.get(int(act_idx), 1.0))

    def _has_direct_biosphere(
        self,
        scenario_year: int,
        act_idx: int,
        cache: dict[tuple[int, int], bool],
    ) -> bool:
        return False

    def _expand_temporal_child_demands_fast(
        self,
        *,
        year: int,
        act_idx: int,
        amount: float,
        use_temporal_distributions: bool = True,
        debug: bool = False,
    ) -> dict[tuple[int, int], float]:
        self.expand_calls.append((int(year), int(act_idx), float(amount)))
        return dict(self.expansions.get((int(year), int(act_idx)), {}))


def _install_synthetic_static_scores(
    monkeypatch: pytest.MonkeyPatch,
    unit_scores: dict[int, float],
) -> StaticActivityScores:
    max_activity = max(unit_scores, default=0)
    scores = np.zeros((1, 1, max_activity + 1), dtype=np.float64)
    for activity, score in unit_scores.items():
        scores[0, 0, int(activity)] = float(score)

    result = StaticActivityScores(
        methods=(SYNTHETIC_METHOD,),
        years=np.asarray([2000], dtype=np.int64),
        scores=scores,
        year_index={2000: 0},
        method_index={SYNTHETIC_METHOD: 0},
    )

    monkeypatch.setattr(
        trails_module,
        "_ensure_static_activity_scores",
        lambda *args, **kwargs: result,
    )
    return result


def _graph_nodes_by_depth_and_activity(trails: Trails) -> dict[tuple[int, int], dict]:
    return {
        (int(data["depth"]), int(data["act_idx"])): data
        for _node, data in trails.graph.nodes(data=True)
    }


def _relative_cutoff_for_effective_score(
    trails: Trails,
    *,
    effective_cutoff: float,
    start_year: int,
    start_act_idx: int,
    amount: float,
    methods: list[str],
) -> float:
    scores = _ensure_static_activity_scores(
        trails,
        methods=methods,
        years=[int(start_year)],
        use_cache=False,
    )
    root_potential, _ = _activity_score_potential(
        scores,
        year=int(start_year),
        activity=int(start_act_idx),
        amount=float(amount),
    )
    return float(effective_cutoff) / float(root_potential)


def test_static_activity_scores_match_static_lca() -> None:
    """Adjoint activity intensities should match a conventional static score."""
    trails = _load_example_trails()
    scores = _ensure_static_activity_scores(
        trails,
        methods=[GWP_METHOD],
        years=[2005],
        use_cache=False,
    )

    trails.static_lca(2005, 2, [GWP_METHOD])

    assert scores.scores.shape == (1, 1, 17)
    assert float(scores.scores[0, 0, 2]) == pytest.approx(
        float(trails.static_score[0]),
        rel=1e-12,
    )


def test_static_activity_scores_cache_roundtrip_and_matrix_invalidation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The cache should reload matching scores and miss changed matrices."""
    monkeypatch.setattr(
        static_activity_scores,
        "cache_dir_for_package",
        lambda *args, **kwargs: tmp_path,
    )

    trails = _load_example_trails()
    first = _ensure_static_activity_scores(
        trails,
        methods=[GWP_METHOD],
        years=[2005],
        use_cache=True,
    )
    assert first.cache_path is not None
    assert first.cache_path.exists()
    assert not first.loaded_from_cache

    trails._static_activity_score_cache.clear()
    second = _ensure_static_activity_scores(
        trails,
        methods=[GWP_METHOD],
        years=[2005],
        use_cache=True,
    )
    assert second.loaded_from_cache
    np.testing.assert_allclose(second.scores, first.scores)

    trails.B = trails.B * 2.0
    trails._static_activity_score_cache.clear()
    trails._static_activity_score_fingerprint = None
    changed = _ensure_static_activity_scores(
        trails,
        methods=[GWP_METHOD],
        years=[2005],
        use_cache=True,
    )

    assert changed.cache_path != first.cache_path
    assert not changed.loaded_from_cache
    assert float(changed.scores[0, 0, 2]) == pytest.approx(
        2.0 * float(first.scores[0, 0, 2]),
        rel=1e-12,
    )


def test_activity_scores_respect_negative_production_sign() -> None:
    """Activity scores should follow the sign convention used for production."""
    intensities = np.array([[10.0], [20.0]])
    product_indices = np.array([0, 1], dtype=np.int64)
    production_values = np.array([1.0, -1.0], dtype=np.float64)

    scores = _activity_scores_from_product_intensities(
        intensities,
        product_indices,
        production_values,
    )

    assert scores[0, 0] == pytest.approx(10.0)
    assert scores[0, 1] == pytest.approx(-20.0)


def test_activity_score_potential_reuses_finite_absolute_arrays() -> None:
    """Potential lookup should sanitize scores once and reuse derived arrays."""
    scores = StaticActivityScores(
        methods=("m1", "m2"),
        years=np.asarray([2030], dtype=np.int64),
        scores=np.asarray(
            [
                [[np.nan, -2.0, np.inf]],
                [[3.0, -np.inf, -4.0]],
            ],
            dtype=np.float64,
        ),
        year_index={2030: 0},
        method_index={"m1": 0, "m2": 1},
    )

    potential, by_method = _activity_score_potential(
        scores,
        year=2030,
        activity=0,
        amount=-2.0,
    )
    assert potential == pytest.approx(6.0)
    assert by_method == {"m1": pytest.approx(0.0), "m2": pytest.approx(6.0)}

    abs_scores = scores._abs_scores
    max_abs_scores = scores._max_abs_scores
    assert abs_scores is not None
    assert max_abs_scores is not None

    potential, by_method = _activity_score_potential(
        scores,
        year=2030,
        activity=1,
        amount=5.0,
    )
    assert potential == pytest.approx(10.0)
    assert by_method == {"m1": pytest.approx(10.0), "m2": pytest.approx(0.0)}
    assert scores._abs_scores is abs_scores
    assert scores._max_abs_scores is max_abs_scores

    potential, by_method = _activity_score_potential(
        scores,
        year=2030,
        activity=2,
        amount=0.5,
    )
    assert potential == pytest.approx(2.0)
    assert by_method == {"m1": pytest.approx(0.0), "m2": pytest.approx(2.0)}


def test_adaptive_routing_expands_only_branches_above_relative_cutoff(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Adaptive routing should not expand children once potential is at cutoff."""
    routing = _SyntheticRoutingTrails(
        {
            (2000, 0): {
                (2000, 1): 0.02,
                (2000, 2): 0.02,
                (2000, 7): 0.01,
            },
            (2000, 1): {(2000, 3): 0.004, (2000, 4): 0.004},
            (2000, 2): {(2000, 5): 100.0},
            (2000, 3): {},
            (2000, 4): {(2000, 6): 100.0},
            (2000, 7): {(2000, 8): 100.0},
        }
    )
    _install_synthetic_static_scores(
        monkeypatch,
        {
            0: 100.0,
            1: 100.0,
            2: 10.0,
            3: 300.0,
            4: 100.0,
            5: 10_000.0,
            6: 10_000.0,
            7: 100.0,
            8: 10_000.0,
        },
    )

    routing.temporal_routing(
        start_year=2000,
        start_act_idx=0,
        amount=1.0,
        max_depth=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[SYNTHETIC_METHOD],
        adaptive_relative_score_cutoff=0.01,
        adaptive_use_cache=False,
    )

    nodes = _graph_nodes_by_depth_and_activity(routing)

    assert [(year, act) for year, act, _amount in routing.expand_calls] == [
        (2000, 0),
        (2000, 1),
        (2000, 3),
    ]
    assert set(nodes) == {(0, 0), (1, 1), (1, 2), (1, 7), (2, 3), (2, 4)}
    assert routing.graph.number_of_edges() == 5

    assert nodes[(1, 1)]["score_potential"] == pytest.approx(2.0)
    assert nodes[(1, 1)]["adaptive_cutoff_reason"] is None
    assert nodes[(2, 3)]["score_potential"] == pytest.approx(1.2)
    assert nodes[(2, 3)]["frontier_reasons"] == {"leaf": pytest.approx(0.004)}

    assert nodes[(1, 2)]["score_potential"] == pytest.approx(0.2)
    assert nodes[(1, 2)]["adaptive_cutoff_reason"] == ADAPTIVE_CUTOFF_REASON
    assert nodes[(1, 2)]["adaptive_cutoff_potential"] == pytest.approx(0.2)
    assert nodes[(1, 2)]["frontier_reasons"] == {
        ADAPTIVE_CUTOFF_REASON: pytest.approx(0.02)
    }

    assert nodes[(2, 4)]["score_potential"] == pytest.approx(0.4)
    assert nodes[(2, 4)]["adaptive_cutoff_reason"] == ADAPTIVE_CUTOFF_REASON
    assert nodes[(2, 4)]["adaptive_cutoff_potential"] == pytest.approx(0.4)

    assert nodes[(1, 7)]["score_potential"] == pytest.approx(1.0)
    assert nodes[(1, 7)]["adaptive_cutoff_reason"] == ADAPTIVE_CUTOFF_REASON
    assert nodes[(1, 7)]["adaptive_cutoff_potential"] == pytest.approx(1.0)


def test_adaptive_relative_cutoff_scales_with_root_score_potential(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A relative cutoff should be multiplied by the functional unit potential."""
    routing = _SyntheticRoutingTrails(
        {
            (2000, 0): {(2000, 1): 1.0, (2000, 2): 1.0},
            (2000, 2): {},
        }
    )
    _install_synthetic_static_scores(
        monkeypatch,
        {
            0: 200.0,
            1: 0.15,
            2: 0.25,
        },
    )

    routing.temporal_routing(
        start_year=2000,
        start_act_idx=0,
        amount=1.0,
        max_depth=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[SYNTHETIC_METHOD],
        adaptive_relative_score_cutoff=1e-3,
        adaptive_use_cache=False,
    )

    nodes = _graph_nodes_by_depth_and_activity(routing)

    assert routing._routing_params["adaptive_root_score_potential"] == pytest.approx(
        200.0
    )
    assert routing._routing_params["adaptive_effective_score_cutoff"] == pytest.approx(
        0.2
    )
    assert [(year, act) for year, act, _amount in routing.expand_calls] == [
        (2000, 0),
        (2000, 2),
    ]
    assert nodes[(1, 1)]["score_potential"] == pytest.approx(0.15)
    assert nodes[(1, 1)]["adaptive_cutoff_reason"] == ADAPTIVE_CUTOFF_REASON
    assert nodes[(1, 2)]["score_potential"] == pytest.approx(0.25)
    assert nodes[(1, 2)]["adaptive_cutoff_reason"] is None


def test_adaptive_relative_cutoff_uses_functional_unit_product_score(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Non-unit production exchanges must not shrink the relative threshold."""
    routing = _SyntheticRoutingTrails(
        {
            (2000, 0): {(2000, 1): 0.1, (2000, 2): 0.1},
            (2000, 2): {},
        },
        production_amounts={0: 5.0, 1: 10.0, 2: 10.0},
    )
    _install_synthetic_static_scores(
        monkeypatch,
        {
            0: 100.0,
            1: 0.5,
            2: 1.5,
        },
    )

    routing.temporal_routing(
        start_year=2000,
        start_act_idx=0,
        amount=10.0,
        max_depth=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[SYNTHETIC_METHOD],
        adaptive_relative_score_cutoff=1e-3,
        adaptive_use_cache=False,
    )

    nodes = _graph_nodes_by_depth_and_activity(routing)

    assert routing._routing_params["adaptive_root_score_potential"] == pytest.approx(
        1000.0
    )
    assert routing._routing_params["adaptive_effective_score_cutoff"] == pytest.approx(
        1.0
    )
    assert [(year, act) for year, act, _amount in routing.expand_calls] == [
        (2000, 0),
        (2000, 2),
    ]
    assert nodes[(1, 1)]["score_potential"] == pytest.approx(0.5)
    assert nodes[(1, 1)]["adaptive_cutoff_reason"] == ADAPTIVE_CUTOFF_REASON
    assert nodes[(1, 2)]["score_potential"] == pytest.approx(1.5)
    assert nodes[(1, 2)]["adaptive_cutoff_reason"] is None


def test_fast_temporal_child_expansion_matches_public_expansion() -> None:
    """Routing's flat expansion helper should match the public dict API."""
    trails = _load_example_trails()

    for use_temporal_distributions in (False, True):
        public = trails.expand_temporal_exchanges(
            year=2005,
            act_idx=2,
            amount=1.0,
            use_temporal_distributions=use_temporal_distributions,
            debug=False,
        )
        expected = {
            (int(year), int(act)): float(amount)
            for year, mapping in public.items()
            for act, amount in mapping.items()
        }
        fast = trails._expand_temporal_child_demands_fast(
            year=2005,
            act_idx=2,
            amount=1.0,
            use_temporal_distributions=use_temporal_distributions,
            debug=False,
        )

        assert set(fast) == set(expected)
        for key, value in expected.items():
            assert fast[key] == pytest.approx(value, rel=1e-12, abs=1e-15)


def test_production_amount_vector_matches_sparse_diagonal() -> None:
    """Production amount cache should match finite nonzero diagonal entries."""
    trails = _load_example_trails()
    t = trails.scenario_index["2005"]
    vector = trails._production_amount_vector(t)

    assert trails._production_amount_vector(t) is vector
    assert vector.shape == (trails.A.shape[1],)

    for act_idx in range(min(8, int(trails.A.shape[1]))):
        raw = float(trails.A[int(t), act_idx, act_idx])
        expected = 1.0 if (not np.isfinite(raw) or abs(raw) < 1e-30) else abs(raw)

        assert vector[act_idx] == pytest.approx(expected)
        assert trails._production_amount(t, act_idx) == pytest.approx(expected)


def test_regular_depth_assessments_remain_available() -> None:
    """The example case should still run with the fixed-depth routing mode."""
    assert GWP_METHOD in get_lcia_method_names()
    expected = {
        1: (6, 5, 0.15496498346328735),
        2: (13, 14, 0.15498079359531403),
        3: (23, 26, 0.15505138039588928),
        4: (28, 33, 0.15513959527015686),
        5: (28, 33, 0.15513959527015686),
    }

    for depth, (nodes, edges, score) in expected.items():
        trails = _load_example_trails()
        trails.temporal_routing(
            start_year=2005,
            start_act_idx=2,
            amount=1.0,
            max_depth=depth,
            min_amount=0.0,
            show_progress=False,
            attribute_to_roots=False,
        )
        lca(
            trails,
            methods=[GWP_METHOD],
            show_progress=False,
            compute_score=True,
            store_inventory=False,
            attribute_to_roots=False,
            solver_mode="direct",
        )

        assert trails.graph.number_of_nodes() == nodes
        assert trails.graph.number_of_edges() == edges
        assert float(trails.scores.sum()) == pytest.approx(score, rel=1e-7)


def test_adaptive_routing_prunes_low_potential_branches() -> None:
    """Adaptive routing should coexist with max-depth routing and mark pruned nodes."""
    effective_cutoff = 0.001
    fixed = _load_example_trails()
    fixed.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=5,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
    )

    adaptive = _load_example_trails()
    relative_cutoff = _relative_cutoff_for_effective_score(
        adaptive,
        effective_cutoff=effective_cutoff,
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        methods=[GWP_METHOD],
    )
    adaptive.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=5,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[GWP_METHOD],
        adaptive_relative_score_cutoff=relative_cutoff,
        adaptive_use_cache=False,
    )
    lca(
        adaptive,
        methods=[GWP_METHOD],
        show_progress=False,
        compute_score=True,
        store_inventory=False,
        attribute_to_roots=False,
        solver_mode="direct",
    )

    pruned_nodes = [
        data
        for _node, data in adaptive.graph.nodes(data=True)
        if data.get("adaptive_cutoff_reason") == ADAPTIVE_CUTOFF_REASON
    ]

    assert adaptive._routing_params["adaptive_enabled"] is True
    assert adaptive._routing_params["adaptive_relative_score_cutoff"] == pytest.approx(
        relative_cutoff
    )
    assert adaptive._routing_params["adaptive_effective_score_cutoff"] == pytest.approx(
        effective_cutoff
    )
    assert adaptive.graph.number_of_nodes() == 15
    assert adaptive.graph.number_of_edges() == 17
    assert adaptive.graph.number_of_nodes() < fixed.graph.number_of_nodes()
    assert pruned_nodes
    assert all(
        float(data["adaptive_cutoff_potential"]) <= effective_cutoff
        for data in pruned_nodes
    )
    assert float(adaptive.scores.sum()) == pytest.approx(
        0.15497121214866638,
        rel=1e-7,
    )


def test_adaptive_routing_can_run_without_depth_cap() -> None:
    """Adaptive routing can let score potential define the stopping depth."""
    effective_cutoff = 0.1
    fixed_mode = _load_example_trails()
    with pytest.raises(ValueError, match="Trails\\(\\.\\.\\., methods=\\.\\.\\.\\)"):
        fixed_mode.temporal_routing(
            start_year=2005,
            start_act_idx=2,
            amount=1.0,
            min_amount=0.0,
            show_progress=False,
            attribute_to_roots=False,
        )

    adaptive = _load_example_trails()
    relative_cutoff = _relative_cutoff_for_effective_score(
        adaptive,
        effective_cutoff=effective_cutoff,
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        methods=[GWP_METHOD],
    )
    adaptive.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[GWP_METHOD],
        adaptive_relative_score_cutoff=relative_cutoff,
        adaptive_use_cache=False,
    )

    assert adaptive._routing_params["max_depth"] is None
    assert adaptive.graph.number_of_nodes() == 6
    assert adaptive.graph.number_of_edges() == 5
    assert {
        data.get("adaptive_cutoff_reason")
        for _node, data in adaptive.graph.nodes(data=True)
    } == {None, ADAPTIVE_CUTOFF_REASON}


def test_adaptive_routing_uses_constructor_method_defaults() -> None:
    """Default routing should use constructor methods and a relative cutoff."""
    adaptive = Trails(
        Package(str(EXAMPLE_PACKAGE)),
        interpolate_annual=False,
        methods=[GWP_METHOD],
        ei_version="3.11",
    )
    adaptive.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_use_cache=False,
    )

    assert adaptive._routing_params["adaptive_methods"] == [GWP_METHOD]
    assert adaptive._routing_params["adaptive_ei_version"] == "3.11"
    assert adaptive._routing_params["max_depth"] is None
    assert adaptive._routing_params["adaptive_relative_score_cutoff"] == pytest.approx(
        1e-4
    )
    assert adaptive._routing_params["adaptive_enabled"] is True
    assert adaptive.graph.number_of_nodes() == 25
    assert adaptive.graph.number_of_edges() == 29


def test_explicit_max_depth_uses_fixed_depth_by_default() -> None:
    """Passing max_depth should remain the explicit fixed-depth routing mode."""
    fixed = Trails(
        Package(str(EXAMPLE_PACKAGE)),
        interpolate_annual=False,
        methods=[GWP_METHOD],
        ei_version="3.11",
    )
    fixed.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=1,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
    )

    assert fixed._routing_params["adaptive_enabled"] is False
    assert fixed._routing_params["adaptive_relative_score_cutoff"] is None
    assert fixed._routing_params["max_depth"] == 1
    assert fixed.graph.number_of_nodes() == 6


def test_adaptive_routing_rejects_edges_only_defaults() -> None:
    """EDGES methods cannot currently provide adaptive screening scores."""
    trails = Trails(
        Package(str(EXAMPLE_PACKAGE)),
        interpolate_annual=False,
        edges_methods=["edge-method"],
    )

    with pytest.raises(ValueError, match="requires regular LCIA methods"):
        trails.temporal_routing(
            start_year=2005,
            start_act_idx=2,
            amount=1.0,
            max_depth=None,
            min_amount=0.0,
            show_progress=False,
            attribute_to_roots=False,
            adaptive_relative_score_cutoff=0.1,
            adaptive_use_cache=False,
        )
