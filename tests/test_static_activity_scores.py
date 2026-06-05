from __future__ import annotations

from pathlib import Path

from datapackage import Package
import numpy as np
import pytest

from trails.lca import lca
from trails.lcia import get_lcia_method_names
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


def _load_example_trails() -> Trails:
    return Trails(Package(str(EXAMPLE_PACKAGE)), interpolate_annual=False)


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
    adaptive.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=5,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[GWP_METHOD],
        adaptive_score_cutoff=0.001,
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
        if data.get("adaptive_cutoff_reason") == "adaptive_score_cutoff"
    ]

    assert adaptive._routing_params["adaptive_enabled"] is True
    assert adaptive._routing_params["adaptive_score_cutoff"] == pytest.approx(0.001)
    assert adaptive._routing_params["adaptive_relative_score_cutoff"] is None
    assert adaptive.graph.number_of_nodes() == 15
    assert adaptive.graph.number_of_edges() == 17
    assert adaptive.graph.number_of_nodes() < fixed.graph.number_of_nodes()
    assert pruned_nodes
    assert all(
        float(data["adaptive_cutoff_potential"]) <= 0.001 for data in pruned_nodes
    )
    assert float(adaptive.scores.sum()) == pytest.approx(
        0.15497121214866638,
        rel=1e-7,
    )


def test_adaptive_routing_can_run_without_depth_cap() -> None:
    """Adaptive routing can let score potential define the stopping depth."""
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
    adaptive.temporal_routing(
        start_year=2005,
        start_act_idx=2,
        amount=1.0,
        max_depth=None,
        min_amount=0.0,
        show_progress=False,
        attribute_to_roots=False,
        adaptive_methods=[GWP_METHOD],
        adaptive_score_cutoff=0.1,
        adaptive_use_cache=False,
    )

    assert adaptive._routing_params["max_depth"] is None
    assert adaptive.graph.number_of_nodes() == 6
    assert adaptive.graph.number_of_edges() == 5
    assert {
        data.get("adaptive_cutoff_reason")
        for _node, data in adaptive.graph.nodes(data=True)
    } == {None, "adaptive_score_cutoff"}


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
            adaptive_score_cutoff=0.1,
            adaptive_use_cache=False,
        )
