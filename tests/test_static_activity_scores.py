from __future__ import annotations

from pathlib import Path

from datapackage import Package
import numpy as np
import pytest

from trails.lca import lca
from trails.lcia import get_lcia_method_names
from trails.static_activity_scores import (
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
