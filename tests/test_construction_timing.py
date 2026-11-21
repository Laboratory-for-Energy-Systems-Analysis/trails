"""A shared supplier keeps fallback timing on routes without a build year."""

import csv

from datapackage import Package
import pytest

from trails import Trails
from trails.datapackage import BASE_COLS_B


def make_construction_package(tmp_path, snapshots, *, timing="dated"):
    """Build a route with a shared market, an untimed link and genuine lifecycle events."""
    activities = (
        "operating plant",
        "other buyer",
        "market for asset",
        "asset wrapper",
        "asset construction",
        "maintenance",
        "retirement",
        "independent asset",
    )
    cols = (
        "index of activity",
        "index of product",
        "value",
        "flip",
        "temporal_distribution",
        "temporal_offsets",
        "temporal_weights",
        "construction_timing",
    )
    common_rows = [(i, i, 1, 0, "", "", "", "") for i in range(len(activities))]
    common_rows += [
        (1, 2, 1, 1, "", "", "", ""),
        (2, 3, 1, 1, "", "", "", "pass_through"),
        (3, 4, 1, 1, 6, "[-10]", "[1]", "pass_through"),
        (3, 5, 1, 1, 6, "[2]", "[1]", ""),
        (3, 6, 1, 1, 6, "[10]", "[1]", ""),
        (3, 7, 1, 1, 6, "[-1]", "[1]", "dated"),
    ]
    resources = []
    for year, (amount, offsets, weights) in snapshots.items():
        path = f"inventories/remind/test/{year}"
        dest = tmp_path / path
        dest.mkdir(parents=True)
        with (dest / "A_matrix.csv").open("w", newline="") as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerow(cols)
            writer.writerows(common_rows)
            writer.writerow((0, 2, amount, 1, 6, str(offsets), str(weights), timing))
        with (dest / "A_matrix_index.csv").open("w", newline="") as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerow(("name", "reference product", "unit", "location", "index"))
            writer.writerows(
                (name, name, "unit", "GLO", i) for i, name in enumerate(activities)
            )
        with (dest / "B_matrix.csv").open("w", newline="") as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerow(BASE_COLS_B)
            bio = {"index of activity": 4, "value": 1, "loc": 1}
            writer.writerow(bio.get(field, 0) for field in BASE_COLS_B)
        resources.extend(
            {
                "name": f"{name.lower()}_{year}",
                "path": f"{path}/{name}.csv",
                "format": "csv",
            }
            for name in ("A_matrix", "A_matrix_index", "B_matrix")
        )
    return Trails(
        Package(
            {"name": "construction-test", "resources": resources},
            base_path=str(tmp_path),
        ),
        cache_interpolation=False,
        interpolation_start_year_offset=-30,
        interpolation_end_year_offset=15,
    )


@pytest.fixture
def construction_package(tmp_path):
    return make_construction_package(tmp_path, {2050: (2, [-15, -5], [0.25, 0.75])})


@pytest.mark.parametrize(
    "start,expected",
    [
        (
            0,
            {
                4: {2035: 0.5, 2045: 1.5},
                5: {2037: 0.5, 2047: 1.5},
                6: {2045: 0.5, 2055: 1.5},
                7: {2034: 0.5, 2044: 1.5},
            },
        ),
        (1, {4: {2040: 1}, 5: {2052: 1}, 6: {2060: 1}, 7: {2049: 1}}),
    ],
)
def test_construction_date_only_passes_to_same_asset(
    construction_package, start, expected
):
    model = construction_package
    model.temporal_routing(
        start_year=2050, start_act_idx=start, max_depth=3, show_progress=False
    )
    routed = {i: {} for i in expected}
    for _, dst, data in model.graph.edges(data=True):
        if dst[-1] in routed:
            routed[dst[-1]][dst[0]] = data["amount"]
    assert routed == expected
    frontier, _ = model.temporal_traversal(2050, start, max_depth=3)
    assert frontier == {
        (y, i): amount for i, dates in expected.items() for y, amount in dates.items()
    }


def test_interpolated_inventory_preserves_a_fixed_build_year(tmp_path):
    # Snapshot resources can appear in either order in the data package.
    model = make_construction_package(
        tmp_path, {2050: (4, [-40], [1]), 2030: (2, [-20], [1])}
    )
    for year, amount in ((2035, 2.5), (2040, 3), (2030, 2), (2050, 4)):
        model.temporal_routing(
            start_year=year, start_act_idx=0, max_depth=3, show_progress=False
        )
        # Both the market and its supplier receive the full amount in 2010.
        for activity in (2, 4):
            dates = {
                dst[0]: data["amount"]
                for _, dst, data in model.graph.edges(data=True)
                if dst[-1] == activity
            }
            assert dates == pytest.approx({2010: amount})
        frontier, _ = model.temporal_traversal(year, 0, max_depth=3)
        assert {(date, i) for date, i in frontier if i == 4} == {(2010, 4)}
        assert frontier[2010, 4] == pytest.approx(amount)


def test_interpolation_blends_build_years_and_excludes_future_construction(tmp_path):
    model = make_construction_package(
        tmp_path,
        {
            2030: (2, [-20, 0], [0.6, 0.4]),
            2050: (4, [-40, -10, 0], [0.2, 0.4, 0.4]),
        },
    )
    expected = [
        (2040, {2010: 0.5, 2030: 0.25, 2040: 0.25}, 3),
        (2035, {2010: 0.625, 2030: 0.375}, 2.5),
        (2045, {2010: 3 / 7, 2030: 1 / 7, 2040: 3 / 7}, 3.5),
    ]
    for year, weights, amount in expected:
        distribution = model.get_temporal_distribution(year, 0, 2)
        assert {
            year + offset: share
            for offset, share in distribution.iter_offsets_and_weights()
        } == pytest.approx(weights)
        model.temporal_routing(
            start_year=year, start_act_idx=0, max_depth=1, show_progress=False
        )
        assert {
            dst[0]: data["amount"]
            for _, dst, data in model.graph.edges(data=True)
            if dst[-1] == 2
        } == pytest.approx({date: amount * weight for date, weight in weights.items()})


def test_ordinary_empirical_timing_keeps_its_relative_offsets(tmp_path):
    model = make_construction_package(
        tmp_path, {2030: (2, [-20], [1]), 2050: (4, [-40], [1])}, timing=""
    )
    model.temporal_routing(
        start_year=2035, start_act_idx=0, max_depth=1, show_progress=False
    )
    assert {
        dst[0]: data["amount"]
        for _, dst, data in model.graph.edges(data=True)
        if dst[-1] == 2
    } == pytest.approx({2015: 2.5})
