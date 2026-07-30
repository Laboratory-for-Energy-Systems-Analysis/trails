from __future__ import annotations

import pytest

from trails.search import search_activity


class DummySearchTrails:
    activity_indices = {
        "2030": {
            0: {
                "name": "market for electricity, high voltage",
                "reference product": "electricity, high voltage",
                "location": "DE",
            },
            1: {
                "name": "market for electricity, high voltage",
                "reference product": "electricity, high voltage",
                "location": "GLO",
            },
            2: {
                "name": "electricity production, wind",
                "reference product": "electricity, high voltage",
                "location": "DE",
            },
        },
        # Repeated metadata across scenario years must remain deduplicated by index.
        "2040": {
            0: {
                "name": "market for electricity, high voltage",
                "reference product": "electricity, high voltage",
                "location": "DE",
            }
        },
    }
    biosphere_indices = {
        "2030": {
            0: {
                "name": "Carbon dioxide, fossil",
                "compartment": "air",
                "subcompartment": "urban air close to ground",
            }
        }
    }


def _indices(table: object) -> list[int]:
    return [int(row[0]) for row in table.rows]


def test_search_activity_by_location_only() -> None:
    result = search_activity(DummySearchTrails(), location="de", match="exact")

    assert _indices(result) == [0, 2]


def test_search_activity_combines_identity_filters() -> None:
    result = search_activity(
        DummySearchTrails(),
        name="market for electricity",
        reference_product="high voltage",
        location="DE",
    )

    assert _indices(result) == [0]


def test_search_activity_combines_reference_product_and_location() -> None:
    result = search_activity(
        DummySearchTrails(),
        reference_product="electricity, high voltage",
        location="GLO",
        match="exact",
    )

    assert _indices(result) == [1]


def test_search_activity_rejects_location_for_biosphere() -> None:
    with pytest.raises(ValueError, match="only valid for technosphere"):
        search_activity(DummySearchTrails(), location="DE", kind="biosphere")
