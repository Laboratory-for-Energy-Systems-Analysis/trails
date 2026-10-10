"""Calendar and signed LCI regression using the actual premise exporter."""

from pathlib import Path

import pytest

from dev.stock_vintage.check_lifecycle_roundtrip import validate


@pytest.mark.parametrize(
    "solver_mode,attribute_to_roots,interpolate_annual",
    [
        ("direct", True, True),
        ("direct", False, True),
        ("direct", True, False),
        ("bw2calc", True, True),
    ],
)
def test_manufacture_and_retirement_keep_calendar_beyond_matrix_horizon(
    solver_mode, attribute_to_roots, interpolate_annual
):
    checks = validate(
        Path(__file__).parent / "data/stock_lifecycle",
        solver_mode=solver_mode,
        attribute_to_roots=attribute_to_roots,
        interpolate_annual=interpolate_annual,
    )
    assert len(checks) == 8


@pytest.mark.parametrize("adaptive_cutoff", [1e-4, 10.0])
@pytest.mark.parametrize("interpolate_annual", [False, True])
def test_adaptive_lifecycle_keeps_past_and_future_event_dates(
    adaptive_cutoff, interpolate_annual
):
    checks = validate(
        Path(__file__).parent / "data/stock_lifecycle",
        interpolate_annual=interpolate_annual,
        adaptive_cutoff=adaptive_cutoff,
    )
    assert len(checks) == 8
