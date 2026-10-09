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
