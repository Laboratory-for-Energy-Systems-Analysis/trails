import math
import numpy as np
import pytest

from trails.temporal_distributions import (
    TemporalDistribution,
    TemporalExchange,
    resolve_temporal_offset_bounds,
)


def test_default_sigma() -> None:
    """Verify default sigma calculation.

    :returns: None.
    :rtype: None
    """
    offsets = np.array([0, 4])
    assert TemporalDistribution._default_sigma(offsets, None) == 1.0
    assert TemporalDistribution._default_sigma(offsets, 2.0) == 2.0


def test_triangular_weights_centered() -> None:
    """Verify centered triangular weights.

    :returns: None.
    :rtype: None
    """
    offsets = np.array([0, 1, 2])
    weights = TemporalDistribution._triangular_weights(offsets, loc=1.0)
    assert weights.tolist() == [0.0, 1.0, 0.0]


def test_normal_weights_truncated() -> None:
    """Verify truncated normal weights.

    :returns: None.
    :rtype: None
    """
    offsets = np.array([-1, 0, 1])
    weights = TemporalDistribution._normal_weights(
        offsets=offsets, loc=0.0, scale=1.0, offset_min=-1, offset_max=1
    )
    z = 0.5 * (1 + math.erf(1.5 / math.sqrt(2))) - 0.5 * (
        1 + math.erf(-1.5 / math.sqrt(2))
    )
    assert np.all(weights > 0)
    assert np.isclose(weights.sum(), 1.0 / z, atol=1e-6)


def test_lognormal_weights_fallback() -> None:
    """Verify lognormal weights fallback behavior.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=2, loc=None, scale=None, offset_min=-1, offset_max=1
    )
    dist = TemporalDistribution(tex)
    weights = dist._lognormal_weights(np.array([-1, 0, 1]), loc=None, scale=None)
    assert weights[0] == 0.0
    assert weights[1] == 0.0
    assert weights[2] > 0.0


def test_discrete_weights_default_zero() -> None:
    """Verify discrete weights default to zero offset.

    :returns: None.
    :rtype: None
    """
    offsets = np.array([-1, 0, 1])
    weights = TemporalDistribution._discrete_weights(offsets, loc=None)
    assert weights.tolist() == [0.0, 1.0, 0.0]


def test_discrete_bounds_default_to_loc_when_support_missing() -> None:
    """Verify discrete pulses use loc when min/max are omitted."""
    off_min, off_max = resolve_temporal_offset_bounds(
        distribution=1,
        loc=16.0,
        offset_min=None,
        offset_max=None,
    )
    tex = TemporalExchange(
        distribution=1,
        loc=16.0,
        scale=None,
        offset_min=off_min,
        offset_max=off_max,
    )
    assert list(TemporalDistribution(tex).iter_offsets_and_weights()) == [
        (16, pytest.approx(1.0))
    ]


def test_iter_offsets_and_weights_with_scaling() -> None:
    """Verify offsets and weights for uniform distribution.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=4,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=1,
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [(0, pytest.approx(0.5)), (1, pytest.approx(0.5))]


def test_iter_offsets_and_weights_unknown_scale_mode() -> None:
    """Verify unknown distributions fall back to uniform.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=99,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=1,
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [(0, pytest.approx(0.5)), (1, pytest.approx(0.5))]


def test_scale_factor_modes_and_clip() -> None:
    """Verify discrete distribution returns a single pulse.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=1,
        loc=0,
        scale=None,
        offset_min=0,
        offset_max=2,
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [(0, pytest.approx(1.0))]


def test_discrete_empirical_distribution() -> None:
    """Verify explicit discrete pulses with normalization.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=3,
        offsets=[0, 2, 2, 3],
        weights=[0.2, 0.3, 0.1, 0.4],
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [
        (0, pytest.approx(0.2)),
        (2, pytest.approx(0.4)),
        (3, pytest.approx(0.4)),
    ]


def test_discrete_empirical_uses_explicit_offsets_not_min_max() -> None:
    """Verify explicit pulse offsets are used even if min/max differ.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=6,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        offsets=[-1, 9],
        weights=[0.5, 0.5],
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [(-1, pytest.approx(0.5)), (9, pytest.approx(0.5))]
