import math
import numpy as np
import pytest

from trails.temporal_distributions import TemporalDistribution, TemporalExchange


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


def test_iter_offsets_and_weights_with_scaling() -> None:
    """Verify offsets and weights with scaling.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=4,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=1,
        scale_mode="linear",
        scale_base=1.0,
        scale_rate=1.0,
    )
    td = TemporalDistribution(tex)
    results = list(td.iter_offsets_and_weights())
    assert results == [(0, pytest.approx(1 / 3)), (1, pytest.approx(2 / 3))]


def test_iter_offsets_and_weights_unknown_scale_mode() -> None:
    """Verify unknown scale modes raise errors.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=4,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        scale_mode="mystery",
    )
    td = TemporalDistribution(tex)
    with pytest.raises(ValueError):
        list(td.iter_offsets_and_weights())


def test_scale_factor_modes_and_clip() -> None:
    """Verify scale factor modes and clipping.

    :returns: None.
    :rtype: None
    """
    tex = TemporalExchange(
        distribution=1,
        loc=None,
        scale=None,
        offset_min=0,
        offset_max=0,
        scale_mode="compound",
        scale_base=2.0,
        scale_rate=0.1,
    )
    td = TemporalDistribution(tex)
    assert td.scale_factor(2) == pytest.approx(2.42)
    assert td.scale_factor(2, clip=(0.0, 2.0)) == 2.0
