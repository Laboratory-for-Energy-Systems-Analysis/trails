# temporal_distributions.py

from dataclasses import dataclass
from typing import Optional, Iterable, Tuple, Sequence
import numpy as np

import numpy as np
from math import erf, sqrt, exp

import logging

logger = logging.getLogger(__name__)


@dataclass
class TemporalExchange:
    """Metadata for a single temporally-distributed exchange.


    - :param distribution: Integer code for the distribution shape.
    - 1: discrete (all mass at loc)
    - 2: lognormal
    - 3: normal
    - 4: uniform
    - 5: triangular
    - 6: discrete empirical (explicit offsets + weights)
    :param loc: Location parameter (mean, median, or mode depending on distribution).
    :param scale: Scale parameter (stddev for normal, sigma for lognormal).
    :param offset_min: Minimum integer offset (inclusive).
    :param offset_max: Maximum integer offset (inclusive).
    :param amount_source: Source of the amount for each pulse year."""

    distribution: int
    loc: Optional[float]
    scale: Optional[float]
    offset_min: int
    offset_max: int
    # "port" (default): use CSV row value scaled by weights
    # "matrix": ignore CSV row value; at each pulse year use interpolated A/B[t, i, j] or B[t, i, f]
    amount_source: str = "port"
    offsets: Optional[Sequence[int]] = None
    weights: Optional[Sequence[float]] = None


class TemporalDistribution:
    """Turn a TemporalExchange into discrete (offset, weight) pairs."""

    def __init__(self, tex: TemporalExchange) -> None:
        """init  .

        :param tex: Value for `tex`.
        :type tex: TemporalExchange"""
        self.tex = tex

    # ------------------------------------------------------------------
    # Individual distribution helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _default_sigma(offsets: np.ndarray, scale: Optional[float]) -> float:
        """default sigma.

        :param offsets: Value for `offsets`.
        :type offsets: np.ndarray
        :param scale: Value for `scale`.
        :type scale: Optional[float]
        :returns: Return value.
        :rtype: float"""
        if scale is not None and scale > 0:
            return float(scale)

        # Heuristic: spread over about half the range
        rng = float(offsets.max() - offsets.min())
        if rng <= 0:
            return 1.0
        return max(rng / 4.0, 1e-6)

    @staticmethod
    def _triangular_weights(offsets: np.ndarray, loc: Optional[float]) -> np.ndarray:
        """triangular weights.

        :param offsets: Value for `offsets`.
        :type offsets: np.ndarray
        :param loc: Value for `loc`.
        :type loc: Optional[float]
        :returns: Return value.
        :rtype: np.ndarray"""
        if loc is None:
            # If no loc given, fall back to symmetric around the midpoint
            loc = 0.0

        offsets = offsets.astype(float)
        distances = np.abs(offsets - float(loc))

        max_d = distances.max()
        if max_d == 0.0:
            # all offsets equal to loc -> put all mass there
            w = np.zeros_like(offsets, dtype=float)
            w[:] = 1.0
            return w

        # Linear "tent" function: max_d - distance
        w = max_d - distances
        w[w < 0.0] = 0.0
        return w

    @staticmethod
    def _normal_weights(
        offsets: np.ndarray,
        loc: float,
        scale: float,
        offset_min: int,
        offset_max: int,
    ) -> np.ndarray:
        """normal weights.

        :param offsets: Value for `offsets`.
        :type offsets: np.ndarray
        :param loc: Value for `loc`.
        :type loc: float
        :param scale: Value for `scale`.
        :type scale: float
        :param offset_min: Value for `offset_min`.
        :type offset_min: int
        :param offset_max: Value for `offset_max`.
        :type offset_max: int
        :returns: Return value.
        :rtype: np.ndarray"""
        if scale is None or scale <= 0:
            scale = 1.0

        offsets = np.asarray(offsets, dtype=float)

        # Normal PDF evaluated at integer offsets
        pdf_vals = np.exp(-0.5 * ((offsets - loc) / scale) ** 2)

        # --- Normalization for TRUNCATED NORMAL BETWEEN [offset_min, offset_max] ---

        # CDF helper
        def normal_cdf(x: float) -> float:
            """Normal cdf.

            :param x: Value for `x`.
            :type x: float
            :returns: Return value.
            :rtype: float"""
            return 0.5 * (1 + erf((x - loc) / (scale * sqrt(2))))

        # Continuous probability mass inside the allowed range
        Z = normal_cdf(offset_max + 0.5) - normal_cdf(offset_min - 0.5)

        if Z <= 0:
            # fallback to uniform if something is wrong
            return np.ones_like(offsets, dtype=float)

        # Discrete normalization to preserve total mass = 1
        pdf_vals /= pdf_vals.sum()  # sum to 1 over sampled offsets
        pdf_vals /= Z  # rescale to match truncated continuous total mass

        return pdf_vals

    def _lognormal_weights(
        self, offsets: np.ndarray, loc: Optional[float], scale: Optional[float]
    ) -> np.ndarray:
        """lognormal weights.

        :param offsets: Value for `offsets`.
        :type offsets: np.ndarray
        :param loc: Value for `loc`.
        :type loc: Optional[float]
        :param scale: Value for `scale`.
        :type scale: Optional[float]
        :returns: Return value.
        :rtype: np.ndarray"""
        x = offsets.astype(float)
        mask = x > 0
        if not mask.any():
            # No positive offsets to assign lognormal mass -> uniform fallback
            return np.ones_like(offsets, dtype=float)

        if loc is None or loc <= 0:
            # Fallback: use median at geometric mid of positive range
            pos_vals = x[mask]
            loc = float(np.sqrt(pos_vals.min() * pos_vals.max()))

        sigma = self._default_sigma(x[mask], scale)
        mu = float(np.log(loc))

        w = np.zeros_like(x, dtype=float)
        # lognormal pdf, up to a constant factor:
        z = (np.log(x[mask]) - mu) / sigma
        # divide by x to get the right qualitative shape; constant factor ignored
        w[mask] = np.exp(-0.5 * z * z) / x[mask]

        if not np.isfinite(w).any() or w.sum() <= 0:
            return np.ones_like(offsets, dtype=float)
        return w

    @staticmethod
    def _discrete_weights(offsets: np.ndarray, loc: Optional[float]) -> np.ndarray:
        """discrete weights.

        :param offsets: Value for `offsets`.
        :type offsets: np.ndarray
        :param loc: Value for `loc`.
        :type loc: Optional[float]
        :returns: Return value.
        :rtype: np.ndarray"""
        if loc is None:
            # try to anchor at 0 if possible
            if offsets.min() <= 0 <= offsets.max():
                target = 0
            else:
                # nearest boundary
                target = offsets[np.argmin(np.abs(offsets))]
        else:
            target = int(round(loc))
            # clip to available range
            target = int(max(offsets.min(), min(offsets.max(), target)))

        w = np.zeros_like(offsets, dtype=float)
        idx = int(np.argmin(np.abs(offsets - target)))
        w[idx] = 1.0
        return w

    @staticmethod
    def _discrete_empirical_weights(
        offsets: np.ndarray,
        pulse_offsets: Optional[Sequence[int]],
        pulse_weights: Optional[Sequence[float]],
        *,
        fallback_loc: Optional[float],
    ) -> np.ndarray:
        """Build explicit discrete pulse weights over integer offsets."""
        if pulse_offsets is None or pulse_weights is None:
            return TemporalDistribution._discrete_weights(offsets, fallback_loc)
        if len(pulse_offsets) == 0 or len(pulse_weights) == 0:
            return TemporalDistribution._discrete_weights(offsets, fallback_loc)

        if len(pulse_offsets) != len(pulse_weights):
            return TemporalDistribution._discrete_weights(offsets, fallback_loc)

        w = np.zeros_like(offsets, dtype=float)
        offset_to_idx = {int(off): idx for idx, off in enumerate(offsets.tolist())}

        for off, wei in zip(pulse_offsets, pulse_weights):
            try:
                off_i = int(off)
                wei_f = float(wei)
            except Exception:
                continue
            if not np.isfinite(wei_f) or wei_f < 0.0:
                continue
            idx = offset_to_idx.get(off_i)
            if idx is None:
                continue
            w[idx] += wei_f

        if not np.isfinite(w).any() or w.sum() <= 0:
            return TemporalDistribution._discrete_weights(offsets, fallback_loc)
        return w

    def iter_offsets_and_weights(
        self, debug: bool = False
    ) -> Iterable[Tuple[int, float]]:
        """Iter offsets and weights.

        :param debug: Value for `debug`.
        :type debug: bool
        :yields: Yielded values.
        :returns: Return value.
        :rtype: Iterable[Tuple[int, float]]"""
        t = self.tex

        dist = int(t.distribution)

        # For discrete empirical distributions, explicit pulse offsets define support.
        if dist == 6 and getattr(t, "offsets", None):
            try:
                explicit_offsets = [int(round(float(o))) for o in t.offsets]  # type: ignore[arg-type]
            except Exception:
                explicit_offsets = []
            if explicit_offsets:
                offsets = np.array(sorted(set(explicit_offsets)), dtype=int)
            else:
                offsets = np.arange(int(t.offset_min), int(t.offset_max) + 1, dtype=int)
        else:
            offsets = np.arange(int(t.offset_min), int(t.offset_max) + 1, dtype=int)

        if offsets.size == 0:
            return iter(())

        # 1) Base (unnormalized) weights
        if dist == 5:
            weights = self._triangular_weights(offsets, t.loc)
        elif dist == 2:
            weights = self._lognormal_weights(offsets, t.loc, t.scale)
        elif dist == 3:
            weights = self._normal_weights(
                offsets=offsets,
                loc=float(t.loc) if t.loc is not None else 0.0,
                scale=float(t.scale) if (t.scale is not None and t.scale > 0) else 1.0,
                offset_min=int(t.offset_min),
                offset_max=int(t.offset_max),
            )
        elif dist == 4:
            weights = np.ones_like(offsets, dtype=float)
        elif dist == 1:
            weights = self._discrete_weights(offsets, t.loc)
        elif dist == 6:
            weights = self._discrete_empirical_weights(
                offsets,
                getattr(t, "offsets", None),
                getattr(t, "weights", None),
                fallback_loc=t.loc,
            )
        else:
            # Unknown distribution -> uniform fallback
            weights = np.ones_like(offsets, dtype=float)

        weights = np.asarray(weights, dtype=float)
        total = float(weights.sum())
        if not np.isfinite(total) or total <= 0.0:
            if debug:
                logger.warning(
                    "TemporalDistribution: invalid total weight (dist=%s, offsets=[%d..%d]) -> empty",
                    dist,
                    int(t.offset_min),
                    int(t.offset_max),
                )
            return iter(())

        # 2) Normalize
        weights /= total

        # 3) Yield
        for k, w in zip(offsets, weights):
            w = float(w)
            if w != 0.0:
                yield int(k), w


def resolve_temporal_offset_bounds(
    *,
    distribution: int,
    loc: Optional[float],
    offset_min: Optional[int],
    offset_max: Optional[int],
) -> tuple[int, int]:
    """Resolve optional temporal offset bounds to integer support.

    Discrete pulse distributions use ``loc`` as the pulse offset. If no explicit
    ``temporal_min``/``temporal_max`` support is provided, the support must
    therefore be the pulse offset itself, not the generic 0-year fallback.
    """
    if int(distribution) == 1 and offset_min is None and offset_max is None:
        if loc is not None and np.isfinite(float(loc)):
            target = int(round(float(loc)))
            return target, target

    off_min = 0 if offset_min is None else int(offset_min)
    off_max = 0 if offset_max is None else int(offset_max)
    return off_min, off_max
