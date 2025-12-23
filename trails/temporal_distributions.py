# temporal_distributions.py

from dataclasses import dataclass
from typing import Optional, Iterable, Tuple
import numpy as np

import numpy as np
from math import erf, sqrt, exp, pi

import logging

logger = logging.getLogger(__name__)


@dataclass
class TemporalExchange:
    """
    Metadata for a single temporally-distributed exchange.

    :param distribution: Integer code for the distribution shape.
        - 1: discrete (all mass at loc)
        - 2: lognormal
        - 3: normal
        - 4: uniform
        - 5: triangular
    :param loc: Location parameter (mean, median, or mode depending on distribution).
    :param scale: Scale parameter (stddev for normal, sigma for lognormal).
    :param offset_min: Minimum integer offset (inclusive).
    :param offset_max: Maximum integer offset (inclusive).
    :param amount_source: Source of the amount for each pulse year.
    """

    distribution: int
    loc: Optional[float]
    scale: Optional[float]
    offset_min: int
    offset_max: int
    # "port" (default): use CSV row value scaled by weights
    # "matrix": ignore CSV row value; at each pulse year use interpolated A/B[t, i, j] or B[t, i, f]
    amount_source: str = "port"


class TemporalDistribution:
    """Turn a TemporalExchange into discrete (offset, weight) pairs."""

    def __init__(self, tex: TemporalExchange) -> None:
        """Initialize the distribution wrapper.

        :param tex: Temporal exchange metadata to interpret.
        :type tex: TemporalExchange
        """
        self.tex = tex


    # ------------------------------------------------------------------
    # Individual distribution helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _default_sigma(offsets: np.ndarray, scale: Optional[float]) -> float:
        """Choose a reasonable sigma if scale is None or invalid.

        :param offsets: Offset values to inspect.
        :type offsets: numpy.ndarray
        :param scale: Scale parameter if provided.
        :type scale: float | None
        :returns: Fallback sigma value.
        :rtype: float
        """
        if scale is not None and scale > 0:
            return float(scale)

        # Heuristic: spread over about half the range
        rng = float(offsets.max() - offsets.min())
        if rng <= 0:
            return 1.0
        return max(rng / 4.0, 1e-6)

    @staticmethod
    def _triangular_weights(offsets: np.ndarray, loc: Optional[float]) -> np.ndarray:
        """Build a discrete triangular shape over given integer offsets.

        :param offsets: Offset values to weight.
        :type offsets: numpy.ndarray
        :param loc: Mode location for the distribution.
        :type loc: float | None
        :returns: Unnormalized triangular weights.
        :rtype: numpy.ndarray
        """
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
        """Compute truncated normal weights over the integer offsets.

        :param offsets: Offset values to weight.
        :type offsets: numpy.ndarray
        :param loc: Mean of the normal distribution.
        :type loc: float
        :param scale: Standard deviation for the normal distribution.
        :type scale: float
        :param offset_min: Minimum offset bound.
        :type offset_min: int
        :param offset_max: Maximum offset bound.
        :type offset_max: int
        :returns: Truncated normal weights over offsets.
        :rtype: numpy.ndarray
        """
        if scale is None or scale <= 0:
            scale = 1.0

        offsets = np.asarray(offsets, dtype=float)

        # Normal PDF evaluated at integer offsets
        pdf_vals = np.exp(-0.5 * ((offsets - loc) / scale) ** 2)

        # --- Normalization for TRUNCATED NORMAL BETWEEN [offset_min, offset_max] ---

        # CDF helper
        def normal_cdf(x: float) -> float:
            """Return the standard normal CDF at ``x`` for the local parameters.

            :param x: Input value.
            :type x: float
            :returns: CDF value for the truncated normal helper.
            :rtype: float
            """
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
        """Compute discrete lognormal weights for integer offsets.

        Only positive offsets receive mass. ``loc`` is interpreted as the
        median in offset space when provided.

        :param offsets: Offset values to weight.
        :type offsets: numpy.ndarray
        :param loc: Median location in offset space.
        :type loc: float | None
        :param scale: Sigma in log-space.
        :type scale: float | None
        :returns: Lognormal weights over offsets.
        :rtype: numpy.ndarray
        """
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
        """Build a discrete (Dirac-like) distribution over offsets.

        :param offsets: Offset values to weight.
        :type offsets: numpy.ndarray
        :param loc: Location to concentrate weight around.
        :type loc: float | None
        :returns: Discrete weights over offsets.
        :rtype: numpy.ndarray
        """
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

    def iter_offsets_and_weights(self, debug: bool = False) -> Iterable[Tuple[int, float]]:
        """
        Yield (offset, weight) pairs for the temporal distribution.

        Scaling modes have been removed: weights are determined purely by the
        selected distribution and then normalized to sum to 1 over the integer
        offsets [offset_min, offset_max].
        """
        t = self.tex

        offsets = np.arange(int(t.offset_min), int(t.offset_max) + 1, dtype=int)
        if offsets.size == 0:
            return iter(())

        dist = int(t.distribution)

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
