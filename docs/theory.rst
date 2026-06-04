Theory
======

TRAILS models environmental impacts as *time-aware* exchanges moving through a
supply-chain graph. Instead of collapsing all exchanges into a single static
inventory, it carries time-indexed technosphere and biosphere matrices and
allows temporal distributions to shift when impacts are booked.

Temporal LCIA in TRAILS
-----------------------

At a high level, TRAILS works with time-indexed matrices:

* :math:`A_t`: technosphere exchanges by scenario/year :math:`t`
* :math:`B_t`: biosphere exchanges by scenario/year :math:`t`

For a given functional unit in year :math:`y_0`, TRAILS traverses the
technosphere graph over time, collecting a **frontier** of activity demands in
future years. Each frontier slice is solved with the corresponding
:math:`A_t` and combined with :math:`B_t` to construct inventories for impact
years.

This traversal is performed by ``trails.temporal_routing(...)`` before running
``trails.lca(...)``, which consumes the stored routing graph to build the
time-resolved inventory and impact scores.

The final output is a mapping from impact year to characterized impact scores
(plus attribution to first-level suppliers when available).

Adaptive score-potential routing
--------------------------------

Fixed-depth routing expands every branch until it reaches the same graph depth
or an absolute demand cutoff. This is simple and reproducible, but it can spend
substantial time expanding low-impact branches while stopping high-impact
branches at the same depth.

Adaptive routing adds an impact-potential screening step. Before routing,
TRAILS computes static LCIA activity scores for the selected adaptive methods
and scenario years. For a candidate child demand, the routing step estimates
the branch potential as:

.. math::

   p(a, y, q) = |q| \max_m |s_m(a, y)|

where ``a`` is the activity, ``y`` is the target year, ``q`` is the routed
demand amount, and ``s_m(a, y)`` is the static score for one unit of the
activity's reference product under method ``m``.

A branch is routed explicitly only while this potential remains above the
chosen cutoff. Branches below the cutoff become frontier nodes and are still
included in the final calculation through the year-wise matrix solve. The
adaptive cutoff therefore controls how much of the supply-chain graph is made
explicit; it does not discard the remaining demand.

Two cutoff forms are supported:

* ``adaptive_score_cutoff``: an absolute score-potential threshold in the LCIA
  unit of the selected method.
* ``adaptive_relative_score_cutoff``: a dimensionless fraction of the root
  activity's static score potential. For example, ``1e-4`` stops branches whose
  estimated static potential is at most 0.01% of the root potential.

``max_depth`` can still be used with adaptive routing as a hard cap. Setting
``max_depth=None`` removes this cap and lets the adaptive cutoff, minimum
amount, and graph leaves determine the explicit routing depth. The static
activity scores are cached using matrix and LCIA-data fingerprints so repeated
adaptive runs over the same interpolated data package avoid recomputing the
upfront screening intensities.

In the recommended workflow, regular LCIA methods are configured once with
``Trails(..., methods=[...], ei_version="...")``. Adaptive routing uses those
methods when a cutoff is provided and ``adaptive_methods`` is omitted, and
``lca(trails)`` reuses the same methods for final scoring. Call-level methods
can still be supplied to override this default. EDGES methods are currently
limited to final scoring and cannot provide adaptive routing potentials.

Temporal exchange distributions
-------------------------------

Exchanges can carry temporal distributions rather than single-point events.
A temporal distribution defines how a flow is spread over offsets relative to
its parent activity year. TRAILS supports deterministic and distributional
forms (e.g., fixed offset, uniform spread, normal-like kernels), allowing:

* **delayed emissions** (e.g., use-phase emissions)
* **staged production** (e.g., infrastructure build-outs)
* **temporal aggregation** for long-lived supply chains

During traversal, these distributions shift demand and emissions into impact
years. When desired, TRAILS can collapse distributions into scalar multipliers
for faster static approximations.

Scenario-aware inventories
--------------------------

TRAILS is designed for prospective LCA with scenario data packages (e.g.,
from ``premise``). Each scenario year contains its own technosphere and
biosphere slices. TRAILS can interpolate to annual resolution, or snap to the
nearest available year, ensuring that the inventory and characterization logic
remains consistent across time horizons.

Impact attribution across time
------------------------------

The core LCA routine stores its outputs on the Trails instance. Use
``trails.scores`` for impact time series (when compute_score=True) and
``trails.inventory`` / ``trails.characterized_inventory`` for time-resolved
inventories and attribution.

Because impacts are booked in impact years, TRAILS provides a direct answer to
questions such as:

* *When do impacts occur?*
* *Which upstream suppliers contribute most over time?*
* *How do scenario transitions shift the impact profile?*


FaIR_ climate model integration
-------------------------------

TRAILS can translate time-resolved inventories into radiative forcing and
temperature anomalies using the FaIR_ climate model. The method runs a baseline
FaIR_ scenario from the bundled IAMC emissions data, then applies per-species
perturbations derived from the Trails inventory. Positive and negative
emissions are treated separately to preserve long-lived CO2 tails for both
uptake and release. Results are allocated to root activities using cumulative
signed emissions for each (flow, root) pair and stored as
``trails.instant_radiative_forcing`` and ``trails.delta_temperature``. Both
arrays are stored across FaIR_ configurations as quantiles (2.5, 25, 50, 75,
97.5) with dims ``(quantile, year, flow, root activity)``.

.. _FaIR: https://github.com/OMS-NetZero/FAIR
