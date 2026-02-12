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
