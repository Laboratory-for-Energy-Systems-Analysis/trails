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

The core LCA routine produces two complementary views:

* **results_by_solve_year**: diagnostic information about each solved year
  (demand vectors, injected supply pulses, and metadata for debugging).
* **results_by_impact_year**: the impact time series that you plot and analyze.

Because impacts are booked in impact years, TRAILS provides a direct answer to
questions such as:

* *When do impacts occur?*
* *Which upstream suppliers contribute most over time?*
* *How do scenario transitions shift the impact profile?*
