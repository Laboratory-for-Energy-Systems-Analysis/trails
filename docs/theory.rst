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

This traversal is performed by ``trails.temporal_routing(...)``. ``trails.lci()``
then consumes the routing graph to build one time-resolved inventory, and
``trails.lcia(...)`` characterizes that inventory without repeating the solves.

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
reference-product demand amount, and ``s_m(a, y)`` is the static score for one
unit of the activity's reference product under method ``m``.

A branch is routed explicitly only while this potential remains above the
chosen cutoff. Branches below the cutoff become frontier nodes and are still
included in the final calculation through the year-wise matrix solve. The
adaptive cutoff therefore controls how much of the supply-chain graph is made
explicit; it does not discard the remaining demand.

The cutoff is configured with ``adaptive_relative_score_cutoff``, a
dimensionless fraction of the functional unit's static score potential. For
example, ``1e-4`` stops branches whose estimated static potential is at most
0.01% of the functional unit potential.

When ``max_depth`` is omitted, ``temporal_routing()`` uses adaptive routing with
``max_depth=None`` and ``adaptive_relative_score_cutoff=1e-4``; explicit regular
``adaptive_methods`` are required. Passing an integer ``max_depth`` without an
adaptive cutoff selects fixed-depth routing. ``max_depth`` can also be combined
with the adaptive relative cutoff as a hard cap. The static activity scores are
cached using matrix and LCIA-data fingerprints so repeated adaptive runs over
the same interpolated data package avoid recomputing the upfront screening
intensities.

In the recommended workflow, regular screening methods are passed explicitly as
``adaptive_methods``. Final regular or EDGES methods are supplied independently
to ``lcia()``. EDGES methods are currently limited to final characterization and
cannot provide adaptive routing potentials.

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

``lci()`` stores the temporal inventory on ``trails.inventory``. ``lcia()``
stores the current impact time series on ``trails.scores`` and regular methods
also expose ``trails.characterized_inventory``. Multiple LCIA results are kept
in ``trails.lcia_results``.

Inventory storage is selected independently from characterization. With
``inventory_backend="auto"``, small inventories remain eager sparse COO arrays.
Eligible large root-attributed direct or iterative inventories use a factorized
representation of annual supply matrices, biosphere coefficients, and temporal
kernels; other large inventories use bounded chunked sparse storage. Both are
presented through the same xarray dimensions, and repeated ``lcia()`` calls do
not rerun the annual linear systems.

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

Fixed-window CO2 pulse equivalence
----------------------------------

``run_fair_co2_pulse_equivalents`` compares the complete inventory perturbation
with a reference CO2 pulse under the same FaIR background scenario. For every
FaIR configuration, the integrated-radiative-forcing equivalent is

.. math::

   M_{\mathrm{CO2,eq}} = M_{\mathrm{ref}}
   \frac{\int_{t_0}^{t_1} \Delta RF_{\mathrm{LCA}}(t)\,dt}
        {\int_{t_0}^{t_1} \Delta RF_{\mathrm{CO2\ pulse}}(t)\,dt}.

The configuration-level ratios are summarized only after the three FaIR runs
(baseline, inventory perturbation, and reference pulse). A negative equivalent
means that the inventory produces net cooling relative to the background over
the selected window. The result is conditional on the background scenario,
pulse year, and window and is therefore distinct from a GWP metric or a physical
carbon-storage efficiency.

.. _FaIR: https://github.com/OMS-NetZero/FAIR
