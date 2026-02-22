User Guide
==========

This guide covers the main workflows in TRAILS: loading data packages, running
temporal LCA, and interpreting results.

Data packages
-------------

TRAILS expects a Frictionless data package (``datapackage.json``) with
scenario-indexed technosphere and biosphere matrices plus index metadata.
Packages produced by ``premise`` are supported out of the box.

Key components include:

* **Technosphere matrix** (``A``): exchanges between activities/products
* **Biosphere matrix** (``B``): biosphere flows per activity
* **Activity indices**: metadata for activities (name, location, unit)
* **Biosphere indices**: metadata for flows (name, compartment, unit)
* **Temporal exchanges** (optional): distribution metadata for delayed flows

Loading TRAILS
--------------

.. code-block:: python

    from datapackage import Package
    from trails import Trails

    package = Package("path/to/datapackage.json")

    # interpolate_annual=True expands scenario slices to full annual resolution
    trails = Trails(package, interpolate_annual=True)

After initialization, you can access scenario labels and metadata:

.. code-block:: python

    print(trails.scenario_labels)
    activity_indices = next(iter(trails.activity_indices.values()))
    print(list(activity_indices.items())[:5])

Selecting activities
--------------------

Activities are referenced by integer indices from the metadata. A typical
workflow is to select an activity and store the index for repeated calls:

.. code-block:: python

    activity_indices = next(iter(trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))

Running temporal LCA
--------------------

The primary entry point is ``trails.lca.lca``:

.. code-block:: python

    from trails import lca, get_lcia_method_names

    method = get_lcia_method_names(ei_version="3.11")[0]

    # Run temporal routing (builds the traversal graph)
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        amount=1.0,
        max_depth=3,
        min_amount=1e-18,
        show_progress=True,
        debug=False,
        attribute_to_roots=True,
    )

    # Run temporal LCA (stores results on the Trails instance)
    lca(
        trails=trails,
        methods=[method],
    )

Temporal LCA results are stored on the Trails instance. Use ``trails.scores``
for impact scores (when compute_score=True) and ``trails.inventory`` or
``trails.characterized_inventory`` for time-resolved inventories.

Importing Excel inventories
---------------------------

You can import user-provided inventories from Excel using ``bw2io``. When you
omit ``year`` and ``scenario_label``, the exchanges are applied to all template
years and interpolated across annual years.

.. code-block:: python

    from trails import Trails

    trails = Trails(package)
    trails.import_excel_inventory("path/to/inventory.xlsx")

    # Target a single scenario slice instead
    trails.import_excel_inventory("path/to/inventory.xlsx", year=2020)

Example inventory: ``examples/lci-pass_cars.xlsx``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file demonstrates the ``bw2io`` Excel format that TRAILS expects:

* **Activity section** (top of the sheet): key/value pairs that define the
  activity being imported.
* **Exchanges section** (table starting after the ``Exchanges`` row): each row
  is an exchange linked to the activity.

Activity section fields (column A = field name, column B = value):

* ``Activity``: activity name (string).
* ``reference product``: reference product (string).
* ``location``: location code (string).
* ``unit``: activity unit (string).
* Additional fields (e.g., ``lifetime [km]``, ``average age [year]``) are
  allowed and preserved as metadata.

Exchange table columns (from the header row):

* ``name``: exchange name (string).
* ``amount``: base exchange amount (float). Used when no year-specific columns
  are provided.
* **Year-specific columns** (e.g., ``2020``, ``2040``, ``2060``):
  optional numeric columns used as **year-specific amounts**. TRAILS writes
  these values into the corresponding years in the A/B matrices and then
  linearly interpolates between them for intermediate years. For years outside
  the provided range, the nearest endpoint value is used (clamped).
* ``location``: exchange location (string). Required for technosphere/production
  linking.
* ``categories``: biosphere categories (tuple or ``compartment/subcompartment``).
* ``unit``: exchange unit (string).
* ``type``: exchange type: ``production``, ``technosphere``, or ``biosphere``.
* ``reference product``: required for technosphere/production exchanges.
* ``temporal_distribution``: integer code for the temporal distribution.
* ``temporal_loc``: location parameter (mean/median/mode depending on distribution).
* ``temporal_scale``: scale parameter (e.g., stddev/sigma).
* ``temporal_min`` / ``temporal_max``: integer offsets defining the support.
* ``temporal_offsets``: JSON list of integer offsets (used by discrete empirical).
* ``temporal_weights``: JSON list of pulse weights (used by discrete empirical).
* ``temporal_amount_source``: either ``port`` or ``matrix``:
  ``port`` applies the **ported exchange amount** over time; ``matrix`` uses
  the **matrix-stored values** in other years.
* ``comment``: optional notes.

TRAILS-specific temporal fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The temporal columns correspond directly to ``TemporalExchange``:

* ``temporal_distribution``: distribution code (1=discrete, 2=lognormal,
  3=normal, 4=uniform, 5=triangular, 6=discrete empirical).
* ``temporal_loc`` and ``temporal_scale``: distribution parameters.
* ``temporal_min`` and ``temporal_max``: inclusive integer offsets.
* ``temporal_offsets`` and ``temporal_weights``:
  JSON-list pulse definitions for ``distribution=6``. Example:
  ``[0, 5, 12]`` and ``[0.5, 0.3, 0.2]``.
* ``temporal_amount_source``:
  - ``port``: uses the exchange amount for all years (ported value).
  - ``matrix``: uses the matrix values for the selected years.

Temporal distributions
----------------------

Temporal distributions control how exchanges are spread across impact years.
TRAILS can use the distribution data included in the package, or collapse the
effects into scalar multipliers for a static approximation:

.. code-block:: python

    # Use full temporal distributions (default)
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        max_depth=3,
        min_amount=1e-18,
        use_temporal_distributions=True,
    )
    lca(
        trails=trails,
        methods=[method],
    )

    # Collapse temporal distributions for a faster static view
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        max_depth=3,
        min_amount=1e-18,
        use_temporal_distributions=False,
    )
    lca(
        trails=trails,
        methods=[method],
    )

Plotting results
----------------

Use ``plot_temporal_scores`` to visualize the time series and contribution by
first-level suppliers:

.. code-block:: python

    from trails import plot_temporal_scores

    fig = plot_temporal_scores(
        trails,
        method_label=method,
        stacked=True,
        show_cumulative_axis=True,
    )
    fig.show()

Interpreting outputs
--------------------

Impact time series can be accessed from ``trails.scores`` (if computed) or
from ``trails.characterized_inventory`` after characterization.

Troubleshooting and diagnostics
-------------------------------

* If a requested year is not available in the data package, TRAILS will snap to
  the nearest available year and emit a warning.
* Set ``debug=True`` to enable detailed logging and retain additional
  diagnostics on the Trails instance.


FaIR_ radiative forcing
-----------------------

TRAILS integrates with the FaIR_ climate model to convert time-resolved
inventories into radiative forcing and temperature anomalies. The workflow runs
a baseline FaIR_ scenario from the bundled REMIND/FaIR_ emissions data, then
performs per-species perturbation runs derived from the Trails inventory.
Positive and negative emissions are treated separately to preserve long-lived
CO2 tails for both uptake and release. Results are allocated to root activities
using cumulative signed emissions for each (flow, root) pair, and summarized
across all FaIR_ configurations as quantiles (2.5, 25, 50, 75, 97.5).

.. _FaIR: https://github.com/OMS-NetZero/FAIR

.. code-block:: python

    from trails.fair_rf import run_fair_delta_rf

    rf = run_fair_delta_rf(
        trails,
        scenario="high-extension",
    )

The outputs are stored on:

* ``trails.instant_radiative_forcing`` with dims ``(quantile, year, flow, root activity)``
* ``trails.delta_temperature`` with dims ``(quantile, year, flow, root activity)``

Visualization helpers default to the 50th quantile:

.. code-block:: python

    from trails import plot_rf, plot_temp

    plot_rf(trails, year_range=(2000, 2100))
    plot_temp(trails, year_range=(2000, 2100))
