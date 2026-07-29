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
    from trails import Trails, get_lcia_method_names

    package = Package("path/to/datapackage.json")
    method = get_lcia_method_names(ei_version="3.11")[0]

    # interpolate_annual=True expands scenario slices to annual resolution.
    # Default annual bounds are [min_year-1, max_year+1].
    trails = Trails(
        package,
        interpolate_annual=True,
        ei_version="3.11",
    )

    # Optional: widen interpolation bounds for endpoint duplication
    # trails = Trails(
    #     package,
    #     interpolate_annual=True,
    #     interpolation_start_year_offset=-20,
    #     interpolation_end_year_offset=20,
    #     ei_version="3.11",
    # )

After initialization, you can access scenario labels and metadata:

.. code-block:: python

    print(trails.scenario_labels)
    activity_indices = next(iter(trails.activity_indices.values()))
    print(list(activity_indices.items())[:5])

Interpolation boundary behavior:

* Between provided inventory years, TRAILS uses linear interpolation.
* Before the earliest year and after the latest year, TRAILS duplicates the
  nearest endpoint year according to
  ``interpolation_start_year_offset`` / ``interpolation_end_year_offset``.

Selecting activities
--------------------

Activities are referenced by integer indices from the metadata. A typical
workflow is to select an activity and store the index for repeated calls:

.. code-block:: python

    activity_indices = next(iter(trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))

Running temporal LCA
--------------------

The primary workflow follows bw2calc: build the routing graph, calculate one
temporal inventory with ``lci()``, and characterize it with one or more
``lcia()`` calls.

.. code-block:: python

    # Run temporal routing (builds the traversal graph).
    # By default this uses adaptive routing with a relative cutoff of 1e-4.
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        amount=1.0,
        min_amount=1e-18,
        show_progress=True,
        debug=False,
        attribute_to_roots=True,
        adaptive_methods=[method],
    )

    # Build and retain the uncharacterized temporal inventory
    trails.lci(
        solver_mode="iterative",
        iterative_rtol=1e-3,
    )

    # Characterize without rerunning the linear systems
    scores = trails.lcia(methods=[method])

Constructor ``methods`` are optional defaults for ``lcia()``. Call-level
methods override them. Adaptive routing methods are configured independently
with ``adaptive_methods``.

For large root-attributed direct or iterative calculations,
``inventory_backend="auto"`` selects factorized storage when the predicted
eager allocation exceeds ``inventory_memory_budget`` (256 MiB by default).
Small calculations remain eager COO; other large calculations dynamically
promote to disk-backed sparse Dask blocks. The public results remain
``xarray.DataArray`` objects with the same dimensions and coordinates. TRAILS'
FaIR, EDGES, and plotting helpers process the lazy storage sequentially.

Use ``inventory_backend="chunked"`` to force this representation or
``inventory_backend="coo"`` to request the legacy eager form; an unsafe eager
request raises before allocation. ``inventory_store`` selects a user-managed
backing directory. Otherwise, call ``trails.close()`` to remove the managed
temporary store. Guarded ``materialize_inventory()`` and
``materialize_characterized_inventory()`` helpers are available when an eager
COO is genuinely required.

The automatically selected factorized backend writes annual activity-by-root
supply matrices and keeps ordinary biosphere multiplication lazy. Ported
temporal exchanges are stored as compact offset/weight kernels and applied
lazily as well. Only matrix-sourced temporal exchanges and direct biosphere
corrections use the explicit chunked sparse store. Set
``inventory_backend="factorized"`` to force this backend regardless of the
size estimate.

Adaptive score-potential routing
--------------------------------

For deep temporalized systems, a fixed ``max_depth`` can expand many branches
whose eventual contribution is negligible. Adaptive routing lets TRAILS use
static LCIA activity scores as a screening estimate before deciding whether a
child branch should be routed explicitly. This is the default routing mode:
when ``max_depth`` is omitted, TRAILS uses ``max_depth=None`` and
``adaptive_relative_score_cutoff=1e-4``.

.. code-block:: python

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        amount=1.0,
        max_depth=None,
        min_amount=1e-18,
        adaptive_methods=[method],
        adaptive_relative_score_cutoff=1e-4,
        adaptive_min_depth=1,
    )
    trails.lci()
    trails.lcia(methods=[method])

The relative cutoff is multiplied by the functional unit's static score
potential. In the example above, branches are stopped once their estimated
static potential is at most ``1e-4`` of the functional unit potential.

Routing modes
~~~~~~~~~~~~~

There are four common routing configurations:

**1. Default adaptive routing**
    Use this for normal regular-LCIA workflows. ``temporal_routing()`` defaults
    to ``max_depth=None`` and ``adaptive_relative_score_cutoff=1e-4``.

    .. code-block:: python

        trails.temporal_routing(
            start_year=2030,
            start_act_idx=start_act_idx,
        )

**2. Adaptive routing with another relative cutoff**
    Use a smaller value for stricter routing and a larger value for looser
    routing.

    .. code-block:: python

        trails.temporal_routing(
            start_year=2030,
            start_act_idx=start_act_idx,
            adaptive_relative_score_cutoff=1e-5,
        )

**3. Adaptive routing with a hard depth cap**
    Use this when branches should be pruned by score potential but never routed
    beyond a fixed graph depth.

    .. code-block:: python

        trails.temporal_routing(
            start_year=2030,
            start_act_idx=start_act_idx,
            max_depth=5,
            adaptive_relative_score_cutoff=1e-4,
        )

**4. Fixed-depth routing**
    Use this when you want the older depth-based behavior. Passing an integer
    ``max_depth`` without an adaptive relative cutoff disables adaptive pruning.

    .. code-block:: python

        trails.temporal_routing(
            start_year=2030,
            start_act_idx=start_act_idx,
            max_depth=3,
        )

Important behavior:

* Adaptive pruning only changes graph expansion. Stopped branches are stored as
  frontier demands and still enter the year-wise matrix solve in ``lci()``.
* Passing an integer ``max_depth`` without an adaptive cutoff selects
  fixed-depth routing. You can also combine an adaptive relative cutoff with a
  finite ``max_depth`` to keep a hard cap.
* Adaptive routing uses explicit regular ``adaptive_methods``. EDGES methods
  cannot currently be used for adaptive screening.
* Static activity scores are cached by default. Set
  ``adaptive_use_cache=False`` for tests or diagnostics that should recompute
  the screening intensities.
* For multi-method screening, TRAILS uses the maximum absolute potential across
  the selected methods, so a branch is retained if it is relevant for any
  screening indicator.

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
In the main temporal workflow, ``temporal_routing`` and ``lci`` use temporal
exchange distributions by default:

.. code-block:: python

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        min_amount=1e-18,
        adaptive_methods=[method],
    )
    trails.lci()
    trails.lcia(methods=[method])

For a non-temporal baseline in a single year, use ``trails.static_lca(...)``.

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

Use ``plot_adaptive_sankey`` to inspect the explicit graph created by adaptive
``temporal_routing``:

.. code-block:: python

    from trails import plot_adaptive_sankey

    fig = plot_adaptive_sankey(
        trails,
        method=method,
        branch_visual_cutoff=0.001,
        max_sankey_links=0,
        output_path="adaptive_sankey.html",
    )
    fig.show()

The Sankey links represent routed graph edges weighted by the child node's
adaptive score potential. Nodes are placed by routing depth and year, labels
are shown on hover, and optional right/bottom density panels summarize impact
density over time and depth. Matrix-solved frontier demands remain part of the
final LCA results, but only explicitly routed edges are displayed in the
Sankey.

The same helper is available on a ``Trails`` instance:

.. code-block:: python

    fig = trails.plot_adaptive_sankey(method=method)

Interpreting outputs
--------------------

Impact time series are available from ``trails.scores`` after ``lcia()``.
``trails.inventory`` is available after ``lci()``, and regular LCIA also exposes
``trails.characterized_inventory``.

Troubleshooting and diagnostics
-------------------------------

* If a requested year is not available in the data package, TRAILS will snap to
  the nearest available year and emit a warning.
* Set ``debug=True`` to enable detailed logging and retain additional
  diagnostics on the Trails instance.
* Advanced performance tuning:
  - For FaIR perturbation parallelism, use
    ``run_fair_delta_rf(..., per_species_workers=<int>)``.
  - For TD biosphere accumulation buffering, set
    ``trails._bio_inventory_flush_nnz`` before calling ``lca``.
    Default is ``2_000_000``. This is an advanced/private knob.


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

    # Ensure an inventory with root attribution is available
    trails.lci()

    rf = run_fair_delta_rf(
        trails,
        scenario="REMIND|SSP2-PkBudg650",
        # defaults shown explicitly:
        per_species_runs=True,
        per_species_workers=None,  # auto: min(4, cpu_count, n_work_items)
    )

The outputs are stored on:

* ``trails.instant_radiative_forcing`` with dims ``(quantile, year, flow, root activity)``
* ``trails.delta_temperature`` with dims ``(quantile, year, flow, root activity)``

Notes:

* ``run_fair_delta_rf`` requires ``trails.inventory`` with a
  ``root activity`` dimension. Run ``lci()`` first.
* ``scenario`` must match a scenario label present in the emissions CSV used by
  ``run_fair_delta_rf``.
* If ``config_name`` and ``config_names`` are omitted, TRAILS evaluates all
  available FaIR configurations and stores quantiles across the ensemble.

Visualization helpers default to the 50th quantile:

.. code-block:: python

    from trails import plot_rf, plot_temp

    plot_rf(trails, year_range=(2000, 2100))
    plot_temp(trails, year_range=(2000, 2100))
