Quickstart
==========

This quickstart shows how to load a Frictionless data package, run a temporal
LCA, and plot impact scores over time.

Install
-------

.. code-block:: bash

    pip install trails

You can also install from conda:

.. code-block:: bash

    conda install romainsacchi::trails

Solver note: TRAILS relies on fast sparse solvers for large LCAs. Install the
recommended solver for your platform (see the README) to avoid slow SciPy
fallbacks.

Run a temporal LCA
------------------

.. code-block:: python

    from datapackage import Package

    from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores

    # Load a Frictionless data package exported by premise (or compatible tooling)
    package = Package("path/to/datapackage.json")

    # Choose an LCIA method bundled with TRAILS
    method = get_lcia_method_names(ei_version="3.11")[0]

    # Initialize TRAILS (annual interpolation is optional).
    # Default interpolation bounds are [min_year-1, max_year+1], where
    # out-of-range years duplicate the nearest endpoint inventory.
    trails = Trails(
        package,
        interpolate_annual=True,
        methods=[method],
        ei_version="3.11",
    )

    # Optional: widen interpolation bounds (e.g., +/-20 years)
    # trails = Trails(
    #     package,
    #     interpolate_annual=True,
    #     interpolation_start_year_offset=-20,
    #     interpolation_end_year_offset=20,
    #     methods=[method],
    #     ei_version="3.11",
    # )

    # Pick an activity index from the metadata
    activity_indices = next(iter(trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))

    # Run temporal routing (builds the traversal graph).
    # By default this uses adaptive routing with a relative cutoff of 1e-4.
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        min_amount=1e-18,
    )

    # Run temporal LCA (stores scores on trails.scores)
    lca(
        trails=trails,
        # defaults shown explicitly:
        solver_mode="iterative",
        iterative_rtol=1e-3,
    )

    # Plot temporal impact scores
    fig = plot_temporal_scores(trails, method_label=method)
    fig.show()

Adaptive routing is the default. To make the default criterion explicit:

.. code-block:: python

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        max_depth=None,
        adaptive_relative_score_cutoff=1e-4,
    )

Branches below the cutoff remain in the matrix-solved frontier, so they are
still included in the final LCA.

To run fixed-depth routing instead, pass an integer ``max_depth`` and omit
adaptive cutoffs.

What you get
------------

Temporal LCA results are stored on the Trails instance. Use ``trails.scores``
for impact scores (when ``compute_score=True``). If you run
``lca(..., store_inventory=True)``, TRAILS also stores ``trails.inventory``;
with ``compute_score=True`` and ``store_inventory=True``, it also stores
``trails.characterized_inventory``.


Importing Excel inventories
---------------------------

You can import user-provided inventories from Excel using ``bw2io``. When you
omit ``year`` and ``scenario_label``, the exchanges are applied to all template
years and interpolated across annual years.

Sign conventions for Excel imports:
- Technosphere exchanges are sign-flipped on import (positive becomes negative, and vice versa).
- Production and biosphere exchanges are stored as-is.

.. code-block:: python

    from trails import Trails

    trails = Trails(package)
    trails.import_excel_inventory("path/to/inventory.xlsx")

    # Target a single scenario slice instead
    trails.import_excel_inventory("path/to/inventory.xlsx", year=2020)

Excel column meanings (from the exchanges table):

* ``name``: exchange name.
* ``amount``: base exchange amount (used when no year-specific columns exist).
* **Year-specific columns** (e.g., ``2020``, ``2040``): optional numeric columns
  that set year-specific amounts. TRAILS writes these values into the matching
  years and linearly interpolates between them; years outside the provided range
  are clamped to the nearest endpoint.
* ``location``: exchange location (required for technosphere/production).
* ``categories``: biosphere categories (tuple or ``compartment/subcompartment``).
* ``unit``: exchange unit.
* ``type``: ``production``, ``technosphere``, or ``biosphere``.
* ``reference product``: required for technosphere/production exchanges.
* ``temporal_distribution``: temporal distribution code.
* ``temporal_loc`` / ``temporal_scale``: distribution parameters.
* ``temporal_min`` / ``temporal_max``: integer offset bounds.
* ``temporal_offsets`` / ``temporal_weights``: JSON lists used for
  ``temporal_distribution=6`` (discrete empirical), e.g.
  ``[0, 5, 12]`` and ``[0.5, 0.3, 0.2]``.
* ``temporal_amount_source``: ``port`` (ported amount) or ``matrix`` (use matrix values).
* ``comment``: optional notes.


FaIR_ radiative forcing
-----------------------

After running a temporal LCA, you can translate the inventory into radiative
forcing and temperature anomalies using the FaIR_ climate model. This uses a
baseline IAMC scenario and per-species perturbations derived from the Trails
inventory. Outputs are aggregated across all FaIR_ configurations and stored as
quantiles (2.5, 25, 50, 75, 97.5).

.. code-block:: python

    from trails import lca
    from trails.fair_rf import run_fair_delta_rf

    # Ensure an inventory with root attribution is available
    lca(
        trails=trails,
        store_inventory=True,
    )

    rf = run_fair_delta_rf(
        trails,
        scenario="REMIND|SSP2-PkBudg650",
        # defaults shown explicitly:
        per_species_runs=True,
        per_species_workers=None,  # auto: min(4, cpu_count, n_work_items)
    )

The resulting outputs are stored on the Trails instance:

* ``trails.instant_radiative_forcing`` with dims ``(quantile, year, flow, root activity)``
* ``trails.delta_temperature`` with dims ``(quantile, year, flow, root activity)``

Notes:

* ``run_fair_delta_rf`` requires ``trails.inventory`` with a
  ``root activity`` dimension. Run ``lca(..., store_inventory=True)`` first.
* ``scenario`` must match a scenario label present in the emissions CSV used by
  ``run_fair_delta_rf``.
* If ``config_name`` and ``config_names`` are omitted, TRAILS evaluates all
  available FaIR configurations and stores quantiles across the ensemble.

You can visualize these with the built-in plotting helpers (defaults to the 50th
quantile):

.. _FaIR: https://github.com/OMS-NetZero/FAIR

.. code-block:: python

    from trails import plot_rf, plot_temp

    plot_rf(trails, year_range=(2000, 2100))
    plot_temp(trails, year_range=(2000, 2100))
