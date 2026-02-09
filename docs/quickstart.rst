Quickstart
==========

This quickstart shows how to load a Frictionless data package, run a temporal
LCA, and plot impact scores over time.

Install
-------

.. code-block:: bash

    pip install trails

Run a temporal LCA
------------------

.. code-block:: python

    from datapackage import Package

    from trails import Trails, lca, get_lcia_method_names, plot_temporal_scores

    # Load a Frictionless data package exported by premise (or compatible tooling)
    package = Package("path/to/datapackage.json")

    # Initialize the TRAILS wrapper (annual interpolation is optional)
    trails = Trails(package, interpolate_annual=True)

    # Pick an activity index from the metadata
    activity_indices = next(iter(trails.activity_indices.values()))
    start_act_idx = next(iter(activity_indices.keys()))

    # Choose an LCIA method bundled with TRAILS
    method = get_lcia_method_names(ei_version="3.11")[0]

    # Run temporal routing (builds the traversal graph)
    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        max_depth=2,
        min_amount=1e-18,
    )

    # Run temporal LCA (stores scores on trails.scores)
    lca(
        trails=trails,
        methods=[method],
    )

    # Plot temporal impact scores
    fig = plot_temporal_scores(trails, method_label=method)
    fig.show()

What you get
------------

Temporal LCA results are stored on the Trails instance. Use ``trails.scores``
for impact scores (when compute_score=True) and ``trails.inventory`` or
``trails.characterized_inventory`` for time-resolved inventories.


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


FaIR radiative forcing
----------------------

After running a temporal LCA, you can translate the inventory into radiative
forcing using the FaIR climate model. This uses a baseline IAMC scenario and
per-species perturbations derived from the Trails inventory.

.. code-block:: python

    from trails.fair_rf import run_fair_delta_rf

    rf = run_fair_delta_rf(
        trails,
        scenario="high-extension",
    )

The resulting ``rf`` is stored on ``trails.instant_radiative_forcing`` with
coords (year, flow, root activity).
