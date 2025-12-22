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

    results = lca(
        trails=trails,
        start_year=2030,
        start_act_idx=start_act_idx,
        methods=[method],
        amount=1.0,
        max_depth=3,
        min_amount=1e-18,
        show_progress=True,
        debug=False,
        return_provenance=False,
        use_temporal_distributions=True,
    )

The ``results`` structure provides:

* ``results_by_solve_year``: demand vectors and diagnostics per solved year
* ``results_by_impact_year``: scores and attribution by impact year

Temporal distributions
----------------------

Temporal distributions control how exchanges are spread across impact years.
TRAILS can use the distribution data included in the package, or collapse the
effects into scalar multipliers for a static approximation:

.. code-block:: python

    # Use full temporal distributions (default)
    results = lca(
        trails=trails,
        start_year=2030,
        start_act_idx=start_act_idx,
        methods=[method],
        use_temporal_distributions=True,
    )

    # Collapse temporal distributions for a faster static view
    results_static = lca(
        trails=trails,
        start_year=2030,
        start_act_idx=start_act_idx,
        methods=[method],
        use_temporal_distributions=False,
    )

Plotting results
----------------

Use ``plot_temporal_scores`` to visualize the time series and contribution by
first-level suppliers:

.. code-block:: python

    from trails import plot_temporal_scores

    fig = plot_temporal_scores(
        results,
        trails,
        method_label=method,
        stacked=True,
        show_cumulative_axis=True,
    )
    fig.show()

Interpreting outputs
--------------------

``results_by_impact_year`` is structured as:

.. code-block:: python

    {
        2030: {
            "scores": 1.23,
            "scores_by_first_level_child": {
                42: 0.8,
                77: 0.4,
            },
        },
        2031: {...},
    }

Use ``scores`` for the total impact in a year, and
``scores_by_first_level_child`` for attribution to first-level suppliers.

Troubleshooting and diagnostics
-------------------------------

* If a requested year is not available in the data package, TRAILS will snap to
  the nearest available year and emit a warning.
* Set ``debug=True`` to retain more diagnostic information in
  ``results_by_solve_year`` and enable detailed logging.
