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

    # Run a temporal LCA
    results = lca(
        trails=trails,
        start_year=2030,
        start_act_idx=start_act_idx,
        methods=[method],
        max_depth=2,
    )

    # Plot temporal impact scores
    fig = plot_temporal_scores(results, trails, method_label=method)
    fig.show()

What you get
------------

The returned ``results`` dictionary contains:

* ``results_by_solve_year``: diagnostics for each solved year
* ``results_by_impact_year``: impact scores aggregated by impact year (ready for
  plotting or further analysis)
