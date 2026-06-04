LCIA Methods
============

TRAILS ships built-in LCIA method sets for ecoinvent versions ``3.10`` and
``3.11``. These methods are bundled with the package data and exposed through
``trails.lcia`` utilities.

List available method names
---------------------------

.. code-block:: python

    from trails import get_lcia_method_names

    methods_311 = get_lcia_method_names(ei_version="3.11")
    methods_310 = get_lcia_method_names(ei_version="3.10")

    print(methods_311[:5])

Method names are returned as strings joined with ``" - "`` from the original
LCIA tuple labels.

Use methods in temporal LCA
---------------------------

.. code-block:: python

    from trails import Trails, lca, get_lcia_method_names

    method = get_lcia_method_names(ei_version="3.11")[0]
    trails = Trails(
        package,
        methods=[method],
        ei_version="3.11",
    )

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
    )

    lca(
        trails=trails,
        # defaults shown explicitly:
        solver_mode="iterative",
        iterative_rtol=1e-3,
    )

Call-level ``methods=...`` and ``ei_version=...`` still override constructor
defaults when needed.

Use EDGES methods
-----------------

TRAILS can also score the finalized temporal inventory with EDGES
edge-level characterization factors. This requires the optional ``edges``
package and currently supports EDGES methods whose supplier matrix is
``"biosphere"``.

.. code-block:: python

    from trails import Trails, lca, get_edges_lcia_method_names

    method = get_edges_lcia_method_names()[0]
    trails = Trails(
        package,
        edges_methods=[method],
    )

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        # EDGES methods are final-score methods only, so this example uses
        # explicit fixed-depth routing. Use Trails(..., methods=...) as well
        # if you want adaptive routing before EDGES final scoring.
        max_depth=2,
    )

    lca(
        trails=trails,
        # edges_additional_topologies=topology,  # optional
        edges_reuse_cached_cfs=True,
    )

Final scoring with ``edges_methods`` is mutually exclusive with regular
``methods`` in a single ``lca()`` call. In EDGES mode, TRAILS builds
``trails.inventory`` internally because edge-level CFs are applied after the
temporalized biosphere inventory has been assembled. Constructor
``edges_methods`` are used only for final EDGES scoring. Adaptive routing still
needs regular LCIA methods, passed with ``Trails(..., methods=...)`` or
``temporal_routing(..., adaptive_methods=...)``.

By default, ``edges_reuse_cached_cfs=True`` reuses EDGES matched CF templates
across scenario years when the supplier and consumer metadata signatures are
identical. The numeric CF values are still evaluated for each scenario year,
which keeps yearly parameterized methods such as AWARE yearly factors
year-specific while avoiding repeated exchange matching. Set
``edges_reuse_cached_cfs=False`` if an EDGES method has year-specific matching
rules or year-specific CF definitions that should be resolved independently for
each year.

Inspect LCIA flow factors
-------------------------

Use ``get_lcia_methods`` if you want the underlying mapping from biosphere flow
keys to characterization factor values:

.. code-block:: python

    from trails import get_lcia_method_names
    from trails.lcia import get_lcia_methods

    method_name = get_lcia_method_names(ei_version="3.11")[0]
    lcia_data = get_lcia_methods(methods=[method_name], ei_version="3.11")
    flow_to_cf = lcia_data[method_name]

    # Key format: (name, compartment, subcompartment)
    print(next(iter(flow_to_cf.items())))

Matching behavior
-----------------

TRAILS matches LCIA factors to biosphere flows using exact tuple keys:

* ``(name, compartment, subcompartment)``

If a method exchange has only one category level, TRAILS uses
``"unspecified"`` as the subcompartment placeholder.

Version selection
-----------------

Set ``ei_version`` consistently in:

* ``get_lcia_method_names``
* ``Trails(..., ei_version=...)`` or ``trails.lca(..., ei_version=...)``
* ``trails.lcia.get_lcia_methods``

so method names and factors come from the same bundled dataset.
