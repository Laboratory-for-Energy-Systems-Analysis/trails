LCIA Methods
============

TRAILS ships built-in LCIA method sets for ecoinvent versions ``3.10``,
``3.11``, and ``3.12``. These methods are bundled with the package data and
exposed through ``trails.lcia`` utilities. ``Trails`` uses ``3.12`` as its
default regular-LCIA and adaptive-screening version.

List available method names
---------------------------

.. code-block:: python

    from trails import get_lcia_method_names

    methods_312 = get_lcia_method_names(ei_version="3.12")
    methods_311 = get_lcia_method_names(ei_version="3.11")
    methods_310 = get_lcia_method_names(ei_version="3.10")

    print(methods_312[:5])

Method names are returned as strings joined with ``" - "`` from the original
LCIA tuple labels.

Use methods in temporal LCA
---------------------------

.. code-block:: python

    from trails import Trails, get_lcia_method_names

    method = get_lcia_method_names(ei_version="3.12")[0]
    trails = Trails(
        package,
        ei_version="3.12",
    )

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        adaptive_methods=[method],
    )

    trails.lci(
        solver_mode="iterative",
        iterative_rtol=1e-3,
    )
    scores = trails.lcia(methods=[method])

``Trails(..., methods=[...])`` remains an optional default. Methods supplied to
``lcia()`` override constructor defaults.

Use EDGES methods
-----------------

TRAILS can also score the finalized temporal inventory with EDGES
edge-level characterization factors. This requires the optional ``edges``
package and currently supports EDGES methods whose supplier matrix is
``"biosphere"``.

.. code-block:: python

    from trails import Trails, get_edges_lcia_method_names

    method = get_edges_lcia_method_names()[0]
    trails = Trails(package)

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        # EDGES methods are final-score methods only, so this example uses
        # explicit fixed-depth routing.
        max_depth=2,
    )

    trails.lci()
    scores = trails.lcia(
        methods=[method],
        method_backend="edges",
        # additional_topologies=topology,  # optional
        reuse_mappings=True,
    )

Regular and EDGES methods characterize the same finalized inventory in separate
``lcia()`` calls. Adaptive routing still needs regular LCIA methods passed via
``temporal_routing(..., adaptive_methods=...)``.

By default, ``reuse_mappings=True`` reuses EDGES matched CF templates
across temporal inventory years when the supplier and consumer metadata
signatures are identical. For every inventory year, TRAILS calls EDGES CF
evaluation with that year as ``scenario_idx``. This keeps prospective methods
such as ``("AWARE 2.0 prospective", "Country", "unspecified", "yearly")``
year-specific, including EDGES interpolation between source years, while
avoiding repeated exchange matching. The actual inventory year is passed even
when TRAILS maps that year to a nearby database scenario year for its A/B
matrices. Set
``reuse_mappings=False`` if an EDGES method has year-specific matching
rules that change which CF row matches an exchange.

Inspect LCIA flow factors
-------------------------

Use ``get_lcia_methods`` if you want the underlying mapping from biosphere flow
keys to characterization factor values:

.. code-block:: python

    from trails import get_lcia_method_names
    from trails.lcia import get_lcia_methods

    method_name = get_lcia_method_names(ei_version="3.12")[0]
    lcia_data = get_lcia_methods(methods=[method_name], ei_version="3.12")
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
* ``Trails(..., ei_version=...)`` or ``trails.lcia(..., ei_version=...)``
* ``trails.lcia.get_lcia_methods``

so method names and factors come from the same bundled dataset.

When ``ei_version`` is omitted from ``Trails(...)``, ``temporal_routing(...)``,
``lcia(...)``, and ``static_lca(...)``, the instance-level default is ``3.12``.
Pass an explicit version to the standalone method-listing/loading helpers when
you need to keep them aligned with the instance.
