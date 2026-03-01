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

    from trails import lca, get_lcia_method_names

    method = get_lcia_method_names(ei_version="3.11")[0]

    trails.temporal_routing(
        start_year=2030,
        start_act_idx=start_act_idx,
        max_depth=2,
    )

    lca(
        trails=trails,
        methods=[method],
        ei_version="3.11",
    )

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
* ``trails.lca(..., ei_version=...)``
* ``trails.lcia.get_lcia_methods``

so method names and factors come from the same bundled dataset.
