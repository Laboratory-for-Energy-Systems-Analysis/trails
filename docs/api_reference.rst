API Reference
=============

This section documents the public API of TRAILS.

Core classes
------------

.. autoclass:: trails.Trails
   :members:
   :undoc-members:
   :show-inheritance:

Core workflows
--------------

.. autofunction:: trails.lca

.. autofunction:: trails.lca_static_simple

LCIA utilities
--------------

.. autofunction:: trails.get_lcia_method_names

.. autofunction:: trails.lcia.get_lcia_methods

.. autofunction:: trails.lcia.fill_characterization_factors_matrices

Temporal distributions
----------------------

.. autoclass:: trails.temporal_distributions.TemporalExchange
   :members:
   :undoc-members:

.. autoclass:: trails.temporal_distributions.TemporalDistribution
   :members:
   :undoc-members:

Plotting
--------

.. autofunction:: trails.plot_temporal_scores
