Nomenclature
============

This page summarizes key naming conventions used by TRAILS in data packages,
inventory outputs, and LCIA matching.

Scenario labels
---------------

Scenario labels are stored as strings (for example ``"2020"``, ``"2050"``) and
interpreted as calendar years. Internally, TRAILS maps requested years to the
nearest available scenario year when needed.

Activity metadata
-----------------

Activity metadata is stored per scenario label under ``trails.activity_indices``.
Each activity entry is keyed by integer activity id and typically includes:

* ``name``
* ``reference product``
* ``location``
* ``unit``

Biosphere flow keys
-------------------

Biosphere metadata is stored under ``trails.biosphere_indices``. LCIA matching
and several downstream workflows use biosphere flow tuple keys in this form:

* ``(name, compartment, subcompartment)``

Examples:

* ``("Carbon dioxide, fossil", "air", "urban air close to ground")``
* ``("Water, river", "natural resource", "in water")``

When a subcompartment is unavailable in source method data, TRAILS uses the
placeholder ``"unspecified"``.

Temporal exchange fields
------------------------

Temporal exchange metadata (when present) uses these fields:

* ``temporal_distribution``
* ``temporal_loc``
* ``temporal_scale``
* ``temporal_min`` / ``temporal_max``
* ``temporal_offsets`` / ``temporal_weights`` (for discrete empirical pulses)
* ``temporal_amount_source`` (``port`` or ``matrix``)

These fields map to ``TemporalExchange`` objects during routing and inventory
accumulation.

FaIR mapping keys
-----------------

FaIR species mapping uses biosphere flow keys compatible with the tuple format
above. Mapping can also fall back to flow name-only keys when tuple matching is
not available in the mapping file.
