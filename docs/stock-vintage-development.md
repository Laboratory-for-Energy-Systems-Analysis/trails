# Annual stock-vintage profiles (experimental)

This opt-in implementation is being developed on `feat/stock-vintage-pilots`.
It is not a reviewed stock dataset or a change to public defaults. The scientific
plan lives in the premise planning worktree under
`docs/development/stock-asset-temporal/`.

## Purpose and allocation convention

For a service in calendar year `t`, split the existing capital exchange amount
`q(t)` over the manufacture/commissioning years of the assets providing that
service. Each pulse is `q(t) * w(c, t)`, where the nonnegative weights sum to one.
The profile does not introduce a second lifetime denominator or change the
exchange total. `common_amortisation` is the only implemented allocation basis.
Stock shares equal service shares only when utilisation and capacity per asset
are equal; empirical curation must state these assumptions.

Calendar cohort years are independent of inventory anchors. A request for 2021
can use interpolated 2020/2022 inventories while resolving a separate annual
2021 stock profile. Changing the matrix template must not move a cohort's
manufacture date. Historical technology outside the inventory range still uses
the existing inventory boundary behaviour; precise cohort dates do not create
historical LCI evidence.

For these packages, the routing graph and final demand vectors retain calendar
years outside the matrix horizon. Matrix selection alone uses the nearest
available scenario. Inventory and score storage cover the routed dates and
biosphere support, including successive component shifts. Legacy packages keep
their previous year-mapping behaviour.

## Consumer and resource contract

Load an explicitly corrected package with:

```python
from datapackage import Package
from trails import Trails

model = Trails(Package("datapackage.json"), stock_vintage=True)
```

New readers reject a declared stock resource unless the flag is enabled. Legacy
packages use the legacy path, including when the flag is enabled. Older TRAILS
versions do not know this resource and can ignore it: they are **unsupported for
corrected packages**, even if they can open the matrix files. A released minimum
consumer version has not yet been assigned.

The package descriptor declares a JSON resource and its SHA-256:

```json
{
  "stock_vintage": {
    "schema_version": 1,
    "resource": "stock-vintage",
    "sha256": "<hash of the exact resource bytes>"
  }
}
```

The JSON resource contains `schema_version: 1`, a `profiles` list and a `bindings`
list. Each profile includes:

- `id`, `event_role`, `allocation_basis: "common_amortisation"`;
- `asset_unit` (supplier unit), `service_unit` (caller unit), and `provenance`;
- annual `years` records containing integer `service_year`, integer calendar
  `event_years`, and corresponding `weights` summing to one.

Supported roles are `existing_asset_service`, `new_asset_component`,
`replacement_component`, and `lifecycle_service`. These are explicit assignments,
not automatic classifications. Only the lifecycle role permits events after the
caller year. Binding entries associate `profile_id` with exact `caller` and
`supplier` identities: `name`, `reference product`, `unit`, `location`. Matrix
indices are resolved from the package and are never maintained in source data.

Missing service years, conflicting bindings, unresolved identities, malformed
weights, future manufacture dates, unit mismatches and changed resource bytes
are errors. Missing years do not reuse the nearest profile. Lookup emits a type
6 discrete temporal exchange with `amount_source="port"`.

The premise producer also requires `inventory_context` (source database/version,
system model, model and pathway), consecutive annual records and actual matched
inventory exchanges. Parser validation is not scientific approval: provenance,
source rights, stock balance, service assumptions and lifecycle boundaries need
separate evidence.

## Caches and foreground overrides

The JSON resource is reloaded and hashed even when inventory interpolation is
loaded from cache. Temporal expansion caches use the requested calendar year for
packages with annual profiles. Existing matrix caches retain their own identity.

Importing a foreground activity replaces its package stock bindings for the
affected templates. Its explicit temporal metadata then takes precedence.
For stock-vintage packages, `cache_import=True` does not save modified matrices
over the base package cache: removed bindings are instance state and cannot be
restored with the original package descriptor. Reapply the foreground import to
each new model instance. Legacy import caching is unchanged.

The integration test also covers two importer corrections: index zero is a
valid supplier, and replacing all anchor years replaces intermediate annual
rows instead of adding the new exchanges to the old ones.

## Validation and unresolved release work

The synthetic tests cover annual lookup, warm caches, foreground overrides,
signed amount conservation, explicit opt-in and failure cases. A routed static
limiting case compares summed biosphere results with an independent sparse
matrix solve. The producer tests exercise compact/legacy stores, checkpoint
reloads, matrix writing, global index reordering and package/resource hashes.

`tests/data/stock_lifecycle` contains a synthetic fixture made by the actual
premise exporter. Its manufacturing years are 2005/2015/2021/2022 and retirement
years 2025/2035/2041/2042, with inventory anchors only in 2020/2022. Tests check
the exact biosphere calendar and a negative recycling credit in 32 cases across
signed functional units, repeated/reversed years, root attribution, annual
interpolation and direct/Brightway solvers. A separate cascade test checks that
direct biosphere flows and successive component offsets survive outside the
matrix horizon. These tests caught both graph-year clamping and lost frontier
demands caused by clamping before storage was initialised.

Premise's new `stock_cohorts` and `stock_lifecycle` modules provide conditional
survival and scoped separation of embedded disposal. The companion document
`docs/development/stock-asset-temporal/lifecycle-and-cohort-method.md` specifies
end-year stocks, retirement conditional on service-year survival, deterministic
coefficient preservation and the condition for using lifecycle marginals:
lifted lifecycle quantities per capital unit must be invariant across anchors.
Variable quantities require cohort-specific attribution. The implementation
does not claim to preserve Monte Carlo correlations.

`dev/stock_vintage/check_lifecycle_inventory.py` verifies the car-chassis and
lorry rewrites against the complete local ecoinvent 3.12 cut-off matrix. All
3,297 biosphere flows agree within numerical precision (maximum absolute errors
`2.22e-15` and `2.22e-16`). This establishes static structural equivalence, not
empirical lifetime validity or a completed temporal pilot. The check consistently
uses processed matrix coefficients rather than mixing dataset metadata with
their lower-precision processed values.

These checks establish software behaviour only. The five empirical pilots still
need reviewed caller/supplier boundaries, cohort and service reconciliation,
real-database legacy/corrected comparisons, sensitivity analysis and lifecycle
checks. A shifted parent asset plus an unreviewed legacy component profile can
still duplicate the vintage shift. Likewise, a flat manufacture profile alone
does not retain the service-year context needed to condition future retirement
on an asset having survived to that service year. No pilot may be promoted until
its relevant boundaries and event chronology have been demonstrated.

Local helper `dev/stock_vintage/audit_inventory.py` reads an existing Brightway
database and writes only candidate identities and selected capital-exchange
metadata. Raw licensed inventories and IAM inputs remain outside version control.
