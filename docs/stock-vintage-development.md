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

Stock packages also convert signed product requirements into activity scaling
using the supplier's signed reference production. For example, a waste input of
`-0.5` into a supplier with reference output `-2` requires `+0.25` units of
treatment activity. A credit to a positive-reference supplier retains its
negative sign. This corrects an error exposed by the real CCGT cooling-tower
residue market; the old absolute denominator reversed its activity and transport
burdens. Legacy packages retain the historical convention for reproducibility.
The root functional-unit convention remains the same.

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

The new signed-production regression has 64 cases, including 16 legacy controls.
It checks non-unit negative reference outputs, signed functional units, port and
matrix temporal amounts, both solver paths, root attribution and two routing
depths. An independent analytical result is `0.66` for the corrected fixture;
the legacy discrepancies are retained explicitly. Additional focused engine,
inventory-storage, importer, lifecycle and distribution suites pass.

## PV parent and component evidence

The companion premise branch now constructs a bounded PV population of 363
complete US-WECC EIA plants, preserving observed AC/DC capacity and annual net
generation separately. Its fixed-crystalline classification is an explicit
proxy for the source multi-Si ground-mounted inventory. Both pinned local REMIND
scenarios generate eight cases over 2022–2030; these are regional reconstructions,
not native IAM vintage outputs.

Regular replacement generates the inverter actually serving the requested year.
Manufacture is no later than service; disposal is the earlier of its next renewal
or conditional parent retirement. The producer retains the already amortised
quantity. Analytical integration over thirty equal annual services recovers one
initial and one replacement inverter for a fifteen-year interval. Component
scope tests also prevent lifted manufacture from reintroducing embedded disposal.

`dev/stock_vintage/report_pv_sensitivity.py` writes a 144-row local CSV, input-hash
manifest and PNG/SVG/PDF figure from the two cohort outputs. The figure explicitly
distinguishes observed parent dates from assumed inverter dates. IAM-derived
figures and tables remain restricted local evidence.

The actual PV export now separates 47 waste ports, inverter and panel manufacture,
and lifetime land occupation. Fifty lifted ports plus twelve copied capital nodes
and one service caller add 63 activities. The complete matrix proof passes with
relative A/B residuals below `2.62e-13`/`1.83e-15`; all 510 signed annual port
checks and 58 internal zero shifts pass. `check_pv_pilot.py` reuses the matrix and
streamed LCI verifier, adding a direct occupation check at the requested service
year. Its analytical tests include both functional-unit signs and caller losses.

The producer documents expected service-year panel maintenance, installation
handling loss, and explicit factory/end-of-life polymer and mounting-packaging
timing endpoints. Full temporal LCI, cache/no-interpolation, legacy and endpoint
comparisons remain release work. There are 89 passing focused profile/verifier
tests; the three PV requirements are not yet claimed complete.

## Real CCGT evidence

Premise commit `e2989248` provides the EIA block audit, exact REMIND leaf-variable
reader, seven annual reconstruction sensitivities and the scoped inventory
exporter. Its companion `ccgt-pilot-method.md` documents the 91-block US-WECC
boundary, observed 2022 output weights, USA relative-trajectory proxy, initial
old-stock extension, future capacity closure, service weighting and unchanged
common amortisation. Both IAM pathways have nine annual records per sensitivity;
the local table contains 126 scenario/case/year records.

`dev/stock_vintage/extract_local_inventory.py` reads the existing complete
26,533-activity ecoinvent 3.12 cut-off system into a restricted local input.
The actual premise exporter creates legacy and corrected packages at
2022/2025/2030, with three scoped CCGT activities. Technology remains constant
across these anchors to isolate timing and conservation. This is not an
IAM-transformed background LCA or historical technology reconstruction.

`dev/stock_vintage/check_ccgt_pilot.py` checks actual annual capital dates,
internal zero shifts, repeated/reversed service requests and signed quantities.
Real LCI checks pass for service year 2023 with annual interpolation and a warm
cache, and 2030 with interpolation disabled. Quantities are proved independently
through the graph identity `A x_expanded + f_frontier = f_service`; componentwise
relative residuals are below `1.2e-16` in the corrected runs.

All 9,847 exported biosphere rows agree with an independently assembled
year-wise static calculation (maximum absolute discrepancy `2.22e-16` across
the tested signed functional units). This comparison uses the same numerical
backend while assembling demands and biosphere multiplication independently.
A separate one-shot SuperLU reference is retained in the report: its differences
are about `5.74e-9` per positive FU. The complete matrix is poorly conditioned,
so comparing different sparse decompositions alone would mix numerical error
with changed inventory quantities. The analytical regression also checks the
independent verifier, and the graph identity prevents a wrong route from passing
merely because both numerical accumulations use it.

The complete inventory uses factorized storage. Its bounded native flow reducer
avoids materialising the activity/flow/year/root tensor. Generic Xarray
skip-NaN reduction attempted an impractically large sparse broadcast during this
check; the verifier uses the same native reduction path as TRAILS' FaIR adapter.

For a 2030 kWh, the capital-attributed fossil-CO2 total is approximately
`3.83763e-4 kg` in both variants. Its corrected mean age is 23.4754 years;
the legacy no-interpolation router confines it to 2022/2025/2030 (mean age
7.8052 years). The latter includes calendar clamping and must not be interpreted
as the raw legacy distribution's mean: the separately checked, unclamped
construction pulses have a mean age of 16.9338 years. The legacy waste-market sign error is
reported separately from this timing comparison.

`dev/stock_vintage/report_ccgt_sensitivity.py` produces a local CSV and PNG/SVG/PDF
figure from the two generated cohort reports. It separates observed capacity
and output weighting and shows the lifetime, old-stock-tail and future-additions
sensitivities. Raw inventory, IAM-derived series and generated packages stay
outside version control.

The CCGT bounded reference case now has real evidence for the three requirements;
the wider pilot campaign and final promotion audit remain unfinished. Passenger
cars, trucks, PV and water still need their empirical temporal comparisons and
remaining lifecycle/component checks. A shifted parent asset plus an unreviewed legacy component profile can
still duplicate the vintage shift. Likewise, a flat manufacture profile alone
does not retain the service-year context needed to condition future retirement
on an asset having survived to that service year. No pilot may be promoted until
its relevant boundaries and event chronology have been demonstrated.

Local helper `dev/stock_vintage/audit_inventory.py` reads an existing Brightway
database and writes only candidate identities and selected capital-exchange
metadata. Raw licensed inventories and IAM inputs remain outside version control.

## Water-network evidence in progress

Premise commit `ef2b1720` adds the Quebec water pilot with separate pipe and
storage cohorts. Its `water-pilot-method.md` records the source bins, missing
dates, 13 sensitivities, provisional survival and replacement-only future stock
scenario. Water storage has 25.56% unknown dates, so its sensitivity is material.
The primary scenario does not claim an IAM water forecast.

The actual export adds 36 scoped activities to the complete source inventory.
It separates nine tank-material inputs and 24 disposal ports from pipe/network
construction, preserving their original amounts and geographic suppliers.
The caller's own-product loss remains 0.191, leaving 0.809 kg net output.
The full matrix rewrite proof passes, with maximum relative technosphere
residual `5.18e-14` and biosphere residual `2.60e-17`. There are 340 signed annual
pulse checks and 34 internal zero-shift bindings. The old raw network-construction
distribution has mean age 35 years, independently of legacy calendar clamping.

`check_water_pilot.py` records cache requests separately from confirmed cache
loads. Full corrected temporal LCI checks remain in progress at this checkpoint.
The first complete solve returned, but its all-flow activity reduction was
killed during verification: even a reduced sparse cube can require enormous
intermediate arrays across many retirement years and roots. The verifier now
streams native inventory entries in bounded windows, accumulating every flow
total and only four dated fossil-CO2 summaries. Its independent reference also
keeps year/root demands sparse until each solve. Six analytical checks cover
the operator proof and streaming with both storage backends and both signs;
85 focused profile/verifier tests pass. The model and exported profiles did not
change as a result of this verifier fix.

`report_water_sensitivity.py` writes 234 scenario/asset/year records and a
PNG/SVG/PDF figure. Mean ages are labelled as service weighted, since one case
tests age-dependent utilisation. The generated source and package artifacts
remain outside the repository.
