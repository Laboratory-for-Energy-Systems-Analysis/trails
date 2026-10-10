# Minimum data specification for robust REMIND vintage profiles

**Checkpoint:** retain this exact-run extraction specification while further
implementation is paused. Its small first-delivery scope reflects the earlier
four-technology pilots; the subsequent producer covers eleven groups. Native
stock and service cohorts remain the highest-value request.

Prepared 10 October 2026. Companion to the
[robustness plan](remind-stock-vintage-robustness-plan.md). This is an unsent
technical specification: availability of the requested exact-run objects has
not been established, and no provider has been contacted.

## Purpose and minimum first delivery

Establish whether the locally held REMIND scenario exports can support exact
stock/service cohort reproduction. Request a small reference extract first:

- Scenarios: `SSP2-NPi2025` and `SSP2-PkBudg650`, corresponding to the exact
  local `REMIND_generic_*.mif` inputs in the `remind 3.5.2` directory.
- Transport: EUR four-wheeler BEVs and EUR heavy 40-tonne liquid-fuel trucks,
  retaining native class/powertrain leaves required to form those aggregates.
- Power: USA PV and unabated combined-cycle gas, retaining native grades and
  distinctions needed to match the published capacity and generation leaves.
- Minimum evaluation years: native states/reporting periods covering 2020–2035,
  including all historical cohorts alive in those states. Prefer the full
  available historical/future horizon for later 2050/2100 tests.
- Granularity: preserve actual native time periods. Transport annual records
  are preferable where natively available. Do not fabricate annual power
  commissioning dates merely to satisfy this request.

A compact CSV/Parquet export with a dictionary is acceptable. RDS/GDX is also
acceptable if accompanied by the variable/dimension definitions and extraction
instructions. Do not require a costly full-model rerun if the necessary
intermediate outputs are already retained.

## 1. Run and reporting manifest — required

| Field | Needed detail |
|---|---|
| Run identity | Unique run ID, scenario label, model version/commit, configuration and input-data versions/hashes. |
| Reporting lineage | remind2, EDGE-Transport and reporttransport versions/commits, harmonisation options, export date and relevant post-processing. |
| Region definition | Country-to-region map for that run; atomic regions versus EUR/EU27/World aggregates; regional changes between configurations. |
| Technology definition | Native class, powertrain, fuel, plant technology and grade IDs; mapping to published MIF leaves; overlapping aggregates. |
| Native time | Time-point/interval definitions, model step widths, integration weights, stock snapshot conventions and any construction lead time. |
| Units and basis | Vehicles, pkm, vehicle-km, tonne-km, GW, GW/year and generation basis; power AC/DC convention and net/gross/loss treatment. |
| Rights and identity | Permitted local/derived use, source-file hashes, missing-value conventions and any available export precision. |

Keep the two scenarios' manifests distinct. Matching filenames or a folder
version do not establish common run/reporting lineage.

## 2. Transport cohorts — highest priority

The archived upstream fleet routine contains objects named
`fleetVehNumbersConstrYears`, `fleetESdemandConstrYears` and
`reportEarlyRetirement`. These are **leads for extraction**, not confirmed object
names retained by the local scenario runs. Equivalent exact-run exports are
acceptable. The archived default depreciation parameters are not substitutes
for the inputs actually used in those runs.

Request these concepts on compatible dimensions:

| Required concept | Dimensions / interpretation |
|---|---|
| Surviving vehicles | Run, region, vehicle class, powertrain, evaluation year and original construction/first-service year. |
| Annual service | Same cohort dimensions, retaining native pkm/tkm and vehicle-km where represented; identify pre/post harmonisation values. |
| New entries / sales | Annual new vehicles; distinguish annual entrants from period totals or samples of annual values. |
| Natural and additional exits | Cohort/year breakdown where available; native meaning of retirement, export, service reallocation or scrappage. Mark concepts not represented. |
| Mileage and loading | Annual mileage, occupancy/payload and empty-running convention, with the actual age, technology, class, region and time dimensions used. |
| Survival / depreciation | Exact parameter table, age indexing, terminal-age rule, time dependence and early-exit algorithm. |
| Initial state | Historical cohort construction rule, any initialisation before the model start and initial stock/service totals. |

If service is available only by technology/year, state that limitation. A
technology-wide utilisation factor cannot identify service by vintage. If
retirement is not physical disposal, retain the original meaning and provide
no invented dismantling dates.

## 3. Power capacity and service — highest priority

Prefer surviving capacity contributions by original investment cohort, along
with the corresponding native stock totals. If those contributions are not
saved, request the ingredients needed to reproduce the exact capacity equation.

The previously inspected REMIND source identifies `vm_deltaCap`, `vm_cap`,
`pm_omeg`, `pm_ts` and `vm_capEarlyReti` as relevant. Confirm their names,
dimensions, units and meaning in the exact run. Do not treat the names as a
complete, verified extraction schema; request additional initial-state,
technology/grade and active/idle-capacity fields as required by that equation.

| Required concept | Detail |
|---|---|
| Surviving capacity by cohort | Region, technology, grade, evaluation period, investment/commissioning cohort and capacity contribution. |
| Gross additions | Native rate or period total, its represented interval, and the corresponding published New Cap transformation. |
| Initial/historical capacity | Every initial cohort or aggregate initial condition needed for the native stock identity. |
| Depreciation and lifetime | Remaining-capacity weights, exact indexing and whether changing parameters attach to construction or evaluation time. |
| Early exit / idle / reactivation | Distinguish operating removal, idle capacity, re-entry and physical decommissioning. Include cohort information where represented. |
| Electricity service | Native output by cohort where represented; otherwise technology/grade output with an explicit statement that vintage dispatch is unavailable. |
| Output accounting | Capacity factors, gross/net generation, storage losses and curtailment definitions; specify what is already subtracted in the published leaf. |
| Crosswalk | Native capacity and generation categories that reproduce the selected published leaves, including the PV AC/DC basis. |

Power generation may not be modelled by vintage. In that case, native capacity
can be reproduced exactly while service-weighted vintages still require an
additional, documented allocation assumption. A data request cannot create
native detail that the model does not represent.

## 4. Canonical exchange format for our adapter

Use separate stock, service, flow and metadata tables. Avoid pretending that
all rows share annual or physical-retirement semantics.

```text
cohort_stock:
  run_id, scenario, region_id, technology_id, grade_or_vehicle_class,
  evaluation_time, stock_date_convention, cohort_time, cohort_date_convention,
  state [installed/operating/idle/...], quantity, unit, source_kind

cohort_service:
  same identity and cohort keys,
  service_period_start, service_period_end, quantity, unit,
  exposure_basis, harmonisation_stage, source_kind

stock_flows:
  same identity keys, cohort_time if represented,
  flow_kind, period_start, period_end,
  value, value_kind [rate/period_total/annual_count], unit,
  exit_meaning, source_kind

metadata:
  native-to-reporting period map, variable definitions, units,
  run/input/reporting hashes, mappings, precision and missing-value rules
```

`source_kind` distinguishes native model output, a reconstruction, an observation
and a proxy. If a native cohort is a multi-year interval or an undated initial
pool, preserve that representation. Our annualisation layer may split it later
with explicit assumptions and sensitivity bounds. Zero and unavailable values
must remain distinguishable.

## 5. Acceptance checks for the first extract

1. Verify run identity and hashes; unique keys, complete declared dimensions,
   compatible units and an explicit missing-data mask.
2. Sum native cohort stocks to native stock totals for every selected period;
   sum service cohorts to the corresponding service aggregate where available.
   Check both sides of any independent harmonisation step.
3. Reproduce the native stock-transition equation using initial state,
   additions, depreciation and additional exits. Derive numerical/rounding
   tolerances from the exported precision before inspecting discrepancies.
4. Reproduce the matching MIF rows using the documented reporting transformation,
   including period integrals. Report discrepancies separately from internal
   cohort sums: the two checks validate different stages.
5. Retain historic cohorts and preserve sums under any annual split. Compare
   resulting annual profiles with the current labelled reconstruction and
   report full-distribution distances, not mean age alone.
6. State exactly what can be claimed: native capacity reproduction, native
   service reproduction, or reconstructed service. Failure to obtain one field
   must not quietly upgrade a proxy into an observation.

## Fallback if native outputs are unavailable

Obtain the exact parameter/reporting manifest plus stocks, additions, service,
available exits and initialisation information. Build a constrained cohort
reconstruction with explicit residuals and independent historical validation.
Native-model fidelity remains unverified. Continue public-data validation and
period/region tests while waiting; neither a provider response nor permission
to contact one is assumed by this specification.

## Evidence underlying the request

The implementation worktree's
`docs/development/stock-asset-temporal/iam-scenario-assessment.md` records the
local export audit and pinned model/reporting revisions. The private evidence
archive retains the inspected EDGE and reporttransport source files. Relevant
upstream identifiers are:

- [EDGE fleet routine](https://github.com/pik-piam/edgeTransport/blob/1124c9967e9644b41977c4750640629d51b5494b/R/toolCalculateFleetComposition.R).
- [REMIND capacity equations](https://github.com/remindmodel/remind/blob/01d9312ca05c8ec35c0f28f5a177050e7e4b6768/core/equations.gms).
- [remind2 investment-period conversion](https://github.com/pik-piam/remind2/blob/e414aecde1eb0ce334c2424367617d711374ea09/R/modifyInvestmentVariables.R).

These source revisions motivate the request. Their identity with the exact
local runs remains to be demonstrated.
