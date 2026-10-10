# Strengthening the REMIND stock-vintage approach

**Checkpoint:** this document records the earlier bounded-pilot roadmap.
The subsequent eleven-group implementation and validation are documented in
[iam-preprocessor](https://github.com/romainsacchi/iam-preprocessor/blob/feat/stock-evidence-validation/docs/remind-evidence-results.md).
Further implementation is paused pending additional exact-run REMIND/EDGE-T data.
The original scope and findings below are retained as historical context.

Plan dated 10 October 2026, based on the completed bounded pilots. This is the
next research and implementation programme, not a claim that the work below is
already implemented. The [pilot review](stock-vintage-pilot-review.md) remains
the record of the existing release evidence. The earlier work is preserved on
`feat/stock-vintage-pilots`, including its report generator and results.

## Recommendation

Prioritise exact-run provenance, native stock **and service** cohorts, and
explaining the stock-versus-additions discrepancies. Then validate the method
against held-out observations and extend its scenario, geographic and temporal
coverage. More lifetime curves alone will not resolve these issues.

The current pilots establish quantity conservation and correct event routing
under declared assumptions. A stronger REMIND application needs separate
evidence for model fidelity, empirical validity and prospective LCI behaviour.
Passing one of those checks does not establish the others.

Keep common amortisation for this programme. Native cohorts improve the dates
and weights; they do not automatically establish individual lifetime-service
denominators or physical material demand. Water infrastructure remains a
separate constant-stock scenario because the inspected REMIND exports do not
provide its stock dynamics.

## Findings that determine the priorities

1. **The exact runs are not identified fully.** File hashes and scenario labels
   are pinned, but `3.5.2` is a directory label and the MIF model field says
   `REMIND`. Model configuration, input revisions and the matching versions of
   remind2, EDGE-Transport and reporttransport remain unconfirmed.
2. **Native vintages are missing from our exports.** The archived EDGE code
   exposes construction-year vehicle and service outputs. Power reconstruction
   currently approximates the native capacity equation. Neither establishes
   that the exact scenario run's intermediate outputs are available.
3. **Aggregate service reconciliation does not identify vintage use.** If every
   cohort's utilisation is multiplied by the same annual factor, that factor
   cancels from normalised weights. Consequently, the current ES/generation
   target reconciles the total without adding information about relative use
   of old and new cohorts.
4. **The alternative constraints disagree materially.** Fresh diagnostics of
   the eight archived vehicle/power cases show stock-target and reported-sales
   or additions reconstructions diverging. Those differences can arise from
   regional proxies, reinitialisation, survival, early exits or reporting
   conventions. They are not, by themselves, proof of erroneous REMIND output.
5. **The horizon is too short to assess lifetime robustness.** Profiles currently
   cover 2022–2030. Young PV cohorts have little retirement within that interval;
   close parent-age sensitivities therefore provide limited evidence about
   2050/2100 turnover or component disposal.
6. **The implementation is deliberately narrow.** Power additions assume centred
   five-year periods and reviewed constant lifetime values. Geography and
   technology mappings are pilot-specific. Actual background technology was
   held constant across anchors to isolate timing and conservation.

The reproducible diagnostic is
`dev/stock_vintage/diagnose_remind_reconstruction.py`. Its local CSV and source
hashes are saved under `build/reports/remind-robustness/` and remain outside Git.
It reports absolute and relative flow/stock discrepancies, utilisation effects
and age ranges. A small difference in mean age can coexist with a large stock
discrepancy; both need checking.

## Priority work packages

| Priority / package | Improvement and data needed | Deliverable and acceptance condition |
|---|---|---|
| **P0 — A. Exact-run identity** | Recover run IDs, model/configuration/input hashes, regional definitions, native periods, and reporting/EDGE revisions for the two pilot scenarios. Preserve required variables before premise filtering. | Machine-readable run manifest and variable dictionary. Every field has units, scope and date/period meaning; unknown lineage is visible. |
| **P0 — B. Native benchmark** | Obtain native transport stock/service cohorts and native power capacity contributions, early-retirement and initial-stock information. Start with the existing four technologies and regions. | Reproduce one scenario/region/technology at native resolution, then the second scenario. Cohort sums match stock and service outputs within documented precision. |
| **P0 — C. Constraint reconciliation** | Resolve Sales versus New Cap period semantics, initial-stock replacement, survival, idle/service exits and category/region differences. | A residual ledger with absolute and relative discrepancies and explanations; no unexplained balancing term is called physical retirement. |
| **P0 — D. Service by cohort** | Native vehicle service plus mileage/load; power output by vintage if represented. Otherwise measured age-specific use and explicit dispatch assumptions. | Cohort service weights respond to cohort-specific information; totals, units and exposure reconcile. Sources of inferred service remain labelled. |
| **P1 — E. External validation and geography** | Multi-year EIA/DfT observations, stable cohort identifiers and compatible service data; native-region coverage and technology crosswalks. | Held-out predictions and coverage audits, with performance against transparent baselines. Country/subfleet proxies stay separate from native-region claims. |
| **P1 — F. Longer horizons and transitions** | Native time-step widths, changing lifetime parameters, declining fleets, idle/reactivated plants and retirement allocation. | Period-preserving annualisation through 2050 and 2100; contraction and period-transition tests; uncertainty reported beyond the observed horizon. |
| **P1 — G. Prospective LCI integration** | Scenario-transformed inventory anchors, historical technology assumptions and changing lifecycle coefficients. | Real varying-background integration, separation of timing and technology effects, and continued signed conservation under the revised contract. |
| **P2 — H. Components and wider deployment** | Battery/inverter replacement histories, PV repowering, physical dismantling delays and additional technologies/regions. | Evidence-based renewal cases, documented residual uncertainty and performance/compatibility checks before default promotion. |

P0 means necessary for stronger claims about REMIND consistency, not a reason
to retract the bounded, explicitly reconstructed timing demonstrations.

## A–B. Obtain a small native reference before expanding the adapter

Prepare the [data specification](remind-stock-vintage-data-specification.md) for
the owners of the exact local runs. First seek files in the existing authorised
local scenario holdings; record what is present and missing. A compact export
is sufficient if it preserves the required dimensions and definitions. A full
GDX is useful but not mandatory. External requests are a later action; this
plan does not send them.

Start with SSP2-NPi2025 and SSP2-PkBudg650, EUR four-wheeler BEVs, EUR 40-tonne
liquid-fuel trucks, USA PV and USA non-CHP/unabated CCGT as far as the native
technology definitions permit. Export the necessary complete initial vintage
state, not only additions since 2022. Preserve native periods; require annual
data only where the model or reporter actually produces them.

Maintain three explicit result origins:

- **Native model result:** retains the IAM's own initial cohorts, regions,
  technology aggregation, survival and stock/service outputs.
- **Reconstruction initialised from observations:** replaces initial stocks
  with empirical cohorts and reconciles subsequent additions/exits. This is a
  different experiment, even when it follows the same macro trajectory.
- **Regional/technology proxy:** applies a larger region's relative path to a
  selected population, as the existing UK and WECC pilots do.

Compare these modes rather than silently blending them. Native reproduction
can establish consistency with REMIND but still requires external validation
of its historical assumptions. Conversely, observed reinitialisation changes
retirement needs and need not reproduce native gross additions.

## C. Reconcile stocks, additions and exits

Use an explicit accounting identity for each consistent population:

```text
closing stock = opening stock + gross new entries + used imports
                - natural exits - additional exits - used exports
                + documented reclassifications
```

For power, distinguish installed, operating and idle capacity and represent
reactivation where the model does. Additional operating exits are not a
physical demolition measurement. For vehicles, registration disappearance can
include exports or classification changes. A stock series alone cannot identify
these terms.

Resolve discrepancies in this order:

1. Check units, leaf/aggregate membership, AC/DC basis, model regions and the
   timing conventions of stocks versus annual service.
2. Confirm reporting transformations against the exact run. Vehicle Sales is
   sampled annual entrants in the inspected reporter; power New Cap is a rate
   with a different period convention. Do not use the same annualiser for both.
3. Reproduce native initialisation, native survival indexing and early exits.
4. Quantify what changes when observed cohorts replace the native initial state.
5. Only then choose a constrained reconstruction for genuinely missing terms.

A joint constrained fit may be useful, with nonnegative additions/exits and
documented uncertainty weights on conflicting observations. It is not uniquely
identified by stock and additions totals. Keep alternative solutions, report
which constraints were relaxed and avoid allowing arbitrary balancing flows to
hide a wrong crosswalk. The stock-target and additions-target modes remain
diagnostic endpoints.

## D–E. Strengthen service weights and independent validation

| Asset | Highest-value data gap | Data route and limitations | Validation experiment |
|---|---|---|---|
| Passenger BEV | Vehicle-km by vehicle age, plus native service/stock cohorts | Exact-run EDGE outputs first. DfT/DVSA anonymised MOT odometers are a public acquisition lead; NTS mileage/occupancy tables provide aggregate checks. Verify identifiers, age/fuel fields, coverage and revisions before claiming a joined BEV panel. | Predict later first-use distributions from an earlier stock; validate age-specific mileage separately. Include young cars absent from inspection samples and England/GB/UK coverage differences. |
| Heavy trucks | Mileage, payload and empty running by age/type; stock-sales discrepancy | Exact-run EDGE mileage/load/service data, DfT VEH0520 and road-freight statistics. Published type/weight tables do not establish an age-by-payload cross-tab; microdata access is unresolved. | Test vehicle counts and tonne-km independently. Keep diesel, liquid-fuel, weight-class and EURO-class mappings explicit. |
| PV | Vintage output, commissioning exposure, module/system classification and AC/DC mapping | Annual EIA-860 plus EIA-923; monthly records where sufficiently complete. Expand beyond the current single-cohort fixed-crystalline plant subset only with auditable allocation. | Predict capacity and cohort output in held-out years; assess degradation, curtailment, new entrants and selection bias. Distinguish utility-scale and distributed PV. |
| CCGT | Cohort-specific operating exits, idle capacity and dispatch | Native REMIND capacity/retirement state; longitudinal EIA plant/block status and output. Resolve CT/steam-block joins and refurbishment. | Stress the declining fossil scenario. Compare proportional, age- and efficiency-informed exit/dispatch assumptions against historical evidence where possible. |

The [EIA-860 documentation](https://www.eia.gov/electricity/data/eia860/) provides
operating/standby and retired/cancelled categories but warns that a single
annual retired tab is not a complete retirement history. Build a longitudinal
panel and assess EIA-860M's historical retirement records. Pair with
[EIA-923](https://www.eia.gov/electricity/data/eia923/) output at its actual
reporting level; complete generator-level service must not be presumed.

The [DVSA directory](https://open.data.dvsa.gov.uk/) and
[MOT archive](https://www.data.gov.uk/dataset/e3939ef8-30c7-4ca8-9c7c-ad9475cc9b2f/anonymised_mot_tests_and_results)
advertise anonymised tests and odometers; the archive records revisions to
historical files. The [NTS tables](https://www.gov.uk/government/statistical-data-sets/nts09-vehicle-mileage-and-occupancy)
and [road-freight tables](https://www.gov.uk/government/statistical-data-sets/domestic-road-freight-activity-rfs01)
are aggregate validation leads, not verified full cohort-use datasets. These
pages were checked for this plan; their underlying files have not been newly
curated here.

Backtesting should have two tracks. A **conditional cohort test** can use known
later totals/additions to isolate errors in survival and vintage allocation. A
**full historical replay** withholds later totals as well. Report them separately;
fitting a year and then validating against that same constraint is circular.
Use earlier years for fitting and genuinely held-out years for assessment.

Report total stock/service errors, additions and exits, age-bin/CDF errors,
Wasserstein distance in years, and mean/median/tail ages. Weight results by
population/service and retain exclusions. Include a simple equal-use reference
and the present reconstruction as baselines. Predeclare acceptance criteria
using source precision and sampling uncertainty; there is no defensible
universal empirical percentage threshold yet.

## F. Make the temporal representation general

Replace fixed five-year assumptions with an explicit mapping from native
investment intervals to reporting intervals and annual events. Preserve
integrated additions across interval boundaries, including partial initial and
terminal intervals and the switch to longer native steps. Annual detail within
a multi-year investment interval is an assumption unless separately observed.

Handle time-varying lifetime inputs according to whether they apply at build
year, evaluation year or another native index. Preserve initial cohorts older
than an assumed support through an explained treatment rather than deleting
them. Make zero stock/service, technology entry, complete exit, missing data,
CCS conversions, retrofits and reactivation explicit cases. Missing is not zero.

For annual service, move from an end-year portfolio to time-integrated exposure
when data permit:

```text
cohort annual service = integral over the year of
                       surviving stock × service rate
```

Use commissioning/exit months for observations, or a declared within-year
entry rule. Respect leap-year/exposure hours and technology-specific feasible
utilisation. An aggregate electricity path alone does not identify future
cohort dispatch; retain uncertainty when native vintage output is absent.

Test 2030, 2050 and 2100, both existing scenarios, and then one clearly separate
EU21 geography and one additional region. Verify the exact region membership
before adding either. Sensitivity dimensions should cover initialisation,
reporting periods, survival, early-exit allocation, service weighting and
geographic transfer; scenario differences are not confidence intervals.

## G–H. Couple robust cohorts to changing inventory technology

The next real-database comparison should use actual premise scenario changes
across inventory anchors, retaining the constant-background case as a control.
Separate changes caused by dates, cohort service, operating technology and
historical production technology. An old asset's manufacture uses historical
production conditions only where those inventories exist; the present nearest
inventory-year lookup is not evidence of historical technology reconstruction.

Review `port` versus `matrix` amount sourcing and the existing lifecycle
rewrite's requirement that lifted coefficients are invariant across anchors.
Changing background technology or component ratios can invalidate that
condition. Design and test a compatible extension before relaxing the guard;
keep common amortisation and signed quantity accounting explicit.

Avoid applying a technology mix or fleet-efficiency correction twice when
premise already represents it. Record whether an inventory describes a new
asset, a fleet-average service or a maintenance event. Component chronology
needs separate observations even if the parent IAM cohorts are exact.

Retain the current three release requirements, adding scope-specific checks:

- **Boundary:** reviewed units, caller/supplier roles, region/technology
  crosswalks, service year and manufacturing technology year.
- **Cohort evidence:** native-model reconciliation or an explicitly labelled
  constrained reconstruction, plus held-out empirical validation and disclosed
  residuals/coverage.
- **Integration:** annual lookup, changing anchors, both demand signs, cache,
  past manufacture/future disposal and full-LCI/operator consistency. Larger
  horizons also need runtime and memory evidence and explicit reader support.

## Order of execution and decisions

1. **First tranche: diagnostics and data packaging.** Complete the residual
   dashboard, annotate run/variable lineage, and prepare a minimal exact-run
   export request. Run the historical validation design against existing data.
   Estimated effort: 2–4 analyst-days; no external response is assumed.
2. **Second tranche: native consistency and service.** Once exact data arrive,
   build the smallest native benchmark and resolve period/survival/exit
   differences. Develop empirical-use evidence in parallel with that work.
   Estimated effort: 5–10 analyst-days after usable files are available.
3. **Third tranche: generalisation and prospective integration.** Extend
   periods/horizons, backtest across populations, then run varying-background
   LCI checks and sensitivity attribution. Estimated effort: 5–10 analyst-days
   after the preceding methods are settled; data curation may add time.

These are planning estimates, not delivery commitments. Native data availability
and quality determine whether tranche two can demonstrate exact fidelity. If
the files cannot be obtained, finish a transparent constrained-reconstruction
route with bounded claims; do not relabel it as native reproduction.

**First decision point:** can the exact run export stock/service cohorts and
their reporting lineage? **Second:** can discrepancies be explained without
unobserved arbitrary balancing? **Third:** do held-out distributions and
varying-background integrations support the intended scope? Only then consider
additional technologies or default promotion.
