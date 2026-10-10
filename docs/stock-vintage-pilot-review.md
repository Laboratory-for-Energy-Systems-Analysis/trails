# Review of the five bounded stock-vintage pilots

Review in progress, 10 October 2026. Both implementation repositories use
`feat/stock-vintage-pilots`. CCGT, PV and trucks have evidence for all three requirements.
Passenger cars and water have passed their complete matrix rewrites and
annual event checks; their corrected full-LCI comparisons are still running.
The autonomous goal is not complete. No public defaults have changed.

This review concerns opt-in **timing under the original common amortisation**.
For each service-year exchange, the original coefficient is distributed over
the serving stock's dates. Service shares require explicit utilisation weights;
they equal stock shares only under the stated equal-service assumptions. No
second lifetime denominator or replacement multiplier is applied. This is not
a physical material-demand forecast or heterogeneous individual-life allocation.

## Requirement 1: exact boundaries and conserved coefficients

| Pilot | Service and capital boundary | Components and lifecycle treatment |
|---|---|---|
| UK passenger BEV | GLO electric-car service, kg chassis per vehicle-km; source 918.22 kg and 150,000 km retained | Battery manufacture and end are separate; 17 parent disposal ports are lifted; four glider production scraps stay at manufacture. Six event markets have reviewed instantaneous internal links. Source battery coefficient already contains 1.5 equivalent packs. |
| UK articulated diesel truck | RER EURO6 service, 40-tonne lorry units per tonne-km; `9.65e-8` retained | Used-lorry disposal is separate; three manufacturing wastewaters stay at manufacture. Maintenance is assigned to service. Exact maintenance/disposal providers are scoped to prevent repeated delays. |
| US-WECC PV | Source multi-Si ground-mounted 570 kWp inventory per kWh; `4.42343137603572e-8` installation units retained | Parent, active inverter, panels, handling/maintenance losses and lifecycle ports are separated. Occupation is assigned to service; land transformation stays at construction. Ambiguous polymer/packaging dates have tested endpoints. |
| US-WECC gas CCGT | Non-CHP combined-cycle generation, 400 MW plant units per kWh; `1.38888888888889e-11` retained | Construction is dated by serving stock; independent factory capital keeps its own role. No absent end-of-life inventory is invented. |
| Quebec potable water | Pipe/network capital per kg water; own-product loss 0.191 and net production 0.809 retained | Nine tank-material inputs, two tank disposal ports and 22 pipe disposal ports are separated. Pipe and tank dates remain distinct. |

Scoped copies leave unrelated callers and the original datasets intact. The
complete A/B operator rewrite proof passes for all five pilots. Annual checks
cover both functional-unit signs, repeated/reversed requests and exact internal
zero shifts. Signed product requirements are divided by signed reference
production in the opt-in consumer. Legacy negative-reference errors are retained
as separate diagnostics; they are not consequences of the stock method.

## Requirement 2: observations, cohort balance and uncertainties

| Pilot | Observed reference stock in 2022 | Primary weighting and future reconstruction | Material limitations and sensitivity treatment |
|---|---|---|---|
| Passenger BEV | 628,318 UK cars with known first-use years 2010–2022; 99.8941% of the BEV count | Equal service per vehicle; relative REMIND EUR stock trajectory, annual entrants closing stock balance | First use is not necessarily manufacture or battery installation. Pre-2010 and unknown dates are explicit exclusions. Unknown dates, entrant exposure, age-dependent service, survival support, sales constraints and battery interval are varied. |
| Heavy truck | 14,231 UK articulated diesel vehicles, 32–40 tonnes | Equal service per vehicle; relative REMIND EUR liquid-fuel 40-tonne stock trajectory | 370 unknown dates are distributed proportionally; one pre-1980 vehicle has an explicit tail. Source EURO6 inventory is a technology proxy, not an observed fleet classification. Missing mileage/loading, survival and stock/sales conflicts are exposed. |
| PV | 363 complete EIA plants: 1,294.8 MW AC, 1,629.2 MW DC; observed net output | Output-weighted serving stock, including zero-output capacity; relative REMIND USA PV trajectory | Fixed crystalline plant classification proxies multi-Si ground-mounted technology. Parent life, inverter replacement, capacity weighting and additions constraints are varied; source waste ambiguities have real-export endpoints. |
| CCGT | 91 complete EIA blocks: 41,476.3 MW AC | Observed output weights; relative REMIND USA CCGT trajectory | The 410.5 MW beyond the survival prior receives an explicit conditional residual-life extension. Service/capacity weights, survival, old-stock extension and additions constraints are varied. |
| Water | 45,114 km Quebec pipes and 974 storage assets | Pipe length and storage counts; constant stock with replacement through 2030 | Unknown construction dates are 5.1603% and 25.5647%. Tank count is a proxy for volume/type. Source bins, missing dates, survival, service weighting and retirement timing have explicit bounds. No IAM water forecast is claimed. |

The annual outputs cover 2022–2030. The two pinned REMIND scenarios generate
126 CCGT, 144 PV and 468 vehicle case/year records; water supplies 234
asset/case/year records. Balance residuals are numerical round-off. Initial
observed survivors are conditioned on already being alive, not multiplied by
unconditional survival a second time. Future stock and sales constraints retain
their discrepancies. Territorial vehicle exits are separated from the assumed
physical retirement used to date disposal.

Exact native IAM vintage arrays and exact-run transport parameters were not
available. Pinned model/reporting code establishes variable and survival
semantics, but these remain labelled reconstructions. The power and vehicle
kernels differ. Their support parameters are not silently called mean lives.
Observed component replacement histories and age-specific service data remain
missing. Declared timing assumptions and sensitivities make bounded pilots
reviewable; they do not turn those assumptions into measurements.

## Requirement 3: actual producer and consumer evidence

Every package contains the complete local ecoinvent 3.12 cut-off source, exported
by premise at 2022/2025/2030. Technology is held constant to isolate timing and
quantity conservation. These runs do not demonstrate historical technology
reconstruction or an IAM-transformed background LCA.

| Pilot | Matrix rewrite and annual signed ports | Full 2023/2030 signed LCI | Cache/lifecycle evidence |
|---|---|---|---|
| CCGT | Pass | Pass, all 9,847 biosphere rows | Annual interpolation and warm cache; exact parent dates |
| PV | Pass; 510 port checks, 58 zero shifts | Pass, plus both 2030 waste-timing endpoints | Confirmed warm cache; direct occupation at service; component manufacture/end chronology |
| Passenger BEV | Pass; 210 port checks, 30 zero shifts | Running | Confirmed warm cache; parent, battery, disposal and maintenance calendars checked |
| Heavy truck | Pass; 30 port checks, six zero shifts | Pass, all 9,847 biosphere rows | Confirmed warm cache; parent/disposal and service-year maintenance checked |
| Water | Pass; 340 port checks, 34 zero shifts | Running | Separate confirmed warm-cache pulse check; pipe/tank and disposal calendars checked |

Full-LCI checks compare the summed biosphere inventory against an independently
assembled year-wise calculation using the same numerical backend, and separately
prove the graph operator identity. This prevents a wrong route from passing
merely because the numerical reference follows that route. Different sparse
factorisations of the poorly conditioned full database are kept as numerical
diagnostics, not mixed with routing errors. Both +1 and −2 demands are checked.

The private evidence records package hashes, scalar residuals, every signed
port calendar and role-attributed fossil-CO2 series. Synthetic integration tests
also cover varying anchors, resource hashes, absent annual years, opt-in guards,
foreground overrides, multiple solvers and legacy behaviour. See
[the implementation contract](stock-vintage-development.md).

The dated emission illustrations select the biosphere flow **Carbon dioxide,
fossil / air / unspecified**. They are not a sum across every CO2 compartment,
a GWP result or a dynamic climate score. Conservation checks cover all exported
biosphere rows independently of that illustrative selection.

## Actual parent-calendar comparison

`dev/stock_vintage/report_parent_calendars.py` loads all ten real packages and
compares unclamped expansion with actual routing for one unit of 2030 service.
It checks conserved totals and exact corrected dates. The restricted local
CSV and PNG/SVG/PDF figure distinguish these three means:

| Parent | Legacy distribution age | Legacy routed age | Corrected age |
|---|---:|---:|---:|
| BEV chassis | 8.0000 | 6.1333 | 3.6107 |
| Heavy truck | 0.0000 | 0.0000 | 7.4611 |
| PV installation | 13.5018 | 7.3190 | 5.6975 |
| CCGT plant | 16.9338 | 7.8052 | 23.4754 |
| Water pipe/network | 35.0000 | 7.9363 | 36.7757 |

All values are years. The selected legacy truck manufacturer is classified as a
throughput process, so it has no stock-age spread. The corrected result therefore
addresses missing stock timing as well as unsuitable distributions. The other
legacy results also exhibit calendar clamping. These parent-only means exclude
active battery/inverter timing, tank timing, disposal and background factory
capital. The full lifecycle and sign checks remain necessary.

## Remaining completion work

1. Finish and review the running corrected water and BEV full-LCI reports,
   including both demand signs and legacy quantity/timing diagnostics.
2. Update this decision table from the completed evidence; do not infer a pass
   from a running process or a successful parent-calendar check.
3. Verify and create the durable private evidence archive with both clean code
   revisions, restricted inputs, outputs and environment versions.
4. Complete the final cross-repository review and record the exact evidence
   locations and remaining applicability limits. Promotion here means a bounded
   opt-in pilot, not a public stock dataset or a change to defaults.
