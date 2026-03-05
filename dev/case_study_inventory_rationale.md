# Temporal Case-Study Inventory Rationale

This note summarizes the modeling logic used to revise the Excel case-study
inventories in `dev/lci-case-study-*.xlsx`.

## Design Principles

- Keep exchanges linkable to existing datapackage datasets/flows.
- Add ancillary exchanges that are plausible for operations and maintenance.
- Use temporal distributions that match process behavior:
  - `6` (discrete empirical) for event-based pulses and delayed tails.
  - `4` (uniform) for continuous operation over service life.
  - `5` (triangular) for concentrated one-off phases (e.g., construction).
- Put rationale in each exchange `comment` cell so assumptions are inspectable
  directly in the workbook.

## File-by-File Intent

## `lci-case-study-daccs_storage_risk.xlsx`

- Added ancillary O&M inputs (diesel logistics, low-voltage auxiliaries, backup
  natural gas) and associated methane/N2O biosphere emissions.
- Modelled operation-related exchanges (including atmospheric CO2 uptake) with
  a uniform profile over the 20-year plant life to represent continuous
  operation.
- Kept a long-tail post-closure leakage risk profile for fossil CO2 to preserve
  the original storage-risk focus.

## `lci-case-study-fertilizer_n2o_timing.xlsx`

- Added ancillary field-operation exchanges (diesel, water, direct fossil CO2,
  methane from fuel chain).
- Reworked fertilizer timing to be front-loaded, then reduced over later
  applications.
- Updated N2O and NH3 timing to include pronounced short-term peaks plus
  decaying multi-year tails.

## `lci-case-study-marine_fuel_switch.xlsx`

- Strengthened matrix-driven transition dynamics with clearer year-specific
  trajectories for diesel, methanol, electricity, and emissions.
- Added transitional natural-gas demand and small N2O component.
- Updated methane slip to peak during transition years and decline with
  technology maturation.

## `lci-case-study-biomass_growth_vs_gas_heat.xlsx`

- Added ancillary exchanges across forestry operations and heat-service stages
  (water, diesel maintenance, methane tails, additional direct emissions).
- Refined growth/harvest temporal profiles to better distinguish short vs long
  rotation behavior.
- Added/updated methane and N2O tails where decomposition/combustion timing is
  climate-relevant.

## New Example Inventories

## `lci-case-study-refrigerant_leakage_cooling.xlsx`

- Cooling service with annual HFC leakage and a strong end-of-life release
  pulse to highlight timing effects of high-GWP refrigerants.

## `lci-case-study-ccu_polyol_delayed_release.xlsx`

- CCU chain with captured CO2 input and delayed end-of-life non-fossil CO2
  release to demonstrate temporary carbon storage dynamics.
