# Numerical example assumptions

This synthetic datapackage supports `1. simple numerical example.ipynb`.
The values illustrate temporal routing and changing background inventories;
they are not an empirically calibrated vehicle inventory.
The four scenario anchors are 2005, 2020, 2050, and 2100.

## Revised gasoline-car model

Activity 13 supplies passenger-car transport. A demand of 200,000 vehicle-km
requires one gasoline vehicle and one end-of-life treatment. The revised
inventory uses the following assumptions at every scenario anchor:

- Vehicle assembly occurs at the transport start year. It requires one chassis,
  one engine, and one gearbox (activities 17–19).
- Chassis production follows a normal distribution centred two years before
  assembly, with scale 0.5 and support from year −3 to −1. Engine and gearbox
  production occur at years −4 and −3, respectively.
- The chassis, engine, and gearbox require 1,000, 250, and 150 kg of steel.
  Steel production (activity 20) consumes 2 kWh per kg. No direct biosphere
  emissions are assigned to these new activities; their impacts arise through
  the electricity supply chain.
- Direct fossil CO₂ from vehicle production is set to zero, replacing the
  previous 10,000/9,000/8,000/6,000 kg per vehicle across the four anchors.
  The new upstream inventory is not calibrated to preserve those totals.
- Fuel requirements and driving emissions are distributed uniformly over
  offsets 0–16. Their `temporal_amount_source` remains `matrix`, so amounts
  follow the scenario coefficients at each pulse year. This replaces fuel
  offsets −9 to +7 and driving-emission offsets −8 to +8.
- End-of-life treatment is split into 30% at offset +17 and 70% at +18.
  Its electricity requirement is reduced from 5,000 to 2,500 kWh per vehicle,
  consumed at treatment time. This is a revised illustrative assumption,
  not an energy-conserving redistribution of the previous inventory.
- Wind-turbine construction is distributed before wind electricity production,
  using a normal distribution centred at −2 with support −5 to 0. The scale
  is unspecified and uses the library's default. The exchange amount is unchanged.

For a transport start in 2050, these offsets place assembly in 2050,
component production in 2046–2049, operation in 2050–2066, and treatment
in 2067 and 2068. Upstream supply-chain offsets can shift impacts further.

## Numerical checks

Checked with the repository source after commit `a61f01d`, using Python 3.11:

- TRAILS loads four 21 × 21 technosphere matrices and four 21 × 2 biosphere
  matrices with annual interpolation and interpolation caching disabled.
- Independent NumPy static solves at every anchor have residuals below 1e-9
  and finite, nonnegative supply quantities.
- A demand of 200,000 vehicle-km requires one vehicle, one treatment,
  1,400 kg of steel, and 5,300 kWh of electricity in the static solve.
- All 16 temporal distributions per anchor yield nonnegative weights summing
  to one within numerical tolerance. The end-of-life pulses are exactly
  `(17, 0.3)` and `(18, 0.7)`.

Static inventories for 200,000 vehicle-km, without temporal routing:

| Scenario year | Fossil CO₂ (kg) | Net non-fossil CO₂ (kg) |
|---|---:|---:|
| 2005 | 41,921.80640 | 0 |
| 2020 | 31,395.80986 | 0 |
| 2050 | 19,949.60500 | 0 |
| 2100 | 936.09355 | 0 |

The static net non-fossil CO₂ balance includes both uptake and release;
it does not imply that their temporally resolved climate effects cancel.
These checks establish numerical consistency, not empirical validity.
The full notebook, annual temporal-routing results, and FaIR calculations
have not been rerun as part of this datapackage review.
