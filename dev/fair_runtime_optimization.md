# FaIR runtime optimization — BrightCon 2026

The unchanged notebook FaIR cell completed in **36.59 seconds**, compared with **113.24 seconds** in the original full-notebook timing: about **68% faster**. This is a fresh-process run with no persisted FaIR response cache. The installed code reproduced this in **36.56 seconds** in a fresh Jupyter kernel. The entire notebook completed in **159.80 seconds (2 min 39.8 s)**, down from 239.29 seconds in the original rehearsal; all 37 nonempty code cells passed. Both interpolation caches were warm in the final rehearsal, whereas the original car cache was built during its run.

## Calculation preserved

The notebook still calls `run_fair_delta_rf(daccs, scenario="REMIND|SSP2-PkBudg1000")` with its defaults. It retains all **841 calibrated configurations**, **52 signed species perturbations plus the baseline**, **751 output timebounds (1750.5–2500.5)**, five quantiles, and attribution by biosphere flow and root activity. The DACCS amount remains 20 billion kg, with the intervention in 2035. No ensemble subset, shorter output horizon, relaxed precision, or skipped emissions is used.

Complete sparse-output comparison against the original implementation gives:

- Radiative forcing: identical sparse coordinates and values; maximum absolute difference **0**.
- Temperature: maximum absolute difference **4.6871e-17 K** across every quantile, year, flow and root. Floating-point regrouping creates some additional values near zero; they were retained in the comparison.
- Native and vectorized calibration loading produce identical climate and species configuration datasets for every ensemble member.

## Implementation

1. Load calibrated parameters by column across the ensemble, replacing repeated scalar assignments.
2. Save a shared pre-perturbation FaIR state, including gas partitions, cumulative and airborne emissions, temperature layers, and stochastic forcing. Most perturbations begin in 2004; earlier land-use CO₂ perturbations retain the complete historical calculation. All historical output labels are preserved.
3. Reuse gas concentrations and forcing channels only when their drivers and feedback coefficients prove them independent of the perturbation. CO₂ feedbacks and major-gas forcing overlaps remain active.
4. Combine eligible unchanged minor gases into a prescribed background-forcing channel inside FaIR. Their climate forcing remains included. Eligibility requires no relevant temperature or downstream chemistry response and unit forcing efficacy.
5. Use two workers by default for this optimized path. On this 24 GiB machine, additional workers increased contention and memory use. Explicit `per_species_workers` remains supported.

State reuse and compaction are limited to the validated **FaIR 2.2.4** engine, non-prescribed temperature, the Leach methane method, and configurations without inverse GHG or EESC coupling. Unsupported modes retain the full model path. This is a compatibility guard, not a dependency downgrade.

## Validation and reproducibility

- **125 focused tests passed**: signed uptake/emissions and root attribution; native FaIR comparison for eight perturbation species, one/three configurations, and all four supported GHG forcing methods; early perturbation fallback; seeded stochastic continuation; full-ensemble parameter equality; existing pulse-equivalence and regression tests.
- Original notebook SHA-256: `0f17f885eee9a9342f5e2a5a862be80c7d1ace5c9a9811539625cdf638dfd55a`.
- Optimized `fair_rf.py` SHA-256: `0a5340f0a0049e4bd37e63f4ea35149fe4992c1fd9995498c09cf92b5f5c5f4c`.
- Environment: Python 3.11.14, FaIR 2.2.4, TRAILS source 1.1.0, EDGES source 1.4.1, macOS, 24 GiB RAM. The repository base revision was `603f57a315d4b6fa53070ac48020f65e83f6f1dc`; pre-existing local workflow changes were retained in both baseline and optimized runs.
- The independent profiled reference cell took 113.36 seconds. The 113.24-second reference above comes from the original unprofiled full-notebook execution.

The conference workspace contains `dev/fair_profile/run_fair_from_notebook.py`, which executes every preceding notebook code cell and then the exact FaIR cell. Its `--workers` option is for experiments only; the accepted 36.59-second run did not use it. Run it from the conference folder with the `trails` Python environment and a new `--out` directory. The snapshot used for development replay excludes inventory reduction and was not used to claim the notebook result.

Raw timing summaries and sparse-output hashes are saved in `results/fair_optimization/{baseline,optimized_cell,parity}.json`. The baseline full-notebook run is in `results/notebook_timing/run_20260913_full`. The final installed-code rehearsal is in `results/notebook_timing/run_20260913_fair_optimized`.

No notebook code or input data was changed for this optimization. The original full-notebook timing included existing interpolation caches; the final rehearsal also uses available disk caches. Human pauses, presentation time and browser rendering are outside the timing.
