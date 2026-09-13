# Changelog

## Unreleased

## 1.1.1 - 2026-09-13

### Changed

- Reduce cold datapackage initialization overhead by avoiding unused CSV fields,
  reusing parsed temporal pulses, and sharing sparse interpolation work between
  annual slices. The BrightCon DACCS benchmark improves from 112.3 s to 45.9 s
  (59%) with identical matrices and cache metadata; see
  [benchmark details](dev/cold_initialization.md).
- Accelerate FaIR climate-response calculations by loading calibrated
  parameters in bulk, reusing the shared pre-perturbation model state, and
  avoiding repeated calculations for independent gases and forcing channels.
  The BrightCon FaIR cell improves from 113.2 s to 36.6 s (68%) while
  retaining all 841 configurations, signed emissions and removals, the full
  output horizon, and flow/root attribution; see
  [benchmark and numerical validation](dev/fair_runtime_optimization.md).
  These optimizations apply across inventory types; the runtime benefit
  depends on the emission mix and timing.
- Use up to two workers by default for the optimized FaIR 2.2.4 path;
  explicit `per_species_workers` values remain supported. Other FaIR
  versions and unsupported model configurations retain the full
  calculation path.

### Fixed

- Handle empty and disjoint sparse anchor supports during annual interpolation.

### Added

- Regression tests for cold-cache round trips, sparse interpolation and
  independent temporal pulse lists, plus FaIR baseline reuse, unsupported
  versions, early perturbations and seeded stochastic continuation.

## 1.1.0 - 2026-09-10

### Added

- Factorized and chunked inventory backends for large temporal LCI and LCIA
  calculations, with bounded temporal inventory memory.
- Ecoinvent 3.12 defaults, versioned LCIA method data, and improved
  characterization-cache scoping.
- Signed routing and biosphere-credit preservation through LCI, LCIA and FaIR
  calculations.
- Location-aware activity search and expanded publication examples, figures,
  and supplementary-information scripts.
- Warn once per Python session on import when the interpolation cache exceeds
  5 GiB, showing its size and instructions for using `trails.clear_cache()`.
  The check reads file metadata only and never removes data automatically.

### Changed

- EDGES characterization factors are evaluated by inventory year.
- FaIR CO2 pulse-equivalent and per-species calculations were optimized.
- Static LCA now preserves temporal inventory state, and interpolation caches
  are stable for ZIP datapackages.
- Documentation and numerical examples were aligned with the current API and
  workflows.

### Fixed

- Report filesystem errors from `clear_cache()` instead of silently ignoring
  incomplete deletion.

## 1.0.1 - 2026-06-14

### Changed

- `Trails.temporal_routing()` is now adaptive by default. Omitting `max_depth`
  uses `max_depth=None` with `adaptive_relative_score_cutoff=1e-4`.
- Passing an integer `max_depth` without an adaptive relative cutoff keeps
  fixed-depth routing available for workflows that need the previous
  depth-based behavior.
- Documented the four public routing modes: default adaptive routing, adaptive
  routing with a custom relative cutoff, adaptive routing with a hard depth cap,
  and fixed-depth routing.
