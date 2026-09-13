# Changelog

## Unreleased

### Changed

- Reduce cold datapackage initialization overhead by avoiding unused CSV fields,
  reusing parsed temporal pulses, and sharing sparse interpolation work between
  annual slices. The BrightCon DACCS benchmark improves from 112.3 s to 45.9 s
  (59%) with identical matrices and cache metadata; see
  [benchmark details](dev/cold_initialization.md).

### Fixed

- Handle empty and disjoint sparse anchor supports during annual interpolation.

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
