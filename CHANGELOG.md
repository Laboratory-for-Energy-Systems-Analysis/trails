# Changelog

## Unreleased v1.0 backport

### Fixed

- Sparse annual interpolation now handles exchanges that appear or disappear
  between snapshots, including an empty anchor matrix. Check lookup positions
  before indexing so missing coefficients remain zero. The weighted linear
  interpolation rule and temporal routing behavior are unchanged.

## 1.0.0 - 2026-06-14

### Changed

- `Trails.temporal_routing()` is now adaptive by default. Omitting `max_depth`
  uses `max_depth=None` with `adaptive_relative_score_cutoff=1e-4`.
- Passing an integer `max_depth` without an adaptive relative cutoff keeps
  fixed-depth routing available for workflows that need the previous
  depth-based behavior.
- Documented the four public routing modes: default adaptive routing, adaptive
  routing with a custom relative cutoff, adaptive routing with a hard depth cap,
  and fixed-depth routing.
