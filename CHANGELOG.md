# Changelog

## Unreleased

### Added

- Warn once per Python session on import when the interpolation cache exceeds
  5 GiB, showing its size and instructions for using `trails.clear_cache()`.
  The check reads file metadata only and never removes data automatically.

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
