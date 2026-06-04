# Changelog

## Unreleased

### Changed

- `Trails.temporal_routing()` is now adaptive by default. Omitting `max_depth`
  uses `max_depth=None` with `adaptive_relative_score_cutoff=1e-4`.
- Passing an integer `max_depth` without an adaptive cutoff keeps fixed-depth
  routing available for workflows that need the previous depth-based behavior.
- Documented the four public routing modes: default adaptive routing, adaptive
  routing with a custom relative cutoff, adaptive routing with a hard depth cap,
  and fixed-depth routing.

