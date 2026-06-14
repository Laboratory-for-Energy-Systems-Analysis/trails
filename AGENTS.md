# Repository Guidelines

## Project Structure & Module Organization
`trails/` contains the library source (core logic, LCA workflows, plotting, and IO).
`trails/data/` holds packaged datasets shipped with the library (notably
scenario inputs under `trails/data/scenarios/`).
`tests/` contains pytest suites and fixtures; test data lives under `tests/data/`.
`assets/` stores documentation and branding assets used by the README/docs.
`docs/` holds Sphinx documentation sources.
`examples/` and `dev/` include usage examples and local development utilities.

## Manuscript Location
The manuscript working files are stored outside the repository at
`/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/trails/manuscript`.
Check this folder when asked to inspect or edit the manuscript.

## Module Map (Key Files)
- `trails/lca.py` main LCA workflow entry points.
- `trails/lcia.py` LCIA data handling and characterization logic.
- `trails/datapackage.py` Frictionless datapackage loading and inventory parsing.
- `trails/characterization.py` characterization factor generation and utilities.
- `trails/temporal_distributions.py` temporal distribution models and helpers.
- `trails/importer.py` Excel inventory import (bw2io) and matrix updates.
- `trails/fair_io.py` FaIR IO helpers and scenario parsing.
- `trails/fair_rf.py` radiative forcing integration using FaIR.
- `trails/plotting.py` plotting utilities (primarily Plotly-based).
- `trails/cache.py` and `trails/cache_interpolation.py` interpolation cache logic.

## Build, Test, and Development Commands
- `pip install -e .` installs the package in editable mode from this repo.
- `pip install -e .[testing]` adds the optional testing dependencies.
- `pytest` runs the full test suite.
- `pytest tests/test_lca.py` runs a focused subset of tests.

## Coding Style & Naming Conventions
Use 4-space indentation and standard Python naming (snake_case for functions, CapWords for classes).
Line length follows the flake8 configuration in `pyproject.toml` (max 88 characters, with E203/W503 ignored).
Module names in `trails/` are lowercase and descriptive (e.g., `temporal_distributions.py`).
Docstrings generally follow Sphinx-style fields (``:param:``, ``:type:``, ``:returns:``), and type hints
use PEP 604 unions (e.g., ``int | None``).

## Testing Guidelines
Tests use `pytest` with the conventions in `pytest.ini`:
- Files: `tests/test_*.py`
- Classes: `Test*`
- Functions: `test_*`
Markers in `pytest.ini` include `forked`. `pyproject.toml` also declares `ecoinvent` and `serial`;
keep marker declarations in sync when adding or updating markers.
Use `pytest -m ecoinvent` only if the required data is available locally.

## Commit & Pull Request Guidelines
Commit messages in this repo are short, imperative summaries (e.g., "Refine LCA defaults and reset behavior").
When opening a PR, include:
- A brief description of the change and its motivation.
- Links to relevant issues or discussions, if applicable.
- Screenshots or plots for any changes affecting outputs or visualizations.

## Configuration & Environment Notes
Python versions supported are `>=3.10` and `<3.12` (see `pyproject.toml`).
Project dependencies are defined in `requirements.txt`, loaded dynamically at build time.
Runtime dependencies include Brightway (`bw2calc`, `bw_processing`), `sparse`, `xarray`,
`plotly`, and `datapackage`. `scikit-umfpack` is optional (`.[umfpack]`) and enables
UMFPACK; otherwise the solver falls back to SciPy.

## Data Packages & Caching
`trails` loads Frictionless datapackages with inventories under
`inventories/<model>/<pathway>/<year>/...` (see `dev/example/` for a sample
datapackage layout).
Annual interpolation caches are stored under the platformdirs user data path
(`appname="trails"`, `appauthor="pylca"`) in a `cache/` subdirectory.

## Debugging & Logging
Core modules use standard library logging (`logging.getLogger(__name__)`).
`TRAILS_DEBUG_FLOW_ID`, `TRAILS_DEBUG_YEAR`, `TRAILS_DEBUG_ACTIVITY`,
`TRAILS_DEBUG_MAX_PULSES`, and `TRAILS_DEBUG_MAX_MATCHES` influence detailed flow-level
debug output.
Tests write debug artifacts to `.pytest-debug` by default, or `PYTEST_DEBUG_DIR` if set.
