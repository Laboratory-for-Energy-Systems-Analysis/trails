# Repository Guidelines

## Project Structure & Module Organization
`trails/` contains the library source (core logic, LCA workflows, plotting, and IO).
`trails/data/` holds packaged datasets shipped with the library.
`tests/` contains pytest suites and fixtures; test data lives under `tests/data/`.
`assets/` stores documentation and branding assets used by the README/docs.
`docs/` holds Sphinx documentation sources.
`examples/` and `dev/` include usage examples and local development utilities.

## Build, Test, and Development Commands
- `pip install -e .` installs the package in editable mode from this repo.
- `pip install -e .[testing]` adds the optional testing dependencies.
- `pytest` runs the full test suite.
- `pytest tests/test_lca.py` runs a focused subset of tests.

## Coding Style & Naming Conventions
Use 4-space indentation and standard Python naming (snake_case for functions, CapWords for classes).
Line length follows the flake8 configuration in `pyproject.toml` (max 88 characters, with E203/W503 ignored).
Module names in `trails/` are lowercase and descriptive (e.g., `temporal_distributions.py`).

## Testing Guidelines
Tests use `pytest` with the conventions in `pytest.ini`:
- Files: `tests/test_*.py`
- Classes: `Test*`
- Functions: `test_*`
Markers include `ecoinvent`, `serial`, and `forked`; use `pytest -m ecoinvent` only if the required data is available locally.

## Commit & Pull Request Guidelines
Commit messages in this repo are short, imperative summaries (e.g., "Refine LCA defaults and reset behavior").
When opening a PR, include:
- A brief description of the change and its motivation.
- Links to relevant issues or discussions, if applicable.
- Screenshots or plots for any changes affecting outputs or visualizations.

## Configuration & Environment Notes
Python versions supported are `>=3.10` and `<3.13` (see `pyproject.toml`).
Project dependencies are defined in `requirements.txt`, loaded dynamically at build time.
