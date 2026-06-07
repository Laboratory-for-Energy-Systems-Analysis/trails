# TRAILS Knowledge File

Date of repository analysis: 2026-05-13

This file summarizes how the `trails` project operates, what it implements, what
it assumes, and where its current boundaries are. It is based on the repository
source code, tests, documentation, README, and example data package available in
this checkout.

## 1. Project Identity

`TRAILS` stands for `Temporal Routing and Aggregation of Impacts across
Life-cycle Systems`. It is a Python package for temporal life cycle assessment
that routes demands and emissions through time-resolved supply chains.

The package is designed around matrix-based LCA data and is especially aligned
with `premise`-style Frictionless datapackages containing scenario-specific
technosphere and biosphere matrices. Its central object is `trails.Trails`,
which loads sparse matrices, interpolates them to an annual scenario grid,
performs temporal routing from a functional unit, solves year-specific LCA
systems for the remaining frontier demands, and accumulates either characterized
scores or inventories.

The public API exposed by `trails/__init__.py` includes:

- `Trails`
- `lca`
- `get_lcia_method_names`
- `plot_temporal_scores`
- `plot_adaptive_sankey`
- `plot_rf`
- `plot_temp`
- `clear_cache`
- `search_activity`

The repository metadata currently describes the project as
`Country-specific characterization factors for Brightway`; this appears stale
relative to the README, docs, tests, and package code, which implement temporal
LCA workflows.

## 2. What TRAILS Does

TRAILS implements a dynamic or temporal LCA workflow with these main
capabilities:

- Load scenario-year technosphere and biosphere matrices from a Frictionless
  datapackage.
- Represent matrices as sparse three-dimensional arrays:
  - `A`: `(scenario, activity, product)`
  - `B`: `(scenario, activity, biosphere_flow)`
- Interpolate scenario matrices to an annual grid, with optional persistent
  cache.
- Attach temporal distributions to technosphere and biosphere exchanges.
- Route a functional unit through a time-dependent supply chain graph.
- Solve year-specific static LCA systems for unresolved frontier demands.
- Accumulate time-resolved biosphere inventories and/or characterized scores.
- Attribute results to first-tier root activities when requested.
- Use bundled LCIA characterization data for ecoinvent 3.10 and 3.11.
- Import additional inventories from bw2io-compatible Excel workbooks.
- Interface with Brightway/bw2calc through `bw_processing` datapackages.
- Propagate inventory perturbations through FaIR to estimate delta radiative
  forcing and delta temperature.
- Generate Plotly visualizations of temporal scores, traversal graphs,
  radiative forcing, and temperature responses.

In methodological terms, TRAILS combines temporal traversal of explicitly
temporalized exchanges with conventional matrix inversion or iterative solves
for parts of the system that remain static within a given year.

## 3. What TRAILS Does Not Do

The current implementation has clear boundaries:

- It does not solve one global time-expanded technosphere matrix. It performs
  routing and then solves conventional year-specific systems.
- It does not provide sub-annual temporal resolution. The temporal grid is
  calendar-year based, and temporal offsets are integer years.
- It does not perform stochastic Monte Carlo uncertainty sampling from the
  uncertainty fields in the matrix CSV files. Those fields are parsed as part of
  the data schema, but the implemented workflow is deterministic.
- It does not generate `premise` scenarios. It consumes already prepared
  datapackages.
- It does not automatically regionalize characterization factors beyond the
  matching embedded in the bundled LCIA JSON files.
- It does not infer semantic matches between biosphere flows and LCIA flows
  beyond exact metadata matching on flow name and compartments.
- It does not guarantee exact preservation of temporal metadata during annual
  interpolation when incompatible temporal exchange definitions bracket the
  requested year. In such cases, it falls back to nearest temporal metadata.
- It does not offer a graphical user interface; interaction is through Python
  APIs and generated Plotly figures.
- It does not make unavailable scenario years physically meaningful. Years
  outside the matrix year range are clipped to the nearest available scenario
  year for matrix solving, although inventories and emissions can be accumulated
  on a wider year axis.

## 4. Repository Map

Core package:

- `trails/trails.py`: main `Trails` class, matrix access, temporal routing,
  traversal helpers, inventory and score accumulation.
- `trails/lca.py`: LCA orchestration, solver selection, static LCA helper,
  node scoring and temporal sankey tree utilities.
- `trails/datapackage.py`: Frictionless package loading, matrix construction,
  annual interpolation, temporal exchange parsing.
- `trails/temporal_distributions.py`: temporal exchange dataclass and temporal
  distribution weight generation.
- `trails/bw_interface.py`: conversion between Trails sparse matrices and
  `bw_processing` datapackages for `bw2calc`.
- `trails/lcia.py` and `trails/characterization.py`: LCIA method loading,
  characterization vector construction, and characterized inventory building.
- `trails/importer.py`: bw2io Excel inventory import into existing Trails
  matrices.
- `trails/fair_io.py` and `trails/fair_rf.py`: FaIR input loading, emissions
  perturbation construction, radiative forcing and temperature calculations.
- `trails/plotting.py`: Plotly-based visualizations.
- `trails/search.py`: activity and biosphere flow search.
- `trails/cache_interpolation.py` and `trails/cache.py`: annual interpolation
  cache keying, serialization, and clearing.
- `trails/iterative_solver.py`: multi-RHS GMRES solver with preconditioners.
- `trails/logging.py`: logging configuration and contextual log fields.

Packaged data:

- `trails/data/lcia_ei310.json`
- `trails/data/lcia_ei311.json`
- `trails/data/scenarios/fair_species_map.yaml`
- `trails/data/scenarios/ssp_emissions_2024-11-26.csv`
- `trails/data/scenarios/calibrated_constrained_parameters_calibration1.4.1.csv`
- `trails/data/scenarios/species_configs_properties_calibration1.4.1.csv`

Project support:

- `docs/`: Sphinx documentation, including user guide, data format,
  interpolation, Excel import, FaIR, and plotting pages.
- `tests/`: pytest suite with unit, regression, plotting, importer, FaIR, and
  LCA behavior checks.
- `examples/example data package/`: minimal Frictionless datapackage used as a
  concrete example of the expected resource layout.
- `dev/`: development utilities, publication material, profiles, local reports,
  and data experiments.

## 5. Core Conceptual Model

### 5.1 Scenario Years

Input datapackages provide one matrix set per scenario year, such as `2005`,
`2020`, `2050`, and `2100`. TRAILS stores these labels in `scenario_labels` and
maps labels to integer positions with `scenario_index`.

When annual interpolation is enabled, the matrix grid is expanded to annual
labels. By default, initialization uses:

```python
Trails(
    package,
    interpolate_annual=True,
    cache_interpolation=True,
    interpolation_start_year_offset=-1,
    interpolation_end_year_offset=1,
)
```

For an input range of 2005 to 2100, the default annual grid therefore extends
from 2004 to 2101. Tests confirm this behavior. Custom offsets can extend the
annual grid further.

TRAILS also retains `template_labels`, corresponding to the original input
scenario years. These are used for temporal metadata and for importer behavior.

### 5.2 Matrix Orientation

Inside the `Trails` object:

- `A[t, activity, product]` stores technosphere exchanges for scenario index
  `t`.
- `B[t, activity, biosphere_flow]` stores biosphere exchanges for scenario
  index `t`.

The row axis of `A` and `B` is the producing or consuming activity being
evaluated. The column axis of `A` is the product/activity required by that
activity. During Brightway conversion, the technosphere matrix is transformed to
Brightway orientation `(product, activity)`.

Production self-exchanges in `A` are identified by `product_index == act_idx`
and skipped during temporal traversal expansion. Consumption exchanges are
converted into positive child requirements by taking the negative of negative
technosphere values.

### 5.3 Temporal Exchanges

Technosphere and biosphere exchanges can have temporal distributions. A temporal
exchange describes how an exchange anchored in one year is distributed over
other years.

The parsed `TemporalExchange` fields are:

- `distribution`
- `loc`
- `scale`
- `offset_min`
- `offset_max`
- `amount_source`
- `offsets`
- `weights`

`amount_source` controls which matrix value is used to calculate pulse amounts:

- `port`: use the exchange amount at the anchor year, then distribute that
  amount across pulse years according to the temporal weights.
- `matrix`: for each pulse year, look up the exchange amount from the matrix
  corresponding to that pulse year, then apply the temporal weight.

This distinction matters when the exchange amount itself changes over time.
`port` preserves the anchor-year amount, while `matrix` samples the target-year
matrix values.

### 5.4 Temporal Distribution Types

The distribution codes implemented in `trails/temporal_distributions.py` are:

| Code | Meaning | Behavior |
| --- | --- | --- |
| 1 | discrete | Places mass according to `loc` over integer offsets. |
| 2 | lognormal | Computes lognormal weights over positive offsets. |
| 3 | normal | Computes normal weights over the integer offset support. |
| 4 | uniform | Equal weight over all integer offsets. |
| 5 | triangular | Triangular weights with mode controlled by `loc`. |
| 6 | discrete empirical | Uses explicit `offsets` and `weights` when provided. |

For all distributions, weights are normalized to sum to 1 before being yielded.
Unknown distribution codes fall back to uniform weighting. Invalid or empty
weights produce an empty distribution.

For code 6, explicit offsets define the support when present. Tests verify that
explicit pulse vectors are used directly and are included in temporal
distribution cache keys.

### 5.5 Temporal Metadata Interpolation

Matrix values are linearly interpolated between template years when annual
interpolation is active. Temporal metadata is handled separately.

For a requested year, TRAILS first checks whether exact temporal exchange
metadata exists. If not, it attempts to interpolate temporal exchange fields
between neighboring template years when the temporal distributions are
compatible. Compatibility requires the relevant structural fields to match, such
as distribution family, support, amount source, and explicit pulse structure.

If temporal metadata is incompatible across bracketing template years, the
implementation falls back to nearest available temporal metadata rather than
constructing an invalid hybrid distribution.

## 6. Data Package Contract

TRAILS expects Frictionless datapackages with resources organized as:

```text
inventories/<model>/<pathway>/<year>/A_matrix.csv
inventories/<model>/<pathway>/<year>/B_matrix.csv
inventories/<model>/<pathway>/<year>/A_matrix_index.csv
inventories/<model>/<pathway>/<year>/B_matrix_index.csv
```

The example package uses paths such as:

```text
inventories/some model/some pathway/2005/A_matrix.csv
```

### 6.1 A Matrix CSV

The expected technosphere matrix columns include:

- `index of activity`
- `index of product`
- `value`
- `uncertainty type`
- `loc`
- `scale`
- `shape`
- `minimum`
- `maximum`
- `negative`
- `flip`

Temporal columns can also be present:

- `temporal distribution`
- `temporal loc`
- `temporal scale`
- `temporal offset min`
- `temporal offset max`
- `temporal amount source`
- explicit offset and weight fields for empirical pulse distributions

The `flip` flag is applied while loading. When set, the exchange value is
multiplied by `-1`. This is how imported or premise-style technosphere
consumption exchanges are represented as negative requirements.

### 6.2 B Matrix CSV

The biosphere matrix CSV has analogous numeric and uncertainty fields:

- `index of activity`
- `index of biosphere flow`
- `value`
- `uncertainty type`
- `loc`
- `scale`
- `shape`
- `minimum`
- `maximum`
- `negative`
- `flip`

It can also contain the same temporal distribution fields. Biosphere temporal
metadata is stored separately from technosphere temporal metadata.

### 6.3 A Matrix Index

The technosphere index file stores activity metadata. The loader expects, and
downstream tools use, fields such as:

- activity index
- name
- reference product
- unit
- location

Activity search and display functions depend on this metadata.

### 6.4 B Matrix Index

The biosphere index file stores flow metadata. The relevant fields include:

- biosphere flow index
- name
- compartment
- subcompartment
- unit

LCIA characterization and FaIR flow mapping depend on these metadata fields.

## 7. Initialization, Interpolation, and Cache

`Trails.__init__` loads `A`, `B`, temporal exchange dictionaries, and index
metadata. The primary signature is:

```python
Trails(
    package,
    interpolate_annual=True,
    cache_interpolation=True,
    interpolation_start_year_offset=-1,
    interpolation_end_year_offset=1,
    value_dtype=np.float32,
    index_dtype=np.int32,
    debug=False,
)
```

If the package is a zip archive, Frictionless unpacks it and TRAILS reports the
temporary unpacking path.

When annual interpolation is enabled:

1. Template years are sorted.
2. Sparse coordinates from neighboring years are unioned.
3. Matrix values are linearly interpolated between template years.
4. Years before the first template year receive the earliest matrix values.
5. Years after the last template year receive the latest matrix values.

The interpolation cache is stored under the platformdirs user data path:

```text
<platformdirs user data path for appname="trails", appauthor="pylca">/cache
```

The cache key includes the package descriptor, matrix file sizes and mtimes,
dtypes, and interpolation offsets. Cached entries include:

- `A.npz`
- `B.npz`
- `meta.json`
- `temporal.pkl`
- `indices.pkl`

`clear_cache()` deletes and recreates the cache directory.

## 8. Routing Workflow

Temporal routing is explicit and must be run before `lca`.

The main method is:

```python
trails.temporal_routing(
    start_year=<year>,
    start_act_idx=<activity index>,
    amount=1.0,
    max_depth=2,
    min_amount=1e-18,
    show_progress=True,
    attribute_to_roots=True,
)
```

It builds a `networkx.DiGraph` and stores it as `trails.graph`.

Each node records:

- calendar year
- activity index
- depth
- amount
- activity metadata
- frontier amount, if traversal stops at that node
- direct biosphere amount, if direct emissions from an expanded node must be
  injected later
- root attribution information, when root attribution is active

The first-tier children of the functional unit become root activities for
attribution. Deeper nodes inherit their root activity. This allows outputs to
answer questions such as: "which immediate supply-chain branch is responsible
for impacts occurring in later years?"

### 8.1 Expansion Rules

For each expanded node, TRAILS reads the year-specific row of the technosphere
matrix. For each non-production exchange:

- If no temporal distribution is defined, the child requirement occurs in the
  same year.
- If a temporal distribution is defined and `amount_source == "port"`, the
  anchor-year matrix amount is distributed across pulse years.
- If a temporal distribution is defined and `amount_source == "matrix"`, each
  pulse year samples the corresponding matrix amount from that pulse year.

Years used for matrix lookup are mapped to available scenario years:

- annual grids use identity mapping after clipping to available bounds;
- non-annual grids use the nearest available scenario year;
- out-of-range years are clipped to `min_year` or `max_year`.

Small child amounts below `min_amount` are not simply discarded. They are
recorded as frontier amounts so their residual contribution can still be solved
by the year-specific static system.

### 8.2 Frontier and Direct Biosphere

Routing terminates at frontier nodes when:

- `max_depth` is reached;
- a node has no further children;
- a child amount falls below `min_amount`;
- traversal cannot meaningfully expand further.

The frontier is later converted into per-year demand vectors. Direct biosphere
flows from expanded nodes are handled separately from frontier solves. This
prevents direct emissions of already-expanded nodes from being lost when the
remaining upstream frontier is solved by matrix methods.

### 8.3 Legacy Traversal

The project also retains `temporal_traversal`, an older dictionary-based
traversal routine. Tests still cover behavior related to traversal totals and
provenance. The graph-based `temporal_routing` workflow is the required entry
point for the current `lca` function.

## 9. LCA Execution

The high-level LCA function is exposed as both `trails.lca(...)` and
`trails.lca.lca(trails, ...)`.

The core signature is:

```python
lca(
    trails,
    methods=None,
    show_progress=True,
    attribute_to_roots=None,
    *,
    store_inventory=False,
    compute_score=True,
    ei_version="3.11",
    solver_mode="iterative",
    iterative_rtol=1e-3,
    iterative_atol=0.0,
    iterative_restart=50,
    iterative_maxiter=300,
    iterative_use_guess=True,
    iterative_preconditioner="jacobi",
    iterative_ilu_drop_tol=1e-4,
    iterative_ilu_fill_factor=10.0,
    inventory_workers=None,
)
```

`lca` requires a graph produced by `temporal_routing`; otherwise it raises a
runtime error.

If `compute_score=True`, `methods` must be provided. If `store_inventory=True`,
the function stores a full time-resolved biosphere inventory in addition to, or
instead of, characterized scores.

`attribute_to_roots=None` means: reuse the attribution mode from
`temporal_routing`. If the stored routing metadata conflicts with the requested
attribution mode, TRAILS raises and asks the caller to rerun routing.

### 9.1 Year-by-Year Solving

The LCA procedure:

1. Reads frontier and direct-biosphere nodes from `trails.graph`.
2. Aggregates frontier demands by `(year, activity)`.
3. Maps solve years to scenario years.
4. Solves the static technosphere system for each relevant year.
5. Applies biosphere exchanges from `B`.
6. Applies temporal distributions on biosphere exchanges.
7. Accumulates inventory entries and/or characterized scores on the output year
   axis.

The inventory year axis is larger than the matrix year axis. It starts at
`min_year + min_offset` and extends to `max_year + max_offset + 500`. The extra
500-year tail gives delayed emissions and long-lived response calculations room
to appear in outputs even when matrices are only available over a shorter
scenario horizon.

### 9.2 Solver Modes

`solver_mode` can be:

- `iterative`: default. Uses multi-RHS GMRES from `trails/iterative_solver.py`.
- `direct`: uses sparse direct factorization; UMFPACK when available, otherwise
  SciPy `splu`.
- `bw2calc`: builds a `bw_processing` datapackage and calls Brightway's
  `bw2calc.LCA`.

The iterative solver supports:

- relative tolerance, default `1e-3`;
- absolute tolerance;
- restart and max-iteration parameters;
- optional reuse of previous solution as initial guess;
- preconditioners: `jacobi`, `ilu`, or `none`.

Tests verify that direct and iterative modes agree with the `bw2calc` mode
within expected tolerances on the example package.

### 9.3 Multi-Root and Multi-Method Behavior

When root attribution is active, TRAILS can solve multiple right-hand sides for
the same year. This lets it preserve attribution by first-tier root activity
without repeatedly solving each root independently.

When multiple LCIA methods are requested and `store_inventory=False`, scores can
be accumulated with a method dimension directly. When inventory is stored, the
inventory is method-independent and characterization can be built afterward.

Tests confirm that:

- total scores are invariant to whether root attribution is enabled;
- multi-method scores keep a `method` dimension when inventory is not stored;
- root-attributed scoring avoids unnecessary per-root supply extraction when
  possible.

## 10. Outputs and Data Structures

### 10.1 Inventory

`trails.inventory` is an `xarray.DataArray` backed by sparse COO data.

Without root attribution, the dimensions are:

```text
activity, flow, year
```

With root attribution, the dimensions are:

```text
activity, flow, year, root activity
```

Coordinates are stored with integer-safe dtypes to avoid sparse reduction
problems on large root-attributed arrays.

### 10.2 Scores

`trails.characterized_inventory` can represent characterized scores with
different dimension combinations:

- no root attribution and no method dimension: `(activity, year)`
- root attribution: `(activity, year, root activity)`
- multi-method scoring: `(method, activity, year)`
- multi-method plus root attribution:
  `(method, activity, year, root activity)`

The exact shape depends on whether scores were computed directly, whether
inventory was stored, and whether multiple methods were requested.

### 10.3 Static Score

`trails.static_lca(...)` computes a conventional static LCA score for a given
year and activity and stores it in `trails.static_score`. It temporarily uses
inventory and characterization routines but restores any previous dynamic
inventory and characterized inventory afterward.

## 11. LCIA Characterization

LCIA methods are bundled as JSON files for ecoinvent 3.10 and 3.11.

`get_lcia_method_names(ei_version="3.11")` returns method names. Method tuple
labels are joined with `" - "` for display and selection.

Characterization matching uses biosphere metadata:

- flow name
- compartment
- subcompartment

Missing subcompartments are represented with `"unspecified"` where needed.
Matching is exact on these metadata fields. This makes the workflow
transparent, but it also means that inconsistent flow naming or compartment
metadata can lead to missing characterization factors.

`build_characterized_inventory` multiplies sparse inventory data by
characterization factors and stores the result as an xarray object. `get_cf_vector`
builds method-specific characterization vectors for direct score accumulation.

## 12. Excel Inventory Import

The importer is:

```python
import_excel_inventory(
    trails,
    path,
    *,
    year=None,
    scenario_label=None,
    cache_import=False,
)
```

It uses `bw2io.importers.excel.ExcelImporter` and intentionally skips the
`csv_drop_unknown` strategy so year-specific amount columns, such as `2030` or
`2050`, are retained.

The importer can:

- import one Excel file or a list of files;
- add or update technosphere activities;
- add or update biosphere flows;
- apply imports to all template years when no `year` or `scenario_label` is
  provided;
- apply imports to a specific mapped year or scenario label;
- interpolate year-specific amounts across the annual grid;
- clear existing temporal exchange metadata for affected rows before inserting
  imported data;
- optionally save the updated state to the interpolation cache.

Important import semantics:

- Technosphere exchanges are sign-flipped on import and stored as negative
  requirements.
- Production and biosphere exchanges are stored as given.
- Activity matching uses `(name, reference product, location)`.
- Biosphere matching uses `(name, compartment, subcompartment)`.
- Unlinked technosphere references or invalid biosphere metadata raise errors.
- Duplicate new sparse coordinates raise errors.
- `temporal_amount_source` must be `port` or `matrix`.

This importer is useful for injecting foreground inventories into an existing
scenario matrix system, but it is not a general-purpose semantic linker. The
input workbook still needs valid bw2io-style structure and consistent metadata.

## 13. Brightway Interface

`trails/bw_interface.py` converts a Trails matrix slice into a
`bw_processing` datapackage for Brightway.

Key details:

- The internal `A` orientation is converted from `(activity, product)` to
  Brightway's `(product, activity)` orientation.
- The biosphere matrix is converted to `(flow, activity)`.
- Technosphere sign handling supports:
  - `abs_flip`: nonnegative matrix values plus flip arrays;
  - `signed`: signed matrix vectors, used in fast zero-biosphere paths.
- Activity demands are mapped to reference products by finding production
  exchanges in the technosphere column.
- If the exact activity-product production exchange is unavailable, the code
  chooses the production-like entry closest to `+1` or `-1`.
- Metadata labels use exact year labels when available and nearest labels
  otherwise.

The direct and iterative solver paths avoid bw2calc overhead for repeated
year-specific solves, while the `bw2calc` path remains available as a reference
backend and compatibility mode.

## 14. FaIR Climate Response Integration

`run_fair_delta_rf` computes marginal radiative forcing and temperature
responses from a stored Trails inventory.

The public signature includes:

```python
run_fair_delta_rf(
    trails,
    *,
    scenario,
    emissions_csv=DEFAULT_EMISSIONS_CSV,
    mapping_yaml=DEFAULT_MAPPING_YAML,
    config_csv=DEFAULT_CONFIGS_CSV,
    properties_csv=DEFAULT_PROPERTIES_CSV,
    config_name=None,
    config_names=None,
    ghg_method="myhre1998",
    temperature_prescribed=False,
    scale_factor=None,
    scale_target_fraction=0.01,
    validate_emissions_delta=True,
    per_species_runs=True,
    per_species_workers=None,
    quantiles=None,
)
```

FaIR integration requires `trails.inventory` with root attribution, i.e. an
inventory with dimensions including:

```text
activity, flow, year, root activity
```

The workflow:

1. Reduces the inventory over activity to `flow, year, root activity`.
2. Maps biosphere flows to FaIR species using
   `trails/data/scenarios/fair_species_map.yaml`.
3. Filters to atmospheric flows.
4. Loads baseline IAMC-style emissions data.
5. Runs a baseline FaIR scenario.
6. Applies scaled Trails perturbations.
7. Runs perturbed FaIR simulations.
8. Stores quantiles of delta radiative forcing and delta temperature.

Default quantiles are:

```text
2.5, 25, 50, 75, 97.5
```

Outputs are stored as:

- `trails.instant_radiative_forcing`
- `trails.delta_temperature`

Both are xarray DataArrays with dimensions:

```text
quantile, year, flow, root activity
```

Per-species runs are the default. They separate positive and negative emissions
where needed so releases and removals can produce distinct response tails.

FaIR execution is protected by a global lock around `FAIR.run`, which tests
verify. This avoids unsafe concurrent calls into the FaIR model while still
allowing surrounding preparation work to use worker pools.

## 15. Plotting and Reporting

Plotting is implemented with Plotly.

Important plotting functions include:

- `plot_temporal_scores`
- `plot_adaptive_sankey`
- `plot_traversal_grid_flow`
- `plot_rf`
- `plot_temp`

`plot_temporal_scores` can plot:

- characterized inventory or direct scores;
- cumulative or annual values;
- stacked positive and negative contributions;
- root activity contributions;
- activity contributions;
- flow contributions;
- multi-method results;
- static score reference lines;
- log-scaled axes;
- trimmed year ranges.

Tests verify that the plotting code trims empty year windows and separates
positive and negative stack groups when needed.

The graph-like temporal sankey visualization can write interactive HTML and
uses a slider-based interface. Tests confirm that the generated HTML includes
the slider dependency.

`plot_rf` and `plot_temp` plot FaIR outputs by flow or root activity and support
quantile selection and uncertainty bands.

## 16. Search and Inspection Utilities

`search_activity` returns a `PrettyTable` of matching activities or biosphere
flows.

The signature is:

```python
search_activity(
    trails,
    query=None,
    *,
    name=None,
    reference_product=None,
    kind="technosphere",
    scenario_label=None,
    match="contains",
    case_sensitive=False,
)
```

For technosphere searches, results include:

- index
- name
- reference product
- location

For biosphere searches, results include flow metadata. `reference_product` is
only valid for technosphere searches.

The `Trails` object also contains table-printing helpers for exchange
inspection and traversal edge collection utilities for visualization.

## 17. Logging and Debugging

`trails/logging.py` provides `configure_trails_logging` and
`trails_log_context`.

Configured logs include contextual fields:

- `run_id`
- `case`
- `year`
- `depth`

The `Trails` object accepts `debug=True`, and many internal routines pass debug
state through to matrix loading, temporal distribution generation, traversal,
LCA, and FaIR functions.

Environment variables used for detailed flow-level debugging include:

- `TRAILS_DEBUG_FLOW_ID`
- `TRAILS_DEBUG_YEAR`
- `TRAILS_DEBUG_ACTIVITY`
- `TRAILS_DEBUG_MAX_PULSES`
- `TRAILS_DEBUG_MAX_MATCHES`

Core modules use Python's standard `logging` package.

## 18. Performance Design

The implementation is sparse-first and cache-heavy.

Major performance choices:

- Store matrices as `sparse.COO`.
- Use annual interpolation cache to avoid repeated interpolation.
- Cache temporal distribution offsets and weights.
- Cache matrix rows for repeated routing and accumulation.
- Batch inventory accumulation into sparse chunks.
- Use direct sparse factorization or iterative multi-RHS solves for repeated
  year-specific systems.
- Use root-attributed multi-RHS matrices instead of solving every root
  independently when possible.
- Use fast paths for zero-biosphere or no-temporal-distribution cases.
- Optionally use numba kernels for some accumulation paths when available.

Important implementation constraint:

- Direct and iterative solvers require square technosphere matrices for the
  solve. Invalid solver modes or preconditioners raise `ValueError`.

## 19. Test-Backed Behavioral Guarantees

The test suite covers several important invariants:

- Annual interpolation with default offsets expands an example 2005-2100 package
  to 2004-2101.
- Custom interpolation offsets expand the annual range accordingly.
- `lca` defaults to `solver_mode="iterative"` and `iterative_rtol=1e-3`.
- Direct and iterative solver outputs match the `bw2calc` backend within tested
  tolerances.
- Root-attributed totals equal non-root-attributed totals for tested cases.
- Multi-method scoring preserves method dimensions when inventory is not
  stored.
- Explicit temporal pulse vectors are used and included in cache keys.
- Imported technosphere exchanges are sign-flipped.
- Year-specific imported amounts are interpolated across annual years.
- Unlinked importer exchanges raise errors.
- Invalid solver modes and preconditioners raise errors.
- FaIR outputs have expected quantile, year, flow, and root activity dimensions.
- FaIR handles negative emissions, aggregate CO2 splits, precursor response
  species, atmospheric filtering, and name-preferred mappings in tested cases.
- FaIR suppresses all-NaN quantile warnings by treating all-NaN slices as zero.
- Plotting functions handle empty years, positive/negative stacked traces, and
  generated HTML dependencies.

## 20. Methodological Interpretation

TRAILS is best understood as a hybrid temporal LCA engine:

1. It temporalizes selected supply-chain exchanges explicitly through traversal.
2. It preserves conventional static LCA solving for unresolved frontier demand.
3. It can attribute delayed impacts back to first-tier branches of the demand.
4. It keeps matrix values scenario-year specific and can interpolate them
   annually.
5. It can pass time-resolved inventories to climate response modeling through
   FaIR.

This design avoids constructing an enormous full time-expanded technosphere
matrix. Instead, it routes temporally significant exchanges through a graph and
uses repeated static solves where conventional matrix algebra is efficient.

The tradeoff is that results depend on:

- the chosen routing depth;
- the `min_amount` cutoff;
- which exchanges have temporal distributions;
- whether temporal exchange amounts use `port` or `matrix` semantics;
- the scenario-year interpolation assumptions;
- exact metadata matching for LCIA and climate species mapping.

These choices should be explicit in any manuscript methods section.

## 21. Suggested Methodology Wording

A concise description of the implemented method could be:

> TRAILS represents scenario-dependent life cycle inventories as sparse
> technosphere and biosphere matrix tensors indexed by scenario year. Scenario
> matrices can be linearly interpolated to annual resolution. Starting from a
> functional unit, TRAILS builds a temporal traversal graph by applying
> exchange-specific temporal distributions to technosphere links. Traversal
> stops at a configurable depth or cutoff threshold, producing year-specific
> frontier demands. These frontier demands are solved with conventional
> year-specific LCA systems, while direct biosphere flows encountered during
> traversal are injected into the inventory. Biosphere exchanges can themselves
> be temporally distributed. Results are accumulated as sparse time-resolved
> inventories or characterized scores, optionally attributed to first-tier root
> activities.

The limitations should also be stated:

> The implementation uses discrete annual time steps, deterministic matrix
> values, and exact metadata matching for characterization. It does not solve a
> global time-expanded technosphere matrix and does not perform uncertainty
> sampling from exchange uncertainty fields.

## 22. Practical Usage Pattern

A typical workflow is:

```python
from frictionless import Package
from trails import Trails, get_lcia_method_names

package = Package("path/to/datapackage.json")
tr = Trails(package, interpolate_annual=True)

tr.temporal_routing(
    start_year=2020,
    start_act_idx=123,
    amount=1.0,
    max_depth=5,
    attribute_to_roots=True,
)

methods = [get_lcia_method_names("3.11")[0]]
tr.lca(
    methods=methods,
    store_inventory=True,
    compute_score=True,
    solver_mode="iterative",
)
```

For FaIR response calculations:

```python
from trails.fair_rf import run_fair_delta_rf

run_fair_delta_rf(
    tr,
    scenario="SSP2-4.5",
)
```

This second step requires `tr.inventory` with root attribution.

## 23. Key Assumptions to Document in Analyses

When reporting results produced by TRAILS, record:

- input datapackage source and scenario years;
- whether annual interpolation was enabled;
- interpolation start and end offsets;
- functional unit activity index, amount, and start year;
- `max_depth` and `min_amount`;
- whether root attribution was enabled;
- solver mode and solver tolerances;
- LCIA method names and ecoinvent characterization version;
- whether full inventory was stored or only scores were computed;
- whether Excel inventories were imported after package loading;
- whether FaIR was run, including baseline scenario, configs, scaling, and
  quantiles;
- any custom flow-to-species mapping or characterization data.

## 24. Current Implementation Cautions

The following points are useful when developing or interpreting results:

- Because temporal routing is graph-based and cutoff-dependent, increasing
  `max_depth` or lowering `min_amount` can change results and runtime.
- Temporally distributed biosphere exchanges can shift emissions outside the
  original matrix-year range; the inventory axis is intentionally extended.
- Matrix lookup years are clipped to available scenario years, so very late
  emissions use the last available matrix values where matrix lookup is needed.
- LCIA and FaIR mappings are metadata-sensitive. Flow naming and compartment
  consistency matter.
- Annual interpolation is linear for matrix values; this is an implementation
  assumption, not a physical model by itself.
- Imported Excel inventories mutate the in-memory matrices and temporal
  exchange dictionaries.
- Cached interpolation can speed up repeated work, but cache keys depend on
  file metadata and interpolation settings. Use `clear_cache()` when validating
  changes to package loading or interpolation logic.
- Debug logging can be verbose and should be configured intentionally for large
  runs.

## 25. Glossary

- `A matrix`: technosphere matrix tensor storing production and inter-activity
  exchanges by scenario year.
- `B matrix`: biosphere matrix tensor storing elementary flows by scenario
  year.
- `template year`: an original year present in the input datapackage.
- `scenario year`: a year label available in the active Trails matrix grid,
  either original or interpolated.
- `temporal exchange`: metadata attached to an exchange describing how its
  amount is distributed across year offsets.
- `port amount`: temporal amount calculated from the anchor-year matrix value.
- `matrix amount`: temporal amount calculated separately from each pulse-year
  matrix value.
- `frontier`: unresolved demand at the edge of the explicit traversal graph,
  later solved with a static LCA system.
- `root activity`: first-tier child activity of the functional unit used for
  attribution.
- `direct biosphere`: biosphere exchanges from expanded traversal nodes that
  are injected directly rather than recovered only through frontier solves.
- `characterized inventory`: inventory multiplied by LCIA characterization
  factors.
- `delta radiative forcing`: FaIR-estimated radiative forcing perturbation
  caused by the Trails inventory relative to a baseline scenario.
- `delta temperature`: FaIR-estimated temperature perturbation caused by the
  Trails inventory relative to a baseline scenario.
