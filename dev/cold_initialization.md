# Cold datapackage initialization

`Trails(Package(...))` can spend more than a minute preparing annual matrices before it writes its first interpolation cache. For the BrightCon REMIND SSP2-PkBudg1000 DACCS background package, this change reduces median initialization time from **112.33 s to 45.94 s (59.1%)**. Matrix data, coordinates and cached temporal/activity metadata remain identical.

## What changed

- Read the biosphere columns consumed by inventory loading, while retaining required-header validation and the full reader's default behavior.
- Parse numeric strings directly instead of stripping every cell first. Retry failed conversions with whitespace stripped, including Unicode whitespace.
- Exclude zero and NaN temporal-distribution entries before parsing exchanges.
- Reuse parsed explicit pulse definitions within a single load. Bound this temporary lookup to 4,096 entries and return independent lists for each exchange.
- Align the sparse coordinate union once per pair of scenario anchors and reuse it across their intervening annual slices. Preserve float64 intermediate arithmetic, output dtypes and exact-zero pruning. Searching the union for input keys also handles empty and disjoint supports safely.

The constructor defaults, annual resolution, cache format and cache identity are unchanged. A running notebook needs a kernel restart to import the updated code.

## Measurements

Two baseline and two final runs used fresh Python processes and separate empty interpolation-cache directories. Imports, output hashing and tests are outside the measured expression.

| Stage | Baseline median | Final median |
|---|---:|---:|
| Open and extract ZIP | 3.32 s | 3.31 s |
| Load matrices and temporal exchanges | 61.24 s | 23.88 s |
| Load activity/flow metadata | 2.33 s | 2.34 s |
| Interpolate annual matrices | 41.29 s | 12.42 s |
| Write cache | 4.16 s | 3.99 s |
| **Trails constructor** | **109.02 s** | **42.63 s** |
| **Complete expression, including Package** | **112.33 s** | **45.94 s** |

The constructor includes the loading, interpolation and writing stages. Individual complete-expression times were **111.88 / 112.78 s** before and **46.89 / 44.99 s** after the change. The cache remains **4,117,310,945 bytes**. The speedup comes mainly from preparing the matrices, rather than from serialization.

An initial cProfile run attributed 101.0 s to loading, 44.5 s to interpolation and 4.2 s to serialization. That instrumented run took 161.7 s in total; profiler overhead in the original string-processing loops makes it unsuitable as the speedup baseline.

Measurements used TRAILS 1.1.0, Python 3.11.14 and macOS 26.6.2 on ARM64. The complete input ZIP has SHA-256 `90ba5671a2d1c7e08614095b5eb2f7721d021fc1bde02099b96ee7c90b062078`. The original report generator ran under a separate Python interpreter; the version above is the benchmark interpreter's version.

“Cold” means no TRAILS interpolation cache. Operating-system file caches and installed JIT caches were not flushed. Baseline and final benchmarks ran serially; brief focused tests overlapped portions of baseline runs (about 2–3 seconds). Final runs had no concurrent test job. Results describe this package and machine, rather than a guaranteed improvement for every input.

## Correctness and regression coverage

[The benchmark evidence](cold_initialization_benchmark.json) records individual timings and complete matrix hashes. Across all four cold runs:

- Coordinates, values, shapes, dtypes, nonzero counts and year labels match for both complete 98-year matrices.
- `temporal.pkl`, `indices.pkl` and `meta.json` are byte-for-byte identical.
- The optimized code reads an old-format baseline cache without rebuilding it. That separate warm-cache check took 6.38 s including opening the ZIP, and produced matching matrices and labels.

The warm run's original `cache_bytes` field is zero because the early harness counted only its new output directory; the existing cache was elsewhere. This field is not used to infer a size change. The committed harness counts the selected cache directory for both cold and warm runs.

The 52 focused tests cover datapackage loading, interpolation/cache identity and round trips, temporal distributions and Trails accessors. New regressions include float32/float64 interpolation, empty/changing supports, cancellation, unordered anchor years and boundary padding, malformed/blank inputs, Unicode whitespace, required headers and independent mutable pulse lists. Full downstream numerical calculations are not required for the performance comparison because the complete initialization outputs match directly.

```sh
python -m pytest tests/test_datapackage.py tests/test_cache_interpolation.py tests/test_cold_initialization.py tests/test_trails.py tests/test_temporal_distributions.py -q
```

## Reproduce

Use [the profiling helper](profile_cold_initialization.py) from the repository root with a compatible environment and a locally available datapackage. The large conference input is supplied separately; no inventory data or generated caches are committed. Each `--out` directory must be new.

```sh
PYTHONPATH=. python dev/profile_cold_initialization.py --package /path/to/trails_remind_SSP2-PkBudg1000.zip --out /tmp/trails-final
```

To compare with the original loader without changing your checkout, extract its module from the baseline commit:

```sh
git show df831a4:trails/datapackage.py > /tmp/trails-datapackage-baseline.py
PYTHONPATH=. python dev/profile_cold_initialization.py --package /path/to/trails_remind_SSP2-PkBudg1000.zip --out /tmp/trails-baseline --datapackage-source /tmp/trails-datapackage-baseline.py
```

Add `--profile` for cProfile data and a sorted text report. Use `--reuse-cache /tmp/trails-baseline/cache` with a new output directory for the compatibility check. This warm run is separate from the cold performance comparison.

The helper isolates interpolation caches, records stage timings and full matrix hashes, and removes its own ZIP extraction after a successful run. It retains generated cache files for comparisons. Remove these disposable benchmark caches after saving the evidence; do not clear an existing user cache to obtain a cold measurement.
