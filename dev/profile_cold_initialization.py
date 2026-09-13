"""Profile Trails initialization with an isolated interpolation cache.

Run from the repository root with PYTHONPATH=. to select that checkout.
Input packages are supplied separately; the benchmark never removes user caches.
"""

import argparse
import cProfile
import hashlib
import importlib
import json
from pathlib import Path
import pstats
import shutil
import time

parser = argparse.ArgumentParser()
parser.add_argument("--package", type=Path, required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--profile", action="store_true")
parser.add_argument("--reuse-cache", type=Path)
parser.add_argument("--datapackage-source", type=Path)
args = parser.parse_args()
out = Path(args.out).resolve()
out.mkdir(parents=True, exist_ok=False)
from datapackage import Package
import trails

core = importlib.import_module("trails.trails")
cache = importlib.import_module("trails.cache_interpolation")
original_cache_dir = cache.cache_dir_for_package


def isolated_cache_dir(*a, **kw):
    return (args.reuse_cache or out / "cache") / original_cache_dir(*a, **kw).name


cache.cache_dir_for_package = isolated_cache_dir
if args.datapackage_source:
    spec = importlib.util.spec_from_file_location(
        "trails._profile_datapackage", args.datapackage_source
    )
    source_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source_module)
    for name in [
        "load_matrices_from_package",
        "load_indices_from_package",
        "interpolate_to_annual",
    ]:
        setattr(core, name, getattr(source_module, name))
else:
    source_module = importlib.import_module("trails.datapackage")
stages = {}
for name in [
    "load_cached_interpolation",
    "load_matrices_from_package",
    "load_indices_from_package",
    "interpolate_to_annual",
    "save_cached_interpolation",
]:
    original = getattr(core, name)

    def timed(*a, _fn=original, _name=name, **kw):
        start = time.perf_counter()
        result = _fn(*a, **kw)
        stages[_name] = time.perf_counter() - start
        print("STAGE", _name, stages[_name], flush=True)
        return result

    setattr(core, name, timed)
package_path = args.package.resolve()
prof = cProfile.Profile()
if args.profile:
    prof.enable()
start = time.perf_counter()
package = Package(str(package_path))
stages["Package"] = time.perf_counter() - start
print("STAGE Package", stages["Package"], flush=True)
init_start = time.perf_counter()
model = trails.Trails(package)
stages["Trails"] = time.perf_counter() - init_start
stages["total"] = time.perf_counter() - start
if args.profile:
    prof.disable()
    prof.dump_stats(str(out / "init.prof"))
    with (out / "profile.txt").open("w") as f:
        pstats.Stats(prof, stream=f).strip_dirs().sort_stats("cumulative").print_stats(
            70
        )


def digest(x):
    return hashlib.sha256(x.tobytes()).hexdigest()


summary = dict(
    package_base_path=package.base_path,
    datapackage_source=source_module.__file__,
    datapackage_sha256=hashlib.sha256(
        Path(source_module.__file__).read_bytes()
    ).hexdigest(),
    profiled=args.profile,
    cold=args.reuse_cache is None,
    stages=stages,
    version=trails.__version__,
    source=trails.__file__,
    matrices={
        name: dict(
            shape=list(m.shape),
            nnz=m.nnz,
            dtype=str(m.dtype),
            coords=digest(m.coords),
            data=digest(m.data),
        )
        for name, m in [("A", model.A), ("B", model.B)]
    },
    labels=model.scenario_labels,
    cache_bytes=sum(
        p.stat().st_size
        for p in isolated_cache_dir(
            package,
            value_dtype=str(model.value_dtype),
            index_dtype=str(model.index_dtype),
        ).rglob("*")
        if p.is_file()
    ),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
model.close()
# This Package was created from a ZIP by this benchmark; release its extraction.
extracted = getattr(package, "_Package__tempdir", None)
if extracted:
    shutil.rmtree(extracted)
