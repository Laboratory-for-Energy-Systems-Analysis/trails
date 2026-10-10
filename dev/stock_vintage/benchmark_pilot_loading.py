"""Measure isolated cold/warm loading and routing of real pilot packages.

Full-LCI verifier wall times are reported separately from existing log/report
timestamps. They include verification and concurrent workload, not just solving.
This benchmark does not rerun the expensive all-flow comparisons.
"""

import argparse
from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import zipfile

PILOTS = [
    "ccgt-real-npi",
    "pv-real-primary",
    "water-real-primary",
    "trucks-real-primary",
    "bev-real-primary",
]


def worker(root, pilot, variant, mode):
    from datapackage import Package
    import numpy as np
    from trails import Trails

    export = json.loads((root / pilot / "export-evidence.json").read_text())
    package = export["packages"][variant]
    path = Path(package["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != package["sha256"]:
        raise ValueError("Changed pilot package")
    log = io.StringIO()
    start = time.perf_counter()
    with redirect_stdout(log):
        model = Trails(
            Package(str(path)),
            stock_vintage=variant == "corrected",
            interpolate_annual=True,
            cache_interpolation=mode == "warm",
            value_dtype=np.float64,
        )
    loaded = time.perf_counter()
    caller = export["rewrite"][
        "caller" if variant == "corrected" else "original_caller"
    ]
    with model:
        indices = model.activity_indices[model.template_labels[0]]
        activity = next(
            i
            for i, d in indices.items()
            if all(
                d[k] == caller[k]
                for k in ("name", "reference product", "unit", "location")
            )
        )
        model.temporal_routing(
            start_year=2030,
            start_act_idx=activity,
            amount=1,
            max_depth=2,
            adaptive_relative_score_cutoff=None,
            min_amount=0,
            attribute_to_roots=True,
            show_progress=False,
        )
        routed = time.perf_counter()
        with zipfile.ZipFile(path) as archive:
            profile_bytes = sum(
                i.file_size
                for i in archive.infolist()
                if i.filename == "stock_vintage.json"
            )
        return {
            "pilot": pilot,
            "variant": variant,
            "requested_cache_mode": mode,
            "cache_loaded": "Loaded interpolated matrices from cache:"
            in log.getvalue(),
            "load_seconds": loaded - start,
            "route_seconds": routed - loaded,
            "process_peak_rss_mib_through_routing": resource.getrusage(
                resource.RUSAGE_SELF
            ).ru_maxrss
            / (1048576 if sys.platform == "darwin" else 1024),
            "activities": model.A.shape[1],
            "graph_nodes": model.graph.number_of_nodes(),
            "graph_edges": model.graph.number_of_edges(),
            "zip_bytes": path.stat().st_size,
            "annual_resource_uncompressed_bytes": profile_bytes,
            "package_sha256": package["sha256"],
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("research_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--worker", choices=PILOTS)
    parser.add_argument("--variant", choices=["legacy", "corrected"])
    parser.add_argument("--mode", choices=["cold", "warm"])
    args = parser.parse_args()
    if args.worker:
        print(
            json.dumps(worker(args.research_dir, args.worker, args.variant, args.mode))
        )
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for pilot in PILOTS:
        for variant in ("legacy", "corrected"):
            for mode in ("cold", "warm"):
                command = [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    str(args.research_dir),
                    "--output-dir",
                    str(args.output_dir),
                    "--worker",
                    pilot,
                    "--variant",
                    variant,
                    "--mode",
                    mode,
                ]
                result = json.loads(subprocess.check_output(command, text=True))
                if mode == "warm" and not result["cache_loaded"]:
                    # A separate process builds a missing cache before timing
                    # an actual warm load in another fresh process.
                    result = json.loads(subprocess.check_output(command, text=True))
                    if not result["cache_loaded"]:
                        raise AssertionError(
                            "Warm-cache benchmark did not reload cache"
                        )
                rows.append(result)
                print(json.dumps(result), flush=True)
    timings = []
    for pilot in PILOTS:
        for path in sorted((args.research_dir / pilot).glob("checks*.json")):
            log = path.with_suffix(".log")
            if not log.exists() or not hasattr(log.stat(), "st_birthtime"):
                continue
            report = json.loads(path.read_text())
            nested = report.get("reports", {"single": report})
            if not any(v.get("lci_checks") for v in nested.values()):
                continue
            timings.append(
                {
                    "pilot": pilot,
                    "report": path.name,
                    "report_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "log_birth_to_report_write_seconds": path.stat().st_mtime
                    - log.stat().st_birthtime,
                }
            )
    with (args.output_dir / "load-route-benchmark.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "resource-evidence.json").write_text(
        json.dumps(
            {
                "scope": "Fresh processes, annual interpolation, service year 2030, demand +1, fixed depth 2. RSS includes loading/routing, not full LCI. Cold means matrix cache disabled, not flushed OS disk cache.",
                "rights": "Restricted local package-derived evidence",
                "runs": rows,
                "full_verifier_observed_wall_times": timings,
                "wall_time_limitations": "Filesystem timestamp intervals include full signed LCI, independent verification and concurrent system load; no isolated native solve or full-LCI peak RSS claim.",
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
