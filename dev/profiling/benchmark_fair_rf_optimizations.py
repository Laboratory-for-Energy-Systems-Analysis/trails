from __future__ import annotations

import argparse
import json
import resource
import statistics
import time
from pathlib import Path
from typing import Any

import numpy as np
from datapackage import Package

from trails import Trails
from trails.fair_rf import run_fair_delta_rf

HERE = Path(__file__).resolve().parent
DEFAULT_DP = "/Users/romain/GitHub/premise/dev/trails_2026-02-22.zip"
DEFAULT_INVENTORY = str(HERE.parent / "lci-case-study-daccs_storage_risk.xlsx")
DEFAULT_METHOD = (
    "IPCC 2021 - climate change: total (excl. biogenic CO2) - "
    "global warming potential (GWP100)"
)
DEFAULT_SCENARIO = "REMIND|SSP2-PkBudg650"


def _q50_vector(rf: Any) -> np.ndarray:
    q50 = rf.sel(quantile=50, method="nearest").sum(dim=["flow", "root activity"])
    data = q50.data
    if hasattr(data, "todense"):
        return np.asarray(data.todense(), dtype=float).ravel()
    return np.asarray(data, dtype=float).ravel()


def _run_once(args: argparse.Namespace) -> dict[str, Any]:
    dp = Package(args.datapackage)
    trails = Trails(
        dp,
        interpolate_annual=True,
        debug=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )
    trails.import_excel_inventory(args.inventory)

    t0 = time.perf_counter()
    trails.temporal_routing(
        start_year=2025,
        start_act_idx=41792,
        amount=1,
        max_depth=args.max_depth,
        show_progress=False,
        attribute_to_roots=True,
    )
    routing_s = time.perf_counter() - t0

    t1 = time.perf_counter()
    trails.lca(
        methods=[DEFAULT_METHOD],
        show_progress=False,
        compute_score=True,
        store_inventory=True,
        solver_mode="iterative",
        iterative_rtol=1e-4,
    )
    lca_s = time.perf_counter() - t1

    rf_kwargs: dict[str, Any] = {
        "scenario": args.scenario,
        "per_species_runs": args.per_species_runs,
    }
    if args.per_species_workers is not None:
        rf_kwargs["per_species_workers"] = args.per_species_workers
    if args.config_name is not None:
        rf_kwargs["config_name"] = args.config_name

    t2 = time.perf_counter()
    rf = run_fair_delta_rf(trails, **rf_kwargs)
    fair_s = time.perf_counter() - t2

    q50 = _q50_vector(rf)
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "routing_s": routing_s,
        "lca_s": lca_s,
        "fair_s": fair_s,
        "shape": list(map(int, rf.shape)),
        "nnz": int(getattr(rf.data, "nnz", 0)),
        "sum": float(rf.sum().item()),
        "q50_sum": float(np.sum(q50)),
        "q50_l2": float(np.sqrt(np.sum(q50 * q50))),
        "q50": q50.tolist(),
        "maxrss": int(maxrss),
    }


def _median(values: list[float]) -> float:
    return float(statistics.median(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default="run")
    parser.add_argument("--datapackage", default=DEFAULT_DP)
    parser.add_argument("--inventory", default=DEFAULT_INVENTORY)
    parser.add_argument("--scenario", default=DEFAULT_SCENARIO)
    parser.add_argument("--max-depth", type=int, default=2)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--per-species-runs", action="store_true", default=True)
    parser.add_argument(
        "--no-per-species-runs", dest="per_species_runs", action="store_false"
    )
    parser.add_argument("--per-species-workers", type=int, default=None)
    parser.add_argument("--config-name", default=None)
    args = parser.parse_args()

    runs: list[dict[str, Any]] = []
    for _ in range(args.repeat):
        runs.append(_run_once(args))

    summary = {
        "label": args.label,
        "depth": args.max_depth,
        "repeat": args.repeat,
        "scenario": args.scenario,
        "per_species_runs": args.per_species_runs,
        "per_species_workers": args.per_species_workers,
        "config_name": args.config_name,
        "median_routing_s": _median([r["routing_s"] for r in runs]),
        "median_lca_s": _median([r["lca_s"] for r in runs]),
        "median_fair_s": _median([r["fair_s"] for r in runs]),
        "median_maxrss": int(statistics.median([r["maxrss"] for r in runs])),
        "runs": runs,
    }
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
