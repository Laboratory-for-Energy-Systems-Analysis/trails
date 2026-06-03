from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from datapackage import Package
from scipy.sparse.linalg import spsolve

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from trails import Trails, get_lcia_method_names
from trails.characterization import get_cf_vector
from trails.lca import (
    _build_direct_technosphere_for_year,
    _reference_product_from_activity_direct,
)


DEFAULT_DATAPACKAGE = (
    REPO_ROOT / "dev" / "trails_remind_SSP2-PkBudg1000.zip"
)
DEFAULT_LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")
DEFAULT_METHOD = "EF v3.1 - particulate matter formation - impact on human health"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "dev" / "notebook_runs" / "static_activity_intensities"


def _activity_metadata(trails: Trails, activity_index: int, year: int) -> dict[str, Any]:
    labels: list[str] = []
    try:
        labels.append(str(trails._map_year_to_scenario_year(int(year))))
    except Exception:
        pass
    labels.extend(str(label) for label in getattr(trails, "scenario_labels", []))

    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        mapping = getattr(trails, "activity_indices", {}).get(label, {})
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)
    return {}


def _load_trails(args: argparse.Namespace) -> tuple[Trails, float, float]:
    start = time.perf_counter()
    trails = Trails(
        Package(str(Path(args.datapackage).expanduser().resolve())),
        interpolate_annual=True,
        cache_interpolation=not bool(args.no_cache_interpolation),
        interpolation_start_year_offset=int(args.interpolation_start_year_offset),
        interpolation_end_year_offset=int(args.interpolation_end_year_offset),
    )
    load_seconds = time.perf_counter() - start

    import_seconds = 0.0
    inventory_paths = [Path(path).expanduser().resolve() for path in args.inventories]
    if inventory_paths:
        missing = [path for path in inventory_paths if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing inventory file(s):\n- "
                + "\n- ".join(str(path) for path in missing)
            )
        start = time.perf_counter()
        trails.import_excel_inventory([str(path) for path in inventory_paths])
        import_seconds = time.perf_counter() - start

    return trails, load_seconds, import_seconds


def _reference_product_arrays(A_csc: Any) -> tuple[np.ndarray, np.ndarray]:
    n_activities = int(A_csc.shape[1])
    product_indices = np.full(n_activities, -1, dtype=np.int64)
    production_values = np.full(n_activities, np.nan, dtype=np.float64)

    indptr = A_csc.indptr
    indices = A_csc.indices
    data = np.asarray(A_csc.data, dtype=np.float64)
    for activity_index in range(n_activities):
        start = int(indptr[activity_index])
        end = int(indptr[activity_index + 1])
        rows = indices[start:end]
        vals = data[start:end]
        if rows.size == 0:
            continue

        product_index: int | None = None
        production_value = 0.0
        if activity_index < int(A_csc.shape[0]):
            matches = np.where(rows == activity_index)[0]
            if matches.size:
                pos = int(matches[0])
                product_index = int(activity_index)
                production_value = float(vals[pos])

        if product_index is None or production_value == 0.0:
            pos = int(np.argmin(np.abs(np.abs(vals) - 1.0)))
            product_index = int(rows[pos])
            production_value = float(vals[pos])

        product_indices[activity_index] = int(product_index)
        production_values[activity_index] = float(production_value)

    return product_indices, production_values


def _solve_adjoint(
    *,
    trails: Trails,
    year: int,
    method: str,
    ei_version: str,
) -> tuple[np.ndarray, np.ndarray, Any, int, dict[str, float]]:
    timings: dict[str, float] = {}

    start = time.perf_counter()
    context = trails._get_scenario_context(int(year))
    if context is None:
        raise RuntimeError(f"No scenario context available for year={int(year)}.")
    scenario_year, _label, t = context
    timings["map_year_seconds"] = time.perf_counter() - start

    start = time.perf_counter()
    A_csc, _product_dict, _ref_cache = _build_direct_technosphere_for_year(
        trails=trails,
        year=int(scenario_year),
        cache={},
    )
    timings["build_technosphere_seconds"] = time.perf_counter() - start

    start = time.perf_counter()
    cf = get_cf_vector(
        trails=trails,
        methods=[method],
        char_cache={},
        ei_version=str(ei_version),
    )
    timings["build_cf_seconds"] = time.perf_counter() - start

    start = time.perf_counter()
    B_t = trails.B[int(t), :, :]
    direct_scores = np.asarray(B_t @ cf, dtype=np.float64).reshape(-1)
    timings["characterize_direct_biosphere_seconds"] = time.perf_counter() - start

    start = time.perf_counter()
    intensities = np.asarray(spsolve(A_csc.T.tocsc(), direct_scores), dtype=np.float64)
    timings["adjoint_solve_seconds"] = time.perf_counter() - start

    return intensities, direct_scores, A_csc, int(scenario_year), timings


def _activity_scores_from_product_intensities(
    intensities: np.ndarray,
    product_indices: np.ndarray,
    production_values: np.ndarray,
) -> np.ndarray:
    scores = np.full(product_indices.shape, np.nan, dtype=np.float64)
    valid = (product_indices >= 0) & (product_indices < intensities.size)
    signs = np.where(production_values[valid] < 0.0, -1.0, 1.0)
    scores[valid] = signs * intensities[product_indices[valid]]
    return scores


def _write_scores(
    *,
    trails: Trails,
    year: int,
    method: str,
    output_csv: Path,
    intensities: np.ndarray,
    direct_scores: np.ndarray,
    A_csc: Any,
) -> float:
    start = time.perf_counter()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    ref_cache: dict[int, tuple[int, float]] = {}
    n_activities = int(A_csc.shape[1])

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "activity_index",
                "name",
                "reference_product",
                "location",
                "unit",
                "year",
                "method",
                "reference_product_index",
                "production_exchange_value",
                "score_per_reference_product_demand",
                "score_per_product_index_demand",
                "direct_biosphere_score_per_activity_supply",
            ],
        )
        writer.writeheader()
        for activity_index in range(n_activities):
            meta = _activity_metadata(trails, activity_index, int(year))
            try:
                product_index, production_value = _reference_product_from_activity_direct(
                    A_csc=A_csc,
                    activity_id=activity_index,
                    cache=ref_cache,
                )
                sign = -1.0 if production_value < 0.0 else 1.0
                product_score = float(intensities[int(product_index)])
                score = sign * product_score
            except Exception:
                product_index = ""
                production_value = ""
                product_score = float("nan")
                score = float("nan")

            writer.writerow(
                {
                    "activity_index": activity_index,
                    "name": meta.get("name", ""),
                    "reference_product": meta.get("reference product", ""),
                    "location": meta.get("location", ""),
                    "unit": meta.get("unit", ""),
                    "year": int(year),
                    "method": method,
                    "reference_product_index": product_index,
                    "production_exchange_value": production_value,
                    "score_per_reference_product_demand": score,
                    "score_per_product_index_demand": product_score,
                    "direct_biosphere_score_per_activity_supply": float(
                        direct_scores[activity_index]
                    )
                    if activity_index < direct_scores.size
                    else "",
                }
            )
    return time.perf_counter() - start


def _write_activity_metadata(trails: Trails, *, year: int, output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    n_activities = int(trails.A.shape[1])
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "activity_index",
                "name",
                "reference_product",
                "location",
                "unit",
            ],
        )
        writer.writeheader()
        for activity_index in range(n_activities):
            meta = _activity_metadata(trails, activity_index, int(year))
            writer.writerow(
                {
                    "activity_index": int(activity_index),
                    "name": meta.get("name", ""),
                    "reference_product": meta.get("reference product", ""),
                    "location": meta.get("location", ""),
                    "unit": meta.get("unit", ""),
                }
            )


def _all_year_output_paths(args: argparse.Namespace) -> dict[str, Path]:
    output_csv = Path(args.output_csv).expanduser().resolve()
    output_json = Path(args.output_json).expanduser().resolve()
    stem = output_csv.with_suffix("")
    return {
        "scores_csv": output_csv,
        "summary_json": output_json,
        "matrix_npz": Path(args.output_npz).expanduser().resolve()
        if args.output_npz is not None
        else stem.with_suffix(".npz"),
        "timings_csv": Path(args.timings_csv).expanduser().resolve()
        if args.timings_csv is not None
        else Path(str(stem) + "_timings.csv"),
        "metadata_csv": Path(args.metadata_csv).expanduser().resolve()
        if args.metadata_csv is not None
        else Path(str(stem) + "_activity_metadata.csv"),
    }


def _run_all_years(args: argparse.Namespace) -> int:
    datapackage = Path(args.datapackage).expanduser().resolve()
    if not datapackage.exists():
        raise FileNotFoundError(datapackage)

    lcia_json = None if args.lcia_json is None else Path(args.lcia_json).expanduser()
    if lcia_json is not None:
        if not lcia_json.exists():
            raise FileNotFoundError(lcia_json)
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json.resolve())

    available = get_lcia_method_names(ei_version=str(args.ei_version))
    if args.method not in available:
        raise ValueError(
            f"LCIA method not found for ecoinvent {args.ei_version}: {args.method}"
        )

    paths = _all_year_output_paths(args)
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Datapackage: {datapackage}", flush=True)
    print("Years: all interpolated scenario labels", flush=True)
    print(f"Method: {args.method}", flush=True)

    total_start = time.perf_counter()
    trails, load_seconds, import_seconds = _load_trails(args)

    years = np.asarray(
        sorted(int(label) for label in trails.scenario_labels),
        dtype=np.int64,
    )
    n_years = int(years.size)
    n_activities = int(trails.A.shape[1])
    scores_matrix = np.empty((n_years, n_activities), dtype=np.float64)
    direct_scores_matrix = np.empty((n_years, n_activities), dtype=np.float64)

    cf_start = time.perf_counter()
    cf = get_cf_vector(
        trails=trails,
        methods=[args.method],
        char_cache={},
        ei_version=str(args.ei_version),
    )
    build_cf_seconds = time.perf_counter() - cf_start

    _write_activity_metadata(
        trails,
        year=int(years[0]),
        output_csv=paths["metadata_csv"],
    )

    timings_rows: list[dict[str, float | int]] = []
    output_rows = 0
    year_loop_start = time.perf_counter()
    with paths["scores_csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "year",
                "activity_index",
                "score_per_reference_product_demand",
                "direct_biosphere_score_per_activity_supply",
            ]
        )

        for year_position, year in enumerate(years):
            year_total_start = time.perf_counter()

            start = time.perf_counter()
            A_csc, _product_dict, _ref_cache = _build_direct_technosphere_for_year(
                trails=trails,
                year=int(year),
                cache={},
            )
            build_technosphere_seconds = time.perf_counter() - start

            start = time.perf_counter()
            context = trails._get_scenario_context(int(year))
            if context is None:
                raise RuntimeError(f"No scenario context for year {int(year)}.")
            _scenario_year, _label, t = context
            B_t = trails.B[int(t), :, :]
            direct_scores = np.asarray(B_t @ cf, dtype=np.float64).reshape(-1)
            characterize_direct_biosphere_seconds = time.perf_counter() - start

            start = time.perf_counter()
            intensities = np.asarray(
                spsolve(A_csc.T.tocsc(), direct_scores),
                dtype=np.float64,
            )
            adjoint_solve_seconds = time.perf_counter() - start

            start = time.perf_counter()
            product_indices, production_values = _reference_product_arrays(A_csc)
            scores = _activity_scores_from_product_intensities(
                intensities,
                product_indices,
                production_values,
            )
            map_reference_products_seconds = time.perf_counter() - start

            start = time.perf_counter()
            scores_matrix[year_position, :] = scores
            direct_scores_matrix[year_position, :] = direct_scores
            for activity_index in range(n_activities):
                writer.writerow(
                    [
                        int(year),
                        int(activity_index),
                        repr(float(scores[activity_index])),
                        repr(float(direct_scores[activity_index])),
                    ]
                )
            output_rows += n_activities
            write_year_seconds = time.perf_counter() - start
            year_total_seconds = time.perf_counter() - year_total_start

            timings_rows.append(
                {
                    "year": int(year),
                    "technosphere_nnz": int(A_csc.nnz),
                    "biosphere_nnz": int(getattr(B_t, "nnz", 0)),
                    "build_technosphere_seconds": build_technosphere_seconds,
                    "characterize_direct_biosphere_seconds": (
                        characterize_direct_biosphere_seconds
                    ),
                    "adjoint_solve_seconds": adjoint_solve_seconds,
                    "map_reference_products_seconds": (
                        map_reference_products_seconds
                    ),
                    "write_year_csv_seconds": write_year_seconds,
                    "year_total_seconds": year_total_seconds,
                }
            )
            print(
                f"{int(year)}: total={year_total_seconds:.2f}s, "
                f"solve={adjoint_solve_seconds:.2f}s",
                flush=True,
            )

    year_loop_and_scores_csv_seconds = time.perf_counter() - year_loop_start

    with paths["timings_csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timings_rows[0].keys()))
        writer.writeheader()
        writer.writerows(timings_rows)

    np.savez_compressed(
        paths["matrix_npz"],
        years=years,
        activity_indices=np.arange(n_activities, dtype=np.int64),
        score_per_reference_product_demand=scores_matrix,
        direct_biosphere_score_per_activity_supply=direct_scores_matrix,
    )
    total_seconds = time.perf_counter() - total_start

    solve_times = np.asarray(
        [row["adjoint_solve_seconds"] for row in timings_rows],
        dtype=np.float64,
    )
    year_times = np.asarray(
        [row["year_total_seconds"] for row in timings_rows],
        dtype=np.float64,
    )
    write_year_times = np.asarray(
        [row["write_year_csv_seconds"] for row in timings_rows],
        dtype=np.float64,
    )
    summary = {
        "datapackage": str(datapackage),
        "inventories": [
            str(Path(path).expanduser().resolve()) for path in args.inventories
        ],
        "method": args.method,
        "ei_version": str(args.ei_version),
        "n_years": n_years,
        "year_min": int(years.min()),
        "year_max": int(years.max()),
        "n_activities": n_activities,
        "output_rows": int(output_rows),
        "load_seconds": load_seconds,
        "import_inventories_seconds": import_seconds,
        "build_cf_seconds": build_cf_seconds,
        "year_loop_and_scores_csv_seconds": year_loop_and_scores_csv_seconds,
        "write_scores_csv_seconds_sum": float(write_year_times.sum()),
        "total_seconds": total_seconds,
        "adjoint_solve_seconds_sum": float(solve_times.sum()),
        "adjoint_solve_seconds_mean": float(solve_times.mean()),
        "adjoint_solve_seconds_max": float(solve_times.max()),
        "year_total_seconds_mean": float(year_times.mean()),
        "year_total_seconds_max": float(year_times.max()),
        "scores_csv": str(paths["scores_csv"]),
        "matrix_npz": str(paths["matrix_npz"]),
        "timings_csv": str(paths["timings_csv"]),
        "metadata_csv": str(paths["metadata_csv"]),
    }
    with paths["summary_json"].open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2), flush=True)
    print(f"Wrote scores CSV: {paths['scores_csv']}", flush=True)
    print(f"Wrote score matrix: {paths['matrix_npz']}", flush=True)
    print(f"Wrote timings: {paths['timings_csv']}", flush=True)
    print(f"Wrote metadata: {paths['metadata_csv']}", flush=True)
    print(f"Wrote summary: {paths['summary_json']}", flush=True)
    return 0


def run(args: argparse.Namespace) -> int:
    if args.all_years:
        return _run_all_years(args)

    datapackage = Path(args.datapackage).expanduser().resolve()
    if not datapackage.exists():
        raise FileNotFoundError(datapackage)

    lcia_json = None if args.lcia_json is None else Path(args.lcia_json).expanduser()
    if lcia_json is not None:
        if not lcia_json.exists():
            raise FileNotFoundError(lcia_json)
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json.resolve())

    available = get_lcia_method_names(ei_version=str(args.ei_version))
    if args.method not in available:
        raise ValueError(
            f"LCIA method not found for ecoinvent {args.ei_version}: {args.method}"
        )

    output_csv = Path(args.output_csv).expanduser().resolve()
    output_json = Path(args.output_json).expanduser().resolve()

    print(f"Datapackage: {datapackage}", flush=True)
    print(f"Year: {int(args.year)}", flush=True)
    print(f"Method: {args.method}", flush=True)

    total_start = time.perf_counter()
    trails, load_seconds, import_seconds = _load_trails(args)
    intensities, direct_scores, A_csc, scenario_year, timings = _solve_adjoint(
        trails=trails,
        year=int(args.year),
        method=args.method,
        ei_version=str(args.ei_version),
    )
    write_seconds = _write_scores(
        trails=trails,
        year=int(args.year),
        method=args.method,
        output_csv=output_csv,
        intensities=intensities,
        direct_scores=direct_scores,
        A_csc=A_csc,
    )
    total_seconds = time.perf_counter() - total_start

    finite = np.isfinite(intensities)
    summary = {
        "datapackage": str(datapackage),
        "inventories": [
            str(Path(path).expanduser().resolve()) for path in args.inventories
        ],
        "year_requested": int(args.year),
        "year_used": int(scenario_year),
        "method": args.method,
        "ei_version": str(args.ei_version),
        "n_activities": int(A_csc.shape[1]),
        "n_products": int(A_csc.shape[0]),
        "technosphere_nnz": int(A_csc.nnz),
        "biosphere_nnz_for_year": int(
            getattr(
                trails.B[trails.scenario_index[str(scenario_year)], :, :],
                "nnz",
                0,
            )
        ),
        "finite_intensities": int(finite.sum()),
        "min_intensity": float(np.nanmin(intensities)),
        "max_intensity": float(np.nanmax(intensities)),
        "load_seconds": load_seconds,
        "import_inventories_seconds": import_seconds,
        "write_csv_seconds": write_seconds,
        "total_seconds": total_seconds,
        **timings,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2), flush=True)
    print(f"Wrote scores: {output_csv}", flush=True)
    print(f"Wrote summary: {output_json}", flush=True)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute static LCIA score intensities for every activity in one "
            "Trails scenario year and one LCIA indicator."
        )
    )
    parser.add_argument("--datapackage", type=Path, default=DEFAULT_DATAPACKAGE)
    parser.add_argument("--year", type=int, default=2025)
    parser.add_argument("--method", default=DEFAULT_METHOD)
    parser.add_argument("--ei-version", default="3.12")
    parser.add_argument("--lcia-json", type=Path, default=DEFAULT_LCIA_JSON)
    parser.add_argument(
        "--inventories",
        type=Path,
        nargs="*",
        default=[],
        help="Optional Excel inventories to import after annual interpolation.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT_DIR / "activity_static_scores_2025_pm.csv",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=DEFAULT_OUTPUT_DIR / "activity_static_scores_2025_pm_summary.json",
    )
    parser.add_argument("--output-npz", type=Path, default=None)
    parser.add_argument("--timings-csv", type=Path, default=None)
    parser.add_argument("--metadata-csv", type=Path, default=None)
    parser.add_argument(
        "--all-years",
        action="store_true",
        help="Compute activity score intensities for every interpolated time slice.",
    )
    parser.add_argument("--no-cache-interpolation", action="store_true")
    parser.add_argument("--interpolation-start-year-offset", type=int, default=-20)
    parser.add_argument("--interpolation-end-year-offset", type=int, default=20)
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
