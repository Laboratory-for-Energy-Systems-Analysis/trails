from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
from datapackage import Package

from trails import Trails

DEFAULT_METHODS = [
    (
        "IPCC 2021 (incl. biogenic CO2) - climate change: total "
        "(incl. biogenic CO2) - global warming potential (GWP100)"
    )
]
DEFAULT_REFERENCE_YEARS = (2025, 2026)
DEFAULT_ACTIVITY_INDEX = 31387
DEFAULT_PACKAGE_NAME = "trails_2026-03-18.zip"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resolve_package_path(raw_path: str | None) -> Path:
    if raw_path:
        candidate = Path(raw_path).expanduser()
        if candidate.exists():
            return candidate.resolve()
        raise FileNotFoundError(f"Package path does not exist: {candidate}")

    repo_root = _repo_root()
    candidates = [
        repo_root / "dev" / DEFAULT_PACKAGE_NAME,
        repo_root / DEFAULT_PACKAGE_NAME,
        Path.cwd() / DEFAULT_PACKAGE_NAME,
        repo_root / "dev" / DEFAULT_PACKAGE_NAME.removesuffix(".zip"),
        repo_root / DEFAULT_PACKAGE_NAME.removesuffix(".zip"),
        Path.cwd() / DEFAULT_PACKAGE_NAME.removesuffix(".zip"),
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not find the default package. Checked:\n"
        + "\n".join(f"  - {path}" for path in candidates)
    )


def _load_trails(package_path: Path, *, debug: bool) -> Trails:
    package = Package(str(package_path))
    return Trails(package, interpolate_annual=True, debug=debug)


def _to_1d_float_array(value: object) -> np.ndarray:
    raw = getattr(value, "data", value)
    if hasattr(raw, "todense"):
        raw = raw.todense()
    elif hasattr(value, "values"):
        raw = getattr(value, "values")
    arr = np.asarray(raw, dtype=float)
    if arr.ndim == 0:
        return np.array([float(arr)], dtype=float)
    return arr.reshape(-1).astype(float, copy=False)


def _extract_dynamic_score_vector(trails: Trails) -> np.ndarray:
    if trails.scores is None:
        raise RuntimeError("trails.scores is None after temporal LCA.")

    score_da = trails.scores
    if "method" in score_da.dims:
        reduce_dims = [dim for dim in score_da.dims if dim != "method"]
        return _to_1d_float_array(score_da.sum(dim=reduce_dims))

    return _to_1d_float_array(score_da.sum())


def _extract_characterized_inventory_vector(trails: Trails) -> np.ndarray:
    if trails.characterized_inventory is None:
        raise RuntimeError(
            "trails.characterized_inventory is None after temporal LCA."
        )

    char_da = trails.characterized_inventory
    if "method" in char_da.dims:
        reduce_dims = [dim for dim in char_da.dims if dim != "method"]
        return _to_1d_float_array(char_da.sum(dim=reduce_dims))

    return _to_1d_float_array(char_da.sum())


def _extract_static_score_vector(trails: Trails) -> np.ndarray:
    if trails.static_score is None:
        raise RuntimeError("trails.static_score is None after static LCA.")
    return _to_1d_float_array(trails.static_score)


def _activity_metadata(trails: Trails, act_idx: int) -> tuple[str | None, dict]:
    for label, mapping in trails.activity_indices.items():
        if act_idx in mapping:
            meta = mapping[act_idx]
            if isinstance(meta, dict):
                return str(label), meta
            break
    return None, {}


def _format_rel_delta(dynamic_score: float, static_score: float) -> str:
    if static_score == 0.0:
        return "n/a"
    return f"{100.0 * (dynamic_score - static_score) / static_score:+.3f}%"


def _print_header(
    *,
    package_path: Path,
    act_idx: int,
    metadata_label: str | None,
    metadata: dict,
    years: Iterable[int],
    solver_mode: str,
    fresh_trails_per_year: bool,
) -> None:
    print(f"Package: {package_path}")
    print(f"Activity index: {int(act_idx)}")
    if metadata:
        print(f"Scenario label for metadata lookup: {metadata_label}")
        print(f"Name: {metadata.get('name', '')}")
        print(f"Reference product: {metadata.get('reference product', '')}")
        print(f"Location: {metadata.get('location', '')}")
        print(f"Unit: {metadata.get('unit', '')}")
    print(f"Reference years: {', '.join(str(int(y)) for y in years)}")
    print(f"Solver mode: {solver_mode}")
    print(f"Fresh Trails per year: {bool(fresh_trails_per_year)}")


def _run_case(
    *,
    trails: Trails,
    year: int,
    act_idx: int,
    methods: list[str],
    amount: float,
    max_depth: int,
    show_progress: bool,
    solver_mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    trails.temporal_routing(
        start_year=int(year),
        start_act_idx=int(act_idx),
        amount=float(amount),
        max_depth=int(max_depth),
        show_progress=show_progress,
        attribute_to_roots=True,
    )

    trails.lca(
        methods=methods,
        show_progress=show_progress,
        compute_score=True,
        store_inventory=True,
        solver_mode=solver_mode,
    )

    dynamic_scores = _extract_dynamic_score_vector(trails)
    characterized_scores = _extract_characterized_inventory_vector(trails)

    trails.static_lca(
        year=int(year),
        act_idx=int(act_idx),
        methods=methods,
        amount=float(amount),
    )
    static_scores = _extract_static_score_vector(trails)
    return dynamic_scores, characterized_scores, static_scores


def _print_case_results(
    *,
    year: int,
    methods: list[str],
    dynamic_scores: np.ndarray,
    characterized_scores: np.ndarray,
    static_scores: np.ndarray,
    trails: Trails,
) -> None:
    graph = getattr(trails, "graph", None)
    if graph is not None:
        print(
            f"\nReference year {int(year)}"
            f" | routing nodes={graph.number_of_nodes()}"
            f" edges={graph.number_of_edges()}"
        )
    else:
        print(f"\nReference year {int(year)}")

    for method_idx, method_name in enumerate(methods):
        dyn = float(dynamic_scores[min(method_idx, dynamic_scores.size - 1)])
        char = float(
            characterized_scores[min(method_idx, characterized_scores.size - 1)]
        )
        sta = float(static_scores[min(method_idx, static_scores.size - 1)])

        print(f"Method: {method_name}")
        print(f"  Dynamic score:               {dyn:.12g}")
        print(f"  Characterized inventory sum: {char:.12g}")
        print(f"  Static score:                {sta:.12g}")
        print(f"  Dynamic - static:            {dyn - sta:+.12g}")
        print(f"  Relative delta:              {_format_rel_delta(dyn, sta)}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the core workflow from dev/tom.ipynb for activity index 31387 "
            "and compare temporal vs static scores for 2025 and 2026."
        )
    )
    parser.add_argument(
        "--package",
        default=None,
        help=(
            "Path to the TRAILS datapackage zip or directory. Defaults to "
            "dev/trails_2026-03-18.zip when present."
        ),
    )
    parser.add_argument(
        "--activity-index",
        type=int,
        default=DEFAULT_ACTIVITY_INDEX,
        help="Activity index to reproduce.",
    )
    parser.add_argument(
        "--reference-years",
        type=int,
        nargs="+",
        default=list(DEFAULT_REFERENCE_YEARS),
        help="Reference years to run sequentially.",
    )
    parser.add_argument(
        "--amount",
        type=float,
        default=1.0,
        help="Functional unit amount.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=4,
        help="Temporal routing max depth.",
    )
    parser.add_argument(
        "--method",
        action="append",
        default=None,
        help="LCIA method name. Repeat to pass multiple methods.",
    )
    parser.add_argument(
        "--solver-mode",
        choices=("iterative", "direct", "bw2calc"),
        default="iterative",
        help="Solver mode passed to trails.lca().",
    )
    parser.add_argument(
        "--fresh-trails-per-year",
        action="store_true",
        help=(
            "Instantiate a fresh Trails object for each reference year instead of "
            "reusing the notebook-style single instance."
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable Trails debug mode.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable routing and solve progress bars.",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    package_path = _resolve_package_path(args.package)
    methods = list(args.method) if args.method else list(DEFAULT_METHODS)
    reference_years = [int(year) for year in args.reference_years]
    show_progress = not bool(args.quiet)

    trails = None
    if not args.fresh_trails_per_year:
        trails = _load_trails(package_path, debug=bool(args.debug))

    header_trails = trails if trails is not None else _load_trails(
        package_path, debug=bool(args.debug)
    )
    metadata_label, metadata = _activity_metadata(
        header_trails, int(args.activity_index)
    )
    _print_header(
        package_path=package_path,
        act_idx=int(args.activity_index),
        metadata_label=metadata_label,
        metadata=metadata,
        years=reference_years,
        solver_mode=str(args.solver_mode),
        fresh_trails_per_year=bool(args.fresh_trails_per_year),
    )

    for year in reference_years:
        if args.fresh_trails_per_year:
            trails = _load_trails(package_path, debug=bool(args.debug))

        assert trails is not None
        dynamic_scores, characterized_scores, static_scores = _run_case(
            trails=trails,
            year=int(year),
            act_idx=int(args.activity_index),
            methods=methods,
            amount=float(args.amount),
            max_depth=int(args.max_depth),
            show_progress=show_progress,
            solver_mode=str(args.solver_mode),
        )
        _print_case_results(
            year=int(year),
            methods=methods,
            dynamic_scores=dynamic_scores,
            characterized_scores=characterized_scores,
            static_scores=static_scores,
            trails=trails,
        )


if __name__ == "__main__":
    main()
