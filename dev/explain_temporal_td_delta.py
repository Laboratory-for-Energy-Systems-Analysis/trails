from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DEFAULT_METHOD = (
    "EF v3.1 - human toxicity: carcinogenic - comparative toxic unit for human "
    "(CTUh)"
)
DEFAULT_OUTPUT_DIR = REPO_ROOT / "dev" / "temporal_td_delta_diagnostics"
DEFAULT_LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")


@dataclass(frozen=True)
class ScoreEntry:
    root: int
    activity: int
    year: int
    value: float


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "plot_terminal_lci_td_comparison",
        RUNNER_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load runner script: {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


runner = _load_runner()
helpers = runner.helpers


def _activity_label(trails: Any, activity_index: int) -> str:
    for mapping in trails.activity_indices.values():
        meta = mapping.get(int(activity_index))
        if not isinstance(meta, dict):
            continue
        parts = [
            str(meta.get("name") or f"Activity {activity_index}").strip(),
            str(meta.get("reference product") or "").strip(),
            str(meta.get("location") or "").strip(),
        ]
        return " | ".join(part for part in parts if part)
    return f"Activity {int(activity_index)}"


def _select_method_scores(scores: Any, method: str) -> Any:
    if "method" not in scores.dims:
        return scores
    methods = [str(value) for value in scores.coords["method"].values.tolist()]
    try:
        method_index = methods.index(str(method))
    except ValueError as exc:
        raise ValueError(
            f"Method not found in score tensor: {method!r}. Available: {methods}"
        ) from exc
    return scores.isel(method=method_index, drop=True)


def _iter_score_entries(trails: Any, method: str) -> list[ScoreEntry]:
    scores = getattr(trails, "scores", None)
    if scores is None:
        raise RuntimeError(
            "trails.scores is None; run LCA with compute_score=True first."
        )

    scores = _select_method_scores(scores, method)
    if "root activity" not in scores.dims:
        raise RuntimeError(
            "Score tensor has no 'root activity' dimension. Rerun both passes with "
            "root attribution enabled; for the existing runner, set "
            "--foreground-attribute-to-roots."
        )

    required = ("activity", "year", "root activity")
    missing = [dim for dim in required if dim not in scores.dims]
    if missing:
        raise RuntimeError(f"Score tensor is missing dimension(s): {missing}")
    scores = scores.transpose(*required)

    activities = np.asarray(scores.coords["activity"].values, dtype=np.int64)
    years = np.asarray(scores.coords["year"].values, dtype=np.int64)
    roots = np.asarray(scores.coords["root activity"].values, dtype=np.int64)
    data = scores.data

    entries: list[ScoreEntry] = []
    if hasattr(data, "coords") and hasattr(data, "data"):
        coords = np.asarray(data.coords, dtype=np.int64)
        values = np.asarray(data.data, dtype=float)
        for activity_i, year_i, root_i, value in zip(
            coords[0], coords[1], coords[2], values
        ):
            if value == 0.0:
                continue
            entries.append(
                ScoreEntry(
                    root=int(roots[int(root_i)]),
                    activity=int(activities[int(activity_i)]),
                    year=int(years[int(year_i)]),
                    value=float(value),
                )
            )
        return entries

    dense = np.asarray(scores.values, dtype=float)
    activity_i, year_i, root_i = np.nonzero(dense)
    for ai, yi, ri in zip(activity_i, year_i, root_i):
        value = float(dense[int(ai), int(yi), int(ri)])
        if value == 0.0:
            continue
        entries.append(
            ScoreEntry(
                root=int(roots[int(ri)]),
                activity=int(activities[int(ai)]),
                year=int(years[int(yi)]),
                value=value,
            )
        )
    return entries


def _aggregate(
    entries: list[ScoreEntry],
    fields: tuple[str, ...],
) -> dict[tuple[int, ...], float]:
    out: dict[tuple[int, ...], float] = defaultdict(float)
    for entry in entries:
        key = tuple(int(getattr(entry, field)) for field in fields)
        out[key] += float(entry.value)
    return dict(out)


def _delta_rows(
    all_entries: list[ScoreEntry],
    foreground_entries: list[ScoreEntry],
    fields: tuple[str, ...],
) -> list[dict[str, Any]]:
    all_values = _aggregate(all_entries, fields)
    foreground_values = _aggregate(foreground_entries, fields)
    rows: list[dict[str, Any]] = []
    for key in sorted(set(all_values) | set(foreground_values)):
        all_score = float(all_values.get(key, 0.0))
        foreground_score = float(foreground_values.get(key, 0.0))
        delta = all_score - foreground_score
        row = {field: key[pos] for pos, field in enumerate(fields)}
        row.update(
            {
                "score_all_tds": all_score,
                "score_foreground_tds_only": foreground_score,
                "delta_all_minus_foreground": delta,
                "abs_delta": abs(delta),
            }
        )
        rows.append(row)
    rows.sort(key=lambda row: float(row["abs_delta"]), reverse=True)
    return rows


def _add_labels(rows: list[dict[str, Any]], trails: Any) -> None:
    for row in rows:
        if "root" in row:
            row["root_label"] = _activity_label(trails, int(row["root"]))
        if "activity" in row:
            row["activity_label"] = _activity_label(trails, int(row["activity"]))


def _write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    *,
    top_n: int | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_to_write = rows if top_n is None else rows[: int(top_n)]
    fieldnames: list[str] = []
    for row in rows_to_write:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_to_write)


def _write_summary(
    path: Path,
    *,
    method: str,
    all_total: float,
    foreground_total: float,
) -> None:
    delta = float(all_total) - float(foreground_total)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"method: {method}",
                f"foreground + background TDs: {all_total:.12g}",
                f"foreground TDs only: {foreground_total:.12g}",
                f"delta: {delta:.12g}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _selected_activity(case_key: str) -> Any:
    try:
        return runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS[str(case_key)]
    except KeyError as exc:
        choices = ", ".join(sorted(runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS))
        raise ValueError(
            f"Unknown case {case_key!r}; choose one of: {choices}"
        ) from exc


def run(args: argparse.Namespace) -> int:
    if args.lcia_json is not None:
        lcia_json = Path(args.lcia_json).expanduser().resolve()
        if not lcia_json.exists():
            raise FileNotFoundError(f"LCIA JSON not found: {lcia_json}")
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json)

    datapackage = Path(args.datapackage).expanduser().resolve()
    inventories = [
        Path(path).expanduser().resolve()
        for path in (args.inventories or runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS)
    ]
    runner._validate_paths(datapackage, inventories)

    activity = _selected_activity(str(args.case))
    method = str(args.method)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Datapackage: {datapackage}", flush=True)
    print("Excel inventories:", flush=True)
    for path in inventories:
        print(f"  {path}", flush=True)
    print(
        f"Case: {activity.name} | {activity.reference_product} | {activity.location}",
        flush=True,
    )
    print(f"Method: {method}", flush=True)

    contexts: dict[str, Any] = {}
    activity_maps: dict[str, dict[Any, int]] = {}
    results: dict[str, Any] = {}
    try:
        for config in runner.RUN_CONFIGS:
            print(f"\nLoading run context: {config.label}", flush=True)
            contexts[config.key] = runner._load_trails(
                datapackage=datapackage,
                interpolation_cache_dir=(
                    None
                    if args.interpolation_cache_dir is None
                    else Path(args.interpolation_cache_dir).expanduser().resolve()
                ),
                inventory_paths=inventories,
                import_before_interpolation=bool(args.import_before_interpolation),
                remove_base_temporal_distributions=(
                    config.remove_base_temporal_distributions
                ),
                no_cache_interpolation=bool(args.no_cache_interpolation),
                interpolation_start_year_offset=int(
                    args.interpolation_start_year_offset
                ),
                interpolation_end_year_offset=int(args.interpolation_end_year_offset),
            )
            activity_maps[config.key] = helpers._match_activity_indices(
                contexts[config.key],
                [activity],
            )
            if activity not in activity_maps[config.key]:
                raise ValueError(
                    "Could not match case-study activity after import in "
                    f"{config.label}: {activity}"
                )

        routing_min_amount = float(args.routing_min_amount)
        if bool(args.scale_routing_min_amount_by_activity):
            routing_min_amount *= runner.DEFAULT_CASE_STUDY_ROUTING_AMOUNT_SCALES.get(
                activity,
                1.0,
            )
        amount = float(args.amount)
        if bool(args.scale_demand_amount_by_activity):
            amount *= runner.DEFAULT_CASE_STUDY_ROUTING_AMOUNT_SCALES.get(activity, 1.0)

        for config in runner.RUN_CONFIGS:
            activity_index = int(activity_maps[config.key][activity])
            results[config.key] = runner._run_activity(
                trails=contexts[config.key],
                run_config=config,
                activity_index=activity_index,
                methods=[method],
                depth=int(args.depth),
                reference_year=int(args.reference_year),
                amount=amount,
                show_progress=bool(args.show_progress),
                solver_mode=str(args.solver_mode),
                iterative_rtol=float(args.iterative_rtol),
                iterative_maxiter=args.iterative_maxiter,
                iterative_restart=args.iterative_restart,
                ei_version=str(args.ei_version),
                fallback_solver_mode=(
                    None
                    if str(args.fallback_solver_mode) == "none"
                    else str(args.fallback_solver_mode)
                ),
                routing_min_amount=routing_min_amount,
                attribute_to_roots=True,
                write_graph_html=False,
                graph_output_dir=output_dir,
                graph_run="all_td",
                graph_min_edge_amount=float(args.graph_min_edge_amount),
            )

        all_entries = _iter_score_entries(contexts["all_td"], method)
        foreground_entries = _iter_score_entries(contexts["foreground_td_only"], method)

        all_total = float(results["all_td"].total_by_method[method])
        foreground_total = float(results["foreground_td_only"].total_by_method[method])
        _write_summary(
            output_dir / "summary.txt",
            method=method,
            all_total=all_total,
            foreground_total=foreground_total,
        )

        for name, fields in {
            "delta_by_year.csv": ("year",),
            "delta_by_root.csv": ("root",),
            "delta_by_emitting_activity.csv": ("activity",),
            "delta_by_root_and_emitting_activity.csv": ("root", "activity"),
            "delta_by_root_activity_year.csv": ("root", "activity", "year"),
        }.items():
            rows = _delta_rows(all_entries, foreground_entries, fields)
            _add_labels(rows, contexts["all_td"])
            _write_csv(output_dir / name, rows, top_n=args.top_n)

        print(f"\nWrote diagnostic tables under: {output_dir}", flush=True)
        print(f"Total delta: {all_total - foreground_total:.12g}", flush=True)
        return 0
    finally:
        contexts.clear()
        results.clear()
        gc.collect()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Explain the difference between temporal LCA cumulative scores with "
            "foreground + background TDs and foreground TDs only."
        )
    )
    parser.add_argument(
        "--datapackage",
        type=Path,
        default=REPO_ROOT / "dev" / "trails_2026-05-18.zip",
    )
    parser.add_argument("--interpolation-cache-dir", type=Path, default=None)
    parser.add_argument(
        "--import-before-interpolation",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    parser.add_argument(
        "--inventories",
        nargs="+",
        type=Path,
        default=[],
        help=(
            "Excel inventories to import. Defaults to the four OneDrive case-study "
            "inventories."
        ),
    )
    parser.add_argument(
        "--case",
        choices=tuple(runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS),
        default="bev",
    )
    parser.add_argument("--method", default=DEFAULT_METHOD)
    parser.add_argument("--depth", type=int, default=5)
    parser.add_argument("--reference-year", type=int, default=2026)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--ei-version", type=str, default="3.12")
    parser.add_argument("--lcia-json", type=Path, default=DEFAULT_LCIA_JSON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--solver-mode", type=str, default="direct")
    parser.add_argument("--fallback-solver-mode", type=str, default="direct")
    parser.add_argument("--iterative-rtol", type=float, default=1e-8)
    parser.add_argument("--iterative-maxiter", type=int, default=None)
    parser.add_argument("--iterative-restart", type=int, default=None)
    parser.add_argument("--routing-min-amount", type=float, default=1e-12)
    parser.add_argument(
        "--scale-routing-min-amount-by-activity",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    parser.add_argument(
        "--scale-demand-amount-by-activity",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    parser.add_argument("--graph-min-edge-amount", type=float, default=1e-9)
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--no-cache-interpolation", action="store_true")
    parser.add_argument("--interpolation-start-year-offset", type=int, default=-20)
    parser.add_argument("--interpolation-end-year-offset", type=int, default=20)
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Optional number of largest absolute-delta rows to keep in each table.",
    )
    return parser


if __name__ == "__main__":
    raise SystemExit(run(build_parser().parse_args()))
