from __future__ import annotations

import importlib.util
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DATAPACKAGE = Path("/Users/romain/GitHub/premise/dev/trails_2026-05-31.zip")
LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")
OUTPUT_DIR = REPO_ROOT / "dev" / "notebook_runs" / "temporal_lci_depth_sweep_runner"
DELTA_CSV = OUTPUT_DIR / "depth_delta_explanations.csv"
PATH_CSV = OUTPUT_DIR / "depth_delta_routed_paths.csv"

REFERENCE_YEAR = 2025
DEPTH = 5
MIN_AMOUNT = 1e-3
INTERPOLATION_START_YEAR_OFFSET = -20
INTERPOLATION_END_YEAR_OFFSET = 20
TOP_DELTA_ROWS_PER_CASE = 8
TOP_FRONTIER_PATHS_PER_ROOT = 8

CASE_DEMAND_AMOUNTS_BY_KEY = {
    "bev": 150_000.0,
    "polyol": 50_000_000_000.0,
    "marine": 180_000_000_000.0,
    "daccs": 20_000_000_000.0,
}


def _load_runner() -> Any:
    spec = importlib.util.spec_from_file_location(
        "plot_terminal_lci_td_comparison", RUNNER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load runner script: {RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


def _activity_label(trails: Any, idx: int) -> str:
    for mapping in trails.activity_indices.values():
        meta = mapping.get(int(idx))
        if not meta:
            continue
        return (
            f"{meta.get('name', '')} | "
            f"{meta.get('reference product', '')} | "
            f"{meta.get('location', '')}"
        )
    return f"Activity {int(idx)}"


def _short_label(trails: Any, idx: int) -> str:
    return _activity_label(trails, idx).split("|")[0].strip()


def _case_key_for_activity(activity: Any, case_defs: dict[str, Any]) -> str | None:
    for key, candidate in case_defs.items():
        if candidate == activity:
            return str(key)
    return None


def _demand_amount(activity: Any, case_defs: dict[str, Any]) -> float:
    key = _case_key_for_activity(activity, case_defs)
    if key is None:
        return 1.0
    return float(CASE_DEMAND_AMOUNTS_BY_KEY[key])


def _find_one_path_from_root(
    graph: Any, root_idx: int, target: Any
) -> list[Any] | None:
    target_depth = int(graph.nodes[target].get("depth", -1))
    stack: list[tuple[Any, list[Any]]] = [(target, [target])]
    visited: set[Any] = set()
    while stack:
        node, suffix = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        data = graph.nodes[node]
        depth = int(data.get("depth", -1))
        act_idx = int(data.get("act_idx", -1))
        if depth == 1 and act_idx == int(root_idx):
            return list(reversed(suffix))
        if depth <= 1:
            continue
        for pred in graph.predecessors(node):
            pred_depth = int(graph.nodes[pred].get("depth", -1))
            if pred_depth < depth and pred_depth >= 1:
                stack.append((pred, suffix + [pred]))
    if target_depth == 1 and int(graph.nodes[target].get("act_idx", -1)) == int(
        root_idx
    ):
        return [target]
    return None


def _path_signature(
    graph: Any, path: list[Any]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return (
        tuple(int(graph.nodes[node].get("act_idx")) for node in path),
        tuple(int(graph.nodes[node].get("year")) for node in path),
    )


def _format_path(trails: Any, graph: Any, path: list[Any]) -> tuple[str, str]:
    years = [int(graph.nodes[node].get("year")) for node in path]
    labels = [
        _short_label(trails, int(graph.nodes[node].get("act_idx"))) for node in path
    ]
    return " -> ".join(labels), " -> ".join(str(year) for year in years)


def _exact_emitter_paths(
    trails: Any,
    root_idx: int,
    emitting_idx: int,
    max_paths: int = 6,
) -> list[dict[str, Any]]:
    graph = trails.graph
    targets = [
        node
        for node, data in graph.nodes(data=True)
        if int(data.get("act_idx", -1)) == int(emitting_idx)
        and int(data.get("depth", -1)) >= 1
    ]
    rows: list[dict[str, Any]] = []
    seen: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for target in sorted(
        targets,
        key=lambda n: abs(float(graph.nodes[n].get("amount") or 0.0)),
        reverse=True,
    ):
        path = _find_one_path_from_root(graph, root_idx, target)
        if not path:
            continue
        signature = _path_signature(graph, path)
        if signature in seen:
            continue
        seen.add(signature)
        label_path, year_path = _format_path(trails, graph, path)
        rows.append(
            {
                "path_type": "exact_emitting_activity_in_routing_graph",
                "path_depth": len(path) - 1,
                "path_years": year_path,
                "path_activities": label_path,
                "target_amount": float(graph.nodes[target].get("amount") or 0.0),
                "target_frontier_amount": float(
                    graph.nodes[target].get("frontier_amount") or 0.0
                ),
                "target_direct_bio_amount": float(
                    graph.nodes[target].get("direct_bio_amount") or 0.0
                ),
            }
        )
        if len(rows) >= max_paths:
            break
    return rows


def _dominant_frontier_paths(
    trails: Any,
    root_idx: int,
    max_paths: int = TOP_FRONTIER_PATHS_PER_ROOT,
) -> list[dict[str, Any]]:
    graph = trails.graph
    candidates: list[tuple[float, Any, float]] = []
    for node, data in graph.nodes(data=True):
        roots = data.get("frontier_roots") or {}
        amount = float(
            roots.get(int(root_idx), 0.0) or roots.get(str(root_idx), 0.0) or 0.0
        )
        if amount:
            candidates.append((abs(amount), node, amount))

    rows: list[dict[str, Any]] = []
    seen: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for _, node, root_amount in sorted(candidates, reverse=True):
        path = _find_one_path_from_root(graph, root_idx, node)
        if not path:
            continue
        signature = _path_signature(graph, path)
        if signature in seen:
            continue
        seen.add(signature)
        label_path, year_path = _format_path(trails, graph, path)
        rows.append(
            {
                "path_type": "dominant_frontier_path_for_root",
                "path_depth": len(path) - 1,
                "path_years": year_path,
                "path_activities": label_path,
                "target_amount": float(root_amount),
                "target_frontier_amount": float(
                    graph.nodes[node].get("frontier_amount") or 0.0
                ),
                "target_direct_bio_amount": float(
                    graph.nodes[node].get("direct_bio_amount") or 0.0
                ),
            }
        )
        if len(rows) >= max_paths:
            break
    return rows


def main() -> None:
    if LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(LCIA_JSON)

    runner = _load_runner()
    case_defs = dict(runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS)
    case_defs["polyol"] = runner.helpers.ActivityDef(
        "polyol precursor production from captured CO2",
        "polyol precursor",
        "RER",
    )
    activities = [case_defs[key] for key in ("bev", "polyol", "marine", "daccs")]

    delta = pd.read_csv(DELTA_CSV)
    rows: list[dict[str, Any]] = []

    print("Loading Trails context", flush=True)
    trails = runner._load_trails(
        datapackage=DATAPACKAGE.expanduser().resolve(),
        interpolation_cache_dir=None,
        inventory_paths=[
            path.expanduser().resolve()
            for path in runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS
        ],
        import_before_interpolation=False,
        remove_base_temporal_distributions=False,
        no_cache_interpolation=False,
        interpolation_start_year_offset=INTERPOLATION_START_YEAR_OFFSET,
        interpolation_end_year_offset=INTERPOLATION_END_YEAR_OFFSET,
    )
    activity_maps = runner.helpers._match_activity_indices(trails, activities)

    for activity in activities:
        case_key = _case_key_for_activity(activity, case_defs)
        activity_index = int(activity_maps[activity])
        amount = _demand_amount(activity, case_defs)
        print(
            f"\nRouting {case_key}: idx={activity_index}, amount={amount:g}",
            flush=True,
        )
        start = time.perf_counter()
        trails.temporal_routing(
            start_year=REFERENCE_YEAR,
            start_act_idx=activity_index,
            amount=amount,
            max_depth=DEPTH,
            min_amount=MIN_AMOUNT,
            show_progress=False,
            attribute_to_roots=True,
        )
        graph = trails.graph
        print(
            f"  done in {time.perf_counter() - start:.1f}s "
            f"(nodes={graph.number_of_nodes():,}, edges={graph.number_of_edges():,})",
            flush=True,
        )

        case_delta = (
            delta.loc[delta["case_key"] == case_key]
            .sort_values("share_abs_delta_pre_reference", ascending=False)
            .head(TOP_DELTA_ROWS_PER_CASE)
        )
        roots_seen: Counter[int] = Counter()
        for _, item in case_delta.iterrows():
            root_idx = int(item["root_activity_index"])
            emitting_idx = int(item["emitting_activity_index"])
            exact_paths = _exact_emitter_paths(trails, root_idx, emitting_idx)
            if not exact_paths:
                exact_paths = [
                    {
                        "path_type": "emitting_activity_not_in_routing_graph",
                        "path_depth": "",
                        "path_years": "",
                        "path_activities": "",
                        "target_amount": "",
                        "target_frontier_amount": "",
                        "target_direct_bio_amount": "",
                    }
                ]
            for path_row in exact_paths:
                rows.append(
                    {
                        "case_key": case_key,
                        "activity_index": activity_index,
                        "root_activity_index": root_idx,
                        "root_activity": item["root_activity"],
                        "emitting_activity_index": emitting_idx,
                        "emitting_activity": item["emitting_activity"],
                        "delta_pre_reference_score": item["delta_pre_reference_score"],
                        "share_abs_delta_pre_reference": item[
                            "share_abs_delta_pre_reference"
                        ],
                        "first_pre_reference_delta_year": item[
                            "first_pre_reference_delta_year"
                        ],
                        "peak_pre_reference_delta_year": item[
                            "peak_pre_reference_delta_year"
                        ],
                        **path_row,
                    }
                )
            roots_seen[root_idx] += 1

        for root_idx, _count in roots_seen.most_common():
            for path_row in _dominant_frontier_paths(trails, root_idx):
                rows.append(
                    {
                        "case_key": case_key,
                        "activity_index": activity_index,
                        "root_activity_index": root_idx,
                        "root_activity": _activity_label(trails, root_idx),
                        "emitting_activity_index": "",
                        "emitting_activity": "",
                        "delta_pre_reference_score": "",
                        "share_abs_delta_pre_reference": "",
                        "first_pre_reference_delta_year": "",
                        "peak_pre_reference_delta_year": "",
                        **path_row,
                    }
                )

    out = pd.DataFrame(rows)
    out.to_csv(PATH_CSV, index=False)
    print(f"\nWrote {PATH_CSV}", flush=True)


if __name__ == "__main__":
    main()
