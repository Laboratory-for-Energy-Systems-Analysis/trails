from __future__ import annotations


import importlib.util
import os
import sys
import time
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
OUTPUT_DIR = (
    REPO_ROOT / "dev" / "notebook_runs" / "temporal_lci_depth_sweep_runner" / "sankey"
)

REFERENCE_YEAR = 2025
DEPTHS = [1, 5]
ROUTING_MIN_AMOUNT = 1e-3
INTERPOLATION_START_YEAR_OFFSET = -20
INTERPOLATION_END_YEAR_OFFSET = 20

CASE_KEYS = ["bev", "polyol", "marine", "daccs"]
CASE_DEMAND_AMOUNTS_BY_KEY = {
    "bev": 150_000.0,
    "polyol": 50_000_000_000.0,
    "marine": 180_000_000_000.0,
    "daccs": 20_000_000_000.0,
}

CASE_METHODS_BY_KEY = {
    "polyol": (
        "EF v3.1 - material resources: metals/minerals - abiotic depletion "
        "potential (ADP): elements (ultimate reserves)"
    ),
    "bev": (
        "EF v3.1 - ecotoxicity: freshwater - comparative toxic unit for "
        "ecosystems (CTUe)"
    ),
    "daccs": "EF v3.1 - ozone depletion - ozone depletion potential (ODP)",
    "marine": (
        "EF v3.1 - human toxicity: carcinogenic - comparative toxic unit "
        "for human (CTUh)"
    ),
}
ADAPTIVE_RELATIVE_SCORE_CUTOFF = 1e-4
BRANCH_VISUAL_CUTOFF = 0.001
DISPLAY_SCORE_COVERAGE = 1.0
MAX_SANKEY_LINKS = 0


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


def _slug(value: str) -> str:
    return (
        str(value)
        .lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("|", "_")
        .replace(",", "")
        .replace(":", "")
    )


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


def _graph_summary(trails: Any) -> dict[str, Any]:
    graph = trails.graph
    depths: dict[int, int] = {}
    years: list[int] = []
    frontier_nodes = 0
    direct_bio_nodes = 0
    score_potential_nodes = 0
    score_potential_edges = 0
    total_score_potential = 0.0
    depth1_edges: dict[str, float] = {}

    for node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        depths[depth] = depths.get(depth, 0) + 1
        years.append(int(data.get("year", 0)))
        if float(data.get("frontier_amount") or 0.0):
            frontier_nodes += 1
        if float(data.get("direct_bio_amount") or 0.0):
            direct_bio_nodes += 1
        score_potential = abs(float(data.get("score_potential") or 0.0))
        if score_potential:
            score_potential_nodes += 1
            total_score_potential += score_potential

    for u, v, data in graph.edges(data=True):
        amount = abs(float(data.get("amount", 0.0)))
        child_score = abs(float(graph.nodes[v].get("score_potential") or 0.0))
        if child_score:
            score_potential_edges += 1
        if int(graph.nodes[u].get("depth", -1)) == 0:
            label = _activity_label(trails, int(graph.nodes[v].get("act_idx")))
            depth1_edges[label] = depth1_edges.get(label, 0.0) + amount

    out: dict[str, Any] = {
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
        "year_min": min(years) if years else "",
        "year_max": max(years) if years else "",
        "score_potential_nodes": int(score_potential_nodes),
        "score_potential_edges": int(score_potential_edges),
        "total_score_potential": float(total_score_potential),
    }
    for depth, count in sorted(depths.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    for i, (label, amount) in enumerate(
        sorted(depth1_edges.items(), key=lambda item: item[1], reverse=True)[:12],
        start=1,
    ):
        out[f"depth1_branch_{i}"] = label
        out[f"depth1_branch_{i}_amount"] = float(amount)
    return out


def main() -> None:
    if LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(LCIA_JSON)

    from trails.plotting import plot_adaptive_sankey

    runner = _load_runner()
    case_defs = dict(runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS)
    case_defs["polyol"] = runner.helpers.ActivityDef(
        "polyol precursor production from captured CO2",
        "polyol precursor",
        "RER",
    )
    activities = [case_defs[key] for key in CASE_KEYS]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
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

    rows: list[dict[str, Any]] = []
    for key, activity in zip(CASE_KEYS, activities):
        activity_index = int(activity_maps[activity])
        amount = float(CASE_DEMAND_AMOUNTS_BY_KEY[key])
        label = _activity_label(trails, activity_index)
        method = str(CASE_METHODS_BY_KEY[key])
        case_dir = OUTPUT_DIR / f"{key}_{activity_index}"
        case_dir.mkdir(parents=True, exist_ok=True)

        for depth in DEPTHS:
            print(
                f"\nRouting {key}, depth={depth}, idx={activity_index}, "
                f"amount={amount:g}",
                flush=True,
            )
            t0 = time.perf_counter()
            trails.temporal_routing(
                start_year=REFERENCE_YEAR,
                start_act_idx=activity_index,
                amount=amount,
                max_depth=int(depth),
                min_amount=float(ROUTING_MIN_AMOUNT),
                show_progress=False,
                attribute_to_roots=True,
                adaptive_methods=[method],
                adaptive_relative_score_cutoff=float(ADAPTIVE_RELATIVE_SCORE_CUTOFF),
            )
            routing_seconds = time.perf_counter() - t0
            summary = _graph_summary(trails)
            print(
                f"  routing done in {routing_seconds:.1f}s "
                f"(nodes={summary['graph_nodes']:,}, edges={summary['graph_edges']:,}; "
                f"score-potential edges={summary['score_potential_edges']:,})",
                flush=True,
            )

            html_path = case_dir / f"{key}_depth_{depth}_sankey.html"
            title = (
                f"{label}<br>Adaptive routed Sankey, depth {depth} "
                f"(cutoff={ADAPTIVE_RELATIVE_SCORE_CUTOFF:g})"
            )
            t1 = time.perf_counter()
            plot_adaptive_sankey(
                trails,
                method=method,
                title=title,
                adaptive_relative_score_cutoff=float(ADAPTIVE_RELATIVE_SCORE_CUTOFF),
                branch_visual_cutoff=float(BRANCH_VISUAL_CUTOFF),
                display_score_coverage=float(DISPLAY_SCORE_COVERAGE),
                max_sankey_links=int(MAX_SANKEY_LINKS),
                width=1500,
                height=1100 if depth == 1 else 1400,
                output_path=html_path,
            )
            plot_seconds = time.perf_counter() - t1
            print(f"  wrote {html_path} in {plot_seconds:.1f}s", flush=True)
            rows.append(
                {
                    "case_key": key,
                    "activity_index": activity_index,
                    "activity": label,
                    "amount": amount,
                    "method": method,
                    "depth": int(depth),
                    "routing_min_amount": ROUTING_MIN_AMOUNT,
                    "adaptive_relative_score_cutoff": (ADAPTIVE_RELATIVE_SCORE_CUTOFF),
                    "branch_visual_cutoff": BRANCH_VISUAL_CUTOFF,
                    "routing_seconds": routing_seconds,
                    "plot_seconds": plot_seconds,
                    "html_path": str(html_path),
                    **summary,
                }
            )

    summary_path = OUTPUT_DIR / "sankey_summary.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    print(f"\nWrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
