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

# Absolute plotting thresholds. They do not affect routing; they only keep
# the Sankey diagrams legible by showing major routed flows.
SANKEY_MIN_EDGE_BY_KEY = {
    "bev": 50.0,
    "polyol": 50_000_000.0,
    "marine": 100_000_000.0,
    "daccs": 50_000_000.0,
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


def _graph_summary(trails: Any, *, min_edge_amount: float) -> dict[str, Any]:
    graph = trails.graph
    depths: dict[int, int] = {}
    years: list[int] = []
    frontier_nodes = 0
    direct_bio_nodes = 0
    plotted_edges = 0
    plotted_nodes: set[Any] = set()
    depth1_edges: dict[str, float] = {}

    for node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        depths[depth] = depths.get(depth, 0) + 1
        years.append(int(data.get("year", 0)))
        if float(data.get("frontier_amount") or 0.0):
            frontier_nodes += 1
        if float(data.get("direct_bio_amount") or 0.0):
            direct_bio_nodes += 1

    for u, v, data in graph.edges(data=True):
        amount = abs(float(data.get("amount", 0.0)))
        if amount >= float(min_edge_amount):
            plotted_edges += 1
            plotted_nodes.add(u)
            plotted_nodes.add(v)
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
        "plotted_nodes": int(len(plotted_nodes)),
        "plotted_edges": int(plotted_edges),
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

    from trails.plotting import plot_temporal_sankey_graphlike

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
        min_edge = float(SANKEY_MIN_EDGE_BY_KEY[key])
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
            )
            routing_seconds = time.perf_counter() - t0
            summary = _graph_summary(trails, min_edge_amount=min_edge)
            print(
                f"  routing done in {routing_seconds:.1f}s "
                f"(nodes={summary['graph_nodes']:,}, edges={summary['graph_edges']:,}; "
                f"plotted_edges={summary['plotted_edges']:,})",
                flush=True,
            )

            html_path = case_dir / f"{key}_depth_{depth}_sankey.html"
            title = (
                f"{label}<br>Temporal routed Sankey, depth {depth} "
                f"(shown edges >= {min_edge:g})"
            )
            t1 = time.perf_counter()
            plot_temporal_sankey_graphlike(
                trails,
                min_edge_amount=min_edge,
                edge_weight="amount",
                title=title,
                amount_label="Routed amount",
                fig_width=1500,
                fig_height=1100 if depth == 1 else 1400,
                node_thickness=12,
                node_pad=6,
                font_size=11,
                max_label_chars=42,
                layout_by_year_depth=True,
                orientation="year_x_depth_y",
                branch_dropdown=True,
                depth_dropdown=True,
                default_depth_level=min(int(depth), 2),
                year_slider=True,
                filename=str(html_path),
            )
            plot_seconds = time.perf_counter() - t1
            print(f"  wrote {html_path} in {plot_seconds:.1f}s", flush=True)
            rows.append(
                {
                    "case_key": key,
                    "activity_index": activity_index,
                    "activity": label,
                    "amount": amount,
                    "depth": int(depth),
                    "routing_min_amount": ROUTING_MIN_AMOUNT,
                    "sankey_min_edge_amount": min_edge,
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
