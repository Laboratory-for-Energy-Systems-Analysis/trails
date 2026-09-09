from __future__ import annotations

from dataclasses import dataclass
import argparse
import time
from pathlib import Path
from typing import Any

import pandas as pd
from datapackage import Package

from trails import Trails
from trails.plotting import plot_adaptive_sankey


REPO_ROOT = Path(__file__).resolve().parents[2]
PUBLICATION_DIR = Path(__file__).resolve().parent
DATAPACKAGE = PUBLICATION_DIR / "trails_remind_SSP2-PkBudg1000.zip"
LCI_DIR = PUBLICATION_DIR / "LCIs"
INVENTORY_PATHS = [
    LCI_DIR / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    LCI_DIR / "lci-case-study-daccs_storage_risk.xlsx",
    LCI_DIR / "lci-case-study-marine_fuel_switch.xlsx",
    LCI_DIR / "lci-pass_cars.xlsx",
]
OUTPUT_DIR = (
    REPO_ROOT / "dev" / "notebook_runs" / "temporal_lci_selected_indicator_sankeys"
)

REFERENCE_YEAR = 2025
EI_VERSION = "3.12"
INTERPOLATION_START_YEAR_OFFSET = -20
INTERPOLATION_END_YEAR_OFFSET = 20
ROUTING_MIN_AMOUNT = 0.0
ADAPTIVE_RELATIVE_SCORE_CUTOFF = 1e-4
ADAPTIVE_MIN_DEPTH = 1
BRANCH_VISUAL_CUTOFF = 0.001
DISPLAY_SCORE_COVERAGE = 1.0
ADAPTIVE_MAX_SANKEY_LINKS = 5000
SANKEY_FIGURE_WIDTH = 1320
SANKEY_FIGURE_HEIGHT = 660
SANKEY_DEPTH_AXIS_MAX = 17
SANKEY_PNG_SCALE = 3


@dataclass(frozen=True)
class ActivityDef:
    name: str
    reference_product: str
    location: str


@dataclass(frozen=True)
class CaseDef:
    key: str
    activity: ActivityDef
    amount: float
    method: str


CASE_DEFS = [
    CaseDef(
        key="bev",
        activity=ActivityDef(
            "transport, passenger, car, battery electric",
            "transport, passenger, car",
            "RER",
        ),
        amount=150_000.0,
        method=(
            "EF v3.1 - ecotoxicity: freshwater, inorganics - comparative "
            "toxic unit for ecosystems (CTUe)"
        ),
    ),
    CaseDef(
        key="polyol",
        activity=ActivityDef(
            "polyol precursor production from captured CO2",
            "polyol precursor",
            "RER",
        ),
        amount=50_000_000_000.0,
        method=(
            "EF v3.1 - material resources: metals/minerals - abiotic "
            "depletion potential (ADP): elements (ultimate reserves)"
        ),
    ),
    CaseDef(
        key="marine",
        activity=ActivityDef(
            "marine freight service, temporal fuel transition",
            "transport service",
            "RER",
        ),
        amount=180_000_000_000.0,
        method="EF v3.1 - ozone depletion - ozone depletion potential (ODP)",
    ),
    CaseDef(
        key="daccs",
        activity=ActivityDef(
            "carbon dioxide, captured, with a solvent-based direct air capture "
            "system, 1MtCO2",
            "carbon dioxide, captured",
            "Europe",
        ),
        amount=20_000_000_000.0,
        method=(
            "EF v3.1 - human toxicity: carcinogenic, inorganics - "
            "comparative toxic unit for human (CTUh)"
        ),
    ),
]


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _metadata_by_idx(trails: Trails) -> dict[int, dict[str, Any]]:
    if not getattr(trails, "activity_indices", None):
        return {}
    first_label = next(iter(trails.activity_indices))
    return {
        int(key): value
        for key, value in trails.activity_indices[first_label].items()
        if isinstance(value, dict)
    }


def _match_activity_index(trails: Trails, target: ActivityDef) -> int:
    matches = [
        index
        for index, metadata in _metadata_by_idx(trails).items()
        if _clean(metadata.get("name")) == target.name
        and _clean(metadata.get("reference product")) == target.reference_product
        and _clean(metadata.get("location")) == target.location
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one match for {target}, found {matches}")
    return int(matches[0])


def _activity_label(trails: Trails, activity_index: int) -> str:
    for mapping in trails.activity_indices.values():
        metadata = mapping.get(int(activity_index))
        if not isinstance(metadata, dict):
            continue
        return " | ".join(
            part
            for part in (
                _clean(metadata.get("name")),
                _clean(metadata.get("reference product")),
                _clean(metadata.get("location")),
            )
            if part
        )
    return f"Activity {int(activity_index)}"


def _method_label(method: str) -> str:
    return method.removeprefix("EF v3.1 - ")


def _graph_summary(trails: Trails) -> dict[str, Any]:
    graph = trails.graph
    depths: dict[int, int] = {}
    years: list[int] = []
    frontier_nodes = 0
    score_potential_edges = 0
    total_score_potential = 0.0

    for _node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        depths[depth] = depths.get(depth, 0) + 1
        years.append(int(data.get("year", 0)))
        if float(data.get("frontier_amount") or 0.0):
            frontier_nodes += 1
        score_potential = abs(float(data.get("score_potential") or 0.0))
        total_score_potential += score_potential

    for _source, target, _data in graph.edges(data=True):
        child_score = abs(float(graph.nodes[target].get("score_potential") or 0.0))
        if child_score:
            score_potential_edges += 1

    out: dict[str, Any] = {
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "frontier_nodes": int(frontier_nodes),
        "year_min": min(years) if years else "",
        "year_max": max(years) if years else "",
        "score_potential_edges": int(score_potential_edges),
        "total_score_potential": float(total_score_potential),
    }
    for depth, count in sorted(depths.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    return out


def _load_trails() -> Trails:
    missing = [path for path in [DATAPACKAGE, *INVENTORY_PATHS] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing input files:\n" + "\n".join(map(str, missing)))

    print(f"Loading datapackage: {DATAPACKAGE}", flush=True)
    trails = Trails(
        Package(str(DATAPACKAGE)),
        interpolate_annual=True,
        cache_interpolation=True,
        interpolation_start_year_offset=INTERPOLATION_START_YEAR_OFFSET,
        interpolation_end_year_offset=INTERPOLATION_END_YEAR_OFFSET,
        ei_version=EI_VERSION,
    )
    print(
        "Importing foreground inventories: "
        + ", ".join(path.name for path in INVENTORY_PATHS),
        flush=True,
    )
    trails.import_excel_inventory([str(path) for path in INVENTORY_PATHS])
    return trails


def _route_and_plot(
    trails: Trails,
    *,
    case: CaseDef,
    activity_index: int,
    mode: str,
    max_depth: int | None,
    max_sankey_links: int,
    height_scale: float,
    show_depth_density: bool,
    year_axis_min: int | None = None,
    year_axis_max: int | None = None,
) -> dict[str, Any]:
    case_dir = OUTPUT_DIR / f"{case.key}_{activity_index}"
    case_dir.mkdir(parents=True, exist_ok=True)
    suffix = "depth_1" if mode == "depth_1" else "adaptive_cutoff_1e-4"
    html_path = case_dir / f"{case.key}_{suffix}_sankey.html"
    png_path = case_dir / f"{case.key}_{suffix}_sankey.png"

    print(
        f"\nRouting {case.key}, mode={mode}, idx={activity_index}, "
        f"amount={case.amount:g}",
        flush=True,
    )
    t0 = time.perf_counter()
    trails.temporal_routing(
        start_year=REFERENCE_YEAR,
        start_act_idx=int(activity_index),
        amount=float(case.amount),
        max_depth=max_depth,
        min_amount=float(ROUTING_MIN_AMOUNT),
        show_progress=False,
        attribute_to_roots=True,
        adaptive_methods=[case.method],
        adaptive_relative_score_cutoff=float(ADAPTIVE_RELATIVE_SCORE_CUTOFF),
        adaptive_min_depth=int(ADAPTIVE_MIN_DEPTH),
    )
    routing_seconds = time.perf_counter() - t0
    summary = _graph_summary(trails)
    print(
        f"  routing done in {routing_seconds:.1f}s "
        f"(nodes={summary['graph_nodes']:,}, edges={summary['graph_edges']:,}; "
        f"score-potential edges={summary['score_potential_edges']:,})",
        flush=True,
    )

    title = (
        f"{_activity_label(trails, activity_index)}<br>"
        f"{_method_label(case.method)}<br>"
        f"{'Fixed depth 1' if mode == 'depth_1' else 'Adaptive routing'}; "
        f"cutoff={ADAPTIVE_RELATIVE_SCORE_CUTOFF:.0e}"
    )
    t1 = time.perf_counter()
    figure_height = max(400, int(SANKEY_FIGURE_HEIGHT * float(height_scale)))
    plot_year_axis_min = int(
        year_axis_min if year_axis_min is not None else summary["year_min"]
    )
    plot_year_axis_max = int(
        year_axis_max if year_axis_max is not None else summary["year_max"]
    )
    fig = plot_adaptive_sankey(
        trails,
        method=case.method,
        title=title,
        adaptive_relative_score_cutoff=float(ADAPTIVE_RELATIVE_SCORE_CUTOFF),
        branch_visual_cutoff=float(BRANCH_VISUAL_CUTOFF),
        display_score_coverage=float(DISPLAY_SCORE_COVERAGE),
        max_sankey_links=int(max_sankey_links),
        width=SANKEY_FIGURE_WIDTH,
        height=figure_height,
        depth_axis_max=SANKEY_DEPTH_AXIS_MAX,
        year_axis_min=plot_year_axis_min,
        year_axis_max=plot_year_axis_max,
        show_depth_density=bool(show_depth_density),
        output_path=html_path,
        png_path=png_path,
        png_scale=SANKEY_PNG_SCALE,
    )
    plot_seconds = time.perf_counter() - t1
    print(f"  wrote {html_path} and {png_path} in {plot_seconds:.1f}s", flush=True)
    stats = dict(fig.layout.meta or {})
    return {
        "case_key": case.key,
        "activity_index": int(activity_index),
        "activity": _activity_label(trails, activity_index),
        "amount": float(case.amount),
        "method": case.method,
        "mode": mode,
        "max_depth": "" if max_depth is None else int(max_depth),
        "adaptive_relative_score_cutoff": float(ADAPTIVE_RELATIVE_SCORE_CUTOFF),
        "adaptive_min_depth": int(ADAPTIVE_MIN_DEPTH),
        "routing_min_amount": float(ROUTING_MIN_AMOUNT),
        "branch_visual_cutoff": float(BRANCH_VISUAL_CUTOFF),
        "max_sankey_links": int(max_sankey_links),
        "figure_width": int(SANKEY_FIGURE_WIDTH),
        "figure_height": int(figure_height),
        "depth_axis_max": int(SANKEY_DEPTH_AXIS_MAX),
        "requested_year_axis_min": "" if year_axis_min is None else int(year_axis_min),
        "requested_year_axis_max": "" if year_axis_max is None else int(year_axis_max),
        "png_scale": int(SANKEY_PNG_SCALE),
        "height_scale": float(height_scale),
        "show_depth_density": bool(show_depth_density),
        "routing_seconds": float(routing_seconds),
        "plot_seconds": float(plot_seconds),
        "html_path": str(html_path),
        "png_path": str(png_path),
        **summary,
        **stats,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate selected case-indicator Sankey diagrams."
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=[case.key for case in CASE_DEFS],
        help="Case key to run. May be supplied more than once. Defaults to all.",
    )
    parser.add_argument(
        "--height-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to the default Sankey figure heights.",
    )
    parser.add_argument(
        "--no-depth-density",
        action="store_true",
        help="Hide the bottom depth-density panel.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    selected_cases = set(args.case or [case.key for case in CASE_DEFS])
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    trails = _load_trails()

    rows: list[dict[str, Any]] = []
    for case in CASE_DEFS:
        if case.key not in selected_cases:
            continue
        activity_index = _match_activity_index(trails, case.activity)
        adaptive_row = _route_and_plot(
            trails,
            case=case,
            activity_index=activity_index,
            mode="adaptive",
            max_depth=None,
            max_sankey_links=ADAPTIVE_MAX_SANKEY_LINKS,
            height_scale=float(args.height_scale),
            show_depth_density=not bool(args.no_depth_density),
        )
        depth_1_row = _route_and_plot(
            trails,
            case=case,
            activity_index=activity_index,
            mode="depth_1",
            max_depth=1,
            max_sankey_links=0,
            height_scale=float(args.height_scale),
            show_depth_density=not bool(args.no_depth_density),
            year_axis_min=int(adaptive_row["year_min"]),
            year_axis_max=int(adaptive_row["year_max"]),
        )
        rows.extend([depth_1_row, adaptive_row])

    suffix = (
        "_".join(sorted(selected_cases))
        if selected_cases != {case.key for case in CASE_DEFS}
        else ""
    )
    summary_name = (
        f"selected_indicator_sankey_summary_{suffix}.csv"
        if suffix
        else "selected_indicator_sankey_summary.csv"
    )
    summary_path = OUTPUT_DIR / summary_name
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    print(f"\nWrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
