from __future__ import annotations


from dataclasses import dataclass
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from datapackage import Package
from trails import Trails
from trails.datapackage import interpolate_to_annual

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

PUBLICATION_DIR = Path(__file__).resolve().parent
DEFAULT_DATAPACKAGE = PUBLICATION_DIR / "trails_remind_SSP2-PkBudg1000.zip"
DEFAULT_LCIA_JSON: Path | None = None


@dataclass(frozen=True)
class ActivityDef:
    name: str
    reference_product: str
    location: str


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _metadata_by_idx(trails: Trails) -> dict[int, dict]:
    if not getattr(trails, "activity_indices", None):
        return {}
    first_label = next(iter(trails.activity_indices))
    return {
        int(key): value for key, value in trails.activity_indices[first_label].items()
    }


def _match_activity_indices(
    trails: Trails,
    targets: list[ActivityDef],
) -> dict[ActivityDef, int]:
    by_idx = _metadata_by_idx(trails)
    out: dict[ActivityDef, int] = {}

    for target in targets:
        exact: list[int] = []
        for index, metadata in by_idx.items():
            if not isinstance(metadata, dict):
                continue
            if _clean(metadata.get("name")) != target.name:
                continue
            if _clean(metadata.get("reference product")) != target.reference_product:
                continue
            if _clean(metadata.get("location")) != target.location:
                continue
            exact.append(index)
        if len(exact) == 1:
            out[target] = exact[0]
            continue

        by_name = [
            index
            for index, metadata in by_idx.items()
            if isinstance(metadata, dict)
            and _clean(metadata.get("name")) == target.name
        ]
        if len(by_name) == 1:
            out[target] = by_name[0]

    return out


def _activity_metadata(
    trails: Trails, activity_index: int, reference_year: int
) -> dict[str, Any]:
    labels: list[str] = []
    try:
        labels.append(str(trails._map_year_to_scenario_year(int(reference_year))))
    except Exception:
        pass
    labels.extend(str(label) for label in getattr(trails, "scenario_labels", []))

    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        mapping = trails.activity_indices.get(label, {})
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)

    for mapping in trails.activity_indices.values():
        meta = mapping.get(int(activity_index))
        if isinstance(meta, dict):
            return dict(meta)

    return {}


def _activity_label(
    trails: Trails,
    activity_index: int,
    reference_year: int = 2025,
    *,
    max_length: int | None = None,
) -> str:
    meta = _activity_metadata(trails, activity_index, reference_year)
    name = _clean(meta.get("name")) or f"Activity {activity_index}"
    reference_product = _clean(meta.get("reference product"))
    location = _clean(meta.get("location"))
    parts = [name]
    if reference_product:
        parts.append(reference_product)
    if location:
        parts.append(location)
    label = " | ".join(parts)
    if max_length is not None:
        return textwrap.shorten(label, width=max_length, placeholder="...")
    return label


def _clear_base_temporal_distributions(trails: Trails) -> tuple[int, int]:
    tech_count = len(trails.temporal_technosphere_exchanges or {})
    bio_count = len(trails.temporal_biosphere_exchanges or {})
    trails.temporal_technosphere_exchanges = {}
    trails.temporal_biosphere_exchanges = {}
    for cache_name in (
        "_td_offsets_cache",
        "_tech_td_cache",
        "_tech_td_expanded_cache",
        "_direct_bio_cache_by_year",
        "_bio_td_expanded_cache",
        "_bio_score_row_char_cache",
        "_bio_score_row_char_matrix_cache",
    ):
        cache = getattr(trails, cache_name, None)
        if hasattr(cache, "clear"):
            cache.clear()
    return tech_count, bio_count


def _interpolate_trails_after_import(
    trails: Trails,
    *,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> None:
    print("Interpolating foreground-augmented matrices to annual resolution", flush=True)
    trails.A, trails.B, trails.scenario_labels, trails.scenario_index = (
        interpolate_to_annual(
            trails.A,
            trails.B,
            trails.scenario_labels,
            value_dtype=trails.value_dtype,
            start_year_offset=int(interpolation_start_year_offset),
            end_year_offset=int(interpolation_end_year_offset),
        )
    )
    trails.years_int = np.array(
        [int(label) for label in trails.scenario_labels],
        dtype=int,
    )
    trails.min_year = int(trails.years_int.min())
    trails.max_year = int(trails.years_int.max())
    for cache_name in (
        "_td_offsets_cache",
        "_tech_td_cache",
        "_tech_td_expanded_cache",
        "_direct_bio_cache_by_year",
    ):
        cache = getattr(trails, cache_name, None)
        if hasattr(cache, "clear"):
            cache.clear()


def _load_trails(
    *,
    datapackage: Path,
    interpolation_cache_dir: Path | None,
    inventory_paths: list[Path],
    import_before_interpolation: bool,
    remove_base_temporal_distributions: bool,
    no_cache_interpolation: bool,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> Trails:
    if interpolation_cache_dir is not None:
        raise ValueError("Explicit interpolation cache loading is not supported here.")

    trails = Trails(
        Package(str(datapackage)),
        interpolate_annual=not bool(import_before_interpolation),
        cache_interpolation=not bool(no_cache_interpolation),
        interpolation_start_year_offset=int(interpolation_start_year_offset),
        interpolation_end_year_offset=int(interpolation_end_year_offset),
    )
    if remove_base_temporal_distributions:
        tech_count, bio_count = _clear_base_temporal_distributions(trails)
        print(
            "Cleared base temporal distributions before LCI import: "
            f"technosphere={tech_count:,}, biosphere={bio_count:,}",
            flush=True,
        )

    print(
        "Importing foreground inventories together: "
        + ", ".join(path.name for path in inventory_paths),
        flush=True,
    )
    trails.import_excel_inventory([str(path) for path in inventory_paths])
    if import_before_interpolation:
        _interpolate_trails_after_import(
            trails,
            interpolation_start_year_offset=int(interpolation_start_year_offset),
            interpolation_end_year_offset=int(interpolation_end_year_offset),
        )
    return trails


DEFAULT_CASE_STUDY_ACTIVITY_KEYS = {
    "bev": ActivityDef(
        "transport, passenger, car, battery electric",
        "transport, passenger, car",
        "RER",
    ),
    "polyol": ActivityDef(
        "polyol precursor from captured CO2",
        "polyol precursor",
        "RER",
    ),
    "marine": ActivityDef(
        "marine freight service, temporal fuel transition",
        "transport service",
        "RER",
    ),
    "daccs": ActivityDef(
        "carbon dioxide, captured, with a solvent-based direct air capture system, 1MtCO2",
        "carbon dioxide, captured",
        "Europe",
    ),
}

DATAPACKAGE = DEFAULT_DATAPACKAGE
LCIA_JSON: Path | None = None
INVENTORY_PATHS = [
    PUBLICATION_DIR / "LCIs" / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-case-study-daccs_storage_risk.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-case-study-marine_fuel_switch.xlsx",
    PUBLICATION_DIR / "LCIs" / "lci-pass_cars.xlsx",
]

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
    if LCIA_JSON is not None and LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(LCIA_JSON)

    from trails.plotting import plot_adaptive_sankey

    case_defs = dict(DEFAULT_CASE_STUDY_ACTIVITY_KEYS)
    case_defs["polyol"] = ActivityDef(
        "polyol precursor production from captured CO2",
        "polyol precursor",
        "RER",
    )
    activities = [case_defs[key] for key in CASE_KEYS]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading Trails context", flush=True)
    trails = _load_trails(
        datapackage=DATAPACKAGE.expanduser().resolve(),
        interpolation_cache_dir=None,
        inventory_paths=[path.expanduser().resolve() for path in INVENTORY_PATHS],
        import_before_interpolation=False,
        remove_base_temporal_distributions=False,
        no_cache_interpolation=False,
        interpolation_start_year_offset=INTERPOLATION_START_YEAR_OFFSET,
        interpolation_end_year_offset=INTERPOLATION_END_YEAR_OFFSET,
    )
    activity_maps = _match_activity_indices(trails, activities)

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
