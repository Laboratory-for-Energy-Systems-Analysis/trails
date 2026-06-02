from __future__ import annotations

import argparse
import cProfile
import io
import pstats
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

from datapackage import Package

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from trails import Trails

DATAPACKAGE = REPO_ROOT / "dev" / "trails_2026-05-18.zip"
ONEDRIVE_DATA = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/trails/data"
)

# The polyol inventory has a foreground exchange to DACCS, so import both files
# and leave the other case-study inventories out of this profiling run.
INVENTORY_PATHS = [
    ONEDRIVE_DATA / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    ONEDRIVE_DATA / "lci-case-study-daccs_storage_risk.xlsx",
]

ACTIVITY_NAME = "polyol precursor production from captured CO2"
REFERENCE_PRODUCT = "polyol precursor"
LOCATION = "RER"
REFERENCE_YEAR = 2026
FUNCTIONAL_UNIT_AMOUNT = 50_000_000_000.0
DEFAULT_DEPTHS = (1, 2, 3)
DEFAULT_MIN_AMOUNTS = (1e-3,)


def _normalize(value: object) -> str:
    return "" if value is None else str(value).strip()


def _timed(label: str, func: Any, *args: Any, **kwargs: Any) -> Any:
    print(label, flush=True)
    start = time.perf_counter()
    result = func(*args, **kwargs)
    print(f"  done in {time.perf_counter() - start:.1f}s", flush=True)
    return result


def _match_activity_index(trails: Trails) -> int:
    if not trails.activity_indices:
        raise RuntimeError("No activity metadata available.")
    label = next(iter(trails.activity_indices))
    matches: list[int] = []
    for idx, metadata in trails.activity_indices[label].items():
        if not isinstance(metadata, dict):
            continue
        if _normalize(metadata.get("name")) != ACTIVITY_NAME:
            continue
        if _normalize(metadata.get("reference product")) != REFERENCE_PRODUCT:
            continue
        if _normalize(metadata.get("location")) != LOCATION:
            continue
        matches.append(int(idx))
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one polyol activity match, got "
            f"{len(matches)}: {matches}"
        )
    return matches[0]


def _graph_stats(trails: Trails) -> dict[str, int]:
    graph = trails.graph
    if graph is None:
        return {}
    depths = Counter(int(data.get("depth", 0)) for _, data in graph.nodes(data=True))
    frontier_nodes = sum(
        1
        for _, data in graph.nodes(data=True)
        if float(data.get("frontier_amount") or 0.0) != 0.0
    )
    direct_bio_nodes = sum(
        1
        for _, data in graph.nodes(data=True)
        if float(data.get("direct_bio_amount") or 0.0) != 0.0
    )
    stats = {
        "nodes": int(graph.number_of_nodes()),
        "edges": int(graph.number_of_edges()),
        "max_depth": max(depths) if depths else 0,
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
    }
    for depth, count in sorted(depths.items()):
        stats[f"nodes_depth_{depth}"] = int(count)
    return stats


def _write_profile_outputs(
    profiler: cProfile.Profile,
    output_dir: Path,
    *,
    depth: int,
    sort_by: str = "cumulative",
    limit: int = 60,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    binary_path = output_dir / f"polyol_depth{depth}_routing.prof"
    text_path = output_dir / f"polyol_depth{depth}_routing_cprofile.txt"

    profiler.dump_stats(str(binary_path))
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats(sort_by)
    stats.print_stats(limit)
    text = stream.getvalue()
    text_path.write_text(text)

    print(f"\nWrote binary profile: {binary_path}", flush=True)
    print(f"Wrote text profile: {text_path}", flush=True)
    print(f"\nTop cProfile entries for depth {depth}:", flush=True)
    print(text, flush=True)


def _parse_depths(value: str) -> tuple[int, ...]:
    depths = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not depths:
        raise argparse.ArgumentTypeError("At least one depth is required.")
    if any(depth < 0 for depth in depths):
        raise argparse.ArgumentTypeError("Depth values must be non-negative.")
    return depths


def _parse_min_amounts(value: str) -> tuple[float, ...]:
    amounts = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    if not amounts:
        raise argparse.ArgumentTypeError("At least one min_amount is required.")
    if any(amount < 0.0 for amount in amounts):
        raise argparse.ArgumentTypeError("min_amount values must be non-negative.")
    return amounts


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Profile temporal routing for the polyol case-study activity at "
            "selected depths with a 50 billion kg functional unit and "
            "selected min_amount cutoffs."
        )
    )
    parser.add_argument("--datapackage", type=Path, default=DATAPACKAGE)
    parser.add_argument(
        "--output-dir", type=Path, default=REPO_ROOT / "dev" / "profiling"
    )
    parser.add_argument(
        "--depths",
        type=_parse_depths,
        default=DEFAULT_DEPTHS,
        help="Comma-separated max_depth values to profile. Default: 1,2,3.",
    )
    parser.add_argument(
        "--min-amounts",
        type=_parse_min_amounts,
        default=DEFAULT_MIN_AMOUNTS,
        help="Comma-separated routing min_amount values. Default: 1e-3.",
    )
    parser.add_argument(
        "--cprofile",
        action="store_true",
        help="Collect cProfile output. Disabled by default to avoid heavy overhead.",
    )
    parser.add_argument("--show-progress", action="store_true")
    args = parser.parse_args()

    print(f"Python: {sys.version}", flush=True)
    print(f"Datapackage: {args.datapackage}", flush=True)
    print("Inventories:", flush=True)
    for path in INVENTORY_PATHS:
        print(f"  - {path}", flush=True)
        if not path.exists():
            raise FileNotFoundError(path)

    trails = _timed(
        "Loading Trails datapackage",
        Trails,
        Package(str(args.datapackage)),
        interpolate_annual=True,
        cache_interpolation=True,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
        debug=False,
    )

    import_summary = _timed(
        "Importing isolated foreground inventories",
        trails.import_excel_inventory,
        [str(path) for path in INVENTORY_PATHS],
    )
    print(f"  import summary: {import_summary}", flush=True)

    activity_index = _match_activity_index(trails)
    scenario_year = int(trails._map_year_to_scenario_year(REFERENCE_YEAR))
    context = trails._get_scenario_context(scenario_year)
    if context is None:
        raise RuntimeError(f"No scenario context for {scenario_year}.")
    _, _, t_index = context
    production_amount = trails._production_amount(int(t_index), activity_index)
    activity_amount = trails._activity_amount_from_product_demand(
        int(t_index),
        activity_index,
        FUNCTIONAL_UNIT_AMOUNT,
    )
    print(
        "Matched activity: "
        f"idx={activity_index}, scenario_year={scenario_year}, "
        f"production_amount={production_amount:.12g}, "
        f"activity_amount={activity_amount:.12g}",
        flush=True,
    )

    for min_amount in args.min_amounts:
        print(f"\n=== min_amount={min_amount:g} ===", flush=True)
        for depth in args.depths:
            print(
                "\nProfiling temporal_routing("
                f"depth={depth}, amount={FUNCTIONAL_UNIT_AMOUNT:g}, "
                f"min_amount={min_amount:g}, cprofile={bool(args.cprofile)})",
                flush=True,
            )
            start = time.perf_counter()
            profiler = cProfile.Profile() if args.cprofile else None
            if profiler is not None:
                profiler.enable()
            try:
                trails.temporal_routing(
                    start_year=REFERENCE_YEAR,
                    start_act_idx=activity_index,
                    amount=FUNCTIONAL_UNIT_AMOUNT,
                    max_depth=int(depth),
                    min_amount=float(min_amount),
                    show_progress=bool(args.show_progress),
                    attribute_to_roots=True,
                )
            finally:
                if profiler is not None:
                    profiler.disable()
            elapsed = time.perf_counter() - start

            print(f"Routing elapsed: {elapsed:.1f}s", flush=True)
            print("Graph stats:", flush=True)
            for key, value in _graph_stats(trails).items():
                print(f"  {key}: {value:,}", flush=True)

            if profiler is not None:
                _write_profile_outputs(
                    profiler,
                    args.output_dir,
                    depth=int(depth),
                )


if __name__ == "__main__":
    main()
