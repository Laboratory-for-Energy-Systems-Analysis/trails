from __future__ import annotations

import argparse
import importlib
import importlib.util
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import sparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))


RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DATAPACKAGE = Path("/Users/romain/GitHub/premise/dev/trails_2026-05-30.zip")
LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")
OUTPUT_DIR = REPO_ROOT / "dev" / "notebook_runs" / "polyol_propylene_gwp_diagnostic"

METHOD = (
    "IPCC 2021 (incl. biogenic CO2) - climate change: total "
    "(incl. biogenic CO2) - global warming potential (GWP100)"
)

ACTIVITY_NAME = "polyol precursor production from captured CO2"
REFERENCE_PRODUCT = "polyol precursor"
LOCATION = "RER"

REFERENCE_YEAR = 2025
FUNCTIONAL_UNIT_AMOUNT = 50_000_000_000.0
DEPTH = 5
MIN_AMOUNT = 1e-3


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


def _metadata_by_idx(trails: Any) -> dict[int, dict[str, Any]]:
    if not trails.activity_indices:
        return {}
    first_label = next(iter(trails.activity_indices))
    return {
        int(key): value
        for key, value in trails.activity_indices[first_label].items()
        if isinstance(value, dict)
    }


def _clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def _label(metadata: dict[str, Any]) -> str:
    return (
        f"{metadata.get('name', '')} | "
        f"{metadata.get('reference product', '')} | "
        f"{metadata.get('location', '')}"
    )


def _activity_label(by_idx: dict[int, dict[str, Any]], idx: int) -> str:
    metadata = by_idx.get(int(idx), {})
    if not metadata:
        return f"idx={idx}"
    return f"idx={idx}: {_label(metadata)}"


def _dense_1d(data: Any) -> np.ndarray:
    raw = getattr(data, "data", data)
    if isinstance(raw, sparse.COO):
        dense = raw.todense()
    elif hasattr(raw, "todense"):
        dense = raw.todense()
    else:
        dense = raw
    return np.asarray(dense, dtype=float).reshape(-1)


def _selected_scores(scores: Any, method: str) -> Any:
    data = scores
    if "method" in data.dims:
        data = data.sel(method=method)
    return data


def _root_timeseries(scores: Any, *, root_idx: int, method: str) -> pd.DataFrame:
    data = _selected_scores(scores, method)
    if "root activity" not in data.dims:
        raise RuntimeError(f"Scores have no root attribution dimension: {data.dims}")
    data = data.sel({"root activity": int(root_idx)})
    extra_dims = [dim for dim in data.dims if dim != "year"]
    if extra_dims:
        data = data.sum(dim=extra_dims)
    data = data.transpose("year")
    years = np.asarray(data.coords["year"].values, dtype=int)
    annual = _dense_1d(data)
    return pd.DataFrame(
        {
            "year": years,
            "annual": annual,
            "cumulative": np.cumsum(annual),
        }
    )


def _activity_contributions(
    scores: Any,
    *,
    root_idx: int,
    year: int,
    method: str,
    by_idx: dict[int, dict[str, Any]],
) -> pd.DataFrame:
    data = _selected_scores(scores, method)
    data = data.sel({"root activity": int(root_idx), "year": int(year)})
    extra_dims = [dim for dim in data.dims if dim != "activity"]
    if extra_dims:
        data = data.sum(dim=extra_dims)

    raw = data.data
    if isinstance(raw, sparse.COO):
        if raw.nnz:
            activities = raw.coords[0].astype(int)
            values = raw.data.astype(float)
        else:
            activities = np.array([], dtype=int)
            values = np.array([], dtype=float)
    else:
        dense = _dense_1d(data)
        activities = np.flatnonzero(dense)
        values = dense[activities]

    rows: list[dict[str, Any]] = []
    for activity, value in zip(activities, values):
        metadata = by_idx.get(int(activity), {})
        rows.append(
            {
                "activity": int(activity),
                "annual_score": float(value),
                "name": metadata.get("name", ""),
                "reference_product": metadata.get("reference product", ""),
                "location": metadata.get("location", ""),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["abs_score"] = frame["annual_score"].abs()
    return frame.sort_values("abs_score", ascending=False)


def _routing_graph_summary(trails: Any) -> dict[str, int]:
    graph = trails.graph
    if graph is None:
        return {}
    depths: dict[int, int] = {}
    frontier_nodes = 0
    direct_bio_nodes = 0
    for _, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        depths[depth] = depths.get(depth, 0) + 1
        if float(data.get("frontier_amount") or 0.0):
            frontier_nodes += 1
        if float(data.get("direct_bio_amount") or 0.0):
            direct_bio_nodes += 1
    out = {
        "nodes": int(graph.number_of_nodes()),
        "edges": int(graph.number_of_edges()),
        "max_depth": max(depths) if depths else 0,
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
    }
    for depth, count in sorted(depths.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    return out


def _find_activity(
    by_idx: dict[int, dict[str, Any]],
    *,
    name: str,
    reference_product: str,
    location: str,
) -> int:
    matches = []
    for idx, metadata in by_idx.items():
        if _clean(metadata.get("name")) != name:
            continue
        if _clean(metadata.get("reference product")) != reference_product:
            continue
        if _clean(metadata.get("location")) != location:
            continue
        matches.append(int(idx))
    if len(matches) != 1:
        raise ValueError(f"Expected one activity match, got {matches}")
    return matches[0]


def _find_propylene_root_candidates(
    trails: Any, by_idx: dict[int, dict[str, Any]]
) -> pd.DataFrame:
    root_ids: set[int] = set()
    if trails.graph is not None:
        for _, data in trails.graph.nodes(data=True):
            if int(data.get("depth", -1)) == 1:
                root_ids.add(int(data.get("act_idx")))

    rows: list[dict[str, Any]] = []
    for idx in sorted(root_ids):
        metadata = by_idx.get(idx, {})
        text = " ".join(
            str(metadata.get(key, ""))
            for key in ("name", "reference product", "location")
        ).lower()
        if "propylene oxide" not in text:
            continue
        rows.append(
            {
                "activity": idx,
                "name": metadata.get("name", ""),
                "reference_product": metadata.get("reference product", ""),
                "location": metadata.get("location", ""),
            }
        )
    return pd.DataFrame(rows)


def _direct_propylene_exchange_rows(
    trails: Any,
    *,
    year: int,
    act_idx: int,
    amount: float,
    by_idx: dict[int, dict[str, Any]],
) -> pd.DataFrame:
    context = trails._get_scenario_context(int(year))
    if context is None:
        raise RuntimeError(f"No scenario context for year={year}")
    _, _, t = context

    rows: list[dict[str, Any]] = []
    row = trails.A[int(t), int(act_idx), :]
    for product_index, exchange_value in zip(row.coords[0], row.data):
        product_index = int(product_index)
        if product_index == int(act_idx):
            continue
        metadata = by_idx.get(product_index, {})
        text = " ".join(
            str(metadata.get(key, ""))
            for key in ("name", "reference product", "location")
        ).lower()
        if "propylene oxide" not in text:
            continue

        exchange_value = float(exchange_value)
        child_amount = trails._child_activity_amount(
            t=int(t),
            product_index=product_index,
            parent_amount=float(amount),
            exchange_value=exchange_value,
        )
        tex, offsets_and_weights = trails._get_tech_td_expanded(
            year=int(year),
            act_idx=int(act_idx),
            prod_idx=int(product_index),
            debug=False,
        )
        rows.append(
            {
                "activity": product_index,
                "label": _label(metadata),
                "matrix_value": exchange_value,
                "child_activity_amount": float(child_amount),
                "td_distribution": getattr(tex, "distribution", None),
                "td_loc": getattr(tex, "loc", None),
                "td_scale": getattr(tex, "scale", None),
                "td_min": getattr(tex, "offset_min", None),
                "td_max": getattr(tex, "offset_max", None),
                "td_offsets": getattr(tex, "offsets", None),
                "td_weights": getattr(tex, "weights", None),
                "td_amount_source": getattr(tex, "amount_source", None),
                "expanded_offsets_and_weights": offsets_and_weights,
            }
        )
    return pd.DataFrame(rows)


def _depth1_propylene_nodes(
    trails: Any, by_idx: dict[int, dict[str, Any]]
) -> pd.DataFrame:
    graph = trails.graph
    if graph is None:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for _, data in graph.nodes(data=True):
        if int(data.get("depth", -1)) != 1:
            continue
        idx = int(data.get("act_idx"))
        metadata = by_idx.get(idx, {})
        text = " ".join(
            str(metadata.get(key, ""))
            for key in ("name", "reference product", "location")
        ).lower()
        if "propylene oxide" not in text:
            continue
        rows.append(
            {
                "activity": idx,
                "year": int(data.get("year")),
                "amount": float(data.get("amount") or 0.0),
                "frontier_amount": float(data.get("frontier_amount") or 0.0),
                "direct_bio_amount": float(data.get("direct_bio_amount") or 0.0),
                "label": _label(metadata),
            }
        )
    return pd.DataFrame(rows).sort_values(["activity", "year"])


def _timed(label: str, func: Any, *args: Any, **kwargs: Any) -> Any:
    print(label, flush=True)
    start = time.perf_counter()
    result = func(*args, **kwargs)
    print(f"  done in {time.perf_counter() - start:.1f}s", flush=True)
    return result


def _run_lca_with_fallback(trails: Any, methods: list[str]) -> None:
    try:
        return _timed(
            "Temporal LCA (iterative solver)",
            trails.lca,
            methods=methods,
            show_progress=False,
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode="iterative",
            iterative_rtol=1e-3,
            iterative_atol=0.0,
            iterative_restart=100,
            iterative_maxiter=1000,
            iterative_use_guess=True,
            iterative_preconditioner="jacobi",
            iterative_ilu_drop_tol=1e-4,
            iterative_ilu_fill_factor=10.0,
            ei_version="3.12",
        )
    except RuntimeError as exc:
        print(f"  iterative solver failed: {exc}", flush=True)
        return _timed(
            "Temporal LCA (direct fallback)",
            trails.lca,
            methods=methods,
            show_progress=False,
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode="direct",
            ei_version="3.12",
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run one isolated polyol temporal LCA and diagnose the propylene "
            "oxide root contribution for the IPCC 2021 GWP indicator."
        )
    )
    parser.add_argument("--datapackage", type=Path, default=DATAPACKAGE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--reference-year", type=int, default=REFERENCE_YEAR)
    parser.add_argument("--amount", type=float, default=FUNCTIONAL_UNIT_AMOUNT)
    parser.add_argument("--depth", type=int, default=DEPTH)
    parser.add_argument("--min-amount", type=float, default=MIN_AMOUNT)
    args = parser.parse_args()

    if LCIA_JSON.exists():
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(LCIA_JSON)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    runner = _load_runner()
    runner._validate_paths(args.datapackage, runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS)

    print(f"Datapackage: {args.datapackage}", flush=True)
    print(f"LCIA JSON: {os.environ.get('TRAILS_LCIA_EI312_JSON', '<package default>')}")
    print("Foreground inventories:", flush=True)
    for path in runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS:
        print(f"  - {path}", flush=True)
    print(
        "Run: "
        f"activity={ACTIVITY_NAME!r}, year={args.reference_year}, "
        f"amount={args.amount:g}, depth={args.depth}, "
        f"min_amount={args.min_amount:g}",
        flush=True,
    )

    trails = _timed(
        "Loading Trails and importing OneDrive foreground inventories",
        runner._load_trails,
        datapackage=args.datapackage,
        interpolation_cache_dir=None,
        inventory_paths=[
            path.resolve() for path in runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS
        ],
        import_before_interpolation=False,
        remove_base_temporal_distributions=False,
        no_cache_interpolation=False,
        interpolation_start_year_offset=-20,
        interpolation_end_year_offset=20,
    )

    by_idx = _metadata_by_idx(trails)
    activity_index = _find_activity(
        by_idx,
        name=ACTIVITY_NAME,
        reference_product=REFERENCE_PRODUCT,
        location=LOCATION,
    )
    mapped_year = int(trails._map_year_to_scenario_year(int(args.reference_year)))
    context = trails._get_scenario_context(mapped_year)
    if context is None:
        raise RuntimeError(f"No scenario context for mapped year {mapped_year}")
    _, _, t = context
    production_amount = trails._production_amount(int(t), activity_index)
    activity_amount = trails._activity_amount_from_product_demand(
        int(t), activity_index, float(args.amount)
    )
    print(f"Matched polyol: {_activity_label(by_idx, activity_index)}", flush=True)
    print(
        f"Scenario year={mapped_year}, production_amount={production_amount:.12g}, "
        f"activity_amount={activity_amount:.12g}",
        flush=True,
    )

    static_start = time.perf_counter()
    trails.static_lca(
        year=mapped_year,
        act_idx=activity_index,
        amount=float(args.amount),
        methods=[METHOD],
        ei_version="3.12",
    )
    static_score = trails.static_score
    if isinstance(static_score, list):
        static_score = static_score[0]
    print(
        f"Static score: {float(static_score):.12g} "
        f"(done in {time.perf_counter() - static_start:.1f}s)",
        flush=True,
    )

    direct_exchange_rows = _direct_propylene_exchange_rows(
        trails,
        year=mapped_year,
        act_idx=activity_index,
        amount=activity_amount,
        by_idx=by_idx,
    )
    direct_exchange_rows.to_csv(
        args.output_dir / "direct_polyol_propylene_exchanges.csv", index=False
    )
    print("\nDirect propylene-like exchanges from the polyol foreground:", flush=True)
    if direct_exchange_rows.empty:
        print("  none found", flush=True)
    else:
        print(
            direct_exchange_rows[
                [
                    "activity",
                    "label",
                    "matrix_value",
                    "child_activity_amount",
                    "td_distribution",
                    "td_loc",
                    "expanded_offsets_and_weights",
                ]
            ].to_string(index=False),
            flush=True,
        )

    _timed(
        "\nTemporal routing",
        trails.temporal_routing,
        start_year=int(args.reference_year),
        start_act_idx=activity_index,
        amount=float(args.amount),
        max_depth=int(args.depth),
        min_amount=float(args.min_amount),
        show_progress=False,
        attribute_to_roots=True,
    )
    print("Routing graph summary:", flush=True)
    for key, value in _routing_graph_summary(trails).items():
        print(f"  {key}: {value:,}", flush=True)

    depth1_nodes = _depth1_propylene_nodes(trails, by_idx)
    depth1_nodes.to_csv(args.output_dir / "depth1_propylene_nodes.csv", index=False)
    print("\nDepth-1 propylene-like routed nodes:", flush=True)
    if depth1_nodes.empty:
        print("  none found", flush=True)
    else:
        print(depth1_nodes.to_string(index=False), flush=True)

    _run_lca_with_fallback(trails, [METHOD])
    if trails.scores is None:
        raise RuntimeError("Temporal LCA did not produce scores.")
    print(f"Score dimensions: {trails.scores.dims}", flush=True)

    candidates = _find_propylene_root_candidates(trails, by_idx)
    candidates.to_csv(args.output_dir / "propylene_root_candidates.csv", index=False)
    print("\nPropylene-like root activity candidates:", flush=True)
    if candidates.empty:
        raise RuntimeError("No propylene oxide root activity candidate was found.")
    print(candidates.to_string(index=False), flush=True)

    root_summaries = []
    for root_idx in candidates["activity"].astype(int):
        series = _root_timeseries(trails.scores, root_idx=root_idx, method=METHOD)
        nonzero = series.loc[series["annual"] != 0.0].copy()
        final_cumulative = float(series["annual"].sum())
        min_pre_ref = series.loc[series["year"] < int(args.reference_year), "annual"]
        root_summaries.append(
            {
                "activity": int(root_idx),
                "label": _activity_label(by_idx, int(root_idx)),
                "final_cumulative": final_cumulative,
                "min_pre_reference_annual": (
                    float(min_pre_ref.min()) if not min_pre_ref.empty else 0.0
                ),
                "nonzero_years": int(len(nonzero)),
            }
        )
        series.to_csv(
            args.output_dir / f"root_{int(root_idx)}_timeseries.csv", index=False
        )

    summary = pd.DataFrame(root_summaries).sort_values(
        "final_cumulative", key=lambda s: s.abs(), ascending=False
    )
    summary.to_csv(args.output_dir / "propylene_root_summary.csv", index=False)
    print("\nPropylene root score summary:", flush=True)
    print(summary.to_string(index=False), flush=True)

    selected_root = int(summary.iloc[0]["activity"])
    selected_series = _root_timeseries(
        trails.scores, root_idx=selected_root, method=METHOD
    )
    window = selected_series[
        (selected_series["year"] >= int(args.reference_year) - 10)
        & (selected_series["year"] <= int(args.reference_year) + 10)
        & (selected_series["annual"] != 0.0)
    ]
    print(
        f"\nSelected root for detailed diagnosis: "
        f"{_activity_label(by_idx, selected_root)}",
        flush=True,
    )
    print("\nNon-zero annual values around the reference year:", flush=True)
    if window.empty:
        print("  none", flush=True)
    else:
        print(window.to_string(index=False), flush=True)

    pre_ref = selected_series[
        (selected_series["year"] < int(args.reference_year))
        & (selected_series["annual"] < 0.0)
    ]
    if pre_ref.empty:
        print(
            "\nNo negative pre-reference annual value found for this root.", flush=True
        )
        return

    suspicious = pre_ref.sort_values("annual").iloc[0]
    suspicious_year = int(suspicious["year"])
    suspicious_value = float(suspicious["annual"])
    print(
        f"\nMost negative pre-reference annual value: year={suspicious_year}, "
        f"annual={suspicious_value:.12g}",
        flush=True,
    )

    contributions = _activity_contributions(
        trails.scores,
        root_idx=selected_root,
        year=suspicious_year,
        method=METHOD,
        by_idx=by_idx,
    )
    contributions.to_csv(
        args.output_dir
        / f"root_{selected_root}_{suspicious_year}_activity_contribs.csv",
        index=False,
    )
    if contributions.empty:
        print(
            "No activity-level contributors found for the suspicious year.", flush=True
        )
        return

    print("\nTop negative activity contributors in the suspicious year:", flush=True)
    negative = contributions[contributions["annual_score"] < 0].head(20)
    print(
        negative[
            [
                "activity",
                "annual_score",
                "name",
                "reference_product",
                "location",
            ]
        ].to_string(index=False),
        flush=True,
    )

    print("\nTop positive activity contributors in the same year:", flush=True)
    positive = contributions[contributions["annual_score"] > 0].head(10)
    if positive.empty:
        print("  none", flush=True)
    else:
        print(
            positive[
                [
                    "activity",
                    "annual_score",
                    "name",
                    "reference_product",
                    "location",
                ]
            ].to_string(index=False),
            flush=True,
        )

    print(f"\nWrote diagnostic CSV files to: {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
