from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from datapackage import Package
from openpyxl import load_workbook

from trails import Trails

DEFAULT_METHODS = [
    "IPCC 2021 (incl. biogenic CO2) - climate change: total (incl. biogenic CO2) - global warming potential (GWP100)"
]


@dataclass(frozen=True)
class ActivityDef:
    name: str
    reference_product: str
    location: str


def _normalize(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _read_activity_and_exchanges(ws) -> tuple[ActivityDef | None, list[dict[str, str]]]:
    activity_name = ""
    reference_product = ""
    location = ""

    for r in range(1, 80):
        k = _normalize(ws.cell(r, 1).value).lower()
        v = _normalize(ws.cell(r, 2).value)
        if k == "activity":
            activity_name = v
        elif k == "reference product":
            reference_product = v
        elif k == "location":
            location = v

    if not activity_name:
        return None, []

    header_row = None
    for r in range(1, 120):
        if (
            _normalize(ws.cell(r, 1).value).lower() == "name"
            and _normalize(ws.cell(r, 10).value).lower() == "type"
        ):
            header_row = r
            break

    exchanges: list[dict[str, str]] = []
    if header_row is not None:
        headers = [_normalize(ws.cell(header_row, c).value) for c in range(1, 21)]
        for r in range(header_row + 1, 400):
            row_vals = [_normalize(ws.cell(r, c).value) for c in range(1, 21)]
            if not any(row_vals):
                continue
            row = {headers[i]: row_vals[i] for i in range(len(headers)) if headers[i]}
            exchanges.append(row)

    return ActivityDef(activity_name, reference_product, location), exchanges


def collect_terminal_activities(inventory_paths: Iterable[Path]) -> list[ActivityDef]:
    all_activities: list[ActivityDef] = []
    activity_names: set[str] = set()
    technosphere_supplier_names: list[str] = []

    for path in inventory_paths:
        wb = load_workbook(path, data_only=True)
        for ws in wb.worksheets:
            activity, exchanges = _read_activity_and_exchanges(ws)
            if activity is None:
                continue
            all_activities.append(activity)
            activity_names.add(activity.name)
            for exc in exchanges:
                if _normalize(exc.get("type", "")).lower() != "technosphere":
                    continue
                supplier_name = _normalize(exc.get("name", ""))
                if supplier_name:
                    technosphere_supplier_names.append(supplier_name)

    suppliers_used_by_imported = {
        supplier
        for supplier in technosphere_supplier_names
        if supplier in activity_names
    }

    terminal = [a for a in all_activities if a.name not in suppliers_used_by_imported]

    # Stable order, deduplicated
    seen: set[tuple[str, str, str]] = set()
    out: list[ActivityDef] = []
    for a in terminal:
        k = (a.name, a.reference_product, a.location)
        if k in seen:
            continue
        seen.add(k)
        out.append(a)
    return out


def _metadata_by_idx(trails: Trails) -> dict[int, dict]:
    if not trails.activity_indices:
        return {}
    first_label = next(iter(trails.activity_indices))
    by_idx = trails.activity_indices[first_label]
    return {int(k): v for k, v in by_idx.items()}


def match_activity_indices(
    trails: Trails, targets: list[ActivityDef]
) -> dict[ActivityDef, int]:
    by_idx = _metadata_by_idx(trails)
    out: dict[ActivityDef, int] = {}

    for target in targets:
        exact = []
        for idx, md in by_idx.items():
            if not isinstance(md, dict):
                continue
            if _normalize(md.get("name")) != target.name:
                continue
            if _normalize(md.get("reference product")) != target.reference_product:
                continue
            if _normalize(md.get("location")) != target.location:
                continue
            exact.append(idx)
        if len(exact) == 1:
            out[target] = exact[0]
            continue

        # Fallback to unique name match
        by_name = [
            idx
            for idx, md in by_idx.items()
            if isinstance(md, dict) and _normalize(md.get("name")) == target.name
        ]
        if len(by_name) == 1:
            out[target] = by_name[0]

    return out


def _to_method_vector(score_obj: object, methods: list[str]) -> np.ndarray:
    arr = np.asarray(score_obj, dtype=float)
    if arr.ndim == 0:
        return np.array([float(arr)], dtype=float)
    if arr.ndim == 1:
        return arr.astype(float, copy=False)
    return arr.ravel().astype(float, copy=False)


def _temporal_total_scores_by_method(trails: Trails, methods: list[str]) -> np.ndarray:
    if trails.scores is None:
        raise RuntimeError("trails.scores is None after temporal LCA run.")
    da = trails.scores
    if "method" in da.dims:
        reduce_dims = [d for d in da.dims if d != "method"]
        vals = da.sum(dim=reduce_dims).values
        return np.asarray(vals, dtype=float).reshape(-1)
    return np.array([float(da.sum().values)], dtype=float)


def run(args: argparse.Namespace) -> int:
    package = Package(str(args.datapackage))
    trails = Trails(
        package,
        interpolate_annual=True,
        interpolation_start_year_offset=args.interpolation_start_year_offset,
        interpolation_end_year_offset=args.interpolation_end_year_offset,
        debug=False,
    )

    inventory_paths = [Path(p).resolve() for p in args.inventories]
    trails.import_excel_inventory([str(p) for p in inventory_paths])

    terminals = collect_terminal_activities(inventory_paths)
    idx_map = match_activity_indices(trails, terminals)

    if not idx_map:
        raise RuntimeError(
            "No terminal imported activities could be matched to indices."
        )

    depths = [1, 5]
    methods = list(args.methods)

    rows: list[dict[str, object]] = []
    total_acts = len(idx_map)
    for act_i, (act_def, idx) in enumerate(idx_map.items(), start=1):
        print(
            f"\n[{act_i}/{total_acts}] Activity idx={idx} | "
            f"{act_def.name} ({act_def.reference_product}, {act_def.location})"
        )

        t0_static = time.perf_counter()
        trails.static_lca(
            year=int(args.reference_year),
            act_idx=int(idx),
            methods=methods,
            amount=float(args.amount),
        )
        static_vec = _to_method_vector(trails.static_score, methods)
        dt_static = time.perf_counter() - t0_static
        print(f"  static_lca done in {dt_static:.2f}s")

        temporal_by_depth: dict[int, np.ndarray] = {}
        for depth in depths:
            print(f"  depth={depth}: temporal_routing + lca ...")
            t0_depth = time.perf_counter()
            trails.temporal_routing(
                start_year=int(args.reference_year),
                start_act_idx=int(idx),
                amount=float(args.amount),
                max_depth=int(depth),
                show_progress=bool(args.show_progress),
                attribute_to_roots=True,
            )
            trails.lca(
                methods=methods,
                show_progress=bool(args.show_progress),
                compute_score=True,
                store_inventory=False,
            )
            temporal_by_depth[depth] = _temporal_total_scores_by_method(trails, methods)
            dt_depth = time.perf_counter() - t0_depth
            print(f"  depth={depth}: done in {dt_depth:.2f}s")

        for m_i, method in enumerate(methods):
            static_score = float(static_vec[m_i if m_i < len(static_vec) else 0])
            d1 = float(
                temporal_by_depth[1][m_i if m_i < len(temporal_by_depth[1]) else 0]
            )
            d5 = float(
                temporal_by_depth[5][m_i if m_i < len(temporal_by_depth[5]) else 0]
            )
            rows.append(
                {
                    "activity_index": int(idx),
                    "activity_name": act_def.name,
                    "reference_product": act_def.reference_product,
                    "location": act_def.location,
                    "method": method,
                    "static_score": static_score,
                    "temporal_depth_1": d1,
                    "temporal_depth_5": d5,
                    "delta_d1_vs_static": d1 - static_score,
                    "delta_d5_vs_d1": d5 - d1,
                }
            )

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("\nTerminal imported activities analyzed:")
    for act_def, idx in sorted(idx_map.items(), key=lambda x: x[1]):
        print(
            f"- idx={idx} | {act_def.name} ({act_def.reference_product}, {act_def.location})"
        )

    print(f"\nWrote results: {out_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare static LCA and temporal LCA scores at routing depths 1 and 5 "
            "for terminal imported activities."
        )
    )
    parser.add_argument(
        "--datapackage",
        type=Path,
        default=Path("/Users/romain/GitHub/premise/dev/trails_2026-02-22.zip"),
        help="Path to Frictionless datapackage zip/json.",
    )
    parser.add_argument(
        "--inventories",
        nargs="+",
        default=[
            str(
                Path(
                    "/Users/romain/GitHub/trails/dev/lci-case-study-daccs_storage_risk.xlsx"
                )
            ),
            str(
                Path(
                    "/Users/romain/GitHub/trails/dev/lci-case-study-fertilizer_n2o_timing.xlsx"
                )
            ),
            str(
                Path(
                    "/Users/romain/GitHub/trails/dev/lci-case-study-marine_fuel_switch.xlsx"
                )
            ),
            str(
                Path(
                    "/Users/romain/GitHub/trails/dev/lci-case-study-biomass_growth_vs_gas_heat.xlsx"
                )
            ),
            str(Path("/Users/romain/GitHub/trails/dev/lci-pass_cars.xlsx")),
        ],
        help="Excel inventory files to import.",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=DEFAULT_METHODS,
        help="LCIA method names.",
    )
    parser.add_argument("--reference-year", type=int, default=2035)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--interpolation-start-year-offset", type=int, default=-20)
    parser.add_argument("--interpolation-end-year-offset", type=int, default=20)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/Users/romain/GitHub/trails/dev/temporal_depth_comparison.csv"),
        help="Output CSV path.",
    )
    return parser


if __name__ == "__main__":
    parser = build_parser()
    raise SystemExit(run(parser.parse_args()))
