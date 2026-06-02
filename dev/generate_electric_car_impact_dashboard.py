from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from datapackage import Package

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from trails import Trails, export_impact_dashboard, get_lcia_method_names

DEFAULT_DATAPACKAGE = Path("/Users/romain/GitHub/premise/dev/trails_2026-05-18.zip")
DEFAULT_METHOD = (
    "IPCC 2021 (incl. biogenic CO2) - climate change: total "
    "(incl. biogenic CO2) - global warming potential (GWP100)"
)
ELECTRIC_CAR_NAME = "transport, passenger, car, battery electric"
ELECTRIC_CAR_REFERENCE_PRODUCT = "transport, passenger, car"
ELECTRIC_CAR_LOCATION = "RER"


def _repo_root() -> Path:
    """Return the repository root for this script."""
    return REPO_ROOT


def _resolve_inventory_path(path: str | Path | None) -> Path:
    """Resolve the foreground Excel inventory path."""
    if path is not None:
        candidate = Path(path).expanduser().resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"Inventory workbook not found: {candidate}")
        return candidate

    root = _repo_root()
    candidates = [
        root / "dev" / "lci-pass_cars.xlsx",
        root / "dev" / "lci-pass-cars.xlsx",
        root / "examples" / "lci-pass_cars.xlsx",
        root / "examples" / "lci-pass-cars.xlsx",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        "Could not find lci-pass_cars.xlsx in dev/ or examples/. "
        "Pass --inventory explicitly."
    )


def _norm(value: Any) -> str:
    """Normalize metadata values for matching."""
    return "" if value is None else str(value).strip()


def _find_activity_index(
    trails: Trails,
    *,
    name: str,
    reference_product: str,
    location: str,
) -> int:
    """Find one activity index matching exact metadata."""
    matches: dict[int, dict[str, Any]] = {}
    for _label, mapping in trails.activity_indices.items():
        for idx, meta in (mapping or {}).items():
            if (
                _norm(meta.get("name")) == name
                and _norm(meta.get("reference product")) == reference_product
                and _norm(meta.get("location")) == location
            ):
                matches[int(idx)] = meta

    if len(matches) != 1:
        rows = [
            (
                idx,
                _norm(meta.get("name")),
                _norm(meta.get("reference product")),
                _norm(meta.get("location")),
            )
            for idx, meta in sorted(matches.items())
        ]
        raise ValueError(
            "Expected exactly one battery electric passenger car activity, "
            f"found {len(matches)}: {rows}"
        )
    return next(iter(matches))


def _validate_method(method: str, ei_version: str) -> str:
    """Ensure the requested LCIA method exists in the bundled method list."""
    available = get_lcia_method_names(ei_version=ei_version)
    if method not in available:
        close = [name for name in available if "IPCC 2021" in name][:10]
        raise ValueError(
            f"Method not found for ecoinvent {ei_version}: {method!r}. "
            f"Example available IPCC methods: {close}"
        )
    return method


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate an impact dashboard for the battery electric passenger "
            "car from lci-pass_cars.xlsx."
        )
    )
    parser.add_argument(
        "--datapackage",
        type=Path,
        default=DEFAULT_DATAPACKAGE,
        help=f"TRAILS datapackage zip. Default: {DEFAULT_DATAPACKAGE}",
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        default=None,
        help="Foreground Excel inventory. Defaults to dev/lci-pass_cars.xlsx.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "electric_car_impact_dashboard.html",
        help="Output HTML path.",
    )
    parser.add_argument(
        "--method",
        default=DEFAULT_METHOD,
        help="LCIA method to compute and export.",
    )
    parser.add_argument(
        "--ei-version",
        default="3.11",
        help="Bundled ecoinvent LCIA method version.",
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=2050,
        help="Functional unit start year.",
    )
    parser.add_argument(
        "--amount",
        type=float,
        default=1.0,
        help="Functional unit amount in passenger-kilometers.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=3,
        help="Temporal routing depth.",
    )
    parser.add_argument(
        "--top-locations",
        type=int,
        default=15,
        help="Number of impacting locations retained in the dashboard.",
    )
    parser.add_argument(
        "--top-roots",
        type=int,
        default=8,
        help="Number of root activities retained in the dashboard.",
    )
    parser.add_argument(
        "--top-flows",
        type=int,
        default=10,
        help="Number of flow contributors retained in detail tables.",
    )
    parser.add_argument(
        "--solver-mode",
        choices=("iterative", "direct", "bw2calc"),
        default="iterative",
        help="LCA solver backend.",
    )
    parser.add_argument(
        "--iterative-rtol",
        type=float,
        default=1e-3,
        help="Relative tolerance for iterative solves.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars.",
    )
    return parser.parse_args()


def main() -> None:
    """Generate the electric car impact dashboard."""
    args = parse_args()
    datapackage_path = args.datapackage.expanduser().resolve()
    inventory_path = _resolve_inventory_path(args.inventory)
    output_path = args.output.expanduser().resolve()
    method = _validate_method(args.method, args.ei_version)
    show_progress = not bool(args.no_progress)

    if not datapackage_path.exists():
        raise FileNotFoundError(f"Datapackage not found: {datapackage_path}")

    print(f"Loading datapackage: {datapackage_path}")
    package = Package(str(datapackage_path))
    trails = Trails(package=package, interpolate_annual=True)

    print(f"Importing foreground inventory: {inventory_path}")
    import_summary = trails.import_excel_inventory(str(inventory_path))
    print(f"Import summary: {import_summary}")

    start_act_idx = _find_activity_index(
        trails,
        name=ELECTRIC_CAR_NAME,
        reference_product=ELECTRIC_CAR_REFERENCE_PRODUCT,
        location=ELECTRIC_CAR_LOCATION,
    )
    print(
        "Selected electric car activity: "
        f"{start_act_idx} | {ELECTRIC_CAR_NAME} | {ELECTRIC_CAR_LOCATION}"
    )

    print("Running temporal routing")
    trails.temporal_routing(
        start_year=args.start_year,
        start_act_idx=start_act_idx,
        amount=args.amount,
        max_depth=args.max_depth,
        show_progress=show_progress,
        attribute_to_roots=True,
    )

    print("Running temporal LCA")
    trails.lca(
        methods=[method],
        show_progress=show_progress,
        attribute_to_roots=True,
        store_inventory=True,
        compute_score=True,
        ei_version=args.ei_version,
        solver_mode=args.solver_mode,
        iterative_rtol=args.iterative_rtol,
    )

    print(f"Writing dashboard: {output_path}")
    html_path = export_impact_dashboard(
        trails,
        filename=output_path,
        method=method,
        top_locations=args.top_locations,
        top_roots=args.top_roots,
        top_flows=args.top_flows,
    )
    print(f"Dashboard HTML: {html_path}")
    print(f"Dashboard JSON: {Path(html_path).with_suffix('.json')}")


if __name__ == "__main__":
    main()
