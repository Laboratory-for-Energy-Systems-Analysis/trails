from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import json
import os
import pickle
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from textwrap import shorten, wrap
from typing import Any

import numpy as np
import sparse
from datapackage import Package
import plotly.graph_objects as go

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

from trails import Trails, get_lcia_method_names, plot_temporal_scores
from trails.datapackage import interpolate_to_annual
from trails.plotting import plot_temporal_graph

HELPER_PATH = REPO_ROOT / "dev" / "plot_temporal_lca_depths.py"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "dev" / "terminal_lci_td_comparison_ef31_depth5"
DEFAULT_RESULTS_CSV = DEFAULT_OUTPUT_DIR / "terminal_lci_td_comparison_scores.csv"
DEFAULT_ONEDRIVE_LCI_DIR = Path(
    "/Users/romain/Library/CloudStorage/OneDrive-PaulScherrerInstitut/trails/data"
)
DEFAULT_ONEDRIVE_INVENTORY_PATHS = [
    DEFAULT_ONEDRIVE_LCI_DIR / "lci-case-study-ccu_polyol_delayed_release.xlsx",
    DEFAULT_ONEDRIVE_LCI_DIR / "lci-case-study-daccs_storage_risk.xlsx",
    DEFAULT_ONEDRIVE_LCI_DIR / "lci-case-study-marine_fuel_switch.xlsx",
    DEFAULT_ONEDRIVE_LCI_DIR / "lci-pass_cars.xlsx",
]
EF_V31_HEADLINE_METHODS = [
    "EF v3.1 - acidification - accumulated exceedance (AE)",
    "EF v3.1 - climate change - global warming potential (GWP100)",
    "EF v3.1 - ecotoxicity: freshwater - comparative toxic unit for ecosystems (CTUe)",
    "EF v3.1 - energy resources: non-renewable - abiotic depletion potential (ADP): fossil fuels",
    "EF v3.1 - eutrophication: freshwater - fraction of nutrients reaching freshwater end compartment (P)",
    "EF v3.1 - eutrophication: marine - fraction of nutrients reaching marine end compartment (N)",
    "EF v3.1 - eutrophication: terrestrial - accumulated exceedance (AE)",
    "EF v3.1 - human toxicity: carcinogenic - comparative toxic unit for human (CTUh)",
    "EF v3.1 - human toxicity: non-carcinogenic - comparative toxic unit for human (CTUh)",
    "EF v3.1 - ionising radiation: human health - human exposure efficiency relative to u235",
    "EF v3.1 - land use - soil quality index",
    "EF v3.1 - material resources: metals/minerals - abiotic depletion potential (ADP): elements (ultimate reserves)",
    "EF v3.1 - ozone depletion - ozone depletion potential (ODP)",
    "EF v3.1 - particulate matter formation - impact on human health",
    "EF v3.1 - photochemical oxidant formation: human health - tropospheric ozone concentration increase",
    "EF v3.1 - water use - user deprivation potential (deprivation-weighted water consumption)",
]


def _load_depth_helpers():
    spec = importlib.util.spec_from_file_location(
        "plot_temporal_lca_depths", HELPER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load helper script: {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


helpers = _load_depth_helpers()

DEFAULT_CASE_STUDY_ACTIVITIES = [
    helpers.ActivityDef(
        "transport, passenger, car, battery electric",
        "transport, passenger, car",
        "RER",
    ),
    helpers.ActivityDef(
        "polyol precursor from captured CO2",
        "polyol precursor",
        "RER",
    ),
    helpers.ActivityDef(
        "marine freight service, temporal fuel transition",
        "transport service",
        "RER",
    ),
    helpers.ActivityDef(
        "carbon dioxide, captured, with a solvent-based direct air capture system, 1MtCO2",
        "carbon dioxide, captured",
        "Europe",
    ),
]
DEFAULT_CASE_STUDY_ACTIVITY_KEYS = {
    "bev": DEFAULT_CASE_STUDY_ACTIVITIES[0],
    "polyol": DEFAULT_CASE_STUDY_ACTIVITIES[1],
    "marine": DEFAULT_CASE_STUDY_ACTIVITIES[2],
    "daccs": DEFAULT_CASE_STUDY_ACTIVITIES[3],
}
DEFAULT_CASE_STUDY_ROUTING_AMOUNT_SCALES = {
    DEFAULT_CASE_STUDY_ACTIVITY_KEYS["bev"]: 150_000.0,
    DEFAULT_CASE_STUDY_ACTIVITY_KEYS["polyol"]: 50_000_000_000.0,
    DEFAULT_CASE_STUDY_ACTIVITY_KEYS["daccs"]: 20_000_000_000.0,
    DEFAULT_CASE_STUDY_ACTIVITY_KEYS["marine"]: 180_000_000_000.0,
}


@dataclass(frozen=True)
class RunConfig:
    key: str
    label: str
    remove_base_temporal_distributions: bool


@dataclass
class RunResult:
    run_key: str
    run_label: str
    activity_index: int
    activity_label: str
    graph_html: Path | None
    static_scores: dict[str, float]
    years_by_method: dict[str, np.ndarray]
    annual_by_method: dict[str, np.ndarray]
    cumulative_by_method: dict[str, np.ndarray]
    total_by_method: dict[str, float]


RUN_CONFIGS = [
    RunConfig(
        key="all_td",
        label="All TDs",
        remove_base_temporal_distributions=False,
    ),
    RunConfig(
        key="foreground_td_only",
        label="Foreground TDs only",
        remove_base_temporal_distributions=True,
    ),
]

RUN_COLORS = {
    "all_td": {
        "annual": "#2563eb",
        "cumulative": "#111827",
    },
    "foreground_td_only": {
        "annual": "#f97316",
        "cumulative": "#92400e",
    },
}


class _Heartbeat:
    def __init__(self, label: str, interval: float = 20.0) -> None:
        self.label = label
        self.interval = float(interval)
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._start = time.perf_counter()

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        self._thread.join(timeout=1.0)

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            elapsed = time.perf_counter() - self._start
            print(f"{self.label}: still running ({elapsed:.0f}s)", flush=True)


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


def _dtype_from_cache(value: object, default: str) -> np.dtype:
    text = str(value or default)
    if "float32" in text:
        return np.dtype("float32")
    if "float64" in text:
        return np.dtype("float64")
    if "int32" in text:
        return np.dtype("int32")
    if "int64" in text:
        return np.dtype("int64")
    return np.dtype(text)


def _trails_from_interpolation_cache(
    cache_dir: Path,
    *,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> Trails:
    cache_dir = Path(cache_dir).expanduser().resolve()
    meta_path = cache_dir / "meta.json"
    a_path = cache_dir / "A.npz"
    b_path = cache_dir / "B.npz"
    temporal_path = cache_dir / "temporal.pkl"
    indices_path = cache_dir / "indices.pkl"
    missing = [
        path
        for path in (meta_path, a_path, b_path, temporal_path, indices_path)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Interpolation cache is incomplete:\n- "
            + "\n- ".join(str(path) for path in missing)
        )

    with meta_path.open("r", encoding="utf-8") as handle:
        meta = json.load(handle)
    with temporal_path.open("rb") as handle:
        temporal_payload = pickle.load(handle)
    with indices_path.open("rb") as handle:
        indices_payload = pickle.load(handle)

    trails = object.__new__(Trails)
    trails.package = SimpleNamespace(
        descriptor={"name": f"interpolation-cache:{cache_dir.name}"},
        resources=[],
        path=str(cache_dir),
        basepath=str(cache_dir),
    )
    trails.value_dtype = _dtype_from_cache(meta.get("value_dtype"), "float32")
    trails.index_dtype = _dtype_from_cache(meta.get("index_dtype"), "int32")
    trails.debug = False
    trails.interpolation_start_year_offset = int(interpolation_start_year_offset)
    trails.interpolation_end_year_offset = int(interpolation_end_year_offset)
    trails.scenario_labels = list(meta.get("scenario_labels", []))
    trails.scenario_index = {
        label: pos for pos, label in enumerate(trails.scenario_labels)
    }
    trails.template_labels = list(meta.get("template_labels", trails.scenario_labels))
    trails.A = sparse.load_npz(a_path)
    trails.B = sparse.load_npz(b_path)
    trails.temporal_technosphere_exchanges = temporal_payload.get(
        "temporal_technosphere_exchanges",
        {},
    )
    trails.temporal_biosphere_exchanges = temporal_payload.get(
        "temporal_biosphere_exchanges",
        {},
    )
    trails.activity_indices = indices_payload.get("activity_indices", {})
    trails.biosphere_indices = indices_payload.get("biosphere_indices", {})

    trails.template_years_int = np.array(
        [int(label) for label in trails.template_labels],
        dtype=int,
    )
    trails.years_int = np.array(
        [int(label) for label in trails.scenario_labels],
        dtype=int,
    )
    trails.min_year = int(trails.years_int.min())
    trails.max_year = int(trails.years_int.max())

    trails.inventory = None
    trails.characterized_inventory = None
    trails.static_score = None
    trails._inventory_years = None
    trails._inventory_year_index = {}
    trails._inventory_coords = None
    trails._inventory_data = None
    trails.provenance = None
    trails.scores = None
    trails._score_years = None
    trails._score_year_index = {}
    trails.graph = None
    trails._routing_attribute_to_roots = None
    trails._routing_params = None
    trails._td_offsets_cache = {}
    trails._tech_td_cache = {}
    trails._A_row_cache = {}
    trails._direct_bio_cache_by_year = {}
    trails._tech_td_expanded_cache = {}
    print(f"Loaded interpolated matrices from explicit cache: {cache_dir}", flush=True)
    return trails


def _interpolate_trails_after_import(
    trails: Trails,
    *,
    interpolation_start_year_offset: int,
    interpolation_end_year_offset: int,
) -> None:
    print(
        "Interpolating foreground-augmented matrices to annual resolution", flush=True
    )
    with _Heartbeat("  annual interpolation"):
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
    trails._td_offsets_cache.clear()
    trails._tech_td_cache.clear()
    trails._tech_td_expanded_cache.clear()
    trails._direct_bio_cache_by_year.clear()


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
    if interpolation_cache_dir is None and import_before_interpolation:
        trails = Trails(
            Package(str(datapackage)),
            interpolate_annual=False,
            cache_interpolation=False,
            interpolation_start_year_offset=int(interpolation_start_year_offset),
            interpolation_end_year_offset=int(interpolation_end_year_offset),
        )
    elif interpolation_cache_dir is None:
        trails = Trails(
            Package(str(datapackage)),
            interpolate_annual=True,
            cache_interpolation=not bool(no_cache_interpolation),
            interpolation_start_year_offset=int(interpolation_start_year_offset),
            interpolation_end_year_offset=int(interpolation_end_year_offset),
        )
    else:
        trails = _trails_from_interpolation_cache(
            interpolation_cache_dir,
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
    with _Heartbeat("  foreground inventory import"):
        trails.import_excel_inventory([str(path) for path in inventory_paths])
    if interpolation_cache_dir is None and import_before_interpolation:
        _interpolate_trails_after_import(
            trails,
            interpolation_start_year_offset=int(interpolation_start_year_offset),
            interpolation_end_year_offset=int(interpolation_end_year_offset),
        )
    return trails


def _validate_paths(datapackage: Path, inventory_paths: list[Path]) -> None:
    if not datapackage.exists():
        raise FileNotFoundError(f"Datapackage not found: {datapackage}")
    if not inventory_paths:
        raise ValueError("At least one Excel inventory is required.")
    missing = [path for path in inventory_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Inventory file(s) not found:\n- "
            + "\n- ".join(str(path) for path in missing)
        )


def _validate_activity_index(trails: Trails, activity_index: int) -> None:
    if trails.A is None:
        raise RuntimeError("Trails.A is not initialized.")
    n_activities = int(trails.A.shape[1])
    if int(activity_index) < 0 or int(activity_index) >= n_activities:
        raise ValueError(
            f"Activity index {activity_index} is outside the A matrix activity "
            f"axis with length {n_activities}."
        )


def _method_label(method: str) -> str:
    if method.startswith("EF v3.1 - "):
        return method.removeprefix("EF v3.1 - ")
    return method


def _wrapped(value: str, width: int) -> str:
    lines = wrap(
        str(value), width=width, break_long_words=False, break_on_hyphens=False
    )
    return "<br>".join(lines) if lines else str(value)


def _window_series(
    years: np.ndarray,
    values: np.ndarray,
    *,
    year_start: int,
    year_end: int,
    hold_value: bool,
) -> tuple[list[int], list[float]]:
    years = np.asarray(years, dtype=int)
    values = np.asarray(values, dtype=float)
    if years.size == 0:
        fill = 0.0
        return [int(year_start), int(year_end)], [fill, fill]

    order = np.argsort(years)
    years = years[order]
    values = values[order]

    mask = (years >= int(year_start)) & (years <= int(year_end))
    x = years[mask].astype(int).tolist()
    y = values[mask].astype(float).tolist()

    if hold_value:
        before_start = np.flatnonzero(years <= int(year_start))
        start_value = float(values[before_start[-1]]) if before_start.size else 0.0
        before_end = np.flatnonzero(years <= int(year_end))
        end_value = float(values[before_end[-1]]) if before_end.size else start_value
    else:
        start_value = 0.0
        end_value = 0.0

    if not x or x[0] > int(year_start):
        x.insert(0, int(year_start))
        y.insert(0, start_value)
    elif x[0] == int(year_start):
        y[0] = start_value if hold_value else y[0]

    if x[-1] < int(year_end):
        x.append(int(year_end))
        y.append(end_value)
    elif x[-1] == int(year_end) and hold_value:
        y[-1] = end_value

    return x, y


def _format_score(value: float) -> str:
    if not np.isfinite(float(value)):
        return "nan"
    return f"{float(value):.6g}"


def _hide_stacked_area_lines(fig: go.Figure) -> None:
    for trace in fig.data:
        fill = str(getattr(trace, "fill", "") or "")
        if getattr(trace, "stackgroup", None) or fill in {"tozeroy", "tonexty"}:
            trace.update(mode="none", line=dict(width=0))


def _retitle_all_td_cumulative_trace(fig: go.Figure) -> None:
    for trace in fig.data:
        if str(getattr(trace, "name", "") or "") == "Cumulative total":
            trace.update(
                name="Foreground + background TDs",
                legendgroup="cumulative-all-td",
                line=dict(width=3, color="#111827"),
            )


def _extend_secondary_lines_to_window(
    fig: go.Figure,
    *,
    year_start: int,
    year_end: int,
) -> None:
    for trace in fig.data:
        if getattr(trace, "yaxis", None) != "y2":
            continue

        x_raw = getattr(trace, "x", None)
        y_raw = getattr(trace, "y", None)
        if x_raw is None or y_raw is None:
            continue
        x = np.asarray(list(x_raw), dtype=int)
        y = np.asarray(list(y_raw), dtype=float)
        if x.size == 0 or y.size == 0:
            continue

        name = str(getattr(trace, "name", "") or "")
        if name.startswith("Static score"):
            static_value = float(y[0])
            trace.x = [int(year_start), int(year_end)]
            trace.y = [static_value, static_value]
            continue

        clipped_x, clipped_y = _window_series(
            x,
            y,
            year_start=int(year_start),
            year_end=int(year_end),
            hold_value=True,
        )
        trace.x = clipped_x
        trace.y = clipped_y


def _secondary_axis_values(fig: go.Figure) -> np.ndarray:
    values: list[float] = []
    for trace in fig.data:
        if getattr(trace, "yaxis", None) != "y2":
            continue
        y_raw = getattr(trace, "y", None)
        if y_raw is None:
            continue
        values.extend(float(value) for value in y_raw if np.isfinite(float(value)))
    return np.asarray(values, dtype=float)


def _primary_axis_visible_values(
    fig: go.Figure,
    *,
    year_start: int,
    year_end: int,
) -> np.ndarray:
    values: list[float] = []
    for trace in fig.data:
        yaxis = getattr(trace, "yaxis", None)
        if yaxis not in (None, "y", "y1"):
            continue

        x_raw = getattr(trace, "x", None)
        y_raw = getattr(trace, "y", None)
        if x_raw is None or y_raw is None:
            continue

        for x_value, y_value in zip(x_raw, y_raw):
            try:
                year = int(x_value)
                value = float(y_value)
            except (TypeError, ValueError):
                continue
            if int(year_start) <= year <= int(year_end) and np.isfinite(value):
                values.append(value)

    return np.asarray(values, dtype=float)


def _set_primary_axis_to_visible_value_range(
    fig: go.Figure,
    *,
    year_start: int,
    year_end: int,
) -> None:
    values = _primary_axis_visible_values(
        fig,
        year_start=int(year_start),
        year_end=int(year_end),
    )
    if values.size == 0:
        return

    y_min = float(np.nanmin(values))
    y_max = float(np.nanmax(values))
    if not np.isfinite(y_min) or not np.isfinite(y_max):
        return
    if y_min == y_max:
        pad = max(abs(y_min), 1.0) * 0.05
        y_min -= pad
        y_max += pad

    fig.update_layout(yaxis=dict(range=[y_min, y_max]))


def _align_secondary_axis_to_visible_values(fig: go.Figure) -> None:
    values = _secondary_axis_values(fig)
    if values.size == 0:
        return

    primary_range = getattr(fig.layout.yaxis, "range", None)
    if primary_range is None or len(primary_range) != 2:
        return

    y1_min = float(primary_range[0])
    y1_max = float(primary_range[1])
    if not np.isfinite(y1_min) or not np.isfinite(y1_max) or y1_max == y1_min:
        return

    y2_min_data = min(float(np.nanmin(values)), 0.0)
    y2_max_data = max(float(np.nanmax(values)), 0.0)
    if y2_max_data == y2_min_data:
        y2_max_data = y2_min_data + 1.0

    zero_fraction = (0.0 - y1_min) / (y1_max - y1_min)
    headroom = 1.05

    if zero_fraction <= 0.0:
        y2_range = [0.0, y2_max_data * headroom]
    elif zero_fraction >= 1.0:
        y2_range = [y2_min_data * headroom, 0.0]
    else:
        y2_max_eff = y2_max_data
        if y2_min_data < 0.0:
            y2_max_eff = max(
                y2_max_eff,
                -y2_min_data * (1.0 - zero_fraction) / zero_fraction,
            )
        y2_max_eff *= headroom
        y2_min_eff = -(zero_fraction / (1.0 - zero_fraction)) * y2_max_eff
        y2_min_eff = min(y2_min_eff, y2_min_data)
        y2_range = [y2_min_eff, y2_max_eff]

    fig.update_layout(yaxis2=dict(range=y2_range))


def _run_activity(
    *,
    trails: Trails,
    run_config: RunConfig,
    activity_index: int,
    methods: list[str],
    depth: int,
    reference_year: int,
    amount: float,
    show_progress: bool,
    solver_mode: str,
    iterative_rtol: float,
    iterative_maxiter: int | None,
    iterative_restart: int | None,
    ei_version: str,
    fallback_solver_mode: str | None,
    routing_min_amount: float,
    attribute_to_roots: bool,
    write_graph_html: bool,
    graph_output_dir: Path,
    graph_run: str,
    graph_min_edge_amount: float,
) -> RunResult:
    _validate_activity_index(trails, activity_index)
    label = helpers._activity_label(trails, activity_index, int(reference_year))

    print(
        f"  {run_config.label}: static_lca for activity {activity_index}",
        flush=True,
    )
    static_t0 = time.perf_counter()
    trails.static_lca(
        year=int(reference_year),
        act_idx=int(activity_index),
        methods=methods,
        amount=float(amount),
        ei_version=str(ei_version),
    )
    static_scores = helpers._score_to_method_map(trails.static_score, methods)
    print(f"    static_lca done in {time.perf_counter() - static_t0:.1f}s", flush=True)

    print(
        f"    routing start: depth={depth}, min_amount={routing_min_amount:g}, "
        f"attribute_to_roots={bool(attribute_to_roots)}",
        flush=True,
    )
    routing_t0 = time.perf_counter()
    trails.temporal_routing(
        start_year=int(reference_year),
        start_act_idx=int(activity_index),
        amount=float(amount),
        max_depth=int(depth),
        min_amount=float(routing_min_amount),
        show_progress=bool(show_progress),
        attribute_to_roots=bool(attribute_to_roots),
    )
    print(
        "    routing done in "
        f"{time.perf_counter() - routing_t0:.1f}s "
        f"({helpers._routing_graph_summary(trails)})",
        flush=True,
    )

    graph_html: Path | None = None
    should_write_graph = bool(write_graph_html) and str(graph_run) in {
        str(run_config.key),
        "both",
    }
    if should_write_graph:
        graph_html = _expected_graph_path(
            output_dir=graph_output_dir,
            activity_index=int(activity_index),
            activity_label=label,
            run_key=run_config.key,
            depth=int(depth),
        )
        graph_html.parent.mkdir(parents=True, exist_ok=True)
        print(
            "    writing routing graph HTML: "
            f"{graph_html} (min_edge_amount={graph_min_edge_amount:g})",
            flush=True,
        )
        graph_t0 = time.perf_counter()
        plot_temporal_graph(
            trails,
            min_edge_amount=float(graph_min_edge_amount),
            filename=str(graph_html),
            height="900px",
            width="100%",
            physics=False,
        )
        print(
            f"    graph HTML done in {time.perf_counter() - graph_t0:.1f}s",
            flush=True,
        )

    print(
        f"    lca start: {len(methods)} method(s), solver_mode={solver_mode}, "
        f"attribute_to_roots={bool(attribute_to_roots)}",
        flush=True,
    )
    lca_t0 = time.perf_counter()
    try:
        trails.lca(
            methods=methods,
            show_progress=bool(show_progress),
            attribute_to_roots=bool(attribute_to_roots),
            compute_score=True,
            store_inventory=False,
            solver_mode=str(solver_mode),
            iterative_rtol=float(iterative_rtol),
            iterative_maxiter=iterative_maxiter,
            iterative_restart=iterative_restart,
            ei_version=str(ei_version),
        )
        print(f"    lca done in {time.perf_counter() - lca_t0:.1f}s", flush=True)
    except RuntimeError as exc:
        if (
            str(solver_mode) != "iterative"
            or fallback_solver_mode is None
            or fallback_solver_mode == str(solver_mode)
            or "GMRES failed to converge" not in str(exc)
        ):
            raise
        print(
            "    iterative solver did not converge; retrying with "
            f"solver_mode={fallback_solver_mode}",
            flush=True,
        )
        fallback_t0 = time.perf_counter()
        trails.lca(
            methods=methods,
            show_progress=bool(show_progress),
            attribute_to_roots=bool(attribute_to_roots),
            compute_score=True,
            store_inventory=False,
            solver_mode=str(fallback_solver_mode),
            ei_version=str(ei_version),
        )
        print(
            "    fallback lca done in " f"{time.perf_counter() - fallback_t0:.1f}s",
            flush=True,
        )

    print("    reducing temporal scores", flush=True)
    scores = helpers._scores_reduced_to_method_root_year(trails)
    years_by_method: dict[str, np.ndarray] = {}
    annual_by_method: dict[str, np.ndarray] = {}
    cumulative_by_method: dict[str, np.ndarray] = {}
    total_by_method: dict[str, float] = {}

    for method in methods:
        years, annual = helpers._temporal_score_by_year_from_scores(scores, method)
        annual = np.asarray(annual, dtype=float)
        cumulative = np.cumsum(annual)
        years_by_method[method] = np.asarray(years, dtype=int)
        annual_by_method[method] = annual
        cumulative_by_method[method] = cumulative
        total_by_method[method] = float(cumulative[-1]) if cumulative.size else 0.0

    return RunResult(
        run_key=run_config.key,
        run_label=run_config.label,
        activity_index=int(activity_index),
        activity_label=label,
        graph_html=graph_html,
        static_scores=static_scores,
        years_by_method=years_by_method,
        annual_by_method=annual_by_method,
        cumulative_by_method=cumulative_by_method,
        total_by_method=total_by_method,
    )


def _plot_comparison(
    *,
    all_td_trails: Trails,
    activity_index: int,
    activity_label: str,
    method: str,
    results: dict[str, RunResult],
    reference_year: int,
    year_start: int,
    year_end: int,
    output_path: Path,
    width: int,
    height: int,
) -> None:
    all_td = results["all_td"]
    fg_td = results["foreground_td_only"]
    static_all = float(all_td.static_scores[method])
    total_all = float(all_td.total_by_method[method])
    total_fg = float(fg_td.total_by_method[method])
    delta = total_all - total_fg

    method_title = _wrapped(_method_label(method), 92)
    activity_title = _wrapped(activity_label, 88)
    subtitle = (
        f"Activity {activity_index}<br>"
        f"{activity_title}<br>"
        f"{method_title}<br>"
        "Cumulative: "
        f"foreground + background TDs={_format_score(total_all)}; "
        f"foreground TDs only={_format_score(total_fg)}; "
        f"delta={_format_score(delta)}; "
        f"static={_format_score(static_all)}"
    )
    fig = plot_temporal_scores(
        trails=all_td_trails,
        title="",
        method_label="Annual impact",
        method=method,
        cumulative=False,
        stacked=False,
        legend_top_n=7,
        show_flow_contributions=False,
        width=int(width),
        height=int(height),
        year_tick=5,
        year_range=None,
        reference_year=int(reference_year),
        show_cumulative_axis=True,
        cumulative_axis_label="Cumulative / static impact",
        show_cumulative_in_legend=True,
        static_score=static_all,
        static_score_label="Static score",
        static_score_dash="dash",
        static_score_color="#dc2626",
    )
    if isinstance(fig, list):
        raise TypeError("Expected one Plotly figure for the selected method.")

    _hide_stacked_area_lines(fig)
    _retitle_all_td_cumulative_trace(fig)

    fg_x, fg_y = _window_series(
        fg_td.years_by_method[method],
        fg_td.cumulative_by_method[method],
        year_start=year_start,
        year_end=year_end,
        hold_value=True,
    )
    fig.add_trace(
        go.Scatter(
            x=fg_x,
            y=fg_y,
            mode="lines",
            name="Foreground TDs only",
            legendgroup="cumulative-foreground-td-only",
            yaxis="y2",
            line=dict(color=RUN_COLORS["foreground_td_only"]["cumulative"], width=3),
            hovertemplate=(
                "<b>Foreground TDs only</b><br>"
                "Year: %{x}<br>"
                "Cumulative impact: %{y:.6g}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        annotations=[
            dict(
                text=subtitle,
                x=0.5,
                y=1.17,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="center",
                font=dict(size=18, color="#1f3557"),
            )
        ],
        width=int(width),
        height=int(height),
        title=dict(text="", x=0.5, xanchor="center"),
        font=dict(size=16, color="#1f3557"),
        title_font=dict(size=22, color="#1f3557"),
        legend=dict(
            orientation="h",
            x=0,
            xanchor="left",
            y=1.02,
            yanchor="bottom",
            font=dict(size=12),
            groupclick="togglegroup",
        ),
        margin=dict(l=80, r=85, t=215, b=55),
    )
    fig.update_xaxes(
        title_text="Year",
        range=[int(year_start), int(year_end)],
        tick0=int(year_start),
        dtick=5,
        tickmode="linear",
        tickfont=dict(size=14),
        title_font=dict(size=16),
        zeroline=False,
        showgrid=True,
        gridcolor="#e6eef8",
    )
    fig.update_yaxes(
        title_text="Annual impact",
        tickfont=dict(size=14),
        title_font=dict(size=16),
        zeroline=True,
        zerolinecolor="#64748b",
        showgrid=True,
        gridcolor="#e6eef8",
    )
    fig.update_layout(
        yaxis2=dict(
            title_text="Cumulative / static impact",
            tickfont=dict(size=14),
            title_font=dict(size=16),
            zeroline=True,
            zerolinecolor="#64748b",
            showgrid=False,
        )
    )
    _extend_secondary_lines_to_window(
        fig,
        year_start=int(year_start),
        year_end=int(year_end),
    )
    _set_primary_axis_to_visible_value_range(
        fig,
        year_start=int(year_start),
        year_end=int(year_end),
    )
    _align_secondary_axis_to_visible_values(fig)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(
        fig.to_image(
            format="png",
            width=int(width),
            height=int(height),
        )
    )


def _write_rows(csv_path: Path, rows: list[dict[str, Any]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "activity_index_all_td",
        "activity_index_foreground_td_only",
        "activity",
        "method",
        "static_all_td",
        "static_foreground_td_only",
        "temporal_cumulative_all_td",
        "temporal_cumulative_foreground_td_only",
        "temporal_cumulative_delta",
        "graph_html_all_td",
        "graph_html_foreground_td_only",
        "figure",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(csv_path: Path) -> list[dict[str, Any]]:
    if not csv_path.exists():
        return []
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _expected_figure_path(
    *,
    output_dir: Path,
    activity_index: int,
    activity_label: str,
    method: str,
    depth: int,
) -> Path:
    activity_slug = helpers._slugify(
        f"idx_{int(activity_index)}_{activity_label}",
        max_length=140,
    )
    return (
        output_dir
        / activity_slug
        / f"{helpers._slugify(method, max_length=170)}_td_comparison_depth_{int(depth)}.png"
    )


def _expected_graph_path(
    *,
    output_dir: Path,
    activity_index: int,
    activity_label: str,
    run_key: str,
    depth: int,
) -> Path:
    activity_slug = helpers._slugify(
        f"idx_{int(activity_index)}_{activity_label}",
        max_length=140,
    )
    return (
        output_dir
        / activity_slug
        / "graphs"
        / f"{helpers._slugify(run_key, max_length=60)}_routing_depth_{int(depth)}.html"
    )


def _activity_complete(
    rows: list[dict[str, Any]],
    *,
    output_dir: Path,
    activity_index: int,
    activity_label: str,
    methods: list[str],
    depth: int,
) -> bool:
    row_methods = {
        row.get("method")
        for row in rows
        if str(row.get("activity_index_all_td")) == str(int(activity_index))
    }
    if any(method not in row_methods for method in methods):
        return False
    return all(
        _expected_figure_path(
            output_dir=output_dir,
            activity_index=int(activity_index),
            activity_label=activity_label,
            method=method,
            depth=int(depth),
        ).exists()
        for method in methods
    )


def run(args: argparse.Namespace) -> int:
    if args.lcia_json is not None:
        lcia_json = Path(args.lcia_json).expanduser().resolve()
        if not lcia_json.exists():
            raise FileNotFoundError(f"LCIA JSON not found: {lcia_json}")
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(lcia_json)

    datapackage = Path(args.datapackage).expanduser().resolve()
    selected_inventories = list(args.inventories)
    if bool(args.all_dev_lci_inventories):
        selected_inventories.extend(helpers.DEFAULT_INVENTORY_PATHS)
    if bool(args.all_onedrive_lci_inventories):
        selected_inventories.extend(DEFAULT_ONEDRIVE_INVENTORY_PATHS)
    inventory_paths = [
        Path(path).expanduser().resolve() for path in selected_inventories
    ]
    inventory_paths = list(dict.fromkeys(inventory_paths))
    _validate_paths(datapackage, inventory_paths)

    available_methods = get_lcia_method_names(ei_version=str(args.ei_version))
    if args.methods:
        methods = [str(method) for method in args.methods]
    elif bool(args.headline_ef_v31_only):
        methods = list(EF_V31_HEADLINE_METHODS)
    else:
        methods = helpers._ef_v31_methods(available_methods)
    methods = list(dict.fromkeys(methods))
    missing_methods = [method for method in methods if method not in available_methods]
    if missing_methods:
        raise ValueError(
            "LCIA method(s) not found for ecoinvent "
            f"{args.ei_version}:\n- " + "\n- ".join(missing_methods)
        )
    if not methods:
        raise ValueError("No EF v3.1 methods were found.")

    if int(args.plot_window_years) < 0:
        raise ValueError("--plot-window-years must be non-negative.")
    year_start = int(args.reference_year) - int(args.plot_window_years)
    year_end = int(args.reference_year) + int(args.plot_window_years)

    print(f"Datapackage: {datapackage}", flush=True)
    print("Excel inventories:", flush=True)
    for path in inventory_paths:
        print(f"  {path}", flush=True)
    print(f"EF v3.1 methods: {len(methods)}", flush=True)
    for method in methods:
        print(f"  {method}", flush=True)

    if bool(args.case_study_activities):
        if args.case_study_activity:
            terminal_activities = [
                DEFAULT_CASE_STUDY_ACTIVITY_KEYS[str(key)]
                for key in args.case_study_activity
            ]
        else:
            terminal_activities = list(DEFAULT_CASE_STUDY_ACTIVITIES)
        print(
            f"Explicit case-study activities: {len(terminal_activities)}",
            flush=True,
        )
    else:
        terminal_activities = helpers._collect_terminal_activities(inventory_paths)
        if not terminal_activities:
            raise ValueError("No terminal activities found in the Excel inventories.")
        print(f"Terminal imported activities: {len(terminal_activities)}", flush=True)
    if args.max_activities is not None:
        terminal_activities = terminal_activities[: int(args.max_activities)]
    for activity in terminal_activities:
        print(
            "  "
            f"{activity.name} | {activity.reference_product} | {activity.location}",
            flush=True,
        )

    output_dir = Path(args.output_dir).expanduser().resolve()
    csv_path = Path(args.results_csv).expanduser().resolve()
    rows: list[dict[str, Any]] = _read_rows(csv_path) if bool(args.resume) else []
    if rows:
        print(
            f"Loaded {len(rows)} existing result row(s) from {csv_path}",
            flush=True,
        )

    contexts: dict[str, Trails] = {}
    activity_maps: dict[str, dict[Any, int]] = {}
    for config in RUN_CONFIGS:
        print(f"\nLoading run context: {config.label}", flush=True)
        contexts[config.key] = _load_trails(
            datapackage=datapackage,
            interpolation_cache_dir=(
                None
                if args.interpolation_cache_dir is None
                else Path(args.interpolation_cache_dir).expanduser().resolve()
            ),
            inventory_paths=inventory_paths,
            import_before_interpolation=bool(args.import_before_interpolation),
            remove_base_temporal_distributions=(
                config.remove_base_temporal_distributions
            ),
            no_cache_interpolation=bool(args.no_cache_interpolation),
            interpolation_start_year_offset=int(args.interpolation_start_year_offset),
            interpolation_end_year_offset=int(args.interpolation_end_year_offset),
        )
        activity_maps[config.key] = helpers._match_activity_indices(
            contexts[config.key],
            terminal_activities,
        )
        missing = [
            activity
            for activity in terminal_activities
            if activity not in activity_maps[config.key]
        ]
        if missing:
            raise ValueError(
                f"Could not match terminal activities after import in "
                f"{config.label}:\n- "
                + "\n- ".join(
                    f"{a.name} | {a.reference_product} | {a.location}" for a in missing
                )
            )

    total_activities = len(terminal_activities)
    try:
        for activity_pos, activity in enumerate(terminal_activities, start=1):
            all_td_activity_index = int(activity_maps["all_td"][activity])
            activity_label = helpers._activity_label(
                contexts["all_td"],
                all_td_activity_index,
                int(args.reference_year),
            )
            if bool(args.resume) and _activity_complete(
                rows,
                output_dir=output_dir,
                activity_index=all_td_activity_index,
                activity_label=activity_label,
                methods=methods,
                depth=int(args.depth),
            ):
                print(
                    "\n"
                    f"[{activity_pos}/{total_activities}] "
                    f"{activity.name} | {activity.reference_product} | "
                    f"{activity.location}",
                    flush=True,
                )
                print(
                    f"  skipping complete activity {all_td_activity_index}: "
                    f"{activity_label}",
                    flush=True,
                )
                continue
            print(
                "\n"
                f"[{activity_pos}/{total_activities}] "
                f"{activity.name} | {activity.reference_product} | "
                f"{activity.location}",
                flush=True,
            )
            routing_min_amount = float(args.routing_min_amount)
            if bool(args.scale_routing_min_amount_by_activity):
                routing_min_amount *= DEFAULT_CASE_STUDY_ROUTING_AMOUNT_SCALES.get(
                    activity,
                    1.0,
                )
            activity_amount = float(args.amount)
            if bool(args.scale_demand_amount_by_activity):
                activity_amount *= DEFAULT_CASE_STUDY_ROUTING_AMOUNT_SCALES.get(
                    activity,
                    1.0,
                )
            results: dict[str, RunResult] = {}
            for config in RUN_CONFIGS:
                activity_index = int(activity_maps[config.key][activity])
                results[config.key] = _run_activity(
                    trails=contexts[config.key],
                    run_config=config,
                    activity_index=activity_index,
                    methods=methods,
                    depth=int(args.depth),
                    reference_year=int(args.reference_year),
                    amount=activity_amount,
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
                    attribute_to_roots=(
                        bool(args.attribute_to_roots)
                        if config.key == "all_td"
                        else bool(args.foreground_attribute_to_roots)
                    ),
                    write_graph_html=bool(args.write_graph_html),
                    graph_output_dir=output_dir,
                    graph_run=str(args.graph_run),
                    graph_min_edge_amount=float(args.graph_min_edge_amount),
                )

            activity_label = results["all_td"].activity_label
            activity_slug = helpers._slugify(
                f"idx_{results['all_td'].activity_index}_{activity_label}",
                max_length=140,
            )
            activity_dir = output_dir / activity_slug
            for method_pos, method in enumerate(methods, start=1):
                figure_path = _expected_figure_path(
                    output_dir=output_dir,
                    activity_index=results["all_td"].activity_index,
                    activity_label=activity_label,
                    method=method,
                    depth=int(args.depth),
                )
                print(
                    f"  plotting method {method_pos}/{len(methods)}: "
                    f"{shorten(method, width=88, placeholder='...')}",
                    flush=True,
                )
                _plot_comparison(
                    all_td_trails=contexts["all_td"],
                    activity_index=results["all_td"].activity_index,
                    activity_label=activity_label,
                    method=method,
                    results=results,
                    reference_year=int(args.reference_year),
                    year_start=year_start,
                    year_end=year_end,
                    output_path=figure_path,
                    width=int(args.width),
                    height=int(args.height),
                )
                rows.append(
                    {
                        "activity_index_all_td": results["all_td"].activity_index,
                        "activity_index_foreground_td_only": results[
                            "foreground_td_only"
                        ].activity_index,
                        "activity": activity_label,
                        "method": method,
                        "static_all_td": results["all_td"].static_scores[method],
                        "static_foreground_td_only": results[
                            "foreground_td_only"
                        ].static_scores[method],
                        "temporal_cumulative_all_td": results["all_td"].total_by_method[
                            method
                        ],
                        "temporal_cumulative_foreground_td_only": results[
                            "foreground_td_only"
                        ].total_by_method[method],
                        "temporal_cumulative_delta": (
                            results["all_td"].total_by_method[method]
                            - results["foreground_td_only"].total_by_method[method]
                        ),
                        "graph_html_all_td": (
                            ""
                            if results["all_td"].graph_html is None
                            else str(results["all_td"].graph_html)
                        ),
                        "graph_html_foreground_td_only": (
                            ""
                            if results["foreground_td_only"].graph_html is None
                            else str(results["foreground_td_only"].graph_html)
                        ),
                        "figure": str(figure_path),
                    }
                )
            _write_rows(csv_path, rows)
            print(
                f"  wrote {len(methods)} comparison figure(s) to {activity_dir}",
                flush=True,
            )
    finally:
        contexts.clear()
        gc.collect()

    _write_rows(csv_path, rows)
    print(f"\nWrote summary CSV: {csv_path}")
    print(f"Wrote figures under: {output_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare terminal imported LCI temporal LCA scores with all TDs "
            "against scores with only foreground LCI TDs."
        )
    )
    parser.add_argument(
        "--datapackage",
        type=Path,
        default=REPO_ROOT / "dev" / "trails_2026-05-18.zip",
    )
    parser.add_argument(
        "--interpolation-cache-dir",
        type=Path,
        default=None,
        help="Optional Trails interpolation cache directory to load instead of opening the datapackage.",
    )
    parser.add_argument(
        "--import-before-interpolation",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Import foreground LCIs on template matrices before annual interpolation.",
    )
    parser.add_argument(
        "--inventories",
        nargs="+",
        type=Path,
        default=[],
        help="Excel inventories to import. Defaults to none unless --all-dev-lci-inventories is set.",
    )
    parser.add_argument(
        "--all-dev-lci-inventories",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Import all dev/lci-*.xlsx inventories.",
    )
    parser.add_argument(
        "--all-onedrive-lci-inventories",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=("Import the four LCI workbooks from the OneDrive trails/data " "folder."),
    )
    parser.add_argument(
        "--case-study-activities",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Process the four explicit case-study foreground activities: "
            "battery electric car, polyol, marine fuel switch, and DACCS."
        ),
    )
    parser.add_argument(
        "--case-study-activity",
        nargs="+",
        choices=tuple(DEFAULT_CASE_STUDY_ACTIVITY_KEYS),
        default=[],
        help=(
            "Restrict --case-study-activities to one or more named cases. "
            "Choices: bev, polyol, marine, daccs."
        ),
    )
    parser.add_argument("--depth", type=int, default=5)
    parser.add_argument("--reference-year", type=int, default=2035)
    parser.add_argument("--plot-window-years", type=int, default=40)
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument("--ei-version", type=str, default="3.12")
    parser.add_argument(
        "--lcia-json",
        type=Path,
        default=Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json"),
    )
    parser.add_argument("--methods", nargs="+", default=[])
    parser.add_argument(
        "--headline-ef-v31-only",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Use only the 16 EF v3.1 headline indicators and skip split "
            "sub-indicators."
        ),
    )
    parser.add_argument(
        "--max-activities",
        type=int,
        default=None,
        help="Limit terminal activities processed, useful for diagnostics.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--results-csv", type=Path, default=DEFAULT_RESULTS_CSV)
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip activities whose rows and PNGs already exist in the output directory.",
    )
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=860)
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
        help=(
            "For explicit case-study activities, multiply --routing-min-amount "
            "by the lifetime system size represented by the foreground LCI. "
            "This keeps the routing cutoff relative to the represented system "
            "without normalizing the inventory amounts."
        ),
    )
    parser.add_argument(
        "--scale-demand-amount-by-activity",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "For explicit case-study activities, multiply --amount by the "
            "lifetime system size represented by the foreground LCI. Use with "
            "lifetime-output production amounts to calculate the full system."
        ),
    )
    parser.add_argument(
        "--attribute-to-roots",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep root-activity attribution during routing/LCA so root areas can be plotted.",
    )
    parser.add_argument(
        "--foreground-attribute-to-roots",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Also compute root attribution for the foreground-only pass. Not "
            "needed for the overlaid cumulative-only foreground line."
        ),
    )
    parser.add_argument(
        "--write-graph-html",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Write interactive routing graph HTML files after temporal routing.",
    )
    parser.add_argument(
        "--graph-run",
        choices=("all_td", "foreground_td_only", "both"),
        default="all_td",
        help=(
            "Which run context to export as routing graph HTML when "
            "--write-graph-html is enabled."
        ),
    )
    parser.add_argument(
        "--graph-min-edge-amount",
        type=float,
        default=1e-9,
        help=(
            "Minimum absolute edge amount included in graph HTML exports. "
            "Use 0 for the full routed graph; large routed graphs can produce "
            "very large HTML files."
        ),
    )
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--no-cache-interpolation", action="store_true")
    parser.add_argument("--interpolation-start-year-offset", type=int, default=-20)
    parser.add_argument("--interpolation-end-year-offset", type=int, default=20)
    return parser


if __name__ == "__main__":
    raise SystemExit(run(build_parser().parse_args()))
