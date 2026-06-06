from __future__ import annotations

import argparse
import cProfile
import csv
import importlib.util
import io
import json
import os
import pstats
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

RUNNER_PATH = REPO_ROOT / "dev" / "plot_terminal_lci_td_comparison.py"
DEFAULT_DATAPACKAGE = REPO_ROOT / "dev" / "trails_remind_SSP2-PkBudg1000.zip"
DEFAULT_LCIA_JSON = Path("/Users/romain/GitHub/pathways/pathways/data/lcia_ei312.json")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "dev" / "profiling" / "temporal_routing_modes"
DEFAULT_CASE_KEY = "polyol"
DEFAULT_REFERENCE_YEAR = 2025
DEFAULT_ROUTING_MIN_AMOUNT = 1e-3
DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF = 1e-4
DEFAULT_ADAPTIVE_MIN_DEPTH = 1
DEFAULT_ADAPTIVE_MAX_DEPTH = 16
DEFAULT_ADAPTIVE_MAX_DEPTHS = (2, 4, 6, 8, 10, 12, 14, 16)

CASE_METHODS_BY_KEY = {
    "polyol": "EF v3.1 - eutrophication: terrestrial - accumulated exceedance (AE)",
    "bev": "EF v3.1 - human toxicity: carcinogenic - comparative toxic unit for human (CTUh)",
    "marine": (
        "EF v3.1 - material resources: metals/minerals - abiotic depletion "
        "potential (ADP): elements (ultimate reserves)"
    ),
    "daccs": "EF v3.1 - particulate matter formation - impact on human health",
}

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


runner = _load_runner()
helpers = runner.helpers


def _case_activity_defs() -> dict[str, Any]:
    case_defs = dict(runner.DEFAULT_CASE_STUDY_ACTIVITY_KEYS)
    case_defs["polyol"] = helpers.ActivityDef(
        "polyol precursor production from captured CO2",
        "polyol precursor",
        "RER",
    )
    return case_defs


def _default_inventory_paths(*, use_local: bool = False) -> list[Path]:
    paths = [
        REPO_ROOT / "dev" / "lci-case-study-ccu_polyol_delayed_release.xlsx",
        REPO_ROOT / "dev" / "lci-case-study-daccs_storage_risk.xlsx",
        REPO_ROOT / "dev" / "lci-case-study-marine_fuel_switch.xlsx",
        REPO_ROOT / "dev" / "lci-pass_cars.xlsx",
    ]
    if use_local and all(path.exists() for path in paths):
        return paths
    return list(runner.DEFAULT_ONEDRIVE_INVENTORY_PATHS)


def _validate_paths(datapackage: Path, inventory_paths: list[Path]) -> None:
    if not datapackage.exists():
        raise FileNotFoundError(datapackage)
    missing = [path for path in inventory_paths if not Path(path).exists()]
    if missing:
        raise FileNotFoundError(
            "Missing inventory file(s):\n- "
            + "\n- ".join(str(path) for path in missing)
        )


def _graph_stats(trails: Any) -> dict[str, Any]:
    graph = getattr(trails, "graph", None)
    if graph is None:
        return {}
    routing_params = getattr(trails, "_routing_params", {}) or {}

    nodes_by_depth: dict[int, int] = {}
    frontier_reasons: dict[str, float] = {}
    frontier_nodes = 0
    direct_bio_nodes = 0
    adaptive_pruned_nodes = 0
    adaptive_pruned_amount_abs = 0.0
    adaptive_cutoff_potentials: list[float] = []
    frontier_score_potential_sum = 0.0

    for _node, data in graph.nodes(data=True):
        depth = int(data.get("depth", 0))
        nodes_by_depth[depth] = nodes_by_depth.get(depth, 0) + 1
        frontier_amount = float(data.get("frontier_amount") or 0.0)
        if frontier_amount:
            frontier_nodes += 1
            frontier_score_potential_sum += float(data.get("score_potential") or 0.0)
            for reason, amount in (data.get("frontier_reasons") or {}).items():
                key = f"frontier_reason_{reason}_amount"
                frontier_reasons[key] = frontier_reasons.get(key, 0.0) + float(amount)
        if data.get("direct_bio_amount"):
            direct_bio_nodes += 1
        if data.get("adaptive_cutoff_reason") == "adaptive_relative_score_cutoff":
            adaptive_pruned_nodes += 1
            adaptive_pruned_amount_abs += abs(frontier_amount)
            adaptive_cutoff_potentials.append(
                float(data.get("adaptive_cutoff_potential") or 0.0)
            )

    out: dict[str, Any] = {
        "graph_nodes": int(graph.number_of_nodes()),
        "graph_edges": int(graph.number_of_edges()),
        "graph_max_depth": int(max(nodes_by_depth, default=0)),
        "routing_nodes_processed": int(routing_params.get("nodes_processed") or 0),
        "routing_max_processed_depth": int(
            routing_params.get("max_processed_depth") or 0
        ),
        "frontier_nodes": int(frontier_nodes),
        "direct_bio_nodes": int(direct_bio_nodes),
        "adaptive_pruned_nodes": int(adaptive_pruned_nodes),
        "adaptive_pruned_amount_abs": float(adaptive_pruned_amount_abs),
        "frontier_score_potential_sum": float(frontier_score_potential_sum),
        "adaptive_cutoff_potential_sum": float(sum(adaptive_cutoff_potentials)),
        "adaptive_cutoff_potential_max": (
            float(max(adaptive_cutoff_potentials))
            if adaptive_cutoff_potentials
            else 0.0
        ),
    }
    for depth, count in sorted(nodes_by_depth.items()):
        out[f"nodes_depth_{depth}"] = int(count)
    out.update(frontier_reasons)
    return out


def _adaptive_depth_threshold_proof(
    trails: Any,
    *,
    depth: int,
    output_dir: Path,
    prefix: str,
    run_stats: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write and summarize adaptive threshold checks for one graph depth."""
    graph = getattr(trails, "graph", None)
    if graph is None:
        return {}

    routing_params = getattr(trails, "_routing_params", {}) or {}
    relative_cutoff = routing_params.get("adaptive_relative_score_cutoff")
    root_potential = routing_params.get("adaptive_root_score_potential")
    effective_cutoff = routing_params.get("adaptive_effective_score_cutoff")
    if effective_cutoff is None and relative_cutoff is not None and root_potential:
        effective_cutoff = abs(float(relative_cutoff)) * float(root_potential)
    if effective_cutoff is None:
        return {}

    effective = float(effective_cutoff)
    static_score = None
    static_score_threshold = None
    if run_stats and "static_score" in run_stats and relative_cutoff is not None:
        static_score = float(run_stats["static_score"])
        static_score_threshold = abs(float(relative_cutoff)) * abs(static_score)
    proof_depth = int(depth)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    expanded_violations: list[dict[str, Any]] = []
    processed_violations: list[dict[str, Any]] = []
    static_expanded_violations: list[dict[str, Any]] = []
    static_processed_violations: list[dict[str, Any]] = []
    processed_ratios: list[float] = []
    expanded_ratios: list[float] = []
    static_processed_ratios: list[float] = []
    static_expanded_ratios: list[float] = []
    adaptive_frontier = 0
    min_amount_frontier = 0

    for node, data in graph.nodes(data=True):
        node_depth = int(data.get("depth", 0))
        if node_depth != proof_depth:
            continue

        score_potential = float(data.get("score_potential") or 0.0)
        ratio = float("inf") if effective == 0.0 else score_potential / effective
        static_ratio = (
            None
            if static_score_threshold is None
            else (
                float("inf")
                if static_score_threshold == 0.0
                else score_potential / float(static_score_threshold)
            )
        )
        amount = float(data.get("amount") or 0.0)
        frontier_amount = float(data.get("frontier_amount") or 0.0)
        frontier_reasons = data.get("frontier_reasons") or {}
        adaptive_reason = data.get("adaptive_cutoff_reason")
        out_degree = int(graph.out_degree(node))
        expanded = out_degree > 0
        processed = amount != 0.0 and adaptive_reason is None

        if adaptive_reason == "adaptive_relative_score_cutoff":
            adaptive_frontier += 1
        if "min_amount" in frontier_reasons:
            min_amount_frontier += 1
        if processed:
            processed_ratios.append(ratio)
            if static_ratio is not None:
                static_processed_ratios.append(float(static_ratio))
        if expanded:
            expanded_ratios.append(ratio)
            if static_ratio is not None:
                static_expanded_ratios.append(float(static_ratio))

        row = {
            "node": repr(node),
            "year": int(data.get("year", 0)),
            "depth": node_depth,
            "act_idx": int(data.get("act_idx", -1)),
            "name": str(data.get("name") or ""),
            "reference_product": str(data.get("reference_product") or ""),
            "location": str(data.get("location") or ""),
            "amount": amount,
            "frontier_amount": frontier_amount,
            "score_potential": score_potential,
            "threshold": effective,
            "ratio_to_threshold": ratio,
            "static_score_threshold": (
                "" if static_score_threshold is None else float(static_score_threshold)
            ),
            "ratio_to_static_score_threshold": (
                "" if static_ratio is None else float(static_ratio)
            ),
            "out_degree": out_degree,
            "processed": bool(processed),
            "expanded": bool(expanded),
            "adaptive_cutoff_reason": (
                "" if adaptive_reason is None else str(adaptive_reason)
            ),
            "frontier_reasons": json.dumps(frontier_reasons, sort_keys=True),
        }
        rows.append(row)
        if expanded and score_potential <= effective:
            expanded_violations.append(row)
        if processed and score_potential <= effective:
            processed_violations.append(row)
        if static_score_threshold is not None:
            if expanded and score_potential <= float(static_score_threshold):
                static_expanded_violations.append(row)
            if processed and score_potential <= float(static_score_threshold):
                static_processed_violations.append(row)

    rows.sort(key=lambda row: float(row["ratio_to_threshold"]))
    proof_path = (
        output_dir / f"{prefix}_depth{proof_depth}_adaptive_threshold_proof.csv"
    )
    violation_path = (
        output_dir / f"{prefix}_depth{proof_depth}_adaptive_threshold_violations.csv"
    )
    static_violation_path = (
        output_dir
        / f"{prefix}_depth{proof_depth}_static_score_threshold_violations.csv"
    )
    fieldnames = [
        "node",
        "year",
        "depth",
        "act_idx",
        "name",
        "reference_product",
        "location",
        "amount",
        "frontier_amount",
        "score_potential",
        "threshold",
        "ratio_to_threshold",
        "static_score_threshold",
        "ratio_to_static_score_threshold",
        "out_degree",
        "processed",
        "expanded",
        "adaptive_cutoff_reason",
        "frontier_reasons",
    ]

    def _min_or_nan(values: list[float]) -> float:
        return float(min(values)) if values else float("nan")

    summary: dict[str, Any] = {
        "proof_depth": proof_depth,
        "proof_threshold": effective,
        "proof_relative_cutoff": (
            None if relative_cutoff is None else float(relative_cutoff)
        ),
        "proof_root_score_potential": (
            None if root_potential is None else float(root_potential)
        ),
        "proof_static_score": static_score,
        "proof_static_score_threshold": static_score_threshold,
        "proof_static_score_threshold_to_algorithm_threshold_ratio": (
            None
            if static_score_threshold is None or effective == 0.0
            else float(static_score_threshold) / effective
        ),
        "proof_root_potential_to_abs_static_score_ratio": (
            None
            if static_score is None or static_score == 0.0
            else float(root_potential) / abs(float(static_score))
        ),
        "proof_depth_nodes": int(len(rows)),
        "proof_depth_processed_nodes": int(len(processed_ratios)),
        "proof_depth_expanded_nodes": int(len(expanded_ratios)),
        "proof_depth_adaptive_frontier_nodes": int(adaptive_frontier),
        "proof_depth_min_amount_frontier_nodes": int(min_amount_frontier),
        "proof_depth_processed_min_ratio": _min_or_nan(processed_ratios),
        "proof_depth_expanded_min_ratio": _min_or_nan(expanded_ratios),
        "proof_depth_processed_min_ratio_to_static_score_threshold": _min_or_nan(
            static_processed_ratios
        ),
        "proof_depth_expanded_min_ratio_to_static_score_threshold": _min_or_nan(
            static_expanded_ratios
        ),
        "proof_depth_processed_violations": int(len(processed_violations)),
        "proof_depth_expanded_violations": int(len(expanded_violations)),
        "proof_depth_processed_static_score_threshold_violations": int(
            len(static_processed_violations)
        ),
        "proof_depth_expanded_static_score_threshold_violations": int(
            len(static_expanded_violations)
        ),
        "proof_depth_csv": str(proof_path),
        "proof_depth_violations_csv": str(violation_path),
        "proof_depth_static_score_threshold_violations_csv": str(static_violation_path),
    }
    if run_stats:
        for key in (
            "routing_mode",
            "case_key",
            "activity_index",
            "activity",
            "method",
            "amount",
            "reference_year",
            "load_seconds",
            "static_lca_seconds",
            "static_score",
            "routing_seconds",
            "routing_nodes_processed",
            "routing_max_processed_depth",
            "graph_nodes",
            "graph_edges",
            "graph_max_depth",
            "frontier_nodes",
            "direct_bio_nodes",
            "adaptive_pruned_nodes",
            "adaptive_pruned_amount_abs",
            "frontier_score_potential_sum",
            "adaptive_cutoff_potential_sum",
            "adaptive_cutoff_potential_max",
            "temporal_lca_seconds",
            "temporal_score",
            "score_delta_static",
            "relative_delta_static",
        ):
            if key in run_stats:
                summary[key] = run_stats[key]
        for key, value in run_stats.items():
            if str(key).startswith("nodes_depth_"):
                summary[key] = value

    summary_path = (
        output_dir / f"{prefix}_depth{proof_depth}_adaptive_threshold_proof.json"
    )
    summary["proof_depth_summary_json"] = str(summary_path)
    _write_csv(rows, proof_path, fieldnames=fieldnames)
    _write_csv(
        processed_violations + expanded_violations,
        violation_path,
        fieldnames=fieldnames,
    )
    _write_csv(
        static_processed_violations + static_expanded_violations,
        static_violation_path,
        fieldnames=fieldnames,
    )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote proof summary JSON: {summary_path}", flush=True)
    print(f"Wrote proof node CSV: {proof_path}", flush=True)
    print(f"Wrote proof violation CSV: {violation_path}", flush=True)
    print(
        "Wrote static-score threshold violation CSV: " f"{static_violation_path}",
        flush=True,
    )
    return summary


def _current_rss_mb() -> float | None:
    """Return current process RSS in MiB when available."""
    try:
        import psutil  # type: ignore[import-not-found]

        return float(psutil.Process().memory_info().rss) / 1024.0 / 1024.0
    except Exception:
        pass

    try:
        output = subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(os.getpid())],
            text=True,
        )
        value = output.strip()
        if value:
            return float(value) / 1024.0
    except Exception:
        return None
    return None


def _process_peak_rss_mb() -> float | None:
    """Return max RSS for the process in MiB when the platform exposes it."""
    try:
        import resource

        value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except Exception:
        return None

    if sys.platform == "darwin":
        return value / 1024.0 / 1024.0
    return value / 1024.0


def _round_optional(value: float | None) -> float | None:
    return None if value is None else round(float(value), 3)


class _MemoryMonitor:
    """Sample process RSS while an expensive step is profiled."""

    def __init__(self, *, interval: float = 2.0) -> None:
        self.interval = float(interval)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._samples: list[float] = []
        self.start_rss_mb: float | None = None
        self.end_rss_mb: float | None = None
        self.process_peak_rss_mb: float | None = None

    def __enter__(self) -> "_MemoryMonitor":
        self.start_rss_mb = _current_rss_mb()
        if self.start_rss_mb is not None:
            self._samples.append(float(self.start_rss_mb))
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(0.1, self.interval))
        self.end_rss_mb = _current_rss_mb()
        if self.end_rss_mb is not None:
            self._samples.append(float(self.end_rss_mb))
        self.process_peak_rss_mb = _process_peak_rss_mb()

    def _sample_loop(self) -> None:
        while not self._stop.wait(self.interval):
            value = _current_rss_mb()
            if value is not None:
                self._samples.append(float(value))

    def as_dict(self) -> dict[str, float | None]:
        sampled_peak = max(self._samples) if self._samples else None
        delta = (
            None
            if self.start_rss_mb is None or self.end_rss_mb is None
            else self.end_rss_mb - self.start_rss_mb
        )
        return {
            "rss_start_mb": _round_optional(self.start_rss_mb),
            "rss_end_mb": _round_optional(self.end_rss_mb),
            "rss_delta_mb": _round_optional(delta),
            "rss_peak_sampled_mb": _round_optional(sampled_peak),
            "rss_peak_process_mb": _round_optional(self.process_peak_rss_mb),
        }


def _prefix_metrics(prefix: str, values: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def _format_optional_mb(value: float | None) -> str:
    return "n/a" if value is None else f"{float(value):.1f} MiB"


def _profile_step(
    label: str,
    output_dir: Path,
    func: Callable[[], Any],
    *,
    sort_by: str,
    limit: int,
) -> tuple[Any, float, Path, Path, dict[str, float | None]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / f"{label}.prof"
    text_path = output_dir / f"{label}_cprofile.txt"

    profiler = cProfile.Profile()
    start = time.perf_counter()
    with _MemoryMonitor() as memory_monitor:
        profiler.enable()
        try:
            result = func()
        finally:
            profiler.disable()
    seconds = time.perf_counter() - start
    memory = memory_monitor.as_dict()

    profiler.dump_stats(str(profile_path))
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats(sort_by)
    stats.print_stats(limit)
    text_path.write_text(stream.getvalue(), encoding="utf-8")

    print(f"{label}: {seconds:.1f}s", flush=True)
    print(
        "  rss: "
        f"start={_format_optional_mb(memory['rss_start_mb'])}, "
        f"peak_sampled={_format_optional_mb(memory['rss_peak_sampled_mb'])}, "
        f"end={_format_optional_mb(memory['rss_end_mb'])}, "
        f"process_peak={_format_optional_mb(memory['rss_peak_process_mb'])}",
        flush=True,
    )
    print(f"  profile: {profile_path}", flush=True)
    print(f"  summary: {text_path}", flush=True)
    return result, seconds, profile_path, text_path, memory


def _parse_depth_list(value: str) -> list[int]:
    depths = [int(item.strip()) for item in str(value).split(",") if item.strip()]
    if not depths:
        raise argparse.ArgumentTypeError("At least one depth is required.")
    if any(depth < 0 for depth in depths):
        raise argparse.ArgumentTypeError("Depth values must be non-negative.")
    return depths


def _score_to_float(score: object) -> float:
    arr = np.asarray(score, dtype=float).ravel()
    if arr.size != 1:
        raise ValueError(f"Expected one score, got {arr.size}.")
    return float(arr[0])


def _temporal_total_score(trails: Any, method: str) -> float:
    if getattr(trails, "scores", None) is None:
        raise RuntimeError("trails.scores is missing after temporal LCA.")
    data = trails.scores
    if "method" in data.dims:
        methods = [str(value) for value in data.coords["method"].values.tolist()]
        data = data.isel(method=methods.index(str(method)), drop=True)
    if data.dims:
        data = data.sum(dim=list(data.dims))
    values = data.data
    if hasattr(values, "todense"):
        return float(np.asarray(values.todense(), dtype=float).sum())
    return float(np.asarray(data.values, dtype=float).sum())


def _run_lca(
    trails: Any,
    *,
    method: str,
    solver_mode: str,
    fallback_solver_mode: str,
    iterative_rtol: float,
    iterative_atol: float,
    iterative_restart: int | None,
    iterative_maxiter: int | None,
    iterative_use_guess: bool,
    iterative_preconditioner: str,
    iterative_ilu_drop_tol: float,
    iterative_ilu_fill_factor: float,
    ei_version: str,
) -> str:
    def call(mode: str) -> None:
        trails.lca(
            methods=[method],
            show_progress=False,
            attribute_to_roots=True,
            compute_score=True,
            store_inventory=False,
            solver_mode=str(mode),
            iterative_rtol=float(iterative_rtol),
            iterative_atol=float(iterative_atol),
            iterative_restart=iterative_restart,
            iterative_maxiter=iterative_maxiter,
            iterative_use_guess=bool(iterative_use_guess),
            iterative_preconditioner=str(iterative_preconditioner),
            iterative_ilu_drop_tol=float(iterative_ilu_drop_tol),
            iterative_ilu_fill_factor=float(iterative_ilu_fill_factor),
            ei_version=str(ei_version),
        )

    try:
        call(str(solver_mode))
        return str(solver_mode)
    except RuntimeError as exc:
        fallback = str(fallback_solver_mode).strip().lower()
        if fallback in {"", "none"} or fallback == str(solver_mode).strip().lower():
            raise
        print(
            f"lca solver_mode={solver_mode!r} failed: {exc}; "
            f"retrying with {fallback_solver_mode!r}",
            flush=True,
        )
        call(str(fallback_solver_mode))
        return str(fallback_solver_mode)


def _write_csv(
    rows: list[dict[str, Any]],
    path: Path,
    *,
    fieldnames: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def _write_outputs(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "summary.json"
    csv_path = output_dir / "summary.csv"
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote summary JSON: {json_path}", flush=True)
    print(f"Wrote summary CSV: {csv_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Profile fixed-depth and adaptive Trails routing modes."
    )
    parser.add_argument(
        "--routing-mode",
        choices=("depth1", "adaptive", "adaptive-sweep"),
        required=True,
    )
    parser.add_argument("--case-key", default=DEFAULT_CASE_KEY)
    parser.add_argument("--datapackage", type=Path, default=DEFAULT_DATAPACKAGE)
    parser.add_argument("--lcia-json", type=Path, default=DEFAULT_LCIA_JSON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--use-local-inventories",
        action="store_true",
        help=(
            "Use the local dev/*.xlsx copies instead of the OneDrive inventories. "
            "The default matches the manuscript notebook."
        ),
    )
    parser.add_argument("--reference-year", type=int, default=DEFAULT_REFERENCE_YEAR)
    parser.add_argument("--amount", type=float, default=None)
    parser.add_argument("--method", default=None)
    parser.add_argument("--ei-version", default="3.12")
    parser.add_argument(
        "--routing-min-amount", type=float, default=DEFAULT_ROUTING_MIN_AMOUNT
    )
    parser.add_argument(
        "--adaptive-relative-score-cutoff",
        type=float,
        default=DEFAULT_ADAPTIVE_RELATIVE_SCORE_CUTOFF,
    )
    parser.add_argument(
        "--adaptive-min-depth", type=int, default=DEFAULT_ADAPTIVE_MIN_DEPTH
    )
    parser.add_argument(
        "--adaptive-max-depth", type=int, default=DEFAULT_ADAPTIVE_MAX_DEPTH
    )
    parser.add_argument(
        "--adaptive-max-depths",
        type=_parse_depth_list,
        default=list(DEFAULT_ADAPTIVE_MAX_DEPTHS),
        help="Comma-separated adaptive max_depth values for adaptive-sweep.",
    )
    parser.add_argument("--no-adaptive-max-depth", action="store_true")
    parser.add_argument(
        "--stop-after-routing-seconds",
        type=float,
        default=600.0,
        help=(
            "In adaptive-sweep mode, stop after a cap whose routing time "
            "exceeds this value."
        ),
    )
    parser.add_argument(
        "--skip-temporal-lca",
        action="store_true",
        help=(
            "Profile static_lca and temporal_routing, but skip the temporal "
            "LCA solve."
        ),
    )
    parser.add_argument(
        "--run-temporal-lca-in-sweep",
        action="store_true",
        help="In adaptive-sweep mode, run temporal LCA after each routed cap.",
    )
    parser.add_argument(
        "--proof-depth",
        type=int,
        default=None,
        help=(
            "Write adaptive score-threshold proof diagnostics for this graph "
            "depth. Expanded and processed nodes at this depth must exceed the "
            "effective adaptive threshold."
        ),
    )
    parser.add_argument("--sort-by", default="cumulative")
    parser.add_argument("--profile-limit", type=int, default=80)
    parser.add_argument("--solver-mode", default="iterative")
    parser.add_argument("--fallback-solver-mode", default="direct")
    parser.add_argument("--iterative-rtol", type=float, default=1e-3)
    parser.add_argument("--iterative-atol", type=float, default=0.0)
    parser.add_argument("--iterative-restart", type=int, default=100)
    parser.add_argument("--iterative-maxiter", type=int, default=1000)
    parser.add_argument(
        "--iterative-use-guess", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--iterative-preconditioner", default="jacobi")
    parser.add_argument("--iterative-ilu-drop-tol", type=float, default=1e-4)
    parser.add_argument("--iterative-ilu-fill-factor", type=float, default=10.0)
    args = parser.parse_args()

    case_defs = _case_activity_defs()
    if args.case_key not in case_defs:
        raise ValueError(f"Unknown case key {args.case_key!r}: {sorted(case_defs)}")
    activity = case_defs[str(args.case_key)]
    method = str(args.method or CASE_METHODS_BY_KEY[str(args.case_key)])
    amount = float(
        args.amount
        if args.amount is not None
        else CASE_DEMAND_AMOUNTS_BY_KEY[str(args.case_key)]
    )
    mode_slug = str(args.routing_mode)
    output_dir = Path(args.output_dir) / f"{args.case_key}_{mode_slug}"

    if args.lcia_json is not None:
        os.environ["TRAILS_LCIA_EI312_JSON"] = str(
            args.lcia_json.expanduser().resolve()
        )

    inventory_paths = _default_inventory_paths(
        use_local=bool(args.use_local_inventories)
    )
    _validate_paths(Path(args.datapackage).expanduser().resolve(), inventory_paths)

    print(f"Python: {sys.executable}", flush=True)
    print(f"Routing mode: {args.routing_mode}", flush=True)
    print(f"Case: {args.case_key}", flush=True)
    print(f"Method: {method}", flush=True)
    print(f"Amount: {amount:g}", flush=True)
    print(f"Datapackage: {args.datapackage}", flush=True)
    print("Inventories:", flush=True)
    for path in inventory_paths:
        print(f"  - {path}", flush=True)

    load_start = time.perf_counter()
    with _MemoryMonitor() as load_memory_monitor:
        trails = runner._load_trails(
            datapackage=Path(args.datapackage).expanduser().resolve(),
            interpolation_cache_dir=None,
            inventory_paths=inventory_paths,
            import_before_interpolation=False,
            remove_base_temporal_distributions=False,
            no_cache_interpolation=False,
            interpolation_start_year_offset=-20,
            interpolation_end_year_offset=20,
        )
    load_seconds = time.perf_counter() - load_start
    load_memory = load_memory_monitor.as_dict()
    trails.methods = [method]
    trails.default_methods = [method]
    trails.ei_version = str(args.ei_version)
    trails.default_ei_version = str(args.ei_version)
    print(f"load_trails_context: {load_seconds:.1f}s", flush=True)
    print(
        "  rss: "
        f"start={_format_optional_mb(load_memory['rss_start_mb'])}, "
        f"peak_sampled={_format_optional_mb(load_memory['rss_peak_sampled_mb'])}, "
        f"end={_format_optional_mb(load_memory['rss_end_mb'])}, "
        f"process_peak={_format_optional_mb(load_memory['rss_peak_process_mb'])}",
        flush=True,
    )

    activity_maps = helpers._match_activity_indices(trails, [activity])
    if activity not in activity_maps:
        raise ValueError(
            "Could not match activity: "
            f"{activity.name} | {activity.reference_product} | {activity.location}"
        )
    activity_index = int(activity_maps[activity])
    activity_label = helpers._activity_label(
        trails,
        activity_index,
        int(args.reference_year),
    )
    print(f"Matched activity: {activity_label} (idx={activity_index})", flush=True)

    rows: list[dict[str, Any]] = []
    row: dict[str, Any] = {
        "routing_mode": args.routing_mode,
        "case_key": args.case_key,
        "activity_index": activity_index,
        "activity": activity_label,
        "method": method,
        "amount": amount,
        "reference_year": int(args.reference_year),
        "load_seconds": load_seconds,
    }
    row.update(_prefix_metrics("load", load_memory))

    def run_static() -> None:
        trails.static_lca(
            year=int(args.reference_year),
            act_idx=activity_index,
            amount=amount,
            methods=[method],
            ei_version=str(args.ei_version),
        )

    _, static_seconds, _, static_text, static_memory = _profile_step(
        f"{args.case_key}_{mode_slug}_static_lca",
        output_dir,
        run_static,
        sort_by=str(args.sort_by),
        limit=int(args.profile_limit),
    )
    static_score = _score_to_float(trails.static_score)
    row.update(
        {
            "static_lca_seconds": static_seconds,
            "static_score": static_score,
            "static_profile_text": str(static_text),
        }
    )
    row.update(_prefix_metrics("static_lca", static_memory))

    def add_depth_proof(current_row: dict[str, Any], *, prefix: str) -> None:
        if args.proof_depth is None:
            return
        proof_stats = _adaptive_depth_threshold_proof(
            trails,
            depth=int(args.proof_depth),
            output_dir=output_dir,
            prefix=prefix,
            run_stats=current_row,
        )
        current_row.update(proof_stats)
        root_score_potential = float(proof_stats.get("proof_root_score_potential"))
        static_threshold = float(proof_stats.get("proof_static_score_threshold"))
        processed_static_violations = proof_stats.get(
            "proof_depth_processed_static_score_threshold_violations"
        )
        expanded_static_violations = proof_stats.get(
            "proof_depth_expanded_static_score_threshold_violations"
        )
        print(
            "  proof result: "
            f"depth={proof_stats.get('proof_depth')}, "
            f"static_score={float(proof_stats.get('proof_static_score')):.6g}, "
            f"root_potential={root_score_potential:.6g}, "
            f"threshold={float(proof_stats.get('proof_threshold')):.6g}, "
            f"static_score_threshold={static_threshold:.6g}, "
            f"nodes={proof_stats.get('proof_depth_nodes')}, "
            f"processed={proof_stats.get('proof_depth_processed_nodes')}, "
            f"expanded={proof_stats.get('proof_depth_expanded_nodes')}, "
            "processed_min_ratio="
            f"{float(proof_stats.get('proof_depth_processed_min_ratio')):.6g}, "
            "expanded_min_ratio="
            f"{float(proof_stats.get('proof_depth_expanded_min_ratio')):.6g}, "
            "processed_violations="
            f"{proof_stats.get('proof_depth_processed_violations')}, "
            "expanded_violations="
            f"{proof_stats.get('proof_depth_expanded_violations')}, "
            "processed_static_violations="
            f"{processed_static_violations}, "
            "expanded_static_violations="
            f"{expanded_static_violations}",
            flush=True,
        )
        if (
            int(proof_stats.get("proof_depth_processed_violations") or 0) > 0
            or int(proof_stats.get("proof_depth_expanded_violations") or 0) > 0
            or int(
                proof_stats.get(
                    "proof_depth_processed_static_score_threshold_violations"
                )
                or 0
            )
            > 0
            or int(
                proof_stats.get(
                    "proof_depth_expanded_static_score_threshold_violations"
                )
                or 0
            )
            > 0
        ):
            static_violation_csv = proof_stats.get(
                "proof_depth_static_score_threshold_violations_csv"
            )
            raise AssertionError(
                "Adaptive threshold proof failed; see "
                f"{proof_stats.get('proof_depth_violations_csv')}"
                " and "
                f"{static_violation_csv}"
            )

    def print_run_summary(current_row: dict[str, Any]) -> None:
        print("\nSummary:", flush=True)
        for key in (
            "load_seconds",
            "static_lca_seconds",
            "routing_seconds",
            "temporal_lca_seconds",
            "routing_nodes_processed",
            "routing_max_processed_depth",
            "graph_nodes",
            "graph_edges",
            "graph_max_depth",
            "frontier_nodes",
            "adaptive_pruned_nodes",
            "static_score",
            "temporal_score",
            "relative_delta_static",
        ):
            if key in current_row:
                print(f"  {key}: {current_row[key]}", flush=True)

    if args.routing_mode == "adaptive-sweep":
        for depth in list(args.adaptive_max_depths):
            routing_kwargs = {
                "start_year": int(args.reference_year),
                "start_act_idx": activity_index,
                "amount": amount,
                "max_depth": int(depth),
                "min_amount": float(args.routing_min_amount),
                "show_progress": False,
                "attribute_to_roots": True,
                "adaptive_methods": [method],
                "adaptive_relative_score_cutoff": float(
                    args.adaptive_relative_score_cutoff
                ),
                "adaptive_ei_version": str(args.ei_version),
                "adaptive_min_depth": int(args.adaptive_min_depth),
            }
            print(
                "\nAdaptive sweep cap "
                f"max_depth={depth}, cutoff={args.adaptive_relative_score_cutoff:g}",
                flush=True,
            )

            def run_sweep_routing() -> None:
                trails.temporal_routing(**routing_kwargs)

            _, routing_seconds, _, routing_text, routing_memory = _profile_step(
                f"{args.case_key}_adaptive_maxdepth_{depth}_temporal_routing",
                output_dir,
                run_sweep_routing,
                sort_by=str(args.sort_by),
                limit=int(args.profile_limit),
            )
            depth_row = dict(row)
            depth_row.update(
                {
                    "routing_mode": "adaptive",
                    "adaptive_sweep": True,
                    "adaptive_max_depth": int(depth),
                    "routing_seconds": routing_seconds,
                    "routing_profile_text": str(routing_text),
                    "routing_kwargs": json.dumps(routing_kwargs, default=str),
                }
            )
            depth_row.update(_prefix_metrics("routing", routing_memory))
            depth_row.update(_graph_stats(trails))

            if bool(args.run_temporal_lca_in_sweep) and not bool(
                args.skip_temporal_lca
            ):

                def run_sweep_temporal_lca() -> str:
                    return _run_lca(
                        trails,
                        method=method,
                        solver_mode=str(args.solver_mode),
                        fallback_solver_mode=str(args.fallback_solver_mode),
                        iterative_rtol=float(args.iterative_rtol),
                        iterative_atol=float(args.iterative_atol),
                        iterative_restart=int(args.iterative_restart),
                        iterative_maxiter=int(args.iterative_maxiter),
                        iterative_use_guess=bool(args.iterative_use_guess),
                        iterative_preconditioner=str(args.iterative_preconditioner),
                        iterative_ilu_drop_tol=float(args.iterative_ilu_drop_tol),
                        iterative_ilu_fill_factor=float(args.iterative_ilu_fill_factor),
                        ei_version=str(args.ei_version),
                    )

                actual_solver_mode, lca_seconds, _, lca_text, lca_memory = (
                    _profile_step(
                        f"{args.case_key}_adaptive_maxdepth_{depth}_temporal_lca",
                        output_dir,
                        run_sweep_temporal_lca,
                        sort_by=str(args.sort_by),
                        limit=int(args.profile_limit),
                    )
                )
                temporal_score = _temporal_total_score(trails, method)
                depth_row.update(
                    {
                        "temporal_lca_seconds": lca_seconds,
                        "actual_solver_mode": actual_solver_mode,
                        "temporal_score": temporal_score,
                        "score_delta_static": temporal_score - static_score,
                        "relative_delta_static": (
                            float("nan")
                            if static_score == 0.0
                            else (temporal_score - static_score) / static_score
                        ),
                        "temporal_lca_profile_text": str(lca_text),
                    }
                )
                depth_row.update(_prefix_metrics("temporal_lca", lca_memory))

            add_depth_proof(
                depth_row,
                prefix=f"{args.case_key}_adaptive_maxdepth_{depth}",
            )
            rows.append(depth_row)
            _write_outputs(rows, output_dir)
            print(
                "  sweep result: "
                f"depth={depth}, routing={routing_seconds:.1f}s, "
                f"visited={depth_row.get('routing_nodes_processed')}, "
                "deepest_visited="
                f"{depth_row.get('routing_max_processed_depth')}, "
                f"nodes={depth_row.get('graph_nodes')}, "
                f"edges={depth_row.get('graph_edges')}, "
                f"max_depth={depth_row.get('graph_max_depth')}, "
                f"pruned={depth_row.get('adaptive_pruned_nodes')}, "
                f"static_score={static_score:.6g}, "
                f"score={depth_row.get('temporal_score')}",
                flush=True,
            )
            if routing_seconds > float(args.stop_after_routing_seconds):
                print(
                    "Stopping adaptive sweep because routing time exceeded "
                    f"{args.stop_after_routing_seconds:g}s.",
                    flush=True,
                )
                break
        return

    if args.routing_mode == "depth1":
        routing_kwargs: dict[str, Any] = {
            "start_year": int(args.reference_year),
            "start_act_idx": activity_index,
            "amount": amount,
            "max_depth": 1,
            "min_amount": float(args.routing_min_amount),
            "show_progress": False,
            "attribute_to_roots": True,
        }
    else:
        max_depth = (
            None if bool(args.no_adaptive_max_depth) else int(args.adaptive_max_depth)
        )
        routing_kwargs = {
            "start_year": int(args.reference_year),
            "start_act_idx": activity_index,
            "amount": amount,
            "max_depth": max_depth,
            "min_amount": float(args.routing_min_amount),
            "show_progress": False,
            "attribute_to_roots": True,
            "adaptive_methods": [method],
            "adaptive_relative_score_cutoff": float(
                args.adaptive_relative_score_cutoff
            ),
            "adaptive_ei_version": str(args.ei_version),
            "adaptive_min_depth": int(args.adaptive_min_depth),
        }

    def run_routing() -> None:
        trails.temporal_routing(**routing_kwargs)

    _, routing_seconds, _, routing_text, routing_memory = _profile_step(
        f"{args.case_key}_{mode_slug}_temporal_routing",
        output_dir,
        run_routing,
        sort_by=str(args.sort_by),
        limit=int(args.profile_limit),
    )
    row.update(
        {
            "routing_seconds": routing_seconds,
            "routing_profile_text": str(routing_text),
            "routing_kwargs": json.dumps(routing_kwargs, default=str),
        }
    )
    row.update(_prefix_metrics("routing", routing_memory))
    row.update(_graph_stats(trails))

    if bool(args.skip_temporal_lca):
        add_depth_proof(row, prefix=f"{args.case_key}_{mode_slug}")
        rows.append(row)
        _write_outputs(rows, output_dir)
        print_run_summary(row)
        return

    def run_temporal_lca() -> str:
        return _run_lca(
            trails,
            method=method,
            solver_mode=str(args.solver_mode),
            fallback_solver_mode=str(args.fallback_solver_mode),
            iterative_rtol=float(args.iterative_rtol),
            iterative_atol=float(args.iterative_atol),
            iterative_restart=int(args.iterative_restart),
            iterative_maxiter=int(args.iterative_maxiter),
            iterative_use_guess=bool(args.iterative_use_guess),
            iterative_preconditioner=str(args.iterative_preconditioner),
            iterative_ilu_drop_tol=float(args.iterative_ilu_drop_tol),
            iterative_ilu_fill_factor=float(args.iterative_ilu_fill_factor),
            ei_version=str(args.ei_version),
        )

    actual_solver_mode, lca_seconds, _, lca_text, lca_memory = _profile_step(
        f"{args.case_key}_{mode_slug}_temporal_lca",
        output_dir,
        run_temporal_lca,
        sort_by=str(args.sort_by),
        limit=int(args.profile_limit),
    )
    temporal_score = _temporal_total_score(trails, method)
    row.update(
        {
            "temporal_lca_seconds": lca_seconds,
            "actual_solver_mode": actual_solver_mode,
            "temporal_score": temporal_score,
            "score_delta_static": temporal_score - static_score,
            "relative_delta_static": (
                float("nan")
                if static_score == 0.0
                else (temporal_score - static_score) / static_score
            ),
            "temporal_lca_profile_text": str(lca_text),
        }
    )
    row.update(_prefix_metrics("temporal_lca", lca_memory))
    add_depth_proof(row, prefix=f"{args.case_key}_{mode_slug}")
    rows.append(row)
    _write_outputs(rows, output_dir)

    print_run_summary(row)


if __name__ == "__main__":
    main()
