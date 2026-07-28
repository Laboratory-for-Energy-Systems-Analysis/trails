#!/usr/bin/env python3
"""Profile the BrightCon DACCS workflow with a hard RSS safety ceiling.

The parent process launches an isolated worker, samples the resident set size
(RSS) of the worker and all of its descendants, and terminates it if the
configured ceiling is crossed.  The worker follows the DACCS section of the
BrightCon notebook: load the premise package, import the foreground, route the
functional unit, run temporal LCA, then evaluate FaIR and prospective AWARE.

Run this from the BrightCon conference directory, where ``data/`` contains the
notebook inputs::

    /opt/homebrew/Caskroom/miniforge/base/envs/trails/bin/python \
        /Users/romain/GitHub/trails/dev/profile_daccs_memory.py \
        --data-dir data --output results/daccs_memory_profile.json

The requested ceiling defaults to 12 GiB.  Unless ``--allow-memory-pressure``
is passed, the effective ceiling is reduced to 80% of currently available RAM.
This makes the profiler safe to run while other applications are open.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import gc
from importlib.resources import files
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback
from typing import Any, Iterator

import psutil


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

GWP = (
    "IPCC 2021 (incl. biogenic CO2) - climate change: total "
    "(incl. biogenic CO2) - global warming potential (GWP100)"
)
START_YEAR = 2025
DACCS_AMOUNT_KG = 20_000_000_000.0
DACCS_NAME = (
    "carbon dioxide, captured, with a solvent-based direct air capture "
    "system, 1MtCO2"
)
DACCS_PRODUCT = "carbon dioxide, captured"
DACCS_LOCATION = "Europe"


def _jsonable(value: Any) -> Any:
    """Convert diagnostics containing NumPy/path objects to JSON values."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return str(value)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    os.replace(temporary, path)


class EventLog:
    """Append phase boundaries for the supervising parent process."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.started = time.perf_counter()

    def emit(self, event: str, **payload: Any) -> None:
        record = {
            "event": event,
            "elapsed_seconds": time.perf_counter() - self.started,
            **payload,
        }
        with self.path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(_jsonable(record), sort_keys=True) + "\n")

    @contextmanager
    def phase(self, name: str) -> Iterator[None]:
        process = psutil.Process()
        started = time.perf_counter()
        self.emit("phase_start", phase=name, worker_rss=process.memory_info().rss)
        print(f"\n[{name}] started", flush=True)
        try:
            yield
        except Exception as error:
            self.emit(
                "phase_error",
                phase=name,
                seconds=time.perf_counter() - started,
                error=repr(error),
                worker_rss=process.memory_info().rss,
            )
            raise
        else:
            seconds = time.perf_counter() - started
            self.emit(
                "phase_end",
                phase=name,
                seconds=seconds,
                worker_rss=process.memory_info().rss,
            )
            print(f"[{name}] finished in {seconds:.2f} s", flush=True)


def _prospective_aware_method(ssp: str = "SSP126") -> dict[str, Any]:
    path = files("edges").joinpath(
        "data/AWARE 2.0 prospective_Country_all_yearly.json"
    )
    method = json.loads(path.read_text(encoding="utf-8"))
    if ssp not in method["parameters"]:
        raise ValueError(f"Unknown AWARE scenario {ssp!r}")
    method["parameters"] = {ssp: method["parameters"][ssp]}
    method["name"] = (
        f"AWARE 2.0 prospective | Country | all | yearly | {ssp}"
    )
    return method


def _find_daccs_activity(model: Any) -> tuple[int, dict[str, Any]]:
    from trails import search_activity

    matches = []
    table = search_activity(model, "carbon dioxide, captured")
    for row in table.rows:
        item = dict(zip(table.field_names, row))
        if (
            item.get("name") == DACCS_NAME
            and item.get("reference product") == DACCS_PRODUCT
            and item.get("location") == DACCS_LOCATION
        ):
            matches.append(item)
    if len(matches) != 1:
        raise ValueError(
            f"Expected one exact DACCS activity match; found {len(matches)}"
        )
    return int(matches[0]["index"]), matches[0]


def _sum_dataarray(value: Any) -> float:
    reduced = value.sum()
    if hasattr(reduced, "compute"):
        reduced = reduced.compute(scheduler="synchronous")
    if hasattr(reduced, "item"):
        return float(reduced.item())
    return float(reduced)


def _run_worker(args: argparse.Namespace) -> int:
    import dask.array as da
    from datapackage import Package

    from trails import Trails, plot_temp, plot_temporal_scores
    from trails.edges_matrix import score_inventory_with_edges
    from trails.fair_rf import run_fair_delta_rf

    event_log = EventLog(Path(args.events))
    result_path = Path(args.worker_result)
    data_dir = Path(args.data_dir).expanduser().resolve()
    package_path = data_dir / "trails_remind_SSP2-PkBudg1000.zip"
    inventory_path = data_dir / "lci-case-study-daccs_storage_risk.xlsx"
    missing = [path for path in (package_path, inventory_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing DACCS profiler inputs: {missing}")

    event_log.emit("worker_start", pid=os.getpid())
    phase_times: dict[str, float] = {}
    model = None
    report: dict[str, Any] = {
        "inputs": {
            "package": str(package_path),
            "inventory": str(inventory_path),
        },
        "functional_unit": {
            "activity_name": DACCS_NAME,
            "reference_product": DACCS_PRODUCT,
            "location": DACCS_LOCATION,
            "amount_kg": DACCS_AMOUNT_KG,
            "start_year": START_YEAR,
        },
        "configuration": {
            "inventory_backend": "auto",
            "inventory_memory_budget_mib": args.memory_budget_mib,
            "solver_mode": "iterative",
            "iterative_rtol": 1e-3,
            "fair": not args.skip_fair,
            "aware": not args.skip_aware,
            "plots": not args.skip_plots,
        },
    }

    def checkpoint() -> None:
        _atomic_json(result_path, report)
        if args.checkpoint_output:
            _atomic_json(Path(args.checkpoint_output), report)

    @contextmanager
    def measured(name: str) -> Iterator[None]:
        started = time.perf_counter()
        with event_log.phase(name):
            yield
        phase_times[name] = time.perf_counter() - started
        report["phase_seconds"] = phase_times
        report["worker_status"] = "running"
        checkpoint()

    try:
        with measured("load_prospective_background"):
            model = Trails(
                Package(str(package_path)),
                interpolate_annual=True,
                interpolation_start_year_offset=-20,
                interpolation_end_year_offset=20,
                methods=[GWP],
                ei_version="3.12",
            )

        with measured("import_daccs_foreground"):
            model.import_excel_inventory(str(inventory_path))

        with measured("select_functional_unit"):
            daccs_index, activity = _find_daccs_activity(model)
            report["model"] = {
                "activities": int(model.A.shape[1]),
                "biosphere_flows": int(model.B.shape[2]),
                "selected_activity_index": daccs_index,
                "selected_activity": activity,
            }

        with measured("temporal_routing"):
            model.temporal_routing(
                start_year=START_YEAR,
                start_act_idx=daccs_index,
                amount=DACCS_AMOUNT_KG,
                attribute_to_roots=True,
                show_progress=args.show_progress,
            )

        with measured("temporal_lca"):
            model.lca(
                compute_score=True,
                store_inventory=True,
                show_progress=args.show_progress,
                solver_mode="iterative",
                iterative_rtol=1e-3,
                inventory_backend="auto",
                inventory_memory_budget=int(args.memory_budget_mib * 2**20),
                inventory_store=args.inventory_store,
            )

        report["inventory"] = {
            "shape": list(model.inventory.shape),
            "chunks": (
                [list(axis) for axis in model.inventory.data.chunks]
                if isinstance(model.inventory.data, da.Array)
                else None
            ),
            "lazy": isinstance(model.inventory.data, da.Array),
            "characterized_lazy": isinstance(
                model.characterized_inventory.data, da.Array
            ),
            "diagnostics": model.inventory_diagnostics,
            "lca_diagnostics": model.lca_diagnostics,
            "gwp_total": _sum_dataarray(model.scores),
        }
        checkpoint()

        if not args.skip_plots:
            with measured("plot_gwp"):
                figure = plot_temporal_scores(
                    trails=model,
                    stacked=True,
                    legend_top_n=7,
                    show_flow_contributions=False,
                    title="DACCS: life-cycle GWP through time",
                    method_label="kg CO2-eq",
                    year_range=(2020, 2050),
                    year_tick=5,
                    reference_year=START_YEAR,
                    show_cumulative_axis=True,
                    width=850,
                    height=480,
                )
                report["gwp_plot_traces"] = len(figure.data)
                del figure

        if not args.skip_fair:
            with measured("fair_temperature"):
                run_fair_delta_rf(
                    model,
                    scenario="REMIND|SSP2-PkBudg1000",
                    scale_target_fraction=0.1,
                )
            if not args.skip_plots:
                with measured("plot_fair_temperature"):
                    figure = plot_temp(
                        model,
                        by="root activity",
                        title="DACCS: temperature response",
                        method_label="degC",
                        year_range=(2020, 2200),
                        year_tick=20,
                        width=850,
                        height=480,
                    )
                    report["fair_plot_traces"] = len(figure.data)
                    del figure

        if not args.skip_aware:
            with measured("load_prospective_aware"):
                aware_method = _prospective_aware_method("SSP126")
            with measured("edges_aware"):
                aware_scores = score_inventory_with_edges(
                    model,
                    [aware_method],
                    reuse_cached_cfs=True,
                    show_progress=args.show_progress,
                )
            report["aware_total"] = _sum_dataarray(aware_scores)
            if not args.skip_plots:
                with measured("plot_aware"):
                    figure = plot_temporal_scores(
                        trails=model,
                        stacked=True,
                        legend_top_n=7,
                        show_flow_contributions=False,
                        title="DACCS: prospective regionalized AWARE impact",
                        method_label="m3 deprived water-eq.",
                        year_range=(2020, 2050),
                        year_tick=5,
                        reference_year=START_YEAR,
                        show_cumulative_axis=True,
                        width=850,
                        height=480,
                    )
                    report["aware_plot_traces"] = len(figure.data)
                    del figure

        report["phase_seconds"] = phase_times
        report["worker_status"] = "completed"
        checkpoint()
        event_log.emit("worker_complete")
        return 0
    except Exception as error:
        report["phase_seconds"] = phase_times
        report["worker_status"] = "failed"
        report["error"] = repr(error)
        report["traceback"] = traceback.format_exc()
        checkpoint()
        event_log.emit("worker_failed", error=repr(error))
        traceback.print_exc()
        return 1
    finally:
        if model is not None:
            model.close()
        gc.collect()


def _process_tree_rss(process: psutil.Process) -> int:
    rss = 0
    members = [process]
    try:
        members.extend(process.children(recursive=True))
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass
    for member in members:
        try:
            rss += int(member.memory_info().rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return rss


def _read_events(path: Path, offset: int) -> tuple[list[dict[str, Any]], int]:
    if not path.exists():
        return [], offset
    events: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as stream:
        stream.seek(offset)
        for line in stream:
            try:
                events.append(json.loads(line))
            except json.JSONDecodeError:
                continue
        return events, stream.tell()


def _run_supervisor(args: argparse.Namespace) -> int:
    data_dir = Path(args.data_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    requested_limit = int(args.rss_limit_gib * 2**30)
    available = int(psutil.virtual_memory().available)
    if args.allow_memory_pressure:
        effective_limit = requested_limit
    else:
        effective_limit = min(requested_limit, max(1 * 2**30, int(available * 0.8)))

    print(f"Available RAM: {available / 2**30:.2f} GiB")
    print(f"Requested RSS ceiling: {requested_limit / 2**30:.2f} GiB")
    print(f"Effective RSS ceiling: {effective_limit / 2**30:.2f} GiB")
    if effective_limit < requested_limit:
        print(
            "The ceiling was reduced to 80% of available RAM. Close other "
            "applications or pass --allow-memory-pressure to use the requested limit."
        )

    with tempfile.TemporaryDirectory(prefix="trails-daccs-profile-") as temporary:
        temporary_path = Path(temporary)
        event_path = temporary_path / "events.jsonl"
        worker_result_path = temporary_path / "worker-result.json"
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "--data-dir",
            str(data_dir),
            "--events",
            str(event_path),
            "--worker-result",
            str(worker_result_path),
            "--checkpoint-output",
            str(output_path),
            "--memory-budget-mib",
            str(args.memory_budget_mib),
        ]
        for enabled, flag in (
            (args.skip_fair, "--skip-fair"),
            (args.skip_aware, "--skip-aware"),
            (args.skip_plots, "--skip-plots"),
            (args.show_progress, "--show-progress"),
        ):
            if enabled:
                command.append(flag)
        inventory_store = (
            Path(args.inventory_store)
            if args.inventory_store
            else temporary_path / "inventory-store"
        )
        command.extend(["--inventory-store", str(inventory_store)])

        started = time.perf_counter()
        child = subprocess.Popen(command, env={**os.environ, "PYTHONUNBUFFERED": "1"})
        process = psutil.Process(child.pid)
        peak_rss = 0
        phase_peaks: dict[str, int] = {}
        active_phase = "startup"
        event_offset = 0
        exceeded = False

        minimum_available = int(args.min_available_gib * 2**30)
        low_available = False
        while child.poll() is None:
            events, event_offset = _read_events(event_path, event_offset)
            for record in events:
                if record.get("event") == "phase_start":
                    active_phase = str(record.get("phase"))
                elif record.get("event") in {"phase_end", "phase_error"}:
                    active_phase = "between_phases"

            rss = _process_tree_rss(process)
            available_now = int(psutil.virtual_memory().available)
            peak_rss = max(peak_rss, rss)
            phase_peaks[active_phase] = max(phase_peaks.get(active_phase, 0), rss)
            unsafe_reason = None
            if rss > effective_limit:
                exceeded = True
                unsafe_reason = (
                    f"worker RSS {rss / 2**30:.2f} GiB exceeded its "
                    f"{effective_limit / 2**30:.2f} GiB ceiling"
                )
            elif available_now < minimum_available:
                low_available = True
                unsafe_reason = (
                    f"system available RAM {available_now / 2**30:.2f} GiB fell "
                    f"below the {minimum_available / 2**30:.2f} GiB reserve"
                )
            if unsafe_reason is not None:
                print(
                    f"\nMemory safety guard triggered during {active_phase}: "
                    f"{unsafe_reason}. Terminating worker.",
                    file=sys.stderr,
                    flush=True,
                )
                child.terminate()
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    child.kill()
                break
            time.sleep(args.sample_interval)

        return_code = child.wait()
        elapsed = time.perf_counter() - started
        worker_report: dict[str, Any] = {}
        if worker_result_path.exists():
            worker_report = json.loads(worker_result_path.read_text(encoding="utf-8"))

        inventory = worker_report.get("inventory", {})
        diagnostics = inventory.get("diagnostics", {})
        checks = {
            "worker_completed": return_code == 0 and not exceeded,
            "rss_below_ceiling": (
                not exceeded and not low_available and peak_rss <= effective_limit
            ),
            "inventory_is_lazy": inventory.get("lazy") is True,
            "characterized_inventory_is_lazy": (
                inventory.get("characterized_lazy") is True
            ),
            "chunked_backend_used": diagnostics.get("backend") == "chunked",
            "inventory_nonempty": diagnostics.get("canonical_entries", 0) > 0,
        }
        report = {
            **worker_report,
            "supervisor": {
                "status": "passed" if all(checks.values()) else "failed",
                "return_code": return_code,
                "wall_seconds": elapsed,
                "sample_interval_seconds": args.sample_interval,
                "available_ram_at_start": available,
                "requested_rss_limit": requested_limit,
                "effective_rss_limit": effective_limit,
                "minimum_available_ram": minimum_available,
                "peak_rss": peak_rss,
                "phase_peak_rss": phase_peaks,
                "rss_limit_exceeded": exceeded,
                "available_ram_reserve_crossed": low_available,
                "checks": checks,
            },
        }
        _atomic_json(output_path, report)

    print(f"\nProfile written to {output_path}")
    print(f"Wall time: {elapsed:.2f} s")
    print(f"Peak RSS: {peak_rss / 2**30:.2f} GiB")
    print(f"Status: {report['supervisor']['status'].upper()}")
    for name, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}  {name}")
    return 0 if all(checks.values()) else 2


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data", help="BrightCon data directory")
    parser.add_argument(
        "--output",
        default="results/daccs_memory_profile.json",
        help="Supervisor JSON report",
    )
    parser.add_argument("--rss-limit-gib", type=float, default=12.0)
    parser.add_argument("--memory-budget-mib", type=int, default=256)
    parser.add_argument("--sample-interval", type=float, default=0.25)
    parser.add_argument(
        "--min-available-gib",
        type=float,
        default=1.5,
        help="Terminate if system-wide available RAM falls below this reserve",
    )
    parser.add_argument("--inventory-store", type=Path)
    parser.add_argument("--allow-memory-pressure", action="store_true")
    parser.add_argument("--show-progress", action="store_true")
    parser.add_argument("--skip-fair", action="store_true")
    parser.add_argument("--skip-aware", action="store_true")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--events", help=argparse.SUPPRESS)
    parser.add_argument("--worker-result", help=argparse.SUPPRESS)
    parser.add_argument("--checkpoint-output", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.worker:
        if not args.events or not args.worker_result:
            raise ValueError("Worker mode requires --events and --worker-result")
        return _run_worker(args)
    return _run_supervisor(args)


if __name__ == "__main__":
    raise SystemExit(main())
