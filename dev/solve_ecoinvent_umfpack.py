from __future__ import annotations

import argparse
import os
import time
import warnings
from collections.abc import Mapping
from typing import Any

import bw2data as bd
import numpy as np
from scipy import sparse

try:
    from scikits.umfpack import UMFPACK_A, UmfpackContext, UmfpackWarning
    import scikits.umfpack as umfpack
except ImportError as exc:  # pragma: no cover - environment dependent
    raise RuntimeError(
        "scikit-umfpack is required for this benchmark. Install it with "
        "`pip install -e .[umfpack]` or from conda-forge, then rerun this script."
    ) from exc

try:
    import psutil
except ImportError:  # pragma: no cover - optional diagnostic dependency
    psutil = None


DEFAULT_PROJECT = "ecoinvent-3.12-cutoff"
DEFAULT_DATABASE = "ecoinvent-3.12-cutoff"
_PROCESS = psutil.Process(os.getpid()) if psutil is not None else None


def _rss_bytes() -> int | None:
    if _PROCESS is None:
        return None
    return int(_PROCESS.memory_info().rss)


def _format_bytes(value: int | float | None) -> str:
    if value is None:
        return "n/a"
    value = float(value)
    sign = "-" if value < 0 else ""
    value = abs(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    unit = units[0]
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            break
        value /= 1024.0
    return f"{sign}{value:.3f} {unit}"


def _csc_storage_bytes(matrix: sparse.csc_matrix) -> dict[str, int]:
    return {
        "data": int(matrix.data.nbytes),
        "indices": int(matrix.indices.nbytes),
        "indptr": int(matrix.indptr.nbytes),
    }


def _umfpack_info_bytes(ctx: UmfpackContext, name: str) -> int | None:
    index = getattr(umfpack, name, None)
    if index is None:
        return None
    unit_index = getattr(umfpack, "UMFPACK_SIZE_OF_UNIT", None)
    if unit_index is None:
        return None
    value = float(ctx.info[index])
    unit_size = float(ctx.info[unit_index])
    if value < 0 or unit_size <= 0:
        return None
    return int(round(value * unit_size))


def _project_names() -> set[str]:
    return {project.name for project in bd.projects}


def _activity_dataset_mapping(database_names: set[str] | None = None) -> dict[int, Any]:
    try:
        from bw2data.backends import ActivityDataset
    except ImportError:
        return {}

    query = ActivityDataset.select(
        ActivityDataset.database,
        ActivityDataset.code,
        ActivityDataset.id,
    )
    if database_names:
        query = query.where(ActivityDataset.database << sorted(database_names))

    return {
        int(activity_id): (database, code)
        for database, code, activity_id in query.tuples()
    }


def _legacy_mapping() -> dict[int, Any]:
    mapping = getattr(bd, "mapping", None)
    if mapping is None or not hasattr(mapping, "items"):
        return {}
    return {int(value): key for key, value in mapping.items()}


def _reverse_mapping(database_names: set[str] | None = None) -> dict[int, Any]:
    return _activity_dataset_mapping(database_names) or _legacy_mapping()


def _with_brightway_keys(
    mapping: Mapping[Any, int], database_names: set[str] | None = None
) -> dict[Any, int]:
    reverse = _reverse_mapping(database_names)
    return {reverse.get(int(key), key): int(value) for key, value in mapping.items()}


def _format_key(key: Any) -> str:
    if isinstance(key, tuple):
        return "(" + ", ".join(repr(part) for part in key) + ")"
    return repr(key)


def _select_demand_row(
    product_dict: Mapping[Any, int],
    *,
    database: str,
    demand_code: str | None,
    demand_row: int | None,
) -> tuple[Any, int]:
    pos_to_product = {int(pos): key for key, pos in product_dict.items()}

    if demand_code is not None:
        key = (database, demand_code)
        if key not in product_dict:
            raise KeyError(
                f"Demand product key {_format_key(key)} is not in the product matrix."
            )
        return key, int(product_dict[key])

    if demand_row is not None:
        if demand_row not in pos_to_product:
            raise ValueError(
                f"Demand row {demand_row} is outside the product matrix "
                f"(0..{len(product_dict) - 1})."
            )
        return pos_to_product[demand_row], int(demand_row)

    candidates = [
        (int(pos), key)
        for key, pos in product_dict.items()
        if isinstance(key, tuple) and key[0] == database
    ]
    if not candidates:
        candidates = [(int(pos), key) for key, pos in product_dict.items()]

    row, key = min(candidates, key=lambda item: item[0])
    return key, row


def _load_technosphere_matrix(
    database_name: str,
) -> tuple[sparse.csc_matrix, dict[Any, int], dict[Any, int], float]:
    db = bd.Database(database_name)

    start = time.perf_counter()
    if hasattr(db, "datapackage"):
        import matrix_utils as mu

        technosphere_mm = mu.MappedMatrix(
            packages=[db.datapackage()],
            matrix="technosphere_matrix",
            use_arrays=False,
            use_distributions=False,
        )
        activity_dict = technosphere_mm.col_mapper.to_dict()
        product_dict = technosphere_mm.row_mapper.to_dict()
        technosphere_matrix = technosphere_mm.matrix
    else:
        from bw2calc.matrices import TechnosphereBiosphereMatrixBuilder

        (
            _bio_params,
            _tech_params,
            _biosphere_dict,
            activity_dict,
            product_dict,
            _biosphere_matrix,
            technosphere_matrix,
        ) = TechnosphereBiosphereMatrixBuilder.build([db.filepath_processed()])
    elapsed = time.perf_counter() - start

    if technosphere_matrix.shape[0] != technosphere_matrix.shape[1]:
        raise ValueError(
            "Technosphere matrix is not square: "
            f"{technosphere_matrix.shape[0]} products x "
            f"{technosphere_matrix.shape[1]} activities."
        )

    matrix = technosphere_matrix.tocsc()
    if matrix.dtype != np.float64:
        matrix = matrix.astype(np.float64)
    matrix.sort_indices()

    return (
        matrix,
        _with_brightway_keys(activity_dict, {database_name}),
        _with_brightway_keys(product_dict, {database_name}),
        elapsed,
    )


def _solve_with_umfpack(
    matrix: sparse.csc_matrix,
    demand: np.ndarray,
    *,
    repeat: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    ctx = UmfpackContext()
    memory: dict[str, Any] = {
        "rss_before_symbolic": _rss_bytes(),
    }

    start = time.perf_counter()
    symbolic = ctx.symbolic(matrix)
    symbolic_seconds = time.perf_counter() - start
    memory["rss_after_symbolic"] = _rss_bytes()
    memory["umfpack_symbolic_size"] = _umfpack_info_bytes(ctx, "UMFPACK_SYMBOLIC_SIZE")
    memory["umfpack_symbolic_peak"] = _umfpack_info_bytes(
        ctx, "UMFPACK_SYMBOLIC_PEAK_MEMORY"
    )

    start = time.perf_counter()
    try:
        ctx.numeric(matrix, symbolic)
    except TypeError:
        ctx.numeric(matrix)
    numeric_seconds = time.perf_counter() - start
    memory["rss_after_numeric"] = _rss_bytes()
    memory["umfpack_numeric_size"] = _umfpack_info_bytes(ctx, "UMFPACK_NUMERIC_SIZE")
    memory["umfpack_peak_memory"] = _umfpack_info_bytes(ctx, "UMFPACK_PEAK_MEMORY")
    memory["umfpack_variable_peak"] = _umfpack_info_bytes(ctx, "UMFPACK_VARIABLE_PEAK")

    solve_seconds: list[float] = []
    solution = np.empty_like(demand)
    for _ in range(max(1, repeat)):
        start = time.perf_counter()
        solution = ctx.solve(UMFPACK_A, matrix, demand)
        solve_seconds.append(time.perf_counter() - start)
    memory["rss_after_solve"] = _rss_bytes()

    return solution, {
        "symbolic_seconds": symbolic_seconds,
        "numeric_seconds": numeric_seconds,
        "solve_seconds": solve_seconds,
        "memory": memory,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load an ecoinvent Brightway database into a sparse technosphere "
            "matrix and solve one demand vector with UMFPACK."
        )
    )
    parser.add_argument("--project", default=DEFAULT_PROJECT)
    parser.add_argument("--database", default=DEFAULT_DATABASE)
    parser.add_argument(
        "--demand-code",
        help=(
            "Brightway code for the product demand in the selected database. "
            "Defaults to the first product row from the database."
        ),
    )
    parser.add_argument(
        "--demand-row",
        type=int,
        help=("Matrix product row to demand. Ignored when --demand-code is given."),
    )
    parser.add_argument("--amount", type=float, default=1.0)
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Number of UMFPACK backsolves to time after factorization.",
    )
    parser.add_argument(
        "--skip-residual",
        action="store_true",
        help="Skip the A*x - demand residual check.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    if args.project not in _project_names():
        raise KeyError(f"Brightway project not found: {args.project!r}")

    bd.projects.set_current(args.project)
    if args.database not in bd.databases:
        raise KeyError(
            f"Database {args.database!r} not found in Brightway project "
            f"{args.project!r}."
        )

    if UmfpackWarning is not None:
        warnings.filterwarnings("ignore", category=UmfpackWarning)

    db = bd.Database(args.database)
    rss_before_build = _rss_bytes()
    matrix, activity_dict, product_dict, build_seconds = _load_technosphere_matrix(
        args.database
    )
    rss_after_build = _rss_bytes()
    demand_key, demand_row = _select_demand_row(
        product_dict,
        database=args.database,
        demand_code=args.demand_code,
        demand_row=args.demand_row,
    )

    demand = np.zeros(matrix.shape[0], dtype=np.float64)
    demand[demand_row] = float(args.amount)

    solution, solve_stats = _solve_with_umfpack(
        matrix,
        demand,
        repeat=max(1, int(args.repeat)),
    )
    symbolic_seconds = solve_stats["symbolic_seconds"]
    numeric_seconds = solve_stats["numeric_seconds"]
    solve_seconds = solve_stats["solve_seconds"]
    memory = solve_stats["memory"]
    csc_bytes = _csc_storage_bytes(matrix)
    csc_total = sum(csc_bytes.values())

    print(f"Project: {args.project}")
    print(f"Database: {args.database}")
    print(f"Database entries: {len(db):,}")
    print(
        "Technosphere matrix: "
        f"shape={matrix.shape}, nnz={matrix.nnz:,}, "
        f"density={matrix.nnz / (matrix.shape[0] * matrix.shape[1]):.3e}"
    )
    print(f"Activities in matrix: {len(activity_dict):,}")
    print(f"Products in matrix: {len(product_dict):,}")
    print(f"Demand key: {_format_key(demand_key)}")
    print(f"Demand row: {demand_row:,}")
    print(f"Demand amount: {float(args.amount):.12g}")
    print(
        "Sparse matrix storage: "
        f"total={_format_bytes(csc_total)} "
        f"(data={_format_bytes(csc_bytes['data'])}, "
        f"indices={_format_bytes(csc_bytes['indices'])}, "
        f"indptr={_format_bytes(csc_bytes['indptr'])})"
    )
    if rss_before_build is not None and rss_after_build is not None:
        print(
            "Process RSS after matrix build: "
            f"{_format_bytes(rss_after_build)} "
            f"(delta={_format_bytes(rss_after_build - rss_before_build)})"
        )
    print(f"Matrix build time: {build_seconds:.6f} s")
    print(f"UMFPACK symbolic time: {symbolic_seconds:.6f} s")
    print(f"UMFPACK numeric time: {numeric_seconds:.6f} s")
    solve_times = ", ".join(f"{value:.6f} s" for value in solve_seconds)
    print(f"UMFPACK solve times: {solve_times}")
    print(
        "UMFPACK total factorize+first-solve time: "
        f"{symbolic_seconds + numeric_seconds + solve_seconds[0]:.6f} s"
    )
    if len(solve_seconds) > 1:
        print("UMFPACK mean backsolve time: " f"{float(np.mean(solve_seconds)):.6f} s")
    print(
        "UMFPACK symbolic object size: "
        f"{_format_bytes(memory['umfpack_symbolic_size'])}"
    )
    print(
        "UMFPACK symbolic peak memory: "
        f"{_format_bytes(memory['umfpack_symbolic_peak'])}"
    )
    print(
        "UMFPACK LU numeric object size: "
        f"{_format_bytes(memory['umfpack_numeric_size'])}"
    )
    print(
        "UMFPACK peak memory during factorization: "
        f"{_format_bytes(memory['umfpack_peak_memory'])}"
    )
    print(
        "UMFPACK variable peak memory during factorization: "
        f"{_format_bytes(memory['umfpack_variable_peak'])}"
    )
    if (
        memory["rss_before_symbolic"] is not None
        and memory["rss_after_symbolic"] is not None
    ):
        symbolic_rss_delta = (
            memory["rss_after_symbolic"] - memory["rss_before_symbolic"]
        )
        print("Process RSS symbolic delta: " f"{_format_bytes(symbolic_rss_delta)}")
    if (
        memory["rss_after_symbolic"] is not None
        and memory["rss_after_numeric"] is not None
    ):
        numeric_rss_delta = memory["rss_after_numeric"] - memory["rss_after_symbolic"]
        print(
            "Process RSS numeric factorization delta: "
            f"{_format_bytes(numeric_rss_delta)}"
        )
    if (
        memory["rss_after_numeric"] is not None
        and memory["rss_after_solve"] is not None
    ):
        solve_rss_delta = memory["rss_after_solve"] - memory["rss_after_numeric"]
        print("Process RSS solve delta: " f"{_format_bytes(solve_rss_delta)}")
    print(f"Solution nonzero entries: {int(np.count_nonzero(solution)):,}")
    print(f"Solution L1 norm: {float(np.linalg.norm(solution, ord=1)):.12g}")

    if not args.skip_residual:
        residual = matrix @ solution - demand
        print(
            "Residual infinity norm: "
            f"{float(np.linalg.norm(residual, ord=np.inf)):.3e}"
        )


if __name__ == "__main__":
    main()
