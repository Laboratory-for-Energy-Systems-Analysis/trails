from __future__ import annotations
import gc
import warnings
from typing import Any, Dict, List, Literal, TYPE_CHECKING

import bw2calc as bc
import numpy as np
from tqdm import tqdm

from scipy import sparse as sp
import sparse
import xarray as xr

from .bw_interface import (
    _extract_supply_fast,
    _extract_supply_fast_cached,
    _get_datapackage,
    _reference_product_from_activity_id,
    build_datapackage_for_year_from_trails,
)
from .iterative_solver import solve_many_rhs_jacobi_gmres
from .characterization import build_characterized_inventory
from .characterization import get_cf_matrix
from .characterization import get_cf_vector
from .edges_matrix import score_inventory_with_edges

try:
    from scikits.umfpack import UmfpackContext, UmfpackWarning, UMFPACK_A
except ImportError:  # pragma: no cover - optional dependency
    UmfpackContext = None  # type: ignore[assignment]
    UMFPACK_A = None  # type: ignore[assignment]
    UmfpackWarning = None  # type: ignore[assignment]

# Solver selection: prefer pypardiso (if available), then UMFPACK, else SciPy.
try:
    import pypardiso  # noqa: F401
except Exception:  # pragma: no cover - optional dependency
    _HAS_PYPARDISO = False
else:
    _HAS_PYPARDISO = True

if _HAS_PYPARDISO:
    SOLVER = "pypardiso"
elif UmfpackContext is not None:
    SOLVER = "umfpack"
else:
    SOLVER = "scipy"
    warnings.warn(
        "No accelerated sparse solver detected (pypardiso or scikits.umfpack). "
        "Falling back to SciPy's solver, which can be significantly slower.",
        RuntimeWarning,
        stacklevel=2,
    )

if TYPE_CHECKING:
    from .trails import Trails

if UmfpackWarning is not None:
    warnings.filterwarnings("ignore", category=UmfpackWarning)
    warnings.filterwarnings("ignore", module="scikits")

_CHAR_CACHE: dict = {}


def _default_ei_version(trails: Trails, ei_version: str | None) -> str:
    """Resolve a call-level LCIA version against Trails instance defaults."""
    if ei_version is not None:
        return str(ei_version)
    return str(
        getattr(
            trails,
            "default_ei_version",
            getattr(trails, "ei_version", "3.11"),
        )
    )


def _resolve_lca_method_defaults(
    trails: Trails,
    *,
    methods: List[str] | None,
    edges_methods: List[Any] | None,
    ei_version: str | None,
) -> tuple[List[str] | None, List[Any] | None, str]:
    """Resolve regular/EDGES method defaults for an LCA-style call.

    Explicit call arguments always win. If neither regular nor EDGES methods
    are provided, EDGES defaults are preferred for final scoring; regular
    default methods remain available for adaptive routing and static scoring.
    """
    resolved_methods = methods
    resolved_edges_methods = edges_methods
    if methods is None and edges_methods is None:
        default_edges = getattr(trails, "default_edges_methods", None)
        default_methods = getattr(trails, "default_methods", None)
        if default_edges:
            resolved_edges_methods = list(default_edges)
        elif default_methods:
            resolved_methods = list(default_methods)
    return (
        resolved_methods,
        resolved_edges_methods,
        _default_ei_version(trails, ei_version),
    )


class _IdentityProductMap:
    """Mapping-like object returning row index for integer product IDs."""

    def __init__(self, n: int) -> None:
        self.n = int(n)

    def __getitem__(self, key: int) -> int:
        i = int(key)
        if i < 0 or i >= self.n:
            raise KeyError(i)
        return i


def _get_mapping_arrays(
    mapping: Any,
) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """get mapping arrays.

    :param mapping: Value for `mapping`.
    :type mapping: Any
    :returns: Return value.
    :rtype: tuple[np.ndarray, np.ndarray] | tuple[None, None]"""
    if mapping:
        try:
            ids = np.fromiter(mapping.keys(), dtype=np.int64, count=len(mapping))
            pos = np.fromiter(mapping.values(), dtype=np.int64, count=len(mapping))
            return ids, pos
        except Exception:
            return None, None
    return None, None


def _map_activity_demands_to_products(
    lca_obj: Any,
    activity_demands: dict[int, float],
) -> dict[int, float]:
    """map activity demands to products.

    :param lca_obj: Value for `lca_obj`.
    :type lca_obj: Any
    :param activity_demands: Value for `activity_demands`.
    :type activity_demands: dict[int, float]
    :returns: Return value.
    :rtype: dict[int, float]
    :raises ValueError: If an error occurs."""
    if not activity_demands:
        return {}

    act_map = getattr(lca_obj.dicts, "activity", None)
    prod_map = getattr(lca_obj.dicts, "product", None)
    if not act_map or not prod_map:
        raise ValueError("LCA object missing dicts.activity or dicts.product")

    mapped: dict[int, float] = {}
    prod_map_keys = prod_map.keys()
    cache: dict[int, tuple[int, float]] = {}
    for act_id, amount in activity_demands.items():
        amt = float(amount)
        act_id = int(act_id)
        cached = cache.get(act_id)
        if cached is None:
            prod_id, prod_value = _reference_product_from_activity_id(lca_obj, act_id)
            cache[act_id] = (int(prod_id), float(prod_value))
        else:
            prod_id, prod_value = cached
        sign = -1.0 if prod_value < 0.0 else 1.0
        mapped[int(prod_id)] = mapped.get(int(prod_id), 0.0) + amt * sign

    return mapped


def _build_rhs_matrix_from_root_demands(
    *,
    per_root_demands: dict[int, dict[int, float]],
    product_dict: dict,
    n: int,
    min_amount: float,
) -> tuple[list[int], np.ndarray]:
    """build rhs matrix from root demands.

    :param per_root_demands: Value for `per_root_demands`.
    :type per_root_demands: dict[int, dict[int, float]]
    :param product_dict: Value for `product_dict`.
    :type product_dict: dict
    :param n: Value for `n`.
    :type n: int
    :param min_amount: Value for `min_amount`.
    :type min_amount: float
    :returns: Return value.
    :rtype: tuple[list[int], np.ndarray]
    :raises KeyError: If a product id cannot be resolved."""
    roots = [int(r) for r in per_root_demands.keys()]
    k = len(roots)
    if k == 0:
        return [], np.zeros((n, 0), dtype=np.float64)

    B = np.zeros((n, k), dtype=np.float64)

    # Fill columns by mapping product IDs -> row indices (product_dict maps id -> row index)
    for j, root in enumerate(roots):
        demand = per_root_demands[root]
        for prod_id, v in demand.items():
            v = float(v)
            if v == 0.0:
                continue
            # Assume prod_id keys are valid; skip check_demand to avoid overhead
            try:
                i = product_dict[int(prod_id)]
            except KeyError:
                # Keep behavior strict if desired; otherwise `continue`
                raise
            B[int(i), j] += v

    return roots, B


def _build_direct_technosphere_for_year(
    trails: Trails,
    year: int,
    cache: dict[
        int, tuple[sp.csc_matrix, _IdentityProductMap, dict[int, tuple[int, float]]]
    ],
) -> tuple[sp.csc_matrix, _IdentityProductMap, dict[int, tuple[int, float]]]:
    """Build (or fetch cached) direct technosphere matrix for one solve year."""
    y = int(year)
    cached = cache.get(y)
    if cached is not None:
        return cached

    context = trails._get_scenario_context(y)
    if context is None:
        raise RuntimeError(f"No scenario context available for year={y}")
    _scenario_year, _label, t = context

    if trails.A is None:
        raise RuntimeError("Trails.A is None")
    A_t = trails.A[int(t), :, :]
    coords = A_t.coords
    if coords.shape[0] == 3:
        act_idx = np.asarray(coords[1], dtype=np.int64)
        prod_idx = np.asarray(coords[2], dtype=np.int64)
    elif coords.shape[0] == 2:
        act_idx = np.asarray(coords[0], dtype=np.int64)
        prod_idx = np.asarray(coords[1], dtype=np.int64)
    else:
        raise ValueError(f"Unsupported A coords ndim={coords.shape[0]}")

    data = np.asarray(A_t.data, dtype=np.float64)
    n_acts = int(trails.A.shape[1])
    n_prods = int(trails.A.shape[2])

    A_csc = sp.coo_matrix((data, (prod_idx, act_idx)), shape=(n_prods, n_acts)).tocsc()
    if A_csc.shape[0] != A_csc.shape[1]:
        raise RuntimeError(
            f"Direct solver requires square technosphere; got {A_csc.shape} in year={y}"
        )

    result = (A_csc, _IdentityProductMap(A_csc.shape[0]), {})
    cache[y] = result
    return result


def _reference_product_from_activity_direct(
    A_csc: sp.csc_matrix,
    activity_id: int,
    cache: dict[int, tuple[int, float]],
) -> tuple[int, float]:
    """Resolve reference product row and sign from direct technosphere matrix."""
    act = int(activity_id)
    cached = cache.get(act)
    if cached is not None:
        return cached

    if act < 0 or act >= int(A_csc.shape[1]):
        raise KeyError(f"activity_id={act} outside technosphere columns")

    col = A_csc.getcol(act)
    rows = col.indices
    vals = np.asarray(col.data, dtype=np.float64)
    if rows.size == 0:
        raise ValueError(f"No technosphere entries found for activity_id={act}")

    prod_row: int | None = None
    prod_value = 0.0

    if 0 <= act < int(A_csc.shape[0]):
        mask = rows == act
        if np.any(mask):
            prod_row = act
            prod_value = float(vals[mask][0])

    if prod_row is None or prod_value == 0.0:
        k = int(np.argmin(np.abs(np.abs(vals) - 1.0)))
        prod_row = int(rows[k])
        prod_value = float(vals[k])

    result = (int(prod_row), float(prod_value))
    cache[act] = result
    return result


def _map_activity_demands_to_products_direct(
    A_csc: sp.csc_matrix,
    activity_demands: dict[int, float],
    ref_cache: dict[int, tuple[int, float]],
) -> dict[int, float]:
    """Map activity demand to product demand using direct technosphere columns."""
    if not activity_demands:
        return {}

    mapped: dict[int, float] = {}
    for act_id, amount in activity_demands.items():
        amt = float(amount)
        act_id = int(act_id)
        prod_id, prod_value = _reference_product_from_activity_direct(
            A_csc=A_csc, activity_id=act_id, cache=ref_cache
        )
        sign = -1.0 if prod_value < 0.0 else 1.0
        mapped[int(prod_id)] = mapped.get(int(prod_id), 0.0) + amt * sign
    return mapped


def _extract_supply_fast_direct(
    supply_array: np.ndarray,
    *,
    n_activities: int,
    min_amount: float = 0.0,
) -> dict[int, float]:
    """Extract non-zero activity supplies directly from solved vector."""
    supply = np.asarray(supply_array, dtype=np.float64)
    limit = min(int(n_activities), int(supply.shape[0]))
    if limit <= 0:
        return {}

    vals = supply[:limit]
    threshold = float(abs(min_amount))
    if threshold > 0.0:
        idx = np.flatnonzero(np.abs(vals) > threshold)
    else:
        idx = np.flatnonzero(vals != 0.0)
    if idx.size == 0:
        return {}

    return {int(i): float(vals[i]) for i in idx}


def solve_many_rhs_umfpack_factorized(
    A_csc: sp.csc_matrix,
    B: np.ndarray,
    *,
    cache: dict | None = None,
) -> np.ndarray:
    """Solve many rhs umfpack factorized.

    :param A_csc: Value for `A_csc`.
    :type A_csc: sp.csc_matrix
    :param B: Value for `B`.
    :type B: np.ndarray
    :param cache: Value for `cache`.
    :type cache: dict | None
    :returns: Return value.
    :rtype: np.ndarray
    :raises ValueError: If an error occurs."""
    if not sp.isspmatrix_csc(A_csc):
        A_csc = A_csc.tocsc()

    # UMFPACK is float64-centric; keep this strict to avoid silent slow paths or failures.
    if A_csc.dtype != np.float64:
        A_csc = A_csc.astype(np.float64)

    # UMFPACK can be sensitive to unsorted indices; enforce canonical CSC.
    A_csc.sort_indices()

    B = np.asarray(B)
    if B.ndim != 2:
        raise ValueError("B must be 2D (n, k)")
    if B.dtype != np.float64:
        B = B.astype(np.float64, copy=False)

    n, k = B.shape
    if A_csc.shape != (n, n):
        raise ValueError(f"Shape mismatch: A {A_csc.shape}, B {B.shape}")

    if UmfpackContext is None:
        # Fallback to SciPy when UMFPACK is unavailable.
        lu = sp.linalg.splu(A_csc)
        X = np.empty((n, k), dtype=np.float64)
        for j in range(k):
            X[:, j] = lu.solve(B[:, j])
        return X

    ctx = None
    symbolic = None
    pattern_sig = None
    if cache is not None:
        ctx = cache.get("ctx")
        symbolic = cache.get("symbolic")
        pattern_sig = cache.get("pattern_sig")
    if ctx is None:
        ctx = UmfpackContext()
        if cache is not None:
            cache["ctx"] = ctx

    # Reuse symbolic factorization when sparsity pattern is unchanged.
    sig = (
        A_csc.indptr.size,
        A_csc.indices.size,
        int(np.bitwise_xor.reduce(A_csc.indptr)),
        int(np.bitwise_xor.reduce(A_csc.indices)),
    )
    if cache is None or sig != pattern_sig or symbolic is None:
        symbolic = ctx.symbolic(A_csc)
        if cache is not None:
            cache["symbolic"] = symbolic
            cache["pattern_sig"] = sig

    # These calls frequently return None on success; failures raise exceptions.
    try:
        ctx.numeric(A_csc, symbolic)
    except TypeError:
        ctx.numeric(A_csc)

    X = np.empty((n, k), dtype=np.float64)

    # Reuse a preallocated x to reduce allocations
    for j in range(k):
        # In this scikits.umfpack version: solve(system, A, b, autoTranspose=True)
        X[:, j] = ctx.solve(UMFPACK_A, A_csc, B[:, j])

    return X


def lca_static(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str] | None = None,
    amount: float = 1.0,
    debug: bool = False,
    ei_version: str | None = None,
) -> None:
    """Run a conventional static LCA for one activity in one year.

    The static score uses regular LCIA methods only. If ``methods`` or
    ``ei_version`` are omitted, the corresponding ``Trails`` constructor
    defaults are used.

    :param trails: Trails instance containing the scenario matrices.
    :type trails: Trails
    :param year: Value for `year`.
    :type year: int
    :param fu_act_idx: Value for `fu_act_idx`.
    :type fu_act_idx: int
    :param methods: Regular LCIA methods. If omitted, uses
        ``Trails(..., methods=...)``.
    :type methods: List[str] | None
    :param amount: Value for `amount`.
    :type amount: float
    :param debug: Value for `debug`.
    :type debug: bool
    :param ei_version: LCIA data version. If omitted, uses
        ``Trails(..., ei_version=...)``.
    :type ei_version: str | None
    :raises ValueError: If no regular methods are available.
    :raises RuntimeError: If an error occurs."""
    if methods is None:
        default_methods = getattr(trails, "default_methods", None)
        methods = list(default_methods) if default_methods else None
    ei_version = _default_ei_version(trails, ei_version)
    if not methods:
        raise ValueError(
            "methods or Trails(..., methods=...) must be provided for static_lca()."
        )

    prev_inventory = trails.inventory
    prev_characterized = trails.characterized_inventory

    trails.reset_inventory(reset_scores=False)
    trails.static_score = None

    dp, _, _, _ = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        debug=debug,
    )

    # Build supply using the same activity->product mapping as temporal LCA
    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    inv = lca_obj.inventory  # SciPy sparse (flow_pos x act_pos)

    inv_coo = sparse.COO.from_scipy_sparse(inv)
    flow_pos = inv_coo.coords[0].astype(np.int64, copy=False)
    act_pos = inv_coo.coords[1].astype(np.int64, copy=False)

    act_map = getattr(lca_obj.dicts, "activity", None) or {}
    bio_map = getattr(lca_obj.dicts, "biosphere", None) or {}
    pos_to_act = {int(pos): int(act_id) for act_id, pos in act_map.items()}
    pos_to_flow = {int(pos): int(flow_id) for flow_id, pos in bio_map.items()}

    act_ids = np.array([pos_to_act.get(int(p), -1) for p in act_pos], dtype=np.int64)
    flow_ids = np.array([pos_to_flow.get(int(p), -1) for p in flow_pos], dtype=np.int64)

    valid = (act_ids >= 0) & (flow_ids >= 0)
    act_ids = act_ids[valid]
    flow_ids = flow_ids[valid]
    data = inv_coo.data[valid]

    n_acts = int(trails.A.shape[1]) if trails.A is not None else int(act_ids.max() + 1)
    n_flows = (
        int(trails.B.shape[2]) if trails.B is not None else int(flow_ids.max() + 1)
    )

    inv_coo = sparse.COO(
        coords=np.vstack([act_ids, flow_ids]),
        data=data,
        shape=(n_acts, n_flows),
    )

    trails.inventory = xr.DataArray(
        inv_coo[:, :, None],  # add year axis
        dims=("activity", "flow", "year"),
        coords={
            "activity": np.arange(n_acts, dtype=int),
            "flow": np.arange(n_flows, dtype=int),
            "year": np.array([int(year)], dtype=int),
        },
    )
    trails.demand = lca_obj.demand

    characterized_inventory = build_characterized_inventory(
        trails=trails,
        methods=methods,
        char_cache=_CHAR_CACHE,
        ei_version=ei_version,
    )

    if "method" in characterized_inventory.dims:
        inv = trails.inventory
        if inv is None:
            raise RuntimeError("Static inventory missing while scoring methods.")
        inv_data = inv.data
        if not isinstance(inv_data, sparse.COO):
            inv_data = sparse.COO.from_numpy(np.asarray(inv_data))
        flow_coords = inv_data.coords[1]
        vals = inv_data.data
        scores_list: list[float] = []
        for m in methods:
            cf = get_cf_vector(
                trails=trails,
                methods=[m],
                char_cache=_CHAR_CACHE,
                debug=debug,
                ei_version=ei_version,
            )
            score = float(np.dot(vals, cf[flow_coords]))
            scores_list.append(score)
        trails.static_score = scores_list
    else:
        trails.static_score = float(characterized_inventory.data.sum())

    trails.inventory = prev_inventory
    trails.characterized_inventory = prev_characterized


def lca(
    trails: Trails,
    methods: List[str] | None = None,
    edges_methods: List[Any] | None = None,
    show_progress: bool = True,
    attribute_to_roots: bool | None = None,
    *,
    store_inventory: bool = False,
    compute_score: bool = True,
    ei_version: str | None = None,
    solver_mode: Literal["bw2calc", "direct", "iterative"] = "iterative",
    iterative_rtol: float = 1e-3,
    iterative_atol: float = 0.0,
    iterative_restart: int | None = 50,
    iterative_maxiter: int | None = 300,
    iterative_use_guess: bool = True,
    iterative_preconditioner: Literal["jacobi", "ilu", "none"] = "jacobi",
    iterative_ilu_drop_tol: float = 1e-4,
    iterative_ilu_fill_factor: float = 10.0,
    edges_additional_topologies: dict[str, Any] | None = None,
    edges_strategies: list[str] | None = None,
    edges_reuse_cached_cfs: bool = True,
    inventory_workers: int | None = None,
) -> None:
    """Run temporal LCA from a previously built routing graph.

    ``temporal_routing()`` must be called first. If ``methods`` and
    ``edges_methods`` are both omitted and ``compute_score=True``, this function
    uses constructor defaults from the ``Trails`` instance. EDGES defaults are
    preferred for final scoring when both regular and EDGES defaults are set;
    explicit call arguments always override constructor defaults.

    :param trails: Trails instance containing a temporal routing graph.
    :type trails: Trails
    :param methods: Regular LCIA methods. If omitted together with
        ``edges_methods``, uses ``Trails(..., edges_methods=...)`` when set,
        otherwise ``Trails(..., methods=...)``.
    :type methods: List[str] | None
    :param edges_methods: Optional EDGES method definitions or method names for
        edge-level characterization factors. Mutually exclusive with ``methods``.
        If omitted together with ``methods``, uses
        ``Trails(..., edges_methods=...)`` when set.
    :type edges_methods: List[Any] | None
    :param show_progress: Value for `show_progress`.
    :type show_progress: bool
    :param attribute_to_roots: Value for `attribute_to_roots`. If ``None``,
        reuse the value used in ``trails.temporal_routing(...)``; if routing
        metadata is unavailable, defaults to ``True``.
    :type attribute_to_roots: bool | None
    :param store_inventory: Value for `store_inventory`.
    :type store_inventory: bool
    :param compute_score: Value for `compute_score`.
    :type compute_score: bool
    :param ei_version: LCIA data version for regular methods. If omitted, uses
        ``Trails(..., ei_version=...)``.
    :type ei_version: str | None
    :param solver_mode: Solver backend mode. Defaults to ``"iterative"``.
    :type solver_mode: Literal["bw2calc", "direct", "iterative"]
    :param iterative_rtol: Relative tolerance for iterative solves.
        Defaults to ``1e-3``.
    :type iterative_rtol: float
    :param iterative_atol: Value for `iterative_atol`.
    :type iterative_atol: float
    :param iterative_restart: Value for `iterative_restart`.
    :type iterative_restart: int | None
    :param iterative_maxiter: Value for `iterative_maxiter`.
    :type iterative_maxiter: int | None
    :param iterative_use_guess: Value for `iterative_use_guess`.
    :type iterative_use_guess: bool
    :param iterative_preconditioner: Iterative preconditioner mode.
    :type iterative_preconditioner: Literal["jacobi", "ilu", "none"]
    :param iterative_ilu_drop_tol: ILU drop tolerance (if ILU is selected).
    :type iterative_ilu_drop_tol: float
    :param iterative_ilu_fill_factor: ILU fill factor (if ILU is selected).
    :type iterative_ilu_fill_factor: float
    :param edges_additional_topologies: Optional topology definitions passed to
        EDGES when resolving regionalized CF locations.
    :type edges_additional_topologies: dict[str, Any] | None
    :param edges_strategies: Optional explicit EDGES matching strategy sequence.
    :type edges_strategies: list[str] | None
    :param edges_reuse_cached_cfs: Reuse EDGES matched CF templates across
        scenario years when supplier and consumer metadata signatures are the
        same. Numeric CF values are still evaluated for each scenario year. Set
        to ``False`` to force EDGES matching independently for every year.
    :type edges_reuse_cached_cfs: bool
    :param inventory_workers: Optional worker count for no-TD inventory batching.
    :type inventory_workers: int | None
    :raises RuntimeError: If an error occurs.
    :raises ValueError: If an error occurs."""
    if compute_score or methods is not None or edges_methods is not None:
        methods, edges_methods, ei_version = _resolve_lca_method_defaults(
            trails,
            methods=methods,
            edges_methods=edges_methods,
            ei_version=ei_version,
        )
    else:
        ei_version = _default_ei_version(trails, ei_version)

    if solver_mode not in ("bw2calc", "direct", "iterative"):
        raise ValueError(
            "solver_mode must be one of {'bw2calc', 'direct', 'iterative'}"
        )
    if iterative_preconditioner not in {"jacobi", "ilu", "none"}:
        raise ValueError(
            "iterative_preconditioner must be one of {'jacobi', 'ilu', 'none'}"
        )
    if methods and edges_methods:
        raise ValueError("methods and edges_methods are mutually exclusive.")

    edge_mode = bool(edges_methods)
    store_inventory_effective = bool(store_inventory or (edge_mode and compute_score))

    debug = bool(getattr(trails, "debug", False))

    # Routing must be run explicitly before LCA.
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing not initialized; run trails.temporal_routing(...) "
            "before trails.lca()."
        )

    routing_params = getattr(trails, "_routing_params", {}) or {}
    routing_attr_to_roots = getattr(trails, "_routing_attribute_to_roots", None)

    if attribute_to_roots is None:
        if routing_attr_to_roots is not None:
            attribute_to_roots = bool(routing_attr_to_roots)
        else:
            attribute_to_roots = True

    score_methods = (
        methods
        if (compute_score and not edge_mode and methods and len(methods) > 1)
        else None
    )
    trails.reset_inventory(
        attribute_to_roots=attribute_to_roots,
        score_methods=score_methods,
    )

    umfpack_cache: dict | None = (
        {} if (attribute_to_roots and SOLVER == "umfpack") else None
    )

    cf_matrix: np.ndarray | None = None
    cf_vectors: list[np.ndarray] | None = None
    if compute_score and not edge_mode:
        if not methods:
            raise ValueError("methods must be provided when compute_score=True.")
        cf_matrix = get_cf_matrix(
            trails=trails,
            methods=methods,
            char_cache=_CHAR_CACHE,
            debug=debug,
            ei_version=ei_version,
        )
        cf_vectors = [cf_matrix[i, :] for i in range(cf_matrix.shape[0])]

    # Ensure inventory builders are ready when we intend to store inventory data.
    if store_inventory_effective:
        required = (
            hasattr(trails, "_inventory_years")
            and trails._inventory_years is not None
            and hasattr(trails, "_inv_chunk_flows")
            and hasattr(trails, "_inv_chunk_values")
            and hasattr(trails, "_inv_chunk_len")
        )
        if not required:
            raise RuntimeError(
                "BUG: inventory storage is enabled but reset_inventory() did not "
                "initialize chunk-based inventory builders."
            )

    if routing_attr_to_roots is not None and routing_attr_to_roots != bool(
        attribute_to_roots
    ):
        raise RuntimeError(
            "Temporal routing was computed with attribute_to_roots="
            f"{routing_attr_to_roots}; rerun temporal_routing with "
            f"attribute_to_roots={attribute_to_roots}."
        )

    if not routing_params:
        raise RuntimeError(
            "Temporal routing parameters missing; rerun trails.temporal_routing(...)."
        )

    start_year_int = int(routing_params["start_year"])
    start_activity = int(routing_params["start_act_idx"])
    start_amount = float(routing_params["amount"])
    min_amount = float(routing_params.get("min_amount", 1e-18))
    start_context = trails._get_scenario_context(start_year_int)
    if start_context is None:
        raise RuntimeError(
            f"No scenario context available for start year={start_year_int}."
        )
    _, _, start_t = start_context
    start_activity_amount = trails._activity_amount_from_product_demand(
        int(start_t),
        start_activity,
        start_amount,
    )

    frontier: dict[tuple[int, int], float] = {}
    provenance: dict[tuple[int, int], dict[int, float]] = {}
    injected_supply_by_year_act: dict[tuple[int, int], float] = {}
    injected_supply_prov_by_year_act: dict[tuple[int, int], dict[int, float]] = {}
    functional_unit_is_frontier = False

    for node, data in graph.nodes(data=True):
        year = int(data.get("year"))
        act = int(data.get("act_idx"))
        depth = int(data.get("depth", 0))
        frontier_amt = float(data.get("frontier_amount") or 0.0)
        direct_bio_amt = float(data.get("direct_bio_amount") or 0.0)
        if (
            depth == 0
            and year == start_year_int
            and act == start_activity
            and frontier_amt
        ):
            functional_unit_is_frontier = True
        if frontier_amt:
            key = (year, act)
            frontier[key] = float(frontier.get(key, 0.0)) + frontier_amt
            if attribute_to_roots:
                roots = data.get("frontier_roots") or {}
                bucket = provenance.setdefault(key, {})
                for root_act, amt in roots.items():
                    bucket[int(root_act)] = float(
                        bucket.get(int(root_act), 0.0)
                    ) + float(amt)
        # Inject direct biosphere only for expanded nodes (not frontier nodes).
        if direct_bio_amt and not frontier_amt:
            key = (year, act)
            injected_supply_by_year_act[key] = (
                float(injected_supply_by_year_act.get(key, 0.0)) + direct_bio_amt
            )
            if attribute_to_roots:
                roots = data.get("direct_bio_roots") or {}
                bucket = injected_supply_prov_by_year_act.setdefault(key, {})
                for root_act, amt in roots.items():
                    bucket[int(root_act)] = float(
                        bucket.get(int(root_act), 0.0)
                    ) + float(amt)

    # Inject FU directly unless the depth-0 FU node is already frontier
    # (for example when max_depth=0).
    # Routing stores `amount` as the requested product amount, while direct
    # biosphere accumulation expects activity scaling. Use the same conversion
    # as temporal_routing() so non-unit production exchanges are not applied
    # twice.
    if not functional_unit_is_frontier:
        injected_supply_by_year_act[(start_year_int, start_activity)] = (
            float(
                injected_supply_by_year_act.get((start_year_int, start_activity), 0.0)
            )
            + start_activity_amount
        )
        injected_supply_prov_by_year_act.setdefault(
            (start_year_int, start_activity), {}
        )
        injected_supply_prov_by_year_act[(start_year_int, start_activity)][
            start_activity
        ] = (
            float(
                injected_supply_prov_by_year_act[(start_year_int, start_activity)].get(
                    start_activity, 0.0
                )
            )
            + start_activity_amount
        )

    # Frontier -> demand vectors (calendar years preserved)
    frontier_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(frontier_by_year.keys())

    datapackage_cache: Dict[tuple, Any] = {}

    root_demands_by_year: dict[int, dict[int, dict[int, float]]] = {}
    root_injected_by_year: dict[int, dict[int, dict[int, float]]] = {}

    if attribute_to_roots:
        for (y, a), total_amt in frontier.items():
            root_map = provenance.get((y, a))
            if not root_map:
                root_map = {int(a): float(total_amt)}
            year_bucket = root_demands_by_year.setdefault(int(y), {})
            for root, amt in root_map.items():
                root_bucket = year_bucket.setdefault(int(root), {})
                root_bucket[int(a)] = root_bucket.get(int(a), 0.0) + float(amt)

        for (y, a), total_amt in injected_supply_by_year_act.items():
            root_map = injected_supply_prov_by_year_act.get((y, a))
            if not root_map:
                root_map = {int(a): float(total_amt)}
            year_bucket = root_injected_by_year.setdefault(int(y), {})
            for root, amt in root_map.items():
                root_bucket = year_bucket.setdefault(int(root), {})
                root_bucket[int(a)] = root_bucket.get(int(a), 0.0) + float(amt)

    pbar = None
    if show_progress:
        pbar = tqdm(
            total=len(candidate_years),
            desc="Temporal LCA: solve years",
            unit="year",
            leave=True,
        )

    direct_matrix_cache: dict[
        int, tuple[sp.csc_matrix, _IdentityProductMap, dict[int, tuple[int, float]]]
    ] = {}

    for solve_year in candidate_years:
        root_ids = None
        root_supply_matrix = None

        solve_year = int(solve_year)
        demand_vector = np.asarray(frontier_by_year[solve_year])
        nonzero_indices = np.where(demand_vector != 0.0)[0]
        if nonzero_indices.size == 0:
            if pbar is not None:
                pbar.update(1)
            continue

        activity_demand = {int(i): float(demand_vector[i]) for i in nonzero_indices}
        n_activities = int(trails.A.shape[1]) if trails.A is not None else 0
        supplies: list[tuple[Dict[int, float], int | None]] = []
        need_supply_dicts = bool(store_inventory_effective) or (
            bool(compute_score) and not edge_mode and not bool(attribute_to_roots)
        )
        use_direct_solver = bool(solver_mode in {"direct", "iterative"})
        use_iterative_solver = bool(solver_mode == "iterative")

        if use_direct_solver:
            A_csc, product_dict, ref_prod_cache = _build_direct_technosphere_for_year(
                trails=trails,
                year=solve_year,
                cache=direct_matrix_cache,
            )
            functional_unit_demand = _map_activity_demands_to_products_direct(
                A_csc=A_csc,
                activity_demands=activity_demand,
                ref_cache=ref_prod_cache,
            )
            if not functional_unit_demand:
                if pbar is not None:
                    pbar.update(1)
                continue

            if attribute_to_roots:
                per_root_demands_raw = root_demands_by_year.get(solve_year, {})
                per_root_demands: dict[int, dict[int, float]] = {}
                for root_act, demand in per_root_demands_raw.items():
                    mapped = _map_activity_demands_to_products_direct(
                        A_csc=A_csc,
                        activity_demands=demand,
                        ref_cache=ref_prod_cache,
                    )
                    if mapped:
                        per_root_demands[int(root_act)] = mapped
                roots, rhs_matrix = _build_rhs_matrix_from_root_demands(
                    per_root_demands=per_root_demands,
                    product_dict=product_dict,  # type: ignore[arg-type]
                    n=A_csc.shape[0],
                    min_amount=float(min_amount),
                )

                if roots:
                    if use_iterative_solver:
                        root_supply_matrix = solve_many_rhs_jacobi_gmres(
                            A_csc,
                            rhs_matrix,
                            rtol=float(iterative_rtol),
                            atol=float(iterative_atol),
                            restart=iterative_restart,
                            maxiter=iterative_maxiter,
                            use_guess=bool(iterative_use_guess),
                            preconditioner_mode=iterative_preconditioner,
                            ilu_drop_tol=float(iterative_ilu_drop_tol),
                            ilu_fill_factor=float(iterative_ilu_fill_factor),
                        )
                    else:
                        root_supply_matrix = solve_many_rhs_umfpack_factorized(
                            A_csc, rhs_matrix, cache=umfpack_cache
                        )
                    root_ids = np.asarray(roots, dtype=np.int64)
                    if need_supply_dicts:
                        for j, root_act in enumerate(roots):
                            supply_vec = root_supply_matrix[:, j]
                            supply_total = _extract_supply_fast_direct(
                                supply_vec,
                                n_activities=n_activities,
                                min_amount=min_amount,
                            )
                            if supply_total:
                                supplies.append((supply_total, int(root_act)))
            else:
                rhs = np.zeros((A_csc.shape[0], 1), dtype=np.float64)
                for prod_id, v in functional_unit_demand.items():
                    rhs[int(product_dict[int(prod_id)]), 0] += float(v)
                if use_iterative_solver:
                    X = solve_many_rhs_jacobi_gmres(
                        A_csc,
                        rhs,
                        rtol=float(iterative_rtol),
                        atol=float(iterative_atol),
                        restart=iterative_restart,
                        maxiter=iterative_maxiter,
                        use_guess=bool(iterative_use_guess),
                        preconditioner_mode=iterative_preconditioner,
                        ilu_drop_tol=float(iterative_ilu_drop_tol),
                        ilu_fill_factor=float(iterative_ilu_fill_factor),
                    )
                else:
                    X = solve_many_rhs_umfpack_factorized(A_csc, rhs, cache=None)
                supply_vec = X[:, 0]
                if need_supply_dicts:
                    supply_total = _extract_supply_fast_direct(
                        supply_vec,
                        n_activities=n_activities,
                        min_amount=min_amount,
                    )
                    if supply_total:
                        supplies.append((supply_total, None))

        else:
            dp, _, _, _ = _get_datapackage(
                dp_cache=datapackage_cache,
                trails=trails,
                year=solve_year,
                zero_bio=True,
                debug=debug,
            )

            lca_obj = bc.LCA(demand=activity_demand, data_objs=[dp])
            # The zero-bio fast path intentionally omits biosphere resources.
            # bw2calc warns on empty biosphere in this case; suppress this
            # specific warning while keeping all other warnings visible.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=(
                        r"No valid biosphere flows found\. No inventory results can "
                        r"be calculated, `lcia` will raise an error"
                    ),
                    category=UserWarning,
                    module=r"bw2calc\.lca_base",
                )
                lca_obj.load_lci_data()  # build matrices + dicts

            # Cache mappings once per year
            act_map = getattr(lca_obj.dicts, "activity", None)
            prod_map = getattr(lca_obj.dicts, "product", None)

            act_ids, positions = _get_mapping_arrays(act_map)
            # For RHS building we need the product-id -> row-index dict.
            product_dict = prod_map  # mapping-like

            functional_unit_demand = _map_activity_demands_to_products(
                lca_obj,
                activity_demand,
            )
            if not functional_unit_demand:
                if pbar is not None:
                    pbar.update(1)
                continue

            # Prepare the technosphere matrix once per year (UMFPACK expects CSC).
            A_csc = lca_obj.technosphere_matrix
            if not sp.isspmatrix_csc(A_csc):
                A_csc = A_csc.tocsc()

            if attribute_to_roots:
                # Build one dense RHS matrix for all roots in this year
                per_root_demands_raw = root_demands_by_year.get(solve_year, {})
                per_root_demands: dict[int, dict[int, float]] = {}
                for root_act, demand in per_root_demands_raw.items():
                    mapped = _map_activity_demands_to_products(
                        lca_obj,
                        demand,
                    )
                    if mapped:
                        per_root_demands[int(root_act)] = mapped
                roots, rhs_matrix = _build_rhs_matrix_from_root_demands(
                    per_root_demands=per_root_demands,
                    product_dict=product_dict,
                    n=A_csc.shape[0],
                    min_amount=float(min_amount),
                )

                if roots:
                    if SOLVER == "pypardiso":
                        # Keep the PARDISO path for environments that rely on it.
                        for root_act in roots:
                            root_demand = per_root_demands[root_act]
                            lca_obj.build_demand_array(root_demand)
                            lca_obj.supply_array = lca_obj.solve_linear_system()
                            if need_supply_dicts:
                                if act_ids is None or positions is None:
                                    supply_total = _extract_supply_fast(
                                        lca_obj, min_amount
                                    )
                                else:
                                    supply_total = _extract_supply_fast_cached(
                                        lca_obj.supply_array,
                                        act_ids,
                                        positions,
                                        min_amount,
                                    )
                                if supply_total:
                                    supplies.append((supply_total, int(root_act)))
                    else:
                        # UMFPACK/SciPy: factorize once and solve all RHS vectors.
                        root_supply_matrix = solve_many_rhs_umfpack_factorized(
                            A_csc, rhs_matrix, cache=umfpack_cache
                        )
                        root_ids = np.asarray(
                            roots, dtype=np.int64
                        )  # Column order of RHS and solution.

                        if need_supply_dicts:
                            for j, root_act in enumerate(roots):
                                supply_vec = root_supply_matrix[:, j]
                                if act_ids is None or positions is None:
                                    # Fallback: emulate lca_obj.supply_array for extractor
                                    lca_obj.supply_array = supply_vec
                                    supply_total = _extract_supply_fast(
                                        lca_obj, min_amount
                                    )
                                else:
                                    supply_total = _extract_supply_fast_cached(
                                        supply_vec, act_ids, positions, min_amount
                                    )
                                if supply_total:
                                    supplies.append((supply_total, int(root_act)))

                else:
                    if pbar is not None:
                        pbar.update(1)
                    continue
            else:
                # Original single-demand path (no per-root attribution)
                lca_obj.build_demand_array(functional_unit_demand)
                if SOLVER != "pypardiso":
                    lca_obj.decompose_technosphere()
                lca_obj.supply_array = lca_obj.solve_linear_system()

                if need_supply_dicts:
                    if act_ids is None or positions is None:
                        supply_total = _extract_supply_fast(lca_obj, min_amount)
                    else:
                        supply_total = _extract_supply_fast_cached(
                            lca_obj.supply_array, act_ids, positions, min_amount
                        )
                    if supply_total:
                        supplies.append((supply_total, None))

        # ---- Inventory ----
        if supplies and store_inventory_effective:
            trails.accumulate_temporalized_biosphere_inventory_batch(
                base_year=solve_year,
                supplies=supplies,
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
                workers=inventory_workers,
            )

        # ---- Scores ----
        if compute_score and not edge_mode:
            assert cf_vectors is not None
            # Build dense supply matrix (n_acts, n_roots_this_year).
            # Injected supplies remain dictionary-based because they are usually small.
            if attribute_to_roots:
                # Only call the matrix scorer when we solved a multi-RHS system
                # and the root ordering matches the solution columns.
                if root_ids is not None and root_supply_matrix is not None:
                    if cf_matrix is None:
                        raise RuntimeError("CF matrix missing while scoring methods.")
                    trails.accumulate_temporalized_biosphere_score_matrix_multi(
                        base_year=solve_year,
                        supply_matrix=root_supply_matrix,
                        root_activities=root_ids,
                        cf_matrix=cf_matrix,
                        method_indices=np.arange(cf_matrix.shape[0], dtype=np.int64),
                        min_amount=float(min_amount),
                        use_temporal_distributions=True,
                        debug=debug,
                    )

            else:
                # Non-root mode: score all supply blocks for the year
                # (solved supply + injected direct additions).
                for supply_dict, _ in supplies:
                    if not supply_dict:
                        continue
                    for method_idx, cf_vec in enumerate(cf_vectors):
                        trails.accumulate_temporalized_biosphere_score(
                            base_year=solve_year,
                            supply_by_activity=supply_dict,
                            cf=cf_vec,
                            min_amount=float(min_amount),
                            store_activity=None,
                            debug=debug,
                            method_idx=method_idx,
                        )

        if pbar is not None:
            pbar.update(1)
    if pbar is not None:
        pbar.close()

    # Injected/direct biosphere supplies should keep their raw calendar year,
    # while technosphere solving remains scenario-mapped.
    if attribute_to_roots:
        if store_inventory_effective:
            for raw_year, per_root_injected in root_injected_by_year.items():
                supplies_batch = [
                    (injected_supply, int(root_act))
                    for root_act, injected_supply in per_root_injected.items()
                    if injected_supply
                ]
                if not supplies_batch:
                    continue
                trails.accumulate_temporalized_biosphere_inventory_batch(
                    base_year=int(raw_year),
                    supplies=supplies_batch,
                    min_amount=float(min_amount),
                    use_temporal_distributions=True,
                    debug=debug,
                    workers=inventory_workers,
                )
        if compute_score and not edge_mode:
            assert cf_vectors is not None
            for raw_year, per_root_injected in root_injected_by_year.items():
                for root_act, injected_supply in per_root_injected.items():
                    if not injected_supply:
                        continue
                    for method_idx, cf_vec in enumerate(cf_vectors):
                        trails.accumulate_temporalized_biosphere_score(
                            base_year=int(raw_year),
                            supply_by_activity=injected_supply,
                            cf=cf_vec,
                            min_amount=float(min_amount),
                            store_activity=int(root_act),
                            debug=debug,
                            method_idx=method_idx,
                        )
    else:
        injected_by_raw_year: dict[int, dict[int, float]] = {}
        for (raw_year, act_idx), v in injected_supply_by_year_act.items():
            y = int(raw_year)
            a = int(act_idx)
            bucket = injected_by_raw_year.setdefault(y, {})
            bucket[a] = bucket.get(a, 0.0) + float(v)

        if store_inventory_effective:
            for raw_year, injected_supply in injected_by_raw_year.items():
                if not injected_supply:
                    continue
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=int(raw_year),
                    supply_by_activity=injected_supply,
                    min_amount=float(min_amount),
                    store_activity=None,
                    debug=debug,
                )
        if compute_score and not edge_mode:
            assert cf_vectors is not None
            for raw_year, injected_supply in injected_by_raw_year.items():
                if not injected_supply:
                    continue
                for method_idx, cf_vec in enumerate(cf_vectors):
                    trails.accumulate_temporalized_biosphere_score(
                        base_year=int(raw_year),
                        supply_by_activity=injected_supply,
                        cf=cf_vec,
                        min_amount=float(min_amount),
                        store_activity=None,
                        debug=debug,
                        method_idx=method_idx,
                    )

    # Large local containers are no longer needed past this point; clear them
    # before finalize_* to reduce peak RSS during inventory/scores materialization.
    frontier.clear()
    provenance.clear()
    injected_supply_by_year_act.clear()
    injected_supply_prov_by_year_act.clear()
    frontier_by_year.clear()
    root_demands_by_year.clear()
    root_injected_by_year.clear()
    datapackage_cache.clear()
    direct_matrix_cache.clear()
    if "injected_by_raw_year" in locals():
        injected_by_raw_year.clear()
    gc.collect()

    if store_inventory_effective:
        trails.finalize_inventory()

    if compute_score:
        if edge_mode:
            score_inventory_with_edges(
                trails,
                list(edges_methods or []),
                additional_topologies=edges_additional_topologies,
                strategies=edges_strategies,
                reuse_cached_cfs=edges_reuse_cached_cfs,
                show_progress=show_progress,
            )
            if not store_inventory:
                trails.inventory = None
        else:
            trails.finalize_scores()

        if trails.scores is None:
            raise RuntimeError(
                "compute_score=True but trails.scores is still None. "
                "This indicates lca() did not finalize scores correctly."
            )

    # Characterized inventory is optional and only meaningful when scoring.
    if store_inventory_effective and compute_score and not edge_mode:
        build_characterized_inventory(
            trails=trails, methods=methods, char_cache=_CHAR_CACHE
        )


def build_temporal_sankey_tree(
    trails: Trails,
    *,
    root_year: int | None = None,
    root_act_idx: int | None = None,
    max_depth: int | None = None,
    min_amount: float = 0.0,
    sort_children: bool = True,
) -> dict[str, Any]:
    """Build temporal sankey tree.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param root_year: Value for `root_year`.
    :type root_year: int | None
    :param root_act_idx: Value for `root_act_idx`.
    :type root_act_idx: int | None
    :param max_depth: Value for `max_depth`.
    :type max_depth: int | None
    :param min_amount: Value for `min_amount`.
    :type min_amount: float
    :param sort_children: Value for `sort_children`.
    :type sort_children: bool
    :returns: Return value.
    :rtype: dict[str, Any]
    :raises RuntimeError: If an error occurs.
    :raises ValueError: If an error occurs."""
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing graph missing; run trails.temporal_routing(...) first."
        )

    routing_params = getattr(trails, "_routing_params", {}) or {}
    if root_year is None:
        if "start_year" not in routing_params:
            raise RuntimeError(
                "Root year not provided and routing params missing; "
                "rerun temporal_routing or pass root_year."
            )
        root_year = int(routing_params["start_year"])
    if root_act_idx is None:
        if "start_act_idx" not in routing_params:
            raise RuntimeError(
                "Root activity not provided and routing params missing; "
                "rerun temporal_routing or pass root_act_idx."
            )
        root_act_idx = int(routing_params["start_act_idx"])

    min_amount = float(min_amount)

    root_key = None
    for node, data in graph.nodes(data=True):
        if int(data.get("depth", -1)) != 0:
            continue
        if int(data.get("year", -9999)) != int(root_year):
            continue
        if int(data.get("act_idx", -1)) != int(root_act_idx):
            continue
        root_key = node
        break

    if root_key is None:
        raise ValueError(
            "Root node not found in temporal routing graph for "
            f"(year={int(root_year)}, act_idx={int(root_act_idx)})."
        )

    def _node_payload(node_key: tuple) -> dict[str, Any]:
        """node payload.

        :param node_key: Value for `node_key`.
        :type node_key: tuple
        :returns: Return value.
        :rtype: dict[str, Any]"""
        data = graph.nodes[node_key]
        return {
            "key": node_key,
            "year": int(data.get("year")),
            "depth": int(data.get("depth")),
            "act_idx": int(data.get("act_idx")),
            "name": data.get("name") or "",
            "reference_product": data.get("reference_product") or "",
            "location": data.get("location") or "",
            "amount": float(data.get("amount") or 0.0),
            "frontier_amount": float(data.get("frontier_amount") or 0.0),
            "direct_bio_amount": float(data.get("direct_bio_amount") or 0.0),
        }

    def _build_tree(node_key: tuple) -> dict[str, Any]:
        """build tree.

        :param node_key: Value for `node_key`.
        :type node_key: tuple
        :returns: Return value.
        :rtype: dict[str, Any]"""
        node_data = graph.nodes[node_key]
        depth = int(node_data.get("depth"))
        if max_depth is not None and depth >= int(max_depth):
            return {"node": _node_payload(node_key), "children": []}

        children = []
        for child in graph.successors(node_key):
            edge_data = graph.edges[node_key, child]
            edge_amt = float(edge_data.get("amount") or 0.0)
            if abs(edge_amt) < min_amount:
                continue
            child_tree = _build_tree(child)
            children.append(
                {
                    "edge_amount": edge_amt,
                    "node": child_tree["node"],
                    "children": child_tree["children"],
                }
            )

        if sort_children:
            children.sort(key=lambda c: abs(c["edge_amount"]), reverse=True)

        return {"node": _node_payload(node_key), "children": children}

    return _build_tree(root_key)


def score_temporal_graph_nodes(
    trails: Trails,
    methods: List[str] | None = None,
    *,
    min_amount: float = 0.0,
    show_progress: bool = True,
    ei_version: str | None = None,
) -> dict[tuple, float]:
    """Estimate individual temporal graph node scores.

    This diagnostic helper scores routed graph nodes without mutating the main
    temporal score arrays. If ``methods`` or ``ei_version`` are omitted, regular
    constructor defaults from the ``Trails`` instance are used.

    :param trails: Trails instance containing a temporal routing graph.
    :type trails: Trails
    :param methods: Regular LCIA methods. If omitted, uses
        ``Trails(..., methods=...)``.
    :type methods: List[str] | None
    :param min_amount: Value for `min_amount`.
    :type min_amount: float
    :param show_progress: Value for `show_progress`.
    :type show_progress: bool
    :param ei_version: LCIA data version. If omitted, uses
        ``Trails(..., ei_version=...)``.
    :type ei_version: str | None
    :returns: Return value.
    :rtype: dict[tuple, float]
    :raises RuntimeError: If an error occurs."""
    graph = getattr(trails, "graph", None)
    if graph is None:
        raise RuntimeError(
            "Temporal routing graph missing; run trails.temporal_routing(...) first."
        )

    if trails.A is None or trails.B is None:
        raise RuntimeError("Trails matrices missing; run setup before scoring.")

    if methods is None:
        default_methods = getattr(trails, "default_methods", None)
        methods = list(default_methods) if default_methods else None
    if not methods:
        raise ValueError(
            "methods or Trails(..., methods=...) must be provided to score nodes."
        )
    ei_version = _default_ei_version(trails, ei_version)

    cf = get_cf_vector(
        trails=trails,
        methods=methods,
        char_cache=_CHAR_CACHE,
        debug=bool(getattr(trails, "debug", False)),
        ei_version=ei_version,
    )

    min_amount = float(min_amount)

    # Caches
    dp_cache: Dict[tuple, Any] = {}
    char_row_cache: dict[int, np.ndarray] = {}  # t -> per-activity coeffs
    has_bio_cache: dict[tuple[int, int], bool] = {}
    supply_cache: dict[tuple[int, int], dict[int, float]] = {}

    def _char_row_for_year(year: int) -> np.ndarray:
        """char row for year.

        :param year: Value for `year`.
        :type year: int
        :returns: Return value.
        :rtype: np.ndarray"""
        context = trails._get_scenario_context(int(year))
        if context is None:
            return np.zeros(int(trails.A.shape[1]), dtype=np.float64)
        _scenario_year, _label, t = context
        if t in char_row_cache:
            return char_row_cache[t]
        B_t = trails.B[t, :, :]
        # B_t is (activity, flow); multiply by CF (flow) -> per-activity coefficients
        coeff = B_t @ cf  # type: ignore[operator]
        coeff = np.asarray(coeff, dtype=np.float64)
        char_row_cache[t] = coeff
        return coeff

    def _score_direct_td(year: int, act: int, amount: float) -> float:
        """score direct td.

        :param year: Value for `year`.
        :type year: int
        :param act: Value for `act`.
        :type act: int
        :param amount: Value for `amount`.
        :type amount: float
        :returns: Return value.
        :rtype: float"""
        if amount == 0.0:
            return 0.0

        # Intercept score appends to avoid mutating trails.scores
        total = 0.0
        original_append = trails._append_score_entry
        original_append_bulk = getattr(trails, "_append_scores_bulk", None)

        def _capture_append(
            act_idx: int,
            _year_val: int,
            value: float,
            *,
            root_activity: int | None = None,
        ) -> None:
            """capture append.

            :param act_idx: Value for `act_idx`.
            :type act_idx: int
            :param _year_val: Value for `_year_val`.
            :type _year_val: int
            :param value: Value for `value`.
            :type value: float
            :param root_activity: Value for `root_activity`.
            :type root_activity: int | None"""
            nonlocal total
            total += float(value)

        def _capture_append_bulk(
            act_idx: int,
            year_idx: int,
            value: float,
            *,
            root_activity: int | None = None,
        ) -> None:
            """capture append bulk.

            :param act_idx: Value for `act_idx`.
            :type act_idx: int
            :param year_idx: Value for `year_idx`.
            :type year_idx: int
            :param value: Value for `value`.
            :type value: float
            :param root_activity: Value for `root_activity`.
            :type root_activity: int | None"""
            nonlocal total
            total += float(value)

        trails._append_score_entry = _capture_append  # type: ignore[assignment]
        if original_append_bulk is not None:
            trails._append_scores_bulk = _capture_append_bulk  # type: ignore[assignment]

        try:
            trails.accumulate_temporalized_biosphere_score(
                base_year=int(year),
                supply_by_activity={int(act): float(amount)},
                cf=cf,
                min_amount=min_amount,
                store_activity=None,
                use_temporal_distributions=True,
                debug=bool(getattr(trails, "debug", False)),
            )
        finally:
            trails._append_score_entry = original_append  # type: ignore[assignment]
            if original_append_bulk is not None:
                trails._append_scores_bulk = original_append_bulk  # type: ignore[assignment]

        return float(total)

    def _solve_supply(year: int, act: int, amount: float) -> dict[int, float]:
        """solve supply.

        :param year: Value for `year`.
        :type year: int
        :param act: Value for `act`.
        :type act: int
        :param amount: Value for `amount`.
        :type amount: float
        :returns: Return value.
        :rtype: dict[int, float]"""
        key = (int(year), int(act))
        if key in supply_cache:
            if amount == 1.0:
                return supply_cache[key]
            return {a: float(v) * float(amount) for a, v in supply_cache[key].items()}

        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=int(year),
            zero_bio=False,
            debug=bool(getattr(trails, "debug", False)),
        )

        lca_obj = bc.LCA(demand={int(act): 1.0}, data_objs=[dp])
        lca_obj.load_lci_data()

        activity_demand = {int(act): 1.0}
        functional_unit_demand = _map_activity_demands_to_products(
            lca_obj, activity_demand
        )
        if not functional_unit_demand:
            supply_cache[key] = {}
            return {}

        act_map = getattr(lca_obj.dicts, "activity", None)
        act_ids, positions = _get_mapping_arrays(act_map)

        lca_obj.build_demand_array(functional_unit_demand)
        if not bc.PYPARDISO:
            lca_obj.decompose_technosphere()
        lca_obj.supply_array = lca_obj.solve_linear_system()

        if act_ids is None or positions is None:
            supply_total = _extract_supply_fast(lca_obj, 0.0)
        else:
            supply_total = _extract_supply_fast_cached(
                lca_obj.supply_array, act_ids, positions, 0.0
            )

        supply_cache[key] = supply_total
        if amount == 1.0:
            return supply_total
        return {a: float(v) * float(amount) for a, v in supply_total.items()}

    def _score_frontier(year: int, act: int, amount: float) -> float:
        """score frontier.

        :param year: Value for `year`.
        :type year: int
        :param act: Value for `act`.
        :type act: int
        :param amount: Value for `amount`.
        :type amount: float
        :returns: Return value.
        :rtype: float"""
        if amount == 0.0:
            return 0.0
        supply = _solve_supply(year, act, amount)
        if not supply:
            return 0.0
        coeff = _char_row_for_year(year)
        total = 0.0
        for a, v in supply.items():
            a = int(a)
            if a < 0 or a >= coeff.size:
                continue
            total += float(v) * float(coeff[a])
        return float(total)

    def _is_frontier(node_data: dict) -> bool:
        """is frontier.

        :param node_data: Value for `node_data`.
        :type node_data: dict
        :returns: Return value.
        :rtype: bool"""
        frontier_amt = float(node_data.get("frontier_amount") or 0.0)
        if frontier_amt != 0.0:
            return True
        return False

    def _direct_amount(node_data: dict) -> float:
        """direct amount.

        :param node_data: Value for `node_data`.
        :type node_data: dict
        :returns: Return value.
        :rtype: float"""
        direct_amt = float(node_data.get("direct_bio_amount") or 0.0)
        if direct_amt != 0.0:
            return direct_amt
        depth = int(node_data.get("depth", 0))
        if depth == 0:
            year = int(node_data.get("year"))
            act = int(node_data.get("act_idx"))
            if trails._has_direct_biosphere(
                int(trails._map_year_to_scenario_year(year)), act, has_bio_cache
            ):
                return float(node_data.get("amount") or 0.0)
        return 0.0

    nodes = list(graph.nodes(data=True))
    pbar = None
    if show_progress:
        pbar = tqdm(
            total=len(nodes),
            desc="Score temporal nodes",
            unit="node",
            leave=True,
        )

    node_scores: dict[tuple, float] = {}

    for node_key, data in nodes:
        year = int(data.get("year"))
        act = int(data.get("act_idx"))
        if _is_frontier(data):
            amount = float(data.get("frontier_amount") or 0.0)
            if abs(amount) >= min_amount:
                node_scores[node_key] = _score_frontier(year, act, amount)
            else:
                node_scores[node_key] = 0.0
        else:
            amount = _direct_amount(data)
            if abs(amount) >= min_amount:
                node_scores[node_key] = _score_direct_td(year, act, amount)
            else:
                node_scores[node_key] = 0.0

        if pbar is not None:
            pbar.update(1)

    if pbar is not None:
        pbar.close()

    return node_scores
