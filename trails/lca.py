import numpy as np
import bw_processing as bwp
from typing import Dict, Tuple, List, Any
from collections import defaultdict

import bw2calc as bc
from tqdm.auto import tqdm


from .trails import Trails
from .lcia import fill_characterization_factors_matrices, get_lcia_methods
from .temporal_distributions import TemporalDistribution

import logging
logger = logging.getLogger(__name__)


def _ij_from_coords(X_t):
    coords = X_t.coords
    if coords.shape[0] == 3:
        i_idx = np.asarray(coords[1], dtype=np.int64)
        j_idx = np.asarray(coords[2], dtype=np.int64)
    elif coords.shape[0] == 2:
        i_idx = np.asarray(coords[0], dtype=np.int64)
        j_idx = np.asarray(coords[1], dtype=np.int64)
    else:
        raise ValueError(f"Unsupported coords ndim={coords.shape[0]} for sparse array")
    return i_idx, j_idx


def _resolve_matrix_label(trails: Trails, year: int) -> tuple[str, int]:
    label_for_matrix = str(year)
    if label_for_matrix not in trails.scenario_index:
        years = np.array([int(lbl) for lbl in trails.scenario_labels])
        idx = int(np.argmin(np.abs(years - year)))
        label_for_matrix = trails.scenario_labels[idx]
        print(
            f"⚠️ Matrix slice for year {year} not found in A/B; "
            f"using nearest available matrix year {label_for_matrix}."
        )

    t = trails.scenario_index[label_for_matrix]
    return label_for_matrix, t


def _select_metadata_label(trails: Trails, label_for_matrix: str) -> str:
    if label_for_matrix in trails.activity_indices:
        return label_for_matrix
    return _nearest_metadata_label_for_year(trails, int(label_for_matrix))


def _apply_temporal_scaling_to_A(
    trails: Trails,
    year: int,
    A_t,
    debug: bool,
    collapse_tech_temporal_scaling: bool,
):
    A_signed = A_t.data
    if not collapse_tech_temporal_scaling:
        return A_signed

    A_act_idx, A_prod_idx = _ij_from_coords(A_t)
    multipliers_A = np.ones_like(A_signed, dtype=np.float64)

    for n in range(len(A_signed)):
        act = int(A_act_idx[n])
        prod = int(A_prod_idx[n])

        if prod == act and abs(float(A_signed[n])) == 1.0:
            continue

        tex = trails._get_tech_temporal_exchange(int(year), act, prod)
        if tex is None:
            continue

        td = TemporalDistribution(tex)
        pairs = list(td.iter_offsets_and_weights(debug=debug))
        if not pairs:
            multipliers_A[n] = 0.0
            continue

        m = 0.0
        for offset, w in pairs:
            m += float(w) * float(td.scale_factor(offset))
        multipliers_A[n] = float(m)

    return np.asarray(A_signed, dtype=np.float64) * multipliers_A


def _apply_temporal_scaling_to_B(
    trails: Trails,
    year: int,
    B_t,
    zero_biosphere: bool,
    collapse_bio_temporal_scaling: bool,
    debug: bool,
):
    B_signed = B_t.data.astype(np.float64, copy=False)
    if not collapse_bio_temporal_scaling or zero_biosphere:
        return B_signed

    B_act_idx, B_flow_idx = _ij_from_coords(B_t)
    multipliers = np.ones_like(B_signed, dtype=np.float64)

    for i in range(len(B_signed)):
        act = int(B_act_idx[i])
        flow = int(B_flow_idx[i])

        tex = trails._get_bio_temporal_exchange(int(year), act, flow)
        if tex is None:
            continue

        td = TemporalDistribution(tex)

        m = 0.0
        for offset, weight in td.iter_offsets_and_weights(debug=debug):
            m += float(weight) * float(td.scale_factor(int(offset)))

        multipliers[i] = m if (np.isfinite(m) and m > 0.0) else 1.0

    return B_signed * multipliers


def _make_bw_indices_rowcol(row_idx: np.ndarray, col_idx: np.ndarray) -> np.ndarray:
    idx = np.empty(len(row_idx), dtype=bwp.INDICES_DTYPE)
    idx["row"] = row_idx.astype(np.uint32, copy=False)
    idx["col"] = col_idx.astype(np.uint32, copy=False)
    return idx


def _warn_on_missing_metadata(label_for_matrix, meta_label, A_act_idx, A_prod_idx, B_flow_idx, act_meta, bio_meta):
    meta_act_indices = set(act_meta.keys())
    meta_bio_indices = set(bio_meta.keys())

    A_act_indices = set(A_act_idx.tolist())
    A_prod_indices = set(A_prod_idx.tolist())
    B_flow_indices = set(B_flow_idx.tolist())

    A_all_activities = A_act_indices | A_prod_indices

    missing_activities = A_all_activities - meta_act_indices
    if missing_activities:
        print(
            f"⚠️ WARNING: A[{label_for_matrix}] uses {len(missing_activities)} "
            f"activity indices not present in activity_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_activities))[:10]}"
            f"{' ...' if len(missing_activities) > 10 else ''}"
        )

    missing_flows = B_flow_indices - meta_bio_indices
    if missing_flows:
        print(
            f"⚠️ WARNING: B[{label_for_matrix}] uses {len(missing_flows)} "
            f"biosphere flow indices not present in biosphere_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_flows))[:10]}"
            f"{' ...' if len(missing_flows) > 10 else ''}"
        )


def _build_metadata_indices(act_meta, bio_meta):
    technosphere_indices: Dict[tuple, int] = {}
    for idx, meta in act_meta.items():
        key = (
            meta["name"],
            meta["reference product"],
            meta["unit"],
            meta["location"],
        )
        technosphere_indices[key] = idx

    biosphere_indices: Dict[tuple, int] = {}
    for idx, meta in bio_meta.items():
        key = (
            meta["name"],
            meta["compartment"],
            meta["subcompartment"],
            meta["unit"],
        )
        biosphere_indices[key] = idx

    return technosphere_indices, biosphere_indices


def _nearest_metadata_label_for_year(trails: Trails, year: int) -> str:
    """
    Pick the metadata scenario_label whose numeric year is closest to `year`.
    Used only as a fallback if the exact label is not present.
    """
    if not trails.activity_indices:
        raise ValueError("Trails.activity_indices is empty – no metadata available.")

    labels = list(trails.activity_indices.keys())
    years = np.array([int(lbl) for lbl in labels])
    idx = int(np.argmin(np.abs(years - year)))
    return labels[idx]


def build_datapackage_for_year_from_trails(
    trails: Trails,
    year: int,
    zero_biosphere: bool = False,
    collapse_bio_temporal_scaling: bool = False,
    collapse_tech_temporal_scaling: bool = False,
    debug: bool = False,
):
    """
    Build a bw_processing.Datapackage for a given calendar year
    directly from Trails.A / Trails.B.

    IMPORTANT:
    - Matrix indices (rows/cols) are taken as-is from Trails (which already
      reflect the indices from A_matrix_index.csv / B_matrix_index.csv).
    - Metadata (technosphere_indices, biosphere_indices) is taken from the
      *corresponding* index CSV for that scenario/year if available.
    - If no exact metadata exists for `year`, we fall back to the nearest
      available metadata year and warn.

    Returns
    -------
    dp : bwp.Datapackage
        Datapackage with 'technosphere_matrix' and 'biosphere_matrix'.
    technosphere_indices : dict[(name, ref, unit, loc), int]
    biosphere_indices : dict[(name, compartment, subcompartment, unit), int]
    uncertain_parameters : list[tuple[int, int]]
        Empty for now (no uncertainty).
    """
    # ------------------------------------------------------------------
    # 1) Map `year` to the matrix slice we actually use
    # ------------------------------------------------------------------
    label_for_matrix, t = _resolve_matrix_label(trails, int(year))

    # Extract 2D slices
    A_t = trails.A[t, :, :]  # sparse.COO (activities x products)
    B_t = trails.B[t, :, :]  # sparse.COO (activities x flows)

    # ------------------------------------------------------------------
    # 2) Choose metadata year: prefer exact, else nearest with warning
    # ------------------------------------------------------------------
    meta_label = _select_metadata_label(trails, label_for_matrix)

    act_meta = trails.activity_indices.get(meta_label, {})
    bio_meta = trails.biosphere_indices.get(meta_label, {})

    # ------------------------------------------------------------------
    # 3) Build technosphere entries from A_t
    # ------------------------------------------------------------------
    A_signed = _apply_temporal_scaling_to_A(
        trails=trails,
        year=int(year),
        A_t=A_t,
        debug=debug,
        collapse_tech_temporal_scaling=collapse_tech_temporal_scaling,
    )

    # Brightway convention via bw_processing:
    # - data must be non-negative
    # - flip marks entries that should be negated
    flip_A = (A_signed < 0)
    A_data = np.abs(A_signed)

    # Build indices array for bw_processing
    # (keep your existing indices construction, but ensure it uses A_coords)
    # Example (adapt to your existing code):
    # indices_A = np.vstack([A_coords[1], A_coords[2]]).T.astype(np.uint32)

    # Ensure dtype compatibility
    data_A = np.asarray(A_data, dtype=np.float64)
    flip_A = np.asarray(flip_A, dtype=bool)

    # ------------------------------------------------------------------
    # 4) Build biosphere entries from B_t
    # ------------------------------------------------------------------

    B_signed = _apply_temporal_scaling_to_B(
        trails=trails,
        year=int(year),
        B_t=B_t,
        zero_biosphere=zero_biosphere,
        collapse_bio_temporal_scaling=collapse_bio_temporal_scaling,
        debug=debug,
    )

    # bw_processing convention: biosphere also needs abs+flip (same as technosphere)
    flip_B = (B_signed < 0)
    data_B = np.abs(B_signed)

    if zero_biosphere:
        data_B = np.zeros_like(data_B)
        flip_B = np.zeros_like(flip_B, dtype=bool)

    # ------------------------------------------------------------------
    # 5) SAFETY CHECK: matrix indices compatible with metadata
    # ------------------------------------------------------------------
    A_act_idx, A_prod_idx = _ij_from_coords(A_t)
    B_act_idx, B_flow_idx = _ij_from_coords(B_t)

    indices_A = _make_bw_indices_rowcol(A_prod_idx, A_act_idx)
    indices_B = _make_bw_indices_rowcol(B_flow_idx, B_act_idx)

    _warn_on_missing_metadata(
        label_for_matrix=label_for_matrix,
        meta_label=meta_label,
        A_act_idx=A_act_idx,
        A_prod_idx=A_prod_idx,
        B_flow_idx=B_flow_idx,
        act_meta=act_meta,
        bio_meta=bio_meta,
    )


    # ------------------------------------------------------------------
    # 6) Create bw_processing datapackage
    # ------------------------------------------------------------------
    dp = bwp.create_datapackage()

    assert indices_A.shape == data_A.shape, (indices_A.shape, data_A.shape)
    assert indices_A.dtype == bwp.INDICES_DTYPE, indices_A.dtype
    assert np.all(data_A >= 0), "A data must be non-negative for flip convention"
    assert flip_A.shape == data_A.shape

    dp.add_persistent_vector(
        matrix="technosphere_matrix",
        indices_array=indices_A,
        data_array=data_A,
        flip_array=flip_A,
    )

    dp.add_persistent_vector(
        matrix="biosphere_matrix",
        indices_array=indices_B,
        data_array=data_B,
        flip_array=flip_B,
    )

    # ------------------------------------------------------------------
    # 7) Build technosphere_indices / biosphere_indices dictionaries
    #     (respecting CSV indices exactly)
    # ------------------------------------------------------------------
    technosphere_indices, biosphere_indices = _build_metadata_indices(act_meta, bio_meta)

    uncertain_parameters: List[Tuple[int, int]] = []  # none for now

    return dp, technosphere_indices, biosphere_indices, uncertain_parameters


def lca_static_simple(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    debug: bool = False,
) -> Dict[str, Any]:
    """
    Plain static LCA for a single FU activity in a single year:
      demand = {fu_act_idx: amount}
      solves full upstream supply chain
      includes that activity's direct biosphere and upstream biosphere
    """
    # Build datapackage for that year with biosphere enabled
    dp, _, bio_idx, _ = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        collapse_bio_temporal_scaling=True,
        collapse_tech_temporal_scaling=True,
        debug=debug,
    )

    lca_obj = bc.LCA(demand={int(fu_act_idx): float(amount)}, data_objs=[dp])
    lca_obj.lci()

    # Build characterization matrix in BW order (static)
    bw_bio_map = lca_obj.dicts.biosphere  # flow_id -> row position
    biosphere_matrix_dict = {int(flow_id): int(pos) for flow_id, pos in bw_bio_map.items()}

    biosphere_dict_simple = {
        (name, comp, subcomp): int(flow_id)
        for (name, comp, subcomp, unit), flow_id in bio_idx.items()
    }

    C = fill_characterization_factors_matrices(
        methods=methods,
        biosphere_matrix_dict=biosphere_matrix_dict,
        biosphere_dict=biosphere_dict_simple,
        debug=debug,
    )

    inv_vec = np.asarray(lca_obj.inventory.sum(axis=1)).reshape((-1, 1))
    score = float(np.sum(C.dot(inv_vec)))

    return {
        "results_by_solve_year": {
            int(year): {
                "fu": {int(fu_act_idx): float(amount)},
                "lca": lca_obj,
                "scores": score,
                "scores_per_root": {int(fu_act_idx): score},  # trivial "root"
            }
        },
        "results_by_impact_year": {
            int(year): {
                "scores": score,
                "scores_per_root": {int(fu_act_idx): score},
            }
        },
    }


def _normalize_root(fu0: int, root: int, legacy_root: int) -> int:
    root = int(root)
    return fu0 if root == legacy_root else root


def _run_temporal_traversal(
    trails: Trails,
    y0: int,
    fu0: int,
    amt0: float,
    max_depth: int,
    min_amount: float,
    show_progress: bool,
    debug: bool,
):
    (
        frontier,
        provenance,
        injected_supply_by_year_act,
        injected_supply_prov_by_year_act,
    ) = trails.temporal_traversal(
        start_year=y0,
        start_act_idx=fu0,
        amount=amt0,
        max_depth=int(max_depth),
        min_amount=float(min_amount),
        return_provenance=True,
        show_progress=bool(show_progress),
        use_temporal_distributions=True,
        debug=debug,
    )

    if injected_supply_by_year_act is None:
        injected_supply_by_year_act = {}
    if injected_supply_prov_by_year_act is None:
        injected_supply_prov_by_year_act = {}

    return frontier, provenance, injected_supply_by_year_act, injected_supply_prov_by_year_act


def _apply_fu_direct_injection(
    injected_supply_by_year_act,
    injected_supply_prov_by_year_act,
    y0: int,
    fu0: int,
    amt0: float,
    legacy_root: int,
):
    injected_supply_by_year_act[(y0, fu0)] = float(injected_supply_by_year_act.get((y0, fu0), 0.0)) + amt0

    injected_supply_prov_by_year_act.setdefault((y0, fu0), {})
    injected_supply_prov_by_year_act[(y0, fu0)][fu0] = (
        float(injected_supply_prov_by_year_act[(y0, fu0)].get(fu0, 0.0)) + amt0
    )

    if legacy_root in injected_supply_prov_by_year_act[(y0, fu0)]:
        injected_supply_prov_by_year_act[(y0, fu0)][fu0] = (
            float(injected_supply_prov_by_year_act[(y0, fu0)].get(fu0, 0.0))
            + float(injected_supply_prov_by_year_act[(y0, fu0)].pop(legacy_root))
        )


def _build_rooted_frontier(frontier, provenance, fu0: int, root_tol: float):
    rooted_frontier = defaultdict(float)
    roots_seen = set()

    def _iter_path_nodes(path):
        if path is None:
            return
        if isinstance(path, (int, np.integer)):
            yield int(path)
            return
        try:
            for node in path:
                if isinstance(node, (int, np.integer)):
                    yield int(node)
                elif isinstance(node, (tuple, list)) and len(node) >= 2:
                    yield int(node[1])
        except TypeError:
            return

    def _root_from_path(path, fu_act: int, fallback_root: int) -> int:
        if isinstance(path, (int, np.integer)):
            r = int(path)
            return r if r != fu_act else int(fallback_root)

        for act in _iter_path_nodes(path):
            if act != fu_act:
                return int(act)

        return int(fallback_root)

    for (year, act), total_amt in frontier.items():
        year = int(year)
        act = int(act)
        total_amt = float(total_amt)

        prov = provenance.get((year, act), {})

        if not prov:
            rooted_frontier[(year, act, act)] += total_amt
            roots_seen.add(act)
            continue

        prov_sum = 0.0
        for path, amt_share in prov.items():
            amt_share = float(amt_share)
            prov_sum += amt_share

            root_act = _root_from_path(path, fu_act=fu0, fallback_root=act)
            rooted_frontier[(year, act, int(root_act))] += amt_share
            roots_seen.add(int(root_act))

        residual = total_amt - prov_sum
        if abs(residual) > root_tol:
            rooted_frontier[(year, act, act)] += residual
            roots_seen.add(act)

    return rooted_frontier, roots_seen


def _extend_roots_with_injected_supply(roots_seen, injected_supply_prov_by_year_act, normalize_root):
    for (_y, _a), roots_map in injected_supply_prov_by_year_act.items():
        for r in roots_map.keys():
            roots_seen.add(normalize_root(int(r)))


def _build_root_demand_vectors(trails: Trails, rooted_frontier, roots_seen, normalize_root):
    n_activities = int(trails.A.shape[1])
    dtype = trails.value_dtype
    f_by_year_by_root = {r: {} for r in roots_seen}

    for (year, act, root_act), amt in rooted_frontier.items():
        y = int(year)
        a = int(act)
        r = normalize_root(int(root_act))

        if y not in f_by_year_by_root[r]:
            f_by_year_by_root[r][y] = np.zeros(n_activities, dtype=dtype)
        f_by_year_by_root[r][y][a] += dtype(amt)

    return f_by_year_by_root


def _get_datapackage(dp_cache, trails: Trails, year: int, zero_bio: bool, debug: bool):
    dp_key = (year, zero_bio)
    if dp_key not in dp_cache:
        dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
            trails=trails,
            year=year,
            zero_biosphere=zero_bio,
            debug=debug,
        )
        dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
    return dp_cache[dp_key]


def _build_demand_by_first_level_child(f_by_year_by_root, solve_year, min_amount):
    demand_by_first_level_child: Dict[int, Dict[int, float]] = {}
    for root_idx, by_year in f_by_year_by_root.items():
        vec = by_year.get(solve_year)
        if vec is None:
            continue
        arr_r = np.asarray(vec)
        nz_r = np.where(np.abs(arr_r) > float(min_amount))[0]
        if nz_r.size == 0:
            continue
        demand_by_first_level_child[int(root_idx)] = {int(i): float(arr_r[i]) for i in nz_r}
    return demand_by_first_level_child


def _assert_rooted_closure(arr, demand_by_first_level_child, solve_year):
    summed = np.zeros_like(arr, dtype=float)
    for dct in demand_by_first_level_child.values():
        for i, v in dct.items():
            summed[i] += v
    abs_err = float(np.max(np.abs(summed - arr))) if arr.size else 0.0
    scale = float(np.max(np.abs(arr))) if arr.size else 0.0
    tol = max(1e-9, 1e-6 * max(1.0, scale))
    if abs_err > tol:
        raise ValueError(
            f"Rooted demand does not sum to total demand in year {solve_year} (abs_err={abs_err:g}, tol={tol:g})."
        )


def _extract_supply(lca_obj, n_acts: int, min_amount: float) -> Dict[int, float]:
    supply_total: Dict[int, float] = {}
    for act_idx in range(n_acts):
        try:
            act_pos = lca_obj.dicts.activity[act_idx]
        except KeyError:
            continue
        s = float(lca_obj.supply_array[act_pos])
        if abs(s) > float(min_amount):
            supply_total[int(act_idx)] = s
    return supply_total


def _build_injected_supply(injected_supply_by_year_act, solve_year, min_amount):
    injected_supply: Dict[int, float] = {}
    for (y, a), v in injected_supply_by_year_act.items():
        if int(y) != solve_year:
            continue
        v = float(v)
        if abs(v) <= float(min_amount):
            continue
        injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v
    return injected_supply


def _build_injected_supply_by_root(
    injected_supply_prov_by_year_act,
    solve_year,
    min_amount,
    normalize_root,
):
    injected_supply_by_first_level_child: Dict[int, Dict[int, float]] = {}
    for (y, act), roots_map in injected_supply_prov_by_year_act.items():
        if int(y) != solve_year:
            continue
        for r, share in roots_map.items():
            r = normalize_root(int(r))
            share = float(share)
            if abs(share) <= float(min_amount):
                continue
            injected_supply_by_first_level_child.setdefault(r, {})
            injected_supply_by_first_level_child[r][int(act)] = (
                injected_supply_by_first_level_child[r].get(int(act), 0.0) + share
            )
    return injected_supply_by_first_level_child


def _build_global_biosphere_dict_simple(trails: Trails) -> Dict[tuple, int]:
    """
    Build a stable (name, compartment, subcompartment) -> flow_id mapping.

    Supports two possible internal formats for trails.biosphere_indices[label]:
      A) {flow_id: {"name": ..., "compartment": ..., "subcompartment": ..., "unit": ...}}
      B) {(name, compartment, subcompartment, unit): flow_id}
    """
    out: Dict[tuple, int] = {}

    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        if not meta:
            continue

        # Peek one item to infer structure
        k0 = next(iter(meta.keys()))

        # Case A: keyed by int flow_id -> meta dict
        if isinstance(k0, (int, np.integer)):
            for flow_id, md in meta.items():
                if not isinstance(md, dict):
                    continue
                name = md.get("name")
                comp = md.get("compartment")
                subcomp = md.get("subcompartment")
                if name is None or comp is None or subcomp is None:
                    continue
                key = (name, comp, subcomp)
                out.setdefault(key, int(flow_id))

        # Case B: keyed by tuple -> flow_id
        elif isinstance(k0, tuple) and len(k0) >= 3:
            for k, flow_id in meta.items():
                # k is (name, comp, subcomp, unit) or similar
                name = k[0]
                comp = k[1] if len(k) > 1 else None
                subcomp = k[2] if len(k) > 2 else None
                if name is None or comp is None or subcomp is None:
                    continue
                out.setdefault((name, comp, subcomp), int(flow_id))

        else:
            raise TypeError(
                f"Unsupported biosphere_indices structure for label {_label}: "
                f"key type={type(k0)} example={k0}"
            )

    return out

def _build_flowkey_to_flowid(trails: Trails) -> Dict[tuple, int]:
    """
    Build (name, compartment, subcompartment) -> flow_id mapping,
    robust across labels.
    """
    out: Dict[tuple, int] = {}

    for _label, meta in getattr(trails, "biosphere_indices", {}).items():
        if not meta:
            continue

        k0 = next(iter(meta.keys()))

        # Expected: {flow_id: {"name":..., "compartment":..., "subcompartment":...}}
        if isinstance(k0, (int, np.integer)):
            for flow_id, md in meta.items():
                if not isinstance(md, dict):
                    continue
                name = md.get("name")
                comp = md.get("compartment")
                sub = md.get("subcompartment")
                if name is None or comp is None or sub is None:
                    continue
                out.setdefault((name, comp, sub), int(flow_id))
        else:
            raise TypeError(
                f"Unexpected biosphere_indices structure for label {_label}: "
                f"expected int keys, got {type(k0)}"
            )

    return out

def _build_cf_vector_flowid_space(
    trails: Trails,
    methods: List[str],
    ei_version: str,
    char_cache: dict,
    debug: bool = False,
) -> np.ndarray:
    """
    Return cf vector of length n_flows in Trails flow-id space:
      cf[flow_id] = sum_m CF_m(flow_key_of_flow_id)

    Assumes LCIA method exchanges keyed by (name, compartment, subcompartment).
    """
    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        return np.zeros(0, dtype=np.float64)

    cache_key = ("cf_vector_flowid_space", ei_version, tuple(methods))
    if cache_key in char_cache:
        return char_cache[cache_key]

    flowkey_to_flowid = _build_flowkey_to_flowid(trails)

    methods_dict = get_lcia_methods(methods=methods, ei_version=ei_version)

    cf = np.zeros(n_flows, dtype=np.float64)

    # Sum CFs across methods (your code already supports multiple methods)
    for mname, exc in methods_dict.items():
        for flow_key, val in exc.items():
            # flow_key is (name, comp, subcomp)
            fid = flowkey_to_flowid.get(flow_key)
            if fid is None:
                continue
            cf[fid] += float(val)

    char_cache[cache_key] = cf
    return cf


def _characterize_impact_years(
    trails: Trails,
    inventory_total_by_impact_year,
    inventory_by_root_by_impact_year,
    dp_cache,
    char_cache,
    methods,
    min_amount,
    normalize_root,
    debug,
    ei_version="3.11",
):
    results_by_impact_year: Dict[int, Dict[str, Any]] = {}

    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        return results_by_impact_year

    cf = _build_cf_vector_flowid_space(
        trails=trails,
        methods=methods,
        ei_version=ei_version,
        char_cache=char_cache,
        debug=debug,
    )

    impact_years = sorted(set(inventory_total_by_impact_year.keys()))
    impact_iter = tqdm(impact_years, desc="Temporal LCA: impact years", unit="year")

    for impact_year in impact_iter:
        impact_year = int(impact_year)

        inv_total = inventory_total_by_impact_year.get(
            impact_year, np.zeros(n_flows, dtype=trails.value_dtype)
        )
        total_score = float(np.dot(cf, inv_total.astype(np.float64, copy=False)))

        scores_by_first_level_child: Dict[int, float] = {}
        for root_idx, inv_map in inventory_by_root_by_impact_year.items():
            root_norm = normalize_root(int(root_idx))
            inv_root = inv_map.get(impact_year, np.zeros(n_flows, dtype=trails.value_dtype))
            s = float(np.dot(cf, inv_root.astype(np.float64, copy=False)))
            if abs(s) <= float(min_amount):
                continue
            scores_by_first_level_child[root_norm] = scores_by_first_level_child.get(root_norm, 0.0) + s

        results_by_impact_year[impact_year] = {
            "scores": total_score,
            "scores_by_first_level_child": scores_by_first_level_child,
        }

    return results_by_impact_year


def lca(
    trails: Trails,
    start_year: int,
    start_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    max_depth: int = 2,
    min_amount: float = 1e-18,
    show_progress: bool = True,
    debug: bool = False,
    return_provenance: bool = False,
    use_temporal_distributions: bool = True,
) -> Dict[str, Any]:
    """
    Option A output:
      - results_by_solve_year: diagnostic per-year solves
      - results_by_impact_year: time series of impacts booked in impact years (for plotting)

    Naming conventions:
      - fu_demand: the full demand vector solved in a given solve_year (sparse dict form)
      - demand_by_first_level_child: technosphere demand split by first-level child under the FU
      - injected_supply: extra "supply pulses" to book direct biosphere without re-solving technosphere
      - injected_supply_by_first_level_child: same, but attributed to a root bucket

    IMPORTANT BEHAVIOR (FIX):
      - FU-direct biosphere is attributed to the *FU activity index* (fu0),
        not to a synthetic bucket like -1.
      - Any legacy occurrences of -1 are normalized into fu0 during bookkeeping.
    """

    import numpy as np
    import bw2calc as bc
    from collections import defaultdict

    # Legacy sentinel (may still appear in some provenance artifacts)
    LEGACY_FU_DIRECT_ROOT = -1

    ROOT_CLOSURE_TOL = 1e-12

    fu0 = int(start_act_idx)
    y0 = int(start_year)
    amt0 = float(amount)

    # All FU-direct biosphere will be booked under this root:
    FU_DIRECT_ROOT = fu0

    def _normalize_root_local(r: int) -> int:
        return _normalize_root(fu0, r, LEGACY_FU_DIRECT_ROOT)

    # ---------------------------------------------------------
    # STATIC SHORTCUT
    # ---------------------------------------------------------
    if not use_temporal_distributions:
        return lca_static_simple(
            trails=trails,
            year=y0,
            fu_act_idx=fu0,
            methods=methods,
            amount=amt0,
            debug=debug,
        )

    # -----------------------------
    # 1) Temporal traversal (+ provenance for rooting)
    # -----------------------------
    need_roots = True

    if return_provenance or need_roots:
        (
            frontier,
            provenance,
            injected_supply_by_year_act,
            injected_supply_prov_by_year_act,
        ) = _run_temporal_traversal(
            trails=trails,
            y0=y0,
            fu0=fu0,
            amt0=amt0,
            max_depth=max_depth,
            min_amount=min_amount,
            show_progress=show_progress,
            debug=debug,
        )
    else:
        frontier, injected_supply_by_year_act = trails.temporal_traversal(
            start_year=y0,
            start_act_idx=fu0,
            amount=amt0,
            max_depth=int(max_depth),
            min_amount=float(min_amount),
            return_provenance=False,
            show_progress=bool(show_progress),
            use_temporal_distributions=True,
            debug=debug,
        )
        provenance = {}
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        injected_supply_prov_by_year_act = {}

    # -----------------------------
    # 1a) Inject FU self-supply for FU-direct biosphere only
    # -----------------------------
    # We keep the injected supply pulse on (y0, fu0) ...
    _apply_fu_direct_injection(
        injected_supply_by_year_act=injected_supply_by_year_act,
        injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
        y0=y0,
        fu0=fu0,
        amt0=amt0,
        legacy_root=LEGACY_FU_DIRECT_ROOT,
    )

    # -----------------------------
    # 1b) Frontier -> per-year demand vectors
    # -----------------------------
    f_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(f_by_year.keys())

    # -----------------------------
    # 1c) Rooted frontier using provenance paths
    # -----------------------------
    rooted_frontier, roots_seen = _build_rooted_frontier(
        frontier=frontier,
        provenance=provenance,
        fu0=fu0,
        root_tol=ROOT_CLOSURE_TOL,
    )
    _extend_roots_with_injected_supply(
        roots_seen=roots_seen,
        injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
        normalize_root=_normalize_root_local,
    )

    roots_seen = sorted(roots_seen)

    # Convert rooted frontier to per-year demand vectors by root (technosphere only)
    f_by_year_by_root = _build_root_demand_vectors(
        trails=trails,
        rooted_frontier=rooted_frontier,
        roots_seen=roots_seen,
        normalize_root=_normalize_root_local,
    )

    # -----------------------------
    # Results + caches
    # -----------------------------
    results_by_solve_year: Dict[int, Dict[str, Any]] = {}
    results_by_impact_year: Dict[int, Dict[str, Any]] = {}

    dp_cache: Dict[tuple, Any] = {}
    char_cache: Dict[tuple, Any] = {}

    inventory_total_by_impact_year: Dict[int, np.ndarray] = {}
    inventory_by_root_by_impact_year: Dict[int, Dict[int, np.ndarray]] = defaultdict(dict)

    # -----------------------------
    # 2) Solve-year loop
    # -----------------------------
    solve_iter = candidate_years
    if show_progress:
        solve_iter = tqdm(candidate_years, desc="Temporal LCA: solve years", unit="year")
    for solve_year in solve_iter:
        solve_year = int(solve_year)
        arr = np.asarray(f_by_year[solve_year])

        nz_idx = np.where(np.abs(arr) > float(min_amount))[0]
        if nz_idx.size == 0:
            continue

        fu_demand = {int(i): float(arr[i]) for i in nz_idx}

        zero_bio = True
        dp, _, _, _ = _get_datapackage(
            dp_cache=dp_cache,
            trails=trails,
            year=solve_year,
            zero_bio=zero_bio,
            debug=debug,
        )

        lca_total = bc.LCA(demand=fu_demand, data_objs=[dp])
        lca_total.lci()

        # demand split by first-level child roots for this solve_year
        demand_by_first_level_child = _build_demand_by_first_level_child(
            f_by_year_by_root=f_by_year_by_root,
            solve_year=solve_year,
            min_amount=min_amount,
        )
        _assert_rooted_closure(arr, demand_by_first_level_child, solve_year)

        # supply extraction (total)
        n_acts = int(trails.B.shape[1])
        supply_total = _extract_supply(lca_total, n_acts, min_amount)

        # (1) solved-supply contribution
        trails.accumulate_temporalized_biosphere_inventory(
            base_year=solve_year,
            supply_by_activity=supply_total,
            inventory_by_year=inventory_total_by_impact_year,
            min_amount=float(min_amount),
            use_temporal_distributions=True,
            debug=debug,
        )

        # (2) injected supply for this solve_year
        injected_supply = _build_injected_supply(
            injected_supply_by_year_act=injected_supply_by_year_act,
            solve_year=solve_year,
            min_amount=min_amount,
        )

        if injected_supply:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=injected_supply,
                inventory_by_year=inventory_total_by_impact_year,
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

        # Build injected supply by root, normalizing legacy sentinel -> fu0
        injected_supply_by_first_level_child = _build_injected_supply_by_root(
            injected_supply_prov_by_year_act=injected_supply_prov_by_year_act,
            solve_year=solve_year,
            min_amount=min_amount,
            normalize_root=_normalize_root_local,
        )

        # Book FU-direct inventory under FU activity root (fu0)
        fu_direct_injected = injected_supply_by_first_level_child.get(FU_DIRECT_ROOT, {})
        if fu_direct_injected:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=fu_direct_injected,
                inventory_by_year=inventory_by_root_by_impact_year[FU_DIRECT_ROOT],
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

        # Per-root inventories (upstream technosphere and any injected shares attributed to that root)
        for root_idx, root_demand in demand_by_first_level_child.items():
            lca_root = bc.LCA(demand=root_demand, data_objs=[dp])
            lca_root.lci()

            supply_root: Dict[int, float] = {}
            for act_idx in range(n_acts):
                try:
                    act_pos = lca_root.dicts.activity[act_idx]
                except KeyError:
                    continue
                s = float(lca_root.supply_array[act_pos])
                if abs(s) > float(min_amount):
                    supply_root[int(act_idx)] = s

            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=supply_root,
                inventory_by_year=inventory_by_root_by_impact_year[int(root_idx)],
                min_amount=float(min_amount),
                use_temporal_distributions=True,
                debug=debug,
            )

            root_injected = injected_supply_by_first_level_child.get(int(root_idx), {})
            if root_injected:
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=solve_year,
                    supply_by_activity=root_injected,
                    inventory_by_year=inventory_by_root_by_impact_year[int(root_idx)],
                    min_amount=float(min_amount),
                    use_temporal_distributions=True,
                    debug=debug,
                )

        # For diagnostics: include all roots that exist either via rooted demand or injected supply (incl. FU_DIRECT_ROOT)
        all_diag_roots = set(int(k) for k in demand_by_first_level_child.keys()) | set(
            int(k) for k in injected_supply_by_first_level_child.keys()
        )

        results_by_solve_year[solve_year] = {
            "fu_demand": fu_demand,
            "n_nonzero_demand": int(len(fu_demand)),
            "sum_abs_demand": float(np.sum(np.abs(arr[nz_idx]))),
            "max_abs_demand": float(np.max(np.abs(arr[nz_idx]))),

            "n_injected_supply": int(len(injected_supply)),
            "injected_supply": injected_supply,  # keep if you want detail

            "roots": sorted(int(k) for k in all_diag_roots),
            "demand_by_first_level_child": demand_by_first_level_child,
            "injected_supply_by_first_level_child": injected_supply_by_first_level_child,

            # Optional but often useful:
            "n_supply_total": int(len(supply_total)),
            "supply_total": supply_total if debug else None,  # can be huge; gate it

            "lca": lca_total,
        }

    # -----------------------------
    # 3) Impact-year characterization
    # -----------------------------
    results_by_impact_year = _characterize_impact_years(
        trails=trails,
        inventory_total_by_impact_year=inventory_total_by_impact_year,
        inventory_by_root_by_impact_year=inventory_by_root_by_impact_year,
        dp_cache=dp_cache,
        char_cache=char_cache,
        methods=methods,
        min_amount=min_amount,
        normalize_root=_normalize_root_local,
        debug=debug,
    )

    if not results_by_impact_year:
        out = {"results_by_solve_year": results_by_solve_year, "results_by_impact_year": {}}
        return (out, provenance) if return_provenance else out

    out = {"results_by_solve_year": results_by_solve_year, "results_by_impact_year": results_by_impact_year}
    return (out, provenance) if return_provenance else out

