import numpy as np
import bw_processing as bwp
from typing import Dict, Tuple, List, Any
from collections import defaultdict

import bw2calc as bc
import pyprind

from .trails import Trails
from .lcia import fill_characterization_factors_matrices
from .temporal_distributions import TemporalDistribution

import logging
logger = logging.getLogger(__name__)


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
    collapse_tech_temporal_scaling: bool = False
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
    label_for_matrix = str(year)
    if label_for_matrix not in trails.scenario_index:
        # If you used interpolation, this might be rare; still, handle safely.
        years = np.array([int(lbl) for lbl in trails.scenario_labels])
        idx = int(np.argmin(np.abs(years - year)))
        label_for_matrix = trails.scenario_labels[idx]
        print(
            f"⚠️ Matrix slice for year {year} not found in A/B; "
            f"using nearest available matrix year {label_for_matrix}."
        )

    t = trails.scenario_index[label_for_matrix]

    # Extract 2D slices
    A_t = trails.A[t, :, :]  # sparse.COO (activities x products)
    B_t = trails.B[t, :, :]  # sparse.COO (activities x flows)

    def _ij_from_coords(X_t):
        """
        Return (i_idx, j_idx) from a sparse.COO with either:
          - 3D coords: [t, i, j]
          - 2D coords: [i, j]
        """
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

    # ------------------------------------------------------------------
    # 2) Choose metadata year: prefer exact, else nearest with warning
    # ------------------------------------------------------------------
    if label_for_matrix in trails.activity_indices:
        meta_label = label_for_matrix
    else:
        meta_label = _nearest_metadata_label_for_year(trails, int(label_for_matrix))

    act_meta = trails.activity_indices.get(meta_label, {})
    bio_meta = trails.biosphere_indices.get(meta_label, {})

    # ------------------------------------------------------------------
    # 3) Build technosphere entries from A_t
    # ------------------------------------------------------------------
    A_signed = A_t.data

    # --- Optional: collapse temporal scaling metadata into the static A coefficients ---
    if collapse_tech_temporal_scaling:
        A_act_idx, A_prod_idx = _ij_from_coords(A_t)

        multipliers_A = np.ones_like(A_signed, dtype=np.float64)

        # Use template year for metadata lookup (same policy as Trails)
        template_year = trails._map_year_to_template_year(int(year))

        for n in range(len(A_signed)):
            act = int(A_act_idx[n])
            prod = int(A_prod_idx[n])

            # Skip diagonal production exchange if present (same logic as traversal)
            # NOTE: keep consistent with your A convention; adjust if needed.
            if prod == act and abs(float(A_signed[n])) == 1.0:
                continue

            tex = trails._get_tech_temporal_exchange(int(year), act, prod)
            if tex is None:
                continue

            td = TemporalDistribution(tex)
            pairs = list(td.iter_offsets_and_weights())
            if not pairs:
                multipliers_A[n] = 0.0
                continue

            # Expected scaling under the temporal distribution
            m = 0.0
            for offset, w in pairs:
                m += float(w) * float(td.scale_factor(offset))
            multipliers_A[n] = float(m)

        A_signed = np.asarray(A_signed, dtype=np.float64) * multipliers_A

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

    B_signed = B_t.data.astype(np.float64, copy=False)

    # --- Optional: collapse temporal scaling metadata into the static B coefficients ---
    if collapse_bio_temporal_scaling and (not zero_biosphere):
        # We need (act, flow) coordinates aligned with B_signed
        B_act_idx, B_flow_idx = _ij_from_coords(B_t)

        multipliers = np.ones_like(B_signed, dtype=np.float64)

        # Use template year for metadata lookup, consistent with Trails design
        template_year = trails._map_year_to_template_year(int(year))

        for i in range(len(B_signed)):
            act = int(B_act_idx[i])
            flow = int(B_flow_idx[i])

            tex = trails._get_bio_temporal_exchange(int(year), act, flow)
            if tex is None:
                continue

            td = TemporalDistribution(tex)

            m = 0.0
            for offset, weight in td.iter_offsets_and_weights():
                m += float(weight) * float(td.scale_factor(int(offset)))

            # If something goes wrong, fall back to 1.0 instead of zeroing the exchange
            multipliers[i] = m if (np.isfinite(m) and m > 0.0) else 1.0

        B_signed = B_signed * multipliers

    # bw_processing convention: biosphere also needs abs+flip (same as technosphere)
    flip_B = (B_signed < 0)
    data_B = np.abs(B_signed)

    if zero_biosphere:
        data_B = np.zeros_like(data_B)
        flip_B = np.zeros_like(flip_B, dtype=bool)

    # ------------------------------------------------------------------
    # 5) SAFETY CHECK: matrix indices compatible with metadata
    # ------------------------------------------------------------------
    meta_act_indices = set(act_meta.keys())
    meta_bio_indices = set(bio_meta.keys())

    # ---- Extract index vectors from sparse coords (COO: coords[0]=t, coords[1]=i, coords[2]=j) ----


    A_act_idx, A_prod_idx = _ij_from_coords(A_t)  # act, prod
    B_act_idx, B_flow_idx = _ij_from_coords(B_t)  # act, flow

    A_act_indices = set(A_act_idx.tolist())
    A_prod_indices = set(A_prod_idx.tolist())
    B_flow_indices = set(B_flow_idx.tolist())

    # All activity indices used anywhere in A (row/col in your internal orientation)
    A_all_activities = A_act_indices | A_prod_indices

    def _make_bw_indices_rowcol(row_idx: np.ndarray, col_idx: np.ndarray) -> np.ndarray:
        idx = np.empty(len(row_idx), dtype=bwp.INDICES_DTYPE)
        idx["row"] = row_idx.astype(np.uint32, copy=False)
        idx["col"] = col_idx.astype(np.uint32, copy=False)
        return idx

    # Correct for bw2calc:
    # Technosphere is (product rows, activity cols)
    indices_A = _make_bw_indices_rowcol(A_prod_idx, A_act_idx)

    # Biosphere is (flow rows, activity cols)
    indices_B = _make_bw_indices_rowcol(B_flow_idx, B_act_idx)

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
    technosphere_indices: Dict[tuple, int] = {}
    for idx, meta in act_meta.items():
        key = (
            meta["name"],
            meta["reference product"],
            meta["unit"],
            meta["location"],
        )
        # HERE: we keep `idx` as given in A_matrix_index.csv
        technosphere_indices[key] = idx

    biosphere_indices: Dict[tuple, int] = {}
    for idx, meta in bio_meta.items():
        key = (
            meta["name"],
            meta["compartment"],
            meta["subcompartment"],
            meta["unit"],
        )
        # HERE: we keep `idx` as given in B_matrix_index.csv
        biosphere_indices[key] = idx

    uncertain_parameters: List[Tuple[int, int]] = []  # none for now

    return dp, technosphere_indices, biosphere_indices, uncertain_parameters


def lca_static_simple(
    trails: Trails,
    year: int,
    fu_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
) -> Dict[str, Any]:
    """
    Plain static LCA for a single FU activity in a single year:
      demand = {fu_act_idx: amount}
      solves full upstream supply chain
      includes that activity's direct biosphere and upstream biosphere
    """
    # Build datapackage for that year with biosphere enabled
    dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
        trails=trails,
        year=int(year),
        zero_biosphere=False,
        collapse_bio_temporal_scaling=True,
        collapse_tech_temporal_scaling=True,
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

def lca(
    trails: Trails,
    start_year: int,
    start_act_idx: int,
    methods: List[str],
    amount: float = 1.0,
    max_depth: int = 2,
    min_amount: float = 1e-16,
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

    def _normalize_root(r: int) -> int:
        """Collapse any legacy sentinel to FU activity index."""
        r = int(r)
        return fu0 if r == LEGACY_FU_DIRECT_ROOT else r

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
        ) = trails.temporal_traversal(
            start_year=y0,
            start_act_idx=fu0,
            amount=amt0,
            max_depth=int(max_depth),
            min_amount=float(min_amount),
            return_provenance=True,
            show_progress=bool(show_progress),
            use_temporal_distributions=True,
        )
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        if injected_supply_prov_by_year_act is None:
            injected_supply_prov_by_year_act = {}
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
        )
        provenance = {}
        if injected_supply_by_year_act is None:
            injected_supply_by_year_act = {}
        injected_supply_prov_by_year_act = {}

    # -----------------------------
    # 1a) Inject FU self-supply for FU-direct biosphere only
    # -----------------------------
    # We keep the injected supply pulse on (y0, fu0) ...
    injected_supply_by_year_act[(y0, fu0)] = float(injected_supply_by_year_act.get((y0, fu0), 0.0)) + amt0

    # ... but provenance root MUST be fu0 (not -1).
    injected_supply_prov_by_year_act.setdefault((y0, fu0), {})
    injected_supply_prov_by_year_act[(y0, fu0)][FU_DIRECT_ROOT] = (
        float(injected_supply_prov_by_year_act[(y0, fu0)].get(FU_DIRECT_ROOT, 0.0)) + amt0
    )

    # If any earlier logic injected the legacy sentinel, normalize it immediately.
    if LEGACY_FU_DIRECT_ROOT in injected_supply_prov_by_year_act[(y0, fu0)]:
        injected_supply_prov_by_year_act[(y0, fu0)][FU_DIRECT_ROOT] = (
            float(injected_supply_prov_by_year_act[(y0, fu0)].get(FU_DIRECT_ROOT, 0.0))
            + float(injected_supply_prov_by_year_act[(y0, fu0)].pop(LEGACY_FU_DIRECT_ROOT))
        )

    # -----------------------------
    # 1b) Frontier -> per-year demand vectors
    # -----------------------------
    f_by_year = trails.frontier_to_demand_vectors(frontier)
    candidate_years = sorted(f_by_year.keys())

    # -----------------------------
    # 1c) Rooted frontier using provenance paths
    # -----------------------------
    rooted_frontier = defaultdict(float)
    roots_seen = set()  # <-- do not pre-seed with -1; FU_DIRECT_ROOT is a real activity idx now

    def _iter_path_nodes(path):
        """
        Normalize provenance keys into an iterator of (year, act) or act.
        Supported:
          - int root
          - tuple/list path whose elements are either:
              * (year, act) tuples/lists
              * act ints
        """
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
        """
        FIRST-LEVEL CHILD root:
          return the first act in the path that is NOT the FU activity.
        """
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
        if abs(residual) > ROOT_CLOSURE_TOL:
            rooted_frontier[(year, act, act)] += residual
            roots_seen.add(act)

    # Include roots that appear only in injected supply provenance
    for (_y, _a), roots_map in injected_supply_prov_by_year_act.items():
        for r in roots_map.keys():
            roots_seen.add(_normalize_root(int(r)))

    roots_seen = sorted(roots_seen)

    # Convert rooted frontier to per-year demand vectors by root (technosphere only)
    n_activities = int(trails.A.shape[1])
    dtype = trails.value_dtype
    f_by_year_by_root = {r: {} for r in roots_seen}

    for (year, act, root_act), amt in rooted_frontier.items():
        y = int(year)
        a = int(act)
        r = int(root_act)
        # Technosphere rooting should never use the legacy FU-direct sentinel;
        # but normalize defensively anyway.
        r = _normalize_root(r)

        if y not in f_by_year_by_root[r]:
            f_by_year_by_root[r][y] = np.zeros(n_activities, dtype=dtype)
        f_by_year_by_root[r][y][a] += dtype(amt)

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
    for solve_year in candidate_years:
        solve_year = int(solve_year)
        arr = np.asarray(f_by_year[solve_year])

        nz_idx = np.where(np.abs(arr) > float(min_amount))[0]
        if nz_idx.size == 0:
            continue

        fu_demand = {int(i): float(arr[i]) for i in nz_idx}

        zero_bio = True
        dp_key = (solve_year, zero_bio)
        if dp_key not in dp_cache:
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=solve_year,
                zero_biosphere=zero_bio,
            )
            dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

        lca_total = bc.LCA(demand=fu_demand, data_objs=[dp])
        lca_total.lci()

        # demand split by first-level child roots for this solve_year
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

        # closure check
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

        # supply extraction (total)
        n_acts = int(trails.B.shape[1])
        supply_total: Dict[int, float] = {}
        for act_idx in range(n_acts):
            try:
                act_pos = lca_total.dicts.activity[act_idx]
            except KeyError:
                continue
            s = float(lca_total.supply_array[act_pos])
            if abs(s) > float(min_amount):
                supply_total[int(act_idx)] = s

        # (1) solved-supply contribution
        trails.accumulate_temporalized_biosphere_inventory(
            base_year=solve_year,
            supply_by_activity=supply_total,
            inventory_by_year=inventory_total_by_impact_year,
            min_amount=float(min_amount),
            use_temporal_distributions=True,
        )

        # (2) injected supply for this solve_year
        injected_supply: Dict[int, float] = {}
        for (y, a), v in injected_supply_by_year_act.items():
            if int(y) != solve_year:
                continue
            v = float(v)
            if abs(v) <= float(min_amount):
                continue
            injected_supply[int(a)] = injected_supply.get(int(a), 0.0) + v

        if injected_supply:
            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=injected_supply,
                inventory_by_year=inventory_total_by_impact_year,
                min_amount=float(min_amount),
                use_temporal_distributions=True,
            )

        # Build injected supply by root, normalizing legacy sentinel -> fu0
        injected_supply_by_first_level_child: Dict[int, Dict[int, float]] = {}
        for (y, act), roots_map in injected_supply_prov_by_year_act.items():
            if int(y) != solve_year:
                continue
            for r, share in roots_map.items():
                r = _normalize_root(int(r))
                share = float(share)
                if abs(share) <= float(min_amount):
                    continue
                injected_supply_by_first_level_child.setdefault(r, {})
                injected_supply_by_first_level_child[r][int(act)] = (
                    injected_supply_by_first_level_child[r].get(int(act), 0.0) + share
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
            )

            root_injected = injected_supply_by_first_level_child.get(int(root_idx), {})
            if root_injected:
                trails.accumulate_temporalized_biosphere_inventory(
                    base_year=solve_year,
                    supply_by_activity=root_injected,
                    inventory_by_year=inventory_by_root_by_impact_year[int(root_idx)],
                    min_amount=float(min_amount),
                    use_temporal_distributions=True,
                )

        # For diagnostics: include all roots that exist either via rooted demand or injected supply (incl. FU_DIRECT_ROOT)
        all_diag_roots = set(int(k) for k in demand_by_first_level_child.keys()) | set(
            int(k) for k in injected_supply_by_first_level_child.keys()
        )

        results_by_solve_year[solve_year] = {
            "fu_activity": fu0,
            "FU_DIRECT_ROOT": FU_DIRECT_ROOT,  # now equals fu0
            "fu_demand": fu_demand,
            "demand_by_first_level_child": demand_by_first_level_child,
            "injected_supply": injected_supply,
            "injected_supply_by_first_level_child": injected_supply_by_first_level_child,
            "lca": lca_total,
            "scores": 0.0,
            "scores_by_first_level_child": {int(k): 0.0 for k in sorted(all_diag_roots)},
        }

    # -----------------------------
    # 3) Impact-year characterization
    # -----------------------------
    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        out = {"results_by_solve_year": results_by_solve_year, "results_by_impact_year": {}}
        return (out, provenance) if return_provenance else out

    impact_years = sorted(set(inventory_total_by_impact_year.keys()))

    for impact_year in impact_years:
        impact_year = int(impact_year)

        dp_key = (impact_year, True)
        if dp_key not in dp_cache:
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=impact_year,
                zero_biosphere=True,
            )
            dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

        char_key = (impact_year, True, "temporal_flowid_space")
        if char_key not in char_cache:
            biosphere_matrix_dict = {int(flow_id): int(flow_id) for flow_id in set(bio_idx.values())}
            biosphere_dict_simple = {
                (name, comp, subcomp): int(flow_id)
                for (name, comp, subcomp, unit), flow_id in bio_idx.items()
            }
            C = fill_characterization_factors_matrices(
                methods=methods,
                biosphere_matrix_dict=biosphere_matrix_dict,
                biosphere_dict=biosphere_dict_simple,
            )
            char_cache[char_key] = C
        else:
            C = char_cache[char_key]

        inv_total = inventory_total_by_impact_year.get(
            impact_year, np.zeros(n_flows, dtype=trails.value_dtype)
        )
        total_score = float(np.sum(C.dot(inv_total.reshape((-1, 1)))))

        # Merge scores by root, normalizing any legacy -1 into fu0
        scores_by_first_level_child: Dict[int, float] = {}
        for root_idx, inv_map in inventory_by_root_by_impact_year.items():
            root_norm = _normalize_root(int(root_idx))
            inv_root = inv_map.get(impact_year, np.zeros(n_flows, dtype=trails.value_dtype))
            s = float(np.sum(C.dot(inv_root.reshape((-1, 1)))))
            if abs(s) <= float(min_amount):
                continue
            scores_by_first_level_child[root_norm] = scores_by_first_level_child.get(root_norm, 0.0) + s

        results_by_impact_year[impact_year] = {
            "scores": total_score,
            "scores_by_first_level_child": scores_by_first_level_child,
        }

    out = {"results_by_solve_year": results_by_solve_year, "results_by_impact_year": results_by_impact_year}
    return (out, provenance) if return_provenance else out


def compute_node_impact_intensities(
        trails: Trails,
        nodes: List[Tuple[int, int]],
        methods: List[str],
        debug: bool = False,
) -> Dict[Tuple[int, int], float]:
    """
    Compute LCIA impact intensity for a set of (year, act_idx) nodes.

    Impact intensity = LCIA score for 1 unit of that activity in that year,
    including upstream supply chain (i.e., full LCA for {act_idx: 1}).

    Parameters
    ----------
    trails : Trails
        Trails wrapper with A/B matrices and scenario info.
    nodes : list[(year, act_idx)]
        Nodes for which we want impact intensities.
    methods : list[str]
        LCIA methods (as in fill_characterization_factors_matrices).
        For now, we assume a single method and return a scalar per node.
    debug : bool
        Print debug info if True.

    Returns
    -------
    node_intensity : dict[(year, act_idx), impact_score]
    """
    if not nodes:
        return {}

    if len(methods) != 1:
        raise ValueError(
            "compute_node_impact_intensities currently assumes a single LCIA "
            "method. Got methods=%r" % (methods,)
        )

    # Group nodes by year
    nodes_by_year: Dict[int, set[int]] = {}
    for year, act in nodes:
        year = int(year)
        nodes_by_year.setdefault(year, set()).add(int(act))

    node_intensity: Dict[Tuple[int, int], float] = {}

    # Simple caches for datapackage and C matrix per year
    dp_cache: Dict[int, Any] = {}
    char_cache: Dict[int, Any] = {}

    for year, acts in sorted(nodes_by_year.items()):
        # 1) Datapackage for this year
        cache_key = (year, bool(use_temporal_distributions))
        if cache_key not in dp_cache:
            zero_bio = use_temporal_distributions
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=year,
                zero_biosphere=zero_bio,
            )
            dp_cache[year] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[year]

        # 2) Temporary LCA to get biosphere mapping and C matrix
        if year not in char_cache:
            # Use a dummy FU: all activities at 0 except maybe one,
            # we just need biosphere_matrix_dict from bw2calc.
            # Here we pick an arbitrary activity id if available.
            any_act = next(iter(acts))
            dummy_lca = bc.LCA(demand={any_act: 1.0}, data_objs=[dp])
            dummy_lca.lci()  # build technosphere & biosphere matrices

            biosphere_matrix_dict = dummy_lca.dicts.biosphere

            # bio_idx keys: (name, compartment, subcompartment, unit) -> idx
            # We need a simplified mapping (name, comp, subcomp) -> idx
            _, _, bio_meta = None, None, None

            # We can re-use the biosphere_indices built in dp helper,
            # but we only have bio_idx here if we return it from helper.
            # However, build_datapackage_for_year_from_trails already returned bio_idx:
            #   biosphere_indices: Dict[(name, comp, subcomp, unit), int]
            bio_idx = dp_cache[year][2]

            biosphere_dict_simple = {
                (name, comp, subcomp): idx
                for (name, comp, subcomp, unit), idx in bio_idx.items()
            }

            C = fill_characterization_factors_matrices(
                methods=methods,
                biosphere_matrix_dict=biosphere_matrix_dict,
                biosphere_dict=biosphere_dict_simple,
                debug=debug,
            )
            char_cache[year] = C
        else:
            C = char_cache[year]

        # 3) For each activity in this year, run a 1-unit LCA
        for act in sorted(acts):
            lca_node = bc.LCA(demand={act: 1.0}, data_objs=[dp])
            lca_node.lci()
            inv = lca_node.inventory  # biosphere vector

            scores_vec = C.dot(inv)
            # Single method -> scalar
            score_scalar = float(np.sum(scores_vec))

            node_intensity[(year, act)] = score_scalar

    return node_intensity

