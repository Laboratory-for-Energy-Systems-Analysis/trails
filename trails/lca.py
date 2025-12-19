import numpy as np
import bw_processing as bwp
from typing import Dict, Tuple, List, Any
from collections import defaultdict

import bw2calc as bc
import pyprind

from .trails import Trails
from .lcia import fill_characterization_factors_matrices

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

    B_data = B_t.data
    data_B = np.asarray(B_data, dtype=np.float64)  # or trails.value_dtype if you prefer
    if zero_biosphere:
        data_B = np.zeros_like(data_B)
    flip_B = None

    # ------------------------------------------------------------------
    # 5) SAFETY CHECK: matrix indices compatible with metadata
    # ------------------------------------------------------------------
    meta_act_indices = set(act_meta.keys())
    meta_bio_indices = set(bio_meta.keys())

    # ---- Extract index vectors from sparse coords (COO: coords[0]=t, coords[1]=i, coords[2]=j) ----

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
    """

    import numpy as np
    import bw2calc as bc
    from collections import defaultdict

    # -----------------------------
    # 1) Traversal frontier (+ provenance if we need roots)
    # -----------------------------
    need_roots = True  # because Option A wants first-level children as roots

    if return_provenance or need_roots:
        frontier, provenance = trails.temporal_traversal(
            start_year=start_year,
            start_act_idx=start_act_idx,
            amount=amount,
            max_depth=max_depth,
            min_amount=min_amount,
            return_provenance=True,
            show_progress=show_progress,
            use_temporal_distributions=use_temporal_distributions,
        )
    else:
        frontier, _prov = trails.temporal_traversal(
            start_year=start_year,
            start_act_idx=start_act_idx,
            amount=amount,
            max_depth=max_depth,
            min_amount=min_amount,
            return_provenance=False,
            show_progress=show_progress,
            use_temporal_distributions=use_temporal_distributions,
        )
        provenance = {}

    # Frontier -> per-year demand vectors (calendar years preserved)
    f_by_year = trails.frontier_to_demand_vectors(frontier)

    # -----------------------------
    # 1b) Build per-root frontier using provenance paths
    #     Root = first activity inside FU (first node in path)
    # -----------------------------
    from collections import defaultdict
    import numpy as np

    # frontier_rooted[(year, act, root_act)] = amount share
    frontier_rooted = defaultdict(float)

    # Collect roots seen
    roots_seen = set()

    FALLBACK_ROOT = int(start_act_idx)
    TOL = 1e-12  # numeric closure tolerance for the residual allocation

    for (year_act, total_amt) in frontier.items():
        year, act = year_act
        year = int(year)
        act = int(act)
        total_amt = float(total_amt)

        prov = provenance.get((year, act), {})

        # If no provenance, we cannot split; put everything in fallback bucket
        if not prov:
            frontier_rooted[(year, act, FALLBACK_ROOT)] += total_amt
            roots_seen.add(FALLBACK_ROOT)
            continue

        # Split by paths; each path begins at first-level child
        # provenance[(year,act)] is {path_tuple_or_root: amt}
        def _root_from_path(path, fallback_root: int) -> int:
            if not path:
                return int(fallback_root)

            # If provenance stores the root directly as an int
            if isinstance(path, (int, np.integer)):
                return int(path)

            # If provenance stores a path tuple/list, look at first element
            first = path[0]

            # Case A: first element is (year, act)
            if isinstance(first, (tuple, list)) and len(first) >= 2:
                return int(first[1])

            # Case B: first element is act (int)
            if isinstance(first, (int, np.integer)):
                return int(first)

            return int(fallback_root)

        prov_sum = 0.0
        for path, amt_share in prov.items():
            amt_share = float(amt_share)
            prov_sum += amt_share

            root_act = _root_from_path(path, fallback_root=FALLBACK_ROOT)
            frontier_rooted[(year, act, int(root_act))] += amt_share
            roots_seen.add(int(root_act))

        # Critical: enforce closure by assigning any residual to fallback root
        residual = total_amt - prov_sum
        if abs(residual) > TOL:
            frontier_rooted[(year, act, FALLBACK_ROOT)] += residual
            roots_seen.add(FALLBACK_ROOT)

    roots_seen = sorted(roots_seen)

    # -----------------------------
    # 1c) Convert rooted frontier to per-year demand vectors by root
    #     f_by_year_by_root[root][year] = demand vector
    # -----------------------------
    n_activities = trails.A.shape[1]
    dtype = trails.value_dtype

    f_by_year_by_root = {r: {} for r in roots_seen}

    for (year, act, root_act), amt in frontier_rooted.items():
        y = int(year)
        a = int(act)
        r = int(root_act)

        if y not in f_by_year_by_root[r]:
            f_by_year_by_root[r][y] = np.zeros(n_activities, dtype=dtype)

        f_by_year_by_root[r][y][a] += dtype(amt)

    candidate_years = sorted(f_by_year.keys())

    results_by_solve_year: Dict[int, Dict[str, Any]] = {}

    # -----------------------------
    # Caches
    # -----------------------------
    dp_cache: Dict[tuple, Any] = {}
    char_cache: Dict[tuple, Any] = {}

    # -----------------------------
    # Temporal-mode accumulators (flow-id space)
    # -----------------------------
    inventory_total_by_impact_year: Dict[int, np.ndarray] = {}
    inventory_per_root_by_impact_year: Dict[int, Dict[int, np.ndarray]] = defaultdict(dict)

    # -----------------------------
    # Helpers
    # -----------------------------
    def _get_C_static(year: int, lca_obj, bio_idx: Dict[tuple, int], zero_bio: bool):
        key = (year, zero_bio, "static_bw_order")
        if key in char_cache:
            return char_cache[key]

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
        char_cache[key] = C
        return C

    def _score_from_lca_inventory(lca_obj, C) -> float:
        inv = lca_obj.inventory
        inv_vec = np.asarray(inv.sum(axis=1)).reshape((-1, 1))
        if C.shape[1] != inv_vec.shape[0]:
            raise ValueError(f"C columns ({C.shape[1]}) != inventory rows ({inv_vec.shape[0]})")
        return float(np.sum(C.dot(inv_vec)))

    # -----------------------------
    # 2) Per-solve-year loop (CRITICAL: loop over f_by_year)
    # -----------------------------
    for solve_year in candidate_years:
        f_vec = f_by_year[solve_year]
        arr = np.asarray(f_vec)

        nz_idx = np.where(np.abs(arr) > min_amount)[0]
        if nz_idx.size == 0:
            continue

        # IMPORTANT: solve the whole demand vector for this solve_year
        fu_total = {int(i): float(arr[i]) for i in nz_idx}

        # Build dp for this solve_year
        zero_bio = bool(use_temporal_distributions)
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

        lca_total = bc.LCA(demand=fu_total, data_objs=[dp])
        lca_total.lci()

        # -----------------------------
        # First-level supplier roots (only based on start activity in THIS solve_year)
        # -----------------------------
        # -----------------------------
        # Roots: demand vectors restricted to this solve_year
        # -----------------------------
        fu_per_root: Dict[int, Dict[int, float]] = {}

        for root_idx, by_year in f_by_year_by_root.items():
            vec = by_year.get(int(solve_year), None)
            if vec is None:
                continue

            arr_r = np.asarray(vec)
            nz_r = np.where(np.abs(arr_r) > min_amount)[0]
            if nz_r.size == 0:
                continue

            fu_per_root[int(root_idx)] = {int(i): float(arr_r[i]) for i in nz_r}

        # Ensure a "direct-only / unassigned" bucket exists if you rely on it downstream
        fu_per_root.setdefault(int(start_act_idx), {})

        # Sanity check: sum of rooted FU vectors equals total FU vector (within tolerance)
        summed = np.zeros_like(arr, dtype=float)
        for fu_root in fu_per_root.values():
            for i, v in fu_root.items():
                summed[i] += v
        if np.max(np.abs(summed - arr)) > 1e-9:
            raise ValueError(f"Rooted FU does not sum to total FU in year {solve_year}")


        # -----------------------------
        # STATIC MODE: score immediately and book into impact_year == solve_year
        # -----------------------------
        if not use_temporal_distributions:
            C = _get_C_static(year=solve_year, lca_obj=lca_total, bio_idx=bio_idx, zero_bio=zero_bio)

            total_score = _score_from_lca_inventory(lca_total, C)

            scores_per_root: Dict[int, float] = {}

            # Supplier roots: full upstream burdens
            for root_idx, fu_root in fu_per_root.items():
                if root_idx == int(start_act_idx):
                    # optional: compute direct-only separately; keep 0 if you prefer
                    continue
                lca_root = bc.LCA(demand=fu_root, data_objs=[dp])
                lca_root.lci()
                s = _score_from_lca_inventory(lca_root, C)
                if abs(s) > min_amount:
                    scores_per_root[int(root_idx)] = float(s)

            # If you want start root shown, compute it; otherwise omit it.
            # Here we omit if it is zero/unknown in static mode.

            results_by_solve_year[solve_year] = {
                "fu": fu_total,
                "fu_per_root": fu_per_root,
                "lca": lca_total,
                "scores": float(total_score),
                "scores_per_root": scores_per_root,
            }

            continue

        # -----------------------------
        # TEMPORAL MODE: accumulate temporalized biosphere inventories
        # -----------------------------
        n_acts = trails.B.shape[1]
        supply_by_act_idx: Dict[int, float] = {}

        for act_idx in range(n_acts):
            try:
                prod_pos = lca_total.dicts.product[act_idx]
            except KeyError:
                continue
            supply = float(lca_total.supply_array[prod_pos])
            if abs(supply) > min_amount:
                supply_by_act_idx[int(act_idx)] = supply

        trails.accumulate_temporalized_biosphere_inventory(
            base_year=solve_year,
            supply_by_activity=supply_by_act_idx,
            inventory_by_year=inventory_total_by_impact_year,
            min_amount=min_amount,
            use_temporal_distributions=True,
        )

        # Per-root inventories
        for root_idx, fu_root in fu_per_root.items():
            if not fu_root:
                continue

            lca_root = bc.LCA(demand=fu_root, data_objs=[dp])
            lca_root.lci()

            supply_by_act_root = {}
            for act_idx in range(n_acts):
                try:
                    prod_pos = lca_root.dicts.product[act_idx]
                except KeyError:
                    continue
                supply = float(lca_root.supply_array[prod_pos])
                if abs(supply) > min_amount:
                    supply_by_act_root[int(act_idx)] = supply

            trails.accumulate_temporalized_biosphere_inventory(
                base_year=solve_year,
                supply_by_activity=supply_by_act_root,
                inventory_by_year=inventory_per_root_by_impact_year[int(root_idx)],
                min_amount=min_amount,
                use_temporal_distributions=True,
            )

        # Store solve-year diagnostics; leave scores empty until second pass
        results_by_solve_year[solve_year] = {
            "fu": fu_total,
            "fu_per_root": fu_per_root,
            "lca": lca_total,
            "scores": 0.0,
            "scores_per_root": {int(k): 0.0 for k in fu_per_root.keys()},
        }

    # -----------------------------
    # 3) Build results_by_impact_year
    # -----------------------------
    results_by_impact_year: Dict[int, Dict[str, Any]] = {}

    # STATIC MODE: just mirror solve-year results into impact-year space
    if not use_temporal_distributions:
        for y, r in results_by_solve_year.items():
            results_by_impact_year[int(y)] = {
                "scores": float(r.get("scores", 0.0)),
                "scores_per_root": dict(r.get("scores_per_root", {})),
            }
        out = {
            "results_by_solve_year": results_by_solve_year,
            "results_by_impact_year": results_by_impact_year,
        }
        return (out, provenance) if return_provenance else out

    # TEMPORAL MODE: characterize inventories booked in impact years
    n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
    if n_flows <= 0:
        out = {
            "results_by_solve_year": results_by_solve_year,
            "results_by_impact_year": {},
        }
        return (out, provenance) if return_provenance else out

    impact_years = sorted(set(inventory_total_by_impact_year.keys()))

    for impact_year in impact_years:
        # dp/bio mapping for characterization (flow-id space)
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

        scores_per_root: Dict[int, float] = {}
        for root_idx, inv_map in inventory_per_root_by_impact_year.items():
            inv_root = inv_map.get(
                impact_year, np.zeros(n_flows, dtype=trails.value_dtype)
            )
            s = float(np.sum(C.dot(inv_root.reshape((-1, 1)))))
            if abs(s) > min_amount:
                scores_per_root[int(root_idx)] = s

        results_by_impact_year[int(impact_year)] = {
            "scores": total_score,
            "scores_per_root": scores_per_root,
        }

    out = {
        "results_by_solve_year": results_by_solve_year,
        "results_by_impact_year": results_by_impact_year,
    }
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

            if debug:
                print(f"Impact intensity year={year}, act={act}: {score_scalar:g}")

    return node_intensity

