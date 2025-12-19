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
    return_provenance: bool = False,
    use_temporal_distributions: bool = True,
) -> Dict[int, Dict[str, Any]]:
    """
    Temporal LCA + LCIA for one starting activity, with attribution of
    impacts to *first-level* suppliers.

    Static mode (use_temporal_distributions=False):
      - bw2calc inventory is characterized directly with CF matrix C
      - supplier roots are computed with separate bw2calc solves
      - start activity direct-only score is computed from Trails.B row and characterized

    Temporal mode (use_temporal_distributions=True):
      - biosphere is zeroed in datapackage
      - temporalized biosphere inventories are accumulated from Trails.B with temporal shifting
      - characterization happens in a second pass
      - start activity root is direct-only (only its own biosphere)
    """

    logger.info(
        "LCA start: start_year=%d start_act_idx=%d amount=%g max_depth=%d min_amount=%g methods=%s use_td=%s",
        start_year, start_act_idx, amount, max_depth, min_amount, methods, use_temporal_distributions
    )

    # -----------------------------
    # 1) Temporal traversal (used to build f_by_year; provenance optional)
    # -----------------------------
    frontier, provenance = trails.temporal_traversal(
        start_year=start_year,
        start_act_idx=start_act_idx,
        amount=amount,
        max_depth=max_depth,
        min_amount=min_amount,
        return_provenance=True if return_provenance else False,
        show_progress=show_progress,
        use_temporal_distributions=use_temporal_distributions,
    )
    logger.info("LCA: traversal frontier years=%d", len(frontier))
    logger.debug("LCA: frontier years=%s", sorted(frontier.keys()))
    # -----------------------------

    # 2) Frontier -> per-year demand vectors (used only to get FU per year and diagnostics)
    f_by_year = trails.frontier_to_demand_vectors(frontier)

    nz = {y: v for y, v in f_by_year.items() if np.any(np.asarray(v) != 0)}
    logger.info("LCA: FU-by-year years=%d nonzero=%d", len(f_by_year), len(nz))
    logger.debug("LCA: FU-by-year (first 25 nonzero)=%s", list(sorted(nz.items()))[:25])
    if not nz:
        logger.error("LCA: FU-by-year is all zero -> results will be empty. Check traversal sign/min_amount/expansion.")

    # -----------------------------
    # Caches
    # -----------------------------
    # dp_cache[(year, zero_bio)] -> (dp, tech_idx, bio_idx, uncertain)
    dp_cache: Dict[tuple, Any] = {}
    # char_cache[(year, zero_bio)] -> CF matrix
    char_cache: Dict[tuple, Any] = {}

    results_by_year: Dict[int, Dict[str, Any]] = {}

    # -----------------------------
    # Temporal-mode inventory accumulators
    # -----------------------------
    inventory_by_year_total: Dict[int, np.ndarray] = {}
    inventory_by_year_per_root: Dict[int, Dict[int, np.ndarray]] = defaultdict(dict)

    candidate_years = sorted(f_by_year.keys())
    if not candidate_years:
        return ({}, provenance) if return_provenance else {}

    # -----------------------------
    # Optional progress bar
    # -----------------------------
    bar = None
    if show_progress:
        try:
            import pyprind
            bar = pyprind.ProgBar(len(candidate_years), title="Temporal LCA over years")
        except Exception:
            bar = None

    # -----------------------------
    # Helper: build C in bw2calc biosphere ordering (static mode)
    # -----------------------------
    def _get_C_static(year: int, lca_obj, bio_idx: Dict[tuple, int], zero_bio: bool):
        """
        Build characterization matrix C such that C.shape[1] matches
        lca_obj.inventory.shape[0], using lca_obj.dicts.biosphere ordering.
        """
        key = (year, zero_bio, "static_bw_order")
        if key in char_cache:
            return char_cache[key]

        bw_bio_map = lca_obj.dicts.biosphere  # flow_id -> row position
        biosphere_matrix_dict = {int(flow_id): int(pos) for flow_id, pos in bw_bio_map.items()}

        # bio_idx: (name, comp, subcomp, unit) -> flow_id
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

    # -----------------------------
    # Helper: score a bw2calc LCA inventory with C
    # -----------------------------
    def _score_from_lca_inventory(lca_obj, C) -> float:
        inv = lca_obj.inventory
        inv_vec = np.asarray(inv.sum(axis=1)).reshape((-1, 1))
        if C.shape[1] != inv_vec.shape[0]:
            raise ValueError(f"C columns ({C.shape[1]}) != inventory rows ({inv_vec.shape[0]})")
        return float(np.sum(C.dot(inv_vec)))

    # -----------------------------
    # Main per-year loop
    # -----------------------------

    # Diagnostics: per-year screening and root coverage
    n_years = len(candidate_years)
    n_skip_fu = 0
    n_process = 0
    max_abs_fu = 0.0
    max_abs_fu_year = None

    # Diagnostics: root building and temporal accumulation
    n_years_with_any_root = 0
    n_years_with_supplier_roots = 0  # roots excluding start_act direct-only bucket
    max_roots = 0
    max_roots_year = None

    for year in candidate_years:
        # Demand vector for this calendar year coming from traversal frontier
        f_vec = f_by_year[year]
        arr = np.asarray(f_vec)

        # Process this year if there is ANY demand above threshold
        nz_idx = np.where(np.abs(arr) > min_amount)[0]
        if nz_idx.size == 0:
            n_skip_fu += 1
            if bar:
                bar.update()
            continue

        n_process += 1

        # FU for bw2calc is the whole demand vector for this year (row space)
        fu_total = {int(i): float(arr[i]) for i in nz_idx}

        n_process += 1

        # Datapackage: zero biosphere in temporal mode; keep in static mode
        zero_bio = bool(use_temporal_distributions)
        dp_key = (year, zero_bio)

        if dp_key not in dp_cache:
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=year,
                zero_biosphere=zero_bio,
            )
            dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

        # Total FU for bw2calc is ONLY the real FU
        # Total LCI (always needed)
        lca_total = bc.LCA(demand=fu_total, data_objs=[dp])
        lca_total.lci()


        # Define a scalar amount to expand the start activity for root attribution.
        # Use the actual demand for start_act_idx in this year, if present; else 0.
        # -----------------------------
        # Build first-level supplier roots (single-activity FUs)
        # -----------------------------

        # Scalar demand for the start activity in this year (used only for root attribution)
        start_amt_year = float(arr[start_act_idx]) if start_act_idx < arr.size else 0.0

        first_level = {}
        if abs(start_amt_year) > min_amount:
            first_level = trails.expand_temporal_exchanges(
                year=year,
                act_idx=start_act_idx,
                amount=start_amt_year,
                use_temporal_distributions=use_temporal_distributions,
            )

        fu_per_root: Dict[int, Dict[int, float]] = {}

        # Each first-level supplier root: {supplier_act: supplier_amt}
        for child_year, mapping in first_level.items():
            if child_year != year:
                # If you want cross-year supplier roots, handle separately.
                continue
            for supplier_act, supplier_amt in mapping.items():
                supplier_amt = float(supplier_amt)
                if abs(supplier_amt) <= min_amount:
                    continue
                fu_per_root[int(supplier_act)] = {int(supplier_act): supplier_amt}

        # Ensure start root exists (direct-only bucket)
        fu_per_root.setdefault(int(start_act_idx), {})

        # Diagnostics: root coverage
        n_roots = len(fu_per_root)
        n_years_with_any_root += 1 if n_roots > 0 else 0
        n_suppliers = len([k for k in fu_per_root.keys() if k != int(start_act_idx)])
        if n_suppliers > 0:
            n_years_with_supplier_roots += 1

        if n_roots > max_roots:
            max_roots = n_roots
            max_roots_year = year

        # -----------------------------
        # Branch: STATIC scoring
        # -----------------------------
        if not use_temporal_distributions:
            C = _get_C_static(year=year, lca_obj=lca_total, bio_idx=bio_idx, zero_bio=zero_bio)

            total_score_scalar = _score_from_lca_inventory(lca_total, C)

            scores_per_root: Dict[int, float] = {}

            # (A) Direct-only for start activity (include only if non-zero)
            bw_bio_map = lca_total.dicts.biosphere  # flow_id -> row position
            inv_direct = np.zeros(len(bw_bio_map), dtype=float)

            # Build direct biosphere vector from Trails.B row in flow-id space and map into bw2calc ordering
            scenario_year = trails._map_year_to_scenario_year(year)
            label = str(scenario_year)
            if label in trails.scenario_index:
                t = trails.scenario_index[label]
                B_row = trails.B[t, start_act_idx, :]  # (activities x flows) row slice (COO)

                coords = B_row.coords
                data = B_row.data

                if data.size > 0:
                    if coords.shape[0] == 2:
                        flow_ids = coords[1]
                    else:
                        flow_ids = coords[0]

                    for flow_id, v in zip(flow_ids, data):
                        pos = bw_bio_map.get(int(flow_id))
                        if pos is None:
                            continue
                        inv_direct[int(pos)] += float(v) * fu_amt

            direct_score = float(np.sum(C.dot(inv_direct.reshape((-1, 1)))))
            if abs(direct_score) > min_amount:
                scores_per_root[int(start_act_idx)] = direct_score

            # (B) Supplier roots: full upstream burdens
            for root_idx, fu_root in fu_per_root.items():
                if root_idx == int(start_act_idx):
                    continue
                lca_root = bc.LCA(demand=fu_root, data_objs=[dp])
                lca_root.lci()
                s = _score_from_lca_inventory(lca_root, C)
                if abs(s) > min_amount:
                    scores_per_root[int(root_idx)] = float(s)

            results_by_year[year] = {
                "fu": fu_total,
                "fu_per_root": fu_per_root,
                "lca": lca_total,
                "scores": float(total_score_scalar),
                "scores_per_root": scores_per_root,
            }

            if bar:
                bar.update()
            continue  # static year done

        # -----------------------------
        # Branch: TEMPORAL accumulation (score later in second pass)
        # -----------------------------
        # Total temporalized inventory
        supply_by_act_idx: Dict[int, float] = {}
        n_acts = trails.B.shape[1]

        for act_idx in range(n_acts):
            try:
                prod_pos = lca_total.dicts.product[act_idx]
            except KeyError:
                continue
            supply = float(lca_total.supply_array[prod_pos])
            if abs(supply) > min_amount:
                supply_by_act_idx[int(act_idx)] = supply

        trails.accumulate_temporalized_biosphere_inventory(
            year,
            supply_by_act_idx,
            inventory_by_year_total,
            min_amount=min_amount,
            use_temporal_distributions=True,
        )

        # Diagnostics: total inventory added this iteration (in flow-id space)
        inv_now = inventory_by_year_total.get(year)
        if inv_now is not None:
            l1 = float(np.sum(np.abs(inv_now)))
            nnz = int(np.count_nonzero(inv_now))
            if n_process <= 3 or (n_process % 25 == 0):
                logger.info("LCA: year=%d inv_total nnz=%d L1=%g", year, nnz, l1)

        # Per-root temporalized inventories
        for root_idx, fu_root in fu_per_root.items():
            if root_idx == int(start_act_idx):
                # Direct-only: only the start activity itself
                supply_by_act_root = {int(start_act_idx): float(start_amt_year)}

            else:
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
                year,
                supply_by_act_root,
                inventory_by_year_per_root[int(root_idx)],
                min_amount=min_amount,
                use_temporal_distributions=True,
            )

        # Diagnostics: per-root inventories present for this year
        if n_process <= 3 or (n_process % 25 == 0):
            root_counts = 0
            for r, inv_map in inventory_by_year_per_root.items():
                inv_r = inv_map.get(year)
                if inv_r is None:
                    continue
                if np.any(inv_r != 0):
                    root_counts += 1
            logger.info("LCA: year=%d per-root inventories nonzero_roots=%d", year, root_counts)

        # Placeholders; overwritten in second pass
        results_by_year[year] = {
            "fu": fu_total,
            "fu_per_root": fu_per_root,
            "lca": lca_total,
            "scores": 0.0,
            "scores_per_root": {int(k): 0.0 for k in fu_per_root.keys()},
        }

        if bar:
            bar.update()

    logger.info(
        "LCA: per-year summary: candidate_years=%d processed=%d skipped_fu=%d "
        "max_abs_fu=%g (year=%s) max_roots=%d (year=%s) years_with_supplier_roots=%d",
        n_years, n_process, n_skip_fu,
        float(max_abs_fu), str(max_abs_fu_year),
        int(max_roots), str(max_roots_year),
        int(n_years_with_supplier_roots),
    )

    if bar:
        try:
            bar.stop()
        except Exception:
            pass

    # -----------------------------
    # Second pass: characterize temporalized inventories (only temporal mode)
    # -----------------------------
    if use_temporal_distributions:
        # inventory vectors produced by trails.accumulate_* are assumed indexed by flow_id (position == flow_id)
        n_flows = int(trails.B.shape[2]) if trails.B is not None else 0
        if n_flows <= 0:
            return (results_by_year, provenance) if return_provenance else results_by_year

        # Build C once per effective year (can differ if CFs vary by year)
        all_years = sorted(set(results_by_year.keys()) | set(inventory_by_year_total.keys()))

        for y_eff in all_years:
            dp_key = (y_eff, True)
            if dp_key not in dp_cache:
                dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                    trails=trails,
                    year=y_eff,
                    zero_biosphere=True,
                )
                dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
            else:
                dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

            # Characterization for temporal inventories:
            # inventory vectors are length n_flows and indexed by flow_id,
            # so we map flow_id -> flow_id (identity).
            char_key = (y_eff, True, "temporal_flowid_space")
            if char_key not in char_cache:
                max_flow_id = max(int(v) for v in bio_idx.values()) if bio_idx else -1
                if max_flow_id >= n_flows:
                    raise ValueError(
                        f"Temporal inventory length n_flows={n_flows} but max biosphere flow id={max_flow_id}."
                        " Your temporal inventory is not in pure flow-id space."
                    )

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

            inv_total = inventory_by_year_total.get(
                y_eff,
                np.zeros(n_flows, dtype=trails.value_dtype),
            )
            total_score_scalar = float(np.sum(C.dot(inv_total.reshape((-1, 1)))))

            # Diagnostics: second pass total inventory mass and score
            nnz_total = int(np.count_nonzero(inv_total))
            l1_total = float(np.sum(np.abs(inv_total)))
            logger.info(
                "LCIA second pass: y_eff=%d inv_total nnz=%d L1=%g score=%g",
                y_eff, nnz_total, l1_total, float(total_score_scalar)
            )

            scores_per_root_2: Dict[int, float] = {}
            for root_idx, inv_map in inventory_by_year_per_root.items():
                inv_root = inv_map.get(
                    y_eff,
                    np.zeros(n_flows, dtype=trails.value_dtype),
                )
                score_val = float(np.sum(C.dot(inv_root.reshape((-1, 1)))))

                # Always keep the start activity root (direct-only bucket), even if tiny/zero
                if root_idx == int(start_act_idx):
                    scores_per_root_2[int(root_idx)] = score_val
                    continue

                # Keep supplier roots only if meaningful
                if abs(score_val) > min_amount:
                    scores_per_root_2[int(root_idx)] = score_val

            logger.info(
                "LCIA second pass: y_eff=%d scores_per_root kept=%d dropped_by_threshold=%s",
                y_eff, len(scores_per_root_2),
                "yes (start root dropped if <= min_amount)"  # reminder of your existing filter
            )

            if y_eff not in results_by_year:
                results_by_year[y_eff] = {
                    "fu": {},
                    "fu_per_root": {},
                    "lca": None,
                    "scores": total_score_scalar,
                    "scores_per_root": scores_per_root_2,
                }
            else:
                results_by_year[y_eff]["scores"] = total_score_scalar
                results_by_year[y_eff]["scores_per_root"] = scores_per_root_2

    if return_provenance:
        return results_by_year, provenance

    return results_by_year


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

