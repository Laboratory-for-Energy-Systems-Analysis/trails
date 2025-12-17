import numpy as np
import bw_processing as bwp
from typing import Dict, Tuple, List, Any
from collections import defaultdict

import bw2calc as bc
import pyprind

from .trails import Trails
from .lcia import fill_characterization_factors_matrices


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
) -> Dict[int, Dict[str, Any]]:
    """
    Temporal LCA + LCIA for one starting activity, with attribution of
    impacts to *first-level* suppliers.

    If use_temporal_distributions=False, no temporal shifting is applied and
    scores are computed directly from each bw2calc inventory.

    If use_temporal_distributions=True, biosphere is zeroed in the datapackage,
    biosphere flows are accumulated with temporal shifting, and scoring happens
    in the second pass.
    """
    # ------------------------------------------------------------------
    # 1) Temporal traversal WITH provenance (needed for per-root attribution)
    # ------------------------------------------------------------------
    frontier, provenance = trails.temporal_traversal(
        start_year=start_year,
        start_act_idx=start_act_idx,
        amount=amount,
        max_depth=max_depth,
        min_amount=min_amount,
        return_provenance=True,  # always needed internally for per-root attribution
        show_progress=show_progress,
        use_temporal_distributions=use_temporal_distributions,
    )

    # 2) Frontier -> per-year total demand vectors
    f_by_year = trails.frontier_to_demand_vectors(frontier)

    # ------------------------------------------------------------------
    # 3) Build per-year, per-root demand from provenance:
    #     - Supplier roots: only first-level suppliers (path[0])
    #     - Start activity root is treated as DIRECT-ONLY, not a normal upstream root
    #     - Any remainder goes into a RESIDUAL bucket
    # ------------------------------------------------------------------
    from collections import defaultdict

    ROOT_DIRECT = start_act_idx  # we will treat this bucket as direct-only later
    ROOT_RESIDUAL = -1  # special bucket for any unattributed remainder

    # fu_per_root_by_year[year][root_act][act_idx] noted as "demand contributions"
    fu_per_root_by_year: Dict[int, Dict[int, Dict[int, float]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(float))
    )

    # Populate supplier-root contributions from provenance
    # (root is the first node in the path: path[0] = (child_year, child_act))
    for (year, act_idx), path_map in provenance.items():
        for path, amt in path_map.items():
            if not path:
                continue
            root_year, root_act = path[0]
            fu_per_root_by_year[year][root_act][act_idx] += float(amt)

    # Now build a residual bucket so that, for each year:
    # sum_root_fu (elementwise) + residual == total FU vector (elementwise)
    for year, f_vec in f_by_year.items():
        summed = np.zeros_like(f_vec, dtype=float)

        for root_act, mp in fu_per_root_by_year.get(year, {}).items():
            for aidx, v in mp.items():
                if 0 <= int(aidx) < len(summed):
                    summed[int(aidx)] += float(v)

        residual = f_vec.astype(float) - summed
        residual[np.abs(residual) < min_amount] = 0.0

        # The functional unit itself is not an "unattributed supplier"; never put it into residual
        if start_act_idx < len(residual):
            residual[start_act_idx] = 0.0

        # Store residual only once, and only if anything remains
        if np.any(residual != 0.0):
            for aidx, v in enumerate(residual):
                if v != 0.0:
                    fu_per_root_by_year[year][ROOT_RESIDUAL][int(aidx)] += float(v)

        # The functional unit itself is not an "unattributed supplier"; do not put it into residual
        if start_act_idx < len(residual):
            residual[start_act_idx] = 0.0

        # Put residual into a separate bucket (NOT into the start activity root)
        for aidx, v in enumerate(residual):
            if v != 0.0:
                fu_per_root_by_year[year][ROOT_RESIDUAL][int(aidx)] += float(v)

    # Also ensure the DIRECT bucket exists (even if empty), so it shows up consistently
    # We do NOT populate it here with {start_act_idx: demand}; it will be computed as direct-only later.
    for year in f_by_year.keys():
        fu_per_root_by_year[year].setdefault(ROOT_DIRECT, {})

    results_by_year: Dict[int, Dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Temporalized biosphere accounting (used only when temporal mode is on)
    # ------------------------------------------------------------------
    inventory_by_year_total: Dict[int, np.ndarray] = {}
    inventory_by_year_per_root: Dict[int, Dict[int, np.ndarray]] = defaultdict(dict)

    # ------------------------------------------------------------------
    # Caches (mode-aware)
    # ------------------------------------------------------------------
    # dp_cache[(year, zero_biosphere)] -> (dp, tech_idx, bio_idx, uncertain)
    dp_cache: Dict[tuple, Any] = {}
    # char_cache[(year, zero_biosphere)] -> CF matrix
    char_cache: Dict[tuple, Any] = {}

    candidate_years = sorted(f_by_year.keys())
    total = max(len(candidate_years), 1)

    bar = None
    if show_progress:
        bar = pyprind.ProgBar(total, title="Temporal LCA over years")

    for year in candidate_years:
        f = f_by_year[year]

        # Use ONLY the actual functional unit in bw2calc.
        # The traversal vector `f` is for attribution/diagnostics, not for bw2calc demand.
        fu_amt = float(f[start_act_idx]) if start_act_idx < len(f) else 0.0
        if abs(fu_amt) <= min_amount:
            if bar:
                bar.update()
            continue

        # ------------------------------------------------------------------
        # Build datapackage for this year
        # In temporal mode: zero biosphere (we re-inject via temporalized inventory)
        # In static mode: keep biosphere (we score directly from bw2calc inventory)
        # ------------------------------------------------------------------
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

        # 5) Total functional unit for this year
        fu_total = {int(start_act_idx): fu_amt}

        # 6) LCI for total FU
        lca_total = bc.LCA(demand=fu_total, data_objs=[dp])
        lca_total.lci()

        # ------------------------------------------------------------------
        # Build per-root FU maps for this year (used in BOTH modes)
        # fu_per_root[root_idx] = {act_idx: amount, ...}
        # ------------------------------------------------------------------
        # First-level supplier demands from the FU (one-step expansion)
        first_level = trails.expand_temporal_exchanges(
            year=year,
            act_idx=start_act_idx,
            amount=fu_amt,
            use_temporal_distributions=use_temporal_distributions,
        )

        fu_per_root: Dict[int, Dict[int, float]] = {}

        # Each first-level supplier becomes a root with a single-activity FU: {supplier_act: supplier_amt}
        for child_year, mapping in first_level.items():
            if child_year != year:
                # If you allow temporal shifting in technosphere, you can decide how to handle cross-year roots.
                # Simplest: ignore or accumulate separately; for now, keep only same-year roots
                continue
            for supplier_act, supplier_amt in mapping.items():
                supplier_amt = float(supplier_amt)
                if abs(supplier_amt) <= min_amount:
                    continue
                fu_per_root[int(supplier_act)] = {int(supplier_act): supplier_amt}

        # ------------------------------------------------------------------
        # Scoring branch
        # ------------------------------------------------------------------
        if not use_temporal_distributions:
            # ---------------------------
            # STATIC scoring (no temporal shifting)
            # ---------------------------
            char_key = (year, zero_bio)  # zero_bio is False here

            if char_key not in char_cache:
                biosphere_matrix_dict = {int(idx): int(idx) for idx in set(bio_idx.values())}
                biosphere_dict_simple = {
                    (name, comp, subcomp): idx
                    for (name, comp, subcomp, unit), idx in bio_idx.items()
                }
                C = fill_characterization_factors_matrices(
                    methods=methods,
                    biosphere_matrix_dict=biosphere_matrix_dict,
                    biosphere_dict=biosphere_dict_simple,
                )
                char_cache[char_key] = C
            else:
                C = char_cache[char_key]

            def _score_from_lca_inventory(lca_obj) -> float:
                inv = lca_obj.inventory
                inv_vec = np.asarray(inv.sum(axis=1)).reshape((-1, 1))
                scores_vec = C.dot(inv_vec)
                return float(np.sum(scores_vec))

            # Total score
            total_score_scalar = _score_from_lca_inventory(lca_total)

            # --- Direct-only contribution for the start activity (no upstream) ---
            # Compute direct biosphere inventory for the start activity only, then characterize
            try:
                fu_amt = float(f_by_year[year][start_act_idx])
            except Exception:
                fu_amt = 0.0

            if abs(fu_amt) > min_amount:
                # Build direct inventory vector (biosphere flows) for this activity
                # Your internal B slice is (activities x flows) in Trails
                scenario_year = trails._map_year_to_scenario_year(year)
                label = str(scenario_year)
                if label in trails.scenario_index:
                    t = trails.scenario_index[label]
                    B_row = trails.B[t, start_act_idx, :]  # sparse.COO row slice
                    inv_direct = np.zeros(len(bio_idx), dtype=float)

                    # B_row.coords for 2D slice: [flows] if it collapses; safest: use coords
                    coords = B_row.coords
                    data = B_row.data

                    # Handle possible shapes from sparse slicing
                    if coords.shape[0] == 2:
                        flow_idx = coords[1]
                    else:
                        flow_idx = coords[0]

                    inv_direct[np.asarray(flow_idx, dtype=int)] = np.asarray(data, dtype=float) * fu_amt

                    direct_score = float(np.sum(C.dot(inv_direct.reshape((-1, 1)))))
                    if abs(direct_score) > min_amount:
                        scores_per_root[start_act_idx] = direct_score

            # Per-root scores
            scores_per_root: Dict[int, float] = {}
            for root_idx, fu_root in fu_per_root.items():
                lca_root = bc.LCA(demand=fu_root, data_objs=[dp])
                lca_root.lci()
                scores_per_root[root_idx] = _score_from_lca_inventory(lca_root)

        else:
            # ---------------------------
            # TEMPORAL scoring (existing logic)
            # - build temporally shifted inventories now
            # - score in second pass later
            # ---------------------------
            supply_by_act_idx: Dict[int, float] = {}
            n_acts = trails.B.shape[1]

            fu_per_root.setdefault(start_act_idx, {})

            for act_idx in range(n_acts):
                try:
                    prod_pos = lca_total.dicts.product[act_idx]
                except KeyError:
                    continue
                supply = float(lca_total.supply_array[prod_pos])
                if abs(supply) > min_amount:
                    supply_by_act_idx[act_idx] = supply

            trails.accumulate_temporalized_biosphere_inventory(
                year,
                supply_by_act_idx,
                inventory_by_year_total,
                min_amount=min_amount,
                use_temporal_distributions=True,
            )

        # Per-root temporalized inventories (scored in second pass)
        for root_idx, fu_root in fu_per_root.items():

            if root_idx == start_act_idx:
                # DIRECT-ONLY bucket:
                # Only count biosphere emissions that belong to the start activity itself.
                # Do NOT run an LCA solve, as that would include the full upstream chain.
                supply_start = float(f_by_year[year][start_act_idx]) if start_act_idx < len(
                    f_by_year[year]) else 0.0
                if abs(supply_start) <= min_amount:
                    continue
                supply_by_act_root = {int(start_act_idx): supply_start}

            else:
                # Normal roots (supplier roots + residual bucket): solve and use their supply arrays
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
                        supply_by_act_root[act_idx] = supply

            trails.accumulate_temporalized_biosphere_inventory(
                year,
                supply_by_act_root,
                inventory_by_year_per_root[root_idx],
                min_amount=min_amount,
                use_temporal_distributions=True,
            )

            # Placeholders; overwritten in second pass
            total_score_scalar = 0.0
            scores_per_root = {root_idx: 0.0 for root_idx in fu_per_root.keys()}

        # ------------------------------------------------------------------
        # Store first-pass results
        # ------------------------------------------------------------------
        results_by_year[year] = {
            "fu": fu_total,
            "fu_per_root": fu_per_root,
            "lca": lca_total,
            "scores": total_score_scalar,
            "scores_per_root": scores_per_root,
        }

        if bar:
            bar.update()

    # ------------------------------------------------------------------
    # Second pass: characterize temporally shifted biosphere inventories
    # (only in temporal mode)
    # ------------------------------------------------------------------
    if use_temporal_distributions:
        all_years = sorted(set(candidate_years) | set(inventory_by_year_total.keys()))
        n_flows = int(trails.B.shape[2]) if trails.B is not None else 0

        for y_eff in all_years:
            # In second pass we always need a mapping consistent with the temporal-mode datapackage:
            # zero_biosphere=True
            zero_bio_2 = True
            dp_key = (y_eff, zero_bio_2)

            if dp_key not in dp_cache:
                dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                    trails=trails,
                    year=y_eff,
                    zero_biosphere=True,
                )
                dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
            else:
                dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

            char_key = (y_eff, zero_bio_2)
            if char_key not in char_cache:
                biosphere_matrix_dict = {int(idx): int(idx) for idx in set(bio_idx.values())}
                biosphere_dict_simple = {
                    (name, comp, subcomp): idx
                    for (name, comp, subcomp, unit), idx in bio_idx.items()
                }
                C = fill_characterization_factors_matrices(
                    methods=methods,
                    biosphere_matrix_dict=biosphere_matrix_dict,
                    biosphere_dict=biosphere_dict_simple,
                )
                char_cache[char_key] = C
            else:
                C = char_cache[char_key]

            inv_vec = inventory_by_year_total.get(
                y_eff,
                np.zeros(n_flows, dtype=trails.value_dtype) if n_flows else np.zeros(0),
            )

            total_scores = C.dot(inv_vec.reshape((-1, 1)))
            total_score_scalar = float(np.sum(total_scores))

            # Per-root scores
            scores_per_root_2: Dict[int, float] = {}
            for root_idx, inv_map in inventory_by_year_per_root.items():
                inv_root = inv_map.get(
                    y_eff,
                    np.zeros(n_flows, dtype=trails.value_dtype) if n_flows else np.zeros(0),
                )
                root_scores = C.dot(inv_root.reshape((-1, 1)))
                score_val = float(np.sum(root_scores))

                # Only include start/root activity if it has characterized direct emissions
                if root_idx == start_act_idx and abs(score_val) <= min_amount:
                    continue

                # Optionally also drop near-zero entries for readability
                if abs(score_val) <= min_amount:
                    continue

                scores_per_root_2[root_idx] = score_val

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

