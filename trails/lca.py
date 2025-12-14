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
    remove_uncertainty: bool = True,
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
    A_coords = A_t.coords  # (2, nnz) -> [activities, products]
    A_data = -A_t.data
    nnz_A = A_data.size

    indices_A = np.empty(nnz_A, dtype=bwp.INDICES_DTYPE)

    act_idx = A_coords[0].astype(int)
    prod_idx = A_coords[1].astype(int)

    # Convention consistent with your CSV pipeline:
    #   row = product (input), col = activity (producer)
    indices_A["row"] = prod_idx
    indices_A["col"] = act_idx

    data_A = np.asarray(A_data, dtype=trails.value_dtype)
    flip_A = np.zeros(nnz_A, dtype=bool)

    dist_A = np.zeros(nnz_A, dtype=bwp.UNCERTAINTY_DTYPE)
    if remove_uncertainty:
        dist_A["uncertainty_type"] = 0
        dist_A["negative"] = False
    else:
        # Placeholder for future uncertainty integration
        dist_A["uncertainty_type"] = 0
        dist_A["negative"] = False

    # ------------------------------------------------------------------
    # 4) Build biosphere entries from B_t
    # ------------------------------------------------------------------
    B_coords = B_t.coords  # (2, nnz) -> [activities, flows]
    B_data = B_t.data
    nnz_B = B_data.size

    indices_B = np.empty(nnz_B, dtype=bwp.INDICES_DTYPE)

    act_idx_B = B_coords[0].astype(int)
    flow_idx = B_coords[1].astype(int)

    indices_B["row"] = flow_idx
    indices_B["col"] = act_idx_B

    data_B = np.asarray(B_data, dtype=trails.value_dtype)
    if zero_biosphere:
        data_B = np.zeros_like(data_B)
    flip_B = None

    dist_B = np.zeros(nnz_B, dtype=bwp.UNCERTAINTY_DTYPE)
    if remove_uncertainty:
        dist_B["uncertainty_type"] = 0
        dist_B["negative"] = False
    else:
        dist_B["uncertainty_type"] = 0
        dist_B["negative"] = False

    # ------------------------------------------------------------------
    # 5) SAFETY CHECK: matrix indices compatible with metadata
    # ------------------------------------------------------------------
    meta_act_indices = set(act_meta.keys())
    meta_bio_indices = set(bio_meta.keys())

    A_act_indices = set(act_idx.tolist())
    A_prod_indices = set(prod_idx.tolist())
    B_act_indices = set(act_idx_B.tolist())
    B_flow_indices = set(flow_idx.tolist())

    # All activity indices used anywhere in A (producer or product)
    A_all_activities = A_act_indices | A_prod_indices

    missing_activities = A_all_activities - meta_act_indices
    if missing_activities:
        print(
            f"⚠️ WARNING: A[{label_for_matrix}] uses {len(missing_activities)} "
            f"activity indices not present in activity_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_activities))[:10]}"
            f"{' ...' if len(missing_activities) > 10 else ''}"
        )
        # Optionally, you can raise instead:
        # raise ValueError(
        #     f"Matrix year {label_for_matrix} uses activity indices not in "
        #     f"metadata {meta_label}: {sorted(missing_activities)}"
        # )

    missing_flows = B_flow_indices - meta_bio_indices
    if missing_flows:
        print(
            f"⚠️ WARNING: B[{label_for_matrix}] uses {len(missing_flows)} "
            f"biosphere flow indices not present in biosphere_indices[{meta_label}]. "
            f"Examples: {sorted(list(missing_flows))[:10]}"
            f"{' ...' if len(missing_flows) > 10 else ''}"
        )
        # Optional strict mode:
        # raise ValueError(
        #     f"Matrix year {label_for_matrix} uses biosphere indices not in "
        #     f"metadata {meta_label}: {sorted(missing_flows)}"
        # )

    # ------------------------------------------------------------------
    # 6) Create bw_processing datapackage
    # ------------------------------------------------------------------
    dp = bwp.create_datapackage()

    dp.add_persistent_vector(
        matrix="technosphere_matrix",
        indices_array=indices_A,
        data_array=data_A,
        flip_array=flip_A,
        distributions_array=dist_A,
    )

    dp.add_persistent_vector(
        matrix="biosphere_matrix",
        indices_array=indices_B,
        data_array=data_B,
        flip_array=flip_B,
        distributions_array=dist_B,
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
    remove_uncertainty: bool = True,
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

    # 3) Build per-year, per-root demand from provenance:
    # fu_per_root_by_year[year][root_act][act_idx] = amount
    fu_per_root_by_year: Dict[int, Dict[int, Dict[int, float]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(float))
    )

    for (year, act_idx), path_map in provenance.items():
        for path, amt in path_map.items():
            if not path:
                continue
            root_year, root_act = path[0]
            fu_per_root_by_year[year][root_act][act_idx] += float(amt)

    # Explicitly add the start activity as a root (captures direct emissions and any unattributed share)
    for year, f_vec in f_by_year.items():
        if start_act_idx < len(f_vec) and float(f_vec[start_act_idx]) != 0.0:
            fu_per_root_by_year[year][start_act_idx][start_act_idx] += float(f_vec[start_act_idx])

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
        nonzero_indices = np.where(f != 0)[0]
        if nonzero_indices.size == 0:
            if bar:
                bar.update()
            continue

        if debug:
            print(
                f"Temporal LCIA: year={year}, "
                f"nonzero activities={len(nonzero_indices)}"
            )

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
                remove_uncertainty=remove_uncertainty,
            )
            dp_cache[dp_key] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[dp_key]

        # 5) Total functional unit for this year
        fu_total = {int(i): float(f[i]) for i in nonzero_indices}

        # 6) LCI for total FU
        lca_total = bc.LCA(demand=fu_total, data_objs=[dp])
        lca_total.lci()

        # ------------------------------------------------------------------
        # Build per-root FU maps for this year (used in BOTH modes)
        # fu_per_root[root_idx] = {act_idx: amount, ...}
        # ------------------------------------------------------------------
        roots_for_year = fu_per_root_by_year.get(year, {})
        fu_per_root: Dict[int, Dict[int, float]] = {}
        for root_idx, fu_map in roots_for_year.items():
            fu_root = {int(i): float(v) for i, v in fu_map.items() if v != 0.0}
            if fu_root:
                fu_per_root[root_idx] = fu_root

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
                    # "Direct-only" contribution for the start activity root:
                    # use the demand of the start activity in this year (if present)
                    direct_amt = float(fu_root.get(start_act_idx, amount))
                    supply_by_act_root = {int(start_act_idx): direct_amt}
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
                    remove_uncertainty=remove_uncertainty,
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
                scores_per_root_2[root_idx] = float(np.sum(root_scores))

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
        remove_uncertainty: bool = True,
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
    remove_uncertainty : bool
        Passed through to build_datapackage_for_year_from_trails.
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
                remove_uncertainty=remove_uncertainty,
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

