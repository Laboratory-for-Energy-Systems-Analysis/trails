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
) -> Dict[int, Dict[str, Any]]:

    """
    Temporal LCA + LCIA for one starting activity, with attribution of
    impacts to *first-level* suppliers.

    For each year:
      1. Temporal traversal -> (year, activity) frontier + provenance.
      2. Frontier -> demand vector f_t.
      3. From provenance -> per-root (first-level child) demand vectors.
      4. Build bw_processing datapackage from Trails.A/B.
      5. Run bw2calc LCI:
         - once for the total FU,
         - once per root FU to get scores_per_root.
      6. Build LCIA CF matrix for that year's biosphere mapping.
      7. Compute LCIA scores = C @ inventory.

    Returns
    -------
    results_by_year : dict[int, dict]
        For each year:
            {
                "fu": {act_idx: amount, ...},
                "fu_per_root": {root_idx: {act_idx: amount, ...}, ...},
                "lca": LCA object for total FU,
                "scores": float (total LCIA score for that year),
                "scores_per_root": {root_idx: float, ...},
            }
    """
    # ------------------------------------------------------------------
    # 1) Temporal traversal WITH provenance (path-based)
    frontier, provenance = trails.temporal_traversal(
        start_year=start_year,
        start_act_idx=start_act_idx,
        amount=amount,
        max_depth=max_depth,
        min_amount=min_amount,
        return_provenance=True,
        show_progress=show_progress,
    )

    # 2) Frontier -> per-year total demand vectors
    f_by_year = trails.frontier_to_demand_vectors(frontier)

    # 3) Build per-year, per-root demand from provenance
    #    fu_per_root_by_year[year][root_idx][act_idx] = amount
    # 3) Build per-year, per-root demand from provenance
    #    fu_per_root_by_year[year][root_act][act_idx] = amount
    fu_per_root_by_year: Dict[int, Dict[int, Dict[int, float]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(float))
    )

    for (year, act_idx), path_map in provenance.items():
        for path, amt in path_map.items():
            if not path:
                # No first-level supplier (root itself); skip to keep semantics
                continue

            # path[0] is (root_year, root_act)
            root_year, root_act = path[0]
            fu_per_root_by_year[year][root_act][act_idx] += float(amt)

    results_by_year: Dict[int, Dict[str, Any]] = {}

    # Optional caches
    dp_cache: Dict[int, Any] = {}   # year -> (dp, tech_idx, bio_idx, uncertain)
    char_cache: Dict[int, Any] = {} # year -> CF matrix

    candidate_years = sorted(f_by_year.keys())
    total = max(len(candidate_years), 1)
    bar = pyprind.ProgBar(total, title="Temporal LCA over years")

    for year in candidate_years:
        f = f_by_year[year]
        nonzero_indices = np.where(f != 0)[0]
        if nonzero_indices.size == 0:
            bar.update()
            continue

        if debug:
            print(
                f"Temporal LCIA: year={year}, "
                f"nonzero activities={len(nonzero_indices)}"
            )

        # 4) Build datapackage for this year (from Trails)
        if year not in dp_cache:
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=year,
                remove_uncertainty=remove_uncertainty,
            )
            dp_cache[year] = (dp, tech_idx, bio_idx, uncertain_params)
        else:
            dp, tech_idx, bio_idx, uncertain_params = dp_cache[year]

        # 5) Total functional unit for this year
        fu_total = {int(i): float(f[i]) for i in nonzero_indices}

        # 6) LCI for total FU
        lca_total = bc.LCA(demand=fu_total, data_objs=[dp])
        lca_total.lci()
        inventory_total = lca_total.inventory  # 1D biosphere vector

        # 7) Build LCIA CF matrix for this year (if not cached)
        if year not in char_cache:
            biosphere_matrix_dict = lca_total.dicts.biosphere  # {flow_id: row}

            # Our bio_idx keys: (name, compartment, subcompartment, unit) -> idx
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

        # 8) Total LCIA score
        # C shape: (n_methods, n_bio_rows)
        # inventory_total shape: (n_bio_rows,)
        total_scores = C.dot(inventory_total)
        # For now assume a single method; otherwise you may want the full vector
        total_score_scalar = float(np.sum(total_scores))

        # 9) Per-root FU and scores
        fu_per_root: Dict[int, Dict[int, float]] = {}
        scores_per_root: Dict[int, float] = {}

        roots_for_year = fu_per_root_by_year.get(year, {})

        for root_idx, fu_map in roots_for_year.items():
            # Filter zeros just in case
            fu_root = {int(i): float(v) for i, v in fu_map.items() if v != 0.0}
            if not fu_root:
                continue

            # Independent LCA for this root
            lca_root = bc.LCA(demand=fu_root, data_objs=[dp])
            lca_root.lci()
            inv_root = lca_root.inventory

            scores_root_vec = C.dot(inv_root)
            scores_root = float(np.sum(scores_root_vec))

            fu_per_root[root_idx] = fu_root
            scores_per_root[root_idx] = scores_root

            if debug:
                print(
                    f"  Year {year}, root {root_idx}: "
                    f"FU size={len(fu_root)}, score={scores_root:g}"
                )

        results_by_year[year] = {
            "fu": fu_total,
            "fu_per_root": fu_per_root,
            "lca": lca_total,
            "scores": total_score_scalar,
            "scores_per_root": scores_per_root,
        }

        bar.update()

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
        if year not in dp_cache:
            dp, tech_idx, bio_idx, uncertain_params = build_datapackage_for_year_from_trails(
                trails=trails,
                year=year,
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


