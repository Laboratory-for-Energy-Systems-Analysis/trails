from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
from tqdm.auto import tqdm

from .lcia import get_lcia_methods
from .trails import Trails


def _top_flow_contributions(lca_obj, bio_idx_simple, cf_vector, top_n=30):
    """
    Return top contributions by flow in BW biosphere row space.
    cf_vector is a 1D array aligned with BW biosphere rows.
    """
    inv = np.asarray(lca_obj.inventory.sum(axis=1)).ravel()  # per BW flow row
    contrib = cf_vector * inv

    # Build reverse mapping: BW row position -> flow_id
    bw_bio_map = lca_obj.dicts.biosphere  # flow_id -> row position
    pos_to_flow_id = {int(pos): int(fid) for fid, pos in bw_bio_map.items()}

    # Reverse of your (name,comp,subcomp)->flow_id map (not always 1-1, but good enough)
    flow_id_to_key = {}
    for k, fid in bio_idx_simple.items():
        flow_id_to_key.setdefault(int(fid), k)

    idx = np.argsort(np.abs(contrib))[::-1][:top_n]
    out = []
    for p in idx:
        fid = pos_to_flow_id.get(int(p))
        key = flow_id_to_key.get(int(fid), None)
        out.append(
            {
                "bw_row": int(p),
                "flow_id": None if fid is None else int(fid),
                "flow_key": key,  # (name, comp, subcomp)
                "inventory": float(inv[p]),
                "cf": float(cf_vector[p]),
                "contribution": float(contrib[p]),
            }
        )
    return out


def top_activity_contributions_from_cfvec(
    lca_obj,
    cf_vec_bw_bio_rows: np.ndarray,
    trails=None,
    year_label=None,
    top_n=30,
    min_abs=0.0,
):
    cf = np.asarray(cf_vec_bw_bio_rows, dtype=np.float64).ravel()
    B = lca_obj.biosphere_matrix
    s = np.asarray(lca_obj.supply_array, dtype=np.float64).ravel()

    char_intensity = np.asarray(cf @ B).ravel()
    contrib = char_intensity * s

    act_id_by_pos = {int(pos): int(aid) for aid, pos in lca_obj.dicts.activity.items()}

    abs_contrib = np.abs(contrib)
    idx_all = (
        np.where(abs_contrib >= float(min_abs))[0]
        if min_abs > 0
        else np.arange(contrib.size)
    )
    if idx_all.size == 0:
        return []

    idx_sorted = idx_all[np.argsort(abs_contrib[idx_all])[::-1]][: int(top_n)]

    meta = None
    if trails is not None and year_label is not None:
        meta = trails.activity_indices.get(str(year_label), None)

    out = []
    for pos in idx_sorted:
        aid = act_id_by_pos.get(int(pos))
        if aid is None:
            continue

        row = {
            "activity_id": int(aid),
            "bw_pos": int(pos),
            "supply": float(s[pos]),
            "char_intensity_per_unit_supply": float(char_intensity[pos]),
            "contribution": float(contrib[pos]),
        }

        if meta and int(aid) in meta:
            md = meta[int(aid)]
            row.update(
                {
                    "name": md.get("name"),
                    "reference product": md.get("reference product"),
                    "location": md.get("location"),
                    "unit": md.get("unit"),
                }
            )

        out.append(row)

    return out


def _build_flowkey_to_flowid(trails: Trails) -> Dict[tuple, int]:
    """Build a flow-key to flow-id mapping across labels.

    :param trails: Trails instance with biosphere metadata.
    :type trails: Trails
    :returns: Mapping of ``(name, compartment, subcompartment)`` to flow id.
    :rtype: dict[tuple, int]
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
    char_cache: dict[tuple, np.ndarray],
    debug: bool = False,
) -> np.ndarray:
    """Build a CF vector aligned with Trails flow-id space.

    :param trails: Trails instance with biosphere metadata.
    :type trails: Trails
    :param methods: LCIA methods to include.
    :type methods: list[str]
    :param ei_version: Ecoinvent release identifier.
    :type ei_version: str
    :param char_cache: Cache mapping for characterization vectors.
    :type char_cache: dict
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :returns: CF vector in flow-id space.
    :rtype: numpy.ndarray
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
    for _mname, exc in methods_dict.items():
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
    inventory_total_by_impact_year: dict[int, np.ndarray],
    inventory_by_root_by_impact_year: dict[int, dict[int, np.ndarray]],
    dp_cache: dict,
    char_cache: dict[tuple, np.ndarray],
    methods: list[str],
    min_amount: float,
    normalize_root: Any,
    debug: bool,
    ei_version: str = "3.11",
) -> dict[int, dict[str, Any]]:
    """Characterize inventories into impact scores by impact year.

    :param trails: Trails instance with inventory dimensions.
    :type trails: Trails
    :param inventory_total_by_impact_year: Total inventory by impact year.
    :type inventory_total_by_impact_year: dict[int, numpy.ndarray]
    :param inventory_by_root_by_impact_year: Inventory by root and impact year.
    :type inventory_by_root_by_impact_year: dict[int, dict[int, numpy.ndarray]]
    :param dp_cache: Datapackage cache (unused but kept for parity).
    :type dp_cache: dict
    :param char_cache: Cache for characterization vectors.
    :type char_cache: dict
    :param methods: List of LCIA methods.
    :type methods: list[str]
    :param min_amount: Minimum magnitude to include.
    :type min_amount: float
    :param normalize_root: Function to normalize root identifiers.
    :type normalize_root: callable
    :param debug: Whether to emit debug logging.
    :type debug: bool
    :param ei_version: Ecoinvent release identifier.
    :type ei_version: str
    :returns: Results by impact year.
    :rtype: dict[int, dict[str, typing.Any]]
    """
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
            inv_root = inv_map.get(
                impact_year, np.zeros(n_flows, dtype=trails.value_dtype)
            )
            s = float(np.dot(cf, inv_root.astype(np.float64, copy=False)))
            if abs(s) <= float(min_amount):
                continue
            scores_by_first_level_child[root_norm] = (
                scores_by_first_level_child.get(root_norm, 0.0) + s
            )

        results_by_impact_year[impact_year] = {
            "scores": total_score,
            "scores_by_first_level_child": scores_by_first_level_child,
        }

    return results_by_impact_year
