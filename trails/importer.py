from __future__ import annotations

from pathlib import Path
import logging
from typing import TYPE_CHECKING

import numpy as np
import sparse

from .temporal_distributions import TemporalExchange
from .utils import _parse_float_or_none, _parse_intish_or_none
from .cache_interpolation import save_cached_interpolation

if TYPE_CHECKING:
    from .trails import Trails

logger = logging.getLogger(__name__)


def import_excel_inventory(
    trails: "Trails",
    path: str | Path,
    *,
    year: int | None = None,
    scenario_label: str | None = None,
    cache_import: bool = False,
) -> dict[str, int]:
    """Import a user-provided inventory spreadsheet into A/B tensors.

    Uses ``bw2io.importers.excel.ExcelImporter`` to load inventories and applies
    its strategies before mapping exchanges onto Trails activity/biosphere indices.
    Temporal distribution fields on exchanges are parsed when present.
    Existing matrix entries at the same coordinates are replaced by default.
    When neither ``year`` nor ``scenario_label`` is provided, exchanges are
    applied to all template years and then interpolated across annual years.

    :param trails: Trails instance to update.
    :type trails: Trails
    :param path: Path to the Excel inventory file.
    :type path: str | pathlib.Path
    :param year: Calendar year to map to a scenario slice.
    :type year: int | None
    :param scenario_label: Explicit scenario label to target.
    :type scenario_label: str | None
    :param cache_import: Whether to save interpolated matrices to cache.
    :type cache_import: bool
    :returns: Summary counts of imported exchanges.
    :rtype: dict[str, int]
    """
    if trails.A is None or trails.B is None:
        raise RuntimeError("A/B matrices are not initialized.")

    from bw2io.importers.excel import ExcelImporter

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Excel inventory file not found: {path}")

    importer = ExcelImporter(str(path))
    importer.apply_strategies()

    apply_to_all_template_years = scenario_label is None and year is None

    if apply_to_all_template_years:
        template_labels = [str(lbl) for lbl in trails.template_labels]
        targets: list[tuple[str, int]] = []
        for label in template_labels:
            if label not in trails.scenario_index:
                raise ValueError(
                    f"Scenario label {label!r} not found in Trails scenarios."
                )
            targets.append((label, int(trails.scenario_index[label])))
        template_label = None
        scenario_year = None
    elif scenario_label is None:
        context = trails._get_scenario_context(int(year))
        if context is None:
            raise ValueError(f"Scenario not available for year={year}")
        scenario_year, scenario_label, t = context
        template_label = str(trails._map_year_to_template_year(scenario_year))
        targets = [(template_label, int(t))]
    else:
        if scenario_label not in trails.scenario_index:
            raise ValueError(
                f"Scenario label {scenario_label!r} not found in Trails scenarios."
            )
        t = int(trails.scenario_index[scenario_label])
        scenario_year = int(scenario_label) if scenario_label.isdigit() else None
        template_label = (
            str(trails._map_year_to_template_year(scenario_year))
            if scenario_year is not None
            else scenario_label
        )
        targets = [(template_label, int(t))]

    mapping_label = template_label or targets[0][0]

    activity_mapping = trails.activity_indices.get(mapping_label)
    if activity_mapping is None:
        seed = next(iter(trails.activity_indices.values()), {})
        activity_mapping = dict(seed)
        trails.activity_indices[mapping_label] = activity_mapping

    biosphere_mapping = trails.biosphere_indices.get(mapping_label)
    if biosphere_mapping is None:
        seed = next(iter(trails.biosphere_indices.values()), {})
        biosphere_mapping = dict(seed)
        trails.biosphere_indices[mapping_label] = biosphere_mapping

    if apply_to_all_template_years:
        for label in trails.template_labels:
            label = str(label)
            trails.activity_indices.setdefault(label, activity_mapping)
            trails.biosphere_indices.setdefault(label, biosphere_mapping)

    def _norm(value: object) -> str:
        return "" if value is None else str(value).strip()

    def _parse_amount_source(value: object) -> str:
        if value is None:
            return "port"
        s = str(value).strip().lower()
        if s == "":
            return "port"
        if s in {"port", "matrix"}:
            return s
        raise ValueError(f"Invalid temporal_amount_source: {value!r}")

    def _parse_temporal_exchange(exchange: dict) -> TemporalExchange | None:
        dist_code = _parse_intish_or_none(exchange.get("temporal_distribution"))
        if dist_code in (None, 0):
            return None

        loc = _parse_float_or_none(exchange.get("temporal_loc"))
        scale = _parse_float_or_none(exchange.get("temporal_scale"))
        off_min = _parse_intish_or_none(exchange.get("temporal_min")) or 0
        off_max = _parse_intish_or_none(exchange.get("temporal_max")) or 0
        amount_source = _parse_amount_source(exchange.get("temporal_amount_source"))

        return TemporalExchange(
            distribution=int(dist_code),
            loc=loc,
            scale=scale,
            offset_min=int(off_min),
            offset_max=int(off_max),
            amount_source=amount_source,
        )

    def _activity_key(
        name: object,
        ref_product: object,
        location: object,
    ) -> tuple[str, str, str]:
        return (
            _norm(name),
            _norm(ref_product),
            _norm(location),
        )

    def _biosphere_key(
        name: object,
        compartment: object,
        subcompartment: object,
    ) -> tuple[str, str, str]:
        return (
            _norm(name),
            _norm(compartment),
            _norm(subcompartment),
        )

    activity_lookup: dict[tuple[str, str, str], int] = {}
    for idx, meta in activity_mapping.items():
        key = _activity_key(
            meta.get("name"),
            meta.get("reference product"),
            meta.get("location"),
        )
        if key in activity_lookup and activity_lookup[key] != int(idx):
            raise ValueError(
                "Duplicate activity metadata in indices; cannot disambiguate "
                f"{key!r} (indices {activity_lookup[key]} and {idx})."
            )
        activity_lookup[key] = int(idx)

    biosphere_lookup: dict[tuple[str, str, str], int] = {}
    for idx, meta in biosphere_mapping.items():
        key = _biosphere_key(
            meta.get("name"),
            meta.get("compartment"),
            meta.get("subcompartment"),
        )
        if key in biosphere_lookup and biosphere_lookup[key] != int(idx):
            raise ValueError(
                "Duplicate biosphere metadata in indices; cannot disambiguate "
                f"{key!r} (indices {biosphere_lookup[key]} and {idx})."
            )
        biosphere_lookup[key] = int(idx)

    def _sync_indices_across_labels() -> None:
        """Ensure all label mappings include new indices with consistent ids."""
        for label, mapping in trails.activity_indices.items():
            if mapping is activity_mapping:
                continue
            lookup = {
                _activity_key(
                    meta.get("name"),
                    meta.get("reference product"),
                    meta.get("location"),
                ): int(idx)
                for idx, meta in mapping.items()
            }
            for key, idx in activity_lookup.items():
                if key in lookup and lookup[key] != idx:
                    raise ValueError(
                        "Activity indices are inconsistent across labels; "
                        f"key={key!r} label={label!r} has {lookup[key]} vs {idx}."
                    )
                if key not in lookup:
                    mapping[idx] = activity_mapping[idx]

        for label, mapping in trails.biosphere_indices.items():
            if mapping is biosphere_mapping:
                continue
            lookup = {
                _biosphere_key(
                    meta.get("name"),
                    meta.get("compartment"),
                    meta.get("subcompartment"),
                ): int(idx)
                for idx, meta in mapping.items()
            }
            for key, idx in biosphere_lookup.items():
                if key in lookup and lookup[key] != idx:
                    raise ValueError(
                        "Biosphere indices are inconsistent across labels; "
                        f"key={key!r} label={label!r} has {lookup[key]} vs {idx}."
                    )
                if key not in lookup:
                    mapping[idx] = biosphere_mapping[idx]

    def _next_index(mapping: dict[int, dict]) -> int:
        return (max(mapping) + 1) if mapping else 0

    new_activity_count = 0
    new_biosphere_count = 0
    updated_activity_ids: set[int] = set()
    updated_biosphere_ids: set[int] = set()
    new_activity_rows: list[dict[str, str]] = []
    unlinked_count = 0

    def _ensure_activity_index(
        name: object,
        ref_product: object,
        location: object,
        *,
        unit: object | None = None,
    ) -> tuple[int, bool]:
        nonlocal new_activity_count
        key = _activity_key(name, ref_product, location)
        idx = activity_lookup.get(key)
        if idx is not None:
            return int(idx), False

        idx = _next_index(activity_mapping)
        activity_mapping[idx] = {
            "name": _norm(name),
            "reference product": _norm(ref_product),
            "unit": _norm(unit),
            "location": _norm(location),
        }
        activity_lookup[key] = int(idx)
        new_activity_count += 1
        new_activity_rows.append(
            {
                "index": str(idx),
                "name": _norm(name),
                "reference product": _norm(ref_product),
                "unit": _norm(unit),
                "location": _norm(location),
            }
        )
        _sync_indices_across_labels()
        return int(idx), True

    def _ensure_biosphere_index(
        name: object,
        compartment: object,
        subcompartment: object,
        *,
        unit: object | None = None,
    ) -> tuple[int, bool]:
        nonlocal new_biosphere_count
        key = _biosphere_key(name, compartment, subcompartment)
        idx = biosphere_lookup.get(key)
        if idx is not None:
            return int(idx), False

        idx = _next_index(biosphere_mapping)
        biosphere_mapping[idx] = {
            "name": _norm(name),
            "compartment": _norm(compartment),
            "subcompartment": _norm(subcompartment),
            "unit": _norm(unit),
        }
        biosphere_lookup[key] = int(idx)
        new_biosphere_count += 1
        _sync_indices_across_labels()
        return int(idx), True

    def _activity_index_for_dataset(dataset: dict) -> int:
        idx, created = _ensure_activity_index(
            dataset.get("name"),
            dataset.get("reference product") or dataset.get("product"),
            dataset.get("location"),
            unit=dataset.get("unit"),
        )
        if not created:
            updated_activity_ids.add(int(idx))
        return idx

    def _biosphere_index_for_exchange(exchange: dict) -> int:
        categories = exchange.get("categories")
        compartment = exchange.get("compartment")
        subcompartment = exchange.get("subcompartment")
        if categories:
            if isinstance(categories, (list, tuple)):
                compartment = categories[0] if len(categories) > 0 else compartment
                subcompartment = (
                    categories[1] if len(categories) > 1 else subcompartment
                )
            else:
                compartment = categories

        key = _biosphere_key(
            exchange.get("name"),
            compartment,
            subcompartment,
        )
        idx, created = _ensure_biosphere_index(
            exchange.get("name"),
            compartment,
            subcompartment,
            unit=exchange.get("unit"),
        )
        if not created:
            updated_biosphere_ids.add(int(idx))
        return idx

    a_coords: list[tuple[int, int, int]] = []
    a_data: list[float] = []
    b_coords: list[tuple[int, int, int]] = []
    b_data: list[float] = []
    production_count = 0
    technosphere_count = 0
    biosphere_count = 0

    if trails.temporal_technosphere_exchanges is None:
        trails.temporal_technosphere_exchanges = {}
    if trails.temporal_biosphere_exchanges is None:
        trails.temporal_biosphere_exchanges = {}

    a_template_values: dict[tuple[int, int], dict[int, float]] = {}
    b_template_values: dict[tuple[int, int], dict[int, float]] = {}
    unlinked_rows: list[dict[str, str]] = []

    dataset_act_indices: dict[tuple[str, str, str], int] = {}
    affected_acts: set[int] = set()
    datasets = getattr(importer, "data", []) or []

    for dataset in datasets:
        act_idx = _activity_index_for_dataset(dataset)
        dataset_key = _activity_key(
            dataset.get("name"),
            dataset.get("reference product") or dataset.get("product"),
            dataset.get("location"),
        )
        dataset_act_indices[dataset_key] = int(act_idx)
        affected_acts.add(int(act_idx))

    # Clear existing temporal exchanges for affected activities/targets first.
    target_t = {int(t) for _, t in targets}
    for label, _t in targets:
        for act in affected_acts:
            for key in list(trails.temporal_technosphere_exchanges.keys()):
                if key[0] == label and key[1] == act:
                    del trails.temporal_technosphere_exchanges[key]
            for key in list(trails.temporal_biosphere_exchanges.keys()):
                if key[0] == label and key[1] == act:
                    del trails.temporal_biosphere_exchanges[key]

    for dataset in datasets:
        dataset_key = _activity_key(
            dataset.get("name"),
            dataset.get("reference product") or dataset.get("product"),
            dataset.get("location"),
        )
        act_idx = dataset_act_indices[dataset_key]
        for exchange in dataset.get("exchanges", []) or []:
            ex_type = exchange.get("type")
            try:
                amount = float(exchange.get("amount", 0.0))
            except (TypeError, ValueError):
                continue
            if amount == 0.0:
                continue

            if ex_type == "production":
                key = _activity_key(
                    exchange.get("name"),
                    exchange.get("reference product") or exchange.get("product"),
                    exchange.get("location"),
                )
                prod_idx = activity_lookup.get(key) or dataset_act_indices.get(key)
                if prod_idx is None:
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(exchange.get("name")),
                            "reference product": _norm(
                                exchange.get("reference product")
                                or exchange.get("product")
                            ),
                            "categories": "",
                            "location": _norm(exchange.get("location")),
                        }
                    )
                    unlinked_count += 1
                    continue
                production_count += 1
                tech_pair = (act_idx, prod_idx)
                stored_amount = amount
            elif ex_type == "technosphere":
                key = _activity_key(
                    exchange.get("name"),
                    exchange.get("reference product") or exchange.get("product"),
                    exchange.get("location"),
                )
                prod_idx = activity_lookup.get(key) or dataset_act_indices.get(key)
                if prod_idx is None:
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(exchange.get("name")),
                            "reference product": _norm(
                                exchange.get("reference product")
                                or exchange.get("product")
                            ),
                            "categories": "",
                            "location": _norm(exchange.get("location")),
                        }
                    )
                    unlinked_count += 1
                    continue
                technosphere_count += 1
                tech_pair = (act_idx, prod_idx)
                stored_amount = -amount
            elif ex_type == "biosphere":
                flow_idx = _biosphere_index_for_exchange(exchange)
                biosphere_count += 1
                stored_amount = amount
            else:
                continue

            if ex_type in {"production", "technosphere"}:
                if (
                    _norm(exchange.get("name")) == ""
                    or _norm(
                        exchange.get("reference product") or exchange.get("product")
                    )
                    == ""
                    or _norm(exchange.get("location")) == ""
                ):
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(exchange.get("name")),
                            "reference product": _norm(
                                exchange.get("reference product")
                                or exchange.get("product")
                            ),
                            "categories": "",
                            "location": _norm(exchange.get("location")),
                        }
                    )
                    unlinked_count += 1
            elif ex_type == "biosphere":
                categories = exchange.get("categories")
                if isinstance(categories, (list, tuple)):
                    cat_text = "/".join([_norm(c) for c in categories if _norm(c)])
                elif categories is None:
                    cat_text = ""
                else:
                    cat_text = _norm(categories)
                if _norm(exchange.get("name")) == "" or cat_text == "":
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(exchange.get("name")),
                            "reference product": "",
                            "categories": cat_text,
                            "location": _norm(exchange.get("location")),
                        }
                    )
                    unlinked_count += 1

            for label, t in targets:
                if ex_type in {"production", "technosphere"}:
                    if not apply_to_all_template_years:
                        a_coords.append((t, act_idx, prod_idx))
                        a_data.append(stored_amount)

                    tex = _parse_temporal_exchange(exchange)
                    key = (label, int(act_idx), int(prod_idx))
                    if tex is None:
                        trails.temporal_technosphere_exchanges.pop(key, None)
                    else:
                        trails.temporal_technosphere_exchanges[key] = tex

                    year_int = int(label) if label.isdigit() else None
                    if year_int is not None:
                        a_template_values.setdefault(tech_pair, {})[
                            year_int
                        ] = stored_amount
                elif ex_type == "biosphere":
                    if not apply_to_all_template_years:
                        b_coords.append((t, act_idx, flow_idx))
                        b_data.append(stored_amount)

                    tex = _parse_temporal_exchange(exchange)
                    key = (label, int(act_idx), int(flow_idx))
                    if tex is None:
                        trails.temporal_biosphere_exchanges.pop(key, None)
                    else:
                        trails.temporal_biosphere_exchanges[key] = tex

                    year_int = int(label) if label.isdigit() else None
                    if year_int is not None:
                        b_template_values.setdefault((act_idx, flow_idx), {})[
                            year_int
                        ] = stored_amount

    def _resize_sparse(
        matrix: sparse.COO, new_shape: tuple[int, int, int]
    ) -> sparse.COO:
        if matrix.shape == new_shape:
            return matrix
        return sparse.COO(coords=matrix.coords, data=matrix.data, shape=new_shape)

    def _replace_sparse_entries(
        matrix: sparse.COO,
        coords: list[tuple[int, int, int]],
        data: list[float],
    ) -> tuple[sparse.COO, int]:
        if not coords:
            return matrix, 0

        coords_arr = np.array(coords, dtype=trails.index_dtype).T
        data_arr = np.array(data, dtype=trails.value_dtype)

        flat = np.ravel_multi_index(coords_arr, matrix.shape)
        flat_unique, inverse = np.unique(flat, return_inverse=True)
        data_agg = np.zeros_like(flat_unique, dtype=trails.value_dtype)
        np.add.at(data_agg, inverse, data_arr)
        coords_arr = np.array(
            np.unravel_index(flat_unique, matrix.shape), dtype=trails.index_dtype
        )

        old_flat = np.ravel_multi_index(matrix.coords, matrix.shape)
        keep_mask = ~np.isin(old_flat, flat_unique)
        replaced_count = int(np.isin(flat_unique, old_flat).sum())
        if keep_mask.any():
            coords_kept = matrix.coords[:, keep_mask]
            data_kept = matrix.data[keep_mask]
            coords_arr = np.concatenate([coords_kept, coords_arr], axis=1)
            data_arr = np.concatenate([data_kept, data_agg], axis=0)
        else:
            data_arr = data_agg

        return (
            sparse.COO(coords=coords_arr, data=data_arr, shape=matrix.shape),
            replaced_count,
        )

    max_act_idx = max(activity_mapping) if activity_mapping else -1
    max_flow_idx = max(biosphere_mapping) if biosphere_mapping else -1

    new_a_shape = (
        trails.A.shape[0],
        max(trails.A.shape[1], max_act_idx + 1),
        max(trails.A.shape[2], max_act_idx + 1),
    )
    new_b_shape = (
        trails.B.shape[0],
        max(trails.B.shape[1], max_act_idx + 1),
        max(trails.B.shape[2], max_flow_idx + 1),
    )

    trails.A = _resize_sparse(trails.A, new_a_shape)
    trails.B = _resize_sparse(trails.B, new_b_shape)

    def _apply_interpolation(
        values_by_coord: dict[tuple[int, int], dict[int, float]],
        years_all: np.ndarray,
        add_coords: list[tuple[int, int, int]],
        add_data: list[float],
    ) -> None:
        if not values_by_coord:
            return
        years_sorted = np.array(sorted(set(int(y) for y in years_all)), dtype=int)
        year_to_t = {int(y): int(trails.scenario_index[str(int(y))]) for y in years_all}

        for (act_idx, j_idx), values in values_by_coord.items():
            xs = np.array(sorted(values.keys()), dtype=int)
            ys = np.array([values[x] for x in xs], dtype=float)
            if xs.size == 1:
                interp = np.full(years_sorted.shape, ys[0], dtype=float)
            else:
                interp = np.interp(years_sorted, xs, ys)
            for year, value in zip(years_sorted, interp):
                if value == 0.0:
                    continue
                add_coords.append((year_to_t[int(year)], act_idx, j_idx))
                add_data.append(float(value))

    def _print_unlinked(rows: list[dict[str, str]]) -> None:
        try:
            from prettytable import PrettyTable

            table = PrettyTable()
            table.field_names = [
                "type",
                "name",
                "reference product",
                "categories",
                "location",
            ]
            for row in rows:
                table.add_row(
                    [
                        row.get("type", ""),
                        row.get("name", ""),
                        row.get("reference product", ""),
                        row.get("categories", ""),
                        row.get("location", ""),
                    ]
                )
            print("Unlinked exchanges:")
            print(table)
        except Exception:
            print("Unlinked exchanges detected (unable to format table).")

    if unlinked_rows:
        _print_unlinked(unlinked_rows)
        raise ValueError(
            "Unlinked exchanges detected. Fix the Excel inventory and re-import."
        )

    if apply_to_all_template_years:
        years_all = trails.years_int
        if years_all.size:
            _apply_interpolation(a_template_values, years_all, a_coords, a_data)
            _apply_interpolation(b_template_values, years_all, b_coords, b_data)

    def _drop_rows(
        matrix: sparse.COO, t_indices: set[int], act_indices: set[int]
    ) -> sparse.COO:
        if not act_indices:
            return matrix
        t_coords = matrix.coords[0]
        a_coords = matrix.coords[1]
        mask = ~(
            np.isin(t_coords, np.array(sorted(t_indices), dtype=t_coords.dtype))
            & np.isin(a_coords, np.array(sorted(act_indices), dtype=a_coords.dtype))
        )
        if mask.all():
            return matrix
        return sparse.COO(
            coords=matrix.coords[:, mask],
            data=matrix.data[mask],
            shape=matrix.shape,
        )

    trails.A = _drop_rows(trails.A, target_t, affected_acts)
    trails.B = _drop_rows(trails.B, target_t, affected_acts)

    trails.A, replaced_a = _replace_sparse_entries(trails.A, a_coords, a_data)
    trails.B, replaced_b = _replace_sparse_entries(trails.B, b_coords, b_data)

    trails._A_row_cache.clear()
    trails._direct_bio_cache_by_year.clear()
    trails._tech_td_cache.clear()
    trails._tech_td_expanded_cache.clear()
    trails._td_offsets_cache.clear()
    if hasattr(trails, "_debug_flow_filters"):
        delattr(trails, "_debug_flow_filters")

    logger.info(
        "Imported Excel inventory: prod=%d tech=%d bio=%d new_acts=%d new_flows=%d",
        production_count,
        technosphere_count,
        biosphere_count,
        new_activity_count,
        new_biosphere_count,
    )
    print(
        "Imported Excel inventory: "
        f"new_activities={new_activity_count} matched_existing_activities={len(updated_activity_ids)} "
        f"new_flows={new_biosphere_count} matched_existing_flows={len(updated_biosphere_ids)} "
        f"unlinked={unlinked_count} replaced_A={replaced_a} replaced_B={replaced_b}"
    )

    if new_activity_rows:
        try:
            from prettytable import PrettyTable

            table = PrettyTable()
            table.field_names = [
                "index",
                "name",
                "reference product",
                "unit",
                "location",
            ]
            for row in new_activity_rows:
                table.add_row(
                    [
                        row.get("index", ""),
                        row.get("name", ""),
                        row.get("reference product", ""),
                        row.get("unit", ""),
                        row.get("location", ""),
                    ]
                )
            print("New activities created from import:")
            print(table)
        except Exception:
            print("New activities created from import (unable to format table).")

    if cache_import:
        save_cached_interpolation(
            trails.package,
            value_dtype=str(trails.value_dtype),
            index_dtype=str(trails.index_dtype),
            A=trails.A,
            B=trails.B,
            scenario_labels=trails.scenario_labels,
            template_labels=trails.template_labels,
            temporal_technosphere_exchanges=trails.temporal_technosphere_exchanges,
            temporal_biosphere_exchanges=trails.temporal_biosphere_exchanges,
            activity_indices=trails.activity_indices,
            biosphere_indices=trails.biosphere_indices,
        )

    return {
        "production": production_count,
        "technosphere": technosphere_count,
        "biosphere": biosphere_count,
        "new_activities": new_activity_count,
        "new_flows": new_biosphere_count,
        "matched_existing_activities": len(updated_activity_ids),
        "matched_existing_flows": len(updated_biosphere_ids),
        "unlinked": unlinked_count,
        "scenario_indices": [int(t) for _, t in targets],
    }
