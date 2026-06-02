from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import logging
import json
from typing import TYPE_CHECKING

import numpy as np
import sparse

from .temporal_distributions import TemporalExchange, resolve_temporal_offset_bounds
from .utils import _parse_float_or_none, _parse_intish_or_none
from .cache_interpolation import save_cached_interpolation

if TYPE_CHECKING:
    from .trails import Trails

logger = logging.getLogger(__name__)


def import_excel_inventory(
    trails: "Trails",
    path: str | Path | list[str | Path],
    *,
    year: int | None = None,
    scenario_label: str | None = None,
    cache_import: bool = False,
) -> dict[str, int]:
    """Import excel inventory.

    :param trails: Value for `trails`.
    :type trails: 'Trails'
    :param path: Value for `path`.
    :type path: str | Path | list[str | Path]
    :param year: Value for `year`.
    :type year: int | None
    :param scenario_label: Value for `scenario_label`.
    :type scenario_label: str | None
    :param cache_import: Value for `cache_import`.
    :type cache_import: bool
    :returns: Return value.
    :rtype: dict[str, int]
    :raises FileNotFoundError: If an error occurs.
    :raises RuntimeError: If an error occurs.
    :raises ValueError: If an error occurs."""
    if trails.A is None or trails.B is None:
        raise RuntimeError("A/B matrices are not initialized.")

    from bw2io.importers.excel import ExcelImporter

    paths: list[Path]
    if isinstance(path, (str, Path)):
        paths = [Path(path)]
    elif isinstance(path, Sequence):
        paths = [Path(p) for p in path]
    else:
        raise TypeError("path must be a filepath or list of filepaths")

    if not paths:
        raise ValueError("path list is empty")

    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Excel inventory file not found: {p}")

    datasets: list[dict] = []
    for p in paths:
        importer = ExcelImporter(str(p))
        # Keep year-specific columns (e.g., 2010, 2020) by skipping csv_drop_unknown.
        strategies = getattr(importer, "strategies", None)
        if strategies:
            strategies = [
                s
                for s in strategies
                if getattr(s, "__name__", "") != "csv_drop_unknown"
            ]
            importer.apply_strategies(strategies=strategies)
        else:
            importer.apply_strategies()
        datasets.extend(getattr(importer, "data", []) or [])

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
        """norm.

        :param value: Value for `value`.
        :type value: object
        :returns: Return value.
        :rtype: str"""
        return "" if value is None else str(value).strip()

    def _parse_amount_source(value: object) -> str:
        """parse amount source.

        :param value: Value for `value`.
        :type value: object
        :returns: Return value.
        :rtype: str
        :raises ValueError: If an error occurs."""
        if value is None:
            return "port"
        s = str(value).strip().lower()
        if s == "":
            return "port"
        if s in {"port", "matrix"}:
            return s
        raise ValueError(f"Invalid temporal_amount_source: {value!r}")

    def _parse_temporal_exchange(exchange: dict) -> TemporalExchange | None:
        """parse temporal exchange.

        :param exchange: Value for `exchange`.
        :type exchange: dict
        :returns: Return value.
        :rtype: TemporalExchange | None"""
        dist_code = _parse_intish_or_none(exchange.get("temporal_distribution"))
        if dist_code in (None, 0):
            return None

        loc = _parse_float_or_none(exchange.get("temporal_loc"))
        scale = _parse_float_or_none(exchange.get("temporal_scale"))
        off_min, off_max = resolve_temporal_offset_bounds(
            distribution=int(dist_code),
            loc=loc,
            offset_min=_parse_intish_or_none(exchange.get("temporal_min")),
            offset_max=_parse_intish_or_none(exchange.get("temporal_max")),
        )
        amount_source = _parse_amount_source(exchange.get("temporal_amount_source"))
        offsets = _parse_json_number_list(
            exchange.get("temporal_offsets"), integer=True
        )
        weights = _parse_json_number_list(
            exchange.get("temporal_weights"), integer=False
        )

        return TemporalExchange(
            distribution=int(dist_code),
            loc=loc,
            scale=scale,
            offset_min=int(off_min),
            offset_max=int(off_max),
            amount_source=amount_source,
            offsets=offsets,
            weights=weights,
        )

    def _parse_json_number_list(
        value: object, *, integer: bool
    ) -> list[int] | list[float] | None:
        """parse json number list.

        :param value: Value for `value`.
        :type value: object
        :param integer: Value for `integer`.
        :type integer: bool
        :returns: Return value.
        :rtype: list[int] | list[float] | None"""
        if value is None:
            return None
        if isinstance(value, str):
            s = value.strip()
            if s == "":
                return None
            try:
                parsed = json.loads(s)
            except json.JSONDecodeError:
                return None
        elif isinstance(value, (list, tuple, np.ndarray)):
            parsed = value
        else:
            return None

        if not isinstance(parsed, (list, tuple)):
            return None

        out: list[int] | list[float] = []
        for item in parsed:
            try:
                num = float(item)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(num):
                continue
            if integer:
                out.append(int(round(num)))
            else:
                out.append(float(num))
        return out if out else None

    def _extract_year_amounts(exchange: dict) -> dict[int, float]:
        """extract year amounts.

        :param exchange: Value for `exchange`.
        :type exchange: dict
        :returns: Return value.
        :rtype: dict[int, float]"""
        year_amounts: dict[int, float] = {}
        for key, value in exchange.items():
            year_key = None
            if isinstance(key, (int, np.integer)):
                year_key = int(key)
            elif isinstance(key, (float, np.floating)):
                if float(key).is_integer():
                    year_key = int(round(float(key)))
            elif isinstance(key, str):
                k = key.strip()
                if k.isdigit():
                    year_key = int(k)
                else:
                    try:
                        fk = float(k)
                    except (TypeError, ValueError):
                        fk = None
                    if fk is not None and float(fk).is_integer():
                        year_key = int(round(float(fk)))
            if year_key is None:
                continue
            try:
                val = float(value)
            except (TypeError, ValueError):
                continue
            year_amounts[int(year_key)] = val
        return year_amounts

    def _amount_for_year(
        year_int: int | None, year_amounts: dict[int, float], fallback: float
    ) -> float | None:
        """amount for year.

        :param year_int: Value for `year_int`.
        :type year_int: int | None
        :param year_amounts: Value for `year_amounts`.
        :type year_amounts: dict[int, float]
        :param fallback: Value for `fallback`.
        :type fallback: float
        :returns: Return value.
        :rtype: float | None"""
        if year_amounts:
            if year_int is None:
                return None
            if year_int in year_amounts:
                return float(year_amounts[year_int])
            return None
        return float(fallback)

    def _label_to_year(label: str, t_idx: int) -> int | None:
        """label to year.

        :param label: Value for `label`.
        :type label: str
        :param t_idx: Value for `t_idx`.
        :type t_idx: int
        :returns: Return value.
        :rtype: int | None"""
        if label.isdigit():
            return int(label)
        try:
            scen_label = str(trails.scenario_labels[int(t_idx)])
        except Exception:
            return None
        if scen_label.isdigit():
            return int(scen_label)
        return None

    def _activity_key(
        name: object,
        ref_product: object,
        location: object,
    ) -> tuple[str, str, str]:
        """activity key.

        :param name: Value for `name`.
        :type name: object
        :param ref_product: Value for `ref_product`.
        :type ref_product: object
        :param location: Value for `location`.
        :type location: object
        :returns: Return value.
        :rtype: tuple[str, str, str]"""
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
        """biosphere key.

        :param name: Value for `name`.
        :type name: object
        :param compartment: Value for `compartment`.
        :type compartment: object
        :param subcompartment: Value for `subcompartment`.
        :type subcompartment: object
        :returns: Return value.
        :rtype: tuple[str, str, str]"""
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
        """sync indices across labels.

        :raises ValueError: If an error occurs."""
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
        """next index.

        :param mapping: Value for `mapping`.
        :type mapping: dict[int, dict]
        :returns: Return value.
        :rtype: int"""
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
        """ensure activity index.

        :param name: Value for `name`.
        :type name: object
        :param ref_product: Value for `ref_product`.
        :type ref_product: object
        :param location: Value for `location`.
        :type location: object
        :param unit: Value for `unit`.
        :type unit: object | None
        :returns: Return value.
        :rtype: tuple[int, bool]"""
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
        """ensure biosphere index.

        :param name: Value for `name`.
        :type name: object
        :param compartment: Value for `compartment`.
        :type compartment: object
        :param subcompartment: Value for `subcompartment`.
        :type subcompartment: object
        :param unit: Value for `unit`.
        :type unit: object | None
        :returns: Return value.
        :rtype: tuple[int, bool]"""
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
        """activity index for dataset.

        :param dataset: Value for `dataset`.
        :type dataset: dict
        :returns: Return value.
        :rtype: int"""
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
        """biosphere index for exchange.

        :param exchange: Value for `exchange`.
        :type exchange: dict
        :returns: Return value.
        :rtype: int"""
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

    def _validate_production_exchange_metadata(
        dataset: dict,
        exchange: dict,
    ) -> None:
        """Ensure explicit production exchange metadata matches its dataset."""
        checks = [
            ("name", exchange.get("name"), dataset.get("name")),
            (
                "reference product",
                exchange.get("reference product") or exchange.get("product"),
                dataset.get("reference product") or dataset.get("product"),
            ),
            ("location", exchange.get("location"), dataset.get("location")),
            ("unit", exchange.get("unit"), dataset.get("unit")),
        ]
        mismatches = [
            f"{field}: exchange={_norm(exchange_value)!r} dataset={_norm(dataset_value)!r}"
            for field, exchange_value, dataset_value in checks
            if _norm(exchange_value) and _norm(exchange_value) != _norm(dataset_value)
        ]
        if mismatches:
            raise ValueError(
                "Production exchange metadata must match dataset header fields: "
                + "; ".join(mismatches)
            )

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
            year_amounts = _extract_year_amounts(exchange)
            try:
                amount = float(exchange.get("amount", 0.0))
            except (TypeError, ValueError):
                continue
            if not year_amounts and amount == 0.0:
                continue

            if ex_type == "production":
                _validate_production_exchange_metadata(dataset, exchange)
                prod_name = exchange.get("name") or dataset.get("name")
                prod_ref_product = (
                    exchange.get("reference product")
                    or exchange.get("product")
                    or dataset.get("reference product")
                    or dataset.get("product")
                )
                prod_location = exchange.get("location") or dataset.get("location")
                key = _activity_key(
                    prod_name,
                    prod_ref_product,
                    prod_location,
                )
                prod_idx = activity_lookup.get(key) or dataset_act_indices.get(key)
                if prod_idx is None:
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(prod_name),
                            "reference product": _norm(prod_ref_product),
                            "categories": "",
                            "location": _norm(prod_location),
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

            # If year-specific amounts are provided, store them directly for interpolation.
            if year_amounts:
                if ex_type in {"production", "technosphere"}:
                    for yk, val in year_amounts.items():
                        stored_val = val if ex_type == "production" else -val
                        a_template_values.setdefault(tech_pair, {})[int(yk)] = float(
                            stored_val
                        )
                elif ex_type == "biosphere":
                    for yk, val in year_amounts.items():
                        b_template_values.setdefault((act_idx, flow_idx), {})[
                            int(yk)
                        ] = float(val)

            if ex_type in {"production", "technosphere"}:
                check_name = exchange.get("name")
                check_ref_product = exchange.get("reference product") or exchange.get(
                    "product"
                )
                check_location = exchange.get("location")
                if ex_type == "production":
                    check_name = check_name or dataset.get("name")
                    check_ref_product = (
                        check_ref_product
                        or dataset.get("reference product")
                        or dataset.get("product")
                    )
                    check_location = check_location or dataset.get("location")
                if (
                    _norm(check_name) == ""
                    or _norm(check_ref_product) == ""
                    or _norm(check_location) == ""
                ):
                    unlinked_rows.append(
                        {
                            "type": ex_type,
                            "name": _norm(check_name),
                            "reference product": _norm(check_ref_product),
                            "categories": "",
                            "location": _norm(check_location),
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
                    year_int = _label_to_year(label, t)
                    if not apply_to_all_template_years:
                        amt_for_year = _amount_for_year(year_int, year_amounts, amount)
                        if amt_for_year is not None and amt_for_year != 0.0:
                            stored_val = (
                                amt_for_year
                                if ex_type == "production"
                                else -amt_for_year
                            )
                            a_coords.append((t, act_idx, prod_idx))
                            a_data.append(float(stored_val))

                    tex = _parse_temporal_exchange(exchange)
                    key = (label, int(act_idx), int(prod_idx))
                    if tex is None:
                        trails.temporal_technosphere_exchanges.pop(key, None)
                    else:
                        trails.temporal_technosphere_exchanges[key] = tex

                    year_int = _label_to_year(label, t)
                    if year_int is not None and not year_amounts:
                        amt_for_year = _amount_for_year(year_int, year_amounts, amount)
                        if amt_for_year is not None:
                            stored_val = (
                                amt_for_year
                                if ex_type == "production"
                                else -amt_for_year
                            )
                            a_template_values.setdefault(tech_pair, {})[year_int] = (
                                float(stored_val)
                            )
                elif ex_type == "biosphere":
                    year_int = _label_to_year(label, t)
                    if not apply_to_all_template_years:
                        amt_for_year = _amount_for_year(year_int, year_amounts, amount)
                        if amt_for_year is not None and amt_for_year != 0.0:
                            b_coords.append((t, act_idx, flow_idx))
                            b_data.append(float(amt_for_year))

                    tex = _parse_temporal_exchange(exchange)
                    key = (label, int(act_idx), int(flow_idx))
                    if tex is None:
                        trails.temporal_biosphere_exchanges.pop(key, None)
                    else:
                        trails.temporal_biosphere_exchanges[key] = tex

                    year_int = _label_to_year(label, t)
                    if year_int is not None and not year_amounts:
                        amt_for_year = _amount_for_year(year_int, year_amounts, amount)
                        if amt_for_year is not None:
                            b_template_values.setdefault((act_idx, flow_idx), {})[
                                year_int
                            ] = float(amt_for_year)

    def _resize_sparse(
        matrix: sparse.COO, new_shape: tuple[int, int, int]
    ) -> sparse.COO:
        """resize sparse.

        :param matrix: Value for `matrix`.
        :type matrix: sparse.COO
        :param new_shape: Value for `new_shape`.
        :type new_shape: tuple[int, int, int]
        :returns: Return value.
        :rtype: sparse.COO"""
        if matrix.shape == new_shape:
            return matrix
        return sparse.COO(coords=matrix.coords, data=matrix.data, shape=new_shape)

    def _append_sparse_entries(
        matrix: sparse.COO,
        coords: list[tuple[int, int, int]],
        data: list[float],
    ) -> sparse.COO:
        """append sparse entries.

        :param matrix: Value for `matrix`.
        :type matrix: sparse.COO
        :param coords: Value for `coords`.
        :type coords: list[tuple[int, int, int]]
        :param data: Value for `data`.
        :type data: list[float]
        :returns: Return value.
        :rtype: sparse.COO"""
        if not coords:
            return matrix

        coords_arr = np.array(coords, dtype=trails.index_dtype).T
        data_arr = np.array(data, dtype=trails.value_dtype)

        coords_all = np.concatenate([matrix.coords, coords_arr], axis=1)
        data_all = np.concatenate([matrix.data, data_arr], axis=0)

        return sparse.COO(coords=coords_all, data=data_all, shape=matrix.shape)

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
        """apply interpolation.

        :param values_by_coord: Value for `values_by_coord`.
        :type values_by_coord: dict[tuple[int, int], dict[int, float]]
        :param years_all: Value for `years_all`.
        :type years_all: np.ndarray
        :param add_coords: Value for `add_coords`.
        :type add_coords: list[tuple[int, int, int]]
        :param add_data: Value for `add_data`.
        :type add_data: list[float]"""
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
                interp = np.interp(
                    years_sorted,
                    xs,
                    ys,
                    left=float(ys[0]),
                    right=float(ys[-1]),
                )
            for year, value in zip(years_sorted, interp):
                if value == 0.0:
                    continue
                add_coords.append((year_to_t[int(year)], act_idx, j_idx))
                add_data.append(float(value))

    def _print_unlinked(rows: list[dict[str, str]]) -> None:
        """print unlinked.

        :param rows: Value for `rows`.
        :type rows: list[dict[str, str]]"""
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
        """drop rows.

        :param matrix: Value for `matrix`.
        :type matrix: sparse.COO
        :param t_indices: Value for `t_indices`.
        :type t_indices: set[int]
        :param act_indices: Value for `act_indices`.
        :type act_indices: set[int]
        :returns: Return value.
        :rtype: sparse.COO"""
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

    if a_coords:
        flat = np.ravel_multi_index(
            np.array(a_coords, dtype=trails.index_dtype).T, trails.A.shape
        )
        if np.unique(flat).size != flat.size:
            raise ValueError(
                "Duplicate technosphere exchanges detected in import; "
                "resolve duplicates before importing."
            )
    if b_coords:
        flat = np.ravel_multi_index(
            np.array(b_coords, dtype=trails.index_dtype).T, trails.B.shape
        )
        if np.unique(flat).size != flat.size:
            raise ValueError(
                "Duplicate biosphere exchanges detected in import; "
                "resolve duplicates before importing."
            )

    trails.A = _append_sparse_entries(trails.A, a_coords, a_data)
    trails.B = _append_sparse_entries(trails.B, b_coords, b_data)

    trails._A_row_cache.clear()
    if hasattr(trails, "_production_amount_cache"):
        trails._production_amount_cache.clear()
    trails._direct_bio_cache_by_year.clear()
    trails._tech_td_cache.clear()
    trails._tech_td_expanded_cache.clear()
    trails._td_offsets_cache.clear()
    for cache_name in (
        "_B_csr_cache",
        "_B_row_cache",
        "_B_row_index_map",
        "_B_cf_actvec_cache",
        "_bio_score_row_char_cache",
        "_bio_score_row_char_matrix_cache",
    ):
        cache = getattr(trails, cache_name, None)
        if hasattr(cache, "clear"):
            cache.clear()
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
        f"unlinked={unlinked_count} replaced_A=rows replaced_B=rows"
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
            interpolation_start_year_offset=int(
                getattr(trails, "interpolation_start_year_offset", -1)
            ),
            interpolation_end_year_offset=int(
                getattr(trails, "interpolation_end_year_offset", 1)
            ),
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
        "replaced_A": "rows",
        "replaced_B": "rows",
        "scenario_indices": [int(t) for _, t in targets],
    }
