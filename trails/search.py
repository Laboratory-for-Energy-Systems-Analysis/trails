from __future__ import annotations

from typing import Iterable, Optional

from prettytable import PrettyTable

from .trails import Trails


def _normalize_text(value: Optional[str], case_sensitive: bool) -> str:
    """normalize text.

    :param value: Value for `value`.
    :type value: Optional[str]
    :param case_sensitive: Value for `case_sensitive`.
    :type case_sensitive: bool
    :returns: Return value.
    :rtype: str"""
    text = (value or "").strip()
    return text if case_sensitive else text.lower()


def _match_text(
    value: Optional[str],
    needle: str,
    *,
    match: str,
    case_sensitive: bool,
) -> bool:
    """match text.

    :param value: Value for `value`.
    :type value: Optional[str]
    :param needle: Value for `needle`.
    :type needle: str
    :param match: Value for `match`.
    :type match: str
    :param case_sensitive: Value for `case_sensitive`.
    :type case_sensitive: bool
    :returns: Return value.
    :rtype: bool
    :raises ValueError: If an error occurs."""
    if match not in {"contains", "exact"}:
        raise ValueError("match must be 'contains' or 'exact'")

    value_norm = _normalize_text(value, case_sensitive)
    needle_norm = _normalize_text(needle, case_sensitive)

    if match == "exact":
        return value_norm == needle_norm
    return needle_norm in value_norm


def _iter_metadata(
    trails: Trails,
    *,
    kind: str,
    scenario_label: Optional[str],
) -> Iterable[tuple[int, dict]]:
    """iter metadata.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param kind: Value for `kind`.
    :type kind: str
    :param scenario_label: Value for `scenario_label`.
    :type scenario_label: Optional[str]
    :yields: Yielded values.
    :rtype: Iterable[tuple[int, dict]]
    :raises ValueError: If an error occurs."""
    if kind not in {"technosphere", "biosphere"}:
        raise ValueError("kind must be 'technosphere' or 'biosphere'")

    mapping = (
        trails.activity_indices if kind == "technosphere" else trails.biosphere_indices
    )

    if scenario_label is not None:
        for idx, meta in mapping.get(scenario_label, {}).items():
            yield int(idx), meta
        return

    for _label, meta_map in mapping.items():
        for idx, meta in meta_map.items():
            yield int(idx), meta


def search_activity(
    trails: Trails,
    query: Optional[str] = None,
    *,
    name: Optional[str] = None,
    reference_product: Optional[str] = None,
    location: Optional[str] = None,
    kind: str = "technosphere",
    scenario_label: Optional[str] = None,
    match: str = "contains",
    case_sensitive: bool = False,
) -> PrettyTable:
    """Search activity.

    :param trails: Value for `trails`.
    :type trails: Trails
    :param query: Value for `query`.
    :type query: Optional[str]
    :param name: Value for `name`.
    :type name: Optional[str]
    :param reference_product: Value for `reference_product`.
    :type reference_product: Optional[str]
    :param location: Activity location filter. Can be used alone or together
        with ``query``/``name`` and ``reference_product``.
    :type location: Optional[str]
    :param kind: Value for `kind`.
    :type kind: str
    :param scenario_label: Value for `scenario_label`.
    :type scenario_label: Optional[str]
    :param match: Value for `match`.
    :type match: str
    :param case_sensitive: Value for `case_sensitive`.
    :type case_sensitive: bool
    :returns: Return value.
    :rtype: PrettyTable
    :raises ValueError: If an error occurs."""
    needle = name if name is not None else query
    if needle is None and reference_product is None and location is None:
        raise ValueError("Provide query/name, reference_product, and/or location.")

    if kind == "biosphere" and (reference_product is not None or location is not None):
        raise ValueError(
            "reference_product and location are only valid for technosphere searches."
        )

    results: dict[int, dict] = {}
    for idx, meta in _iter_metadata(trails, kind=kind, scenario_label=scenario_label):
        if needle is not None and not _match_text(
            meta.get("name"), needle, match=match, case_sensitive=case_sensitive
        ):
            continue

        if reference_product is not None and not _match_text(
            meta.get("reference product"),
            reference_product,
            match=match,
            case_sensitive=case_sensitive,
        ):
            continue

        if location is not None and not _match_text(
            meta.get("location"),
            location,
            match=match,
            case_sensitive=case_sensitive,
        ):
            continue

        if int(idx) not in results:
            results[int(idx)] = meta

    table = PrettyTable()
    if kind == "technosphere":
        table.field_names = ["index", "name", "reference product", "location"]
        for idx in sorted(results):
            meta = results[idx]
            table.add_row(
                [
                    idx,
                    meta.get("name") or "",
                    meta.get("reference product") or "",
                    meta.get("location") or "",
                ]
            )
    else:
        table.field_names = ["index", "name", "categories"]
        for idx in sorted(results):
            meta = results[idx]
            compartment = meta.get("compartment") or ""
            subcompartment = meta.get("subcompartment") or ""
            categories = " / ".join([v for v in (compartment, subcompartment) if v])
            table.add_row([idx, meta.get("name") or "", categories])

    return table
