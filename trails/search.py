from __future__ import annotations

from typing import Iterable, Optional

from prettytable import PrettyTable

from .trails import Trails


def _normalize_text(value: Optional[str], case_sensitive: bool) -> str:
    text = (value or "").strip()
    return text if case_sensitive else text.lower()


def _match_text(
    value: Optional[str],
    needle: str,
    *,
    match: str,
    case_sensitive: bool,
) -> bool:
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
    kind: str = "technosphere",
    scenario_label: Optional[str] = None,
    match: str = "contains",
    case_sensitive: bool = False,
) -> PrettyTable:
    """Search activity or biosphere indices by metadata fields.

    :param trails: Trails instance with indices loaded.
    :param query: Default text to match against the activity name.
    :param name: Optional name string to match (overrides query if provided).
    :param reference_product: Optional reference product filter (technosphere only).
    :param kind: Either ``technosphere`` or ``biosphere``.
    :param scenario_label: Optional metadata label to restrict search.
    :param match: ``contains`` or ``exact`` string matching.
    :param case_sensitive: Whether matching is case-sensitive.
    :returns: PrettyTable of matching indices and metadata.
    """
    needle = name if name is not None else query
    if needle is None and reference_product is None:
        raise ValueError("Provide query/name and/or reference_product.")

    if kind == "biosphere" and reference_product is not None:
        raise ValueError("reference_product is only valid for technosphere searches.")

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
