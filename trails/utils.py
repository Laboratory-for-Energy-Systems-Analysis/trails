# utils.py
from __future__ import annotations

from typing import Iterable, Mapping


def _parse_float_or_none(value: object) -> float | None:
    """parse float or none.

    :param value: Value for `value`.
    :type value: object
    :returns: Return value.
    :rtype: float | None"""
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return float(s)


def _parse_int_or_none(value: object) -> int | None:
    """parse int or none.

    :param value: Value for `value`.
    :type value: object
    :returns: Return value.
    :rtype: int | None"""
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return int(float(s))


def _format_path_label(path: Iterable[int], idx_to_label: Mapping[int, str]) -> str:
    """format path label.

    :param path: Value for `path`.
    :type path: Iterable[int]
    :param idx_to_label: Value for `idx_to_label`.
    :type idx_to_label: Mapping[int, str]
    :returns: Return value.
    :rtype: str"""
    parts = [idx_to_label.get(idx, f"Activity {idx}") for idx in path]
    return " → ".join(parts)


def _format_path_label_with_years(
    path: Iterable[tuple[int, int]], idx_to_label: Mapping[int, str]
) -> str:
    """format path label with years.

    :param path: Value for `path`.
    :type path: Iterable[tuple[int, int]]
    :param idx_to_label: Value for `idx_to_label`.
    :type idx_to_label: Mapping[int, str]
    :returns: Return value.
    :rtype: str"""
    parts = []
    for year, act in path:
        base = idx_to_label.get(act, f"Activity {act}")
        parts.append(f"{year}: {base}")
    return " → ".join(parts)


def _parse_intish_or_none(v: object) -> int | None:
    """parse intish or none.

    :param v: Value for `v`.
    :type v: object
    :returns: Return value.
    :rtype: int | None"""
    f = _parse_float_or_none(v)
    return None if f is None else int(f)
