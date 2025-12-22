# utils.py
from __future__ import annotations

from typing import Iterable, Mapping


def _parse_float_or_none(value: object) -> float | None:
    """Parse a value into a float or return None for empty inputs.

    :param value: Input value that may be convertible to a float.
    :type value: object
    :returns: Parsed float value or None when input is blank.
    :rtype: float | None
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return float(s)


def _parse_int_or_none(value: object) -> int | None:
    """Parse a value into an int or return None for empty inputs.

    :param value: Input value that may be convertible to an integer.
    :type value: object
    :returns: Parsed integer value or None when input is blank.
    :rtype: int | None
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return int(float(s))


def _format_path_label(path: Iterable[int], idx_to_label: Mapping[int, str]) -> str:
    """Format a path of activity indices into a label string.

    :param path: Iterable of activity indices.
    :type path: Iterable[int]
    :param idx_to_label: Mapping of activity index to display label.
    :type idx_to_label: Mapping[int, str]
    :returns: Human-readable label with arrow separators.
    :rtype: str
    """
    parts = [idx_to_label.get(idx, f"Activity {idx}") for idx in path]
    return " → ".join(parts)


def _format_path_label_with_years(
    path: Iterable[tuple[int, int]], idx_to_label: Mapping[int, str]
) -> str:
    """Format a path of (year, activity) pairs into a label string.

    :param path: Iterable of (year, activity index) pairs.
    :type path: Iterable[tuple[int, int]]
    :param idx_to_label: Mapping of activity index to display label.
    :type idx_to_label: Mapping[int, str]
    :returns: Human-readable label with year prefixes.
    :rtype: str
    """
    parts = []
    for year, act in path:
        base = idx_to_label.get(act, f"Activity {act}")
        parts.append(f"{year}: {base}")
    return " → ".join(parts)


def _parse_intish_or_none(v: object) -> int | None:
    """Parse a value into an int, allowing float-like strings.

    :param v: Input value that may represent an integer.
    :type v: object
    :returns: Parsed integer value or None when input is blank.
    :rtype: int | None
    """
    f = _parse_float_or_none(v)
    return None if f is None else int(f)
