# utils.py


def _parse_float_or_none(value):
    """Parse a value into a float or return None for empty inputs.

    Args:
        value: Input value that may be convertible to a float.

    Returns:
        float | None: Parsed float value or None when input is blank.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return float(s)


def _parse_int_or_none(value):
    """Parse a value into an int or return None for empty inputs.

    Args:
        value: Input value that may be convertible to an integer.

    Returns:
        int | None: Parsed integer value or None when input is blank.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return int(float(s))


def _format_path_label(path, idx_to_label):
    """Format a path of activity indices into a label string.

    Args:
        path: Iterable of activity indices.
        idx_to_label: Mapping of activity index to display label.

    Returns:
        str: Human-readable label with arrow separators.
    """
    parts = [idx_to_label.get(idx, f"Activity {idx}") for idx in path]
    return " → ".join(parts)


def _format_path_label_with_years(path, idx_to_label):
    """Format a path of (year, activity) pairs into a label string.

    Args:
        path: Iterable of (year, activity index) pairs.
        idx_to_label: Mapping of activity index to display label.

    Returns:
        str: Human-readable label with year prefixes.
    """
    parts = []
    for year, act in path:
        base = idx_to_label.get(act, f"Activity {act}")
        parts.append(f"{year}: {base}")
    return " → ".join(parts)


def _parse_intish_or_none(v):
    """Parse a value into an int, allowing float-like strings.

    Args:
        v: Input value that may represent an integer.

    Returns:
        int | None: Parsed integer value or None when input is blank.
    """
    f = _parse_float_or_none(v)
    return None if f is None else int(f)
