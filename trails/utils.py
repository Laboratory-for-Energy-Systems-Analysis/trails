# utils.py

def _parse_float_or_none(value):
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return float(s)


def _parse_int_or_none(value):
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    return int(float(s))


def _format_path_label(path, idx_to_label):
    """Turn a tuple of activity indices into a nice string."""
    parts = [idx_to_label.get(idx, f"Activity {idx}") for idx in path]
    return " → ".join(parts)

def _format_path_label_with_years(path, idx_to_label):
    """
    Turn a tuple of (year, act_idx) pairs into a nice string.
    """
    parts = []
    for year, act in path:
        base = idx_to_label.get(act, f"Activity {act}")
        parts.append(f"{year}: {base}")
    return " → ".join(parts)

def _parse_intish_or_none(v):
    f = _parse_float_or_none(v)
    return None if f is None else int(f)
