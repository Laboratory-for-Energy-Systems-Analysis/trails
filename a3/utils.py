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
    return int(s)
