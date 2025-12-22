import trails.utils as utils


def test_parse_float_or_none() -> None:
    """Verify parsing of float values and empty inputs.

    :returns: None.
    :rtype: None
    """
    assert utils._parse_float_or_none(None) is None
    assert utils._parse_float_or_none("") is None
    assert utils._parse_float_or_none("  ") is None
    assert utils._parse_float_or_none("3.5") == 3.5


def test_parse_int_or_none() -> None:
    """Verify parsing of integer values and empty inputs.

    :returns: None.
    :rtype: None
    """
    assert utils._parse_int_or_none(None) is None
    assert utils._parse_int_or_none("") is None
    assert utils._parse_int_or_none("  ") is None
    assert utils._parse_int_or_none("3.0") == 3


def test_format_path_label() -> None:
    """Verify formatting of activity index paths.

    :returns: None.
    :rtype: None
    """
    mapping = {1: "A", 2: "B"}
    assert utils._format_path_label((1, 3, 2), mapping) == "A → Activity 3 → B"


def test_format_path_label_with_years() -> None:
    """Verify formatting of year-tagged activity paths.

    :returns: None.
    :rtype: None
    """
    mapping = {1: "Act A"}
    assert utils._format_path_label_with_years(((2005, 1), (2020, 2)), mapping) == (
        "2005: Act A → 2020: Activity 2"
    )


def test_parse_intish_or_none() -> None:
    """Verify parsing of int-like inputs.

    :returns: None.
    :rtype: None
    """
    assert utils._parse_intish_or_none(None) is None
    assert utils._parse_intish_or_none("") is None
    assert utils._parse_intish_or_none("4.0") == 4
