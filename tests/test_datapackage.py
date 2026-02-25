import numpy as np
from datapackage import Package

import trails.datapackage as datapackage


def test_parse_intish_or_none() -> None:
    """Verify int-like parsing helper.

    :returns: None.
    :rtype: None
    """
    assert datapackage._parse_intish_or_none(None) is None
    assert datapackage._parse_intish_or_none("") is None
    assert datapackage._parse_intish_or_none("3.0") == 3


def test_iter_inventory_resources(example_package: Package) -> None:
    """Verify inventory resource iteration.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: None.
    :rtype: None
    """
    resources = list(
        datapackage._iter_inventory_resources(example_package, "A_matrix.csv")
    )
    years = sorted({year for year, _ in resources})
    assert years == ["2005", "2020", "2050", "2100"]


def test_label_to_year_and_sorting() -> None:
    """Verify label-to-year parsing and sorting.

    :returns: None.
    :rtype: None
    """
    assert datapackage._label_to_year("2050") == 2050
    assert datapackage._label_to_year("model/pathway/2005") == 2005

    years, order = datapackage._years_and_sorted_indices(["2050", "2005", "2020"])
    assert years.tolist() == [2005, 2020, 2050]
    assert order.tolist() == [1, 2, 0]


def test_load_matrices_from_package(example_package: Package) -> None:
    """Verify matrix loading from datapackage.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: None.
    :rtype: None
    """
    A, B, labels, scenario_index, temporal_exchanges, temporal_bio = (
        datapackage.load_matrices_from_package(example_package)
    )
    assert labels == ["2005", "2020", "2050", "2100"]
    assert A.shape[0] == 4
    assert B.shape[0] == 4
    assert np.isclose(float(A[scenario_index["2005"], 1, 5]), -0.7)
    assert ("2005", 2, 0) in temporal_exchanges
    assert isinstance(temporal_bio, dict)
    assert all(hasattr(v, "distribution") for v in temporal_bio.values())


def test_interpolate_to_annual(example_package: Package) -> None:
    """Verify annual interpolation of matrices.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: None.
    :rtype: None
    """
    A, B, labels, _, _, _ = datapackage.load_matrices_from_package(example_package)
    A_interp, B_interp, new_labels, new_index = datapackage.interpolate_to_annual(
        A, B, labels
    )

    assert new_labels[0] == "2004"
    assert new_labels[-1] == "2101"
    assert len(new_labels) == (2100 - 2005 + 1) + 2
    assert np.isclose(float(A_interp[new_index["2004"], 1, 5]), -0.7)
    value_2010 = float(A_interp[new_index["2010"], 1, 5])
    assert np.isclose(value_2010, -0.6333333, atol=1e-6)
    assert np.isclose(
        float(A_interp[new_index["2101"], 1, 5]),
        float(A_interp[new_index["2100"], 1, 5]),
    )
    assert np.isclose(
        float(A_interp[new_index["2101"], 1, 5]),
        float(A_interp[new_index["2100"], 1, 5]),
    )
    assert B_interp.shape[0] == len(new_labels)


def test_interpolate_to_annual_with_custom_offsets(example_package: Package) -> None:
    """Verify annual interpolation with custom boundary offsets.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: None.
    :rtype: None
    """
    A, B, labels, _, _, _ = datapackage.load_matrices_from_package(example_package)
    A_interp, _, new_labels, new_index = datapackage.interpolate_to_annual(
        A, B, labels, start_year_offset=-20, end_year_offset=20
    )

    assert new_labels[0] == "1985"
    assert new_labels[-1] == "2120"
    assert len(new_labels) == (2100 - 2005 + 1) + 40
    assert np.isclose(float(A_interp[new_index["1985"], 1, 5]), -0.7)
    assert np.isclose(
        float(A_interp[new_index["2120"], 1, 5]),
        float(A_interp[new_index["2100"], 1, 5]),
    )


def test_load_indices_from_package(example_package: Package) -> None:
    """Verify index loading from datapackage.

    :param example_package: Example datapackage fixture.
    :type example_package: datapackage.Package
    :returns: None.
    :rtype: None
    """
    activity_indices, biosphere_indices = datapackage.load_indices_from_package(
        example_package
    )
    assert activity_indices["2005"][0]["name"] == "battery electric vehicle, production"
    assert biosphere_indices["2005"][0]["name"] == "Carbon dioxide, fossil"


def test_parse_temporal_exchange_row_with_explicit_pulses() -> None:
    """Verify parsing of JSON pulse lists for temporal distributions.

    :returns: None.
    :rtype: None
    """
    row = {
        "temporal_distribution": "6",
        "temporal_loc": "",
        "temporal_scale": "",
        "temporal_min": "0",
        "temporal_max": "10",
        "temporal_amount_source": "port",
        "temporal_offsets": "[0, 5, 10]",
        "temporal_weights": "[0.5, 0.25, 0.25]",
    }
    tex = datapackage._parse_temporal_exchange_row(row)
    assert tex is not None
    assert tex.distribution == 6
    assert tex.offsets == [0, 5, 10]
    assert tex.weights == [0.5, 0.25, 0.25]
