import numpy as np

import trails.datapackage as datapackage


def test_parse_intish_or_none():
    assert datapackage._parse_intish_or_none(None) is None
    assert datapackage._parse_intish_or_none("") is None
    assert datapackage._parse_intish_or_none("3.0") == 3


def test_iter_inventory_resources(example_package):
    resources = list(
        datapackage._iter_inventory_resources(example_package, "A_matrix.csv")
    )
    years = sorted({year for year, _ in resources})
    assert years == ["2005", "2020", "2050", "2100"]


def test_label_to_year_and_sorting():
    assert datapackage._label_to_year("2050") == 2050
    assert datapackage._label_to_year("model/pathway/2005") == 2005

    years, order = datapackage._years_and_sorted_indices(["2050", "2005", "2020"])
    assert years.tolist() == [2005, 2020, 2050]
    assert order.tolist() == [1, 2, 0]


def test_load_matrices_from_package(example_package):
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


def test_interpolate_to_annual(example_package):
    A, B, labels, _, _, _ = datapackage.load_matrices_from_package(example_package)
    A_interp, B_interp, new_labels, new_index = datapackage.interpolate_to_annual(
        A, B, labels
    )

    assert new_labels[0] == "2005"
    assert new_labels[-1] == "2100"
    assert len(new_labels) == (2100 - 2005 + 1)
    value_2010 = float(A_interp[new_index["2010"], 1, 5])
    assert np.isclose(value_2010, -0.6333333, atol=1e-6)
    assert B_interp.shape[0] == len(new_labels)


def test_load_indices_from_package(example_package):
    activity_indices, biosphere_indices = datapackage.load_indices_from_package(
        example_package
    )
    assert activity_indices["2005"][0]["name"] == "battery electric vehicle, production"
    assert biosphere_indices["2005"][0]["name"] == "Carbon dioxide, fossil"
