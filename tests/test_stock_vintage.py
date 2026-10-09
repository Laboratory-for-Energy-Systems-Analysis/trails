"""Scientific and routing checks for opt-in calendar-year stock profiles."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import sys
import types

from datapackage import Package
import numpy as np
import pytest
import sparse
from scipy.sparse.linalg import spsolve

from trails import Trails
from trails.stock_vintage import StockVintageProfiles
from trails.importer import import_excel_inventory

CALLER = {
    "name": "transport, passenger car, BEV",
    "reference product": "transport, passenger car, BEV",
    "unit": "vkm",
    "location": "GLO",
}
SUPPLIER = {
    "name": "battery electric vehicle, production",
    "reference product": "battery electric vehicle",
    "unit": "unit",
    "location": "GLO",
}


def payload():
    return {
        "schema_version": 1,
        "profiles": [
            {
                "id": "synthetic-fleet",
                "event_role": "existing_asset_service",
                "allocation_basis": "common_amortisation",
                "asset_unit": "unit",
                "service_unit": "vkm",
                "provenance": {"evidence_tier": "synthetic", "equal_utilisation": True},
                "years": [
                    {
                        "service_year": 2020,
                        "event_years": [2010, 2015],
                        "weights": [0.75, 0.25],
                    },
                    {
                        "service_year": 2021,
                        "event_years": [2010, 2021],
                        "weights": [0.4, 0.6],
                    },
                    {
                        "service_year": 2022,
                        "event_years": [2010, 2021],
                        "weights": [0.2, 0.8],
                    },
                ],
            }
        ],
        "bindings": [
            {"profile_id": "synthetic-fleet", "caller": CALLER, "supplier": SUPPLIER}
        ],
    }


@pytest.fixture
def stock_package(tmp_path):
    source = Path(__file__).parents[1] / "examples" / "example data package"
    root = tmp_path / "package"
    shutil.copytree(source, root)
    content = json.dumps(payload()).encode()
    (root / "stock_vintage.json").write_bytes(content)
    descriptor = json.loads((root / "datapackage.json").read_text())
    descriptor["resources"].append(
        {
            "name": "stock-vintage",
            "path": "stock_vintage.json",
            "profile": "data-resource",
            "format": "json",
            "mediatype": "application/json",
        }
    )
    descriptor["stock_vintage"] = {
        "schema_version": 1,
        "resource": "stock-vintage",
        "sha256": hashlib.sha256(content).hexdigest(),
    }
    (root / "datapackage.json").write_text(json.dumps(descriptor))
    return Package(str(root / "datapackage.json"))


def make_model(package, **kwargs):
    return Trails(
        package,
        stock_vintage=True,
        cache_interpolation=False,
        value_dtype=np.float64,
        **kwargs,
    )


def test_package_requires_explicit_opt_in(stock_package):
    with pytest.raises(ValueError, match="stock_vintage=True"):
        Trails(stock_package, cache_interpolation=False)


def test_annual_lookup_preserves_calendar_cohorts_and_signed_total(stock_package):
    with make_model(stock_package) as model:
        # These years share the same matrix template. Reverse/repeat order also
        # exposes cache pollution between successive routing requests.
        for year in [2022, 2020, 2021, 2020]:
            record = payload()["profiles"][0]["years"][year - 2020]
            for amount in [1.0, -2.0]:
                demand = model.expand_temporal_exchanges(
                    year=year, act_idx=2, amount=amount
                )
                actual = {
                    event: mapping[0]
                    for event, mapping in demand.items()
                    if 0 in mapping
                }
                coefficient = -float(model.A[model.scenario_index[str(year)], 2, 0])
                expected = dict(
                    zip(
                        record["event_years"],
                        [amount * coefficient * w for w in record["weights"]],
                    )
                )
                assert actual == pytest.approx(expected)
                assert sum(actual.values()) == pytest.approx(amount * coefficient)
            tex = model.get_temporal_exchange(year, 2, 0)
            assert [year + offset for offset in tex.offsets] == record["event_years"]


def test_missing_year_fails_instead_of_reusing_nearest_profile(stock_package):
    with make_model(stock_package) as model:
        with pytest.raises(ValueError, match="no service year 2023"):
            model.expand_temporal_exchanges(year=2023, act_idx=2, amount=1.0)


def test_profile_year_is_independent_of_matrix_year(stock_package):
    with make_model(stock_package, interpolate_annual=False) as model:
        assert model._map_year_to_scenario_year(2021) == 2020
        demand = model.expand_temporal_exchanges(year=2021, act_idx=2, amount=1.0)
        assert {year for year, row in demand.items() if 0 in row} == {2010, 2021}


def test_router_and_lci_conserve_static_inventory_with_annual_profiles(stock_package):
    with make_model(stock_package) as model:
        # Hold technologies and biosphere coefficients constant to isolate timing.
        a = model.A[0]
        b = model.B[0]
        model.A = sparse.stack([a] * len(model.scenario_labels))
        model.B = sparse.stack([b] * len(model.scenario_labels))
        functional_unit = np.zeros(a.shape[0])
        functional_unit[2] = 1.0
        expected_supply = spsolve(a.to_scipy_sparse().T.tocsc(), functional_unit)
        expected_biosphere = b.to_scipy_sparse().T @ expected_supply
        for year in (2020, 2021):
            model.temporal_routing(
                start_year=year,
                start_act_idx=2,
                amount=1.0,
                max_depth=1,
                adaptive_relative_score_cutoff=None,
                min_amount=0.0,
                show_progress=False,
            )
            car_nodes = [node for node in model.graph if node[-1] == 0]
            expected_years = set(
                payload()["profiles"][0]["years"][year - 2020]["event_years"]
            )
            assert {node[0] for node in car_nodes} == expected_years
            inventory = model.lci(show_progress=False, solver_mode="direct")
            totals = inventory.sum(dim=[d for d in inventory.dims if d != "flow"]).data
            assert np.allclose(
                totals.todense(), expected_biosphere, rtol=1e-9, atol=1e-10
            )


def test_stock_resource_hash_is_checked_on_warm_cache(
    stock_package, monkeypatch, tmp_path
):
    monkeypatch.setattr(
        "trails.cache_interpolation.platformdirs.user_data_path",
        lambda **kwargs: tmp_path / "cache",
    )
    with Trails(stock_package, stock_vintage=True, value_dtype=np.float64) as model:
        assert model.get_temporal_exchange(2021, 2, 0).weights == (0.4, 0.6)
    with Trails(stock_package, stock_vintage=True, value_dtype=np.float64) as model:
        assert model.get_temporal_exchange(2022, 2, 0).weights == (0.2, 0.8)
    resource = next(r for r in stock_package.resources if r.name == "stock-vintage")
    Path(resource.source).write_text(json.dumps(payload()) + " ")
    with pytest.raises(ValueError, match="hash mismatch"):
        Trails(stock_package, stock_vintage=True, value_dtype=np.float64)


@pytest.mark.parametrize(
    "weights", [[0.2, 0.3], [-0.1, 1.1], [float("nan"), 1.0], [0, 0]]
)
def test_malformed_weights_are_rejected(weights):
    data = payload()
    data["profiles"][0]["years"][0]["weights"] = weights
    with pytest.raises(ValueError, match="weights"):
        StockVintageProfiles(data, {"2020": {0: SUPPLIER, 2: CALLER}}, sha256="test")


def test_identity_units_and_unique_binding_are_required():
    indices = {"2020": {0: SUPPLIER, 2: CALLER}}
    data = payload()
    data["bindings"].append(deepcopy(data["bindings"][0]))
    with pytest.raises(ValueError, match="Conflicting"):
        StockVintageProfiles(data, indices, sha256="test")
    data = payload()
    data["profiles"][0]["service_unit"] = "passenger-kilometre"
    with pytest.raises(ValueError, match="units"):
        StockVintageProfiles(data, indices, sha256="test")
    data = payload()
    data["bindings"][0]["caller"] = dict(CALLER, location="UNKNOWN")
    with pytest.raises(ValueError, match="Unmatched"):
        StockVintageProfiles(data, indices, sha256="test")


def test_legacy_package_keeps_legacy_lookup(example_package):
    with (
        Trails(example_package, cache_interpolation=False) as legacy,
        make_model(example_package) as enabled,
    ):
        assert enabled.stock_vintage_profiles is None
        assert enabled._get_tech_temporal_exchange(
            2021, 2, 0
        ) == legacy._get_tech_temporal_exchange(2021, 2, 0)


def test_foreground_caller_can_supersede_package_profile():
    profiles = StockVintageProfiles(
        payload(), {"2020": {0: SUPPLIER, 2: CALLER}}, sha256="test"
    )
    assert profiles.get(2021, 2020, 2, 0) is not None
    profiles.remove_callers(["2020"], {2})
    assert profiles.get(2021, 2020, 2, 0) is None


def test_foreground_override_does_not_poison_base_package_cache(
    stock_package, monkeypatch, tmp_path
):
    monkeypatch.setattr(
        "trails.cache_interpolation.platformdirs.user_data_path",
        lambda **kwargs: tmp_path / "cache",
    )
    dataset = dict(
        CALLER,
        database="foreground",
        code="caller",
        exchanges=[
            dict(CALLER, type="production", amount=1.0),
            dict(
                SUPPLIER,
                type="technosphere",
                amount=0.1,
                temporal_distribution=6,
                temporal_offsets=[-1],
                temporal_weights=[1.0],
                temporal_amount_source="port",
            ),
        ],
    )

    class ExcelImporter:
        def __init__(self, path):
            self.data = [deepcopy(dataset)]

        def apply_strategies(self):
            pass

    # Stub only the workbook reader; exercise the real importer, matrices,
    # temporal lookup and persistent interpolation-cache implementation.
    for name in ("bw2io", "bw2io.importers", "bw2io.importers.excel"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    sys.modules["bw2io.importers.excel"].ExcelImporter = ExcelImporter
    workbook = tmp_path / "foreground.xlsx"
    workbook.touch()
    with Trails(stock_package, stock_vintage=True, value_dtype=np.float64) as model:
        original = float(model.A[model.scenario_index["2021"], 2, 0])
        model.expand_temporal_exchanges(year=2021, act_idx=2, amount=1.0)
        import_excel_inventory(model, workbook, cache_import=True)
        actual = model.expand_temporal_exchanges(year=2021, act_idx=2, amount=1.0)
        assert {year: values[0] for year, values in actual.items() if 0 in values} == {
            2020: pytest.approx(0.1)
        }
    with Trails(stock_package, stock_vintage=True, value_dtype=np.float64) as fresh:
        assert float(fresh.A[fresh.scenario_index["2021"], 2, 0]) == original
        assert fresh.get_temporal_exchange(2021, 2, 0).weights == (0.4, 0.6)
