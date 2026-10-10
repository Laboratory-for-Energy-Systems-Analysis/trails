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
from trails.temporal_distributions import TemporalExchange
from dev.stock_vintage.check_ccgt_pilot import independent_graph_reference

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


@pytest.mark.parametrize("solver_mode", ["direct", "bw2calc"])
@pytest.mark.parametrize("attribute_to_roots", [False, True])
@pytest.mark.parametrize("depth", [2, 3])
@pytest.mark.parametrize("amount", [1.0, -2.0])
@pytest.mark.parametrize(
    "legacy,amount_source",
    [(False, "none"), (False, "port"), (False, "matrix"), (True, "none")],
)
def test_negative_reference_waste_chain_under_stock_and_legacy_routing(
    stock_package, solver_mode, attribute_to_roots, depth, amount, legacy, amount_source
):
    """Waste product signs must become positive treatment activity supplies.

    A capital supplier generates waste through a market with reference output
    -2 and a treatment with reference output -1. Both contribute biosphere
    flows. Check the waste market both at and before the routing frontier.
    Retain the historical legacy discrepancy as an explicit compatibility check.
    """
    if legacy:
        stock_package.descriptor.pop("stock_vintage")
    with make_model(stock_package) as model:
        a = np.eye(4)
        a[2, 0] = -0.1
        a[0, 1] = 0.5
        a[1, 1], a[1, 3], a[3, 3] = -2, 2, -1
        b = np.zeros((4, model.B.shape[2]))
        b[:, 0] = [2, 0.4, 0.3, 3]
        model.A = sparse.stack([sparse.COO.from_numpy(a)] * len(model.scenario_labels))
        model.B = sparse.stack([sparse.COO.from_numpy(b)] * len(model.scenario_labels))
        for label in model.activity_indices:
            model.activity_indices[label][3] = {
                "name": "synthetic negative-reference treatment",
                "reference product": "synthetic waste",
                "unit": "unit",
                "location": "GLO",
            }
        model.temporal_technosphere_exchanges = {}
        model.temporal_biosphere_exchanges = {}
        if amount_source != "none":
            for label in model.template_labels:
                for parent, child in ((0, 1), (1, 3)):
                    model.temporal_technosphere_exchanges[(label, parent, child)] = (
                        TemporalExchange(
                            distribution=6,
                            loc=None,
                            scale=None,
                            offset_min=0,
                            offset_max=0,
                            offsets=(0,),
                            weights=(1.0,),
                            amount_source=amount_source,
                        )
                    )
        direct = model.expand_temporal_exchanges(
            year=2021, act_idx=0, amount=0.1 * amount
        )
        expected_market_supply = (-0.025 if legacy else 0.025) * amount
        assert sum(row.get(1, 0) for row in direct.values()) == pytest.approx(
            expected_market_supply
        )
        model.temporal_routing(
            start_year=2021,
            start_act_idx=2,
            amount=amount,
            max_depth=depth,
            adaptive_relative_score_cutoff=None,
            min_amount=0,
            show_progress=False,
            attribute_to_roots=attribute_to_roots,
        )
        # q=.1, waste-market supply=.025 and treatment supply=.05 per FU.
        inventory = model.lci(show_progress=False, solver_mode=solver_mode)
        total = float(inventory.sum().item())
        expected = (0.34 if depth == 2 else 0.64) if legacy else 0.66
        assert total == pytest.approx(expected * amount, rel=1e-12, abs=1e-14)
        reference, audit = independent_graph_reference(
            model, model.A[0].to_scipy_sparse().T.tocsc(), 2, amount
        )
        assert reference.sum() == pytest.approx(expected * amount, rel=1e-12, abs=1e-14)
        assert audit["operator_conservation_passed"] is not legacy
        waste_nodes = [d for n, d in model.graph.nodes(data=True) if n[-1] == 1]
        assert sum(d["amount"] for d in waste_nodes) == pytest.approx(
            expected_market_supply
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


def test_successive_offsets_extend_storage_and_keep_expanded_direct_biosphere(
    stock_package,
):
    with make_model(stock_package) as model:
        data = payload()
        for record in data["profiles"][0]["years"]:
            record.update(event_years=[1900], weights=[1.0])
        model.stock_vintage_profiles = StockVintageProfiles(
            data, model.activity_indices, sha256="synthetic-cascade"
        )
        a = np.eye(model.A.shape[1])
        a[2, 0], a[0, 1] = -0.1, -0.5
        model.A = sparse.stack([sparse.COO.from_numpy(a)] * len(model.scenario_labels))
        b = np.zeros(model.B.shape)
        b[:, 2, 0], b[:, 0, 0], b[:, 1, 0] = 0.3, 2.0, 4.0
        model.B = sparse.COO.from_numpy(b)
        # A deliberately long component lead time makes the combined route
        # older than either single offset. It is a numerical boundary test.
        model.temporal_technosphere_exchanges = {
            (label, 0, 1): TemporalExchange(
                distribution=6,
                loc=None,
                scale=None,
                offset_min=-200,
                offset_max=-200,
                offsets=[-200],
                weights=[1.0],
                amount_source="port",
            )
            for label in model.template_labels
        }
        model.temporal_biosphere_exchanges = {}
        assert model.min_year + model._inventory_offset_bounds()[0] > 1700
        for year in (2022, 2020):
            model.temporal_routing(
                start_year=year,
                start_act_idx=2,
                max_depth=3,
                adaptive_relative_score_cutoff=None,
                min_amount=0,
                show_progress=False,
            )
            inventory = model.lci(show_progress=False, solver_mode="direct")
            annual = inventory.sum(dim=[d for d in inventory.dims if d != "year"])
            result = {
                int(y): v
                for y, v in zip(annual.year.values, annual.data.todense())
                if v != 0
            }
            assert result == pytest.approx({year: 0.3, 1900: 0.2, 1700: 0.2})
            assert model._score_years[0] <= 1700


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
