"""Signed biosphere amounts must survive the inventory/FaIR boundary."""

import types

import numpy as np
import pandas as pd
import pytest
import sparse
import xarray as xr

import trails.fair_rf as rf
from trails.fair_io import load_species_mapping


@pytest.mark.parametrize(
    "name",
    [
        "Carbon dioxide, in air",
        "Carbon dioxide, to soil or biomass stock",
        "Carbon dioxide, non-fossil, resource correction",
        "Carbon dioxide, fossil",
        "Carbon dioxide, non-fossil",
        "Methane, fossil",
    ],
)
@pytest.mark.parametrize("amount", [5.0, -5.0])
def test_mapping_multiplies_direction_without_discarding_supply_sign(name, amount):
    mapping, signs = load_species_mapping()
    inv = sparse.COO(np.array([[0], [0], [0]]), [amount], shape=(1, 1, 1))
    result = rf._inventory_emissions_by_fair_species(
        inv, [2035], 1, {0: name}, mapping, signs
    )
    assert result.iloc[0, 0] == amount * signs.get(name, 1.0)


def test_uptake_and_avoided_uptake_cancel_independently_of_root_grouping():
    mapping, signs = load_species_mapping()
    separate = sparse.COO(
        np.array([[0, 0], [0, 0], [0, 1]]), [5.0, -5.0], shape=(1, 1, 2)
    )
    collapsed = separate.sum(axis=2, keepdims=True)
    for inventory in (separate, collapsed):
        result = rf._inventory_emissions_by_fair_species(
            inventory, [2035], 1, {0: "Carbon dioxide, in air"}, mapping, signs
        )
        assert result.iloc[0, 0] == 0.0


@pytest.fixture
def linear_fair(monkeypatch):
    """A transparent climate-response stub to isolate mapping and attribution."""
    baseline = pd.DataFrame(
        dict(
            scenario=["s", "s"],
            region=["World", "World"],
            variable=["CO2 FFI", "CO2 AFOLU"],
            unit=["Gt CO2/yr", "Gt CO2/yr"],
            **{"2000.5": [0.0, 0.0], "2001.5": [0.0, 0.0], "2002.5": [0.0, 0.0]},
        )
    )
    monkeypatch.setattr(rf, "load_emissions_csv", lambda *a, **k: baseline.copy())

    def run(emissions, scenario, **kwargs):
        configs = kwargs.get("config_names") or ["test"]
        values = emissions[["2000.5", "2001.5", "2002.5"]].sum().to_numpy(float)
        response = np.repeat(np.cumsum(values)[None, None, :], len(configs), axis=1)
        coords = dict(
            scenario=[scenario], config=configs, timebounds=[2000.5, 2001.5, 2002.5]
        )
        return types.SimpleNamespace(
            forcing=xr.DataArray(
                response[..., None],
                dims=(*coords, "specie"),
                coords={**coords, "specie": ["CO2"]},
            ),
            temperature=xr.DataArray(response * 2, dims=tuple(coords), coords=coords),
            configs=configs,
        )

    monkeypatch.setattr(rf, "_run_fair_emissions", run)


@pytest.mark.parametrize("per_species", [False, True])
@pytest.mark.parametrize("name", ["Carbon dioxide, in air", "Carbon dioxide, fossil"])
@pytest.mark.parametrize("amounts", [(5.0, 2.0), (-5.0, -2.0), (-5.0, 2.0)])
def test_public_rf_and_temperature_paths_preserve_signed_inventory(
    linear_fair, per_species, name, amounts
):
    uptake = name == "Carbon dioxide, in air"
    direction = -1.0 if uptake else 1.0
    inventory = sparse.COO(
        np.array([[0, 0], [0, 0], [0, 1], [0, 0]]), amounts, shape=(1, 1, 3, 1)
    )
    model = types.SimpleNamespace(
        debug=False,
        inventory=xr.DataArray(
            inventory,
            dims=("activity", "flow", "year", "root activity"),
            coords={
                "activity": [0],
                "flow": [0],
                "year": [2000, 2001, 2002],
                "root activity": [0],
            },
        ),
        biosphere_indices={
            "2000": {
                0: dict(
                    name=name,
                    compartment="natural resource" if uptake else "air",
                    subcompartment="in air" if uptake else "",
                )
            }
        },
    )
    rf.run_fair_delta_rf(
        model,
        scenario="s",
        config_names=["test"],
        per_species_runs=per_species,
        scale_factor=1.0,
        quantiles=[50.0],
        validate_emissions_delta=True,
        validate_raise=True,
        show_progress=False,
    )
    expected = np.cumsum([*amounts, 0.0]) * direction / 1e12
    for attribute, scale in [
        ("instant_radiative_forcing", 1.0),
        ("delta_temperature", 2.0),
    ]:
        result = (
            getattr(model, attribute)
            .sum(dim=["flow", "root activity"])
            .data.todense()[0]
        )
        np.testing.assert_allclose(result, expected * scale, rtol=1e-10, atol=1e-24)
