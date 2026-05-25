from __future__ import annotations

from typing import Any
import types

import numpy as np
import pandas as pd
import pytest
import sparse
import xarray as xr

from trails.fair_rf import (
    calculate_co2_pulse_equivalents,
    integrate_window,
    make_reference_co2_pulse_emissions,
    run_fair_co2_pulse_equivalents,
)


def test_co2_pulse_equivalents_equal_reference_response() -> None:
    years = np.array([2050.0, 2051.0, 2052.0])
    d_rf_ref = np.array([1.0, 2.0, 3.0])
    d_t_ref = np.array([0.1, 0.2, 0.3])

    out = calculate_co2_pulse_equivalents(
        d_rf_ref,
        d_t_ref,
        d_rf_ref,
        d_t_ref,
        years,
        reference_pulse_mass=1.0e9,
        start_year=2050,
        end_year=2052,
    )

    assert out["co2_pulse_equivalent_integrated_rf"] == pytest.approx(1.0e9)
    assert out["co2_pulse_equivalent_integrated_temperature"] == pytest.approx(1.0e9)


def test_co2_pulse_equivalents_scale_linearly() -> None:
    years = np.array([2050.0, 2051.0, 2052.0])
    d_rf_ref = np.array([1.0, 2.0, 3.0])
    d_t_ref = np.array([0.1, 0.2, 0.3])

    out = calculate_co2_pulse_equivalents(
        2.0 * d_rf_ref,
        2.0 * d_t_ref,
        d_rf_ref,
        d_t_ref,
        years,
        reference_pulse_mass=1.0e9,
        start_year=2050,
        end_year=2052,
    )

    assert out["co2_pulse_equivalent_integrated_rf"] == pytest.approx(2.0e9)
    assert out["co2_pulse_equivalent_integrated_temperature"] == pytest.approx(2.0e9)


def test_co2_pulse_equivalents_rf_and_temperature_can_differ() -> None:
    years = np.array([2050.0, 2051.0, 2052.0])

    out = calculate_co2_pulse_equivalents(
        np.array([1.0, 1.0, 1.0]),
        np.array([2.0, 2.0, 2.0]),
        np.array([1.0, 1.0, 1.0]),
        np.array([1.0, 1.0, 1.0]),
        years,
        reference_pulse_mass=1.0,
        start_year=2050,
        end_year=2052,
    )

    assert out["co2_pulse_equivalent_integrated_rf"] == pytest.approx(1.0)
    assert out["co2_pulse_equivalent_integrated_temperature"] == pytest.approx(2.0)


def test_integrate_window_truncates_to_selected_years() -> None:
    years = np.array([2049.0, 2050.0, 2090.0, 2100.0, 2101.0])
    values = np.array([100.0, 1.0, 1.0, 3.0, 100.0])

    out = integrate_window(values, years, start_year=2050, end_year=2100)

    assert out == pytest.approx(60.0)


def test_co2_pulse_equivalents_zero_reference_integral_raises() -> None:
    years = np.array([2050.0, 2051.0, 2052.0])

    with pytest.raises(ZeroDivisionError, match="reference CO2 pulse"):
        calculate_co2_pulse_equivalents(
            np.ones(3),
            np.ones(3),
            np.zeros(3),
            np.ones(3),
            years,
            reference_pulse_mass=1.0,
        )


def test_co2_pulse_equivalents_reference_pulse_size_scaling() -> None:
    years = np.array([2050.0, 2051.0, 2052.0])
    kg_ref_rf = np.array([1.0, 2.0, 3.0])
    kg_ref_t = np.array([0.1, 0.2, 0.3])
    lca_rf = 5_000.0 * kg_ref_rf
    lca_t = 5_000.0 * kg_ref_t

    out_kg = calculate_co2_pulse_equivalents(
        lca_rf,
        lca_t,
        kg_ref_rf,
        kg_ref_t,
        years,
        reference_pulse_mass=1.0,
    )
    out_t = calculate_co2_pulse_equivalents(
        lca_rf,
        lca_t,
        1_000.0 * kg_ref_rf,
        1_000.0 * kg_ref_t,
        years,
        reference_pulse_mass=1_000.0,
    )
    out_mt = calculate_co2_pulse_equivalents(
        lca_rf,
        lca_t,
        1.0e9 * kg_ref_rf,
        1.0e9 * kg_ref_t,
        years,
        reference_pulse_mass=1.0e9,
    )

    kg_value = out_kg["co2_pulse_equivalent_integrated_rf"]
    tonne_reference_value = out_t["co2_pulse_equivalent_integrated_rf"]
    mt_reference_value = out_mt["co2_pulse_equivalent_integrated_rf"]

    assert kg_value == pytest.approx(5_000.0)
    assert tonne_reference_value == pytest.approx(kg_value)
    assert mt_reference_value == pytest.approx(kg_value)


def test_make_reference_co2_pulse_emissions_adds_pulse_in_native_unit() -> None:
    df = pd.DataFrame(
        {
            "scenario": ["s", "s"],
            "region": ["World", "World"],
            "variable": ["CO2 FFI", "CH4"],
            "unit": ["Gt CO2/yr", "Mt CH4/yr"],
            "2050.5": [10.0, 1.0],
            "2051.5": [10.0, 1.0],
        }
    )

    out = make_reference_co2_pulse_emissions(
        df,
        scenario="s",
        pulse_year=2050,
        pulse_mass_kg=1.0e9,
    )

    co2_row = out[out["variable"] == "CO2 FFI"].iloc[0]
    ch4_row = out[out["variable"] == "CH4"].iloc[0]
    assert co2_row["2050.5"] == pytest.approx(10.001)
    assert co2_row["2051.5"] == pytest.approx(10.0)
    assert ch4_row["2050.5"] == pytest.approx(1.0)


class DummyPulseTrails:
    def __init__(self) -> None:
        coords = np.array([[0], [0], [0], [0]], dtype=int)
        data = np.array([1.0e9], dtype=float)
        inv = sparse.COO(coords=coords, data=data, shape=(1, 1, 1, 1))
        self.inventory = xr.DataArray(
            inv,
            dims=("activity", "flow", "year", "root activity"),
            coords={
                "activity": [0],
                "flow": [0],
                "year": [2050],
                "root activity": [0],
            },
        )
        self.biosphere_indices = {
            "2050": {
                0: {
                    "name": "Carbon dioxide, fossil",
                    "compartment": "air",
                    "subcompartment": "",
                }
            }
        }
        self.debug = False


def test_run_fair_co2_pulse_equivalents_summarizes_by_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_load_emissions_csv(*args: Any, **kwargs: Any) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "scenario": ["s"],
                "region": ["World"],
                "variable": ["CO2 FFI"],
                "unit": ["Gt CO2/yr"],
                "2050.5": [0.0],
                "2051.5": [0.0],
            }
        )

    def fake_load_species_mapping(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[object, str], dict[object, float]]:
        return {("Carbon dioxide, fossil", "air", ""): "CO2 FFI"}, {}

    def fake_run_fair_emissions(*args: Any, **kwargs: Any) -> Any:
        emissions_df = args[0]
        scenario = args[1]
        config_names = kwargs.get("config_names") or ["c1", "c2"]
        local = emissions_df[
            (emissions_df["scenario"] == scenario)
            & (emissions_df["region"].str.lower() == "world")
        ]
        pulse = float(local.loc[local["variable"] == "CO2 FFI", "2050.5"].iloc[0])
        config_multiplier = np.array([1.0, 2.0], dtype=float)
        years = np.array([2050.5, 2051.5], dtype=float)
        rf_vals = np.zeros((1, len(config_names), len(years), 1), dtype=float)
        rf_vals[0, :, :, 0] = pulse * config_multiplier[:, None]
        temp_vals = np.zeros((1, len(config_names), len(years)), dtype=float)
        temp_vals[0, :, :] = 2.0 * pulse * config_multiplier[:, None]

        return types.SimpleNamespace(
            forcing=xr.DataArray(
                rf_vals,
                dims=("scenario", "config", "timebounds", "specie"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": years,
                    "specie": ["CO2"],
                },
            ),
            temperature=xr.DataArray(
                temp_vals,
                dims=("scenario", "config", "timebounds"),
                coords={
                    "scenario": [scenario],
                    "config": config_names,
                    "timebounds": years,
                },
            ),
            configs=config_names,
        )

    monkeypatch.setattr("trails.fair_rf.load_emissions_csv", fake_load_emissions_csv)
    monkeypatch.setattr(
        "trails.fair_rf.load_species_mapping", fake_load_species_mapping
    )
    monkeypatch.setattr("trails.fair_rf._run_fair_emissions", fake_run_fair_emissions)

    trails = DummyPulseTrails()
    result = run_fair_co2_pulse_equivalents(
        trails,
        scenario="s",
        config_names=["c1", "c2"],
        reference_pulse_year=2050,
        window_start=2050,
        window_end=2052,
        reference_pulse_mass_kg=1.0e9,
        scale_factor=1.0,
    )

    indicator = result["co2_pulse_equivalent"]
    assert indicator["integrated_rf"]["by_config"] == pytest.approx([1.0e9, 1.0e9])
    assert indicator["integrated_temperature"]["by_config"] == pytest.approx(
        [1.0e9, 1.0e9]
    )
    assert indicator["integrated_rf"]["median"] == pytest.approx(1.0e9)
    assert trails.co2_pulse_equivalent is indicator
