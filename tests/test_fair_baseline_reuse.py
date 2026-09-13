"""Compare optimized perturbations with independent, native FaIR runs."""

import copy
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import sparse
import xarray as xr
from trails import fair_rf as rf

SCENARIO = "REMIND|SSP2-PkBudg1000"


@pytest.mark.parametrize(
    "ghg_method", ["myhre1998", "etminan2016", "meinshausen2020", "leach2021"]
)
@pytest.mark.parametrize("count", [1, 3])
@pytest.mark.parametrize(
    "species", ["CO2 FFI", "CO2 AFOLU", "CH4", "N2O", "CF4", "Sulfur", "NOx", "VOC"]
)
def test_native_reuse_matches_full_history(count, species, ghg_method):
    frame = rf.load_emissions_csv(rf.DEFAULT_EMISSIONS_CSV)
    frame = rf._ensure_response_species_rows(
        frame, scenario=SCENARIO, drivers=["NOx", "Sulfur", "BC", "OC"]
    )
    columns, years = rf._extract_year_columns(frame)
    keep = [c for c, y in zip(columns, years) if 2000 <= y <= 2020]
    frame = frame[["scenario", "region", "variable", "unit", *keep]]
    configs = pd.read_csv(rf.DEFAULT_CONFIGS_CSV, index_col=0).index[:count].tolist()
    kw = dict(config_names=configs, ghg_method=ghg_method)
    baseline = rf._run_fair_emissions(frame, SCENARIO, **kw)
    reference_forcing = baseline.forcing.values.copy()
    reference_temperature = baseline.temperature.values.copy()
    rf._prepare_fair_checkpoint(baseline, 2008)
    perturbed = frame.copy()
    rows = perturbed.variable == species
    perturbed.loc[rows, keep[8]] += 0.2
    perturbed.loc[rows, keep[10]] -= 0.1
    expected = rf._run_fair_emissions(perturbed, SCENARIO, **kw)
    actual = rf._run_fair_emissions(
        perturbed,
        SCENARIO,
        _baseline=baseline,
        _changed_species=(species,),
        _perturbation_start_year=2008,
        **kw,
    )
    assert actual._n_timepoints < expected._n_timepoints
    for name in ["forcing", "temperature"]:
        a = getattr(actual, name)
        b = getattr(expected, name).sel(timebounds=a.timebounds)
        if name == "forcing":
            np.testing.assert_allclose(
                a.sum("specie").values, b.sum("specie").values, rtol=1e-11, atol=2e-14
            )
            assert actual._n_species < expected._n_species
            retained = [
                n for n in actual.species if n != actual._trails_background_species
            ]
            a = a.sel(specie=retained)
            b = b.sel(specie=retained)
        np.testing.assert_allclose(
            a.values, b.values, rtol=1e-11, atol=2e-14, equal_nan=True
        )
    np.testing.assert_array_equal(baseline.forcing.values, reference_forcing)
    np.testing.assert_array_equal(baseline.temperature.values, reference_temperature)


@pytest.mark.parametrize("count", [1, 3])
def test_public_signed_quantiles_match_independent_runs(monkeypatch, count):
    inv = sparse.COO(
        np.array([[0, 0, 0, 0], [0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 0, 1]]),
        [1e9, -2e8, 3e9, -4e8],
        shape=(1, 2, 3, 2),
    )
    model = SimpleNamespace(
        debug=False,
        inventory=xr.DataArray(
            inv,
            dims=("activity", "flow", "year", "root activity"),
            coords={
                "activity": [0],
                "flow": [0, 1],
                "year": [2035, 2036, 2037],
                "root activity": [0, 1],
            },
        ),
        biosphere_indices={
            "2035": {
                0: dict(
                    name="Carbon dioxide, fossil", compartment="air", subcompartment=""
                ),
                1: dict(
                    name="Carbon dioxide, in air",
                    compartment="natural resource",
                    subcompartment="in air",
                ),
            }
        },
    )
    reference = copy.deepcopy(model)
    configs = pd.read_csv(rf.DEFAULT_CONFIGS_CSV, index_col=0).index[:count].tolist()
    kw = dict(scenario=SCENARIO, config_names=configs, show_progress=False)
    with monkeypatch.context() as patch:
        patch.setattr(
            rf, "_compact_perturbation", lambda t, b, *a, **kw: (copy.deepcopy(t), b)
        )
        patch.setattr(rf, "_prepare_fair_checkpoint", lambda *a, **kw: None)
        patch.setattr(rf, "_reuse_independent_baseline", lambda f, *a, **kw: f)
        rf.run_fair_delta_rf(reference, per_species_workers=1, **kw)
    rf.run_fair_delta_rf(model, **kw)
    for name in ["instant_radiative_forcing", "delta_temperature"]:
        a, b = getattr(reference, name), getattr(model, name)
        assert a.dims == b.dims
        for dim in a.dims:
            np.testing.assert_array_equal(a.coords[dim], b.coords[dim])
        np.testing.assert_allclose(
            a.data.todense(), b.data.todense(), rtol=1e-7, atol=2e-14
        )


def test_vectorized_calibration_matches_all_native_ensemble_parameters():
    species, properties = rf.read_properties(filename=str(rf.DEFAULT_PROPERTIES_CSV))
    configs = pd.read_csv(rf.DEFAULT_CONFIGS_CSV, index_col=0).index.tolist()
    native = rf.fair.FAIR()
    native.define_time(2000, 2001, 1)
    native.define_scenarios([SCENARIO])
    native.define_configs(configs)
    native.define_species(species, properties)
    native.allocate()
    native.fill_species_configs(filename=str(rf.DEFAULT_PROPERTIES_CSV))
    optimized = copy.deepcopy(native)
    native.override_defaults(str(rf.DEFAULT_CONFIGS_CSV))
    rf._override_fair_defaults(optimized, str(rf.DEFAULT_CONFIGS_CSV))
    xr.testing.assert_identical(optimized.species_configs, native.species_configs)
    xr.testing.assert_identical(optimized.climate_configs, native.climate_configs)


def test_unvalidated_fair_version_uses_full_model(monkeypatch):
    frame = rf.load_emissions_csv(rf.DEFAULT_EMISSIONS_CSV)
    columns, years = rf._extract_year_columns(frame)
    keep = [c for c, y in zip(columns, years) if 2000 <= y <= 2010]
    frame = frame[["scenario", "region", "variable", "unit", *keep]]
    baseline = rf._run_fair_emissions(frame, SCENARIO)
    monkeypatch.setattr(rf.fair, "__version__", "2.3.0")
    rf._prepare_fair_checkpoint(baseline, 2005)
    assert not hasattr(baseline, "_trails_perturbation_template")
    result, reference = rf._compact_perturbation(
        baseline._trails_template, baseline, ("CF4",)
    )
    assert result is not baseline._trails_template
    assert result.species == baseline.species
    assert reference is baseline


def test_early_perturbation_preserves_full_history():
    frame = rf.load_emissions_csv(rf.DEFAULT_EMISSIONS_CSV)
    columns, years = rf._extract_year_columns(frame)
    keep = [c for c, y in zip(columns, years) if 2000 <= y <= 2020]
    frame = frame[["scenario", "region", "variable", "unit", *keep]]
    base = rf._run_fair_emissions(frame, SCENARIO)
    rf._prepare_fair_checkpoint(base, 2010)
    frame.loc[frame.variable == "CO2 AFOLU", keep[2]] -= 0.1
    expected = rf._run_fair_emissions(frame, SCENARIO)
    actual = rf._run_fair_emissions(
        frame,
        SCENARIO,
        _baseline=base,
        _changed_species=("CO2 AFOLU",),
        _perturbation_start_year=2002,
    )
    np.testing.assert_array_equal(actual.timebounds, expected.timebounds)
    np.testing.assert_allclose(
        actual.temperature, expected.temperature, rtol=1e-11, atol=2e-14
    )
    np.testing.assert_allclose(
        actual.forcing.sum("specie"),
        expected.forcing.sum("specie"),
        rtol=1e-11,
        atol=2e-14,
    )


def test_checkpoint_preserves_seeded_stochastic_temperature(tmp_path):
    frame = rf.load_emissions_csv(rf.DEFAULT_EMISSIONS_CSV)
    columns, years = rf._extract_year_columns(frame)
    keep = [c for c, y in zip(columns, years) if 2000 <= y <= 2020]
    frame = frame[["scenario", "region", "variable", "unit", *keep]]
    calibration = pd.read_csv(rf.DEFAULT_CONFIGS_CSV, index_col=0).iloc[:2].copy()
    calibration["stochastic_run"] = True
    calibration["use_seed"] = True
    calibration["seed"] = [42, 43]
    path = tmp_path / "seeded.csv"
    calibration.to_csv(path)
    kw = dict(config_csv=path, config_names=calibration.index.tolist())
    base = rf._run_fair_emissions(frame, SCENARIO, **kw)
    rf._prepare_fair_checkpoint(base, 2010)
    frame.loc[frame.variable == "CO2 FFI", keep[12]] += 0.1
    expected = rf._run_fair_emissions(frame, SCENARIO, **kw)
    actual = rf._run_fair_emissions(
        frame,
        SCENARIO,
        _baseline=base,
        _changed_species=("CO2 FFI",),
        _perturbation_start_year=2012,
        **kw,
    )
    expected_temperature = expected.temperature.sel(timebounds=actual.timebounds)
    np.testing.assert_allclose(
        actual.temperature, expected_temperature, rtol=1e-11, atol=2e-14
    )
    np.testing.assert_allclose(
        actual.stochastic_forcing,
        expected.stochastic_forcing.sel(timebounds=actual.timebounds),
        rtol=1e-13,
        atol=2e-14,
    )
