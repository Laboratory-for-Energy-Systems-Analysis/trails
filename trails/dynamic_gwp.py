import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Union, Literal


# -----------------------------
# Parameters
# -----------------------------
@dataclass(frozen=True)
class WMGHGParams:
    """Well-mixed GHG treated as a linear, first-order box model (single lifetime)."""

    lifetime_yr: float  # tau
    rad_eff_Wm2_per_ppb: float  # epsilon (W m-2 ppb-1)


@dataclass(frozen=True)
class CO2Params:
    """
    CO2 handled with:
      - multi-exponential impulse response for concentration (airborne fraction)
      - logarithmic forcing formula
    """

    # Airborne fraction IRF: a0 + sum(ai * exp(-t/taui))
    a0: float
    a: Tuple[float, ...]
    tau: Tuple[float, ...]  # years

    # Forcing formula selection + baseline
    forcing_formula: Literal["myhre1998", "etminan2016", "meinshausen2020"] = (
        "meinshausen2020"
    )
    C0_ppm: float = 278.3  # baseline CO2 concentration (ppm)


GasParams = Union[WMGHGParams, CO2Params]


# ---- Constants for kg <-> mixing ratio conversion (global mean) ----
M_ATM_KG = 5.1480e18
M_AIR_KG_PER_MOL = 28.97e-3
N_AIR_MOL = M_ATM_KG / M_AIR_KG_PER_MOL

# Mean molar masses (g/mol)
MOLAR_MASS_AIR = 28.97  # g/mol


# -----------------------------
# Unit conversions
# -----------------------------
def kg_to_ppb(mass_kg: np.ndarray, molar_mass_g_per_mol: float) -> np.ndarray:
    """kg -> ppb (well-mixed), assuming instantaneous uniform mixing."""
    molar_mass_kg_per_mol = molar_mass_g_per_mol / 1000.0
    kg_per_ppb = N_AIR_MOL * 1e-9 * molar_mass_kg_per_mol
    return mass_kg / kg_per_ppb


def kg_to_ppm_co2(mass_kg: np.ndarray) -> np.ndarray:
    """kg CO2 -> ppm CO2 (well-mixed), assuming instantaneous uniform mixing."""
    molar_mass_co2_kg_per_mol = 44.0095 / 1000.0  # g/mol -> kg/mol
    kg_per_ppm = N_AIR_MOL * 1e-6 * molar_mass_co2_kg_per_mol
    return mass_kg / kg_per_ppm


# -----------------------------
# RF for WMGHGs (your current method)
# -----------------------------
def rf_wmghg_from_annual_emissions(
    years: np.ndarray,
    emissions_kg_per_yr: np.ndarray,
    params: WMGHGParams,
    molar_mass_g_per_mol: float,
    *,
    pulse_placement: str = "midyear",
) -> np.ndarray:
    years = np.asarray(years, dtype=int)
    E = np.asarray(emissions_kg_per_yr, dtype=float)
    if E.ndim == 1:
        if E.shape != years.shape:
            raise ValueError("years and emissions_kg_per_yr must have the same shape")
    elif E.ndim == 2:
        if E.shape[0] != years.shape[0]:
            raise ValueError("emissions_kg_per_yr first dimension must match years")
    else:
        raise ValueError("emissions_kg_per_yr must be 1D or 2D")

    dX_ppb = kg_to_ppb(E, molar_mass_g_per_mol)

    t = years[:, None].astype(float)
    y = years[None, :].astype(float)

    if pulse_placement == "startyear":
        dt = t - y
    elif pulse_placement == "midyear":
        dt = t - y + 0.5
    else:
        raise ValueError("pulse_placement must be 'startyear' or 'midyear'")

    mask = dt >= 0.0
    dt = np.maximum(dt, 0.0)
    kernel = np.exp(-dt / params.lifetime_yr) * mask
    X_ppb = kernel @ dX_ppb
    return params.rad_eff_Wm2_per_ppb * X_ppb


# -----------------------------
# RF for CO2 (new)
# -----------------------------
def co2_airborne_fraction(dt: np.ndarray, p: CO2Params) -> np.ndarray:
    """
    dt: array of years since pulse (>=0)
    returns fraction remaining in atmosphere.
    """
    af = np.full_like(dt, p.a0, dtype=float)
    for ai, taui in zip(p.a, p.tau):
        af += ai * np.exp(-dt / taui)
    return af


def rf_co2_from_annual_emissions(
    years: np.ndarray,
    emissions_kg_per_yr: np.ndarray,
    params: CO2Params,
    *,
    pulse_placement: str = "midyear",
) -> np.ndarray:
    """
    CO2:
      - convert each annual emission pulse to an instantaneous ppm perturbation
      - propagate with multi-exponential airborne fraction
      - compute forcing using a log relationship:
          ΔF = 5.35 * ln(C/C0)   (classic Myhre form)
        (You can swap to another fit later if you want.)
    """
    years = np.asarray(years, dtype=int)
    E = np.asarray(emissions_kg_per_yr, dtype=float)
    if E.ndim == 1:
        if E.shape != years.shape:
            raise ValueError("years and emissions_kg_per_yr must have the same shape")
    elif E.ndim == 2:
        if E.shape[0] != years.shape[0]:
            raise ValueError("emissions_kg_per_yr first dimension must match years")
    else:
        raise ValueError("emissions_kg_per_yr must be 1D or 2D")

    # ppm pulse amplitude per year
    dC_ppm = kg_to_ppm_co2(E)

    t = years[:, None].astype(float)
    y = years[None, :].astype(float)

    if pulse_placement == "startyear":
        dt = t - y
    elif pulse_placement == "midyear":
        dt = t - y + 0.5
    else:
        raise ValueError("pulse_placement must be 'startyear' or 'midyear'")

    mask = dt >= 0.0
    dt = np.maximum(dt, 0.0)

    # Concentration perturbation from pulse train
    af = co2_airborne_fraction(dt, params)  # (T,T)
    af = af * mask
    dC_t_ppm = af @ dC_ppm  # (T,) or (T,N)

    # Forcing (simple, robust default)
    # Treat perturbation as added on top of baseline C0
    C0 = params.C0_ppm
    if params.forcing_formula in ("myhre1998", "etminan2016", "meinshausen2020"):
        # Use log1p to preserve tiny perturbations (avoids precision loss when dC << C0)
        frac = dC_t_ppm / C0
        frac = np.clip(frac, -0.999999999, None)
        return 5.35 * np.log1p(frac)
    else:
        raise ValueError(f"Unknown CO2 forcing formula: {params.forcing_formula}")


# -----------------------------
# Suite driver (modified)
# -----------------------------
def rf_suite(
    years: np.ndarray,
    emissions_by_gas: Dict[str, np.ndarray],
    gas_params: Dict[str, GasParams],
    molar_masses_g_per_mol: Dict[str, float],
    *,
    pulse_placement: str = "midyear",
) -> Tuple[Dict[str, np.ndarray], np.ndarray]:

    rf_by_gas: Dict[str, np.ndarray] = {}
    total = np.zeros_like(np.asarray(years, dtype=float), dtype=float)

    for gas, E in emissions_by_gas.items():
        if gas not in gas_params:
            raise KeyError(f"Missing gas_params for '{gas}'")

        p = gas_params[gas]

        if isinstance(p, CO2Params):
            RF = rf_co2_from_annual_emissions(
                years=years,
                emissions_kg_per_yr=E,
                params=p,
                pulse_placement=pulse_placement,
            )
        elif isinstance(p, WMGHGParams):
            if gas not in molar_masses_g_per_mol:
                raise KeyError(f"Missing molar mass for '{gas}'")
            RF = rf_wmghg_from_annual_emissions(
                years=years,
                emissions_kg_per_yr=E,
                params=p,
                molar_mass_g_per_mol=molar_masses_g_per_mol[gas],
                pulse_placement=pulse_placement,
            )
        else:
            raise TypeError(f"Unsupported params type for gas '{gas}': {type(p)}")

        rf_by_gas[gas] = RF
        total += RF

    return rf_by_gas, total
