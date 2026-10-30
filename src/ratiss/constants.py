"""
Fiducial cosmological parameters.

All values are frozen here, once, so that every number produced by this
repository is traceable to a single source.

Source: Planck 2018 TT,TE,EE+lowE+lensing (Planck Collaboration 2020,
A&A 641, A6), flat Lambda-CDM column. The values below are the ones used
for every figure and every table in ``results/``.
"""

from __future__ import annotations

# --- Planck 2018 (Table 2, last column: TT,TE,EE+lowE+lensing) -------------
H0_KM_S_MPC = 67.36
H0_OVER_100 = H0_KM_S_MPC / 100.0          # h = 0.6736
OMEGA_B = 0.02237                           # physical baryon density
OMEGA_C = 0.1200                            # physical cold dark matter density
OMEGA_M = 0.3153                            # total matter, Omega_b + Omega_c
OMEGA_LAMBDA = 1.0 - OMEGA_M
SIGMA_8 = 0.8111
N_S = 0.9649
T_CMB = 2.7255                              # K
OMEGA_R = 9.15e-5                           # photons + 3 massless neutrinos

# --- derived ---------------------------------------------------------------
H0_SI = H0_KM_S_MPC * 1.0e3 / 3.0856775814913673e22   # s^-1
C_KM_S = 299792.458
MPC_M = 3.0856775814913673e22
AU_M = 1.495978707e11
PC_M = 3.0856775814913673e16
YR_S = 3.15576e7

# --- derived times ---------------------------------------------------------
AGE_UNIVERSE_GYR = (2.0 / (3.0 * H0_SI * OMEGA_M ** 0.5)) * 0.9644 / 1.0e9
# (approximate Lambda-CDM age; the exact value is computed in background.py)


def hubble_time_gyr() -> float:
    """1/H0 in Gyr."""
    return 1.0 / H0_SI / YR_S / 1e9
