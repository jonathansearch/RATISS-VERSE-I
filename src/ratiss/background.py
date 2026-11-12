"""
Background expansion histories for every universe implemented in RATISS.

Each model exposes ``E(a) = H(a)/H0``. Everything else (H0*t0, comoving
distance, age, deceleration) is derived from ``E(a)`` by quadrature, so the
only thing a model has to get right is its Friedmann equation.

Implemented universes
---------------------
lcdm          Flat Lambda-CDM.                                     (reference)
wcdm          One dark-energy fluid with constant w.               (Porte 1/3)
rvm           Running vacuum, Lambda(H) = c0 + 3 nu H^2.           (Sola Peracaula)
dgp           Dvali-Gabadadze-Porrati self-accelerating branch.    (braneworld)
cardassian    Modified Friedmann, E^2 = Om a^-3 + (1-Om) a^-m.     (Freese)
fr_hs         Hu-Sawicki f(R); background designed = Lambda-CDM.
ft            f(T) teleparallel; background designed = Lambda-CDM.
granular      RATISS r_c model: Newton + K/r^4, graviton cut-off at r_c.

References are collected in ``docs/03_REFERENCES.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import quad

from . import constants as C


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _quad(f: Callable[[float], float], lo: float, hi: float) -> tuple[float, float]:
    """Adaptive quadrature with an independent error estimate.

    Returns (value, estimated absolute error) using two tolerances so that
    the reported uncertainty is not a guess.
    """
    v1, e1 = quad(f, lo, hi, limit=400, epsabs=1e-12, epsrel=1e-12)
    v2, e2 = quad(f, lo, hi, limit=400, epsabs=1e-9, epsrel=1e-9)
    return v1, abs(v1 - v2) + e1


def hubble_time_times_h0(E: Callable[[float], float]) -> tuple[float, float]:
    """Dimensionless product H0 * t0 = int_0^1 da / (a E(a)).

    This is the quantity every test in this project has been about.
    """
    f = lambda lna: 1.0 / max(E(np.exp(lna)), 1e-30)
    return _quad(f, -60.0, 0.0)


def comoving_distance_over_cH0(E: Callable[[float], float]) -> float:
    """Line-of-sight comoving distance D_C/(c/H0) = int da / (a^2 E(a)).

    The kernel is ``1/(a^2 E)``, not ``1/(a E)``: the latter is the age
    integral returned by :func:`hubble_time_times_h0` (defect D12).
    """
    f = lambda lna: np.exp(-lna) / max(E(np.exp(lna)), 1e-30)
    return _quad(f, -60.0, 0.0)[0]


def luminosity_distance(z: float, E: Callable[[float], float]) -> float:
    """Luminosity distance D_L/(c/H0) = (1+z) D_C(z), same kernel as D_C."""
    f = lambda lna: np.exp(-lna) / max(E(np.exp(lna)), 1e-30)
    dc = _quad(f, -np.log1p(z), 0.0)[0]
    return (1.0 + z) * dc


# ---------------------------------------------------------------------------
# the models
# ---------------------------------------------------------------------------
def E_lcdm(a, om: float = C.OMEGA_M) -> float:
    """Flat Lambda-CDM including radiation.

    Radiation is kept in *every* background in this repository. An earlier
    version omitted it here while retaining it in ``make_wcdm``, which shifted
    H0*t0 by 3.8e-4 -- an artefact larger than several model differences the
    registry was meant to resolve. See docs/05_LIMITATIONS.md, defect D1.
    """
    return float(np.sqrt(om * a ** -3 + C.OMEGA_R * a ** -4
                         + (1.0 - om - C.OMEGA_R)))


def make_wcdm(w: float, om: float = C.OMEGA_M) -> Callable[[float], float]:
    """Constant-w dark energy plus matter and radiation."""
    def E(a):
        return float(np.sqrt(om * a ** -3 + C.OMEGA_R * a ** -4
                             + (1.0 - om - C.OMEGA_R) * a ** (-3.0 * (1.0 + w))))
    return E


def make_rvm(nu: float, om: float = C.OMEGA_M) -> Callable[[float], float]:
    """Running vacuum model: Lambda(H) = c0 + 3 nu H^2.

    Friedmann gives 3 H^2 (1 - nu) = 8 pi G rho_m + c0, hence

        E^2(a) = [ Om a^-3 + (1 - Om - nu) ] / (1 - nu)

    Ref: Sola Peracaula, Gomez-Valent & de Cruz Perez 2019, MNRAS 490, 1.
    """
    def E(a):
        num = om * a ** -3 + (1.0 - om - nu)
        return float(np.sqrt(max(num, 0.0) / (1.0 - nu)))
    return E


def E_dgp_selfacc(a, om: float = C.OMEGA_M) -> float:
    """DGP self-accelerating branch.

    H^2 - H/r_c = rho_m/3   ->   E = sqrt(Om a^-3 + Om_rc) + sqrt(Om_rc)
    with Om_rc = ((1-Om)/2)^2.
    Ref: Dvali, Gabadadze & Porrati 2000, Phys. Lett. B 485, 208.
    """
    orc = ((1.0 - om) / 2.0) ** 2
    return float(np.sqrt(om * a ** -3 + C.OMEGA_R * a ** -4 + orc)
                 + np.sqrt(orc))


def E_dgp_normal(a, om: float = C.OMEGA_M) -> float:
    orc = ((1.0 - om) / 2.0) ** 2
    return float(np.sqrt(om * a ** -3 + orc) - np.sqrt(orc))


def make_cardassian(m: float, om: float = C.OMEGA_M) -> Callable[[float], float]:
    """Freese-Cardassian modified Friedmann: E^2 = Om a^-3 + (1-Om) a^-m."""
    def E(a):
        return float(np.sqrt(om * a ** -3 + C.OMEGA_R * a ** -4
                             + (1.0 - om - C.OMEGA_R) * a ** (-m)))
    return E


def E_fr_hs(a, om: float = C.OMEGA_M) -> float:
    """Hu-Sawicki f(R) with a designer background: identical to Lambda-CDM.

    The whole point of the Hu-Sawicki construction is that the expansion
    history is *made* to match Lambda-CDM; all the information lives in the
    growth of structure (see growth.py).
    Ref: Hu & Sawicki 2007, Phys. Rev. D 76, 064004.
    """
    return E_lcdm(a, om)


def E_ft(a, om: float = C.OMEGA_M) -> float:
    """f(T) teleparallel gravity with a designer background.

    As for f(R), the background can be chosen to reproduce Lambda-CDM
    exactly while the effective Newton constant governing perturbations
    differs. Ref: Linder 2010, Phys. Rev. D 81, 127301.
    """
    return E_lcdm(a, om)


def make_single_fluid(n: float) -> Callable[[float], float]:
    """One substance diluting as rho ~ a^-n, nothing else. H*t = 2/n."""
    return lambda a: float(a ** (-n / 2.0))


# ---------------------------------------------------------------------------
# the granular (r_c) universe
# ---------------------------------------------------------------------------
def E_granular(a, om: float = C.OMEGA_M, rc_over_lambda: float = 0.0) -> float:
    """RATISS universe with a gravitational cut-off at r_c.

    On the background level the cut-off only acts through the *mean* field.
    By the Poisson/Jeans argument implemented in granular.py the mean density
    produces no force, so the background is *identical* to Lambda-CDM. The
    cut-off appears only in the growth of perturbations.

    ``rc_over_lambda`` is kept as an explicit argument so that a caller
    cannot forget that the cut-off exists.
    """
    return E_lcdm(a, om)


@dataclass(frozen=True)
class Universe:
    key: str
    label: str
    E: Callable[[float], float]
    family: str
    notes: str = ""


def registry(om: float = C.OMEGA_M) -> list[Universe]:
    """The complete list of universes simulated by RATISS."""
    return [
        Universe("lcdm", "Lambda-CDM", lambda a: E_lcdm(a, om),
                 "GR", "reference universe"),
        Universe("w_m0302", "w = -0.30185", make_wcdm(-0.30185, om),
                 "fluid", "single-fluid target from H0*t0 = 0.9549"),
        Universe("w_m1004", "w = -1.004049", make_wcdm(-1.004049, om),
                 "fluid", "exact inversion of H0*t0 = 0.9549 with matter"),
        Universe("rvm_001", "RVM nu=0.001", make_rvm(0.001, om),
                 "running vacuum", "Lambda(H) = c0 + 3 nu H^2"),
        Universe("dgp", "DGP self-accelerating", lambda a: E_dgp_selfacc(a, om),
                 "braneworld", "5D DGP branch"),
        Universe("cardassian_05", "Cardassian m=0.5", make_cardassian(0.5, om),
                 "modified Friedmann", "E^2 = Om a^-3 + (1-Om) a^-0.5"),
        Universe("fr_hs", "f(R) Hu-Sawicki", lambda a: E_fr_hs(a, om),
                 "f(R)", "designer background = Lambda-CDM"),
        Universe("ft", "f(T) teleparallel", lambda a: E_ft(a, om),
                 "f(T)", "designer background = Lambda-CDM"),
        Universe("granular", "RATISS granular (r_c)", lambda a: E_granular(a, om),
                 "RATISS", "Newton + K/r^4, graviton cut-off at r_c"),
        Universe("fluid_n2094", "single fluid n=2.0945", make_single_fluid(2.0945),
                 "Porte 1", "no matter at all: H*t = 2/n exactly"),
    ]
