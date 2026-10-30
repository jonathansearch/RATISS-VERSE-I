"""
Linear growth of structure.

The single equation solved here is the standard one:

    D'' + (2 + dlnH/dN) D' - (3/2) Omega_m(N) mu(k,a) D = 0

with N = ln a, ' = d/dN, and mu(k,a) the ratio of the effective
gravitational coupling to Newton's constant.

``mu = 1`` is General Relativity. Every deviation from GR that RATISS
supports is expressed as a function ``mu(k, a)``, which is why the same
solver covers Lambda-CDM, the granular universe and every modified-gravity
variant.

The growth index is

    gamma = ln f / ln Omega_m(a),      f = dlnD/dlna
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from . import constants as C


OM0 = C.OMEGA_M
N_INIT = -9.0          # a_init = e^-9 ~ 1.2e-4 : deep in the matter era


def growth_factor(
    E: Callable[[float], float],
    mu: Callable[[float, float], float] | float = 1.0,
    k: float = 0.0,
    N_fin: float = 0.0,
    om0: float = OM0,
) -> tuple[Callable, float, float, float]:
    """Solve the growth equation.

    Parameters
    ----------
    E : callable
        H(a)/H0.
    mu : callable or float
        Effective gravitational coupling, ``mu(k, a)``.
    k : float
        Wavenumber in h/Mpc passed to ``mu``.
    N_fin : float
        Final e-fold (0 = today).

    Returns
    -------
    solver, D(N_fin), D'(N_fin), Omega_m(N_fin)
    """
    mu_c = mu if callable(mu) else (lambda kk, aa: float(mu))

    def H2(N):
        return max(E(np.exp(N)) ** 2, 1e-30)

    def rhs(N, y):
        D, dD = y
        a = np.exp(N)
        h2 = H2(N)
        eps = 1e-4
        dlnH = 0.5 * (np.log(H2(N + eps)) - np.log(H2(N - eps))) / (2 * eps)
        omN = om0 * a ** -3 / h2
        m = mu_c(k, a)
        return [dD, -((2.0 + dlnH) * dD - 1.5 * omN * m * D)]

    a_i = np.exp(N_INIT)
    # matter-era growing mode: D ~ a, D' ~ a
    y0 = [a_i, a_i]
    sol = solve_ivp(rhs, (N_INIT, N_fin), y0, method="LSODA",
                    rtol=1e-11, atol=1e-15, dense_output=True,
                    max_step=0.05)
    D, dD = sol.sol(N_fin)
    om_fin = om0 * np.exp(-3 * N_fin) / H2(N_fin)
    return sol, float(D), float(dD), float(om_fin)


def f_sigma8(
    E: Callable[[float], float],
    mu: Callable[[float, float], float] | float = 1.0,
    k: float = 0.0,
    redshifts=(0.0, 0.5, 1.0, 1.5, 2.0),
    om0: float = OM0,
) -> list[dict]:
    """f sigma8(z) for a given universe.

    Sigma8 is normalised to the Planck value today, so what is predicted is
    the *shape* of f sigma8(z), which is what DESI and Euclid measure.
    """
    out = []
    _, D0, dD0, _ = growth_factor(E, mu, k, N_fin=0.0, om0=om0)
    for z in redshifts:
        N = -np.log1p(z)
        _, Dz, dDz, _ = growth_factor(E, mu, k, N_fin=N, om0=om0)
        Dz_rel = Dz / D0
        fz = (dDz / Dz) if Dz != 0 else 0.0
        out.append({
            "z": float(z),
            "f": float(fz),
            "D_rel": float(Dz_rel),
            "f_sigma8": float(fz * Dz_rel * C.SIGMA_8),
        })
    return out


def growth_index(
    E: Callable[[float], float],
    mu: Callable[[float, float], float] | float = 1.0,
    k: float = 0.0,
) -> tuple[float, float]:
    """(f(z=0), gamma) with gamma = ln f / ln Omega_m."""
    _, D, dD, om = growth_factor(E, mu, k, N_fin=0.0)
    if D == 0 or dD / D <= 0 or om <= 0 or om >= 1:
        return float("nan"), float("nan")
    f = dD / D
    return float(f), float(np.log(f) / np.log(om))


def growth_rate_curve(
    E: Callable[[float], float],
    mu: Callable[[float, float], float] | float = 1.0,
    k: float = 0.0,
    n: int = 40,
) -> dict:
    """D(a)/D(1) and f(a) on a log-spaced grid — used by the figures."""
    a = np.linspace(0.02, 1.0, n)
    Ds, fs = [], []
    _, D0, _, _ = growth_factor(E, mu, k, N_fin=0.0)
    for av in a:
        _, Dv, dDv, _ = growth_factor(E, mu, k, N_fin=np.log(av))
        Ds.append(Dv / D0)
        fs.append(dDv / Dv if Dv else np.nan)
    return {"a": a.tolist(), "D": Ds, "f": fs}
