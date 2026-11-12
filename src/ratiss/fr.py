"""
Hu--Sawicki f(R) gravity with chameleon screening.

The v1.0.2 registry carried f(R) and f(T) rows with ``mu = 1`` and
``E = E_LambdaCDM``. Their agreement with Lambda-CDM was therefore an
identity, not a result (see docs/05_LIMITATIONS.md, defect D6, and the
negative result N3). This module supplies the missing physics: a
*scale-dependent* effective coupling ``mu(k, a)`` so that the same growth
solver produces a genuine, falsifiable f(R) prediction.

Quasi-static effective coupling
-------------------------------
On sub-horizon scales the scalaron mediates a fifth force that is screened
by the chameleon mechanism. The standard quasi-static result is

    mu(k, a) = G_eff/G = 1 + (1/3) * x / (1 + x),     x = (k / (a M(a)))^2

with ``M(a)`` the scalaron mass. The two exact limits are the physics:

    x -> 0    mu -> 1     (large scales, unscreened fifth force is absent)
    x -> inf  mu -> 4/3   (small scales, fully unscreened)

Scalaron mass
-------------
For Hu--Sawicki with n = 1 the background field scales as
``f_R(a) ~ -|f_R0| (R0/R(a))^2``, and ``M^2 = 1/(3 f_RR)`` gives

    M(a)^2 = M0^2 * (R(a)/R0)^3,
    M0^2/H0^2 = (Omega_m + 4 Omega_Lambda) / (2 |f_R0|),
    R(a)/R0 = (a^-3 + 4 Omega_Lambda/Omega_m) / (1 + 4 Omega_Lambda/Omega_m).

The single free parameter is ``|f_R0|``. It is scanned over
``1e-4 ... 1e-7`` in ``models.fr_sector``, spanning Compton wavelengths of a
few Mpc/h down to sub-Mpc -- inside the reach of DESI and Euclid.

References: Hu & Sawicki 2007, Phys. Rev. D 76, 064004; Tsujikawa 2008,
Phys. Rev. D 77, 023507; Lombriser 2014, Ann. Phys. 529, 1600338.
"""

from __future__ import annotations

import numpy as np

from . import constants as C


# (100 h km/s/Mpc) / c, i.e. H0 in units of h/Mpc with h factored out.
K_H0 = 100.0 / C.C_KM_S


def _m_over_H0(a, f_R0: float, om: float = C.OMEGA_M) -> np.ndarray:
    """Scalaron mass M(a) in units of H0 (h factors cancel in k/(aM))."""
    ol = 1.0 - om
    m0 = np.sqrt((om + 4.0 * ol) / (2.0 * abs(f_R0)))
    a = np.asarray(a, dtype=float)
    scale = ((a ** -3 + 4.0 * ol / om) / (1.0 + 4.0 * ol / om)) ** 1.5
    return m0 * scale


def compton_wavelength_today(f_R0: float, om: float = C.OMEGA_M) -> float:
    """Scalaron Compton wavelength today, in Mpc/h."""
    return 1.0 / (K_H0 * _m_over_H0(1.0, f_R0, om))


def mu_fR(k_h_mpc, a, f_R0: float, om: float = C.OMEGA_M):
    """Effective coupling mu(k, a) for Hu--Sawicki f(R), k in h/Mpc."""
    a = np.asarray(a, dtype=float)
    M_hmpc = K_H0 * _m_over_H0(a, f_R0, om)
    x = (np.asarray(k_h_mpc, dtype=float) / (a * M_hmpc)) ** 2
    return 1.0 + (1.0 / 3.0) * x / (1.0 + x)


def mu_fR_callable(f_R0: float, om: float = C.OMEGA_M):
    """Callable mu(k, a) in the form expected by ``growth.growth_factor``."""
    def mu(k, a):
        return float(mu_fR(k, a, f_R0, om))
    return mu
