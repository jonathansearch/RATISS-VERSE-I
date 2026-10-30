"""
The RATISS granular universe: SYNTORON as a cosmological model.

The toy model that seeded this repository used the pair force

    F(r) = -G/r^2 + K/r^4                        (SYNTORON 24, S03b)

which factorises exactly as

    F(r) = -(G/r^2) * [ 1 - (r_c/r)^2 ],         r_c = sqrt(K/G)

so that three regimes exist:

    r >> r_c    Newton, unmodified
    r  = r_c    the force vanishes exactly
    r << r_c    the force is repulsive

``r_c`` is the only length in the model, and ``eta = r_c/EPS`` is the only
dimensionless number. By the theorem proved in ``docs/02_THEOREMS.md``
(and in the parent project), ``r_c`` carries no absolute value on its own:
the simulation is invariant under a simultaneous rescaling of G, K and the
softening.

What this module does is therefore *not* to invent a value of r_c. It
translates the force law into a scale-dependent gravitational coupling
``mu(k)`` and propagates it into the growth of structure, which is where a
cut-off leaves an observable trace.

Mapping
-------
The relative force deficit is ``1 - mu = (r_c/r)^2``: a two-derivative
(higher-order) correction. Fourier transforming and regularising at the
cut-off itself gives the bounded, monotone interpolation

    mu(x) = (1 - x^2) / (1 + x^2),      x = k * r_c

which reproduces all three regimes:

    mu(0) = +1   (Newton)
    mu(1) =  0   (cut-off)
    mu(->inf) = -1 (repulsion, bounded)

The exact Fourier transform of a 1/r^4 force is not a tempered
distribution (it diverges logarithmically in the infrared), which is
precisely why the model *needs* a cut-off at r_c and cannot exist without
one. The interpolation above is the minimal regularisation that respects
the exact limits.
"""

from __future__ import annotations

import numpy as np

from . import constants as C


# ---------------------------------------------------------------------------
# force law
# ---------------------------------------------------------------------------
def force(r, G: float = 1.0, K: float = 0.01):
    """F(r) for unit mass. Negative = attractive."""
    r = np.asarray(r, dtype=float)
    return -G / r ** 2 + K / r ** 4


def r_c(G: float = 1.0, K: float = 0.01) -> float:
    """The equilibrium radius where the two terms cancel: r_c = sqrt(K/G)."""
    return float(np.sqrt(K / G))


def eta(G: float = 1.0, K: float = 0.01, eps: float = 0.1) -> float:
    return r_c(G, K) / eps


def force_ratio(r_over_rc):
    """F / F_Newton = 1 - (r_c/r)^2."""
    x = np.asarray(r_over_rc, dtype=float)
    return 1.0 - 1.0 / x ** 2


# ---------------------------------------------------------------------------
# effective coupling
# ---------------------------------------------------------------------------
def mu_effective(x):
    """mu(x) with x = k * r_c. See the module docstring for the derivation."""
    x = np.asarray(x, dtype=float)
    return (1.0 - x ** 2) / (1.0 + x ** 2)


def mu_of_k(k_h_mpc, r_c_mpc: float):
    """Callable mu(k, a) in the form expected by ``growth.growth_factor``."""
    def mu(kk, a):
        return float(mu_effective(np.asarray(kk) * r_c_mpc))
    return mu


# ---------------------------------------------------------------------------
# consequences
# ---------------------------------------------------------------------------
def growth_suppression(x_values, D_of_x, D_gr):
    """D(x)/D_GR for a list of x values."""
    return [d / D_gr for d in D_of_x]


def smallest_surviving_scale(r_c_mpc: float, x_cut: float = 0.1) -> float:
    """Largest comoving wavelength that survives a cut-off r_c.

    Structure on scale L corresponds to k = 2 pi / L. Requiring
    x = k r_c <= x_cut gives L >= 2 pi r_c / x_cut.
    """
    return 2.0 * np.pi * r_c_mpc / x_cut


def rc_bound_from(L_min_mpc: float, x_cut: float = 0.1) -> float:
    """Upper bound on r_c implied by the existence of structure of size L_min.

    Inverting the relation above:  r_c <= x_cut * L_min / (2 pi).
    """
    return x_cut * L_min_mpc / (2.0 * np.pi)


# ---------------------------------------------------------------------------
# reference objects used for the bound
# ---------------------------------------------------------------------------
REFERENCE_STRUCTURES = [
    # (name, comoving size in Mpc, status)
    ("galaxy clusters", 2.0e0, "observed"),
    ("galaxies", 1.0e-2, "observed"),
    ("ultra-faint dwarfs", 1.0e-4, "observed"),
    ("dark matter microhalos", 1.0e-8, "predicted (not yet observed)"),
]


def bounds_table(x_cut: float = 0.1) -> list[dict]:
    rows = []
    for name, L, status in REFERENCE_STRUCTURES:
        rc_max = rc_bound_from(L, x_cut)
        rows.append({
            "structure": name,
            "size_mpc": L,
            "status": status,
            "r_c_max_mpc": rc_max,
            "r_c_max_m": rc_max * C.MPC_M,
            "r_c_max_au": rc_max * C.MPC_M / C.AU_M,
        })
    return rows
