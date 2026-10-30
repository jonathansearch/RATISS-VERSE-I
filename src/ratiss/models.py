"""
The RATISS universe registry and its predictions.

A RATISS *reference universe* is a tuple

    (background expansion, linear growth, effective coupling mu(k,a))

from which every observable quantity is derived by the same two solvers
(background.py and growth.py). No quantity is ever hand-typed into a table:
everything in ``results/`` is recomputed from the registry at build time.

The registry is deliberately redundant. Several universes give the same
``H0*t0`` and are therefore *indistinguishable* by expansion alone. That
degeneracy is itself one of the results this repository reports.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import background as bg
from . import constants as C
from . import granular as gr
from . import growth as gw


OM = C.OMEGA_M
# The observed expansion history is the best-fit Lambda-CDM one. Every model
# is therefore compared with Lambda-CDM computed at the *same* Omega_m, so
# that the comparison isolates the model and not the choice of parameters.
_HT_REF, _ = bg.hubble_time_times_h0(bg.E_lcdm)
HT_MEASURED = _HT_REF

# The growth index of the reference universe is *derived*, never typed.
# An earlier version compared against the literal 0.5544, which is a magic
# number (see docs/05_LIMITATIONS.md, defect D4).
_GAMMA_REF = None


def gamma_reference() -> float:
    """gamma of the Lambda-CDM reference, recomputed on first use."""
    global _GAMMA_REF
    if _GAMMA_REF is None:
        _GAMMA_REF = gw.growth_index(bg.E_lcdm, 1.0)[1]
    return _GAMMA_REF


@dataclass
class Prediction:
    key: str
    label: str
    family: str
    H0t0: float
    H0t0_err: float
    gamma: float
    f_z0: float
    fsigma8_z1: float
    status: str
    note: str


def _status(ht: float, gamma: float) -> tuple[str, str]:
    d_ht = abs(ht - HT_MEASURED)
    if d_ht > 0.02:
        return "EXCLUDED", f"H0*t0 off by {d_ht:.3f}"
    if np.isnan(gamma):
        return "UNDEFINED", "growth index undefined"
    gref = gamma_reference()
    if abs(gamma - gref) > 0.03:
        return "TESTABLE", f"gamma deviates by {gamma - gref:+.3f}"
    return "DEGENERATE", "indistinguishable from Lambda-CDM at background level"


def build_registry() -> list[Prediction]:
    rows: list[Prediction] = []
    for u in bg.registry(OM):
        ht, ht_err = bg.hubble_time_times_h0(u.E)
        if u.key == "granular":
            # the cut-off lives in the growth sector only (see granular.py)
            f, gam = gw.growth_index(u.E, 1.0)
        else:
            f, gam = gw.growth_index(u.E, 1.0)
        fs8 = gw.f_sigma8(u.E, 1.0)
        fs8_z1 = next((r["f_sigma8"] for r in fs8 if r["z"] == 1.0), float("nan"))
        st, note = _status(ht, gam)
        if u.key == "fluid_n2094":
            st, note = "EXCLUDED", "contains no matter: no galaxies, no structures"
        if u.key == "granular":
            st, note = "ALIVE-BUT-INVISIBLE", "cut-off below any probed scale (see docs)"
        rows.append(Prediction(u.key, u.label, u.family, ht, ht_err,
                               gam, f, fs8_z1, st, note))
    return rows


# ---------------------------------------------------------------------------
# the granular sector
# ---------------------------------------------------------------------------
def granular_sector(rc_over_box_values=(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)) -> dict:
    """Growth index and suppression versus the cut-off scale, in units of r_c."""
    xs = np.logspace(-3, 0.0, 40)
    gamma_x, D_x = [], []
    _, D_gr, _, _ = gw.growth_factor(bg.E_lcdm, 1.0, 0.0)
    for x in xs:
        mu = (lambda xx: (lambda k, a: float(gr.mu_effective(xx))))(x)
        _, D, _, om = gw.growth_factor(bg.E_lcdm, mu, 0.0)
        f = gw.growth_factor(bg.E_lcdm, mu, 0.0)[2] / D
        gamma_x.append(float(np.log(f) / np.log(om)) if f > 0 else float("nan"))
        D_x.append(float(D / D_gr))
    return {
        "x": xs.tolist(),
        "gamma": gamma_x,
        "D_rel": D_x,
        "bounds": gr.bounds_table(),
    }


def observable_window() -> dict:
    """Where the granular signature would have to appear to be observable.

    Euclid / DESI probe roughly 0.01 < k < 0.3 h/Mpc. The cut-off signature
    peaks near x = k r_c ~ 0.3, so an observable deviation requires
    r_c ~ 0.3 / k, i.e. of order a few Mpc: which is excluded by the
    existence of dwarf galaxies. This is the central tension of the
    granular universe and it is reported, not hidden.
    """
    k_obs_lo, k_obs_hi = 0.01, 0.3
    rc_needed_lo = 0.3 / k_obs_hi
    rc_needed_hi = 0.3 / k_obs_lo
    rc_bound_dwarf = gr.rc_bound_from(1.0e-4)          # ultra-faint dwarfs, 100 pc
    rc_bound_micro = gr.rc_bound_from(1.0e-8)          # microhalos, 1e-2 pc
    return {
        "k_observed_range_h_mpc": [k_obs_lo, k_obs_hi],
        "r_c_needed_for_observable_gamma_mpc": [rc_needed_lo, rc_needed_hi],
        "r_c_bound_from_dwarfs_mpc": rc_bound_dwarf,
        "r_c_bound_from_microhalos_mpc": rc_bound_micro,
        "gap_dex_dwarfs": float(np.log10(rc_needed_lo / rc_bound_dwarf)),
        "gap_dex_microhalos": float(np.log10(rc_needed_lo / rc_bound_micro)),
    }
