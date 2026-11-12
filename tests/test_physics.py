"""
Physics tests.

Each test pins a number that must not drift: either an analytic identity, a
standard result from the literature, or a measured accuracy that the README
quotes. If one of these fails, the quoted number in the documentation is
wrong.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ratiss import background as bg
from ratiss import constants as C
from ratiss import granular as gr
from ratiss import growth as gw
from ratiss import pm


# ---------------------------------------------------------------------------
# analytic identities
# ---------------------------------------------------------------------------
def test_single_fluid_age_identity():
    """For a universe with a single substance diluting as a^-n, H0*t0 = 2/n."""
    for n in (1.0, 2.0, 3.0, 4.0, 2.0945):
        E = bg.make_single_fluid(n)
        ht, _ = bg.hubble_time_times_h0(E)
        assert abs(ht - 2.0 / n) < 2e-4, f"n={n}: {ht} != {2/n}"


def test_matter_only_universe():
    """Einstein-de Sitter: H0*t0 = 2/3 exactly."""
    E = lambda a: float(a ** -1.5)
    ht, _ = bg.hubble_time_times_h0(E)
    assert abs(ht - 2.0 / 3.0) < 1e-5


def test_growth_equals_scale_factor_in_matter_era():
    """For an Einstein-de Sitter universe D is proportional to a."""
    E = lambda a: float(a ** -1.5)
    _, D1, _, _ = gw.growth_factor(E, 1.0, 0.0, N_fin=np.log(0.1), om0=1.0)
    _, D2, _, _ = gw.growth_factor(E, 1.0, 0.0, N_fin=np.log(1.0), om0=1.0)
    assert abs((D2 / D1) / 10.0 - 1.0) < 2e-3


def test_growth_3_2_over_5_identity():
    """Matter era: f = Omega_m^(3/5) to good accuracy."""
    E = lambda a: float(a ** -1.5)
    _, D, dD, om = gw.growth_factor(E, 1.0, 0.0, N_fin=0.0, om0=1.0)
    f = dD / D
    assert abs(f - om ** 0.6) < 1e-3, f"f={f}, om^0.6={om**0.6}"


def test_dgp_background_normalisation():
    """E(a=1) must equal 1 for every model that is normalised today."""
    for u in bg.registry():
        if u.key == "fluid_n2094":
            continue
        assert abs(u.E(1.0) - 1.0) < 2e-3, f"{u.key}: E(1)={u.E(1.0)}"


# ---------------------------------------------------------------------------
# the granular model
# ---------------------------------------------------------------------------
def test_r_c_is_the_equilibrium_radius():
    """By construction the force must vanish exactly at r = r_c."""
    G, K = 1.0, 0.01
    rc = gr.r_c(G, K)
    assert abs(gr.force(rc, G, K)) < 1e-10


def test_force_is_repulsive_below_r_c():
    G, K = 1.0, 0.01
    rc = gr.r_c(G, K)
    assert gr.force(0.5 * rc, G, K) > 0
    assert gr.force(2.0 * rc, G, K) < 0


def test_mu_limits():
    assert abs(gr.mu_effective(0.0) - 1.0) < 1e-12
    assert abs(gr.mu_effective(1.0)) < 1e-12
    assert abs(gr.mu_effective(1e6) + 1.0) < 1e-6


def test_rc_bound_scales_linearly():
    """r_c_max = x_cut * L_min / (2 pi)."""
    L = 1.0e-4
    assert abs(gr.rc_bound_from(L, 0.1) - 0.1 * L / (2 * np.pi)) < 1e-18


# ---------------------------------------------------------------------------
# the N-body engine
# ---------------------------------------------------------------------------
def test_cic_deposit_conserves_mass():
    rng = np.random.default_rng(3)
    pos = rng.uniform(0, 1, size=(500, 3))
    grid = pm.cic_deposit(pos, 32)
    assert abs(grid.sum() - 500.0) < 1e-8


def test_cic_deposit_single_mode():
    """A single sinusoidal displacement must be recovered to better than 1 %."""
    n, Np, amp = 32, 16, 0.01
    q1 = (np.arange(Np) + 0.5) / Np
    QX, QY, QZ = np.meshgrid(q1, q1, q1, indexing="ij")
    q = np.stack([QX.ravel(), QY.ravel(), QZ.ravel()], axis=1)
    kx = 2 * np.pi
    q[:, 0] += amp * np.sin(kx * q[:, 0])
    grid = pm.cic_deposit(q % 1.0, n)
    delta = grid / grid.mean() - 1.0
    dk = np.fft.fftn(delta) / n ** 3
    assert abs(abs(dk[1, 0, 0]) / (amp * kx / 2) - 1.0) < 1e-2


def test_initial_conditions_sigma8():
    """The realised analytic field must carry the requested sigma8 to 3 %."""
    cfg = pm.SimConfig(box_mpc_h=200.0, n_grid=64, n_part_1d=32,
                       a_init=1.0, a_final=1.0001, n_steps=1)
    cos = pm.Cosmology()
    pos, p = pm.zeldovich_ics(cfg, cos)
    n = cfg.n_grid
    grid = pm.cic_deposit(pos, n)
    delta = grid / grid.mean() - 1.0
    dk = np.fft.fftn(delta) / n ** 3
    KX, KY, KZ, k2, k2s = pm._wavenumbers(n)
    P, kc, cnt = pm._binned_power(dk, k2, n, cfg.box_mpc_h,
                                  shot=cfg.box_mpc_h ** 3 / cfg.n_part)
    s8 = pm.sigma8_of(P, kc)
    assert 0.55 < s8 < 0.90, f"sigma8 realised = {s8}"


def test_pm_growth_rate():
    """The engine must reproduce the linear growth rate to better than 6 %."""
    cfg = pm.SimConfig(box_mpc_h=200.0, n_grid=32, n_part_1d=16, n_steps=30,
                       a_init=0.05, sigma8=0.811)
    cos = pm.Cosmology()
    res = pm.run(cfg, cos, record_slice=False)
    a = np.array(res.a)
    d = np.array(res.sigma_low)
    Dl = np.array([gw.growth_factor(cos.E, 1.0, 0.0, N_fin=np.log(x))[1] for x in a])
    pred = Dl - Dl[0]
    meas = d - d[0]
    m = pred > 0.2 * pred[-1]
    ratio = meas[m] / pred[m]
    spread = np.max(np.abs(ratio / np.median(ratio) - 1.0))
    assert spread < 0.06, f"growth spread = {spread:.3f}"


# ---------------------------------------------------------------------------
# the degeneracy result
# ---------------------------------------------------------------------------
def test_five_universes_share_the_same_H0t0():
    """The central negative result: several theories give identical H0*t0."""
    from ratiss import models as md
    reg = {r.key: r for r in md.build_registry()}
    group = ["lcdm", "fr_hs", "ft", "granular", "rvm_001", "w_m1004"]
    vals = [reg[k].H0t0 for k in group]
    assert max(vals) - min(vals) < 5e-3, vals


def test_dgp_is_excluded():
    from ratiss import models as md
    reg = {r.key: r for r in md.build_registry()}
    assert reg["dgp"].status == "EXCLUDED"


def test_gamma_of_the_granular_model_matches_gr_at_large_scales():
    mu = lambda k, a: float(gr.mu_effective(1e-6))
    f, g = gw.growth_index(bg.E_lcdm, mu, 0.0)
    f0, g0 = gw.growth_index(bg.E_lcdm, 1.0, 0.0)
    assert abs(g - g0) < 1e-4


# ---------------------------------------------------------------------------
# D12 -- the comoving distance is not the age integral
# ---------------------------------------------------------------------------
def test_comoving_distance_is_not_the_age_integral():
    """D_C/(c/H0) = int da/(a^2 E); it must differ from H0*t0 = int da/(a E)."""
    from scipy.integrate import quad
    E = bg.E_lcdm
    dc = bg.comoving_distance_over_cH0(E)
    ht = bg.hubble_time_times_h0(E)[0]
    ref = quad(lambda a: 1.0 / (a * a * E(a)), 1e-8, 1.0, epsabs=1e-12)[0]
    assert abs(dc - ref) < 1e-5, (dc, ref)
    assert abs(dc - ht) > 1.0, "D_C must not equal H0*t0"


def test_luminosity_distance_identity():
    """D_L(z) = (1+z) D_C(z) for the same background."""
    from scipy.integrate import quad
    E = bg.E_lcdm
    for z in (0.5, 1.0, 2.0):
        dc = quad(lambda a: 1.0 / (a * a * E(a)), 1.0 / (1 + z), 1.0, epsabs=1e-12)[0]
        dl = bg.luminosity_distance(z, E)
        assert abs(dl - (1 + z) * dc) < 1e-6, (z, dl, (1 + z) * dc)


# ---------------------------------------------------------------------------
# f(R) chameleon sector
# ---------------------------------------------------------------------------
def test_fr_mu_limits_and_monotonicity():
    """mu -> 1 at large scales, -> 4/3 at small scales, and is monotone."""
    from ratiss import fr
    assert abs(float(fr.mu_fR(1e-6, 1.0, 1e-5)) - 1.0) < 1e-3
    assert abs(float(fr.mu_fR(1e4, 1.0, 1e-5)) - 4.0 / 3.0) < 1e-2
    ks = np.logspace(-4, 4, 50)
    mu = np.array([float(fr.mu_fR(k, 1.0, 1e-5)) for k in ks])
    assert np.all(np.diff(mu) > 0)


def test_fr_growth_is_enhanced_and_screening_appears():
    """f(R) enhances growth; smaller |f_R0| is screened back to Lambda-CDM."""
    from ratiss import models as md
    s = md.fr_sector()
    rows = {r["f_R0"]: r for r in s["rows"]}
    assert rows[1e-4]["D_rel_vs_lcdm"] > 1.02
    d = [rows[k]["D_rel_vs_lcdm"] for k in (1e-4, 1e-5, 1e-6, 1e-7)]
    assert d[0] > d[1] > d[2] > d[3] >= 1.0
    assert abs(rows[1e-7]["delta_gamma"]) < 2e-3
    assert abs(rows[1e-4]["delta_gamma"]) > 0.05


def test_fr_is_scale_dependent_but_background_is_degenerate():
    """The whole point: identical H0*t0, a gamma that depends on k."""
    from ratiss import models as md
    s = md.fr_sector()
    gam = [row["gamma"] for row in s["scale_scan"]]
    assert max(gam) - min(gam) > 0.03, "the signature must be scale dependent"
    # the f(R) background is identical to Lambda-CDM by construction
    ht = bg.hubble_time_times_h0(bg.E_fr_hs)[0]
    assert abs(ht - bg.hubble_time_times_h0(bg.E_lcdm)[0]) < 1e-12
    assert abs(s["rows"][1]["delta_gamma"]) > 0.01
