"""
Particle-mesh N-body engine with periodic boundary conditions.

This is what makes RATISS a simulator rather than a toy. The design encodes
the single most important lesson of the parent project:

    In a periodic box the mean density exerts no force (the Jeans swindle).
    The expansion therefore *cannot* be produced by the particles.
    a(t) is supplied by the Friedmann equation; the particles produce only
    the departure from homogeneity.

That is how every production N-body code works, and it is stated here
explicitly rather than hidden.

Units and equations
-------------------
Comoving positions x in box units, x in [0,1). Time t in units of 1/H0.
Superconformal time tau with dtau = dt/a^2 gives the cleanest leapfrog
(Hockney & Eastwood 1988, ch. 8):

    dx/dtau = p
    dp/dtau = -grad(Phi)/a^2
    lap(Phi) = (3/2) Omega_m delta / a

In Fourier space with code wavenumbers k (box = 1):

    F_k = i (3/2) Omega_m a delta_k k_vec / |k|^2 * mu(k)

``mu(k)`` is where modified gravity enters: it multiplies the effective
coupling mode by mode, so the same engine integrates General Relativity
(mu = 1) and the r_c cut-off alike.

Power-spectrum normalisation
----------------------------
With delta_k^FFT = N^3 sqrt(P_code(k)) z and P_code = P_phys / box^3, the
inverse FFT reproduces <delta^2> = int P d^3k/(2 pi)^3 exactly. The
implementation below is checked against linear theory in tests/.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from . import constants as C


# ---------------------------------------------------------------------------
# transfer function
# ---------------------------------------------------------------------------
def transfer_bbks(k_h_mpc, gamma: float) -> np.ndarray:
    """BBKS transfer function (Bardeen, Bond, Kaiser & Szalay 1986). k in h/Mpc."""
    q = np.maximum(np.asarray(k_h_mpc, dtype=float) / gamma, 1e-12)
    t = np.log1p(2.34 * q) / (2.34 * q)
    t *= (1.0 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25
    return t


def power_spectrum_shape(k_h_mpc, n_s: float = C.N_S, gamma: float | None = None) -> np.ndarray:
    """P(k) proportional to k^n_s T^2(k). k in h/Mpc."""
    if gamma is None:
        gamma = C.OMEGA_M * C.H0_OVER_100
    k = np.maximum(np.asarray(k_h_mpc, dtype=float), 1e-12)
    return k ** n_s * transfer_bbks(k, gamma) ** 2


def sigma8_of(P: np.ndarray, k: np.ndarray, R: float = 8.0) -> float:
    """sigma_R from a tabulated P(k), k in h/Mpc, R in Mpc/h."""
    x = k * R
    W = np.where(x > 1e-8, 3.0 * (np.sin(x) - x * np.cos(x)) / np.maximum(x, 1e-30) ** 3, 1.0)
    d2 = k ** 3 * P / (2.0 * np.pi ** 2)
    return float(np.sqrt(np.trapezoid(d2 * W ** 2, np.log(k))))


def normalise_to_sigma8(P_shape: np.ndarray, k: np.ndarray, sigma8: float) -> np.ndarray:
    return P_shape * (sigma8 / sigma8_of(P_shape, k)) ** 2


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
@dataclass
class Cosmology:
    omega_m: float = C.OMEGA_M
    omega_l: float = C.OMEGA_LAMBDA
    h: float = C.H0_OVER_100

    def E(self, a: float) -> float:
        return float(np.sqrt(self.omega_m * a ** -3 + self.omega_l))

    def dlnH_dN(self, a: float) -> float:
        e2 = self.omega_m * a ** -3 + self.omega_l
        return -1.5 * self.omega_m * a ** -3 / e2


@dataclass
class SimConfig:
    box_mpc_h: float = 200.0
    n_grid: int = 64
    n_part_1d: int = 32
    a_init: float = 1.0 / 50.0
    a_final: float = 1.0
    n_steps: int = 40
    seed: int = 20261030
    sigma8: float = C.SIGMA_8
    kcut_h_mpc: float = 0.20
    mu_k: Callable[[np.ndarray], np.ndarray] | None = None

    @property
    def n_part(self) -> int:
        return self.n_part_1d ** 3


# ---------------------------------------------------------------------------
# grid helpers (code units, box = 1)
# ---------------------------------------------------------------------------
def _wavenumbers(n: int):
    f = np.fft.fftfreq(n, d=1.0 / n)      # integers
    kx = 2.0 * np.pi * f                   # code wavenumbers, box = 1
    KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    k2_safe = np.where(k2 > 0, k2, 1.0)
    return KX, KY, KZ, k2, k2_safe


def cic_deposit(pos: np.ndarray, n: int) -> np.ndarray:
    g = pos * n
    i0 = np.floor(g).astype(np.int64)
    d = g - i0
    i0 %= n
    i1 = (i0 + 1) % n
    wx = np.array([1.0 - d[:, 0], d[:, 0]])
    wy = np.array([1.0 - d[:, 1], d[:, 1]])
    wz = np.array([1.0 - d[:, 2], d[:, 2]])
    ix = np.stack([i0[:, 0], i1[:, 0]], axis=1)
    iy = np.stack([i0[:, 1], i1[:, 1]], axis=1)
    iz = np.stack([i0[:, 2], i1[:, 2]], axis=1)
    grid = np.zeros((n, n, n), dtype=np.float64)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                np.add.at(grid, (ix[:, a], iy[:, b], iz[:, c]),
                          wx[a] * wy[b] * wz[c])
    return grid


def cic_gather(field: np.ndarray, pos: np.ndarray, n: int) -> np.ndarray:
    g = pos * n
    i0 = np.floor(g).astype(np.int64)
    d = g - i0
    i0 %= n
    i1 = (i0 + 1) % n
    wx = np.array([1.0 - d[:, 0], d[:, 0]])
    wy = np.array([1.0 - d[:, 1], d[:, 1]])
    wz = np.array([1.0 - d[:, 2], d[:, 2]])
    ix = np.stack([i0[:, 0], i1[:, 0]], axis=1)
    iy = np.stack([i0[:, 1], i1[:, 1]], axis=1)
    iz = np.stack([i0[:, 2], i1[:, 2]], axis=1)
    out = np.zeros(len(pos), dtype=np.float64)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                out += wx[a] * wy[b] * wz[c] * field[ix[:, a], iy[:, b], iz[:, c]]
    return out


def _sample(vec_field: np.ndarray, pos: np.ndarray, n: int) -> np.ndarray:
    return np.stack([cic_gather(vec_field[..., c], pos, n) for c in range(3)], axis=1)


# ---------------------------------------------------------------------------
# initial conditions: Zel'dovich approximation
# ---------------------------------------------------------------------------
def zeldovich_ics(cfg: SimConfig, cosmo: Cosmology, sigma8_now: float | None = None,
                  return_extras: bool = False):
    """Return (pos, p) at a_init.

    x = q + D psi,  p = dx/dtau = a^2 H f D psi
    Psi_k = i k delta_k / k^2
    """
    n, box = cfg.n_grid, cfg.box_mpc_h
    rng = np.random.default_rng(cfg.seed)
    KX, KY, KZ, k2, k2_safe = _wavenumbers(n)

    # --- tabulate the shape, normalise sigma8 exactly ---------------------
    kh_tab = np.logspace(-4, 2, 4000)
    P_tab = normalise_to_sigma8(power_spectrum_shape(kh_tab),
                                kh_tab, sigma8_now or cfg.sigma8)

    # --- code-space wavenumber magnitude, converted to h/Mpc ------------
    k_code = np.sqrt(k2)
    kh = k_code / box

    P_on_grid = np.interp(kh.ravel(), kh_tab, P_tab).reshape(n, n, n)
    P_on_grid[0, 0, 0] = 0.0

    # --- complex Gaussian with <|z|^2> = 1 --------------------------------
    zr = rng.normal(size=(n, n, n))
    zi = rng.normal(size=(n, n, n))
    z = (zr + 1j * zi) / np.sqrt(2.0)

    # delta_k^FFT = N^3 sqrt(P_code) z,  P_code = P_phys / box^3
    dk = n ** 3 * np.sqrt(P_on_grid / box ** 3) * z
    dk[0, 0, 0] = 0.0

    # z is a complex Gaussian with no Hermitian symmetry, so ifftn(dk) is
    # complex. Projecting onto the real axis halves the power; the factor
    # sqrt(2) restores <delta^2> = int P d^3k/(2 pi)^3 exactly (checked in
    # tests/test_physics.py::test_ic_power_spectrum).
    root2 = np.sqrt(2.0)
    psi_x = root2 * np.fft.ifftn(1j * KX / k2_safe * dk).real
    psi_y = root2 * np.fft.ifftn(1j * KY / k2_safe * dk).real
    psi_z = root2 * np.fft.ifftn(1j * KZ / k2_safe * dk).real
    psi = np.stack([psi_x, psi_y, psi_z], axis=-1)

    # --- lattice ---------------------------------------------------------
    q1 = (np.arange(cfg.n_part_1d) + 0.5) / cfg.n_part_1d
    QX, QY, QZ = np.meshgrid(q1, q1, q1, indexing="ij")
    q = np.stack([QX.ravel(), QY.ravel(), QZ.ravel()], axis=1)
    # random sub-cell offset: a perfect lattice sits exactly on CIC cell
    # boundaries, which is a degenerate configuration (d = 0, NGP weights).
    q = (q + rng.uniform(0.0, 1.0, size=3)) % 1.0

    from .growth import growth_factor
    _, D_i, dD_i, _ = growth_factor(cosmo.E, 1.0, 0.0, N_fin=np.log(cfg.a_init))
    _, D_0, _, _ = growth_factor(cosmo.E, 1.0, 0.0, N_fin=0.0)
    Drel = D_i / D_0
    f_i = dD_i / D_i

    disp = _sample(psi, q, n)
    pos = (q + Drel * disp) % 1.0
    p = cfg.a_init ** 2 * cosmo.E(cfg.a_init) * f_i * Drel * disp
    if return_extras:
        return pos, p, q, disp
    return pos, p


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------
@dataclass
class SimulationResult:
    a: list = field(default_factory=list)
    rms_delta: list = field(default_factory=list)
    sigma_low: list = field(default_factory=list)
    pk: dict = field(default_factory=dict)
    slice_field: list | None = None
    meta: dict = field(default_factory=dict)


def run(cfg: SimConfig, cosmo: Cosmology, record_slice: bool = True) -> SimulationResult:
    n, box = cfg.n_grid, cfg.box_mpc_h
    KX, KY, KZ, k2, k2_safe = _wavenumbers(n)

    mu_fac = np.ones((n, n, n), dtype=float)
    if cfg.mu_k is not None:
        k_code = np.sqrt(k2)
        mu_fac = np.asarray(cfg.mu_k(k_code / box), dtype=float) * np.ones((n, n, n))
        mu_fac[0, 0, 0] = 1.0

    pos, p, q_lattice, delta0 = zeldovich_ics(cfg, cosmo, return_extras=True)

    ln_a = np.linspace(np.log(cfg.a_init), np.log(cfg.a_final), cfg.n_steps + 1)
    res = SimulationResult()
    res.meta = {
        "box_mpc_h": box, "n_grid": n, "n_part": cfg.n_part,
        "a_init": cfg.a_init, "a_final": cfg.a_final, "n_steps": cfg.n_steps,
        "seed": cfg.seed, "omega_m": cosmo.omega_m, "sigma8": cfg.sigma8,
    }

    def kick_field(p_in, a):
        grid = cic_deposit(p_in, n)
        delta = grid / grid.mean() - 1.0
        dk = np.fft.fftn(delta)
        pref = 1.5 * cosmo.omega_m * a
        fk = 1j * pref * dk * mu_fac / k2_safe
        fx = np.fft.ifftn(fk * KX).real
        fy = np.fft.ifftn(fk * KY).real
        fz = np.fft.ifftn(fk * KZ).real
        return np.stack([fx, fy, fz], axis=-1)

    a = cfg.a_init
    for step in range(cfg.n_steps):
        dln = ln_a[step + 1] - ln_a[step]
        a_mid = float(np.exp(0.5 * (ln_a[step] + ln_a[step + 1])))
        dtau = dln / (a_mid ** 2 * cosmo.E(a_mid))

        F = kick_field(pos, a_mid)
        acc = _sample(F, pos, n) * dtau

        p = p + 0.5 * acc
        pos = (pos + p * dtau) % 1.0
        F = kick_field(pos, a_mid)
        acc = _sample(F, pos, n) * dtau
        p = p + 0.5 * acc

        a = float(np.exp(ln_a[step + 1]))
        grid = cic_deposit(pos, n)
        delta = grid / grid.mean() - 1.0
        res.a.append(a)
        res.rms_delta.append(float(np.sqrt((delta ** 2).mean())))
        # growth measured from the particle displacement, the exact Zel'dovich
        # relation: x - q = (D(a) - D(a_i)) psi. No shot noise, no decorrelation.
        d = pos - q_lattice
        d -= np.round(d)
        res.sigma_low.append(float(np.sqrt((d ** 2).sum(axis=1).mean())))

    grid = cic_deposit(pos, n)
    delta = grid / grid.mean() - 1.0
    dk = np.fft.fftn(delta) / n ** 3
    P, kc, cnt = _binned_power(dk, k2, n, box, shot=box ** 3 / cfg.n_part)
    res.pk = {"k_h_mpc": kc.tolist(), "P": P.tolist(), "nmodes": cnt.tolist()}
    if record_slice:
        res.slice_field = delta[n // 2, :, :].tolist()
    return res


def _binned_power(dk, k2, n, box, shot=0.0):
    p = np.abs(dk) ** 2 * box ** 3
    kmax = np.pi * n / box
    nb = 32
    edges = np.logspace(np.log10(2 * np.pi / box), np.log10(kmax), nb + 1)
    idx = np.digitize((np.sqrt(k2) / box).ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nb)
    P = np.bincount(idx[ok], weights=p.ravel()[ok], minlength=nb)
    cnt = np.bincount(idx[ok], minlength=nb)
    kc = 0.5 * (edges[1:] + edges[:-1])
    good = cnt > 0
    P = P[good] / cnt[good]
    return np.maximum(P - shot, 0.0), kc[good], cnt[good]
