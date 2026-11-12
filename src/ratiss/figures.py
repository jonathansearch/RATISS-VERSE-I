"""
Figure generation for the RATISS reference universe.

Every figure is produced from numeric arrays computed at build time; no
value is ever hand-typed. Matplotlib only renders.
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import background as bg
from . import constants as C
from . import granular as gr
from . import growth as gw

plt.rcParams.update({
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "figure.dpi": 130,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.28,
})

C_BG = "#1b2a41"
C_ACC = "#c9a227"
C_RED = "#b03a2e"
C_BLUE = "#2c6e9b"
C_GREEN = "#2e7d5b"
C_PUR = "#6d4c8a"
C_ORG = "#c9702a"
C_GREY = "#7b8794"


def _save(fig, path):
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
def fig01_background(reg, path):
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    a = np.linspace(0.05, 1.0, 400)
    keys = ["lcdm", "w_m0302", "dgp", "cardassian_05", "fluid_n2094"]
    for u in bg.registry():
        if u.key not in keys:
            continue
        e = [u.E(x) for x in a]
        ax[0].plot(a, e, lw=2.0, label=u.label)
    ax[0].set_xlabel("scale factor  a")
    ax[0].set_ylabel("E(a) = H(a)/H0")
    ax[0].set_ylim(0, 12)
    ax[0].set_title("Expansion histories: one solver, every universe")
    ax[0].legend()

    labels = [r.label for r in reg]
    vals = [r.H0t0 for r in reg]
    y = np.arange(len(reg))
    cols = [C_GREEN if r.status == "DEGENERATE" else
            (C_ORG if r.status.startswith("ALIVE") else C_RED) for r in reg]
    ax[1].barh(y, vals, color=cols, edgecolor="k", lw=0.6)
    ax[1].axvline(vals[0], ls="--", color="k", lw=1.4)
    ax[1].set_yticks(y)
    ax[1].set_yticklabels(labels, fontsize=7.5)
    ax[1].invert_yaxis()
    ax[1].set_xlabel("H0 * t0  (dimensionless age)")
    ax[1].set_title("Five universes share the same H0*t0")
    for i, v in enumerate(vals):
        ax[1].text(v + 0.01, i, f"{v:.4f}", va="center", fontsize=7)
    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig02_growth(reg, path):
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    a = np.logspace(-2, 0, 60)
    show = [("lcdm", C_BG, "-"), ("w_m0302", C_ORG, "--"),
            ("dgp", C_PUR, "-."), ("cardassian_05", C_BLUE, ":")]
    lookup = {u.key: u for u in bg.registry()}
    for key, col, ls in show:
        u = lookup[key]
        D, f = [], []
        _, D0, _, _ = gw.growth_factor(u.E, 1.0, 0.0, N_fin=0.0)
        for av in a:
            _, Dv, dDv, _ = gw.growth_factor(u.E, 1.0, 0.0, N_fin=np.log(av))
            D.append(Dv / D0)
            f.append(dDv / Dv if Dv else np.nan)
        ax[0].plot(a, D, lw=2.0, color=col, ls=ls, label=u.label)
        ax[1].plot(a, f, lw=2.0, color=col, ls=ls, label=u.label)
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_xlabel("a")
    ax[0].set_ylabel("D(a)/D(1)")
    ax[0].set_title("Linear growth factor")
    ax[1].set_xlabel("a")
    ax[1].set_ylabel("f(a) = dlnD/dlna")
    ax[1].set_title("Growth rate: the models separate here")
    ax[0].legend()
    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig03_discriminator(reg, path):
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    reg = sorted(reg, key=lambda r: r.gamma if np.isfinite(r.gamma) else -1)
    y = np.arange(len(reg))
    g = [r.gamma if np.isfinite(r.gamma) else 0.0 for r in reg]
    cols = [C_GREEN if r.status == "DEGENERATE" else
            (C_ORG if r.status.startswith("ALIVE") else C_RED) for r in reg]
    ax[0].barh(y, g, color=cols, edgecolor="k", lw=0.6)
    ax[0].axvline(g[0] if g else 0.545, ls="--", color="k", lw=1.3)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([r.label for r in reg], fontsize=7.5)
    ax[0].set_xlim(0.50, 0.62)
    ax[0].set_xlabel("growth index  gamma = ln f / ln Omega_m")
    ax[0].set_title("Growth index: the only discriminator found")
    for i, v in enumerate(g):
        ax[0].text(v + 0.002, i, f"{v:.4f}", va="center", fontsize=7)

    # relative f sigma8 (normalised to the reference universe)
    z = np.linspace(0, 2, 60)
    ref = {r["z"]: r["f_sigma8"] for r in gw.f_sigma8(bg.E_lcdm)}
    ax[1].axhline(0, color="k", lw=1.2)
    for key, col, ls in [("w_m0302", C_ORG, "--"), ("dgp", C_PUR, "-."),
                         ("cardassian_05", C_BLUE, ":")]:
        u = {x.key: x for x in bg.registry()}[key]
        vals = []
        for zz in z:
            row = [r for r in gw.f_sigma8(u.E, redshifts=(float(zz),))][0]
            r0 = [r for r in gw.f_sigma8(bg.E_lcdm, redshifts=(float(zz),))][0]
            vals.append(100.0 * (row["f_sigma8"] / r0["f_sigma8"] - 1.0))
        ax[1].plot(z, vals, lw=2.0, color=col, ls=ls, label=u.label)
    ax[1].set_xlabel("redshift z")
    ax[1].set_ylabel("f sigma8 relative to Lambda-CDM  [%]")
    ax[1].set_title("Predicted deviation in f sigma8")
    ax[1].legend()
    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig04_granular(gran, path):
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    x = np.logspace(-3, 1, 400)
    ax[0].plot(x, gr.mu_effective(x), lw=2.6, color=C_BLUE)
    ax[0].axhline(1, ls="--", color=C_GREEN, lw=1.3, label="Newton / GR")
    ax[0].axhline(0, color="k", lw=0.9)
    ax[0].fill_between(x, gr.mu_effective(x), 0,
                       where=gr.mu_effective(x) < 0, color=C_RED, alpha=0.22)
    ax[0].axvline(1.0, ls=":", color=C_RED, lw=1.6)
    ax[0].text(1.25, -0.55, "repulsion", color=C_RED, fontsize=9)
    ax[0].text(0.012, 0.55, "Newton", color=C_GREEN, fontsize=9)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("x = k * r_c")
    ax[0].set_ylabel("mu(x)")
    ax[0].set_title("Effective coupling of the granular universe")
    ax[0].legend(loc="upper right")

    xs = np.array(gran["x"])
    gam = np.array(gran["gamma"])
    D = np.array(gran["D_rel"])
    ax[1].plot(xs, gam, lw=2.4, color=C_PUR, label="gamma(x)")
    ax[1].axhline(0.5543, ls="--", color=C_GREEN, lw=1.4, label="GR: 0.554")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("x = k * r_c")
    ax[1].set_ylabel("growth index gamma", color=C_PUR)
    ax2 = ax[1].twinx()
    ax2.semilogy(xs, D, lw=2.0, color=C_ORG, label="D/D_GR")
    ax2.set_ylabel("growth suppression D/D_GR", color=C_ORG)
    ax2.grid(False)
    ax[1].set_title("Cut-off: growth deviates, then structures vanish")
    ax[1].legend(loc="upper left")
    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig05_observability(window, path):
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.set_xscale("log")
    lo, hi = window["r_c_needed_for_observable_gamma_mpc"]
    ax.axvspan(lo, hi, color=C_BLUE, alpha=0.22,
               label="r_c needed for an observable gamma")
    ax.axvline(window["r_c_bound_from_dwarfs_mpc"], color=C_RED, lw=2.2,
               label="upper bound from ultra-faint dwarfs")
    ax.axvline(window["r_c_bound_from_microhalos_mpc"], color=C_RED, lw=2.2, ls="--",
               label="upper bound from microhalos (predicted)")
    ax.set_xlim(1e-11, 1e3)
    ax.set_ylim(0, 1)
    ax.set_xlabel("r_c  [Mpc]")
    ax.set_yticks([])
    ax.text(lo * 1.15, 0.72, f"observable window\n{lo:.0f} - {hi:.0f} Mpc",
            color=C_BLUE, fontsize=9.5)
    ax.text(window["r_c_bound_from_dwarfs_mpc"] * 1.6, 0.42,
            f"allowed by data\n$r_c < {window['r_c_bound_from_dwarfs_mpc']:.1e}$ Mpc",
            color=C_RED, fontsize=9.5)
    gap = window["gap_dex_dwarfs"]
    ax.annotate("", xy=(lo, 0.16), xytext=(window["r_c_bound_from_dwarfs_mpc"], 0.16),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.6))
    ax.text(np.sqrt(lo * window["r_c_bound_from_dwarfs_mpc"]), 0.20,
            f"{gap:.1f} decades of gap", ha="center", fontsize=10, fontweight="bold")
    ax.set_title("The central tension of the granular universe\n"
                 "alive, but its signature lies far outside any probed scale")
    ax.legend(loc="upper left")
    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig06_simulation(sim, path):
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    field = np.array(sim["slice_field"])
    im = ax[0].imshow(np.log10(np.abs(field) + 1e-2), cmap="inferno",
                      origin="lower")
    a_fin = sim["config"].get("a_final", 1.0)
    ax[0].set_title(f"Particle-mesh density slice, a = {a_fin:.2f}")
    ax[0].set_xlabel("grid cell")
    ax[0].set_ylabel("grid cell")
    ax[0].grid(False)
    fig.colorbar(im, ax=ax[0], label="log10 |delta + 1|")

    pk = sim["pk"]
    k = np.array(pk["k_h_mpc"])
    P = np.array(pk["P"])
    ax[1].loglog(k, P, "o-", ms=3.5, lw=1.4, color=C_RED,
                 label="RATISS PM simulation")
    kk = np.logspace(-3, np.log10(k.max()), 300)
    Plin = pm_shape(kk)
    ax[1].loglog(kk, Plin, "-", lw=1.6, color=C_BG, label="linear theory (BBKS)")
    ax[1].set_xlabel("k  [h/Mpc]")
    ax[1].set_ylabel("P(k)  [(Mpc/h)^3]")
    ax[1].set_title("Matter power spectrum at z = 0")
    ax[1].legend()
    fig.tight_layout()
    _save(fig, path)


def pm_shape(k):
    from . import pm
    return pm.normalise_to_sigma8(pm.power_spectrum_shape(k), k, C.SIGMA_8)


def fig08_fr(sector, path):
    """Hu--Sawicki f(R): a genuine, scale-dependent, testable signature.

    Left  : mu(k) and gamma(k) at fixed |f_R0| -- the shape of the signature.
    Right : the growth factor today versus |f_R0|, with the screening limit.
    """
    rows = sector["rows"]
    fR0 = np.array([r["f_R0"] for r in rows])
    drel = np.array([r["D_rel_vs_lcdm"] for r in rows])
    dgam = np.array([r["delta_gamma"] for r in rows])
    scan = sector["scale_scan"]
    k = np.array([s["k_h_mpc"] for s in scan])
    mu = np.array([s["mu"] for s in scan])
    gam = np.array([s["gamma"] for s in scan])

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))

    ax[0].set_xscale("log")
    ax[0].plot(k, mu, lw=2.4, color=C_PUR, label=r"$\mu(k)=G_{\rm eff}/G$")
    ax[0].axhline(1.0, color="k", lw=1.0, ls=":")
    ax[0].axhline(4.0 / 3.0, color=C_RED, lw=1.0, ls="--")
    ax[0].text(k[0] * 1.2, 4.0 / 3.0 + 0.012, r"$4/3$ (unscreened)",
               color=C_RED, fontsize=8.5)
    ax[0].text(k[0] * 1.2, 1.0 + 0.012, r"$1$ (GR)", color="k", fontsize=8.5)
    ax[0].set_xlabel(r"$k$  [$h$/Mpc]")
    ax[0].set_ylabel(r"$\mu(k)$")
    ax2 = ax[0].twinx()
    ax2.plot(k, gam, lw=2.0, color=C_BLUE, ls="-.", label=r"$\gamma(k)$")
    ax2.axhline(sector["gamma_lcdm_at_k_ref"], color=C_GREY, lw=1.0, ls=":")
    ax2.set_ylabel(r"growth index  $\gamma$", color=C_BLUE)
    ax2.grid(False)
    ax[0].set_title(r"$|f_{R0}|=10^{-5}$: the coupling is scale dependent")
    l1, la1 = ax[0].get_legend_handles_labels()
    l2, la2 = ax2.get_legend_handles_labels()
    ax[0].legend(l1 + l2, la1 + la2, loc="center right")

    ax[1].set_xscale("log")
    ax[1].plot(fR0, drel, "o-", lw=2.2, color=C_GREEN,
               label=r"$D(a{=}1)/D_{\Lambda\rm CDM}$")
    ax[1].axhline(1.0, color="k", lw=1.0, ls=":")
    ax[1].set_xlabel(r"$|f_{R0}|$")
    ax[1].set_ylabel(r"$D/D_{\Lambda\rm CDM}$")
    axb = ax[1].twinx()
    axb.plot(fR0, np.abs(dgam), "s--", lw=1.8, color=C_ORG, label=r"$|\Delta\gamma|$")
    axb.axhline(0.01, color=C_RED, lw=1.0, ls="--")
    axb.text(fR0[0], 0.011, "DESI / Euclid reach", color=C_RED, fontsize=8.5)
    axb.set_yscale("log")
    axb.set_ylabel(r"$|\Delta\gamma|$", color=C_ORG)
    axb.grid(False)
    ax[1].set_title("Screening: large $|f_{R0}|$ is testable, small is invisible")
    l1, la1 = ax[1].get_legend_handles_labels()
    l2, la2 = axb.get_legend_handles_labels()
    ax[1].legend(l1 + l2, la1 + la2, loc="upper right")

    fig.tight_layout()
    _save(fig, path)


# ---------------------------------------------------------------------------
def fig07_architecture(path):
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)

    boxes = [
        (0.2, 4.4, 2.2, 1.2, "constants.py\nPlanck 2018", C_GREY),
        (2.9, 4.4, 2.2, 1.2, "background.py\nE(a), H0*t0", C_BLUE),
        (2.9, 2.7, 2.2, 1.2, "growth.py\ngamma, f sigma8", C_PUR),
        (2.9, 1.0, 2.2, 1.2, "granular.py\nmu(k), r_c", C_ORG),
        (5.6, 2.7, 2.0, 1.2, "models.py\nuniverse registry", C_GREEN),
        (5.6, 4.4, 2.0, 1.2, "pm.py\nN-body engine", C_RED),
        (8.0, 3.55, 1.9, 1.2, "build.py\nresults/ + figures/", C_BG),
    ]
    for (x, y, w, h, txt, col) in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=col, ec="k", lw=1.1, alpha=0.88))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                color="w", fontsize=8.6, fontweight="bold")

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", lw=1.5, color="k"))

    arrow(2.4, 5.0, 2.9, 5.0)
    arrow(4.0, 4.4, 4.0, 3.9)
    arrow(5.1, 3.3, 5.6, 3.3)
    arrow(5.1, 1.6, 5.6, 2.7)
    arrow(5.1, 5.0, 5.6, 5.0)
    arrow(6.6, 3.9, 7.4, 4.2)
    arrow(7.6, 4.6, 8.0, 4.3)
    arrow(7.6, 5.0, 8.0, 5.0)
    ax.text(5.0, 0.35, "one solver per question, every number recomputed at build time",
            ha="center", fontsize=9.5, style="italic")
    ax.set_title("RATISS architecture", fontsize=12)
    fig.tight_layout()
    _save(fig, path)
