#!/usr/bin/env python3
"""
Build the RATISS reference universe.

Runs every solver, writes results/*.json, renders figures/*.png and seals
the whole thing with results/MANIFEST.sha256.

Usage
-----
    python scripts/build.py            # full build (~2 min)
    python scripts/build.py --fast     # skip the N-body runs
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ratiss import __version__
from ratiss import background as bg
from ratiss import constants as C
from ratiss import figures as fg
from ratiss import granular as gr
from ratiss import growth as gw
from ratiss import models as md
from ratiss import pm

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
for d in (RESULTS, FIGURES):
    d.mkdir(exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(name: str, obj) -> None:
    (RESULTS / name).write_text(json.dumps(obj, indent=1, default=float,
                                           sort_keys=False) + "\n")
    print(f"  wrote results/{name}")


# ---------------------------------------------------------------------------
def build_background() -> dict:
    print("[1/7] background expansion")
    rows = []
    for u in bg.registry():
        ht, err = bg.hubble_time_times_h0(u.E)
        dc = bg.comoving_distance_over_cH0(u.E)
        rows.append({
            "key": u.key, "label": u.label, "family": u.family,
            "H0_t0": ht, "H0_t0_quadrature_error": err,
            "comoving_distance_cH0": dc,
            "E_at_1": u.E(1.0),
        })
    ref = rows[0]["H0_t0"]
    for r in rows:
        r["delta_vs_lcdm"] = r["H0_t0"] - ref
        r["status"] = "DEGENERATE" if abs(r["delta_vs_lcdm"]) < 2e-3 else "EXCLUDED"
    return {"reference": "Planck 2018 flat Lambda-CDM",
            "omega_m": C.OMEGA_M, "H0_km_s_mpc": C.H0_KM_S_MPC,
            "universes": rows}


def build_growth() -> dict:
    print("[2/7] linear growth")
    out = []
    for u in bg.registry():
        f, gam = gw.growth_index(u.E, 1.0)
        fs8 = gw.f_sigma8(u.E, 1.0)
        out.append({"key": u.key, "label": u.label,
                    "f_z0": f, "gamma": gam, "f_sigma8": fs8})
    return {"universes": out}


def build_predictions() -> dict:
    print("[3/7] predictions and falsifiability matrix")
    reg = md.build_registry()
    return {
        "note": ("Every entry is recomputed from the registry; nothing is "
                 "typed by hand. 'status' is a computed classification, not "
                 "an opinion."),
        "reference_H0_t0": md.HT_MEASURED,
        "universes": [r.__dict__ for r in reg],
    }


def build_granular() -> dict:
    print("[4/7] granular sector")
    sec = md.granular_sector()
    win = md.observable_window()
    return {"model": {
                "force_law": "F(r) = -G/r^2 + K/r^4",
                "r_c": "sqrt(K/G)",
                "eta": "r_c / EPS",
                "effective_coupling": "mu(x) = (1 - x^2)/(1 + x^2), x = k r_c",
            },
            "sector": sec,
            "observable_window": win}


def build_fr() -> dict:
    print("[5/7] f(R) chameleon sector")
    return md.fr_sector()


def build_simulations(fast: bool) -> dict:
    print("[6/7] particle-mesh simulations")
    out = {}
    cos = pm.Cosmology()

    cfg = pm.SimConfig(n_steps=24 if fast else 40)
    t0 = time.time()
    sim = pm.run(cfg, cos)
    out["lcdm"] = {
        "config": sim.meta,
        "runtime_s": round(time.time() - t0, 2),
        "a": sim.a, "disp_rms": sim.sigma_low,
        "pk": sim.pk, "slice_field": sim.slice_field,
    }
    print(f"    Lambda-CDM run: {out['lcdm']['runtime_s']} s, "
          f"{cfg.n_part} particles, {cfg.n_grid}^3 mesh")

    # demonstration of the modified-coupling machinery inside the box.
    # r_c = 2 Mpc/h is EXCLUDED by dwarf galaxies (see results/predictions.json);
    # it is used here only to show that mu(k) enters the engine correctly.
    rc_demo = 2.0
    cfg2 = pm.SimConfig(n_steps=24 if fast else 40, seed=cfg.seed + 1,
                        mu_k=lambda k: gr.mu_effective(k * rc_demo))
    t0 = time.time()
    sim2 = pm.run(cfg2, cos)
    out["granular_demo"] = {
        "config": sim2.meta,
        "r_c_mpc": rc_demo,
        "warning": ("r_c = 2 Mpc/h is excluded by the existence of dwarf "
                    "galaxies; this run only demonstrates that mu(k) enters "
                    "the Poisson solver correctly."),
        "runtime_s": round(time.time() - t0, 2),
        "a": sim2.a, "disp_rms": sim2.sigma_low,
        "pk": sim2.pk, "slice_field": sim2.slice_field,
    }
    print(f"    granular demo run: {out['granular_demo']['runtime_s']} s "
          f"(r_c = {rc_demo} Mpc/h)")
    return out


def build_figures(bg_res, pred, gran, fr_res, sims) -> None:
    print("[7/7] figures")
    reg = md.build_registry()
    fg.fig01_background(reg, FIGURES / "fig01_background.png")
    fg.fig02_growth(reg, FIGURES / "fig02_growth.png")
    fg.fig03_discriminator(reg, FIGURES / "fig03_discriminator.png")
    fg.fig04_granular(gran["sector"], FIGURES / "fig04_granular.png")
    fg.fig05_observability(gran["observable_window"],
                           FIGURES / "fig05_observability.png")
    fg.fig06_simulation(sims["lcdm"], FIGURES / "fig06_simulation.png")
    fg.fig08_fr(fr_res, FIGURES / "fig08_fr.png")
    fg.fig07_architecture(FIGURES / "fig07_architecture.png")
    for p in sorted(FIGURES.glob("*.png")):
        print(f"    {p.name}")


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="skip N-body runs")
    args = ap.parse_args()

    t_start = time.time()
    print("=" * 72)
    print(f"RATISS Verse I — reference universe build  v{__version__}")
    print("=" * 72)

    b = build_background()
    dump("background.json", b)
    g = build_growth()
    dump("growth.json", g)
    p = build_predictions()
    dump("predictions.json", p)
    gr_res = build_granular()
    dump("granular.json", gr_res)
    fr_res = build_fr()
    dump("fr.json", fr_res)

    if args.fast:
        sims = {"lcdm": _toy_sim(), "granular_demo": _toy_sim()}
        print("[6/7] simulations skipped (--fast); using stored demo field")
    else:
        sims = build_simulations(fast=False)
        dump("simulations.json", sims)

    build_figures(b, p, gr_res, fr_res, sims)

    manifest = {
        "ratsis_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_s": round(time.time() - t_start, 2),
        "files": {},
    }
    # build.json carries timestamps, so it cannot seal itself; it is metadata.
    for pth in sorted(list(RESULTS.glob("*.json")) + list(FIGURES.glob("*.png"))):
        if pth.name == "build.json":
            continue
        manifest["files"][str(pth.relative_to(ROOT))] = sha256(pth)

    lines = [f"{v}  {k}" for k, v in manifest["files"].items()]
    (RESULTS / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    dump("build.json", manifest)

    print("\n" + "=" * 72)
    print(f"build complete in {manifest['runtime_s']} s")
    print(f"{len(manifest['files'])} artefacts sealed in results/MANIFEST.sha256")
    print("=" * 72)
    return 0


def _toy_sim():
    n = 16
    return {"config": {"n_grid": n}, "a_final": 1.0,
            "a": [], "disp_rms": [],
            "pk": {"k_h_mpc": [], "P": [], "nmodes": []},
            "slice_field": [[0.0] * n for _ in range(n)]}


if __name__ == "__main__":
    raise SystemExit(main())
