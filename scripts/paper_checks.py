"""Recompute every derived number quoted in the manuscript.

Run:  PYTHONPATH=src python3 scripts/paper_checks.py
Nothing here is typed by hand; every value printed is a consequence of the
registry in src/ratiss/background.py.
"""
from __future__ import annotations

import json
import numpy as np
from scipy.optimize import brentq

from ratiss import background as bg, growth as gw, models as m, constants as C

OM = C.OMEGA_M


def main() -> None:
    out = {}

    ht_ref, err_ref = bg.hubble_time_times_h0(bg.E_lcdm)
    out["H0t0_lcdm"] = ht_ref
    out["H0t0_quad_err"] = err_ref

    # --- analytic identity: a single fluid rho ~ a^-n has H0t0 = 2/n ---------
    n = 2.0945
    out["H0t0_single_fluid_2_over_n"] = 2.0 / n
    out["H0t0_single_fluid_numeric"] = bg.hubble_time_times_h0(
        bg.make_single_fluid(n))[0]

    # --- what w reproduces the historical target 0.9549 with matter present --
    def f(w):
        return bg.hubble_time_times_h0(bg.make_wcdm(w, OM))[0] - 0.9549
    out["w_for_0p9549_with_matter"] = brentq(f, -1.6, -1.0, xtol=1e-10)

    # --- and what Omega_m would do it in Lambda-CDM -------------------------
    def g(om):
        return bg.hubble_time_times_h0(lambda a: bg.E_lcdm(a, om))[0] - 0.9549
    out["Omega_m_for_0p9549"] = brentq(g, 0.15, 0.3153, xtol=1e-10)

    # --- the registry -------------------------------------------------------
    rows = []
    for p in m.build_registry():
        rows.append({k: getattr(p, k) for k in
                     ("key", "label", "H0t0", "H0t0_err", "gamma",
                      "f_z0", "fsigma8_z1", "status")})
    out["registry"] = rows

    DEG = ("lcdm", "w_m1004", "rvm_001", "fr_hs", "ft", "granular")
    ht = [r["H0t0"] for r in rows if r["key"] in DEG]
    out["spread_degenerate"] = max(ht) - min(ht)
    out["spread_degenerate_rel"] = (max(ht) - min(ht)) / ht_ref
    out["degenerate_group"] = list(DEG)

    # --- physical clock: age of the universe in Gyr -------------------------
    hub = C.hubble_time_gyr()
    out["hubble_time_gyr"] = hub
    out["age_gyr"] = {r["key"]: r["H0t0"] * hub for r in rows}
    out["age_planck_2018_gyr"] = 13.801
    out["age_oldest_star_gyr"] = 13.5          # HD 140283 lower bound, ~14.0
    out["age_excess_lcdm_over_planck_gyr"] = (
        out["age_gyr"]["lcdm"] - out["age_planck_2018_gyr"])

    gam = [r["gamma"] for r in rows]
    out["gamma_min"] = min(gam)
    out["gamma_max"] = max(gam)
    out["gamma_spread"] = max(gam) - min(gam)
    out["gamma_spread_rel"] = (max(gam) - min(gam)) / m.gamma_reference()

    # --- granular sector ----------------------------------------------------
    gr = json.load(open("results/granular.json"))
    out["granular_gap_dex_dwarfs"] = gr["observable_window"]["gap_dex_dwarfs"]
    out["granular_gap_dex_microhalos"] = gr["observable_window"]["gap_dex_microhalos"]
    out["rc_needed_mpc"] = gr["observable_window"]["r_c_needed_for_observable_gamma_mpc"]
    out["rc_bound_dwarfs_mpc"] = gr["observable_window"]["r_c_bound_from_dwarfs_mpc"]

    # --- growth: ratio D_granular / D_lcdm over the scanned x ---------------
    x = np.array(gr["sector"]["x"]); Dr = np.array(gr["sector"]["D_rel"])
    out["granular_Drel_min"] = float(Dr.min())
    out["granular_Drel_at_x0p1"] = float(np.interp(0.1, x, Dr))
    out["granular_Drel_at_x1"] = float(np.interp(1.0, x, Dr))

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
