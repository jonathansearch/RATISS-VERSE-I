# Limitations and defect register

This file is the machine-readable companion to Section 9 of the manuscript
(`paper/paper.tex`). Every defect below was found by us, on our own code, and
is listed whether or not it has been fixed.

The manuscript quotes these identifiers (D1 … D11). Do not renumber them
without updating `paper/paper.tex` and `paper/paper.pdf`.

---

## Fixed defects

### D1 — Radiation was handled inconsistently across background families

`E_lcdm()` omitted radiation while `make_wcdm()` retained it.

| convention | `H0*t0` |
|---|---|
| without radiation (`E_lcdm`, original) | 0.9507304260038759 |
| with radiation (`make_wcdm`)            | 0.950350929653681  |
| difference                              | **3.795e-04**      |

The difference between ΛCDM and the `w = -1.004` model — the very signal the
registry was built to resolve — is 5.64e-04. The convention inconsistency was
therefore 67 % of the signal.

**Fix applied.** Radiation (`OMEGA_R = 9.15e-05`) is now present in every
background of the registry. All numbers in the manuscript are post-fix.
Pre-fix values are preserved in `git` history (commit `6f214ce`) for audit.

### D4 — A magic number in the classifier

`models._status()` compared the growth index against the literal `0.5544`,
which is computed nowhere in the code.

**Fix applied.** `models.gamma_reference()` now derives it from the growth
solver on first use: `gamma = 0.5543609664573533`. Difference from the
literal: 4e-05.

---

## Retained defects

### D2 — A quadrature error presented as a model uncertainty

`H0t0_err` in `results/predictions.json` is of order 1e-14. It is the
discrepancy between two quadrature tolerances — the accuracy of an integral,
not an uncertainty on a cosmological model. No statement in this repository
carries a 1e-14 uncertainty. Keep the field (it proves convergence), never
read it as an error bar.

### D3 — Arbitrary classification thresholds

A universe is labelled `EXCLUDED` when `|ΔH0t0| > 0.02` and `TESTABLE` when
`|Δγ| > 0.03`. Both numbers are chosen. The quantity they partition spans
four orders of magnitude (5.6e-04 to 1.6e-01). The labels are an ordering
device, not a measurement. Not replaced, because replacing an arbitrary
threshold with a differently arbitrary one is not an improvement.

### D5 — The growth initial condition is a pure matter-era mode

Equation (5) is integrated from `N = -9` with `D = D' = a`. This discards
(i) the radiation-era modification of the growing mode, (ii) any early
modified-gravity transient, (iii) the effect of `μ ≠ 1` before `N = -9`.
Expected effect on `γ`: order 1e-03. Not quantified.

### D6 — `f(R)` and `f(T)` rows are degenerate by construction

Both use `μ = 1` and `E = E_ΛCDM`. Their agreement with ΛCDM is an identity,
not a result. A genuine `f(R)` model needs chameleon screening and a
scale-dependent `μ(k, a)`, neither of which is implemented.

### D7 — `μ(x)` is an interpolation, not a transform

The exact Fourier transform of a `1/r^4` force does not exist (logarithmic
infrared divergence). `μ(x) = (1-x²)/(1+x²)` is the minimal bounded
regularisation reproducing the three exact limits. The assignment `x = k·r_c`
is a further modelling choice: a factor of order unity in the mapping would
shift the 5.8-decade gap by a fraction of a decade, and would not close it.

### D8 — "EXCLUDED" does not mean "contradicts data"

The classification compares each universe against the ΛCDM value of `H0*t0`
at the same `Ω_m`. It does not compare against a measurement, because
`H0*t0` is not directly measured. The one genuine observational exclusion in
this work — `w = -0.302`, falsified by stellar ages at 11.51 Gyr — is made
separately, in physical units.

### D9 — The simulations are not cosmological

`32^3 = 32 768` particles in a `(200 Mpc/h)^3` box. Mean interparticle
separation 6.3 Mpc/h; grid Nyquist `k ≈ 1.0 h/Mpc`; 6–12 modes in the lowest
power-spectrum bins. The runs establish that leapfrog, CIC deposition and the
FFT Poisson solver operate and conserve. They are not a measurement of the
matter power spectrum. A production run needs ≥ `512^3` particles over a Gpc.

### D10 — A stale label

`w_m1004` was constructed by inverting `H0*t0 = 0.9549` under the pre-fix
(radiation-free) convention. After the D1 fix, `w = -1.004049` gives
`H0*t0 = 0.950915`, and the value that actually reproduces 0.9549 with matter
present is `w = -1.0332263`. The computed values in that row are correct; the
key name is stale by 0.029 in `w`.

### D11 — The transfer function is a fit

Initial conditions use BBKS (Bardeen et al. 1986) with shape parameter
`Γ = Ω_m h`. Accurate to a few percent against a Boltzmann code in the range
of interest. Contains no baryon acoustic feature, no neutrino suppression, no
radiation-era physics.

---

## Negative results

| ID | Result |
|---|---|
| N1 | `H0*t0` is not a discriminator: six universes agree to 0.059 % |
| N2 | The target `H0*t0 = 0.9549` is an artefact of emptiness, not of any equation of state |
| N3 | `f(R)` and `f(T)` rows reproduce ΛCDM by construction; they constitute no evidence |
| N4 | The granular registry row is evaluated at `μ = 1`, hence vacuous; the cut-off lives only in the dedicated sector |
| N5 | The PM runs do not resolve the nonlinear regime |
| N6 | `fσ8` is normalised, not predicted: `σ8` is an input |

---

## What is *not* wrong

For balance, the following were checked and hold:

- `H0*t0 = 2/n` for a single fluid, exact, verified to 3e-16.
- ΛCDM age `13.7952 Gyr` against the Planck 2018 value `13.801 ± 0.024 Gyr`;
  agreement `-0.006 Gyr` from an independent quadrature.
- The Jeans swindle: the `k = 0` force mode vanishes to machine precision.
- `μ(0) = 1`, `μ(1) = 0`, `μ(∞) = -1`: the three limits hold exactly.
- CIC mass conservation to 1e-12.
- Twelve artefacts, twelve SHA-256 hashes, `verify.py` exit code 0.
- Twenty unit tests, twenty passing.
