# Limitations and defect register

This file is the machine-readable companion to Section 9 of the manuscript
(`paper/paper.tex`). Every defect below was found by us, on our own code, and
is listed whether or not it has been fixed.

The manuscript quotes these identifiers (D1 … D11). Do not renumber them
without updating `paper/paper.tex` and `paper/paper.pdf`. **D12** was added
after publication by a post-v1.0.2 audit; it is not quoted in the manuscript
and does not change any number in it.

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

### D12 — The comoving distance was the age integral

`background.comoving_distance_over_cH0()` integrated `1/(a E(a))`, which is
`H0*t0` — the *age* — not the line-of-sight comoving distance, which is
`∫ da/(a² E(a))`. The two functions therefore returned the same number for
every universe (`0.95035…`), and the comoving-distance column of
`results/background.json` was a duplicate of the age column.

The defect was invisible to the existing tests because none of them
exercised `comoving_distance_over_cH0`; the manuscript quotes the correct
integral (Eq. 3), so the code disagreed with its own paper. `H0*t0` and the
`gamma` result — the actual content of Verse I — were never affected.

**Fix applied.** The kernel is now `1/(a² E)` and `luminosity_distance`
uses the same kernel with the `(1+z)` prefactor. Two regression tests were
added (`test_comoving_distance_is_not_the_age_integral`,
`test_luminosity_distance_identity`). ΛCDM now returns `D_C/(c/H0) =
3.178966`, versus `H0*t0 = 0.950351`.

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

**Fix applied (f(R) only).** `src/ratiss/fr.py` now implements the
Hu–Sawicki quasi-static chameleon coupling

```
μ(k, a) = 1 + (1/3) x/(1+x),     x = (k/(a M(a)))²,
```

with `M(a)` the scalaron mass and `|f_R0|` the single free parameter. The
coupling is exercised in a **dedicated sector** (`models.fr_sector`,
`results/fr.json`, `figures/fig08_fr.png`), exactly as the granular cut-off
is: the sealed registry is left untouched (`μ = 1`), so no number in the
manuscript moves. The sector reports a genuine, *scale-dependent*
prediction — `μ → 1` at large scales, `μ → 4/3` below the Compton
wavelength, and `Δγ ≈ −0.055` at `k = 0.1 h/Mpc` for `|f_R0| = 1e-5`
(Compton wavelength 7.7 Mpc/h) — inside the reach of DESI and Euclid, while
`|f_R0| = 1e-6` is screened back onto ΛCDM. This is the repository's first
non-degenerate, *testable* modified-gravity signature. The `f(T)` row
remains `μ = 1` and is still degenerate by construction.

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
| N3 | The `f(R)` and `f(T)` registry rows reproduce ΛCDM by construction; the sealed registry is unchanged. The genuine, scale-dependent f(R) signature now lives in a dedicated sector (see D6). |
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
- Fourteen artefacts, fourteen SHA-256 hashes, `verify.py` exit code 0.
- Twenty-five unit tests, twenty-five passing.
- `D_C/(c/H0) = 3.178966` for ΛCDM, against `H0*t0 = 0.950351` — distinct
  quantities, checked against an independent quadrature (post-fix, D12).
- f(R) chameleon limits `μ → 1` (large scales) and `μ → 4/3` (small scales),
  with growth enhanced relative to ΛCDM and screened monotonically as
  `|f_R0|` decreases (post-fix, D6).
