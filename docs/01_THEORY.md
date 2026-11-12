# 01 — Theory

Everything in RATISS follows from three functions. This document states them
precisely so that a reader can check the code against the physics without
reading the code.

---

## 1. Background expansion

Every universe is defined by `E(a) = H(a)/H0`, normalised so that `E(1) = 1`.

From `E(a)` every background observable follows by quadrature:

```
H0 t0            = ∫₀¹ da / (a E(a))
D_C / (c/H0)     = ∫₀¹ da / (a² E(a))
D_L(z)           = (1+z) D_C(z)
```

The age integral is the one that has dominated the parent project, and it is
the one that turns out to carry almost no information (see §5).

Implemented forms:

| key | `E²(a)` | family |
|---|---|---|
| `lcdm` | `Ω_m a⁻³ + Ω_R a⁻⁴ + Ω_Λ` | GR |
| `w_m0302` | `Ω_m a⁻³ + Ω_R a⁻⁴ + Ω_DE a^{−3(1+w)}`, `w = −0.30185` | fluid |
| `w_m1004` | same, `w = −1.004049` | fluid |
| `rvm_001` | `[Ω_m a⁻³ + (1−Ω_m−ν)] / (1−ν)`, `ν = 0.001` | running vacuum |
| `dgp` | `(√(Ω_m a⁻³ + Ω_rc) + √Ω_rc)²`, `Ω_rc = ((1−Ω_m)/2)²` | braneworld |
| `cardassian_05` | `Ω_m a⁻³ + (1−Ω_m) a^{−0.5}` | modified Friedmann |
| `fr_hs` | `= lcdm` (designer background) | `f(R)` |
| `ft` | `= lcdm` (designer background) | `f(T)` |
| `granular` | `= lcdm` (cut-off acts only in the growth sector) | RATISS |
| `fluid_n2094` | `a^{−n}`, `n = 2.0945`, nothing else | single fluid |

---

## 2. Linear growth

```
D'' + (2 + dlnH/dN) D' − (3/2) Ω_m(N) μ(k,a) D = 0
```

with `N = ln a`, `' = d/dN`, and

```
Ω_m(N) = Ω_m0 e^{−3N} / E²(N)
```

The **effective coupling** `μ(k,a)` is the ratio of the gravitational coupling
felt by perturbations to Newton's constant. Everything that departs from
General Relativity is expressed as a `μ`, which is why one solver covers every
universe in the registry.

Observables:

```
D(a)        growth factor
f(a)        = dlnD/dlna                    growth rate
γ           = ln f / ln Ω_m(a=1)           growth index
fσ8(z)      = f(z) · D(z)/D(0) · σ8
```

`solve_ivp` with LSODA, `rtol = 1e−11`, started at `a_i = e⁻⁹` in the matter
era where `D ∝ a` is exact to machine precision. Radiation is not included in
the growth equation; at `a_i = e⁻⁹` this is a `10⁻⁶`-level approximation to
the true initial condition and it is documented rather than hidden.

---

## 3. Local gravity: the granular coupling

The RATISS force law is

```
F(r) = −G/r² + K/r⁴ = −(G/r²)[1 − (r_c/r)²],     r_c = √(K/G)
```

Three regimes, exact:

| `r/r_c` | `F/F_Newton` | behaviour |
|---|---|---|
| ≫ 1 | → 1 | Newton |
| = 1 | **0** | force vanishes |
| < 1 | < 0 | repulsive |

The relative deficit is `1 − μ = (r_c/r)²`, i.e. a two-derivative correction.
Fourier transforming and regularising at the cut-off itself gives the bounded
interpolation

```
μ(x) = (1 − x²)/(1 + x²),        x = k r_c
```

which reproduces all three limits: `μ(0) = +1`, `μ(1) = 0`, `μ(∞) = −1`.

### Why a cut-off is not optional

The exact Fourier transform of a `1/r⁴` force is not a tempered distribution:
it diverges logarithmically in the infrared. A universe with a `K/r⁴` term and
no short-distance cut-off has no well-defined background. The cut-off at `r_c`
is therefore not a convenience, it is a consistency requirement — and that is
the strongest theoretical argument for taking `r_c` seriously as a physical
scale.

---

## 4. What the cut-off implies

`μ(x) < 1` for any `x > 0`, so growth is **suppressed** relative to GR, and at
`x = 1` it stops entirely. A gravitational cut-off means a **minimum length**:

> Below `r_c`, no structure can form.

That statement converts a free parameter into an observational bound, which is
the whole reason RATISS exists. The bound is derived in
[`04_FALSIFIABILITY.md`](04_FALSIFIABILITY.md).

---

## 5. The degeneracy

`H0 t0` is an integral of `1/E` over the whole history. It is a single number
and it is insensitive to the shape of `E(a)`. The registry demonstrates this
computationally:

| universe | `H0 t0` |
|---|---|
| ΛCDM | 0.95073 |
| `f(R)` Hu–Sawicki | 0.95073 |
| `f(T)` teleparallel | 0.95073 |
| RATISS granular | 0.95073 |
| RVM `ν = 0.001` | 0.95046 |
| `w = −1.004049` | 0.95091 |

Spread over five physically distinct theories: **4.5 × 10⁻⁴**, i.e. 0.05 %.

This is a *negative* result and it is the most useful one in the repository:
it says that the observable the parent project spent its effort on cannot, in
principle, select between theories. Any further work on `H0 t0` as a
discriminator is unfalsifiable by construction.

The growth index `γ` spreads over `0.554 → 0.592` for the same set, i.e. 7 % —
two orders of magnitude more discriminating, and inside the reach of
stage-IV surveys.

---

## 6. Numerical policy

* Every quadrature reports an independent error estimate (`_quad` runs twice
  at different tolerances and takes the difference).
* Every figure is rendered from an array computed at build time; no number is
  typed into a plotting call.
* Every random draw comes from a seeded `numpy.random.Generator`.
* The build is sealed with SHA-256 (`results/MANIFEST.sha256`).
