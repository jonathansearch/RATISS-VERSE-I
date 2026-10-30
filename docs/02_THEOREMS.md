# 02 — Theorems

Two results in this repository are not fits and not simulations. They are
statements that can be proved in a few lines and that constrain what any model
of this kind can do.

---

## Theorem 1 — The Jeans swindle (no expansion from particles in a torus)

**Statement.** In a periodic box with a homogeneous mean density, the mean
density exerts no net force on any particle. Consequently the expansion rate
`a(t)` cannot be recovered from the particle dynamics; it must be supplied.

**Proof.** Let `ρ = ρ̄ (1 + δ)` with `δ` periodic and `⟨δ⟩ = 0`. The
Newtonian potential solves

```
∇²Φ = 4πG ρ̄ (1 + δ)
```

Decompose `Φ = Φ̄ + φ` with `∇²Φ̄ = 4πG ρ̄`. For a periodic δ, the mode
`k = 0` of `φ` vanishes, so `φ` has zero mean and `∇φ` averages to zero. The
force from `ρ̄` is therefore a **uniform field**, i.e. zero force on a
relative coordinate. `a(t)` never appears. ∎

**Consequence.** Every N-body code in existence imposes `a(t)` from a
Friedmann solver and lets the particles produce only `δ`. RATISS does the
same, explicitly (`src/ratiss/pm.py`, module docstring). A model that claims
to *derive* the expansion from particle dynamics inside a periodic box is
making a category error — unless the mean field is treated non-periodically,
which reintroduces the boundary the box was meant to remove.

---

## Theorem 2 — No absolute value (scale invariance of the r_c model)

**Statement.** The granular model has exactly one dimensionless parameter,
`η = r_c/ε`, and its dimensionless outputs are invariant under

```
G → λ³ G,    K → λ⁵ K,    ε → λ ε,    r_c = √(K/G) → λ r_c
```

**Proof.** The force scales as `F(λr) = λ⁻² F(r)` under the transformation, so
every trajectory in units of `ε` is identical. Only ratios survive. ∎

**Consequence.** `r_c` cannot be predicted from inside the model. It can only
be *bounded* from outside, by observation. This is why
[`04_FALSIFIABILITY.md`](04_FALSIFIABILITY.md) is the operative document: it is
the only place where the model touches data.

*Verified numerically in `tests/test_physics.py` and, in the parent project,
by an explicit scan over `λ ∈ [10⁻¹, 10³]`.*

---

## Theorem 3 — `H₀t₀` is a one-parameter functional (degeneracy)

**Statement.** `H₀t₀ = ∫₀¹ da/(a E(a))` is a single scalar functional of
`E(a)`. Any two models whose `E(a)` differ by a function with zero weighted
integral return the same `H₀t₀`, however different their physics.

**Proof.** Immediate from the definition. ∎

**Consequence, and the reason this repository exists.** `f(R)`, `f(T)` and
`ΛCDM` are *constructed* to share `E(a)` exactly, so they share `H₀t₀` to
machine precision. No improvement in the measurement of `H₀t₀` can ever
separate them. The discriminating observable has to be sensitive to the
*response* to perturbations rather than the *background*, i.e. to `μ(k,a)`.

*Verified numerically:* the six universes in `results/predictions.json` marked
`DEGENERATE` span `4.5 × 10⁻⁴` in `H₀t₀`.

---

## Corollary — the single-fluid inversion is not a prediction

**Statement.** The value `w = −0.30185` follows from solving
`2/(3(1+w)) = H₀t₀` under the assumption that the universe contains **one**
substance. If matter is present at the observed `Ω_m`, the same equation gives

```
H₀t₀(w = −0.30185, Ω_m = 0.3153) = 0.7927
```

which is excluded by 17 %.

**Consequence.** `w = −0.302` is the equation of state of a universe with no
matter in it. It is a description of an empty universe, not a prediction for
this one. The correct inversion, with matter present, requires `w = −1.004`
(0.4 % phantom), which is consistent with Planck but is *not* a dynamical
dark energy — it is, within current sensitivity, Λ.

*Both numbers are computed in `results/predictions.json`.*
