# 04 — Falsifiability

For every universe in the registry: **what observation would kill it**, and
**which instrument can make that observation**.

*Applies to RATISS Verse I. Any universe added in a later verse inherits the
same obligation: a registry entry without a kill condition is not admitted.*

A model that cannot be killed is not a scientific model. This document is
therefore the operative one.

---

## Summary matrix

| Universe | `H₀t₀` | `γ` | Status | Killed by |
|---|---|---|---|---|
| ΛCDM | 0.95073 | 0.5543 | reference | — |
| `f(R)` Hu–Sawicki | 0.95073 | 0.5543 | DEGENERATE | `γ` at `≲0.01` precision (background identical) |
| `f(T)` teleparallel | 0.95073 | 0.5543 | DEGENERATE | same |
| RVM `ν = 0.001` | 0.95046 | 0.5544 | DEGENERATE | `ν` from `H(z)` curvature, future |
| `w = −1.004049` | 0.95091 | 0.5542 | DEGENERATE | `w` at 0.1 % precision; already allowed |
| RATISS granular | 0.95073 | 0.5543 | **ALIVE-BUT-INVISIBLE** | a structure below `r_c`; see §2 |
| DGP self-accelerating | 0.89525 | 0.5663 | **EXCLUDED** | already excluded by `H₀t₀` at 5.8 % |
| Cardassian `m = 0.5` | 0.92503 | 0.5601 | **EXCLUDED** | already excluded at 2.7 % |
| `w = −0.30185` | 0.79267 | 0.5916 | **EXCLUDED** | already excluded at 17 % |
| single fluid `n = 2.0945` | 0.95488 | 0.5706 | **EXCLUDED** | contains no matter; no galaxies |

Status is computed in `src/ratiss/models.py::_status`, not assigned.

---

## 1. The degenerate set: what would break it

`ΛCDM`, `f(R)`, `f(T)` and RATISS-granular share `E(a)` **exactly**, because
`f(R)` and `f(T)` are constructed with a designer background. No measurement
of `H(z)`, `D_A(z)` or `H₀t₀`, at any precision, can separate them.

The separation lives in the **response to perturbations**. Three observables
carry it:

| Observable | What it measures | Sensitivity needed | Instrument |
|---|---|---|---|
| `γ = ln f / ln Ω_m` | how fast growth runs | ~0.01 absolute | DESI, Euclid |
| `fσ₈(z)` | growth × amplitude, redshift-resolved | ~2 % | DESI Y3, Euclid |
| `μ(k)` scale dependence | whether gravity is scale-free | ~5 % in `P(k)` at `k > 0.1` | Euclid, SKA |

For the models in the registry the spread in `γ` is `0.554 → 0.592`, i.e.
**7 %** — an order of magnitude above the projected sensitivity. If the
universe is one of these, the next data release will see it.

---

## 2. The granular universe: the bound, and the gap

The cut-off implies a minimum length: below `r_c`, `μ → 0` and growth stops. A
structure of comoving size `L` corresponds to `k = 2π/L`, so requiring the
cut-off to stay below the observed structures gives

```
r_c  ≤  x_cut · L_min / (2π),     x_cut = 0.1
```

| Constraint source | `L_min` | `r_c` bound |
|---|---|---|
| galaxy clusters | 2 Mpc | 32 kpc |
| galaxies | 10 kpc | 160 pc |
| ultra-faint dwarfs *(observed)* | 100 pc | **1.6 pc** |
| dark matter microhalos *(predicted)* | 10⁻² pc | **33 AU** |

The binding constraint is the smallest thing that is known to exist.

### The gap

An **observable** deviation of `γ` needs `x = k r_c ≳ 0.3` at a scale the
surveys probe, `0.01 < k < 0.3 h/Mpc`. That requires

```
r_c ≈ 0.3 / k  ≈  1 – 30 Mpc
```

while the dwarf-galaxy bound requires `r_c < 1.6 pc`.

```
gap = log10(1 Mpc / 1.6 pc) = 5.8 decades
```

Against predicted microhalos the gap is **9.8 decades**.

**Verdict.** The granular universe is **not excluded**. It is **not
observable** either, at any scale that current or planned instruments can
reach. Its signature sits 6 to 10 orders of magnitude below the smallest
structure we can see.

This is a real result and it is the honest one. A model that is alive but
invisible has exactly one scientific status: it must wait for a probe of
sub-parsec structure, which does not exist and is not planned.

### What would kill it

* **A structure smaller than `r_c`.** If `r_c = 1.6 pc` and a `0.1 pc` dark
  matter subhalo is detected (via millilensing or stellar-stream heating), the
  model dies immediately.
* **An observed `γ` deviating from 0.554 at scales where `r_c` is already
  excluded.** If Euclid measures `γ = 0.59` at `k = 0.1 h/Mpc`, the granular
  explanation would require `r_c ≈ 3 Mpc`, which is already excluded by
  galaxies — so the granular model could not be the cause.
* **A detection of the repulsive core.** At `r < r_c` the force reverses sign.
  Any dynamical probe of the regime would be definitive.

---

## 3. What is *not* falsifiable, and why it is reported anyway

The single-fluid inversion `w = −0.30185` deserves a paragraph because it is
the most attractive-looking number in the parent project and it is not a
prediction.

`w = −0.302` follows from solving `2/(3(1+w)) = H₀t₀` **assuming the universe
contains one substance**. That assumption is false: matter exists at
`Ω_m = 0.3153 ≫ 0`. With matter present the same equation gives `H₀t₀ =
0.7927`, excluded by 17 %.

The number is the equation of state of an empty universe. It is reported in
`results/predictions.json` with the status `EXCLUDED`, and the derivation of
*why* it is excluded is in [`02_THEOREMS.md`](02_THEOREMS.md), Corollary.

---

## 4. Reproducibility as a second axis of falsification

Physical falsifiability says what observation kills a model. Computational
falsifiability says what *disagreement* kills a claim.

* `scripts/build.py` is deterministic: fixed parameters, fixed seeds, no
  network, no unseeded generator.
* `results/MANIFEST.sha256` seals all 12 artefacts.
* `scripts/verify.py` recomputes every hash; the CI workflow runs it on every
  push.
* `tests/test_physics.py` pins 20 numbers, including four analytic identities
  (`H₀t₀ = 2/n` for a single fluid, `H₀t₀ = 2/3` for Einstein–de Sitter,
  `D ∝ a` in matter domination, `f = Ω_m^{3/5}`).

If any of these fails, the number quoted in the README is wrong. That is a
stronger guarantee than peer review can provide for a numerical result.
