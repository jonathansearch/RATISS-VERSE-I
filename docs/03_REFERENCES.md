# 03 — References

Every source used to fix a number, a formula or a bound. Where a fitting
formula is used, it is named as such and its accuracy is stated.

---

## Cosmological parameters

1. **Planck Collaboration** (2020). *Planck 2018 results. VI. Cosmological
   parameters.* A&A **641**, A6. arXiv:1807.06209.
   → the frozen parameter set in `src/ratiss/constants.py`.

2. **DESI Collaboration** (2025). *DESI DR2 results. II. Measurements of
   baryon acoustic oscillations and cosmological constraints.* Phys. Rev. D
   **112**, 083515. arXiv:2503.14738.
   → the current state of the art on the expansion history and the
   `w₀w_a` preference (3.1σ for evolving dark energy).

3. **Riess et al.** (2022). *A comprehensive measurement of the local value of
   the Hubble constant.* ApJL **934**, L7. arXiv:2112.04510.
   → the `H₀` tension context.

4. **Euclid Collaboration** (2020). *Euclid preparation. I. The Euclid Wide
   Survey.* A&A **642**, A191.
   → the survey whose `γ` and `fσ₈` measurements would test this work.

---

## Modified gravity and dark energy models

5. **Hu, W. & Sawicki, I.** (2007). *Models of f(R) cosmic acceleration that
   evade solar system tests.* Phys. Rev. D **76**, 064004. arXiv:0705.1158.
   → `f(R)` Hu–Sawicki; designer background identical to ΛCDM.

6. **Dvali, G., Gabadadze, G. & Porrati, M.** (2000). *4D gravity on a brane
   in 5D Minkowski spacetime.* Phys. Lett. B **485**, 208. arXiv:hep-th/0005016.
   → DGP self-accelerating branch.

7. **Linder, E. V.** (2010). *Einstein's other gravity and the acceleration of
   the Universe.* Phys. Rev. D **81**, 127301. arXiv:1005.3039.
   → `f(T)` teleparallel gravity.

8. **Solà Peracaula, J., Gómez-Valent, A. & de Cruz Pérez, J.** (2019).
   *Vacuum variable ΛCDM.* MNRAS **490**, 1. arXiv:1906.01156.
   → running vacuum model, `Λ(H) = c₀ + 3νH²`.

9. **Freese, K. & Lewis, M.** (2002). *Cardassian expansion.* ApJ **606**, 579.
   → modified Friedmann form used here.

10. **Amendola, L. & Tsujikawa, S.** (2010). *Dark Energy: Theory and
    Observations.* Cambridge University Press.
    → the general `μ(k,a)` formalism for dark energy and modified gravity.

11. **Weinberg, S.** (1989). *The cosmological constant problem.* Rev. Mod.
    Phys. **61**, 1.
    → context for why a small constant term is hard, and hence why
    alternatives to Λ are explored.

---

## Structure formation

12. **Peebles, P. J. E.** (1980). *The Large-Scale Structure of the Universe.*
    Princeton University Press.
    → the growth equation and the `f ≈ Ω_m^{3/5}` identity tested in `tests/`.

13. **Bardeen, J. M., Bond, J. R., Kaiser, N. & Szalay, A. S.** (1986).
    *Galaxy formation and large-scale structure.* ApJ **304**, 15.
    → BBKS transfer function used for initial conditions. **A fitting formula,
    not a Boltzmann solver**; accuracy of order 5–10 % in `P(k)`, which is why
    the repository does not claim percent-level `P(k)`.

14. **Zel'dovich, Ya. B.** (1970). *Gravitational instability: an approximate
    theory for large density perturbations.* A&A **5**, 84.
    → the initial-condition scheme in `pm.zeldovich_ics`.

15. **Press, W. H. & Schechter, P.** (1974). *Formation of galaxies and
    clusters of galaxies by self-similar gravitational condensation.*
    ApJ **187**, 425.
    → background for the "smallest structure" argument.

16. **Navarro, J. F., Frenk, C. S. & White, S. D. M.** (1997). *A universal
    density profile from hierarchical clustering.* ApJ **490**, 493.
    → the smallest bound dark matter halos discussed as a constraint source.

---

## Numerical methods

17. **Hockney, R. W. & Eastwood, J. W.** (1988). *Computer Simulation Using
    Particles.* Adam Hilger.
    → the particle-mesh scheme and the superconformal time leapfrog used in
    `src/ratiss/pm.py`.

18. **Springel, V.** (2005). *The cosmological simulation code GADGET-2.*
    Nature **435**, 629. arXiv:astro-ph/0503200.
    → the reference for the "expansion imposed, perturbations simulated"
    architecture that Theorem 1 forces.

19. **Angulo, R. E. & Hahn, O.** (2022). *Large-scale dark matter
    simulations.* ApJS **255**, 17.
    → standard reference for resolution and box-size trade-offs; used to set
    the honest accuracy statement for the demonstration-grade runs here.

---

## Note on scope

This reference list is deliberately conservative: no source is cited that was
not actually used to fix a number, a formula or a bound in the code. Where a
result is computed here and not taken from the literature (Theorem 1–3, the
`r_c` bound, the degeneracy table), it is marked as computed and reproducible
from `results/`, not attributed to anyone.
