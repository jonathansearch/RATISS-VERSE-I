# results/

Every file here is produced by `scripts/build.py` and sealed in
`MANIFEST.sha256`. Nothing is hand-edited.

| file | content |
|---|---|
| `background.json` | `E(a)`, `H₀t₀` with quadrature error, comoving distance, per universe |
| `growth.json` | `f(z=0)`, `γ`, `fσ₈(z)` for every redshift, per universe |
| `predictions.json` | the falsifiability matrix: status and reason, computed |
| `granular.json` | the `μ(x)` sector, growth suppression, `r_c` bounds, observable window |
| `simulations.json` | particle-mesh runs: configuration, growth history, `P(k)`, density slice |
| `build.json` | build metadata: versions, timestamps, runtime (not sealed: it seals itself) |
| `MANIFEST.sha256` | SHA-256 of every artefact |

Verify with:

```bash
python scripts/verify.py
```

## Provenance of the laboratory mark

`logo/ratiss-labs-mark.png` is the existing RATISS Labs logo, taken unchanged
from `jonathansearch/Ratiss-Jonathan-Labs-` (`docs/assets/logo.png`). Its
SHA-256 is recorded in `logo/CREDITS.md`. It was not redrawn or regenerated.
