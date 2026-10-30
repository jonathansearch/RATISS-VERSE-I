"""
RATISS — Reference Atlas of Testable Inhomogeneous Structure Simulations.

A reference implementation of a family of cosmological universes, each
defined by three ingredients (background expansion, linear growth, local
gravitational coupling) and each propagated through the same solvers to
observable predictions.

The package is intentionally small and dependency-light: numpy, scipy and
matplotlib only. Every number in results/ is recomputed by scripts/build.py
and sealed by a SHA-256 manifest, so the repository is falsifiable on two
levels at once: physically (every prediction states what would kill it)
and computationally (anyone can rerun and compare hashes).
"""

__version__ = "1.0.0"
__all__ = ["background", "growth", "granular", "pm", "models", "figures", "constants"]
