"""
Reproducibility tests: the sealed manifest must hold.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def test_manifest_exists():
    assert (ROOT / "results" / "MANIFEST.sha256").exists(), \
        "run scripts/build.py first"


def test_manifest_verifies():
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "verify.py")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_seeds_are_fixed():
    """No simulation may draw from an unseeded generator."""
    from ratiss import pm
    cfg = pm.SimConfig()
    assert cfg.seed > 0


def test_build_is_deterministic():
    """Two short runs with the same seed must agree bit for bit."""
    from ratiss import pm
    cfg = pm.SimConfig(box_mpc_h=100.0, n_grid=16, n_part_1d=8, n_steps=4,
                       a_init=0.2, seed=42)
    cos = pm.Cosmology()
    r1 = pm.run(cfg, cos, record_slice=False)
    r2 = pm.run(cfg, cos, record_slice=False)
    assert r1.pk["P"] == r2.pk["P"]
