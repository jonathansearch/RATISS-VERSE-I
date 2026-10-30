#!/usr/bin/env python3
"""
Verify the sealed build.

Recomputes the SHA-256 of every artefact listed in results/MANIFEST.sha256 and
compares. Exit code 0 means: this repository is bit-identical to the build
that produced the reported numbers.

    python scripts/verify.py
    python scripts/verify.py --rebuild   # rebuild first, then verify
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "MANIFEST.sha256"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    args = ap.parse_args()

    if args.rebuild:
        subprocess.run([sys.executable, str(ROOT / "scripts" / "build.py")], check=True)

    if not MANIFEST.exists():
        print("no manifest found: run scripts/build.py first")
        return 2

    ok = bad = missing = 0
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(maxsplit=1)
        path = ROOT / name.strip()
        if not path.exists():
            print(f"MISS  {name}")
            missing += 1
            continue
        if sha256(path) == digest:
            print(f"OK    {name}")
            ok += 1
        else:
            print(f"FAIL  {name}")
            bad += 1

    print("-" * 56)
    print(f"{ok} verified, {bad} mismatched, {missing} missing")
    return 0 if (bad == 0 and missing == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
