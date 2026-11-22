#!/usr/bin/env python3
"""Unified end-to-end verification runner for MitoFitness Atlas (Phase 0 -> Stage 4)."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYTHON = ROOT / ".venv" / "bin" / "python"

TEST_SUITES = [
    "tests/test_phase0.py",
    "tests/test_stage1.py",
    "tests/test_stage2.py",
    "tests/test_stage3.py",
    "tests/test_stage4.py",
    "tests/test_stage5_proofs.py",
]


def main() -> int:
    exe = str(PYTHON) if PYTHON.exists() else sys.executable
    print("=" * 72)
    print("MitoFitness Atlas — Unified Verification Suite (Phase 0 -> Stage 4)")
    print("=" * 72)
    for rel in TEST_SUITES:
        path = ROOT / rel
        print(f"\n>>> Running {rel} ...")
        proc = subprocess.run([exe, str(path)], cwd=str(ROOT))
        if proc.returncode != 0:
            print(f"\nFAIL: {rel} exited with code {proc.returncode}")
            return proc.returncode
    print("\n" + "=" * 72)
    print("ALL 6 VERIFICATION SUITES PASSED (Phase 0, Stage 1, Stage 2, Stage 3, Stage 4, Stage 5 Proofs)")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
