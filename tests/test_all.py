"""pytest wrapper of the two test suites (run from the repository root:  python -m pytest -q tests)."""
import os
import subprocess
import sys

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code")


def _run(args):
    r = subprocess.run([sys.executable, "-X", "utf8"] + args, cwd=CODE, capture_output=True, text=True)
    print(r.stdout[-1500:])
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "FAIL" not in r.stdout


def test_original_model():
    _run(["oral_peptide_mc.py", "--test"])


def test_extension():
    _run(["run_extended.py", "--test"])
