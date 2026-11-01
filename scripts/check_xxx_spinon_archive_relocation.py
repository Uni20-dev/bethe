#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Integration check for the declared inputs of retained spinon benchmarks.

Requires the raw archive and reference; no numerical recalculation is needed.
Copies are temporary and placed below --scratch-directory (host-local storage).
"""

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store_directory", type=Path)
    parser.add_argument("--pilot-reference", type=Path, required=True)
    parser.add_argument("--scratch-directory", type=Path, required=True)
    args = parser.parse_args()
    analyzer = ROOT / "scripts/check_xxx_six_spinon_precision_probe.py"
    expected = (ROOT / "docs/validation/six-spinon-precision-floor-20261101.json").read_bytes()

    def run(store, reference):
        return subprocess.check_output(
            [sys.executable, str(analyzer), str(store), "--pilot-reference", str(reference)], timeout=60
        )

    missing = subprocess.run([sys.executable, str(analyzer), str(args.store_directory)],
                             capture_output=True, timeout=10)
    assert missing.returncode == 2 and b"--pilot-reference" in missing.stderr
    assert run(args.store_directory, args.pilot_reference) == expected
    with tempfile.TemporaryDirectory(dir=args.scratch_directory, prefix="spinon-relocation-") as temporary:
        root = Path(temporary)
        store = root / "observations" / args.store_directory.name
        reference = root / "independent-reference" / args.pilot_reference.name
        reference.parent.mkdir(parents=True)
        shutil.copytree(args.store_directory, store)
        shutil.copy2(args.pilot_reference, reference)
        # Neither input has any of the original archive's sibling directories.
        assert run(store, reference) == expected
    print("Required reference argument checked; original and relocated inputs reproduce the committed summary exactly.")
    print("Summary SHA256:", hashlib.sha256(expected).hexdigest())


if __name__ == "__main__":
    main()
