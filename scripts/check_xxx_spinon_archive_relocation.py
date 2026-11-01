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
    references = parser.add_mutually_exclusive_group(required=True)
    references.add_argument("--pilot-reference", type=Path, help="archived pilot file for the PR #93 analyzer")
    references.add_argument("--reference-directory", type=Path, help="frozen PR #93 archive for the stage-1 analyzer")
    parser.add_argument("--scratch-directory", type=Path, required=True)
    args = parser.parse_args()
    if args.pilot_reference is not None:
        analyzer = ROOT / "scripts/check_xxx_six_spinon_precision_probe.py"
        expected_path = ROOT / "docs/validation/six-spinon-precision-floor-20261101.json"
        reference_input, flag = args.pilot_reference, "--pilot-reference"
    else:
        analyzer = ROOT / "scripts/check_xxx_scaled_spinon_endpoints.py"
        expected_path = ROOT / "docs/validation/scaled-spinon-endpoints-20261101.json"
        reference_input, flag = args.reference_directory, "--reference-directory"
    expected = expected_path.read_bytes()

    def run(store, reference):
        return subprocess.check_output(
            [sys.executable, str(analyzer), str(store), flag, str(reference)], timeout=60
        )

    missing = subprocess.run([sys.executable, str(analyzer), str(args.store_directory)],
                             capture_output=True, timeout=10)
    assert missing.returncode == 2 and flag.encode() in missing.stderr
    assert run(args.store_directory, reference_input) == expected
    with tempfile.TemporaryDirectory(dir=args.scratch_directory, prefix="spinon-relocation-") as temporary:
        root = Path(temporary)
        store = root / "observations" / args.store_directory.name
        reference = root / "independent-reference" / reference_input.name
        reference.parent.mkdir(parents=True)
        shutil.copytree(args.store_directory, store)
        if args.reference_directory is None:
            shutil.copy2(reference_input, reference)
        else:
            reference.mkdir()
            for name in ("geometry-bits-512.jsonl", "geometry-receipt.json"):
                shutil.copy2(reference_input / name, reference / name)
        # Neither input has any of the original archive's sibling directories.
        assert run(store, reference) == expected
        if args.reference_directory is not None:
            # Relocation must not weaken validation of the frozen reference.
            with (reference / "geometry-bits-512.jsonl").open("ab") as output:
                output.write(b"\n")
            corrupted = subprocess.run([sys.executable, str(analyzer), str(store), flag, str(reference)],
                                       capture_output=True, timeout=10)
            assert corrupted.returncode != 0 and b"reference SHA256 mismatch" in corrupted.stderr
    print("Required reference argument checked; original and relocated inputs reproduce the committed summary exactly.")
    print("Summary SHA256:", hashlib.sha256(expected).hexdigest())


if __name__ == "__main__":
    main()
