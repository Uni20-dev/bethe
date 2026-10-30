# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional higher-precision W4 oracle for the recovery precision-growth test.

Uses the independent 8F7 contour sum, not the C++ rational-tail recurrence.
The pair factor uses mpmath Barnes G here: the separate physical-density tests
already check its independent oscillatory-integral representation. The older
40-digit quadrature fixture is not an 80-digit reference.
"""
from pathlib import Path
import sys

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import reference_xxx_four_spinon as reference
from reference_xxx_transition_rate import transition_rate


def main():
    reference.rate = transition_rate
    values = []
    for dps in (60, 85):
        with mp.workdps(dps):
            beta = [2 * mp.pi * mp.mpf(s) for s in ("-.6", "-.1", ".3", ".8")]
            value = reference.density(beta)
            values.append(value)
            print(dps, mp.nstr(value, dps - 5), flush=True)
    with mp.workdps(85):
        if abs(values[0] / values[1] - 1) >= mp.mpf("1e-55"):
            raise ArithmeticError("independent contour reference did not converge")


if __name__ == "__main__":
    main()
