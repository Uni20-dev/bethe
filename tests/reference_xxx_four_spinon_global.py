# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent finite-cubature global moments; no thermodynamic claim.

Two Gauss nodes on each gap interval (0,2). Locate P=pi by scalar root solving,
evaluate the independent 8F7 physical density, and integrate both common-shift
moments by direct positive mpmath quadrature. No C++ geometry or A0/A1 formula
is used. The output includes the ordered 24 and both observable 1/(2*pi)
measures. The fixed shift padding of 60 limits this reference to at most 45 digits.
mpmath is only needed to regenerate this small reference.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import reference_xxx_four_spinon as reference
from reference_xxx_transition_rate import transition_rate


def calculate(dps):
    with mp.workdps(dps):
        reference.rate = transition_rate
        nodes = [1-1/mp.sqrt(3), 1+1/mp.sqrt(3)]
        total = [mp.mpf(0), mp.mpf(0)]
        for index, gaps in enumerate(itertools.product(nodes, repeat=3)):
            a, b, c = gaps
            delta = [mp.mpf(0), a, a+b, a+b+c]

            def phase(eta):
                return 2*sum(mp.atan(mp.exp(eta+d)) for d in delta)

            eta = mp.findroot(lambda x: phase(x)-mp.pi,
                              (-delta[-1]-4, 4), solver="bisect",
                              maxsteps=4*mp.mp.prec)
            beta = [eta+d for d in delta]
            shape = reference.density(beta)*mp.fprod(mp.sech(x) for x in beta)
            cuts = [-delta[-1]-60, -delta[-1], -delta[2], -delta[1], 0, 60]
            for moment in range(2):
                def integrand(shift):
                    sine = mp.sin(phase(shift)/2)
                    energy = mp.pi/2*sum(mp.sech(shift+d) for d in delta)
                    return sine*sine*(energy if moment else 1)
                total[moment] += shape*mp.quad(integrand, cuts)*24/(2*mp.pi)**2
            print(json.dumps({"node": index, "partial": [mp.nstr(x,dps-5) for x in total]}), flush=True)
        return {"dps": dps, "order": 2, "gap_cutoffs": [2,2,2],
                "moments": [mp.nstr(x,dps-5) for x in total]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=40)
    args = parser.parse_args()
    if args.dps < 20 or args.dps > 45:
        parser.error("use 20..45 digits for this finite-shift reference")
    print(json.dumps(calculate(args.dps), indent=2), flush=True)
