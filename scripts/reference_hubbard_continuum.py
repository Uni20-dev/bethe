#!/usr/bin/env python3
"""Independent double-precision oracle for half-filled Hubbard charge continua.

Developer tool only: SciPy Fourier-Bessel quadrature and bounded optimization,
not the native library's nonoscillatory dispersion evaluation. No production
dependency. U>=1 avoids weak-coupling cancellation in this reference method.
Grid refinement is a diagnostic, not a proof of global optimality.
"""
import argparse
from functools import lru_cache
import json
import math
import warnings

from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import j0, j1


def wrap(p):
    return (p + math.pi) % (2*math.pi) - math.pi


class Lines:
    def __init__(self, interaction):
        self.u = interaction / 4
        self.cutoff = 60 / self.u
        self.calls = 0

    def integral(self, function):
        self.calls += 1
        value, error = quad(function, 0, self.cutoff, epsabs=2e-12, epsrel=2e-12, limit=500)
        if error > 1e-10:
            raise ArithmeticError(f"unresolved reference quadrature: {error}")
        return value

    def weight(self, w, spin):
        z = math.exp(-self.u*w)
        return 2*z/(1+z*z) if spin else 2*z*z/(1+z*z)

    def spin_momentum(self, rapidity):
        return math.pi/2-self.integral(lambda w: j0(w)*math.sin(w*rapidity)/w*self.weight(w, True))

    def charge_momentum(self, k):
        return math.pi/2-k-self.integral(lambda w: j0(w)*math.sin(w*math.sin(k))/w*self.weight(w, False))

    @lru_cache(maxsize=None)
    def energy(self, branch, p):
        if branch == "spinon":
            p = min(p, math.pi-p)
            if p == 0:
                return 0.0
            bound = 1+2*self.u/math.pi*math.log(4/p)
            rapidity = brentq(lambda x: self.spin_momentum(x)-p, 0, bound, xtol=1e-12)
            return 2*self.integral(lambda w: j1(w)*math.cos(w*rapidity)/w*self.weight(w, True))
        target = wrap(p + (math.pi if branch == "antiholon" else 0))
        if target < -math.pi/2:
            target += 2*math.pi
        k = brentq(lambda x: self.charge_momentum(x)-target, -math.pi, math.pi, xtol=1e-12)
        return 2*self.u+2*math.cos(k)+2*self.integral(
            lambda w: j1(w)*math.cos(w*math.sin(k))/w*self.weight(w, False))


def edges(lines, channel, total, intervals):
    first, second = channel.split("-")
    lo, hi = (0.0, math.pi) if first == "spinon" else (-math.pi, math.pi)

    def objective(p):
        return lines.energy(first, p)+lines.energy(second, wrap(total-p))

    xs = [lo+(hi-lo)*i/intervals for i in range(intervals+1)]
    ys = [objective(x) for x in xs]
    result = {}
    for label, sign in (("lower", 1), ("upper", -1)):
        candidates = [(ys[0], xs[0]), (ys[-1], xs[-1])]
        brackets = [(xs[0], xs[1]), (xs[-2], xs[-1])]
        for i in range(1, intervals):
            if sign*ys[i] <= sign*ys[i-1] and sign*ys[i] <= sign*ys[i+1]:
                brackets.append((xs[i-1], xs[i+1]))
        for bounds in brackets:
            opt = minimize_scalar(lambda x: sign*objective(x), bounds=bounds,
                                  method="bounded", options={"xatol": 1e-9, "maxiter": 200})
            if not opt.success:
                raise ArithmeticError("reference extremum refinement failed")
            candidates.append((objective(opt.x), float(opt.x)))
        energy, p = min(candidates, key=lambda pair: sign*pair[0])
        result[label] = dict(energy=energy, p1=p, p2=wrap(total-p), candidates=len(candidates))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--u", type=float, default=4)
    parser.add_argument("--momenta-over-pi", type=float, nargs="+", default=[0, 0.25, 0.5, 1])
    parser.add_argument("--intervals", type=int, nargs="+", default=[32, 64])
    parser.add_argument("--channels", nargs="+", choices=["spinon-holon", "spinon-antiholon", "holon-antiholon"],
                        default=["spinon-holon", "holon-antiholon"])
    args = parser.parse_args()
    if not math.isfinite(args.u) or args.u < 1 or any(n < 8 for n in args.intervals):
        parser.error("require finite U>=1 and at least 8 intervals")
    if any(not math.isfinite(p) or abs(p) > 1 for p in args.momenta_over_pi):
        parser.error("momenta/pi must be finite and in [-1,1]")
    warnings.simplefilter("error", IntegrationWarning)
    lines = Lines(args.u)
    results = []
    for channel in args.channels:
        for p in args.momenta_over_pi:
            scans = [dict(intervals=n, **edges(lines, channel, p*math.pi, n)) for n in args.intervals]
            results.append(dict(channel=channel, momentum_over_pi=p, scans=scans))
    print(json.dumps(dict(interaction=args.u, convention="symmetric Hamiltonian; t=1; half filling; zero field",
                          method="independent double-precision Fourier-Bessel quadrature; not certified",
                          quadratures=lines.calls, results=results), indent=2))


if __name__ == "__main__":
    main()
