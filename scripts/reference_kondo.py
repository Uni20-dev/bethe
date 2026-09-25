#!/usr/bin/env python3
"""Independent mpmath oracle for the planned universal T=0 Kondo evaluator.

Not a production solver or a build dependency. Equations (124)-(126) of
Barcza et al., arXiv:1911.08279. x=b/T_B with b the full Zeeman splitting
and T_B=2*T1 in that paper. Energy is [E_imp(b)-E_imp(0)]/T_B.
"""
import argparse
import json

import mpmath as mp


def coefficient(n):
    return mp.exp((n-mp.mpf("0.5"))*mp.log(n+mp.mpf("0.5"))
                  - n-mp.mpf("0.5")-mp.loggamma(n+1))/(2*mp.sqrt(mp.pi))


def low(x, energy=False):
    def term(n):
        value = (-1)**n * coefficient(n) * x**(2*n+1)
        return -value*x/(2*n+2) if energy else value
    return mp.nsum(term, [0, mp.inf], method="alternating")


def amplitude(w):
    if not w:
        return mp.sqrt(mp.pi)
    return mp.exp(mp.loggamma(w+mp.mpf("0.5")) - w*mp.log(w) + w)


def high(x, energy_at_one):
    log_x = mp.log(x)
    prefactor = 1/(2*mp.pi**mp.mpf("1.5"))

    def magnetization_integrand(w):
        if not w:
            return mp.pi*mp.sqrt(mp.pi)
        return mp.sinpi(w)/w * amplitude(w)*mp.exp(-2*w*log_x)

    def energy_integrand(w):
        z = 1-2*w
        integral = mp.expm1(z*log_x)/z if z else log_x
        if not w:
            return mp.pi*mp.sqrt(mp.pi)*integral
        return mp.sinpi(w)/w * amplitude(w)*integral

    # Resolve the endpoint's w*log(w) nonanalyticity separately. The tail
    # is oscillatory even for the integrated energy at x>1.
    def integrate(f):
        return mp.quad(f, [0, mp.mpf("0.5"), 1]) + mp.quadosc(f, [1, mp.inf], omega=mp.pi)

    return (mp.mpf("0.5")-prefactor*integrate(magnetization_integrand),
            energy_at_one-(x-1)/2+prefactor*integrate(energy_integrand))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=50)
    parser.add_argument("--points", nargs="+", default=["0.1", "0.5", "1", "2", "10"])
    args = parser.parse_args()
    if args.digits < 20:
        parser.error("use at least 20 digits")
    mp.mp.dps = args.digits + 15
    e1 = low(mp.mpf(1), energy=True)
    rows = []
    for token in args.points:
        x = mp.mpf(token)
        if not mp.isfinite(x) or x < 0:
            parser.error("reference x values must be finite and nonnegative")
        if x <= 1:
            m, e = low(x), low(x, energy=True)
        else:
            m, e = high(x, e1)
        row = dict(x=token, magnetization=mp.nstr(m, args.digits), energy_change_over_TB=mp.nstr(e, args.digits))
        if x == 1:
            independent, _ = high(x, e1)
            if abs(m-independent) > mp.mpf(10)**(-args.digits):
                raise ArithmeticError("low-field series and high-field integral disagree at crossover")
            row["high_integral_magnetization"] = mp.nstr(independent, args.digits)
            row["crossover_difference"] = mp.nstr(abs(m-independent), 5)
        rows.append(row)
    print(json.dumps(dict(digits=args.digits, convention="b=full Zeeman splitting; x=b/T_B; T_B=2*T1",
                          rows=rows), indent=2))


if __name__ == "__main__":
    main()
