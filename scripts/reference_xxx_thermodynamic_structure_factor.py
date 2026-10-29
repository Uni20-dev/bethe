#!/usr/bin/env python3
"""Optional mpmath oracle: different tail subtraction and tanh-sinh meshes.

Print I(rho) for the C++ regression points. Not a runtime dependency.
The C++ code subtracts 2*(1-exp(-2*x))*cos(a*x)/x; here use exp(-x),
whose analytic integral is log(1+1/a^2). The cutoff grows with precision.
For the MPFR regression references use --digits 105 or 130 --points 0.1 1
--panel-width 1. The printed kernel values are I; tests store log(rate)=-I.
"""
import argparse
import mpmath as mp

def transition_integral(rho, cutoff=160, panel_width=mp.mpf("0.5")):
    a = 4 * rho

    def remainder(x):
        t = mp.tanh(x)
        return ((1+t)/mp.tanh(2*x)*(-2*mp.sin(a*x/2)**2)
                + t+t*t + 2*mp.expm1(-x)*mp.cos(a*x))/x

    panels = int(mp.ceil(cutoff / panel_width))
    edges = [mp.mpf(i) * cutoff / panels for i in range(panels + 1)]
    return mp.log1p(1/a**2) + mp.quad(remainder, edges)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=65)
    parser.add_argument("--points", nargs="+", default=['0.0001', '0.1', '0.5', '1', '2', '5'])
    parser.add_argument("--panel-width", default="0.5")
    args = parser.parse_args()
    if args.digits < 30:
        parser.error("at least 30 digits are required")
    mp.mp.dps = args.digits
    width = mp.mpf(args.panel_width)
    if not mp.isfinite(width) or width <= 0:
        parser.error("panel width must be finite and positive")
    # Independent exp(-x) tail decays more slowly than production exp(-2x).
    cutoff = int(mp.ceil((args.digits + 20) * mp.log(10)))
    print(f"# tanh-sinh; exp(-x) subtraction; digits={args.digits}; cutoff={cutoff}; panel-width={width}")
    for text in args.points:
        print(text, mp.nstr(transition_integral(mp.mpf(text), cutoff, width), args.digits-5), flush=True)
    q, omega = mp.mpf('1.5'), mp.mpf('1.8')
    lower, upper = mp.pi/2*mp.sin(q), mp.pi*mp.sin(q/2)
    rho = mp.acosh(mp.sqrt((upper**2-lower**2)/(omega**2-lower**2)))/mp.pi
    density = mp.exp(-transition_integral(rho, cutoff, width))/(2*mp.sqrt(upper**2-omega**2))
    print('Szz(1.5,1.8)', mp.nstr(density, args.digits-5))


if __name__ == "__main__":
    main()
