"""Independent high-precision SG references (requires mpmath).

Uses the closed p=2 kernel in rapidity space, not production Fourier quadrature.
The pair is asymptotic Bethe-Yang, not an exact finite-volume excitation.
"""
import argparse
import json
import mpmath as mp


def phase(theta):
    return mp.quad(lambda x: x / (mp.pi * mp.sinh(x)) if x else 1 / mp.pi, [0, theta])



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=90)
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    if args.digits < 30:
        parser.error('digits must be at least 30')
    mp.mp.dps = args.digits
    theta = mp.findroot(lambda x: 10 * mp.sinh(x) + phase(2 * x) - mp.pi, mp.mpf("0.3"))
    z = mp.mpc('0.7', '0.2')
    kernel = z / (2 * mp.pi**2 * mp.sinh(z))
    values = {
        "phase_p2": phase(mp.mpf('0.7')),
        "breather_energy": mp.hypot(2 * mp.sin(mp.pi * mp.mpf('0.4') / 2), mp.mpf('0.7')),
        "pair_rapidity": theta,
        "pair_energy": 2 * mp.cosh(theta),
        "kernel_p2_real": mp.re(kernel),
        "kernel_p2_imag": mp.im(kernel),
    }
    values = {k: mp.nstr(v, args.digits-6) for k, v in values.items()}
    if args.json:
        print(json.dumps(dict(digits=args.digits, mpmath=mp.__version__, values=values), indent=2))
    else:
        for name, value in values.items():
            print(name, '=', value)


if __name__ == '__main__':
    main()
