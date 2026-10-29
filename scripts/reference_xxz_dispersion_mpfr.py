#!/usr/bin/env python3
"""Independent thermodynamic XXZ references; requires mpmath.

Theta constants use mpmath's elliptic nome directly, without Bethe's modular
reconstruction. Massive bulk energies use the direct Fourier series even where
the solver uses a Poisson-resummed integral; gapless energies use the original
rapidity Fourier integral. Repeat with --dps 128 and 160 to check stability.
"""
import argparse
import json
import mpmath as mp


def evaluate(text):
    delta = mp.mpf(text)
    if delta > 1:
        eta = mp.acosh(delta)
        nome = mp.exp(-eta)
        maximum = mp.sinh(eta) * mp.jtheta(3, 0, nome)**2 / 2
        gap = mp.sinh(eta) * mp.jtheta(4, 0, nome)**2 / 2
        count = int(mp.ceil((mp.mp.dps + 12) * mp.log(10) / (2 * eta)))
        energy = delta / 4 - mp.sinh(eta) * (mp.mpf('0.5') + 2 * mp.fsum(
            1 / (1 + mp.exp(2*n*eta)) for n in range(1, count+1)))
    else:
        gamma = mp.acos(delta)
        maximum = mp.pi * mp.sin(gamma) / (2 * gamma)
        gap = mp.mpf(0)
        def kernel(w):
            if not w:
                return (mp.pi-gamma)/mp.pi
            return mp.sinh((mp.pi-gamma)*w) / (mp.sinh(mp.pi*w)*mp.cosh(gamma*w))
        energy = delta/4 - mp.sin(gamma)*mp.quad(kernel, [0, 1, mp.inf])
    dispersion = mp.sqrt((gap*mp.cos(mp.mpf('0.3')))**2 +
                         (maximum*mp.sin(mp.mpf('0.3')))**2)
    return {key: mp.nstr(value, mp.mp.dps) for key, value in
            [('maximum',maximum),('gap',gap),('energy_at_0.3',dispersion),('bulk',energy)]}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps', type=int, default=128)
    args=parser.parse_args()
    if args.dps < 80:
        parser.error('use at least 80 digits')
    mp.mp.dps=args.dps
    print(json.dumps(dict(dps=args.dps, mpmath=mp.__version__,
        values={delta: evaluate(delta) for delta in ('-0.75', '0.5', '1.0078125', '1.125', '2', '100')}), indent=2))
