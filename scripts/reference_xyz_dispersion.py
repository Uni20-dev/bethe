#!/usr/bin/env python3
"""Independent mpmath XYZ oracle using JKM (1973) Eqs. (7.8) and (7.11).

Uses Jacobi functions, not the production theta-derivative prefactor. Checks
bound E(Q) against the complex-rapidity dispersion before printing references.
Not a build dependency; moderate eta,t with sufficient digits are intended.
"""
import argparse
import json

import mpmath as mp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta", default="0.75")
    parser.add_argument("--t", default="1")
    parser.add_argument("--digits", type=int, default=65)
    args = parser.parse_args()
    if args.digits < 20:
        parser.error("require at least 20 digits")
    mp.mp.dps = args.digits + 30
    eta, t = mp.mpf(args.eta), mp.mpf(args.t)
    if not 0 < eta < 1 or not mp.isfinite(t) or t <= 0:
        parser.error("require 0<eta<1, finite t>0")
    if eta/(1-eta) > 1000:
        parser.error("reference enumeration limited to 1000 bound indices")
    q = mp.exp(-mp.pi*t)
    theta = lambda kind, u: mp.jtheta(kind, mp.pi*u, q)
    jx, jy, jz = [theta(kind, eta)/theta(kind, 0) for kind in (4, 3, 2)]
    ell = theta(2, 0)**2/theta(3, 0)**2
    big_k = mp.pi*theta(3, 0)**2/2
    q1 = mp.exp(-mp.pi*eta/t)
    k1 = mp.jtheta(2, 0, q1)**2/mp.jtheta(3, 0, q1)**2
    big_k1 = mp.pi*mp.jtheta(3, 0, q1)**2/2
    r = mp.jtheta(4, 0, q1)**2/mp.jtheta(3, 0, q1)**2
    maximum = jx*mp.ellipfun("sn", 2*big_k*eta, ell**2)*big_k1/(2*t*big_k)
    tolerance = mp.mpf(10)**(-args.digits)
    assert abs(mp.ellipfun("cn", 2*big_k*eta, ell**2)-jz/jx) < tolerance
    assert abs(mp.ellipfun("dn", 2*big_k*eta, ell**2)-jy/jx) < tolerance
    fmt = lambda x: mp.nstr(x, args.digits)
    result = dict(eta=args.eta, t=args.t, jx=fmt(jx), jy=fmt(jy), jz=fmt(jz),
                  maximum=fmt(maximum), gap=fmt(maximum*r), bound=[])
    s = 1
    while s*(1-eta) < eta:
        a = mp.ellipfun("sn", s*big_k1*(1-eta)/t, r*r)

        def energy(momentum):
            sine, cosine = mp.sin(momentum/2), mp.cos(momentum/2)
            return 2*maximum/a*mp.sqrt(sine*sine+r*r*a*a*cosine*cosine)*mp.sqrt(sine*sine+a*a*cosine*cosine)

        residuals = []
        for phi in (mp.mpf(0), mp.mpf("0.4"), mp.mpf("1.3")):
            z = big_k1/mp.pi*(phi+1j*mp.pi*((s+1)*eta-s)/t)
            dn = mp.ellipfun("dn", z, k1*k1)
            # am(z)=-i log(cn(z)+i sn(z)); Eq. (7.11b), analytically integrated.
            am = -1j*mp.log(mp.ellipfun("cn", z, k1*k1)+1j*mp.ellipfun("sn", z, k1*k1))
            momentum = -mp.pi-2*am.real
            residuals.append(abs(energy(momentum)-2*maximum*dn.real)/maximum)
        assert max(residuals) < tolerance, (s, residuals)
        result["bound"].append(dict(s=s, a=fmt(a), gap=fmt(2*maximum*r*a),
                                    energy_at_0_7=fmt(energy(mp.mpf("0.7"))),
                                    parametric_relative_error=mp.nstr(max(residuals), 5)))
        s += 1
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
