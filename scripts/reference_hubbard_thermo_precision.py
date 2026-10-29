#!/usr/bin/env python3
"""Independent MPFR-test references; requires mpmath, never used by the solvers.

Example: --kind lines --digits 128; --kind doped --digits 112 --nodes 96.
Repeat at larger precision and, for doped backgrounds, more nodes. Agreement
is numerical evidence, not an interval certificate. Output is JSON on stdout.
"""
import argparse
import json
import time
import mpmath as mp


def lines(digits):
    # Original Fourier-Bessel integrals, U=4, spin rapidity=2, bare k=1.
    cutoff = mp.mpf(digits + 12) * mp.log(10)
    panels = [mp.mpf(i) for i in range(0, int(cutoff), 8)] + [cutoff]
    def integral(f):
        return mp.quad(f, panels)
    es = 2 * integral(lambda w: mp.besselj(1, w) * mp.cos(2*w) / (w*mp.cosh(w)))
    ps = mp.pi/2 - integral(lambda w: mp.besselj(0, w)*mp.sin(2*w)/(w*mp.cosh(w)))
    ec = 2 + 2*mp.cos(1) + 2*integral(
        lambda w: mp.besselj(1, w)*mp.cos(w*mp.sin(1))*mp.exp(-w)/(w*mp.cosh(w)))
    pc = mp.pi/2 - 1 - integral(
        lambda w: mp.besselj(0, w)*mp.sin(w*mp.sin(1))*mp.exp(-w)/(w*mp.cosh(w)))
    # U=1/4 minimum holon energy: positive odd K1 series. Choose a truncation
    # beyond the requested decimal accuracy rather than relying on cancellation.
    terms = max(7, int((digits+12)*mp.log(10)/(16*mp.pi)) + 1)
    gap = 8/mp.pi * mp.fsum(mp.besselk(1, 8*mp.pi*n)/n for n in range(1,2*terms,2))
    return dict(spin_energy=es, spin_momentum=ps, charge_energy=ec,
                charge_momentum=pc, weak_charge_gap=gap)


def doped(nodes):
    # Fix U=40 and Q=1 and derive density, rather than fit the Fermi rapidity.
    u, q = mp.mpf(10), mp.mpf(1)
    xs, ws = mp.gauss_quadrature(nodes, 'legendre')
    k = [q*(x+1)/2 for x in xs]
    w = [q*x/2 for x in ws]
    sine, cosine = list(map(mp.sin,k)), list(map(mp.cos,k))
    def kernel(x):
        a = 1 + 1j*x/(2*u)
        return mp.re(mp.digamma((a+1)/2)-mp.digamma(a/2))/(4*mp.pi*u)
    ar, ae = mp.eye(nodes), mp.eye(nodes)
    for i in range(nodes):
        for j in range(nodes):
            value = kernel(sine[i]-sine[j]) + kernel(sine[i]+sine[j])
            ar[i,j] -= cosine[i]*w[j]*value
            ae[i,j] -= cosine[j]*w[j]*value
    rho = mp.lu_solve(ar, mp.matrix([1/(2*mp.pi)]*nodes))
    h = mp.lu_solve(ae, mp.matrix([-2*c for c in cosine]))
    xi = mp.lu_solve(ae, mp.matrix([1]*nodes))
    edge = [w[j]*cosine[j]*(kernel(mp.sin(q)-sine[j])+kernel(mp.sin(q)+sine[j])) for j in range(nodes)]
    mu = (-2*mp.cos(q)+mp.fsum(edge[j]*h[j] for j in range(nodes))) / (1+mp.fsum(edge[j]*xi[j] for j in range(nodes)))
    epsilon = [h[j]-mu*xi[j] for j in range(nodes)]
    return dict(density=2*mp.fsum(w[j]*rho[j] for j in range(nodes)), mu=mu,
                e0=-4*mp.fsum(w[j]*cosine[j]*rho[j] for j in range(nodes)),
                spin_energy_center=-mp.fsum(w[j]*cosine[j]*epsilon[j]/mp.cosh(mp.pi*sine[j]/(2*u))
                                            for j in range(nodes))/(2*u))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=('lines','doped'), default='lines')
    parser.add_argument('--digits', type=int, default=96)
    parser.add_argument('--nodes', type=int, default=64)
    args = parser.parse_args()
    if args.digits < 30 or args.nodes < 4:
        parser.error('require digits>=30 and nodes>=4')
    mp.mp.dps = args.digits
    start = time.monotonic()
    values = lines(args.digits) if args.kind=='lines' else doped(args.nodes)
    print(json.dumps(dict(kind=args.kind, digits=args.digits, nodes=args.nodes if args.kind=='doped' else None,
                          mpmath=mp.__version__, seconds=time.monotonic()-start,
                          values={k:mp.nstr(v,args.digits-6) for k,v in values.items()}), indent=2))


if __name__ == '__main__':
    main()
