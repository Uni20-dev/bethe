# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent finite-q two/four-spinon contour reference.

Independent implementation of Abada--Bougourzi--Si-Lakhal,
hep-th/9702028 Eqs. (2.6),(2.7), with one contour for n=4.
Residues are derived directly from the q-Pochhammer integrand, rather than
using its printed isotropic formula or any withdrawn-preprint input.
Requires mpmath; no Bethe runtime dependency.
"""
import argparse
import mpmath as mp


def qpochhammer(*args):
    # Close to q=-1, mpmath's default 50*prec product budget can expire
    # during its adaptive accurate multiplication. This raises only the
    # iteration ceiling; precision and convergence remain mpmath-controlled.
    return mp.qp(*args, maxterms=100000)


class Vertex:
    def __init__(self, epsilon):
        self.epsilon = mp.mpf(epsilon)
        self.q = -mp.exp(-self.epsilon)
        self.r = self.q**4

    def double_product(self, z):
        # Product (z;r,r)_infinity, evaluated as sum (h+1)*log(1-z*r**h).
        total = mp.mpc(0)
        power = z
        for h in range(100000):
            term = (h + 1) * mp.log1p(-power)
            total += term
            if abs(term) < mp.eps / 100:
                return mp.exp(total)
            power *= self.r
        raise ArithmeticError('double product budget exhausted')

    @staticmethod
    def theta(z, base):
        return qpochhammer(base, base) * qpochhammer(z, base) * qpochhammer(base/z, base)

    def common(self, beta, vacuum, *, remove_pair_zeros=()):
        q, r = self.q, self.r
        n = len(beta)
        removed = set(remove_pair_zeros)
        if len(removed) != len(remove_pair_zeros) or any(not 0 <= j < k < n for j, k in removed):
            raise ValueError('removed zero factors must be distinct ordered spinon pairs')
        xi = [1j * mp.exp(-1j*self.epsilon*b/mp.pi) for b in beta]
        u = [-x*x for x in xi]
        rho = self.double_product(r) / self.double_product(q**6)
        # Eq. (2.7), not its n=2,4-only simplification 1-i*n/2.
        value = (-q)**(-n*n//4+(3-vacuum)*n//2-1) * (1-q**-2)**(n//2-1)
        value *= qpochhammer(q*q, r) * qpochhammer(r, r)**(n-1) * rho**n
        for j in range(n):
            value /= qpochhammer(-q*u[j], r) * qpochhammer(-q**3/u[j], r)
            for k in range(j+1, n):
                a = u[k]/u[j]
                # D(z)/(1-z)=D(r*z)*(r*z;r)_infinity. This evaluates
                # the removable within-pair zero even at coincident roots.
                second = (self.double_product(r/a)*qpochhammer(r/a, r) if (j, k) in removed
                          else self.double_product(1/a))
                value *= (self.double_product(r*a) * second
                          / self.double_product(q**6*a) / self.double_product(q*q/a))
        return value, xi, u

    def polynomial(self, u, ell, v):
        q = self.q
        return (mp.fprod(q/v - 1/(q*u[j]) for j in range(ell))
                * mp.fprod(1/v - 1/u[j] for j in range(ell+1, len(u))))

    def residue(self, u, ell, vacuum, j, m):
        """Direct residue at v=u[j]*r**m, with the dv/v factor included."""
        q, r = self.q, self.r
        v = u[j]*r**m
        if m == 0 and j > ell:
            return mp.mpc(0)  # polynomial cancels this pole
        result = (-1)**m * r**(m*(m+1)//2) / (qpochhammer(r, r, m) * qpochhammer(r, r))
        result *= self.polynomial(u, ell, v)
        for a in range(4):
            if a != j:
                result /= qpochhammer(u[a]/v, r)
            result /= qpochhammer(q**-2*v/u[a], r)
        result *= qpochhammer(-v/q, q*q) * qpochhammer(-q**3/v, q*q) * (-v)**vacuum
        result *= self.theta(-q**(-2+4*vacuum)*v*v/mp.fprod(u), q**8)
        return result

    def contour_integrand(self, u, ell, vacuum, v):
        """Eq. (2.6)'s contour integrand, excluding dv/(2*pi*i*v)."""
        q, r = self.q, self.r
        result = self.polynomial(u, ell, v)
        for a in u:
            result /= qpochhammer(a/v, r)*qpochhammer(q**-2*v/a, r)
        result *= qpochhammer(-v/q, q*q)*qpochhammer(-q**3/v, q*q)*(-v)**vacuum
        return result*self.theta(-q**(-2+4*vacuum)*v*v/mp.fprod(u), q**8)

    def residue_ratio(self, u, ell, vacuum, j, m):
        """Successive-pole ratio from product shifts and theta quasi-periodicity."""
        q, r = self.q, self.r
        v = u[j]*r**m
        z = -q**(-2+4*vacuum)*v*v/mp.fprod(u)
        result = -r**(m+1)/(1-r**(m+1))
        result *= self.polynomial(u, ell, r*v)/self.polynomial(u, ell, v)
        result *= mp.fprod(1-q**-2*v/a for a in u)
        result /= mp.fprod(1-u[a]/(r*v) for a in range(4) if a != j)
        result *= (1+1/(q*v))*(1+q/v)/((1+v/q)*(1+q*v))
        return result * r**vacuum * (-1/z)

    def amplitudes(self, beta, vacuum=0, terms=256):
        n = len(beta)
        if n not in (2, 4):
            raise ValueError('only n=2 and n=4 are implemented')
        common, xi, u = self.common(beta, vacuum)
        if n == 2:
            return [common * mp.fprod(xi[j]**(j+1-vacuum) for j in range(n))
                    * self.theta(-self.q**(4*vacuum)/mp.fprod(u), self.q**8)]
        amplitudes = []
        for ell in range(4):
            integral = mp.mpc(0)
            for j in range(4):
                start = int(j > ell)
                term = self.residue(u, ell, vacuum, j, start)
                subtotal = term
                for m in range(start, terms-1):
                    term *= self.residue_ratio(u, ell, vacuum, j, m)
                    subtotal += term
                integral += subtotal
            monomial = mp.fprod(xi[j]**(j+1-int(j == ell)-vacuum) for j in range(4))
            amplitudes.append(common * monomial * integral)
        return amplitudes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epsilon', nargs='+', default=['0.4', '0.2', '0.1'])
    parser.add_argument('--dps', type=int, default=50)
    parser.add_argument('--terms', type=int, default=256)
    parser.add_argument('--rho', nargs=4, default=['-.6', '-.1', '.3', '.8'])
    args = parser.parse_args()
    mp.mp.dps = args.dps
    if args.dps < 25 or args.terms < 4 or any(not mp.isfinite(mp.mpf(e)) or mp.mpf(e) <= 0 for e in args.epsilon):
        parser.error('require dps>=25, terms>=4, finite epsilon>0')
    beta = [2*mp.pi*mp.mpf(r) for r in args.rho]
    if any(not mp.isfinite(b) for b in beta) or len(set(beta)) != 4:
        parser.error('require four finite distinct rapidities')
    print('Finite-q amplitude reference; no quadrature-convergence certificate', flush=True)
    print('rho', args.rho, 'dps', args.dps, 'residue terms', args.terms, flush=True)
    for epsilon in args.epsilon:
        vertex = Vertex(epsilon)
        for n in (2, 4):
            amplitudes = [vertex.amplitudes(beta[:n], vacuum, args.terms) for vacuum in (0, 1)]
            densities = [vertex.epsilon**n * sum(abs(a)**2 for a in row)
                         * mp.fprod(mp.cosh(b) for b in beta[:n]) for row in amplitudes]
            # This is the raw squared amplitude density, before completeness,
            # Fourier, spin-component or energy-convention factors.
            print('epsilon', epsilon, 'n', n, 'raw densities',
                  *[mp.nstr(d, 25) for d in densities], 'vacuum ratios',
                  [mp.nstr(b/a, 12) for a, b in zip(*amplitudes)], flush=True)


if __name__ == '__main__':
    main()
