#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional mpmath oracle for the massive XXZ transverse two-spinon DSF.

Caux--Mossel--Perez Castillo (2008), PDF Eqs. (44), (51), (56)
(arXiv HTML Eqs. (49), (56), (61)). Physical momentum k = pi-Q.
The theta-A harmonic tail is summed analytically, not truncated.
Original two-spinon vertex amplitudes supply a separate momentum-root control;
a NumPy full-measure integral checks intensity and the first moment. Neither
imports research-only four-spinon contours, data, or compiled Bethe kernels.
"""
import argparse
import mpmath as mp


class TwoSpinon:
    def __init__(self, delta):
        self.delta = mp.mpf(delta)
        if not mp.isfinite(self.delta) or self.delta <= 1:
            raise ValueError('require finite Delta>1')
        self.eta = mp.acosh(self.delta)
        self.nome = mp.exp(-self.eta)
        self.K = mp.pi * mp.jtheta(3, 0, self.nome)**2 / 2
        dual = mp.exp(-mp.pi**2/self.eta)
        self.kp = ((mp.jtheta(2, 0, dual)/mp.jtheta(3, 0, dual))**2 if self.eta < mp.pi else
                   (mp.jtheta(4, 0, self.nome)/mp.jtheta(3, 0, self.nome))**2)
        self.modulus2 = 1-self.kp**2
        self.I = self.K*mp.sinh(self.eta)/mp.pi
        self.kappa = (1-self.kp)/(1+self.kp)

    def theta_ratio(self, beta):
        x = 2*mp.pi*beta/self.K
        total = 2*mp.log(2*abs(mp.sin(x/2)))
        for n in range(1, 1000000):
            t = mp.exp(-2*n*self.eta)
            total -= 2*t*((-1+2*t+t*t)*mp.cos(n*x)-2)/(n*(1+t)*(1-t*t))
            if t < mp.eps**2:
                break
        else:
            raise ArithmeticError('theta-A budget exhausted')
        # Neville theta_d is normalized Jacobi theta_3 (DLMF 20.1).
        d = mp.jtheta(3, mp.pi*beta/(2*self.K), self.nome)/mp.jtheta(3, 0, self.nome)
        return mp.exp(total)/d**2

    def sheets(self, momentum):
        k = min(mp.mpf(momentum), 2*mp.pi-mp.mpf(momentum))
        Q = mp.pi-k
        critical = mp.acos(self.kappa)
        def plus(q):
            if q < critical:
                return None
            upper = 2*self.I*mp.sqrt(mp.sin(q/2)**2+self.kp**2*mp.cos(q/2)**2)
            lower = (self.I*(1+self.kp)*mp.sin(q) if q <= mp.pi-critical else
                     2*self.I*mp.sqrt(mp.cos(q/2)**2+self.kp**2*mp.sin(q/2)**2))
            return lower, upper
        return {1: plus(Q), -1: plus(k)}

    def density(self, momentum, omega):
        momentum, omega = mp.mpf(momentum), mp.mpf(omega)
        Q = mp.pi-min(momentum, 2*mp.pi-momentum)
        w = omega/self.I
        result = mp.mpf(0)
        for sign, interval in self.sheets(momentum).items():
            if interval is None or not interval[0] < omega < interval[1]:
                continue
            if Q in (0, mp.pi):
                z = w/2
                beta = mp.ellipf(mp.asin(mp.sqrt((1-z*z)/self.modulus2)), self.modulus2)
                result += (self.kp if sign == 1 else 1)*self.theta_ratio(beta)*w/mp.sqrt((4-w*w)*(w*w-4*self.kp**2))
                continue
            c, s = mp.cos(Q), mp.sin(Q)
            w0 = (1+self.kp)*s
            T = mp.sqrt((w*w-self.kappa**2*w0*w0)*(w*w-w0*w0))
            z = (1+sign*c)/s*mp.sqrt((w*w-self.kappa*w0*w0+T)/(w*w+self.kappa*w0*w0-T))
            beta = mp.ellipf(mp.asin(mp.sqrt((1-z*z)/self.modulus2)), self.modulus2)
            W = mp.sqrt(self.kappa**2*w0**4/w**4-(T/w**2+sign*c)**2)
            numerator = (self.kp*(w*w+self.kappa*w0*w0+T) if sign == 1 else
                         w*w-self.kappa*w0*w0+T)
            result += (1+sign*c)*self.theta_ratio(beta)*numerator/(2*w*T*W)
        return result/self.I

    def density_from_vertex(self, momentum, omega):
        """Independent momentum-root search and original q-product amplitudes.

        Brackets simple roots on a mesh; for validation at interior points,
        not an exhaustive near-threshold root finder. No theta-A series,
        analytic inverse-dn solution, or sheet Jacobian is used.
        """
        vertex = TwoSpinonVertex(self.eta)
        Q = mp.pi-mp.mpf(momentum)
        omega = mp.mpf(omega)
        energy = lambda p: self.I*mp.sqrt(mp.sin(p)**2+self.kp**2*mp.cos(p)**2)
        velocity = lambda p: self.I**2*self.modulus2*mp.sin(p)*mp.cos(p)/energy(p)
        total = mp.mpf(0)
        for shift, sign in ((0, -1), (mp.pi, 1)):
            P = (Q+shift) % (2*mp.pi)
            low, high = max(0, P-mp.pi), P/2
            function = lambda p: energy(p)+energy(P-p)-omega
            last, f_last = low, function(low)
            for i in range(1, 257):
                end = low+(high-low)*i/256
                f_end = function(end)
                if f_last*f_end < 0:
                    root = mp.findroot(function, (last, end), solver='anderson')
                    momenta = [root, P-root]
                    beta = [mp.ellipf(p-mp.pi/2, self.modulus2) for p in momenta]
                    arguments = [-mp.pi**2*b/(2*self.K*self.eta) for b in beta]
                    x0, x1 = [vertex.amplitude(arguments, vacuum=i) for i in (0, 1)]
                    jacobian = energy(root)*energy(P-root)/self.I**2*abs(velocity(root)-velocity(P-root))
                    total += (mp.pi/(2*self.K))**2*abs(x0-sign*x1)**2/jacobian
                last, f_last = end, f_end
        return total


class TwoSpinonVertex:
    """Original two-spinon q-product amplitude (hep-th/9702028, Eqs. 2.6-2.7).

    This has no contour integral, reduced theta-A factor or energy Jacobian.
    Kept separate from the reduced DSF as a pointwise normalization control.
    """
    def __init__(self, eta):
        self.eta = mp.mpf(eta)
        if not mp.isfinite(self.eta) or self.eta <= 0:
            raise ValueError('require finite eta>0')
        self.q = -mp.exp(-self.eta)
        self.r = self.q**4

    @staticmethod
    def product(z, base):
        return mp.qp(z, base, maxterms=100000)

    def double_product(self, z):
        total, power = mp.mpc(0), z
        for h in range(100000):
            term = (h + 1)*mp.log1p(-power)
            total += term
            if abs(term) < mp.eps/100:
                return mp.exp(total)
            power *= self.r
        raise ArithmeticError('double-product budget exhausted')

    def amplitude(self, beta, vacuum=0):
        if len(beta) != 2 or vacuum not in (0, 1):
            raise ValueError('require two rapidities and vacuum 0 or 1')
        q, r = self.q, self.r
        qp = self.product
        dp = self.double_product
        xi = [1j*mp.exp(-1j*self.eta*b/mp.pi) for b in beta]
        u = [-x*x for x in xi]
        a = u[1]/u[0]
        value = (-q)**(1-vacuum)*qp(q*q, r)*qp(r, r)*(dp(r)/dp(q**6))**2
        for x in u:
            value /= qp(-q*x, r)*qp(-q**3/x, r)
        value *= dp(r*a)*dp(1/a)/(dp(q**6*a)*dp(q*q/a))
        z = -q**(4*vacuum)/(u[0]*u[1])
        theta = qp(q**8, q**8)*qp(z, q**8)*qp(q**8/z, q**8)
        return value*xi[0]**(1-vacuum)*xi[1]**(2-vacuum)*theta


def integrated_vertex_measure(delta, order=64, products=64):
    """Independent binary64 integration before resolving momentum and energy.

    Periodic staggered grids in both spectral angles avoid coincident nodes.
    Includes the two-vacuum average and 1/2! completeness factor. Product and
    grid truncations are explicit controls, not certified error bounds. This
    direct-product control is intended for massive interior/Ising checks.
    """
    import numpy as np
    if not np.isfinite(delta) or delta < 1.2 or order < 4 or products < 4:
        raise ValueError('measure control requires Delta>=1.2, order/products>=4')
    eta = np.arccosh(delta)
    q, r = -np.exp(-eta), np.exp(-4*eta)
    a, b = np.meshgrid((np.arange(order)+.23)*np.pi/order-np.pi/2,
                      (np.arange(order)+.67)*np.pi/order-np.pi/2, indexing='ij')
    angles = np.column_stack((a.ravel(), b.ravel()))
    xi = 1j*np.exp(-1j*angles)
    u = -xi*xi

    def qp(z, base=r, double=False):
        value = np.ones_like(z, dtype=complex)
        for h in range(products):
            value *= (1-z*base**h)**(h+1 if double else 1)
        return value

    dp = lambda z: qp(z, double=True)
    common = qp(q*q)*qp(r)*(dp(r)/dp(q**6))**2
    for j in range(2):
        common = common/(qp(-q*u[:, j])*qp(-q**3/u[:, j]))
    separation = u[:, 1]/u[:, 0]
    common *= dp(r*separation)*dp(1/separation)/(dp(q**6*separation)*dp(q*q/separation))
    density = np.zeros(len(angles))
    for vacuum in (0, 1):
        z = -q**(4*vacuum)/np.prod(u, axis=1)
        theta = qp(q**8, q**8)*qp(z, q**8)*qp(q**8/z, q**8)
        amplitude = (-q)**(1-vacuum)*common*np.prod(xi**np.arange(1-vacuum, 3-vacuum), axis=1)*theta
        density += abs(amplitude)**2/4
    # Poisson form of elliptic dn, independent of the resolved sheet energies.
    beta = ((angles+np.pi/2) % np.pi-np.pi/2)*np.pi/eta
    period = np.pi**2/eta
    images = int(np.ceil(40/period))+1
    energy = np.zeros_like(beta)
    for image in range(-images, images+1):
        small = np.exp(-abs(beta+image*period))
        energy += 2*small/(1+small*small)
    energy *= np.pi*np.sinh(eta)/(2*eta)
    return float(density.mean()), float((density*energy.sum(axis=1)).mean())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--delta', default='1.5')
    parser.add_argument('--momentum', default='1.5')
    parser.add_argument('--omega', default='2')
    parser.add_argument('--digits', type=int, default=65)
    args = parser.parse_args()
    mp.mp.dps = args.digits
    sf = TwoSpinon(args.delta)
    print('density', mp.nstr(sf.density(args.momentum, args.omega), args.digits-5))
    print('sheets', sf.sheets(mp.mpf(args.momentum)))
