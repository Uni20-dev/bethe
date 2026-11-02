# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional binary64 finite-q moments and Gaussian spectral cuts.

SciPy elliptic geometry and batched NumPy amplitudes, with independent
product/residue/outer-grid controls and physical normalization. Not a
production convergence engine or a general isotropic-limit evaluator.
"""
import argparse
import json
from math import ceil, factorial

import mpmath as mp
import numpy as np
from scipy.special import ellipj, ellipkinc
from numpy.polynomial.legendre import leggauss

from vertex_numpy import vertex_amplitudes


def local_spinon_density(angles, t, poles, products):
    """Two-vacuum average including 1/n!, against product dphi/pi.

    The second result is an arithmetic screen, not a truncation error bound.
    """
    n = angles.shape[1]
    values, errors = np.zeros(len(angles)), np.zeros(len(angles))
    for vacuum in (0, 1):
        amp, indicator = vertex_amplitudes(angles, -t, vacuum, poles, products,
                                            recurrence=True, diagnostics=True)
        values += np.sum(abs(amp)**2, axis=1)/(2*factorial(n))
        errors += np.sum(2*abs(amp)*indicator+indicator**2, axis=1)/(2*factorial(n))
    return values, errors


def spinon_energy(angles, epsilon):
    """J=1 energy as a periodic sech sum (Poisson form of elliptic dn)."""
    beta = ((np.asarray(angles)+np.pi/2) % np.pi-np.pi/2)*np.pi/epsilon
    period = np.pi**2/epsilon
    images = ceil(40/period)+1
    total = np.zeros_like(beta)
    for image in range(-images, images+1):
        small = np.exp(-abs(beta+image*period))
        total += 2*small/(1+small*small)
    return np.pi*np.sinh(epsilon)/(2*epsilon)*total


class MassiveKinematics:
    """Physical momenta in [0,pi], spectral angles in [-pi/2,pi/2].

    Binary64 prototype: reject an unresolved complementary elliptic modulus
    instead of quietly replacing the massive dispersion by the XXX limit.
    """
    def __init__(self, delta):
        self.delta = float(delta)
        if not np.isfinite(self.delta) or self.delta <= 1:
            raise ValueError('require finite Delta > 1')
        with mp.workdps(40):
            eta = mp.acosh(self.delta)
            nome = mp.exp(-eta)
            self.K = float(mp.pi*mp.jtheta(3, 0, nome)**2/2)
            self.kp = float((mp.jtheta(4, 0, nome)/mp.jtheta(3, 0, nome))**2)
        self.eta = float(eta)
        self.q = -float(nome)
        self.m = 1-self.kp**2
        if self.m == 1 or abs(self.q) > .975:
            raise ValueError('binary64 elliptic kinematics unresolved; use a higher-precision implementation')
        self.I = self.K*np.sinh(self.eta)/np.pi

    def momenta(self, angles):
        return np.pi/2+ellipj(-2*self.K*np.asarray(angles)/np.pi, self.m)[3]

    def angles(self, momenta):
        return -np.pi*ellipkinc(np.asarray(momenta)-np.pi/2, self.m)/(2*self.K)

    def energy(self, momenta):
        return self.I*np.hypot(np.sin(momenta), self.kp*np.cos(momenta))

    def angular_jacobian(self, momenta):
        """|dphi/dp|/pi = 1/(2 K dn b)."""
        return self.I/(2*self.K*self.energy(momenta))

def periodic_nodes(dimension, order, batch=2048, offset=0.):
    """Batched deterministic trapezoidal grid with distinct angle offsets.

    No Monte Carlo and no discarded collision nodes. A shifted grid provides
    an independent discretization diagnostic; it is not a statistical error.
    """
    if dimension not in (1, 2, 3, 4) or order < 2 or batch < 1 or not np.isfinite(offset):
        raise ValueError('require dimension 1..4, order >=2, positive batch, finite offset')
    shifts = (np.array([.173, .389, .617, .823])[:dimension]+offset) % 1
    for start in range(0, order**dimension, batch):
        indices = np.arange(start, min(start+batch, order**dimension))
        grid = np.column_stack(np.unravel_index(indices, (order,)*dimension))
        yield np.pi*((grid+shifts)/order-.5)

def integrated(model, spinons, order, poles=8, products=24, batch=2048, offset=0.):
    """Direct full-angle integral W_n and M_n; independent of momentum elimination."""
    if spinons not in (2, 4):
        raise ValueError('spinons must be 2 or 4')
    result = np.zeros(4)
    for angles in periodic_nodes(spinons, order, batch, offset):
        values, screen = local_spinon_density(angles, -model.q, poles, products)
        energy = model.energy(model.momenta(angles)).sum(axis=1)
        result += [values.sum(), np.dot(values, energy), screen.sum(), np.dot(screen, energy)]
    result /= order**spinons
    return dict(weight=result[0], first_moment=result[1],
                arithmetic_indicator=result[2], first_moment_arithmetic_indicator=result[3])

def checked_momentum(momentum):
    if not np.isfinite(momentum):
        raise ValueError('physical momentum must be finite')
    folded = momentum % (2*np.pi)
    q = np.pi-folded
    if folded != 0 and (q == np.pi or q == -np.pi):
        raise ArithmeticError('physical momentum shift unresolved in binary64')
    return folded


def momentum_nodes(model, free_angles, momentum):
    """Resolve delta_2pi(Q+shift-sum(p)), Q=pi-k, for shift=0,pi.

    Taking the remainder includes EVERY winding of total four-spinon
    momentum. For fixed first n-1 momenta there is precisely one admissible
    last momentum in [0,pi), belonging to one of the two translation sectors.
    No other preimages occur since this interval is shorter than 2pi.
    """
    momentum = checked_momentum(momentum)
    p = model.momenta(free_angles)
    remainder = (np.pi-momentum-p.sum(axis=1)) % (2*np.pi)
    sigma = np.where(remainder < np.pi, -1, 1)
    last = np.where(remainder < np.pi, remainder, remainder-np.pi)
    angles = np.column_stack((free_angles, model.angles(last)))
    energy = model.energy(p).sum(axis=1)+model.energy(last)
    return angles, energy, model.angular_jacobian(last), sigma

def cut_nodes(model, spinons, momentum, order, poles=8, products=24, batch=2048, offset=0.):
    """Weighted energies for integral S_n(k,omega) domega/(2pi).

    S_n = pi^2/n! integral prod(dphi/pi) delta(omega-E)
          sum_shift delta_2pi(Q+shift-P) sum_pol |X0-sigma X1|^2.
    Thus a resolved node has static weight pi/(2 n!) J_last B_sigma.
    Its Gaussian spectral density is 2pi times this weight times G_sigma.
    """
    if spinons not in (2, 4) or not np.isfinite(momentum):
        raise ValueError('require two/four spinons and finite momentum')
    momentum = checked_momentum(momentum)
    for free in periodic_nodes(spinons-1, order, batch, offset):
        angles, energy, jacobian, sigma = momentum_nodes(model, free, momentum)
        a, ea = vertex_amplitudes(angles, model.q, 0, poles, products, True, True)
        b, eb = vertex_amplitudes(angles, model.q, 1, poles, products, True, True)
        coherent = a-sigma[:, None]*b
        indicator = ea+eb
        factor = np.pi/(2*factorial(spinons))*jacobian/order**(spinons-1)
        weights = factor*np.sum(abs(coherent)**2, axis=1)
        errors = factor*np.sum(2*abs(coherent)*indicator+indicator**2, axis=1)
        if not np.all(np.isfinite(weights)) or not np.all(np.isfinite(errors)):
            raise ArithmeticError('non-finite momentum-resolved node; no nodes discarded')
        yield energy, weights, errors, sigma

def pair_cut_nodes(model, spinons, momentum, order, poles=8, products=24, batch=2048, offset=0.):
    """Resolve the last pair with a balanced angular importance quadrature.

    For each translation sector P_pair=(Q+shift-P_free) mod 2pi. Put
    p in [max(0,P_pair-pi),min(pi,P_pair)] and p'=P_pair-p.
    If j(p)=|dphi/dp|/pi and L=integral j(p)dp on this interval, use
    the mixture [j(p)+j(p')]/(2L). This resolves either soft end without
    a divergent importance weight, and requires no energy root finding.
    Evaluate BOTH maps deterministically; do not assume permutation symmetry.
    """
    if spinons not in (2, 4) or not np.isfinite(momentum):
        raise ValueError('require two/four spinons and finite momentum')
    if order < 2 or order % 2 or batch < 1 or not np.isfinite(offset):
        raise ValueError('pair quadrature requires even order >=2, positive batch, finite offset; '
                         'odd Gauss rules contain a coincident pair needing a confluent formula')
    momentum = checked_momentum(momentum)
    nodes, gauss = leggauss(order)
    nodes, gauss = (nodes+1)/2, gauss/2
    free_count = order**(spinons-2)
    free_batches = (periodic_nodes(spinons-2, order, max(1, batch//order), offset)
                    if spinons == 4 else [np.empty((1, 0))])
    for free in free_batches:
        free = np.repeat(free, order, axis=0)
        unit = np.tile(nodes, len(free)//order)
        quadrature = np.tile(gauss, len(free)//order)/free_count
        p_free = model.momenta(free)
        for shift, sign in ((0., -1), (np.pi, 1)):
            pair = (np.pi-momentum+shift-p_free.sum(axis=1)) % (2*np.pi)
            low, high = np.maximum(0., pair-np.pi), np.minimum(np.pi, pair)
            # A zero-length interval at an isolated periodic seam has zero
            # measure; it is not a failed amplitude node being discarded.
            present = high > low
            if not np.any(present):
                continue
            phi_low, phi_high = model.angles(low[present]), model.angles(high[present])
            length = (phi_low-phi_high)/np.pi
            phi = phi_high+(phi_low-phi_high)*unit[present]
            p = model.momenta(phi)
            other = pair[present]-p
            ja, jb = model.angular_jacobian(p), model.angular_jacobian(other)
            # Half of 2L*j(p)*j(p')/(j(p)+j(p')) for EACH of the two maps.
            factor = np.pi/(2*factorial(spinons))*quadrature[present]*length*ja*jb/(ja+jb)
            energy = model.energy(p_free[present]).sum(axis=1)+model.energy(p)+model.energy(other)
            for exchange_pair in (False, True):
                last = np.column_stack((phi, model.angles(other)))
                if exchange_pair:
                    last = last[:, ::-1]
                angles = np.column_stack((free[present], last))
                a, ea = vertex_amplitudes(angles, model.q, 0, poles, products, True, True)
                b, eb = vertex_amplitudes(angles, model.q, 1, poles, products, True, True)
                coherent, indicator = a-sign*b, ea+eb
                weights = factor*np.sum(abs(coherent)**2, axis=1)
                errors = factor*np.sum(2*abs(coherent)*indicator+indicator**2, axis=1)
                if not np.all(np.isfinite(weights)) or not np.all(np.isfinite(errors)):
                    raise ArithmeticError('non-finite pair node; no nodes discarded')
                yield energy, weights, errors, np.full(len(weights), sign)

def frequency_cut(model, spinons, momentum, omega, width, order, poles=8, products=24,
                  batch=2048, offset=0., method='pair'):
    """Normalized Gaussian convolution, width = standard deviation in J=1 units.

    Do not renormalize a finite frequency window or its negative-frequency
    Gaussian tail. Static weights/moments below use the unbroadened nodes.
    """
    omega = np.asarray(omega, dtype=float)
    if omega.ndim != 1 or not np.all(np.isfinite(omega)) or not np.isfinite(width) or width <= 0:
        raise ValueError('require finite 1D frequencies and positive Gaussian width')
    if method not in ('pair', 'direct'):
        raise ValueError('cut method must be pair or direct')
    density = np.zeros_like(omega)
    density_screen = np.zeros_like(omega)
    static = moment = screen = 0.
    sectors = np.zeros(2)
    engine = pair_cut_nodes if method == 'pair' else cut_nodes
    for energy, weight, error, sigma in engine(model, spinons, momentum, order, poles, products, batch, offset):
        gaussian = np.exp(-.5*((omega[:, None]-energy)/width)**2)/(np.sqrt(2*np.pi)*width)
        density += 2*np.pi*(gaussian @ weight)
        density_screen += 2*np.pi*(gaussian @ error)
        static += weight.sum()
        moment += np.dot(energy, weight)
        screen += error.sum()
        sectors += [weight[sigma == s].sum() for s in (-1, 1)]
    return dict(density=density.tolist(), density_arithmetic_indicator=density_screen.tolist(),
                static_weight=static, first_moment=moment, arithmetic_indicator=screen,
                translation_sector_weights=sectors.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['moments', 'cut'])
    parser.add_argument('--delta', type=float, default=1.5)
    parser.add_argument('--spinons', type=int, choices=[2, 4], default=4)
    parser.add_argument('--order', type=int, default=4)
    parser.add_argument('--products', type=int, default=24)
    parser.add_argument('--residues', type=int, default=8)
    parser.add_argument('--momentum', type=float, default=1.5)
    parser.add_argument('--omega', type=float, nargs='+', default=[1.8, 2.2, 2.6])
    parser.add_argument('--width', type=float, default=.06)
    args = parser.parse_args()
    if not 2 <= args.order <= 128 or not 4 <= args.products <= 1000 or not 1 <= args.residues <= 1000:
        parser.error('require order 2..128, products 4..1000, residues 1..1000')
    model = MassiveKinematics(args.delta)
    if args.mode == 'moments':
        result = integrated(model, args.spinons, args.order, args.residues, args.products)
    else:
        result = frequency_cut(model, args.spinons, args.momentum, args.omega,
                               args.width, args.order, args.residues, args.products)
    print(json.dumps(dict(controls=vars(args), result=result,
                          scope='fixed-order reference, not certified convergence'), indent=2))


if __name__ == '__main__':
    main()
