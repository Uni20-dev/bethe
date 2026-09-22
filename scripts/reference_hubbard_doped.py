#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent doped-Hubbard fixtures for test_hubbard_doped.cpp.

Optional developer tool; requires mpmath, not used by the build or tests.
At U=4, Q=1, solve on the full interval using digamma/log-gamma kernels.
The production code instead uses an even-half mesh and Euler/Stirling kernels.
Default 72/96-node results agree to all 44 printed digits (48-digit arithmetic).
Run: python3 scripts/reference_hubbard_doped.py [--nodes 72 96] [--digits 48]
"""
import argparse
import mpmath as mp


def reference(nodes):
    u, cutoff = mp.mpf(1), mp.mpf(1)

    def kernel(x):
        z = 1j * x / (4 * u)
        return mp.re(mp.digamma(1 + z) - mp.digamma(mp.mpf(".5") + z)) / (4 * mp.pi * u)

    def primitive(x):
        z = 1j * x / (4 * u)
        return mp.im(mp.loggamma(1 + z) - mp.loggamma(mp.mpf(".5") + z)) / mp.pi

    x, weight = mp.gauss_quadrature(nodes, "legendre")
    x, weight = [cutoff * v for v in x], [cutoff * v for v in weight]
    sine, cosine = [mp.sin(v) for v in x], [mp.cos(v) for v in x]
    arho, aenergy = mp.matrix(nodes), mp.matrix(nodes)
    for i in range(nodes):
        for j in range(nodes):
            r = kernel(sine[i] - sine[j])
            arho[i, j] = int(i == j) - cosine[i] * weight[j] * r
            aenergy[i, j] = int(i == j) - cosine[j] * weight[j] * r
    rho = mp.lu_solve(arho, mp.matrix([1 / (2 * mp.pi)] * nodes))
    h = mp.lu_solve(aenergy, mp.matrix([-2 * v for v in cosine]))
    xi = mp.lu_solve(aenergy, mp.matrix([1] * nodes))
    edge = [weight[j] * cosine[j] * kernel(mp.sin(cutoff) - sine[j]) for j in range(nodes)]
    mu = (-2 * mp.cos(cutoff) + mp.fsum(edge[j] * h[j] for j in range(nodes))) / (
        1 + mp.fsum(edge[j] * xi[j] for j in range(nodes)))
    epsilon = [h[j] - mu * xi[j] for j in range(nodes)]
    density = mp.fsum(weight[j] * rho[j] for j in range(nodes))
    energy = -2 * mp.fsum(weight[j] * rho[j] * cosine[j] for j in range(nodes))
    rapidity = mp.mpf(".5")
    spin_energy = -mp.fsum(weight[j] * cosine[j] * epsilon[j] /
                          (4 * u * mp.cosh(mp.pi * (rapidity - sine[j]) / (2 * u))) for j in range(nodes))
    spin_momentum = 2 * mp.fsum(weight[j] * rho[j] *
                              mp.atan(mp.exp(-mp.pi * (rapidity - sine[j]) / (2 * u))) for j in range(nodes))
    print("nodes", nodes, flush=True)
    for name, value in [("density", density), ("mu_unshifted", mu), ("e0_unshifted", energy),
                        ("spin_energy(lambda=.5)", spin_energy), ("spin_momentum(lambda=.5)", spin_momentum)]:
        print(name, mp.nstr(value, mp.mp.dps - 4), flush=True)
    for k in [mp.mpf(".5"), mp.mpf(2)]:
        epsilon_k = -2 * mp.cos(k) - mu + mp.fsum(
            weight[j] * cosine[j] * kernel(mp.sin(k) - sine[j]) * epsilon[j] for j in range(nodes))
        momentum = k + 2 * mp.pi * mp.fsum(
            weight[j] * rho[j] * primitive(mp.sin(k) - sine[j]) for j in range(nodes))
        print("charge", k, "epsilon", mp.nstr(epsilon_k, mp.mp.dps - 4),
              "p_c", mp.nstr(momentum, mp.mp.dps - 4), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=int, nargs="+", default=[72, 96])
    parser.add_argument("--digits", type=int, default=48)
    args = parser.parse_args()
    if args.digits < 16 or min(args.nodes) < 4:
        parser.error("require at least 16 digits and 4 nodes")
    mp.mp.dps = args.digits
    for order in args.nodes:
        reference(order)
