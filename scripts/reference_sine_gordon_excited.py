"""Independent p=2 excited NLIE audit using closed kernels and Gauss meshes.

Requires NumPy/SciPy; use OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1.
No production Fourier tables, uniform trapezoid grid or root solver are used.
Reports bulk-subtracted energies, not vacuum-relative gaps. This is an
independent discretization of the audited two-hole NLIE, not separate physics.
"""
import argparse
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq


def kernel(z):
    z = np.asarray(z, dtype=complex)
    result = np.full(z.shape, 1 / (2 * np.pi**2), dtype=complex)
    np.divide(z, 2 * np.pi**2 * np.sinh(z), out=result, where=z != 0)
    return result


phase_x, phase_w = leggauss(96)


def phase(z):
    z = np.asarray(z)
    return np.pi * z * (kernel(z[..., None] * (1 + phase_x) / 2) @ phase_w)


def level(u, number, n, cutoff, eta):
    nodes, weights = leggauss(n)
    x, w = cutoff * nodes, cutoff * weights
    z = x + 1j * eta
    k0 = kernel(x[:, None] - x[None, :]) * w
    k2 = kernel(x[:, None] - x[None, :] + 2j * eta) * w
    driving = -1j * u * np.sinh(z)

    def solve(hole):
        source = -1j * (phase(z - hole) + phase(z + hole))
        correction = np.zeros(n, dtype=complex)
        for _ in range(1000):
            epsilon = driving + source + correction
            argument = 1 + np.exp(-epsilon)
            assert np.min(argument.real) > 0
            a = np.log1p(np.exp(-epsilon))
            next_correction = -k0 @ a + k2 @ a.conj()
            residual = np.max(np.abs(next_correction - correction))
            if residual < 3e-14:
                counting = 2 * np.imag(np.sum(w * kernel(hole - x - 1j * eta) * a))
                f = u * np.sinh(hole) + phase(2 * hole).real + counting - 2 * np.pi * number
                value = 2 * u * np.cosh(hole) - u / np.pi * np.imag(np.sum(w * np.sinh(z) * a))
                return f, value / u, residual
            correction = (correction + next_correction) / 2
        raise RuntimeError("independent NLIE iteration exhausted")

    h = brentq(lambda h: solve(h)[0], 0, np.arcsinh(2 * np.pi * number / u) + 1, xtol=2e-14)
    residual, energy, nonlinear = solve(h)
    return dict(u=u, number=number, n=n, cutoff=cutoff, eta=eta, hole=h,
                energy=energy, counting_residual=residual, nonlinear_residual=nonlinear)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=float, default=1.)
    parser.add_argument("--number", type=float, default=.5)
    parser.add_argument("--meshes", type=int, nargs="+", default=[128, 256, 512])
    args = parser.parse_args()
    for n in args.meshes:
        for eta, cutoff in [(.6, 9.), (.75, 10.)]:
            print(json.dumps(level(args.length, args.number, n, cutoff, eta)), flush=True)
