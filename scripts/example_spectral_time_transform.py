# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Analytic finite-window example; no model data or plotting dependency."""
import json

import numpy as np

from spectral_time_transform import weights


def calculate():
    # S(omega)=1 on [1,2], zero outside. This declared support is exact for
    # the toy spectrum; no claim about omitted spin-chain weight is made.
    omega = np.array([1., 1.1, 1.6, 2.])
    times = np.array([0., 1., 5., 20., 100.])
    density = np.ones(len(omega))
    value = weights(omega, times, order=64)@density
    exact = np.exp(-1.5j*times)*np.sinc(times/(2*np.pi))/(2*np.pi)
    m1 = (weights(omega, [0.], energy_power=1)@density)[0].real
    return {
        "spectrum": "Analytic example: S(omega)=1 on [1,2], zero outside",
        "convention": "C(t)=integral S(omega)*exp(-i*omega*t)*domega/(2*pi)",
        "first_moment": float(m1),
        "exact_first_moment": 3/(4*np.pi),
        "time_signal": [dict(time=float(t), real=float(z.real), imag=float(z.imag),
                             absolute_error=float(abs(z-reference)))
                        for t, z, reference in zip(times, value, exact)],
    }


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2, allow_nan=False))
