# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional fixed-grid controls; passing does not certify outer convergence."""
import json
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad, trapezoid

from vertex_numpy import vertex_amplitudes
from vertex_quadrature import (MassiveKinematics, frequency_cut, integrated,
                               momentum_nodes, spinon_energy)


class VertexQuadrature(unittest.TestCase):
    def test_geometry_and_independent_sech_energy(self):
        for delta in (1.5, 4, 40):
            model = MassiveKinematics(delta)
            angles = np.linspace(-1.55, 1.55, 57)
            p = model.momenta(angles)
            np.testing.assert_allclose(model.angles(p), angles, atol=2e-13)
            np.testing.assert_allclose(model.energy(p), spinon_energy(angles, model.eta), rtol=2e-12)
            self.assertAlmostEqual(quad(model.angular_jacobian, 0, np.pi, epsabs=1e-11)[0], 1., places=11)

    def test_both_sectors_and_momentum_windings(self):
        model = MassiveKinematics(1.5)
        p = np.array([[.03, .09, .13], [.9, 1.3, 2.4], [3.01, 3.03, 3.09]])
        k = .81
        angles, energy, jacobian, sigma = momentum_nodes(model, model.angles(p), k)
        resolved = model.momenta(angles)
        self.assertEqual(set(sigma), {-1, 1})
        shift = np.where(sigma == -1, 0., np.pi)
        winding = (resolved.sum(axis=1)-(np.pi-k+shift))/(2*np.pi)
        np.testing.assert_allclose(winding, np.round(winding), atol=1e-13)
        self.assertGreater(np.max(np.round(winding)), 0)
        np.testing.assert_allclose(energy, model.energy(resolved).sum(axis=1), rtol=1e-13)
        self.assertTrue(np.all(jacobian > 0))

    def test_complex_amplitudes_against_multiprecision_fixture(self):
        fixture = json.loads(Path(__file__).with_name('vertex-point.json').read_text())
        model = MassiveKinematics(float(fixture['delta']))
        angles = np.array([[float(x) for x in fixture['angles']]])
        for vacuum in (0, 1):
            actual = vertex_amplitudes(angles, model.q, vacuum, 8, 24, True)[0]
            expected = np.array([complex(float(re), float(im)) for re, im in fixture['values'][vacuum]])
            np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-14)

    def test_fixed_order_values_used_by_cpp_controls(self):
        model = MassiveKinematics(1.5)
        full = integrated(model, 4, 12, poles=8, products=24)
        self.assertAlmostEqual(full['weight'], .015307326980575908, delta=1e-13)
        self.assertAlmostEqual(full['first_moment'], .03594894479635669, delta=1e-12)
        two = frequency_cut(model, 2, np.pi/2, [1.8, 2.2, 2.6], .06, 128)
        np.testing.assert_allclose(two['density'], [2.140377108167388, 2.197348192146067, .17242886353960657], rtol=1e-10)
        four = frequency_cut(model, 4, 1.5, [1.8, 2.2, 2.6, 3.2], .06, 12)
        np.testing.assert_allclose(four['density'], [6.606952745684729e-06, .047679974335733656,
                                                  .09266643564633772, .015852815526838926], rtol=1e-10)
        self.assertAlmostEqual(four['static_weight'], .012093538495760878, delta=1e-12)
        self.assertAlmostEqual(four['first_moment'], .031926812239310653, delta=1e-12)
        self.assertAlmostEqual(sum(four['translation_sector_weights']), four['static_weight'], places=14)

    def test_broadening_preserves_full_axis_weight_without_rescaling(self):
        model = MassiveKinematics(4.)
        omega = np.linspace(-4., 4*model.I+4., 2401)
        cut = frequency_cut(model, 4, 1.5, omega, .3, 4)
        self.assertAlmostEqual(trapezoid(cut['density'], omega)/(2*np.pi), cut['static_weight'], places=11)
        window = omega > cut['first_moment']/cut['static_weight']
        self.assertLess(trapezoid(np.asarray(cut['density'])[window], omega[window])/(2*np.pi), .9*cut['static_weight'])
        with self.assertRaises(ValueError):
            frequency_cut(model, 4, 1.5, [2.], .1, 3)
        with self.assertRaises(ArithmeticError):
            frequency_cut(model, 2, 1e-30, [2.], .1, 4)


if __name__ == '__main__':
    unittest.main()
