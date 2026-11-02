# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss

from spectral_spatial import cumulative_weights, spatial_weights
from spectral_time_transform import statistics


class Spatial(unittest.TestCase):
    def test_neighbours_from_static_structure_factor(self):
        x, w = leggauss(24)
        q, measure = np.pi*(x+1)/2, w/2
        static = .25 - .10*np.cos(q) + .04*np.cos(2*q)
        W = spatial_weights(q, [0, 1, 2, 3], measure, even=True)
        np.testing.assert_allclose(W@static, [.25, -.05, .02, 0], atol=3e-16)
        # A distinct first-frequency moment is NOT a neighbour correlator.
        first = .3*(1-np.cos(q))
        self.assertAlmostEqual(measure@first, .3)
        self.assertGreater(abs((W@first)[1] - (W@static)[1]), .09)

    def test_finite_fourier_sign_and_covariance(self):
        q = 2*np.pi*np.arange(8)/8
        W = spatial_weights(q, [-1, 0, 1, 2], np.full(8, 1/8))
        values = np.sin(q)
        np.testing.assert_allclose(W@values, [-.5j, 0, .5j, 0], atol=2e-16)
        samples = np.array([values + a*np.cos(q) + b for a, b in
                            [(-.2, .1), (.1, -.2), (.4, .3), (-.1, -.1)]])
        covariance = np.cov(samples, rowvar=False, ddof=1)/len(samples)
        result = statistics(W, samples.mean(axis=0), covariance)
        transformed = samples@W.T
        np.testing.assert_allclose(result['real'], transformed.real.mean(axis=0), atol=1e-16)
        np.testing.assert_allclose(result['imag'], transformed.imag.mean(axis=0), atol=1e-16)
        np.testing.assert_allclose(result['real_standard_error'],
                                   transformed.real.std(axis=0, ddof=1)/2, atol=2e-9)
        # Contrast rows retain common-sample cancellation.
        contrast = statistics((W[[0]]-W[[2]]).real, samples.mean(axis=0), covariance)
        self.assertEqual(contrast['real_standard_error'], [0])

    def test_cutoff_retains_original_denominator(self):
        x = [.1, .4, .7, 1.]
        W = cumulative_weights(x, [.1, .5, 1.], np.full(4, 1/4))
        # Zero contributions remain observations, not dropped samples.
        np.testing.assert_array_equal(W@np.array([2., 0., 4., 2.]), [.5, .5, 2.])
        samples = np.array([[1., 2., 3., 4.], [2., 1., 5., 0.], [0., 3., 4., 2.], [3., 0., 2., 3.]])
        mean, C = samples.mean(axis=0), np.cov(samples, rowvar=False)/len(samples)
        for operator in (W, W[[2]]-W[[1]], np.vstack([W[0], W[2]])):
            got = statistics(operator, mean, C)
            direct = samples@operator.T
            np.testing.assert_allclose(got['real'], direct.mean(axis=0), atol=1e-16)
            np.testing.assert_allclose(got['real_standard_error'],
                                       direct.std(axis=0, ddof=1)/2, atol=1e-15)
        # Bin means are already integrals: no extra count normalization.
        np.testing.assert_array_equal(cumulative_weights([1, 2, 3], [1, 3], [1, 1, 1]),
                                      [[1, 0, 0], [1, 1, 1]])

    def test_domain_covariance_and_invalid_inputs(self):
        W = cumulative_weights([1, 2], [1, 2], [1, 1])
        # Explicitly independent domains: add means AND covariance matrices.
        a, b = np.array([2., 3.]), np.array([4., 5.])
        ca, cb = np.array([[1., -.5], [-.5, 2.]]), np.array([[4., 1.], [1., 3.]])
        both = statistics(W, a+b, ca+cb)
        self.assertEqual(both['real'], [6, 14])
        np.testing.assert_allclose(np.array(both['real_standard_error'])**2, [5, 11])
        for args in (([0], [True], [1]), ([0], [1.5], [1]), ([0], [10**30+1], [1]),
                     ([0], [1], [1, 2]), ([np.inf], [1], [1])):
            with self.subTest(args=args), self.assertRaises(ValueError):
                spatial_weights(*args)
        for args in (([], [1], []), ([1], [2, 1], [1]), ([1], [1, 1], [1]),
                     ([1], [1], [np.nan]), ([1, 2], [1], [1])):
            with self.subTest(args=args), self.assertRaises(ValueError):
                cumulative_weights(*args)


if __name__ == '__main__':
    unittest.main()
