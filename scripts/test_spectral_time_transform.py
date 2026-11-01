# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent analytic integrals and correlated-replicate checks."""
from decimal import Decimal, localcontext
import math
import unittest

import numpy as np

from spectral_time_transform import coarser_indices, statistics, weights


def analytic_hat_moments(omega, offsets=None):
    """Integrate basis polynomials algebraically, without phase quadrature."""
    x = omega if offsets is None else np.sqrt(offsets)
    jac = np.ones(len(x)) if offsets is None else 2*x
    threshold = None if offsets is None else omega[0]-offsets[0]
    zeroth, first = np.zeros(len(x)), np.zeros(len(x))
    for i in range(len(x)-1):
        a, b = x[i:i+2]
        h = b-a
        zeroth[i:i+2] += jac[i:i+2]*h/2
        if threshold is None:
            pieces = [(2*a+b)/6, (a+2*b)/6]
        else:
            pieces = [threshold/2+a*a/2+a*h/3+h*h/12,
                      threshold/2+a*a/2+2*a*h/3+h*h/4]
        first[i:i+2] += jac[i:i+2]*h*np.asarray(pieces)
    return zeroth/(2*np.pi), first/(2*np.pi)


class SpectralTimeTests(unittest.TestCase):
    def test_constant_density_exact_phase_and_measure(self):
        omega = np.array([1., 1.1, 1.6, 2.])
        t = np.array([0., .1, 2., 15.25, -15.25])
        actual = weights(omega, t)@np.ones(len(omega))
        expected = np.exp(-1.5j*t)*np.sinc(t/(2*np.pi))/(2*np.pi)
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-16)

    def test_zero_time_derivative_and_independent_hat_moments(self):
        delta = np.array([1e-6, .002, .03, .4, 1.])
        omega = 1.57+delta
        for d in (None, delta):
            zero, first = analytic_hat_moments(omega, d)
            W = weights(omega, [0., 1e-6, -1e-6], offsets=d)
            np.testing.assert_allclose(W[0].real, zero, rtol=1e-13)
            np.testing.assert_allclose(-(W[1]-W[2]).imag/(2e-6), first, rtol=2e-12)
            weighted = weights(omega, [0.], offsets=d, energy_power=1)
            np.testing.assert_allclose(weighted[0].real, first, rtol=1e-13)

    def test_first_moment_multiplies_original_interpolant(self):
        # S(w)=w is reproduced by frequency-linear interpolation. Integrating
        # w*S gives (b^3-a^3)/3, unlike interpolating the nodal values of w*S.
        omega = np.array([1., 1.4, 2.])
        lo, edge, hi = 1.13, 1.61, 1.92
        parts = [weights(omega, [0.], window=interval, energy_power=1)[0].real
                 for interval in [(lo, edge), (edge, hi)]]
        np.testing.assert_allclose([p@omega for p in parts],
            [(edge**3-lo**3)/(6*np.pi), (hi**3-edge**3)/(6*np.pi)], rtol=2e-14)
        # Same distinction with threshold coordinates and nonzero curvature
        # of omega(x); G(x)=1+x is represented exactly.
        x = np.array([.01, .2, .7, 1.])
        w = 2+x*x
        density = (1+x)/(2*x)
        def primitive(y):
            return 2*(y+y*y/2)+y**3/3+y**4/4
        actual = weights(w, [0.], offsets=x*x, energy_power=1)@density
        self.assertAlmostEqual(actual[0].real, (primitive(1)-primitive(.01))/(2*np.pi), places=14)

    def test_long_time_requires_phase_order_independent_of_grid(self):
        t = np.array([200.])
        expected = np.exp(-1.5j*t)*np.sinc(t/(2*np.pi))/(2*np.pi)
        low = weights([1., 2.], t, order=24)@np.ones(2)
        high = weights([1., 2.], t, order=96)@np.ones(2)
        self.assertGreater(abs(low[0]-expected[0]), 1e-5)
        np.testing.assert_allclose(high, expected, rtol=1e-11, atol=1e-14)

    def test_grid_error_is_separate_from_phase_order(self):
        errors = []
        for intervals in (4, 8, 16):
            w = np.linspace(0, 1, intervals+1)
            value = (weights(w, [0.])@(w*w))[0].real
            errors.append(value-1/(6*np.pi))
            higher_order = (weights(w, [0.], order=48)@(w*w))[0].real
            self.assertAlmostEqual(value, higher_order, places=15)
        np.testing.assert_allclose(np.array(errors[:-1])/errors[1:], 4, rtol=1e-12)

    def check_statistics(self, operator, samples):
        mean = samples.mean(axis=0)
        cov = np.cov(samples, rowvar=False, ddof=1)/len(samples)
        transformed = samples@operator.T
        result = statistics(operator, mean, cov)
        np.testing.assert_allclose(result['real'], transformed.real.mean(axis=0), atol=1e-15)
        np.testing.assert_allclose(result['imag'], transformed.imag.mean(axis=0), atol=1e-15)
        for label, data in [('real', transformed.real), ('imag', transformed.imag)]:
            np.testing.assert_allclose(result[label+'_standard_error'],
                data.std(axis=0, ddof=1)/np.sqrt(len(samples)), atol=1e-15)
        centered = transformed-transformed.mean(axis=0)
        cross = (centered.real*centered.imag).sum(axis=0)/(len(samples)*(len(samples)-1))
        np.testing.assert_allclose(result['real_imag_covariance'], cross, atol=1e-15)

    def test_correlated_grids_times_and_sector_contrasts(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(size=(19, 8))@rng.normal(size=(8, 8))
        omega = np.array([1.01, 1.1, 1.4, 2.])
        W = weights(omega, [0., .2, 1.3], offsets=omega-1)
        V = weights(omega, [0., .2, 1.3])
        for operator in (W, W-V, W[[2]]-W[[1]]):
            self.check_statistics(operator, samples[:, :4])
        # Two sectors sampled with the SAME random stream: use the whole
        # block covariance, including the nonzero off-diagonal cross block.
        for operator in (np.hstack([W, W]), np.hstack([W, -W])):
            self.check_statistics(operator, samples)
        # Identical sector estimates cancel exactly. Dropping the cross block
        # instead assigns a positive uncertainty to their identically zero
        # difference. Identity blocks keep this control analytically exact.
        cov = np.block([[np.eye(4), np.eye(4)], [np.eye(4), np.eye(4)]])
        exact = statistics(np.hstack([W, -W]), np.zeros(8), cov)
        wrong = statistics(np.hstack([W, -W]), np.zeros(8), np.eye(8))
        np.testing.assert_allclose(exact['real_standard_error'], 0, atol=1e-10)
        self.assertGreater(max(wrong['real_standard_error']), .01)

    def test_singular_density_has_no_endpoint_fill(self):
        delta = np.array([1e-6, 1e-4, .01, 1.])
        value = weights(2+delta, [0.], offsets=delta)@(1/(2*np.sqrt(delta)))
        self.assertAlmostEqual(value[0].real, (1-.001)/(2*np.pi), places=14)

    def test_precise_threshold_offsets_survive_rounded_frequencies(self):
        delta = np.array([1e-14, 1e-12, 1e-10, 1e-8])
        omega = 2+delta
        self.assertFalse(np.array_equal(omega-2, delta))
        actual = weights(omega, [0.], offsets=delta)@(1/(2*np.sqrt(delta)))
        expected = (np.sqrt(delta[-1])-np.sqrt(delta[0]))/(2*np.pi)
        np.testing.assert_allclose(actual.real, expected, rtol=2e-14)

    def test_clipped_windows_and_coarsening(self):
        w = np.array([1.01, 1.1, 1.4, 2.])
        t = np.array([0., 2., 10.])
        for delta in (None, w-1):
            whole = weights(w, t, offsets=delta, window=(1.04, 1.8))
            left = weights(w, t, offsets=delta, window=(1.04, 1.25))
            right = weights(w, t, offsets=delta, window=(1.25, 1.8))
            np.testing.assert_allclose(whole, left+right, atol=1e-15)
        ids = coarser_indices(4, 2)
        self.assertEqual(ids, [0, 2, 3])
        self.assertEqual(coarser_indices(5, 99), [0, 4])
        W = weights(w, t, indices=ids)
        np.testing.assert_array_equal(W[:, 1], 0)
        np.testing.assert_allclose(W[:, ids], weights(w[ids], t), atol=1e-16)
        with self.assertRaisesRegex(ValueError, 'coverage'):
            weights(w, t, window=(1., 1.8))

    def test_bad_indices_are_not_silently_rounded(self):
        for indices in ([0., 1.9], [0., 2.], [False, True], [0, True],
                        [False, 2], (0, np.bool_(True)),
                        np.array([0, True], dtype=object), [0, 0], [2, 1],
                        [-1, 2], [0, 3], [0, 2**80], [1], [[0, 2]], 1):
            with self.subTest(indices=indices), self.assertRaises(ValueError):
                weights([1, 2, 3], [0], indices=indices)
        expected = weights([1, 2, 3], [0], indices=[0, 2])
        for indices in (np.array([0, 2], dtype=np.int32),
                        np.array([0, 2], dtype=np.uint64), (np.int64(0), 2)):
            np.testing.assert_array_equal(weights([1, 2, 3], [0], indices=indices), expected)

    def test_invalid_inputs_and_controls(self):
        for w, t in [([1, 1], [0]), ([2, 1], [0]), ([1, np.inf], [0]),
                     ([1, 2], [np.nan]), ([1, 2], []), ([1, 2], [1j])]:
            with self.subTest(w=w, t=t), self.assertRaises(ValueError):
                weights(w, t)
        for keyword, value in [('energy_power', -1), ('energy_power', 1.5),
                               ('order', 1), ('order', True), ('window', [1, 2, 3]),
                               ('offsets', [0, 1]), ('offsets', [1, 3]),
                               ('window', [1, np.inf])]:
            with self.subTest(keyword=keyword, value=value), self.assertRaises(ValueError):
                weights([1, 2], [0], **{keyword: value})
        for count, stride in [(1, 1), (4, 0), (4, -1), (True, 1), (4, 1.5)]:
            with self.assertRaises(ValueError):
                coarser_indices(count, stride)
        with self.assertRaises(OverflowError):
            weights([10., 11.], [0.], energy_power=10000)
        with self.assertRaises(OverflowError):
            weights([10., 11.], [1e308])

    def test_covariance_symmetrization_preserves_numeric_range(self):
        tiny = np.nextafter(0., 1.)
        expected = np.sqrt(tiny)*1e160
        result = statistics([[1e160]], [0.], [[tiny]])
        np.testing.assert_allclose(result['real_standard_error'], [expected], rtol=1e-15)
        # Normalizing the whole matrix and scaling it back also loses this
        # covariance. It must survive alongside a much larger diagonal.
        C = np.diag([tiny, 1e300])
        saved = C.copy()
        result = statistics([[1e160, 0], [0, 1e-150]], [0, 0], C)
        np.testing.assert_allclose(result['real_standard_error'], [expected, 1.], rtol=1e-15)
        np.testing.assert_array_equal(C, saved)
        result = statistics([[1e160, 1e160]], [0, 0], np.full((2, 2), tiny))
        np.testing.assert_allclose(result['real_standard_error'], [2*expected], rtol=1e-15)
        # Conversely, adding large entries before halving would overflow.
        C = np.array([[1.7e308, 9e307], [np.nextafter(9e307, np.inf), 1.7e308]])
        result = statistics([[1e-154, -1e-154]], [0, 0], C)
        np.testing.assert_allclose(result['real_standard_error'], [np.sqrt(1.6)], rtol=1e-15)

    def test_standard_errors_do_not_require_representable_variances(self):
        for coefficient, variance, expected in [(1e-100, 1e-130, 1e-165),
                                                 (1e-160, 1., 1e-160),
                                                 (1e-200, 1., 1e-200),
                                                 (1e200, 1., 1e200)]:
            with self.subTest(coefficient=coefficient, variance=variance):
                result = statistics([[coefficient]], [0.], [[variance]])
                np.testing.assert_allclose(result['real_standard_error'], [expected], rtol=2e-15)
                self.assertEqual(result['imag_standard_error'], [0.])
                self.assertEqual(result['real_imag_covariance'], [0.])
        # Both marginal variances exceed binary64's useful exponent range,
        # while the two standard errors and their covariance are representable.
        result = statistics([[complex(1e200, 1e-200)]], [0.], [[1.]])
        np.testing.assert_allclose(result['real_standard_error'], [1e200], rtol=2e-15)
        np.testing.assert_allclose(result['imag_standard_error'], [1e-200], rtol=2e-15)
        np.testing.assert_allclose(result['real_imag_covariance'], [1.], rtol=2e-15)

    def test_covariance_cancellation_does_not_overflow_its_allowance(self):
        for scale in (1e154, 1e200):
            result = statistics([[scale, -scale]], [0., 0.], np.ones((2, 2)))
            self.assertEqual(result['real_standard_error'], [0.])
            self.assertEqual(result['real_imag_covariance'], [0.])
        # A global normalization of C followed by rescaling would lose the
        # first independent contribution, even though it matters to this row.
        tiny = np.nextafter(0., 1.)
        result = statistics([[1e160, 1e-150]], [0., 0.], np.diag([tiny, 1e300]))
        expected = np.hypot(np.sqrt(tiny)*1e160, 1.)
        np.testing.assert_allclose(result['real_standard_error'], [expected], rtol=2e-15)
        with self.assertRaises(OverflowError):
            statistics([[complex(1e200, 1e200)]], [0.], [[1.]])
        with self.assertRaises(OverflowError):
            statistics([[1e300]], [0.], [[1e300]])
        # A negative covariance hidden by the initial global eigenvalue
        # allowance must still fail the observable's relative variance check.
        with self.assertRaisesRegex(ValueError, 'negative transformed variance'):
            statistics([[1e200, 0]], [0., 0.], np.diag([-1e-300, 1e300]))
        near = np.array([[1., 1.+1e-15], [1.+1e-15, 1.]])
        result = statistics([[1e200, -1e200]], [0., 0.], near)
        self.assertEqual(result['real_standard_error'], [0.])

    def test_extreme_covariance_against_decimal_oracle(self):
        # Full, correlated covariance at very different input scales. This
        # independent decimal matrix expression checks the fallback's cross
        # covariance as well as its square-root conversion and range rejection.
        d = np.array([math.ldexp(1., 500), math.ldexp(1., -500), 1.])
        C = np.outer(d, d) + np.diag(d*d)
        with localcontext() as context:
            context.prec = 2200
            dc = [[Decimal.from_float(float(v)) for v in row] for row in C]
            for re in (-900, -400, 0, 400):
                for im in (-900, -400, 0, 400):
                    real = np.array([math.ldexp(1., re), -math.ldexp(1., re-2), math.ldexp(1., re-4)])
                    imag = np.array([-math.ldexp(1., im), math.ldexp(1., im-3), math.ldexp(1., im-5)])
                    parts = [list(map(Decimal.from_float, values.tolist())) for values in (real, imag)]
                    forms = [sum((parts[a][i]*dc[i][j]*parts[b][j]
                                  for i in range(3) for j in range(3)), Decimal(0))
                             for a, b in ((0, 0), (1, 1), (0, 1))]
                    expected = [float(forms[0].sqrt()), float(forms[1].sqrt()), float(forms[2])]
                    with self.subTest(real_exponent=re, imaginary_exponent=im):
                        if not all(math.isfinite(x) for x in expected):
                            with self.assertRaises(OverflowError):
                                statistics([real+1j*imag], np.zeros(3), C)
                        else:
                            result = statistics([real+1j*imag], np.zeros(3), C)
                            actual = [result[key][0] for key in ('real_standard_error',
                                      'imag_standard_error', 'real_imag_covariance')]
                            np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0)

    def test_statistics_rejects_invalid_covariance_and_complex_means(self):
        W = np.array([[1, 1j]])
        for covariance in ([[1, 1], [0, 1]], [[1, 2], [2, 1]], [[-1e-30, 0], [0, 0]],
                           [[1, 0], [0, np.inf]], np.eye(2, dtype=complex)):
            with self.subTest(covariance=covariance), self.assertRaises(ValueError):
                statistics(W, [1, 2], covariance)
        with self.assertRaises(ValueError):
            statistics(W, [1+1j, 2], np.eye(2))
        with self.assertRaises(ValueError):
            statistics([], [], [])
        with self.assertRaises(OverflowError):
            statistics([[1e300]], [1e300], [[1.]])
        exact = statistics(W, [1, 2], np.zeros((2, 2)))
        self.assertEqual(exact['real_standard_error'], [0.])
        self.assertEqual(exact['imag_standard_error'], [0.])
        # Rank-deficient covariance from perfectly shared observations is valid.
        cancel = statistics([[1, -1]], [1, 1], np.ones((2, 2)))
        self.assertEqual(cancel['real_standard_error'], [0.])


if __name__ == '__main__':
    unittest.main()
