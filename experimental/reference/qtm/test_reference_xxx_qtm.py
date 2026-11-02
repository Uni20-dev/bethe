#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Pointwise QTM convention, isotropic-limit and independent-route checks."""
import unittest

import mpmath as mp

import reference_xxx_qtm as qtm


def fixture(ell):
    u = [mp.mpc(x, "-.4") for x in ("-.7", ".8", "1.8")[:ell]]
    v = [mp.mpc(x, ".6") for x in (".3", "1.4", "2.4")[:ell]]
    if ell > 1:
        v[0] = mp.mpc("-.2", ".6")
    return u, v


class QtmReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with mp.workdps(32):
            cls.results = {}
            for ell in (1, 2, 3):
                u, v = fixture(ell)
                cls.results[ell] = qtm.jmp(u, v), qtm.prl(u, v)

    def setUp(self):
        self.precision = mp.mp.dps
        mp.mp.dps = 32

    def tearDown(self):
        mp.mp.dps = self.precision

    def assertClose(self, x, y, tolerance="1e-25"):
        self.assertLess(abs(x-y)/max(1, abs(x), abs(y)), mp.mpf(tolerance))

    def test_independent_routes_through_three_pairs(self):
        for ell, (jmp, prl) in self.results.items():
            with self.subTest(ell=ell):
                self.assertClose(jmp["amplitude"], prl["amplitude"])
                self.assertClose(jmp["density"], prl["density"])
                self.assertLess(jmp["max_quad_error"], mp.mpf("1e-25"))
                for result in (jmp, prl):
                    for matrix, determinant in zip(result["matrices"], result["null_determinants"]):
                        self.assertLess(abs(determinant)/max(1, mp.norm(matrix, mp.inf)**ell),
                                        mp.mpf("1e-25"))

    def test_one_pair_gauss_formula_and_first_nonzero_determinant(self):
        u, v = fixture(1)
        a = -mp.j*(v[0]-u[0])/2
        phi = mp.gamma(mp.mpf(".5")+a)/(mp.sqrt(mp.pi)*mp.gamma(1+a))
        bar = mp.gamma(mp.mpf(".5")-a)/(mp.sqrt(mp.pi)*mp.gamma(1-a))
        value, slope = qtm.hyper_value_and_slope(*qtm.hyper_parameters(u, v, 0))
        self.assertClose(value, phi)
        self.assertClose(slope, a*phi)
        matrix, derivative = qtm.prl_matrix(u, v)
        self.assertClose(matrix[0, 0], 0)
        self.assertClose(derivative[0, 0], mp.j*(v[0]-u[0])*bar)
        self.assertGreater(abs(derivative[0, 0]), 1)

    def test_literal_first_column_is_not_the_symmetric_pointwise_weight(self):
        result, prl = self.results[2]
        self.assertGreater(abs(result["literal_first_column_amplitude"]/prl["amplitude"]-1),
                           mp.mpf(".1"))
        # Independently average all distinguished-column choices. Eq. (282)
        # uses ell^2 and the first columns only AFTER exploiting integration
        # symmetry. Its pointwise comparison requires this label average.
        a, b = result["matrices"]
        da, db = result["modified_columns"]
        pieces = []
        for k in range(2):
            for m in range(2):
                ac, bc = a.copy(), b.copy()
                for j in range(2):
                    ac[j, k], bc[j, m] = da[j, k], db[j, m]
                pieces.append(4*mp.det(ac)*mp.det(bc))
        u, v = fixture(2)
        average = -qtm.sine_squared(u, v)*result["pair_factor"]*mp.fsum(pieces)/4
        self.assertClose(average, prl["amplitude"])

    def test_label_permutation_invariance(self):
        u, v = fixture(2)
        expected = self.results[2][0]["amplitude"]
        self.assertClose(qtm.prl(u[::-1], v)["amplitude"], expected)
        self.assertClose(qtm.prl(u, v[::-1])["amplitude"], expected)
        self.assertClose(qtm.jmp(u[::-1], v[::-1])["amplitude"], expected)

    def test_precision_and_quadrature_splits(self):
        with mp.workdps(42):
            u, v = fixture(2)
            jmp = qtm.jmp(u, v, maxdegree=8, split=False)
            prl = qtm.prl(u, v)
            self.assertClose(jmp["amplitude"], prl["amplitude"], "1e-35")
            self.assertClose(prl["amplitude"], self.results[2][1]["amplitude"])

    def test_nonuniform_q_series_tail_coefficient(self):
        u, v = fixture(2)
        a, b = qtm.hyper_parameters(u, v, 0)
        phi, slope = qtm.hyper_value_and_slope(a, b)
        errors = []
        for gamma in map(mp.mpf, (".08", ".04", ".02")):
            base = mp.exp(-4*gamma)
            finite = mp.qhyper([base**z for z in a], [base**z for z in b], base, base)
            errors.append(abs((finite-phi)/gamma/slope-1))
        self.assertLess(errors[0], mp.mpf(".02"))
        for previous, current in zip(errors, errors[1:]):
            self.assertLess(current, previous*mp.mpf(".6"))

    def test_unexpanded_finite_q_determinant_approaches_analytic_limit(self):
        u, v = fixture(2)
        slope = self.results[2][1]["determinant_slopes"][0]
        errors = [abs(mp.det(qtm.finite_q_matrix(u, v, g))/mp.mpf(g)/slope-1)
                  for g in (".08", ".04", ".02")]
        self.assertLess(errors[0], mp.mpf(".01"))
        self.assertLess(errors[-1], mp.mpf(".0005"))
        for previous, current in zip(errors, errors[1:]):
            self.assertLess(current, previous*mp.mpf(".27"))
            self.assertGreater(current, previous*mp.mpf(".23"))

    def test_multilinear_determinant_derivative_needs_no_inverse(self):
        a = mp.matrix([[1, 2], [2, 4]])
        b = mp.matrix([[3, 5], [7, 11]])
        self.assertEqual(mp.det(a), 0)
        self.assertClose(qtm.determinant_derivative(a, b), -1)

    def test_measure_and_input_guards(self):
        self.assertClose(qtm.measure(2), 1/(4*(2*mp.pi)**4))
        u, v = fixture(2)
        for holes, particles in (([], []), (u, v[:1]), ([u[0]]*2, v),
                                  (u, [v[0]]*2), ([mp.nan], v[:1]),
                                  ([mp.mpc(0, ".1")], v[:1]),
                                  ([mp.mpc(".3", "-.4")], [mp.mpc(".3", ".6")])):
            with self.subTest(holes=holes, particles=particles):
                with self.assertRaises(ValueError):
                    qtm.points(holes, particles)
        with self.assertRaisesRegex(ValueError, "total-rapidity"):
            qtm.prl([mp.mpc("-.7", "-.4"), mp.mpc(".8", "-.4")],
                    [mp.mpc("-.2", ".6"), mp.mpc(".3", ".6")])
        with self.assertRaises(ValueError):
            qtm.finite_q_matrix(u, v, 0)
        with self.assertRaises(ValueError):
            qtm.hyper_value_and_slope([mp.mpf(".5")], [mp.mpf("2")])


if __name__ == "__main__":
    unittest.main()
