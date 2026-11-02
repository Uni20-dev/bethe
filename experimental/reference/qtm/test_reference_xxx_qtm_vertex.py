#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Regressions for measures, full partitions and the physical contour boundary."""
import math
import unittest

import mpmath as mp

import reference_xxx_qtm as qtm
import reference_xxx_qtm_vertex as bridge
from reference_xxx_transition_rate import transition_rate


class QtmVertexTests(unittest.TestCase):
    def setUp(self):
        self.digits = mp.mp.dps
        mp.mp.dps = 32

    def tearDown(self):
        mp.mp.dps = self.digits

    def test_all_ordered_assignments_and_factorials(self):
        for n, count in ((2, 2), (4, 6), (6, 20)):
            terms = bridge.assignments(n)
            self.assertEqual(len(terms), count)
            self.assertEqual(len(set(terms)), count)
            for h, p in terms:
                self.assertEqual(set(h) | set(p), set(range(n)))
                self.assertFalse(set(h) & set(p))
                self.assertIn((p, h), terms)
            self.assertEqual(count*math.factorial(n//2)**2, math.factorial(n))
        for n in (True, 0, 3, 8, 2.0):
            with self.assertRaises(ValueError):
                bridge.assignments(n)

    def test_measure_conversion_without_fitted_scale(self):
        for n in (2, 4, 6):
            beta = [mp.mpf(j)/3 for j in range(n)]
            jacobian = mp.fprod(1/mp.cosh(b) for b in beta)
            # Convert an arbitrary W=7 in both directions. Distinguish Pauli
            # from Szz, and labelled momentum from unlabelled QTM amplitude.
            pauli_beta = 4*7/(2*mp.pi)**2*jacobian
            qtm_beta = 7*bridge.vertex_conversion(beta)/(math.factorial(n)*(2*mp.pi**2)**n)
            self.assertLess(abs(qtm_beta/pauli_beta-1), mp.mpf("1e-29"))

    def test_two_spinon_absolute_normalization(self):
        beta = list(map(mp.mpf, ("-1.1", ".7")))
        expected = (2*mp.pi**2*transition_rate((beta[0]-beta[1])/(2*mp.pi))
                    /mp.fprod(mp.cosh(b) for b in beta))
        result = bridge.partition_sum(beta, "1e-10")
        self.assertLess(abs(result["total"].real/expected-1), mp.mpf("1e-18"))
        self.assertLess(abs(result["total"].imag/expected), mp.mpf("1e-9"))
        self.assertGreater(abs(result["pieces"][0]["amplitude"].imag), 1)
        self.assertLess(abs(bridge.vertex_target(beta)["amplitude"]/expected-1), mp.mpf("1e-29"))

    def test_four_spinon_boundary_and_label_invariance(self):
        beta = ["-2.1", "-1.3", "-.4", ".3"]
        target = bridge.vertex_target(beta)["amplitude"]
        a, b = [bridge.partition_sum(beta, offset) for offset in ("1e-5", "1e-9")]
        self.assertLess(abs(b["total"].real/target-1), mp.mpf("1e-16"))
        self.assertLess(abs(bridge.richardson_real(a, b)/target-1), mp.mpf("1e-24"))
        other = bridge.partition_sum(beta[::-1], "1e-9")
        self.assertLess(abs(other["total"]/b["total"]-1), mp.mpf("1e-25"))
        self.assertGreater(b["cancellation_condition"], 100)

    def test_six_spinon_twenty_term_absolute_check(self):
        beta = ["-2.1", "-1.3", "-.4", ".3", "1.2", "2"]
        q = bridge.partition_sum(beta, "1e-10")
        v = bridge.vertex_target(beta)
        self.assertEqual(len(q["pieces"]), 20)
        self.assertLess(abs(q["total"].real/v["amplitude"]-1), mp.mpf("1e-17"))
        self.assertGreater(q["cancellation_condition"], 1000)
        self.assertLess(abs(q["one_per_complement_pair"].real/q["total"].real-mp.mpf(".5")),
                        mp.mpf("1e-5"))
        self.assertLess(q["complement_conjugacy_scaled_defect"], mp.mpf("1e-8"))
        # The two former fixed-scale possibilities are decisively different.
        self.assertGreater(abs(q["total"].real/(2*v["amplitude"])-1), mp.mpf(".49"))

    def test_jmp_and_prl_agree_on_shifted_strip(self):
        _, t = bridge.shifted_points(["-1.1", ".7"], ".15")
        a, b = qtm.jmp([t[0]-mp.j], [t[1]]), qtm.prl([t[0]-mp.j], [t[1]])
        self.assertLess(abs(a["amplitude"]/b["amplitude"]-1), mp.mpf("1e-26"))

    def test_physical_momentum_and_energy_mapping(self):
        # Boundary evaluation only for elementary functions; the full QTM
        # reference deliberately keeps all inputs strictly inside its strips.
        beta = mp.mpf(".7")
        t = beta/mp.pi
        self.assertLess(abs(qtm.momentum(t)-2*mp.atan(mp.exp(beta))), mp.mpf("1e-30"))
        hole = qtm.momentum(t-mp.j)
        self.assertLess(abs(mp.exp(mp.j*(hole+qtm.momentum(t)))-1), mp.mpf("1e-29"))
        # Pauli-Hamiltonian J=1/4 matches spin-Hamiltonian J=1.
        ep = -mp.pi/(2*mp.cosh(mp.pi*t))
        eh = -mp.pi/(2*mp.cosh(mp.pi*(t-mp.j)))
        self.assertLess(abs(eh+ep), mp.mpf("1e-30"))

    def test_invalid_inputs_and_no_silent_boundary_rounding(self):
        for offset in (0, 1, "-.1", "nan", "1e-100"):
            with self.assertRaises(ValueError):
                bridge.shifted_points(["-.7", ".3"], offset)
        with self.assertRaises(ValueError):
            bridge.partition_sum(["-.7", ".3"], ".1", route="wrong")
        with self.assertRaises(ValueError):
            bridge.richardson_real(dict(offset=mp.mpf(".1"), total=mp.mpc(1)),
                                   dict(offset=mp.mpf(".2"), total=mp.mpc(1)))


if __name__ == "__main__":
    unittest.main()
