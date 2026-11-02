#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Regression checks for complex partition comparison and analytic boundary."""
import unittest
import json
from pathlib import Path
import tempfile

import mpmath as mp

import audit_xxx_qtm_partitions as audit
import reference_xxx_qtm as qtm
import reference_xxx_qtm_vertex as bridge


class ComplexPartitionTests(unittest.TestCase):
    def setUp(self):
        self.precision = mp.mp.dps
        mp.mp.dps = 30

    def tearDown(self):
        mp.mp.dps = self.precision

    def assertClose(self, a, b, tolerance="1e-24"):
        self.assertLess(audit.scaled_difference(a, b), mp.mpf(tolerance))

    def test_same_sine_factor_for_every_six_spinon_partition(self):
        beta = audit.FIXTURES["six-bulk"]
        for offset in (".1", ".25", ".6"):
            _, t = bridge.shifted_points(beta, offset)
            expected = audit.kinematic_factor(beta, offset)
            for h, p in bridge.assignments(6):
                self.assertClose(qtm.sine_squared([t[j]-mp.j for j in h], [t[j] for j in p]),
                                 expected)

    def test_complex_two_spinon_routes_and_boundary_without_fit(self):
        beta = audit.FIXTURES["two"]
        a = audit.matched_partition(beta, (0,), ".25")
        b = audit.matched_partition(beta, (1,), ".1", split=False, maxdegree=8)
        for r in (a, b):
            self.assertClose(r["routes"]["jmp"]["amplitude"], r["routes"]["prl"]["amplitude"])
        ca, cb = (r["routes"]["prl"]["reduced_coefficient"] for r in (a, b))
        self.assertClose(ca, mp.conj(cb))
        total = a["routes"]["prl"]["boundary_amplitude"]+b["routes"]["prl"]["boundary_amplitude"]
        self.assertClose(total, bridge.vertex_target(beta)["amplitude"])
        self.assertGreater(abs(ca.imag), 1)  # taking real/abs too early must be detectable

    def test_common_shift_invariance_and_complement_not_same_amplitude(self):
        beta = audit.FIXTURES["four"]
        coefficients = []
        for offset, h in ((".25", (0, 1)), (".1", (0, 1)), (".25", (2, 3))):
            _, t = bridge.shifted_points(beta, offset)
            p = tuple(j for j in range(4) if j not in h)
            result = qtm.prl([t[j]-mp.j for j in h], [t[j] for j in p])
            coefficients.append(audit.reduced_coefficient(result, "prl", 2))
        self.assertClose(coefficients[0], coefficients[1])
        self.assertClose(coefficients[0], mp.conj(coefficients[2]))
        self.assertGreater(abs(coefficients[0]-coefficients[2]), 1)

    def test_removing_sine_does_not_divide_by_a_small_or_zero_number(self):
        beta = ["-2", "-1", "1", "2"]
        self.assertLess(abs(audit.boundary_factor(beta)), mp.mpf("1e-55"))
        # The algebraic helper needs no amplitude or sine field at all.
        fake = dict(pair_factor=mp.mpc(2, 3), scalar=mp.mpf(7),
                    determinant_slopes=(mp.mpc(1, 2), mp.mpc(3, 4)))
        self.assertEqual(audit.reduced_coefficient(fake, "prl", 2),
                         fake["pair_factor"]*7*mp.fprod(fake["determinant_slopes"]))
        with self.assertRaises(ValueError):
            audit.reduced_coefficient(fake, "unknown", 2)

    def test_refined_unsplit_quadrature_resolves_the_four_spinon_coefficient(self):
        # At maxdegree=8 the discrepancy was about 2.25e-32 despite 42-digit
        # arithmetic. Degree 9 resolves it; this test must not accept that
        # under-resolved inner quadrature as a high-precision reference.
        with mp.workdps(42):
            result = audit.matched_partition(audit.FIXTURES["four"], (0, 1), ".1",
                                             maxdegree=9, split=False)
            self.assertLess(result["reduced_scaled_difference"], mp.mpf("1e-35"))

    def test_public_checkpoint_guards_and_fresh_reference(self):
        before = mp.mp.dps
        fresh = audit.fresh_reference("two", 40)
        self.assertEqual(mp.mp.dps, before)
        self.assertEqual(fresh["beta"], audit.FIXTURES["two"])
        self.assertEqual(len(fresh["records"][0]["pieces"]), 2)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"vertex.json"
            good = dict(schema="xxx-qtm-vertex-pointwise-v1", passed=True,
                        source_sha256=audit.public_source_hashes(), fixtures=[fresh])
            path.write_text(json.dumps(good))
            self.assertEqual(audit.verify_checkpoint(path), good)
            for mutate in (
                lambda d: d.update(passed=False),
                lambda d: d.update(schema="other"),
                lambda d: d.update(source_sha256={}),
                lambda d: d["source_sha256"].update({"../outside.py": "0"*64}),
                lambda d: d["source_sha256"].update({
                    "experimental/reference/qtm/reference_xxx_qtm.py": "0"*64}),
            ):
                changed = json.loads(json.dumps(good))
                mutate(changed)
                path.write_text(json.dumps(changed))
                with self.assertRaises(ValueError):
                    audit.verify_checkpoint(path)

    def test_decimal_encoding_and_invalid_partition(self):
        value = mp.mpc(".123456789012345678901234567", "-.98765432109876543210987654")
        self.assertEqual(audit.decode_complex(audit.encode(value)), value)
        for h in ((1, 0), (0,), (0, 4), (0, 0)):
            with self.assertRaises(ValueError):
                audit.matched_partition(audit.FIXTURES["four"], h, ".25")


if __name__ == "__main__":
    unittest.main()
