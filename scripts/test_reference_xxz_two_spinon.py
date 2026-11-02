# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Optional independent scientific checks; requires mpmath and numpy."""
import unittest
import mpmath as mp
from reference_xxz_two_spinon import TwoSpinon


class MassiveTwoSpinon(unittest.TestCase):
    def test_original_vertex_amplitudes_and_momentum_roots(self):
        with mp.workdps(45):
            for delta, k, w in [('1.5', '1.5', '2'), ('2', '1', '2.3'), ('4', '2', '4')]:
                sf = TwoSpinon(delta)
                direct = sf.density_from_vertex(k, w)
                reduced = sf.density(k, w)
                self.assertLess(abs(direct/reduced-1), mp.mpf('1e-35'))

    def test_integrated_vertex_measure(self):
        import numpy as np
        from reference_xxz_two_spinon import integrated_vertex_measure
        for delta, expected in [(1.5, [.4846429124574824, .7491062383225591]),
                                (4., [.4997481118452324, 1.9999994798394005])]:
            coarse = integrated_vertex_measure(delta, 64, 32)
            fine = integrated_vertex_measure(delta, 128, 64)
            np.testing.assert_allclose(coarse, fine, rtol=1e-12)
            np.testing.assert_allclose(fine, expected, rtol=1e-12)



if __name__ == '__main__':
    unittest.main()
