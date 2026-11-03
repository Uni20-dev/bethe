# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
import json
from pathlib import Path
import unittest

import mpmath as mp
from vertex_reference import Vertex


class VertexReference(unittest.TestCase):
    def test_complex_amplitude_fixture(self):
        fixture = json.loads(Path(__file__).with_name('vertex-point.json').read_text())
        with mp.workdps(100):
            eta = mp.acosh(mp.mpf(fixture['delta']))
            beta = [mp.pi * mp.mpf(phi) / eta for phi in fixture['angles']]
            vertex = Vertex(eta)
            for vacuum in (0, 1):
                result = vertex.amplitudes(beta, vacuum, fixture['residue_levels'])
                for actual, expected in zip(result, fixture['values'][vacuum]):
                    self.assertLess(abs(actual - mp.mpc(*expected)), mp.mpf('1e-93'))

    def test_recurrence_against_literal_residues(self):
        with mp.workdps(80):
            eta = mp.acosh(mp.mpf('1.5'))
            vertex = Vertex(eta)
            beta = [mp.pi * mp.mpf(phi) / eta for phi in ('-1.2', '-.2', '.6', '1.3')]
            _, _, u = vertex.common(beta, 0)
            for vacuum in (0, 1):
                for ell, pole, level in ((1, 2, 1), (2, 0, 0), (3, 1, 3)):
                    a = vertex.residue(u, ell, vacuum, pole, level)
                    b = vertex.residue(u, ell, vacuum, pole, level + 1)
                    ratio = vertex.residue_ratio(u, ell, vacuum, pole, level)
                    self.assertLess(abs(b / a - ratio), mp.mpf('1e-70'))


if __name__ == '__main__':
    unittest.main()
