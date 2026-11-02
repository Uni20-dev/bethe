# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
import copy
import math
import unittest

import numpy as np

from example_xxx_two_spinon_time import LOWER, Q, UPPER, read_segment, time_signals


class PhysicalTimeExample(unittest.TestCase):
    def test_distinct_measures_window_and_moments(self):
        # An independently integrable square-root density: G(x)=1.
        x = np.array([.1, .2, .3, .4])
        spectrum = np.column_stack((LOWER+x*x, 1/(2*x), np.zeros(4)))
        lines = [dict(gap=LOWER+.001, weight=.2), dict(gap=LOWER+.04, weight=.3),
                 dict(gap=LOWER+.3, weight=.4)]
        h = 1e-5
        data, moments = time_signals(lines, spectrum, [0., -h, h, 1.])
        self.assertAlmostEqual(moments['finite_all']['m0'], .9, places=14)
        self.assertAlmostEqual(moments['finite_window']['m0'], .3, places=14)
        self.assertAlmostEqual(moments['thermo_window']['m0'], .3/(2*math.pi), places=14)
        self.assertAlmostEqual(moments['thermo_window']['m1'],
                               (LOWER*.3+(.4**3-.1**3)/3)/(2*math.pi), places=14)
        for name in moments:
            signal = data[name+'_real']+1j*data[name+'_imag']
            self.assertAlmostEqual(signal[0].real, moments[name]['m0'], places=14)
            self.assertLess(abs((signal[2]-signal[1])/(2*h)+1j*moments[name]['m1']), 1e-9)
        self.assertAlmostEqual(data['finite_window_real'][3]+1j*data['finite_window_imag'][3],
                               .3*np.exp(-1j*(LOWER+.04)), places=14)
        self.assertLess(max(data['phase_order_change_abs']), 1e-15)
        self.assertLess(max(data['coarse_grid_change_abs']), 1e-15)

    @staticmethod
    def exports():
        header = '\n'.join('# '+key+': '+value for key, value in {
            'Outcome': 'success', 'Status': 'converged', 'Bethe revision': 'test-only',
            'Uni20 revision': 'test-only', 'Command': 'synthetic test fixture',
            'Program': 'bethe-xxx-structure-factor-thermo', 'Precision': 'fp64',
            'Channel': 'zz', 'Exchange J': '1', 'Spectral points': '3',
            'Failed points': '0', 'Divergent threshold points': '0',
            'Family': 'Thermodynamic two-spinon contribution; partial DSF; no rescaling',
            'Density convention': 'S(q,w) itself; integrals use dw/(2*pi); raising=2*zz',
        }.items())+'\n'
        return {
            'spectrum': header+'q,omega,sector,density,estimated_error,evaluations,status\n'+
                        '\n'.join(f'{Q},{omega},2,1,0,100,converged' for omega in [1.6, 1.65, 1.7]),
            'continuum': header+f'q,lower,upper\n{Q},{LOWER},{UPPER}\n',
        }

    def test_failed_missing_wrong_sector_and_mismatched_exports_rejected(self):
        original = self.exports()
        self.assertEqual(read_segment(original, 3, 1.6, 1.7)[1].shape, (3, 3))
        for old, new in [('2,1,0,100,converged', '2,,0,100,converged'),
                         ('2,1,0,100,converged', '2,1,0,100,numerical_failure'),
                         ('2,1,0,100,converged', '4,1,0,100,converged'),
                         ('2,1,0,100,converged', '2,0,0,100,outside_continuum'),
                         ('2,1,0,100,converged', '2,-1,0,100,converged'),
                         ('# Channel: zz', '# Channel: raising'),
                         ('# Bethe revision: test-only', '# Bethe revision: other')]:
            with self.subTest(replacement=new):
                texts = copy.deepcopy(original)
                texts['spectrum'] = texts['spectrum'].replace(old, new, 1)
                with self.assertRaises(ValueError):
                    read_segment(texts, 3, 1.6, 1.7)
        texts = copy.deepcopy(original)
        texts['spectrum'] = '\n'.join(texts['spectrum'].splitlines()[:-1])
        with self.assertRaises(ValueError):
            read_segment(texts, 3, 1.6, 1.7)


if __name__ == '__main__':
    unittest.main()
