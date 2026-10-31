"""Scientific/export guards for the measured four-spinon convergence tutorial."""
from copy import deepcopy
from decimal import Decimal, localcontext
import json
import unittest

import plot_xxx_four_spinon_tutorial as tutorial


class FourSpinonTutorial(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(tutorial.DATA.read_text())

    def test_checked_complete_grid(self):
        self.assertIs(tutorial.validate(self.data), self.data)
        self.assertEqual(len(self.data['cases']), 8)

    def test_wrong_channel_normalization_is_rejected(self):
        bad = deepcopy(self.data)
        bad['controls']['channel'] = 'raising'
        with self.assertRaises(ValueError):
            tutorial.validate(bad)
        bad = deepcopy(self.data)
        with localcontext() as context:
            context.prec = 100
            for row in bad['cases']:
                for key in ('s2', 's4', 'sum'):
                    row[key] = str(2 * Decimal(row[key]))
        with self.assertRaisesRegex(ValueError, 'normalization'):
            tutorial.validate(bad)

    def test_missing_duplicate_or_mixed_run_is_rejected(self):
        for kind in ('missing', 'duplicate', 'revision', 'threads'):
            with self.subTest(kind=kind):
                bad = deepcopy(self.data)
                if kind == 'missing':
                    bad['cases'].pop()
                elif kind == 'duplicate':
                    bad['cases'][0] = bad['cases'][1]
                elif kind == 'revision':
                    bad['cases'][0]['source_commit'] = '0' * 40
                else:
                    bad['cases'][0]['threads'] += 1
                with self.assertRaises(ValueError):
                    tutorial.validate(bad)

    def test_missing_weight_or_nonfinite_values_are_rejected(self):
        for key, value in [('s4', '0'), ('sum', '1'), ('cpu_seconds', 'NaN'),
                           ('arithmetic_diagnostic', '-1'), ('density_attempts', 0)]:
            with self.subTest(key=key):
                bad = deepcopy(self.data)
                bad['cases'][0][key] = value
                with self.assertRaises(ValueError):
                    tutorial.validate(bad)


if __name__ == '__main__':
    unittest.main()
