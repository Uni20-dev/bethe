"""Scientific/export guards for the measured four-spinon convergence tutorial."""
from copy import deepcopy
from decimal import Decimal, localcontext
import json
import re
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

    def test_import_checks_tolerance_beyond_default_decimal_precision(self):
        record = json.loads((tutorial.ROOT / 'docs/tutorials/data/xxx-four-spinon-precision.json').read_text())
        case = next(row for row in record['cases'] if row['precision'] == 'bits-128')
        document = {'tables': {}}
        for name in ('spectrum', 'diagnostics'):
            names = list(case[name][0])
            document['tables'][name] = {
                'metadata': deepcopy(case['metadata']), 'summary': case['summary'],
                'columns': [{'id': name} for name in names],
                'rows': [[row[name] for name in names] for row in case[name]]}
        tutorial.read_case(document, 4, 8)
        with localcontext() as context:
            context.prec = 100
            changed = str(Decimal('1e-12') * (1 + Decimal('1e-30')))
        for table in document['tables'].values():
            table['metadata']['Two-spinon transition-integral tolerance'] = changed
        # A 28-digit Decimal division rounds this change away, but it is much
        # larger than the 128-bit format's permitted metadata-rounding error.
        with self.assertRaisesRegex(ValueError, 'Two-spinon tolerance changed'):
            tutorial.read_case(document, 4, 8)

    def test_refinement_precision_evidence_and_failures(self):
        path = tutorial.ROOT / 'docs/validation/four-spinon-refinement-precision-20261031.json'
        study = json.loads(path.read_text())
        expected = {(p, n, t) for p in ('long-double', 'bits-128')
                    for n in tutorial.ORDERS for t in tutorial.CUTOFFS}
        expected |= {('fp128', 8, t) for t in tutorial.CUTOFFS}
        actual = [(r['precision'], r['order'], r['cutoff']) for r in study['cases']]
        self.assertEqual(set(actual), expected)
        self.assertEqual(len(actual), len(expected))
        references = {(r['order'], r['cutoff']): r for r in study['reference_192']['cases']}
        grid = {(r['order'], r['cutoff']): r for r in self.data['cases']}
        differences = []
        with localcontext() as context:
            context.prec = 100
            for row in study['cases']:
                key = row['order'], row['cutoff']
                failed = row['precision'] != 'bits-128' and row['order'] >= 8
                self.assertEqual(row['exit_code'], 2 if failed else 0)
                self.assertRegex(row['raw_sha256'], re.compile('[0-9a-f]{64}'))
                self.assertGreater(Decimal(row['cpu_seconds']), 0)
                if failed:
                    self.assertEqual(row['status'], 'numerical_failure')
                    self.assertIsNone(row['s4'])
                    self.assertIsNone(row['partial_sum'])
                    self.assertIn('arithmetic screen failed', row['reason'])
                else:
                    self.assertEqual(row['status'], 'estimated')
                    self.assertGreater(Decimal(row['s4']), 0)
                    reference = references[key]
                    difference = abs(Decimal(row['s4']) / Decimal(reference['s4']) - 1)
                    self.assertLess(difference, Decimal('1e-14'))
                    if row['precision'] == 'bits-128':
                        self.assertEqual(row['s4'], grid[key]['s4'])
                        self.assertEqual(row['raw_sha256'], grid[key]['raw_sha256'])
                        self.assertEqual(row['density_attempts'], reference['density_attempts'])
                        self.assertLess(difference, Decimal('1e-35'))
                        differences.append(difference)
            self.assertLess(abs(max(differences) / Decimal(study['maximum_s4_relative_difference_128_vs_192']) - 1), Decimal('1e-60'))


if __name__ == '__main__':
    unittest.main()
