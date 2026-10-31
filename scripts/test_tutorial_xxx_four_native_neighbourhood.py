"""A native grid must retain the difficult point and its matched reference."""
from copy import deepcopy
import json
import unittest

import check_xxx_four_spinon_native_neighbourhood as audit


class NativeFourNeighbourhood(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(audit.DATA.read_text())
        self.original = json.loads(audit.ORIGINAL.read_text())

    def test_complete_grid_reconstructs_the_coverage_table(self):
        strict, ordinary = audit.validate(self.data, self.original)
        self.assertEqual(strict['complete'], 48)
        self.assertEqual(strict['failed_points'], [['1.5', '4.5']])
        self.assertEqual(ordinary['complete'], 49)
        self.assertEqual(ordinary['failed_points'], [])

    def test_failed_point_cannot_disappear_or_become_a_zero(self):
        failed = next(c for c in self.data['cases'] if c['status'] == 'numerical_failure')
        for change in ('drop', 'zero'):
            bad = deepcopy(self.data)
            row = next(c for c in bad['cases'] if c['raw_sha256'] == failed['raw_sha256'])
            if change == 'drop':
                bad['cases'].remove(row)
            else:
                row.update(status='physical_zero', density='0', exit_code=0)
            with self.subTest(change=change), self.assertRaises(ValueError):
                audit.validate(bad, self.original)

    def test_small_nonzero_density_is_compared_relatively(self):
        row = min((c for c in self.data['cases'] if c['status'] == 'estimated'),
                  key=lambda c: audit.number(c['density']))
        row['density'] = str(audit.number(row['density']) * audit.Decimal('1.000001'))
        with self.assertRaisesRegex(ValueError, 'value disagreement'):
            audit.validate(self.data, self.original)

    def test_reference_identity_and_precision_pair_are_checked(self):
        for change in ('reference', 'precision', 'support', 'summary'):
            bad = deepcopy(self.data)
            if change == 'reference':
                bad['cases'][0]['reference_sha256'] = bad['references'][1]['raw_sha256']
            elif change == 'precision':
                bad['native_metadata']['Density comparison bits'] = '53'
            elif change == 'support':
                bad['excluded_physical_zero_points'].pop()
            else:
                bad['summaries'][0]['complete'] = 49
            with self.subTest(change=change), self.assertRaises(ValueError):
                audit.validate(bad, self.original)


if __name__ == '__main__':
    unittest.main()
