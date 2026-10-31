"""Protect the neighbourhood audit's failures, coverage and reference comparisons."""
from copy import deepcopy
import json
import unittest

import check_xxx_four_spinon_neighbourhood as study


class FourSpinonNeighbourhood(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(study.DATA.read_text())

    def test_complete_audit_includes_coverage_and_wide_references(self):
        self.assertIs(study.validate(self.data), self.data)
        for tol, count in [('1e-8', 26), ('1e-6', 47), ('1e-4', 49)]:
            selected = [c for c in self.data['cases'][:324]
                        if c['precision'] == 'long-double' and c['relative_tolerance'] == tol]
            self.assertEqual(sum(c['status'] == 'estimated' for c in selected), count)

    def test_failed_points_cannot_be_replaced_or_dropped(self):
        for change in ('density', 'status', 'missing'):
            with self.subTest(change=change):
                bad = deepcopy(self.data)
                failed = next(c for c in bad['cases'] if c['status'] == 'numerical_failure')
                if change == 'density':
                    failed['density'] = '0'
                elif change == 'status':
                    failed.update(status='physical_zero', density='0', evaluations='0', exit_code=0,
                                  outcome='success', failed_points='0')
                else:
                    bad['cases'].remove(failed)
                with self.assertRaises(ValueError):
                    study.validate(bad)

    def test_success_is_checked_against_matching_reference(self):
        bad = deepcopy(self.data)
        row = next(c for c in bad['cases'] if c['precision'] == 'long-double'
                   and c['relative_tolerance'] == '1e-6' and c['q'] == '1.5' and c['omega'] == '3.0')
        row['density'] = '0.003'
        with self.assertRaisesRegex(ValueError, 'disagreement'):
            study.validate(bad)

    def test_precision_and_controls_are_not_silently_changed(self):
        for change in ('comparison', 'coordinates', 'encoding', 'cutoff', 'cost'):
            with self.subTest(change=change):
                bad = deepcopy(self.data)
                if change == 'comparison':
                    bad['cases'][0]['comparison_bits'] = 53
                elif change == 'coordinates':
                    bad['cases'][0]['computed_omega'] = '2'
                elif change == 'encoding':
                    bad['cases'][0]['computed_omega'] = 1.0
                elif change == 'cutoff':
                    bad['controls']['Four-spinon log-radius cutoff'] = '12'
                else:
                    bad['cases'][0]['cpu_seconds'] = 'NaN'
                with self.assertRaises(ValueError):
                    study.validate(bad)


if __name__ == '__main__':
    unittest.main()
