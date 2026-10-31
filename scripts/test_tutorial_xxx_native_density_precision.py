"""The native comparison audit keeps complete observations and failed queries."""
from copy import deepcopy
import json
import unittest

import check_xxx_native_density_precision as study


class NativeDensityPrecision(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(study.DATA.read_text())

    def test_retained_matrix_matches_wider_references(self):
        self.assertEqual(study.comparisons(self.data), self.data['maximum_relative_differences'])
        self.assertEqual(len(self.data['cases']), 20)
        self.assertEqual([c['name'] for c in self.data['cases'] if c['exit_code'] == 2],
                         ['spectrum-fp64', 'moments-fp64'])

    def test_weaker_pair_or_changed_controls_are_rejected(self):
        for key, value in [('Density comparison bits', '53'), ('Density relative tolerance', '1e-5'),
                           ('Six-spinon seed', '43')]:
            bad = deepcopy(self.data)
            bad['cases'][2]['metadata'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                study.comparisons(bad)

    def test_missing_query_or_failed_value_cannot_be_hidden(self):
        bad = deepcopy(self.data)
        bad['cases'].pop(0)
        with self.assertRaises(ValueError):
            study.comparisons(bad)
        bad = deepcopy(self.data)
        bad['cases'][0]['tables']['spectrum'][0]['density'] = '0'
        with self.assertRaises(ValueError):
            study.comparisons(bad)

    def test_covariance_and_replicates_are_compared(self):
        for table, column in [('covariance', 'covariance'), ('replicates', 'maximum_m0')]:
            bad = deepcopy(self.data)
            bad['cases'][7]['tables'][table][0][column] = '0.1'
            with self.subTest(table=table), self.assertRaisesRegex(ValueError, 'Disagreement'):
                study.comparisons(bad)
        bad = deepcopy(self.data)
        bad['cases'][7]['tables']['replicates'][0]['completed_samples'] = '3'
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            study.comparisons(bad)

    def test_native_grid_checks_the_wider_reference(self):
        bad = deepcopy(self.data)
        bad['cases'][-1]['tables']['spectrum'][1]['density'] = '0.2'
        with self.assertRaisesRegex(ValueError, 'MPFR192'):
            study.comparisons(bad)


if __name__ == '__main__':
    unittest.main()
