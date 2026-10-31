"""The native comparison keeps the difficult seeds and complete observations."""
from copy import deepcopy
import json
import unittest

import check_xxx_six_spinon_native_precision as check


class SixNativeComparison(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(check.DATA.read_text())

    def test_all_three_seeds_reconstruct_against_mpfr(self):
        self.assertEqual(set(check.validate(self.data)), {'42', '2', '35'})

    def test_difficult_seed_cannot_be_dropped(self):
        self.data['cases'].pop()
        with self.assertRaisesRegex(ValueError, 'Missing case'):
            check.validate(self.data)

    def test_looser_screen_or_weaker_pair_is_rejected(self):
        for field, value in [('Density relative tolerance', '1e-3'), ('Density comparison bits', '53'),
                             ('Six-spinon seed', '43')]:
            data = deepcopy(self.data)
            data['cases'][0]['native']['runs'][0]['metadata'][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                check.validate(data)

    def test_changed_replicate_maximum_is_not_hidden_by_the_mean(self):
        self.data['cases'][0]['native']['runs'][0]['replicates'][0]['maximum_absolute_density'] = '0.003'
        with self.assertRaises(ValueError):
            check.validate(self.data)


if __name__ == '__main__':
    unittest.main()
