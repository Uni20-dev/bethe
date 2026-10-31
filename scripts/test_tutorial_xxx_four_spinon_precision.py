"""Check that the precision study preserves controls, failures and arithmetic evidence."""
from copy import deepcopy
from decimal import Decimal, localcontext
import json
import unittest

import check_xxx_four_spinon_precision as study


class FourSpinonPrecision(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(study.DATA.read_text())

    def test_complete_matrix_includes_failure_and_agreeing_successes(self):
        self.assertIs(study.validate(self.data), self.data)
        self.assertEqual(self.data['cases'][0]['spectrum'][1]['status'], 'numerical_failure')
        self.assertIsNone(self.data['cases'][0]['spectrum'][2]['density'])
        native, reference = [self.data['cases'][i]['spectrum'][1]['density'] for i in (1, 4)]
        self.assertTrue(study.close(native, reference, Decimal('2e-16')))

    def test_controls_and_effective_precision_cannot_change_silently(self):
        for key, value in [('Channel', 'raising'), ('Density comparison bits', '52'),
                           ('Geometry and output bits', '53'), ('Density relative tolerance', '1e-5'),
                           ('Four-spinon log-radius cutoff', '4')]:
            with self.subTest(key=key):
                bad = deepcopy(self.data)
                bad['cases'][1]['metadata'][key] = value
                with self.assertRaises(ValueError):
                    study.validate(bad)

    def test_failure_is_not_zero_or_a_successful_prefix(self):
        for alteration in ('zero', 'reason', 'success', 'missing'):
            with self.subTest(alteration=alteration):
                bad = deepcopy(self.data)
                failed = bad['cases'][0]
                if alteration == 'zero':
                    failed['spectrum'][1]['density'] = '0'
                elif alteration == 'reason':
                    failed['diagnostics'][1]['reason'] = ''
                elif alteration == 'success':
                    failed['exit_code'] = 0
                else:
                    bad['cases'].pop(0)
                with self.assertRaises(ValueError):
                    study.validate(bad)

    def test_agreement_is_checked_beyond_a_valid_partial_sum(self):
        bad = deepcopy(self.data)
        rows = bad['cases'][1]['spectrum']
        with localcontext() as context:
            context.prec = 120
            rows[1]['density'] = str(Decimal(rows[1]['density'])*Decimal('1.001'))
            rows[2]['density'] = str(Decimal(rows[0]['density'])+Decimal(rows[1]['density']))
        with self.assertRaisesRegex(ValueError, 'disagreement'):
            study.validate(bad)

    def test_reproducibility_and_cost_are_retained(self):
        for alteration in ('command', 'source', 'cost'):
            with self.subTest(alteration=alteration):
                bad = deepcopy(self.data)
                row = bad['cases'][1]
                if alteration == 'command':
                    row['command'][-1] = 'bits-192'
                elif alteration == 'source':
                    row['metadata']['Bethe revision'] = '0'*40
                else:
                    row['summary']['Run CPU seconds'] = 'NaN'
                with self.assertRaises(ValueError):
                    study.validate(bad)


if __name__ == '__main__':
    unittest.main()
