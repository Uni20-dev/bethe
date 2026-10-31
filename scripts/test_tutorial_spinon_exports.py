"""Prevent plausible plots assembled from incompatible or incomplete exports."""
from copy import deepcopy
import unittest

from spinon_tutorial_exports import CONVENTION, read_spectrum


def fixture():
    metadata = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
                'Density convention': CONVENTION, 'Geometry and output bits': '192',
                'Exchange J': '1', 'Selected spinon sectors': '2,4', 'Spectral points': '1',
                'Bethe revision': 'a' * 40, 'Uni20 revision': 'b' * 40}
    return {'tables': {'spectrum': {'metadata': metadata,
            'summary': {'Outcome': 'success', 'Failed points': '0'},
            'columns': [{'id': x} for x in ('q', 'omega', 'sector', 'density', 'status')],
            'rows': [['1.5', '2', '2', '0.5', 'converged'],
                     ['1.5', '2', '4', '0.125', 'estimated'],
                     ['1.5', '2', 'sum:2+4', '0.625', 'partial_sum']]}}}


class SpinonExports(unittest.TestCase):
    def test_decimal_strings_are_retained(self):
        result = read_spectrum(fixture(), momentum='1.5')
        self.assertEqual(result['spectrum'][1]['density'], '0.125')

    def test_incompatible_or_unidentifiable_runs_are_rejected(self):
        for key, value in [('Channel', 'raising'), ('Exchange J', '2'),
                           ('Bethe revision', 'a' * 40 + '-dirty'),
                           ('Selected spinon sectors', '2,4,6')]:
            with self.subTest(key=key):
                data = fixture()
                data['tables']['spectrum']['metadata'][key] = value
                with self.assertRaises(ValueError):
                    read_spectrum(data, momentum='1.5')
        with self.assertRaisesRegex(ValueError, 'momentum'):
            read_spectrum(fixture(), momentum='1.6')

    def test_missing_failed_or_rescaled_weight_is_not_hidden(self):
        for kind in ('missing', 'failed', 'sum', 'nonfinite', 'duplicate'):
            with self.subTest(kind=kind):
                data = fixture()
                rows = data['tables']['spectrum']['rows']
                if kind == 'missing':
                    rows.pop(1)
                elif kind == 'failed':
                    rows[1][3:] = [None, 'numerical_failure']
                elif kind == 'sum':
                    rows[-1][3] = '0.7'
                elif kind == 'nonfinite':
                    rows[1][3] = 'NaN'
                else:
                    rows.append(deepcopy(rows[1]))
                with self.assertRaises(ValueError):
                    read_spectrum(data, momentum='1.5')

    def test_all_tables_must_share_provenance(self):
        data = fixture()
        other = deepcopy(data['tables']['spectrum'])
        other['metadata']['Exchange J'] = '2'
        data['tables']['diagnostics'] = other
        with self.assertRaisesRegex(ValueError, 'provenance'):
            read_spectrum(data, momentum='1.5')

    def test_faint_partial_sum_cannot_lose_the_higher_sector(self):
        data = fixture()
        rows = data['tables']['spectrum']['rows']
        rows[0][3:] = ['0', 'physical_zero']
        rows[1][3] = rows[2][3] = '1e-100'
        read_spectrum(data, momentum='1.5')
        rows[2][3] = '0'
        with self.assertRaisesRegex(ValueError, 'Partial sum'):
            read_spectrum(data, momentum='1.5')

    def test_typed_precision_must_match(self):
        data = fixture()
        data['tables']['spectrum']['rows'][0][0] = {'value': '1.5', 'precision': 256}
        with self.assertRaisesRegex(ValueError, 'precision'):
            read_spectrum(data, momentum='1.5')


if __name__ == '__main__':
    unittest.main()
