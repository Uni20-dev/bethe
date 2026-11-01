"""Prevent plausible plots assembled from incompatible or incomplete exports."""
from copy import deepcopy
import unittest

from spinon_tutorial_exports import CONVENTION, output_precision_bits, plotting_float, read_spectrum


def fixture():
    metadata = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
                'Hamiltonian': 'H=J sum S_j.S_(j+1); infinite chain; zero field and temperature',
                'Family': 'Partial DSF; physical normalization; no fitted rescaling',
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
    def test_plot_conversion_preserves_zeros_and_rejects_lost_range(self):
        self.assertEqual(plotting_float('0'), 0)
        self.assertEqual(plotting_float('-0'), 0)
        self.assertGreater(plotting_float('1e-310'), 0)
        self.assertEqual(plotting_float('1.25'), 1.25)
        for value in ('1e-1000', '-1e-1000', '1e1000', '-1e1000', 'NaN', 'Infinity'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                plotting_float(value)

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

    def test_mixed_scalars_use_output_precision_not_an_earlier_stage(self):
        data = fixture()
        table = data['tables']['spectrum']
        meta = table['metadata']
        meta['Accumulation and output bits'] = meta.pop('Geometry and output bits')
        meta['Six-spinon precision stages'] = 'fp80,fp128,reference'
        meta['Selected spinon sectors'] = '6'
        table['rows'] = [['1.5', '2', '6', {'value': '0.125', 'precision': 192}, 'estimated']]
        self.assertEqual(read_spectrum(data, momentum='1.5')['spectrum'][0]['density'], '0.125')
        table['rows'][0][3]['precision'] = 64
        with self.assertRaisesRegex(ValueError, 'precision'):
            read_spectrum(data, momentum='1.5')

    def test_mixed_precision_cannot_masquerade_as_fixed_geometry(self):
        good = {'Accumulation and output bits': '192', 'Six-spinon precision stages': 'fp80,fp128,reference',
                'Selected spinon sectors': '6'}
        self.assertEqual(output_precision_bits(good), 192)
        for changes in ({'Geometry and output bits': '192'}, {'Six-spinon precision stages': ''},
                        {'Six-spinon precision stages': 'fp128,fp80,reference'},
                        {'Six-spinon precision stages': 'fp80,fp80,reference'},
                        {'Six-spinon precision stages': 'fp80'}, {'Selected spinon sectors': '2,4'},
                        {'Accumulation and output bits': True}, {'Accumulation and output bits': '192.5'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                output_precision_bits({**good, **changes})
        for key in ('Six-spinon precision stages', 'Accumulation and output bits'):
            bad = dict(good)
            del bad[key]
            with self.subTest(missing=key), self.assertRaises(ValueError):
                output_precision_bits(bad)
        bad = dict(good)
        bad['Geometry and output bits'] = bad.pop('Accumulation and output bits')
        with self.assertRaisesRegex(ValueError, 'geometry'):
            output_precision_bits(bad)


if __name__ == '__main__':
    unittest.main()
