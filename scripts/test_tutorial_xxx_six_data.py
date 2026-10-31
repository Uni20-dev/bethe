"""Compact slices preserve stream evidence, missing nodes and tiny weights."""
from copy import deepcopy
from decimal import Decimal
import hashlib
import json
import unittest

from test_tutorial_xxx_six_spinon import fixture
from xxx_six_spinon_tutorial_data import compact_slice, match_requested_frequency, read_compact_slice


class CompactSixSlice(unittest.TestCase):
    def packed(self, data=None, **options):
        raw = json.dumps(fixture() if data is None else data).encode()
        return compact_slice([raw], momentum='1.5', rapidity_cutoff='64',
                             frequencies=options.pop('frequencies', ['2']), **options)

    def test_roundtrip_preserves_exact_evidence_and_recomputes_estimates(self):
        data = fixture()
        for table in data['tables'].values():
            table['metadata']['Command'] = 'program --json /local/raw.json'
        packed = self.packed(data)
        result = read_compact_slice(json.loads(json.dumps(packed)), momentum='1.5')
        point = result['points'][0]
        self.assertEqual(point['s6'], '2')
        self.assertEqual(point['sampling_se'], '1')
        self.assertEqual(len(point['replicate_records']), 2)
        self.assertEqual(point['raw_sha256'], hashlib.sha256(json.dumps(data).encode()).hexdigest())
        self.assertNotIn('Command', point['metadata'])
        packed['runs'][0]['diagnostics'][0]['sampling_error'] = '0'
        with self.assertRaisesRegex(ValueError, 'mean/SE'):
            read_compact_slice(packed, momentum='1.5')

    def test_missing_node_requires_explicit_null_not_zero_or_removal(self):
        with self.assertRaisesRegex(ValueError, 'Missing or undeclared'):
            self.packed(frequencies=['2', '3'])
        unavailable = [{'omega': '3', 'status': 'failed', 'reason': 'density recovery exhausted',
                        'raw_sha256': 'c' * 64}]
        packed = self.packed(frequencies=['2', '3'], unavailable=unavailable)
        points = read_compact_slice(packed, momentum='1.5')['points']
        self.assertEqual([point['omega'] for point in points], ['2', '3'])
        self.assertIsNone(points[1]['s6'])
        self.assertIsNone(points[1]['sampling_se'])
        self.assertEqual(points[1]['status'], 'failed')
        self.assertEqual(points[1]['raw_sha256'], 'c' * 64)

    def test_duplicate_and_malformed_grids_are_rejected(self):
        for grid in (['2', '2.0'], ['3', '2'], [], ['2', '4'], [2]):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                self.packed(frequencies=grid)
        packed = self.packed()
        packed['runs'].append(deepcopy(packed['runs'][0]))
        with self.assertRaisesRegex(ValueError, 'Duplicate frequency'):
            read_compact_slice(packed, momentum='1.5')

    def test_wrong_identity_cutoff_or_hash_is_rejected(self):
        for field, value in [('schema', 'future'), ('momentum', '1.6'),
                             ('exchange', '2'), ('rapidity_cutoff', '128'),
                             ('density_convention', 'unit-normalized')]:
            with self.subTest(field=field):
                packed = self.packed()
                packed[field] = value
                with self.assertRaises(ValueError):
                    read_compact_slice(packed, momentum='1.5')
        packed = self.packed()
        packed['runs'][0]['raw_sha256'] = 'missing'
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            read_compact_slice(packed, momentum='1.5')

    def test_faint_weight_stays_decimal_and_zeroed_error_is_rejected(self):
        data = fixture()
        for row in data['tables']['replicates']['rows']:
            for column in (6, 7):
                row[column] = str(Decimal(row[column]) * Decimal('1e-1000'))
        data['tables']['spectrum']['rows'][0][3] = '2e-1000'
        data['tables']['diagnostics']['rows'][0][3] = '1e-1000'
        packed = self.packed(data)
        point = read_compact_slice(packed, momentum='1.5')['points'][0]
        self.assertEqual(Decimal(point['s6']), Decimal('2e-1000'))
        self.assertEqual(Decimal(point['sampling_se']), Decimal('1e-1000'))
        packed['runs'][0]['diagnostics'][0]['sampling_error'] = '0'
        with self.assertRaisesRegex(ValueError, 'mean/SE'):
            read_compact_slice(packed, momentum='1.5')

    def test_separate_frequency_runs_keep_controls_and_reject_mixed_sources(self):
        first, second = fixture(), fixture()
        for table in second['tables'].values():
            table['metadata']['Six-spinon seed'] = '43'
            for row in table['rows']:
                row[1] = '3'
        packed = compact_slice([json.dumps(x).encode() for x in (first, second)],
                               momentum='1.5', rapidity_cutoff='64', frequencies=['2', '3'])
        points = read_compact_slice(packed, momentum='1.5')['points']
        self.assertEqual([point['metadata']['Six-spinon seed'] for point in points], ['42', '43'])
        self.assertEqual([point['run_index'] for point in points], [0, 1])
        packed['runs'][1]['metadata']['Bethe revision'] = 'd' * 40
        with self.assertRaisesRegex(ValueError, 'Mixed source revisions'):
            read_compact_slice(packed, momentum='1.5')

    def test_declared_decimal_node_preserves_rounded_export_coordinate(self):
        # Actual binary80 and binary128 serializations of decimal 1.8.
        for bits, actual in [(64, '1.79999999999999999996'),
                             (113, '1.80000000000000000000000000000000004')]:
            with self.subTest(bits=bits):
                data = fixture()
                for table in data['tables'].values():
                    table['metadata']['Geometry and output bits'] = str(bits)
                    for row in table['rows']:
                        row[1] = actual
                packed = self.packed(data, frequencies=['1.8'])
                point = read_compact_slice(packed, momentum='1.5')['points'][0]
                self.assertEqual(point['omega'], '1.8')
                self.assertEqual(point['computed_omega'], actual)
                self.assertEqual(point['replicate_records'][0]['omega'], actual)
                self.assertEqual(packed['runs'][0]['spectrum'][0]['omega'], actual)

    def test_rounding_match_cannot_interpolate_or_collapse_grid(self):
        self.assertEqual(match_requested_frequency('1.79999999999999999996', ['1.8'], 64), Decimal('1.8'))
        for actual, requested, bits, message in [
                ('1.81', ['1.8'], 64, 'Undeclared'),
                ('1.8', ['1.8', '1.800000000000000001'], 64, 'Ambiguous'),
                ('1.8', ['1.8', '1.80'], 64, 'Duplicate'),
                ('1e-100', ['0'], 64, 'Undeclared'),
                ('1e-100', ['2e-100'], 64, 'Undeclared')]:
            with self.subTest(actual=actual, requested=requested), self.assertRaisesRegex(ValueError, message):
                match_requested_frequency(actual, requested, bits)
        data = self.packed()
        data['frequencies'] = ['2', '2.000000000000000000000000000000000000000000000000000000001']
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            read_compact_slice(data, momentum='1.5')
        first, second = fixture(), fixture()
        for document, actual in [(first, '1.79999999999999999996'), (second, '1.8')]:
            for table in document['tables'].values():
                table['metadata']['Geometry and output bits'] = '64'
                for row in table['rows']:
                    row[1] = actual
        with self.assertRaisesRegex(ValueError, 'Duplicate frequency'):
            compact_slice([json.dumps(d).encode() for d in (first, second)], momentum='1.5',
                          rapidity_cutoff='64', frequencies=['1.8'])


if __name__ == '__main__':
    unittest.main()
