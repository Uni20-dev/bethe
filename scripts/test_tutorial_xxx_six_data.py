"""Compact slices preserve stream evidence, missing nodes and tiny weights."""
from copy import deepcopy
from decimal import Decimal
import hashlib
import json
import unittest

from test_tutorial_xxx_six_spinon import fixture
from xxx_six_spinon_tutorial_data import compact_slice, read_compact_slice


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


if __name__ == '__main__':
    unittest.main()
