"""Do not lose failed work or accept survivor means in mixed-precision exports."""
from copy import deepcopy
from decimal import Decimal
import json
import unittest

from plot_xxx_six_spinon_tutorial import paired_cutoffs, sampling_points
from spinon_recovery_exports import COLUMNS, recovery_policy
from test_tutorial_xxx_six_spinon import fixture as fixed_fixture
from xxx_six_spinon_tutorial_data import (RECOVERY_SCHEMA, SCHEMA, compact_slice,
                                         read_compact_slice)


def mixed_fixture():
    data = fixed_fixture()
    metadata = dict(data['tables']['spectrum']['metadata'])
    metadata['Accumulation and output bits'] = metadata.pop('Geometry and output bits')
    metadata.update({'Six-spinon precision stages': 'fp80,fp128,reference',
                     'Reference density comparison bits': '128', 'Reference density guard bits': '192',
                     'Density relative tolerance': '1e-6', 'Density absolute tolerance': '0',
                     'First density residue budget': '32', 'Second density residue budget': '48',
                     'First density tail order': '12', 'Second density tail order': '16',
                     'Recovery seconds per configuration': '60'})
    for name, fmt, geometry, low, high in (('fp80', 'extended80', 64, 53, 64),
                                         ('fp128', 'binary128', 113, 64, 113),
                                         ('reference', 'arbitrary binary', 192, 128, 192)):
        metadata.update({f'Six-spinon {name} arithmetic format': fmt,
                         f'Six-spinon {name} geometry bits': str(geometry),
                         f'Six-spinon {name} density comparison bits': str(low),
                         f'Six-spinon {name} density guard bits': str(high)})
    # PR92's real fixtures also supply generic fixed-precision density fields.
    metadata.pop('Density comparison bits', None)
    metadata.pop('Density guard bits', None)
    columns = ['q', 'omega', 'replicate', 'stage', 'kind', 'attempts', 'accepted', 'retried',
               'supported_roots', 'density_attempts', 'failure', 'density_status', 'density_failure', 'reason']
    def stage(replica, name, attempts, accepted, retried, roots, density):
        return ['1.5', '2', str(replica), name, 'stage', str(attempts), str(accepted), str(retried),
                str(roots), str(density), '', '', '', '']
    def failure(replica, name, kind, roots=0, density=0, status='', density_failure=''):
        return ['1.5', '2', str(replica), name, 'failure', '1', None, None, str(roots), str(density),
                kind, status, density_failure, 'controlled ' + kind]
    data['tables']['recovery'] = {
        'metadata': {}, 'summary': dict(data['tables']['spectrum']['summary']),
        'columns': [{'id': name} for name in columns],
        'rows': [stage(0, 'fp80', 4, 2, 2, 2, 4), failure(0, 'fp80', 'geometry'),
                 failure(0, 'fp80', 'exponent_range'), stage(0, 'fp128', 2, 1, 1, 1, 1),
                 failure(0, 'fp128', 'density', 1, 1, 'exhausted', 'arithmetic_precision'),
                 stage(0, 'reference', 1, 1, 0, 0, 0),
                 stage(1, 'fp80', 4, 3, 1, 2, 4), failure(1, 'fp80', 'geometry'),
                 stage(1, 'fp128', 1, 1, 0, 2, 4), stage(1, 'reference', 0, 0, 0, 0, 0)]}
    for table in data['tables'].values():
        table['metadata'] = dict(metadata)
    return data


def change_metadata(data, key, value):
    for table in data['tables'].values():
        if value is None:
            del table['metadata'][key]
        else:
            table['metadata'][key] = value


def change_row(data, index, key, value):
    table = data['tables']['recovery']
    column = [c['id'] for c in table['columns']].index(key)
    table['rows'][index][column] = value


def compact(data):
    return compact_slice([json.dumps(data).encode()], momentum='1.5', frequencies=['2'], rapidity_cutoff='64')


class SixSpinonRecoveryExports(unittest.TestCase):
    def test_physical_values_and_all_failed_work_are_retained(self):
        data = mixed_fixture()
        result = sampling_points(data, momentum='1.5')
        point, = result['points']
        self.assertEqual(Decimal(point['s6']), 2)
        self.assertEqual(Decimal(point['sampling_se']), 1)
        self.assertEqual(point['supported_roots'], 6)
        self.assertEqual(point['density_attempts'], 12)
        totals = result['recovery_totals']
        self.assertEqual([totals[s]['attempts'] for s in totals], [8, 3, 1])
        self.assertEqual([totals[s]['accepted'] for s in totals], [5, 2, 1])
        self.assertEqual(sum(t['supported_roots'] for t in totals.values()), 7)
        self.assertEqual(sum(t['density_attempts'] for t in totals.values()), 13)
        self.assertEqual(len(point['recovery_records']), 10)
        self.assertEqual(set(point['recovery_records'][0]), COLUMNS)

    def test_compact_v2_roundtrip_and_strict_v1_preservation(self):
        fixed = compact(fixed_fixture())
        self.assertEqual(fixed['schema'], SCHEMA)
        self.assertNotIn('recovery', fixed['runs'][0])
        packed = compact(mixed_fixture())
        self.assertEqual(packed['schema'], RECOVERY_SCHEMA)
        restored = read_compact_slice(json.loads(json.dumps(packed)), momentum='1.5')
        self.assertEqual(restored['points'][0]['recovery_records'], packed['runs'][0]['recovery'])
        for case in ('old-schema', 'drop-table', 'drop-zero-stage', 'fake-fixed'):
            bad = deepcopy(packed)
            if case == 'old-schema': bad['schema'] = SCHEMA
            elif case == 'drop-table': del bad['runs'][0]['recovery']
            elif case == 'drop-zero-stage': bad['runs'][0]['recovery'].pop()
            else: bad['runs'][0]['metadata'].pop('Six-spinon precision stages')
            with self.subTest(case=case), self.assertRaises(ValueError):
                read_compact_slice(bad, momentum='1.5')
        fixed['runs'][0]['recovery'] = []
        with self.assertRaises(ValueError):
            read_compact_slice(fixed, momentum='1.5')

    def test_missing_or_duplicate_records_never_disappear(self):
        for case in ('missing-table', 'missing-unused-stage', 'duplicate-stage', 'duplicate-failure'):
            data = mixed_fixture()
            if case == 'missing-table': del data['tables']['recovery']
            elif case == 'missing-unused-stage': data['tables']['recovery']['rows'].pop()
            else: data['tables']['recovery']['rows'].append(deepcopy(data['tables']['recovery']['rows'][case == 'duplicate-failure']))
            with self.subTest(case=case), self.assertRaises(ValueError):
                sampling_points(data, momentum='1.5')

    def test_corrupt_stage_identity_counts_and_failure_work_are_rejected(self):
        changes = [(0, 'q', '1.6'), (0, 'omega', '3'), (0, 'replicate', '9'), (0, 'stage', 'fp64'),
                   (0, 'kind', 'unknown'), (0, 'attempts', '3.5'), (0, 'accepted', None),
                   (0, 'retried', '1'), (0, 'failure', 'geometry'), (0, 'density_status', None), (0, 'reason', 'invented'),
                   (0, 'supported_roots', '1'), (0, 'density_attempts', '3'),
                   (1, 'accepted', '0'), (1, 'retried', '0'), (1, 'attempts', '0'),
                   (1, 'reason', ''), (1, 'failure', 'cancelled'), (1, 'density_status', 'exhausted'), (1, 'density_failure', None),
                   (1, 'supported_roots', '3'), (1, 'density_attempts', '1'),
                   (4, 'density_status', 'time_limit'), (4, 'density_failure', 'unknown'),
                   (5, 'accepted', '0'), (9, 'attempts', '1')]
        for index, key, value in changes:
            data = mixed_fixture()
            change_row(data, index, key, value)
            with self.subTest(row=index, key=key, value=value), self.assertRaises(ValueError):
                sampling_points(data, momentum='1.5')

    def test_retry_chain_cannot_lose_a_sample_even_when_local_counts_balance(self):
        data = mixed_fixture()
        change_row(data, 3, 'attempts', '1')
        change_row(data, 3, 'accepted', '0')
        with self.assertRaisesRegex(ValueError, 'retry chain'):
            sampling_points(data, momentum='1.5')
        data = mixed_fixture()
        data['tables']['replicates']['rows'][0][4] = '4'
        with self.assertRaisesRegex(ValueError, 'physical replicate'):
            sampling_points(data, momentum='1.5')

    def test_malformed_or_mislabelled_precision_policy_is_rejected(self):
        changes = [('Six-spinon fp80 arithmetic format', 'binary64'),
                   ('Six-spinon reference arithmetic format', 'binary128'),
                   ('Six-spinon fp80 geometry bits', '80'), ('Six-spinon fp80 density comparison bits', '64'),
                   ('Six-spinon fp128 density guard bits', '64'), ('Six-spinon uniform grid bits', '64'),
                   ('Six-spinon reference geometry bits', '193'), ('Six-spinon reference density guard bits', '100'),
                   ('Reference density guard bits', '113'), ('Density guard bits', '192'),
                   ('Six-spinon fp128 arithmetic format', None), ('Two/four-spinon geometry bits', '192'),
                   ('Density relative tolerance', '-1'), ('First density tail order', '17'),
                   ('Recovery seconds per configuration', '0')]
        for key, value in changes:
            data = mixed_fixture()
            change_metadata(data, key, value)
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                sampling_points(data, momentum='1.5')
        # Native binary128 can compare with binary64 when long double is not wider.
        data = mixed_fixture()
        change_metadata(data, 'Six-spinon fp128 density comparison bits', '53')
        sampling_points(data, momentum='1.5')

    def test_reference_only_policy_and_zero_attempt_native_stages(self):
        data = mixed_fixture()
        meta = data['tables']['spectrum']['metadata']
        change_metadata(data, 'Six-spinon precision stages', 'reference')
        for key in list(meta):
            if key.startswith(('Six-spinon fp80 ', 'Six-spinon fp128 ')):
                change_metadata(data, key, None)
        rows = data['tables']['recovery']['rows']
        rows[:] = [rows[5], rows[9]]
        for row, roots in zip(rows, (2, 4)):
            row[5:10] = ['4', '4', '0', str(roots), str(2 * roots)]
        result = sampling_points(data, momentum='1.5')
        self.assertEqual(result['recovery_totals']['reference']['accepted'], 8)

    def test_combined_sector_reader_keeps_the_full_recovery_records(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path
        from test_tutorial_spinon_exports import fixture as four_fixture
        import plot_xxx_spinon_slice_tutorial as combined
        four = four_fixture()
        four['tables']['spectrum']['metadata']['Command'] = 'bethe-xxx-structure-factor-thermo --sectors 2,4'
        six = compact(mixed_fixture())
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'four.json'
            path.write_text(json.dumps(four))
            result = combined.collect([path], six, momentum='1.5')
        self.assertEqual(result['six_slice'], six)
        self.assertEqual(Decimal(result['points'][0]['s2_plus_s4_plus_s6']), Decimal('2.625'))
        self.assertEqual(combined.validate(json.loads(json.dumps(result))), result)
        result['six_slice']['runs'][0]['recovery'].pop()
        with self.assertRaises(ValueError): combined.validate(result)

    def test_reference_precision_is_also_explicit_for_two_four_sectors(self):
        metadata = mixed_fixture()['tables']['spectrum']['metadata']
        metadata['Selected spinon sectors'] = '2,4,6'
        with self.assertRaisesRegex(ValueError, 'Two/four'): recovery_policy(metadata)
        metadata['Two/four-spinon geometry bits'] = '192'
        self.assertEqual(len(recovery_policy(metadata)), 3)
        metadata['Two/four-spinon geometry bits'] = '64'
        with self.assertRaisesRegex(ValueError, 'Two/four'): recovery_policy(metadata)

    def test_failed_or_cancelled_runs_cannot_become_successful_compact_data(self):
        for case in ('cancelled', 'failed-sample'):
            data = mixed_fixture()
            if case == 'cancelled':
                for table in data['tables'].values(): table['summary']['Outcome'] = 'cancelled'
            else: data['tables']['replicates']['rows'][0][-1] = '3'
            with self.subTest(case=case), self.assertRaises(ValueError): compact(data)

    def test_paired_cutoffs_bind_all_stage_and_reference_controls(self):
        lower, upper = mixed_fixture(), mixed_fixture()
        change_metadata(upper, 'Six-spinon absolute rapidity cutoff', '256')
        a, b = [sampling_points(d, momentum='1.5') for d in (lower, upper)]
        paired = paired_cutoffs(a, b)
        self.assertEqual(Decimal(paired['points'][0]['change']), 0)
        for key, value in [('Six-spinon fp128 density comparison bits', '53'),
                           ('Recovery seconds per configuration', '120'),
                           ('Reference density comparison bits', '127')]:
            changed = deepcopy(b)
            changed['metadata'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): paired_cutoffs(a, changed)
        with self.assertRaises(ValueError): paired_cutoffs(sampling_points(fixed_fixture(), momentum='1.5'), b)



if __name__ == '__main__':
    unittest.main()
