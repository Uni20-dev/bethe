"""Reviewed source cohorts cannot hide changed controls, gaps or replacements."""
from copy import deepcopy
from decimal import Decimal
from importlib.util import find_spec
import json
from pathlib import Path
import unittest

from spinon_tutorial_provenance import (PHYSICAL_KEYS, SCHEMA, compatibility_sources,
                                       physical_identity, validate_physics, validate_sources)
from test_tutorial_xxx_six_spinon import fixture
from xxx_six_spinon_tutorial_data import (NODE_CONTROLS, cohort_metadata, compact_slice,
                                        compact_slice_set, read_compact_slice, read_slice)

OLD = {'Bethe revision': 'a' * 40, 'Uni20 revision': 'b' * 40}
NEW = {'Bethe revision': 'c' * 40, 'Uni20 revision': 'b' * 40}


def cohort_fixture(name, source, omega, seed):
    data = fixture()
    for table in data['tables'].values():
        table['metadata'].update(source)
        table['metadata']['Six-spinon seed'] = seed
        for row in table['rows']:
            row[1] = omega
    compact = compact_slice([json.dumps(data).encode()], momentum='1.5',
                            frequencies=[omega], rapidity_cutoff='64')
    return {'id': name, 'source': dict(source), 'slice': compact}


def set_fixture():
    cohorts = [cohort_fixture('replicas', OLD, '2', '42'), cohort_fixture('blocks', NEW, '3', '43')]
    record = {'schema': SCHEMA, 'sources': [dict(OLD), dict(NEW)],
              'reason': 'Reviewed execution-only migration; synthetic test evidence',
              'evidence': [{'description': 'Matched-sample comparison', 'sha256': 'd' * 64}],
              'sector_sources': {'2,4': [dict(OLD)], '6': [dict(OLD), dict(NEW)]},
              'cohort_controls': {}}
    for cohort in cohorts:
        metadata = cohort['slice']['runs'][0]['metadata']
        record['cohort_controls'][cohort['id']] = {
            'source': dict(cohort['source']), 'binary_sha256': 'e' * 64, 'manifest_sha256': 'f' * 64,
            'metadata': cohort_metadata(metadata),
            'nodes': [{'omega': cohort['slice']['frequencies'][0],
                       'metadata': {key: metadata[key] for key in NODE_CONTROLS}}]}
    return compact_slice_set(cohorts, record, momentum='1.5', frequencies=['2', '3'], rapidity_cutoff='64')


class PhysicalIdentity(unittest.TestCase):
    def test_every_identity_field_is_required_and_missing_is_not_an_agreement(self):
        original = fixture()['tables']['spectrum']['metadata']
        for key in PHYSICAL_KEYS:
            for value in (None, '', '   ', 1):
                with self.subTest(key=key, value=value):
                    metadata = dict(original)
                    metadata[key] = value
                    with self.assertRaisesRegex(ValueError, key):
                        validate_physics([metadata, metadata])
        with self.assertRaisesRegex(ValueError, 'No physical'):
            validate_physics([])

    def test_each_physical_field_is_compared_without_source_compatibility(self):
        original = fixture()['tables']['spectrum']['metadata']
        for key in PHYSICAL_KEYS:
            with self.subTest(key=key):
                changed = dict(original)
                changed[key] = '2' if key == 'Exchange J' else 'different ' + changed[key]
                with self.assertRaisesRegex(ValueError, key):
                    validate_physics([original, changed])
        changed = dict(original, **NEW)
        changed.update({'Exchange J': '1.000', 'Density guard bits': '320',
                        'Selected spinon sectors': '2,4'})
        self.assertEqual(validate_physics([original, changed]), physical_identity(original))
        for value in ('NaN', 'sNaN', 'Infinity', '-Infinity', '0', '-1', 'unknown'):
            with self.subTest(exchange=value), self.assertRaisesRegex(ValueError, 'Exchange J'):
                physical_identity(dict(original, **{'Exchange J': value}))


class SourceCompatibility(unittest.TestCase):
    def test_default_join_stays_strict_and_allowlist_is_sector_scoped(self):
        record = set_fixture()['source_compatibility']
        metadata = [{**OLD, 'Selected spinon sectors': '2,4'}, {**NEW, 'Selected spinon sectors': '6'}]
        with self.assertRaisesRegex(ValueError, 'different source'):
            validate_sources(metadata)
        self.assertEqual(len(validate_sources(metadata, record)), 2)
        metadata[1]['Selected spinon sectors'] = '2,4'
        with self.assertRaisesRegex(ValueError, 'these sectors'):
            validate_sources(metadata, record)
        metadata[1]['Bethe revision'] = 'f' * 40
        with self.assertRaisesRegex(ValueError, 'outside'):
            validate_sources(metadata, record)

    def test_missing_malformed_or_duplicate_evidence_fails(self):
        changes = [('reason', ''), ('evidence', []), ('sources', [OLD, OLD]),
                   ('evidence', [{'description': 'claimed', 'sha256': 'not-a-hash'}]),
                   ('evidence', [{'description': '', 'sha256': 'd' * 64}]),
                   ('sector_sources', {'6': [OLD, OLD]}), ('sector_sources', {})]
        for key, value in changes:
            with self.subTest(key=key, value=value):
                record = set_fixture()['source_compatibility']
                record[key] = value
                with self.assertRaises(ValueError):
                    compatibility_sources(record)
        record = set_fixture()['source_compatibility']
        record['evidence'] *= 2
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            compatibility_sources(record)


class SourceCohorts(unittest.TestCase):
    def test_roundtrip_preserves_cohort_local_indices_and_all_original_records(self):
        document = set_fixture()
        before = deepcopy(document)
        checked = read_slice(json.loads(json.dumps(document)), momentum='1.5')
        self.assertEqual(document, before)
        self.assertEqual([p['omega'] for p in checked['points']], ['2', '3'])
        self.assertEqual([p['cohort_id'] for p in checked['points']], ['replicas', 'blocks'])
        self.assertEqual([p['run_index'] for p in checked['points']], [0, 0])
        for cohort, point in zip(document['cohorts'], checked['points']):
            expected = read_compact_slice(cohort['slice'], momentum='1.5')['points'][0]
            self.assertEqual({k: point[k] for k in expected}, expected)
            self.assertEqual(point['source'], cohort['source'])
        self.assertEqual(checked['source_compatibility'], document['source_compatibility'])

    def test_single_source_path_and_its_rejection_are_unchanged(self):
        document = set_fixture()
        first, second = [c['slice'] for c in document['cohorts']]
        self.assertEqual(read_slice(first, momentum='1.5'), read_compact_slice(first, momentum='1.5'))
        first['runs'] += second['runs']
        first['frequencies'] = ['2', '3']
        with self.assertRaisesRegex(ValueError, 'Mixed source revisions'):
            read_slice(first, momentum='1.5')
        # Adding the new compatibility field cannot bypass the strict v1 path.
        first['source_compatibility'] = document['source_compatibility']
        with self.assertRaisesRegex(ValueError, 'Mixed source revisions'):
            read_slice(first, momentum='1.5')

    def test_duplicate_replacement_and_undeclared_gap_fail(self):
        document = set_fixture()
        document['cohorts'][1] = cohort_fixture('blocks', NEW, '2', '43')
        document['source_compatibility']['cohort_controls']['blocks']['nodes'][0]['omega'] = '2'
        with self.assertRaisesRegex(ValueError, 'Duplicate frequency'):
            read_slice(document, momentum='1.5')
        document = set_fixture()
        document['frequencies'] = ['2', '2.5', '3']
        with self.assertRaisesRegex(ValueError, 'Missing frequency'):
            read_slice(document, momentum='1.5')
        document = set_fixture()
        document['frequencies'] = ['2']
        with self.assertRaisesRegex(ValueError, 'not requested'):
            read_slice(document, momentum='1.5')

    def test_declared_missing_point_stays_missing_with_cohort_identity(self):
        document = set_fixture()
        pending = document['cohorts'][1]['slice']
        pending['runs'] = []
        pending['unavailable'] = [{'omega': '3', 'status': 'failed', 'reason': 'density recovery failed', 'raw_sha256': 'e' * 64}]
        point = read_slice(document, momentum='1.5')['points'][1]
        self.assertIsNone(point['s6'])
        self.assertIsNone(point['sampling_se'])
        self.assertIsNone(point['run_index'])
        self.assertEqual(point['cohort_id'], 'blocks')
        self.assertEqual(point['source'], NEW)
        self.assertEqual(point['status'], 'failed')

    def test_exact_control_contract_rejects_otherwise_valid_seed_and_setting_changes(self):
        for key, value in [('Six-spinon seed', '44'), ('Density guard bits', '512'),
                           ('Maximum samples per Async task', '8')]:
            with self.subTest(key=key):
                document = set_fixture()
                changed = document['cohorts'][1]['slice']
                changed['runs'][0]['metadata'][key] = value
                # Ordinary validation alone still accepts these metadata changes.
                read_compact_slice(changed, momentum='1.5')
                with self.assertRaisesRegex(ValueError, 'controls differ'):
                    read_slice(document, momentum='1.5')

    def test_missing_cohort_cannot_declare_different_physics(self):
        for key, value in [('Channel', 'plus-minus'), ('Exchange J', '2'),
                           ('Six-spinon absolute rapidity cutoff', '128')]:
            with self.subTest(key=key):
                document = set_fixture()
                pending = document['cohorts'][1]['slice']
                pending['runs'] = []
                pending['unavailable'] = [{'omega': '3', 'status': 'unavailable',
                                           'reason': 'not yet calculated', 'raw_sha256': None}]
                document['source_compatibility']['cohort_controls']['blocks']['metadata'][key] = value
                with self.assertRaisesRegex(ValueError, 'physic'):
                    read_slice(document, momentum='1.5')

    def test_bound_sampling_counts_reject_consistently_edited_export(self):
        document = set_fixture()
        changed = document['cohorts'][1]['slice']
        changed['runs'][0]['metadata']['Samples per replicate'] = '8'
        for row in changed['runs'][0]['replicates']:
            row['completed_samples'] = '8'
        read_compact_slice(changed, momentum='1.5')
        with self.assertRaisesRegex(ValueError, 'sampling controls differ'):
            read_slice(document, momentum='1.5')

    def test_reviewed_sources_cannot_combine_different_models_or_normalizations(self):
        for key, value in [('Hamiltonian', 'XXZ Delta=2'), ('Family', 'Fitted intensity')]:
            for available in (True, False):
                with self.subTest(key=key, available=available):
                    document = set_fixture()
                    changed = document['cohorts'][1]['slice']
                    document['source_compatibility']['cohort_controls']['blocks']['metadata'][key] = value
                    if available:
                        changed['runs'][0]['metadata'][key] = value
                    else:
                        changed['runs'] = []
                        changed['unavailable'] = [{'omega': '3', 'status': 'unavailable',
                                                   'reason': 'not yet calculated', 'raw_sha256': None}]
                    # The per-cohort contract is consistent; the inter-cohort physics is not.
                    with self.assertRaisesRegex(ValueError, 'physical'):
                        read_slice(document, momentum='1.5')

    def test_missing_physical_identity_cannot_be_hidden_by_matching_contracts(self):
        for key in ('Hamiltonian', 'Family'):
            for value in (None, '', '   '):
                with self.subTest(key=key, value=value):
                    document = set_fixture()
                    for cohort in document['cohorts']:
                        metadata = [cohort['slice']['runs'][0]['metadata'],
                                    document['source_compatibility']['cohort_controls'][cohort['id']]['metadata']]
                        for item in metadata:
                            if value is None:
                                item.pop(key)
                            else:
                                item[key] = value
                    with self.assertRaisesRegex(ValueError, key):
                        read_slice(document, momentum='1.5')

    def test_reviewed_numerical_control_differences_remain_allowed(self):
        document = set_fixture()
        for metadata in (document['cohorts'][1]['slice']['runs'][0]['metadata'],
                         document['source_compatibility']['cohort_controls']['blocks']['metadata']):
            metadata.update({'Density guard bits': '512', 'Maximum samples per Async task': '16',
                             'Scheduler concurrency limit': '24', 'Exchange J': '1.0'})
        checked = read_slice(document, momentum='1.5')
        self.assertEqual([p['s6'] for p in checked['points']], ['2', '2'])

    def test_source_contract_and_evidence_are_required_even_with_valid_inner_slices(self):
        mutations = [lambda d: d['cohorts'][1].update(source=OLD),
                     lambda d: d['cohorts'][1].update(id='replicas'),
                     lambda d: d['source_compatibility']['cohort_controls']['blocks'].update(binary_sha256='unknown'),
                     lambda d: d['source_compatibility']['cohort_controls']['blocks']['nodes'][0].update(omega='4'),
                     lambda d: d['source_compatibility'].update(evidence=[]),
                     lambda d: d['source_compatibility']['sources'].append({'Bethe revision': 'f'*40, 'Uni20 revision': 'b'*40}),
                     lambda d: d['source_compatibility']['cohort_controls'].update(unused={})]
        for mutate in mutations:
            document = set_fixture()
            mutate(document)
            with self.assertRaises(ValueError):
                read_slice(document, momentum='1.5')

    def test_numerical_checks_still_recompute_replicate_error(self):
        document = set_fixture()
        document['cohorts'][1]['slice']['runs'][0]['diagnostics'][0]['sampling_error'] = '0'
        with self.assertRaisesRegex(ValueError, 'mean/SE'):
            read_slice(document, momentum='1.5')

    def test_only_date_and_invocation_are_excluded_from_common_metadata(self):
        document = set_fixture()
        metadata = document['cohorts'][1]['slice']['runs'][0]['metadata']
        metadata.update(Date='another execution date', Command='program --json anywhere.json')
        self.assertEqual(read_slice(document, momentum='1.5')['points'][1]['s6'], '2')

    @unittest.skipUnless(find_spec('matplotlib'), 'Install plotting requirements')
    def test_standalone_plot_consumes_both_cohorts_without_rescaling(self):
        from plot_xxx_six_spinon_tutorial import spectrum_figure
        import matplotlib.pyplot as plt
        figure = spectrum_figure(set_fixture(), momentum='1.5')
        try:
            for axis in figure.axes:
                line = next(x for x in axis.lines if x.get_gid() == 'six-spinon-spectrum')
                self.assertEqual(list(line.get_xdata()), [2, 3])
                self.assertEqual(list(line.get_ydata()), [2, 2])
        finally:
            plt.close(figure)


class RecordedSchedulerControls(unittest.TestCase):
    def test_committed_control_recomputes_identical_nonzero_estimates(self):
        path = Path(__file__).resolve().parents[1] / 'docs/validation/six-spinon-scheduler-controls-20261101.json'
        record = json.loads(path.read_text())
        self.assertEqual([case['samples_per_task'] for case in record['cases']], [32, 16, 7])
        baseline = None
        for case in record['cases']:
            data = read_compact_slice(case['slice'], momentum=record['momentum'])
            point, = data['points']
            self.assertEqual(point['supported_roots'], 40)
            self.assertEqual(point['density_attempts'], 80)
            self.assertEqual(len(point['replicate_records']), 4)
            self.assertEqual(sum(int(row['completed_samples']) for row in point['replicate_records']), 128)
            self.assertGreater(Decimal(point['s6']), 0)
            run, = case['slice']['runs']
            tables = {key: run[key] for key in ('spectrum', 'diagnostics', 'replicates')}
            if baseline is None:
                baseline = tables
            else:
                self.assertEqual(tables, baseline)
            self.assertNotIn('Command', run['metadata'])


if __name__ == '__main__':
    unittest.main()
