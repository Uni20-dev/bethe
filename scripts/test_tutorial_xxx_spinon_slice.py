"""Do not interpolate, rescale, or mismatch the sectors in a comparison plot."""
from copy import deepcopy
from decimal import Decimal, localcontext
from importlib.util import find_spec
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import json
import math
from unittest.mock import patch

import plot_xxx_spinon_slice_tutorial as tutorial
from spinon_tutorial_exports import CONVENTION, read_spectrum
from test_tutorial_spinon_exports import fixture
from test_tutorial_xxx_six_spinon import fixture as six_fixture
from xxx_six_spinon_tutorial_data import compact_slice, read_compact_slice


def six_document():
    data = six_fixture()
    rows = data['tables']['replicates']['rows']
    rows[0][6:8], rows[1][6:8] = ['0.1', '0.05'], ['0.15', '0.075']
    data['tables']['spectrum']['rows'][0][3] = '0.0625'
    data['tables']['diagnostics']['rows'][0][3] = '0.0125'
    return compact_slice([json.dumps(data).encode()], momentum='1.5',
                         frequencies=['2.0'], rapidity_cutoff='64')


def checked_inputs():
    document = fixture()
    document['tables']['spectrum']['metadata']['Command'] = '/build/bethe-xxx-structure-factor-thermo --sectors 2,4'
    four = read_spectrum(document, momentum='1.5')
    six = read_compact_slice(six_document(), momentum='1.5')['points'][0]
    return four, six


def compact_fixture():
    four, six = checked_inputs()
    f = {**tutorial.provenance(four, b'four fixture'), 'spectrum': four['spectrum'], 'diagnostics': []}
    return {'schema': tutorial.SCHEMA, 'channel': 'zz', 'exchange': '1', 'momentum': '1.5',
            'density_convention': CONVENTION, 'four_runs': [f], 'six_slice': six_document(),
            'points': tutorial.join_points([four], [six])}


def cohort_inputs():
    from test_tutorial_spinon_cohorts import set_fixture
    document = set_fixture()
    four, _ = checked_inputs()
    runs = []
    for omega in document['frequencies']:
        run = deepcopy(four)
        for row in run['spectrum']:
            row['omega'] = omega
        runs.append(run)
    return runs, document


def cohort_compact_fixture():
    from xxx_six_spinon_tutorial_data import read_slice
    four, document = cohort_inputs()
    six = read_slice(document, momentum='1.5')
    return {'schema': tutorial.SCHEMA, 'channel': 'zz', 'exchange': '1', 'momentum': '1.5',
            'density_convention': CONVENTION,
            'four_runs': [{**tutorial.provenance(run, b'four cohort fixture'),
                           'spectrum': run['spectrum'], 'diagnostics': []} for run in four],
            'six_slice': document,
            'points': tutorial.join_points(four, six['points'],
                                          compatibility=six['source_compatibility'])}


class CombinedSourceCohorts(unittest.TestCase):
    def test_collect_retains_full_compatibility_and_cohort_local_run_indices(self):
        four, document = cohort_inputs()
        original = deepcopy(document)
        with TemporaryDirectory() as directory:
            paths = []
            for index, run in enumerate(four):
                raw = fixture()
                raw['tables']['spectrum']['metadata']['Command'] = run['metadata']['Command']
                for row in raw['tables']['spectrum']['rows']:
                    row[1] = document['frequencies'][index]
                path = Path(directory) / f'four-{index}.json'
                path.write_text(json.dumps(raw))
                paths.append(path)
            data = tutorial.collect(paths, document, momentum='1.5')
        self.assertEqual(document, original)
        self.assertEqual(data['six_slice'], original)
        self.assertEqual([(p['six_cohort'], p['six_run']) for p in data['points']],
                         [('replicas', 0), ('blocks', 0)])
        for point in data['points']:
            self.assertEqual(Decimal(point['s2_plus_s4']), Decimal('0.625'))
            self.assertEqual(Decimal(point['s6']), Decimal('2'))
            self.assertEqual(Decimal(point['s2_plus_s4_plus_s6']), Decimal('2.625'))
        self.assertEqual(tutorial.validate(json.loads(json.dumps(data))), data)

    def test_different_sources_still_require_explicit_compatibility(self):
        from xxx_six_spinon_tutorial_data import read_slice
        four, document = cohort_inputs()
        six = read_slice(document, momentum='1.5')
        with self.assertRaisesRegex(ValueError, 'different source'):
            tutorial.join_points(four, six['points'])

    def test_compatibility_does_not_admit_new_or_unrelated_four_sources(self):
        from xxx_six_spinon_tutorial_data import read_slice
        for revision in ('c' * 40, 'f' * 40):
            with self.subTest(revision=revision):
                four, document = cohort_inputs()
                six = read_slice(document, momentum='1.5')
                four[0]['metadata']['Bethe revision'] = revision
                with self.assertRaises(ValueError):
                    tutorial.join_points(four, six['points'], compatibility=six['source_compatibility'])

    def test_physical_checks_are_independent_of_source_permission(self):
        from xxx_six_spinon_tutorial_data import read_slice
        for key, value in [('Channel', 'raising'), ('Exchange J', '2'),
                           ('Density convention', 'another measure')]:
            with self.subTest(key=key):
                four, document = cohort_inputs()
                six = read_slice(document, momentum='1.5')
                six['points'][1]['metadata'][key] = value
                with self.assertRaisesRegex(ValueError, 'physical conventions'):
                    tutorial.join_points(four, six['points'], compatibility=six['source_compatibility'])
        four, document = cohort_inputs()
        six = read_slice(document, momentum='1.5')
        six['points'][1]['replicate_records'][0]['q'] = '1.6'
        with self.assertRaisesRegex(ValueError, 'momenta'):
            tutorial.join_points(four, six['points'], compatibility=six['source_compatibility'])

    def test_failed_new_cohort_point_keeps_its_identity_and_gap(self):
        from xxx_six_spinon_tutorial_data import read_slice
        data = cohort_compact_fixture()
        missing = data['six_slice']['cohorts'][1]['slice']
        missing['runs'] = []
        missing['unavailable'] = [{'omega': '3', 'status': 'failed',
                                   'reason': 'density recovery exhausted', 'raw_sha256': 'd' * 64}]
        four = [tutorial.checked_four_run(run, momentum='1.5', exchange='1') for run in data['four_runs']]
        six = read_slice(data['six_slice'], momentum='1.5')
        data['points'] = tutorial.join_points(four, six['points'], compatibility=six['source_compatibility'])
        self.assertIs(tutorial.validate(data), data)
        point = data['points'][1]
        self.assertEqual((point['six_cohort'], point['six_run'], point['s6_status']),
                         ('blocks', None, 'failed'))
        self.assertIsNone(point['s6'])
        self.assertIsNone(point['s2_plus_s4_plus_s6'])
        self.assertEqual(Decimal(point['s2_plus_s4']), Decimal('0.625'))

    def test_embedded_evidence_controls_and_point_provenance_are_revalidated(self):
        mutations = [lambda d: d['six_slice']['source_compatibility'].update(evidence=[]),
                     lambda d: d['six_slice']['cohorts'][1]['slice']['runs'][0]['metadata'].update(
                         {'Six-spinon seed': '99'}),
                     lambda d: d['points'][1].update(six_cohort='replicas'),
                     lambda d: d['points'][1].update(six_run=1),
                     lambda d: d['six_slice'].update(frequencies=['2', '2.5', '3'])]
        for mutate in mutations:
            data = cohort_compact_fixture()
            mutate(data)
            with self.assertRaises(ValueError):
                tutorial.validate(data)


class SpinonSlice(unittest.TestCase):
    def assert_plot_range_rejection(self, data):
        four = [tutorial.checked_four_run(run, momentum=data['momentum'], exchange='1')
                for run in data['four_runs']]
        six = read_compact_slice(data['six_slice'], momentum=data['momentum'])
        data['points'] = tutorial.join_points(four, six['points'])
        self.assertIs(tutorial.validate(data), data)
        with patch.object(tutorial, 'save_svg') as save:
            with self.assertRaisesRegex(ValueError, 'plotting precision'):
                tutorial.plot(data)
            save.assert_not_called()

    def test_valid_compact_density_outside_plotting_range_is_not_hidden(self):
        for factor in ('1e-1000', '1e1000'):
            with self.subTest(factor=factor), localcontext() as context:
                context.prec = 100
                data = compact_fixture()
                for run in data['six_slice']['runs']:
                    for row in run['replicates']:
                        for key in ('maximum_absolute_density', 'density'):
                            row[key] = str(Decimal(row[key]) * Decimal(factor))
                    for row in run['spectrum']:
                        row['density'] = str(Decimal(row['density']) * Decimal(factor))
                    for row in run['diagnostics']:
                        row['sampling_error'] = str(Decimal(row['sampling_error']) * Decimal(factor))
                self.assert_plot_range_rejection(data)

    def test_valid_tiny_sampling_error_is_not_rounded_to_zero_for_plotting(self):
        data = compact_fixture()
        run = data['six_slice']['runs'][0]
        run['metadata']['Geometry and output bits'] = '4096'
        with localcontext() as context:
            context.prec = 1400
            for row, sign in zip(run['replicates'], (-1, 1)):
                row['density'] = str(Decimal('0.0625') + sign * Decimal('1e-400'))
                row['maximum_absolute_density'] = '0.125'
        run['diagnostics'][0]['sampling_error'] = '1e-400'
        self.assert_plot_range_rejection(data)

    def test_valid_compact_coordinates_must_fit_plotting_range(self):
        for coordinate in ('omega', 'q'):
            for value in ('1e-1000', '1e1000'):
                with self.subTest(coordinate=coordinate, value=value):
                    data = compact_fixture()
                    for run in data['four_runs'] + data['six_slice']['runs']:
                        for table in ('spectrum', 'diagnostics', 'replicates'):
                            for row in run.get(table, []):
                                row[coordinate] = value
                    if coordinate == 'q':
                        data['momentum'] = data['six_slice']['momentum'] = value
                    else:
                        data['six_slice']['frequencies'] = [value]
                    self.assert_plot_range_rejection(data)

    def test_sums_physical_weights_at_identical_decimal_frequencies(self):
        four, six = checked_inputs()
        point, = tutorial.join_points([four], [six])
        self.assertEqual(Decimal(point['s2_plus_s4']), Decimal('0.625'))
        self.assertEqual(Decimal(point['s2_plus_s4_plus_s6']), Decimal('0.6875'))
        self.assertEqual(point['s6_sampling_se'], '0.0125')
        self.assertEqual((point['four_run'], point['six_run']), (0, 0))

    def test_declared_decimal_grid_matches_distinct_arithmetic_formats(self):
        four, six = checked_inputs()
        computed_four = '1.799999999999999999999999999999999999998'
        computed_six = '1.79999999999999999996'
        four['metadata']['Geometry and output bits'] = '128'
        for row in four['spectrum']:
            row['omega'] = computed_four
        six['metadata']['Geometry and output bits'] = '64'
        six['omega'], six['computed_omega'] = '1.8', computed_six
        for row in six['replicate_records']:
            row['omega'] = computed_six
        point, = tutorial.join_points([four], [six])
        self.assertEqual(point['omega'], '1.8')
        self.assertEqual(point['four_computed_omega'], computed_four)
        self.assertEqual(point['six_computed_omega'], computed_six)
        self.assertEqual(Decimal(point['s2_plus_s4_plus_s6']), Decimal('0.6875'))

    def test_two_computed_points_cannot_collapse_onto_one_requested_point(self):
        four, six = checked_inputs()
        four['metadata']['Geometry and output bits'] = '64'
        four['metadata']['Spectral points'] = '2'
        original = deepcopy(four['spectrum'])
        four['spectrum'] = []
        for omega in ('1.79999999999999999996', '1.80000000000000000007'):
            for row in original:
                four['spectrum'].append(dict(row, omega=omega))
        six['omega'] = '1.8'
        with self.assertRaisesRegex(ValueError, 'Duplicate two/four frequency'):
            tutorial.join_points([four], [six])

    def test_missing_or_duplicate_frequencies_are_rejected(self):
        for kind in ('four_duplicate', 'six_duplicate', 'grid'):
            with self.subTest(kind=kind):
                four, six = checked_inputs()
                if kind == 'grid':
                    six['omega'] = '2.01'
                with self.assertRaises(ValueError):
                    tutorial.join_points([four] * (2 if kind == 'four_duplicate' else 1),
                                         [six] * (2 if kind == 'six_duplicate' else 1))

    def test_revisions_conventions_and_momenta_must_match(self):
        for key, value in [('Bethe revision', 'c' * 40), ('Channel', 'raising'), ('Exchange J', '2')]:
            with self.subTest(key=key):
                four, six = checked_inputs()
                six['metadata'][key] = value
                with self.assertRaises(ValueError):
                    tutorial.join_points([four], [six])
        four, six = checked_inputs()
        six['replicate_records'][0]['q'] = '1.6'
        with self.assertRaisesRegex(ValueError, 'momenta'):
            tutorial.join_points([four], [six])

    def test_faint_tail_survives_above_two_spinon_support(self):
        four, six = checked_inputs()
        four['spectrum'][0]['density'] = '0'
        four['spectrum'][1]['density'] = '1e-100'
        six['s6'] = '2e-101'
        point, = tutorial.join_points([four], [six])
        self.assertEqual(Decimal(point['s2_plus_s4_plus_s6']), Decimal('1.2e-100'))

    def test_plot_records_remain_attached_to_their_sources(self):
        data = compact_fixture()
        self.assertIs(tutorial.validate(data), data)
        for kind in ('sum', 'sector', 'momentum', 'exchange', 'failed', 'source'):
            with self.subTest(kind=kind):
                bad = deepcopy(data)
                if kind == 'sum': bad['points'][0]['s2_plus_s4_plus_s6'] = '1'
                elif kind == 'sector': bad['points'][0]['s6'] = '0'
                elif kind == 'momentum': bad['momentum'] = '1.6'
                elif kind == 'exchange': bad['four_runs'][0]['metadata']['Exchange J'] = '2'
                elif kind == 'failed': bad['four_runs'][0]['summary']['Outcome'] = 'partial'
                else: bad['four_runs'][0]['metadata']['Bethe revision'] += '-dirty'
                with self.assertRaises(ValueError):
                    tutorial.validate(bad)

    def test_compact_six_statistics_are_rechecked(self):
        data = compact_fixture()
        data['six_slice']['runs'][0]['diagnostics'][0]['sampling_error'] = '0'
        with self.assertRaisesRegex(ValueError, 'mean/SE'):
            tutorial.validate(data)

    def test_unavailable_six_point_leaves_both_six_curves_missing(self):
        four, _ = checked_inputs()
        document = compact_slice([], momentum='1.5', frequencies=['2'], rapidity_cutoff='64',
                                 unavailable=[{'omega': '2', 'status': 'failed',
                                               'reason': 'density recovery exhausted', 'raw_sha256': 'c' * 64}])
        six = read_compact_slice(document, momentum='1.5')['points']
        point, = tutorial.join_points([four], six)
        self.assertIsNone(point['s6'])
        self.assertIsNone(point['s6_sampling_se'])
        self.assertIsNone(point['s2_plus_s4_plus_s6'])
        self.assertEqual(Decimal(point['s2_plus_s4']), Decimal('0.625'))
        self.assertEqual(point['s6_status'], 'failed')

    @unittest.skipIf(find_spec('matplotlib') is None,
                     'Install docs/tutorials/requirements.txt for rendering checks')
    def test_plot_identifies_unavailable_six_results_without_zeroing_them(self):
        for partly_available in (False, True):
            with self.subTest(partly_available=partly_available):
                data = compact_fixture()
                missing_omega = '2.1' if partly_available else '2'
                if partly_available:
                    extra_four = deepcopy(data['four_runs'][0])
                    for table in ('spectrum', 'diagnostics'):
                        for row in extra_four[table]:
                            row['omega'] = missing_omega
                    data['four_runs'].append(extra_four)
                data['six_slice'] = compact_slice(
                    [json.dumps(six_fixture()).encode()] if partly_available else [],
                    momentum='1.5', frequencies=['2', '2.1'] if partly_available else ['2'],
                    rapidity_cutoff='64', unavailable=[{
                        'omega': missing_omega, 'status': 'failed',
                        'reason': 'density recovery exhausted', 'raw_sha256': 'c' * 64}])
                four = [tutorial.checked_four_run(r, momentum='1.5', exchange='1') for r in data['four_runs']]
                six = read_compact_slice(data['six_slice'], momentum='1.5')['points']
                data['points'] = tutorial.join_points(four, six)
                with patch.object(tutorial, 'save_svg') as save:
                    tutorial.plot(data)
                figure = save.call_args.args[0]
                self.assertIn(f'unavailable at 1/{2 if partly_available else 1}', figure._suptitle.get_text())
                for axis in figure.axes:
                    six_container, = axis.containers
                    plotted = six_container.lines[0].get_ydata()
                    self.assertTrue(math.isnan(plotted[-1]))
                    if partly_available:
                        self.assertGreater(plotted[0], 0)
                        self.assertIn('sampling SE', six_container.get_label())
                    else:
                        self.assertIn('unavailable', six_container.get_label())
                        self.assertNotIn('sampling', six_container.get_label())

    @unittest.skipIf(find_spec('matplotlib') is None,
                     'Install docs/tutorials/requirements.txt for rendering checks')
    def test_zero_sectors_and_missing_six_retain_requested_frequency_range(self):
        for has_available_point in (False, True):
            with self.subTest(has_available_point=has_available_point):
                data = compact_fixture()
                zero = deepcopy(data['four_runs'][0])
                for row in zero['spectrum']:
                    row['omega'], row['density'] = '7', '0'
                    if row['sector'] in ('2', '4'):
                        row['status'] = 'physical_zero'
                data['four_runs'] = ([data['four_runs'][0]] if has_available_point else []) + [zero]
                if not has_available_point:
                    data['six_slice']['runs'] = []
                data['six_slice']['frequencies'] = ['2', '7'] if has_available_point else ['7']
                data['six_slice']['unavailable'] = [{'omega': '7', 'status': 'unavailable',
                    'reason': 'not computed', 'raw_sha256': None}]
                four = [tutorial.checked_four_run(r, momentum='1.5', exchange='1')
                        for r in data['four_runs']]
                six = read_compact_slice(data['six_slice'], momentum='1.5')['points']
                data['points'] = tutorial.join_points(four, six)
                with patch.object(tutorial, 'save_svg') as save:
                    tutorial.plot(data)
                figure = save.call_args.args[0]
                for axis in figure.axes:
                    lower, upper = axis.get_xlim()
                    self.assertLess(lower, 7)
                    self.assertGreater(upper, 7)
                    self.assertLessEqual(lower, math.pi / 2 * abs(math.sin(1.5)))
                    self.assertTrue(math.isnan(axis.containers[0].lines[0].get_ydata()[-1]))
                self.assertIsNone(data['points'][-1]['s2_plus_s4_plus_s6'])

    def test_regeneration_commands_use_installed_program_and_local_outputs(self):
        four, _ = checked_inputs()
        four['metadata']['Command'] += ' --json /archive/machine/run.json'
        command = tutorial.provenance(four, b'fixture')['command']
        self.assertEqual(command[0], 'bethe-xxx-structure-factor-thermo')
        self.assertEqual(command[-1], 'run.json')

    def test_check_import_validates_without_creating_or_overwriting_files(self):
        for existing in (False, True):
            for valid in (False, True):
                with self.subTest(existing=existing, valid=valid), TemporaryDirectory() as directory:
                    directory = Path(directory)
                    raw, six, output = [directory / name for name in ('four.json', 'six.json', 'joined.json')]
                    document = fixture()
                    document['tables']['spectrum']['metadata']['Command'] = 'bethe-xxx-structure-factor-thermo --sectors 2,4'
                    if not valid:
                        document['tables']['spectrum']['rows'][-1][3] = '0.7'
                    raw.write_text(json.dumps(document))
                    six.write_text(json.dumps(six_document()))
                    if existing:
                        output.write_text('Existing output must remain untouched.\n')
                    before = {p.name: p.read_bytes() for p in directory.iterdir()}
                    command = ['plot_xxx_spinon_slice_tutorial.py', '--four-json', str(raw),
                               '--six-data', str(six), '--data', str(output),
                               '--figure', str(directory / 'slice.svg'), '--momentum', '1.5', '--check']
                    with patch('sys.argv', command), patch('sys.stdout', new_callable=StringIO) as stdout, \
                            patch.object(tutorial, 'plot') as plot:
                        if valid:
                            tutorial.main()
                            self.assertIn('Checked 1 matched frequencies', stdout.getvalue())
                        else:
                            with self.assertRaisesRegex(ValueError, 'Partial sum'):
                                tutorial.main()
                        plot.assert_not_called()
                    self.assertEqual({p.name: p.read_bytes() for p in directory.iterdir()}, before)


if __name__ == '__main__':
    unittest.main()
