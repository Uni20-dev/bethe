"""Do not interpolate, rescale, or mismatch the sectors in a comparison plot."""
from copy import deepcopy
from decimal import Decimal, localcontext
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

    def test_regeneration_commands_use_installed_program_and_local_outputs(self):
        four, _ = checked_inputs()
        four['metadata']['Command'] += ' --json /archive/machine/run.json'
        command = tutorial.provenance(four, b'fixture')['command']
        self.assertEqual(command[0], 'bethe-xxx-structure-factor-thermo')
        self.assertEqual(command[-1], 'run.json')


if __name__ == '__main__':
    unittest.main()
