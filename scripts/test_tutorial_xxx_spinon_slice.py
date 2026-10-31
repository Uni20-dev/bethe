"""Do not interpolate, rescale, or mismatch the sectors in a comparison plot."""
from copy import deepcopy
from decimal import Decimal
import unittest

import plot_xxx_spinon_slice_tutorial as tutorial
from spinon_tutorial_exports import CONVENTION, read_spectrum
from test_tutorial_spinon_exports import fixture


def checked_inputs():
    document = fixture()
    document['tables']['spectrum']['metadata']['Command'] = '/build/bethe-xxx-structure-factor-thermo --sectors 2,4'
    four = read_spectrum(document, momentum='1.5')
    six = {'metadata': {**four['metadata'], 'Selected spinon sectors': '6'},
           'summary': dict(four['summary']),
           'points': [{'omega': '2.0', 's6': '0.0625', 'sampling_se': '0.0125',
                       'largest_draw_fraction': '0.3', 'largest_replicate_fraction': '0.6',
                       'replicate_records': [{'q': '1.5', 'density': '0.05'},
                                             {'q': '1.5', 'density': '0.075'}]}]}
    return four, six


def compact_fixture():
    four, six = checked_inputs()
    f = {**tutorial.provenance(four, b'four fixture'), 'spectrum': four['spectrum']}
    s = {**tutorial.provenance(six, b'six fixture'), 'sampling_points': six['points']}
    return {'schema': tutorial.SCHEMA, 'channel': 'zz', 'exchange': '1', 'momentum': '1.5',
            'density_convention': CONVENTION, 'source_runs': [f, s],
            'points': tutorial.join_points([four], [six])}


class SpinonSlice(unittest.TestCase):
    def test_sums_physical_weights_at_identical_decimal_frequencies(self):
        four, six = checked_inputs()
        point, = tutorial.join_points([four], [six])
        self.assertEqual(Decimal(point['s2_plus_s4']), Decimal('0.625'))
        self.assertEqual(Decimal(point['s2_plus_s4_plus_s6']), Decimal('0.6875'))
        self.assertEqual(point['s6_sampling_se'], '0.0125')
        self.assertEqual((point['four_run'], point['six_run']), (0, 1))

    def test_missing_or_duplicate_frequencies_are_rejected(self):
        for kind in ('four_duplicate', 'six_duplicate', 'grid'):
            with self.subTest(kind=kind):
                four, six = checked_inputs()
                if kind == 'grid':
                    six['points'][0]['omega'] = '2.01'
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
        six['points'][0]['replicate_records'][0]['q'] = '1.6'
        with self.assertRaisesRegex(ValueError, 'momenta'):
            tutorial.join_points([four], [six])

    def test_faint_tail_survives_above_two_spinon_support(self):
        four, six = checked_inputs()
        four['spectrum'][0]['density'] = '0'
        four['spectrum'][1]['density'] = '1e-100'
        six['points'][0]['s6'] = '2e-101'
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
                elif kind == 'exchange': bad['source_runs'][0]['metadata']['Exchange J'] = '2'
                elif kind == 'failed': bad['source_runs'][0]['summary']['Outcome'] = 'partial'
                else: bad['source_runs'][0]['metadata']['Bethe revision'] += '-dirty'
                with self.assertRaises(ValueError):
                    tutorial.validate(bad)

    def test_regeneration_commands_use_installed_program_and_local_outputs(self):
        four, _ = checked_inputs()
        four['metadata']['Command'] += ' --json /archive/machine/run.json'
        command = tutorial.provenance(four, b'fixture')['command']
        self.assertEqual(command[0], 'bethe-xxx-structure-factor-thermo')
        self.assertEqual(command[-1], 'run.json')


if __name__ == '__main__':
    unittest.main()
