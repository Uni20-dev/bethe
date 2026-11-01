"""Scientific plotting checks: gaps, unit conversion and validated error bars."""
from decimal import Decimal
from importlib.util import find_spec
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

if find_spec('matplotlib') is not None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as ArtistPath
else:
    plt = None

from plot_xxx_six_spinon_tutorial import plot_compact, spectrum_figure
from test_tutorial_xxx_six_spinon import fixture
from xxx_six_spinon_tutorial_data import compact_slice


def packed(data=None, frequencies=('2',), unavailable=()):
    return compact_slice([json.dumps(fixture() if data is None else data).encode()],
                         momentum='1.5', exchange='1' if data is None else data['tables']['spectrum']['metadata']['Exchange J'],
                         rapidity_cutoff='64', frequencies=list(frequencies), unavailable=list(unavailable))


def spectrum_line(axis):
    return next(line for line in axis.lines if line.get_gid() == 'six-spinon-spectrum')


@unittest.skipIf(plt is None, 'Install docs/tutorials/requirements.txt for rendering checks')
class SixPlot(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_gap_does_not_connect_available_estimates(self):
        first, last = fixture(), fixture()
        for table in last['tables'].values():
            for row in table['rows']:
                row[1] = '3'
        data = compact_slice([json.dumps(x).encode() for x in (first, last)], momentum='1.5',
                             rapidity_cutoff='64', frequencies=['2', '2.5', '3'],
                             unavailable=[{'omega': '2.5', 'status': 'failed',
                                           'reason': 'test failure', 'raw_sha256': None}])
        fig = spectrum_figure(data, momentum='1.5')
        for axis in fig.axes:
            segments = list(spectrum_line(axis).get_path().iter_segments(remove_nans=True))
            self.assertEqual([code for _, code in segments], [ArtistPath.MOVETO, ArtistPath.MOVETO])
        self.assertEqual(fig.axes[1].get_yscale(), 'log')

    def test_exchange_scaling_applies_to_coordinates_density_and_se(self):
        data = fixture()
        for table in data['tables'].values():
            table['metadata']['Exchange J'] = '2'
        fig = spectrum_figure(packed(data), momentum='1.5', exchange='2')
        for axis in fig.axes:
            self.assertEqual(list(spectrum_line(axis).get_xdata()), [1])
            self.assertEqual(list(spectrum_line(axis).get_ydata()), [4])
            segments = axis.containers[0].lines[2][0].get_segments()
            self.assertEqual(segments[0].tolist(), [[1, 2], [1, 6]])

    def test_missing_endpoints_retain_requested_domain_and_count(self):
        missing = [{'omega': omega, 'status': 'unavailable',
                    'reason': 'not computed', 'raw_sha256': None} for omega in ('1.8', '3.7')]
        fig = spectrum_figure(packed(frequencies=('1.8', '2', '3.7'), unavailable=missing), momentum='1.5')
        for axis in fig.axes:
            lower, upper = axis.get_xlim()
            self.assertLessEqual(lower, 1.8)
            self.assertGreaterEqual(upper, 3.7)
            segments = list(spectrum_line(axis).get_path().iter_segments(remove_nans=True))
            self.assertEqual([code for _, code in segments], [ArtistPath.MOVETO])
        self.assertIn('unavailable at 2/3 requested frequencies', fig._suptitle.get_text())

    def test_all_missing_draws_no_physical_curve(self):
        data = packed()
        data['runs'] = []
        data['unavailable'] = [{'omega': '2', 'status': 'unavailable',
                                'reason': 'not computed', 'raw_sha256': None}]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'missing.svg'
            plot_compact(data, momentum='1.5', figure=output)
            self.assertTrue(output.exists())
            self.assertTrue(output.with_name('missing-dark.svg').exists())
            self.assertIn('No available estimates', output.read_text())

    def test_svg_exports_ignore_ambient_id_salt(self):
        data = packed()
        with tempfile.TemporaryDirectory() as directory:
            outputs = [Path(directory) / name for name in ('first.svg', 'second.svg')]
            for salt, output in zip(('caller-one', 'caller-two'), outputs):
                with plt.rc_context({'svg.hashsalt': salt}):
                    plot_compact(data, momentum='1.5', figure=output)
                    self.assertEqual(plt.rcParams['svg.hashsalt'], salt)
            for suffix in ('', '-dark'):
                rendered = [p.with_name(p.stem + suffix + '.svg').read_bytes() for p in outputs]
                self.assertEqual(rendered[0], rendered[1])

    def test_single_missing_frequency_retains_scaled_coordinate(self):
        raw = fixture()
        for table in raw['tables'].values():
            table['metadata']['Exchange J'] = '2'
        data = packed(raw)
        data['runs'] = []
        data['unavailable'] = [{'omega': '2', 'status': 'unavailable',
                                'reason': 'not computed', 'raw_sha256': None}]
        fig = spectrum_figure(data, momentum='1.5', exchange='2')
        for axis in fig.axes:
            lower, upper = axis.get_xlim()
            self.assertLess(lower, 1)
            self.assertGreater(upper, 1)
            self.assertAlmostEqual((lower + upper) / 2, 1)

    def test_invalid_se_is_rejected_before_rendering(self):
        data = packed()
        data['runs'][0]['diagnostics'][0]['sampling_error'] = '0'
        with patch('matplotlib.pyplot.subplots') as subplots:
            with self.assertRaisesRegex(ValueError, 'mean/SE'):
                spectrum_figure(data, momentum='1.5')
            subplots.assert_not_called()

    def test_unrepresentable_faint_data_are_not_silently_zeroed(self):
        data = fixture()
        for row in data['tables']['replicates']['rows']:
            for column in (6, 7):
                row[column] = str(Decimal(row[column]) * Decimal('1e-1000'))
        data['tables']['spectrum']['rows'][0][3] = '2e-1000'
        data['tables']['diagnostics']['rows'][0][3] = '1e-1000'
        with self.assertRaisesRegex(ValueError, 'plotting precision'):
            spectrum_figure(packed(data), momentum='1.5')


if __name__ == '__main__':
    unittest.main()
