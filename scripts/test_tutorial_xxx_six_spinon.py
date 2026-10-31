"""Reject incomplete or misleading six-spinon sampling summaries."""
from copy import deepcopy
from decimal import Decimal, localcontext
import unittest

from plot_xxx_six_spinon_tutorial import sampling_points
from spinon_tutorial_exports import CONVENTION


def fixture():
    meta = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
            'Density convention': CONVENTION, 'Geometry and output bits': '192',
            'Exchange J': '1', 'Selected spinon sectors': '6', 'Spectral points': '1',
            'Bethe revision': 'a' * 40, 'Uni20 revision': 'b' * 40,
            'Six-spinon replicates': '2', 'Samples per replicate': '4',
            'Six-spinon stream identities': 'first_replicate=grid_index*replicates; distinct streams at each point'}
    summary = {'Outcome': 'success', 'Failed points': '0'}
    def table(columns, rows):
        return {'metadata': dict(meta), 'summary': dict(summary),
                'columns': [{'id': x} for x in columns], 'rows': rows}
    return {'tables': {
        'spectrum': table(['q', 'omega', 'sector', 'density', 'status', 'evaluations'],
                          [['1.5', '2', '6', '2', 'estimated', '12']]),
        'diagnostics': table(['q', 'omega', 'sector', 'sampling_error', 'contour_diagnostic',
                              'arithmetic_diagnostic', 'evaluations'],
                             [['1.5', '2', '6', '1', '0', '0', '12']]),
        'replicates': table(['q', 'omega', 'replicate', 'completed_samples', 'supported_roots',
                             'density_attempts', 'maximum_absolute_density', 'density',
                             'contour_diagnostic', 'arithmetic_diagnostic', 'failed_sample'],
                            [['1.5', '2', '0', '4', '2', '4', '2', '1', '0', '0', None],
                             ['1.5', '2', '1', '4', '4', '8', '6', '3', '0', '0', None]])}}


class SixSpinonSampling(unittest.TestCase):
    def test_recomputes_se_of_mean_and_weight_concentrations(self):
        point = sampling_points(fixture(), momentum='1.5')['points'][0]
        self.assertEqual(Decimal(point['s6']), 2)
        self.assertEqual(Decimal(point['sampling_se']), 1)
        self.assertEqual(Decimal(point['largest_draw_fraction']), Decimal('0.375'))
        self.assertEqual(Decimal(point['largest_replicate_fraction']), Decimal('0.75'))

    def test_wrong_se_or_mean_is_rejected(self):
        for table, column, value in [('spectrum', 3, '2.1'), ('diagnostics', 3, '1.4142135623730951')]:
            with self.subTest(table=table):
                data = fixture()
                data['tables'][table]['rows'][0][column] = value
                with self.assertRaisesRegex(ValueError, 'mean/SE'):
                    sampling_points(data, momentum='1.5')

    def test_missing_duplicate_or_failed_replicate_is_rejected(self):
        for kind in ('missing', 'duplicate', 'failed', 'count', 'id', 'maximum', 'attempts'):
            with self.subTest(kind=kind):
                data = fixture()
                rows = data['tables']['replicates']['rows']
                if kind == 'missing': rows.pop()
                elif kind == 'duplicate': rows[1] = deepcopy(rows[0])
                elif kind == 'failed': rows[1][-1] = '3'
                elif kind == 'count': rows[1][3] = '3'
                elif kind == 'id': rows[1][2] = '2'
                elif kind == 'maximum': rows[1][6] = '2'
                else: rows[1][5] = '7'
                with self.assertRaises(ValueError):
                    sampling_points(data, momentum='1.5')

    def test_cancellation_sensitive_decimal_means(self):
        data = fixture()
        with localcontext() as context:
            context.prec = 100
            small = Decimal('1e-40')
            rows = data['tables']['replicates']['rows']
            rows[0][7], rows[1][7] = str(1 - small), str(1 + small)
            data['tables']['spectrum']['rows'][0][3] = '1'
            data['tables']['diagnostics']['rows'][0][3] = str(small)
            point = sampling_points(data, momentum='1.5')['points'][0]
            self.assertEqual(Decimal(point['sampling_se']), small)


if __name__ == '__main__':
    unittest.main()
