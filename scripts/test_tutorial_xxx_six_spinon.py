"""Reject incomplete or misleading six-spinon sampling summaries."""
from copy import deepcopy
from decimal import Decimal, localcontext
import unittest

from plot_xxx_six_spinon_tutorial import paired_cutoffs, sampling_points
from spinon_tutorial_exports import CONVENTION


def fixture():
    meta = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
            'Density convention': CONVENTION, 'Geometry and output bits': '192',
            'Exchange J': '1', 'Selected spinon sectors': '6', 'Spectral points': '1',
            'Bethe revision': 'a' * 40, 'Uni20 revision': 'b' * 40,
            'Six-spinon replicates': '2', 'Samples per replicate': '4',
            'Six-spinon seed': '42', 'Six-spinon uniform grid bits': '52',
            'Six-spinon absolute rapidity cutoff': '64',
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
        for kind in ('missing', 'duplicate', 'failed', 'count', 'id', 'maximum', 'overweight', 'attempts'):
            with self.subTest(kind=kind):
                data = fixture()
                rows = data['tables']['replicates']['rows']
                if kind == 'missing': rows.pop()
                elif kind == 'duplicate': rows[1] = deepcopy(rows[0])
                elif kind == 'failed': rows[1][-1] = '3'
                elif kind == 'count': rows[1][3] = '3'
                elif kind == 'id': rows[1][2] = '2'
                elif kind == 'maximum': rows[1][6] = '2'
                elif kind == 'overweight': rows[1][6] = '13'
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
            rows[0][6] = rows[1][6] = '2'
            data['tables']['spectrum']['rows'][0][3] = '1'
            data['tables']['diagnostics']['rows'][0][3] = str(small)
            point = sampling_points(data, momentum='1.5')['points'][0]
            self.assertEqual(Decimal(point['sampling_se']), small)

    def test_faint_tail_cannot_hide_a_zeroed_sampling_error(self):
        data = fixture()
        for row in data['tables']['replicates']['rows']:
            for column in (6, 7):
                row[column] = str(Decimal(row[column]) * Decimal('1e-100'))
        data['tables']['spectrum']['rows'][0][3] = '2e-100'
        data['tables']['diagnostics']['rows'][0][3] = '1e-100'
        point = sampling_points(data, momentum='1.5')['points'][0]
        self.assertEqual(Decimal(point['sampling_se']), Decimal('1e-100'))
        data['tables']['diagnostics']['rows'][0][3] = '0'
        with self.assertRaisesRegex(ValueError, 'mean/SE'):
            sampling_points(data, momentum='1.5')

    def test_paired_cutoff_error_retains_anticorrelation(self):
        lower, upper = fixture(), fixture()
        for table in upper['tables'].values():
            table['metadata']['Six-spinon absolute rapidity cutoff'] = '256'
        rows = upper['tables']['replicates']['rows']
        rows[0][6:8], rows[1][6:8] = ['6', '3'], ['2', '1']
        a, b = [sampling_points(d, momentum='1.5') for d in (lower, upper)]
        point = paired_cutoffs(a, b)['points'][0]
        self.assertEqual(Decimal(point['change']), 0)
        self.assertEqual(Decimal(point['paired_sampling_se']), 2)
        self.assertEqual(Decimal(point['replicate_correlation']), -1)
        self.assertEqual([Decimal(x) for x in point['replicate_differences']], [2, -2])

    def test_pairing_rejects_changed_seed_controls_or_grid(self):
        a = sampling_points(fixture(), momentum='1.5')
        for kind in ('seed', 'precision', 'grid', 'cutoff'):
            with self.subTest(kind=kind):
                b = deepcopy(a)
                b['metadata']['Six-spinon absolute rapidity cutoff'] = '256'
                if kind == 'seed': b['metadata']['Six-spinon seed'] = '43'
                elif kind == 'precision': b['metadata']['Geometry and output bits'] = '256'
                elif kind == 'grid': b['points'][0]['omega'] = '2.1'
                else: b['metadata']['Six-spinon absolute rapidity cutoff'] = '32'
                with self.assertRaises(ValueError):
                    paired_cutoffs(a, b)


if __name__ == '__main__':
    unittest.main()
