import itertools
import math
import unittest

from plot_exclusion_tutorial import CASES, DATA, SIZES, BIASES, load_cases, read_exports


def characteristic_determinant(size, particles, right, left, value):
    """Independent configuration-space generator, not Bethe roots (small tests only)."""
    configs = [sum(1 << i for i in sites) for sites in itertools.combinations(range(size), particles)]
    indices = {state: i for i, state in enumerate(configs)}
    matrix = [[complex(value if i == j else 0) for j in range(len(configs))] for i in range(len(configs))]
    for col, state in enumerate(configs):
        for source in range(size):
            if not state & (1 << source): continue
            for direction, rate in [(1, right), (-1, left)]:
                target = (source+direction) % size
                if not state & (1 << target):
                    destination = state ^ (1 << source) ^ (1 << target)
                    matrix[indices[destination]][col] -= rate
                    matrix[col][col] += rate
    determinant = 1
    for k in range(len(configs)):
        pivot = max(range(k, len(configs)), key=lambda i: abs(matrix[i][k]))
        if matrix[pivot][k] == 0: return 0
        if pivot != k:
            matrix[pivot], matrix[k] = matrix[k], matrix[pivot]
            determinant = -determinant
        determinant *= matrix[k][k]
        for i in range(k+1, len(configs)):
            factor = matrix[i][k]/matrix[k][k]
            for j in range(k+1, len(configs)):
                matrix[i][j] -= factor*matrix[k][j]
    return determinant


class ExclusionTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_independent_ten_configuration_generator(self):
        for name, gap, frequency in [
                ('tasep-small', .71029006912302867249, .32689551305469573645),
                ('asep-small', 1.04039558904641505636, .15916882161018025320)]:
            row = self.cases[name][1]
            self.assertAlmostEqual(row['gap'], gap, places=12)
            self.assertAlmostEqual(row['frequency'], frequency, places=12)
            _, size, particles, right, left = CASES[name]
            lam = complex(row['lambda_real'], row['lambda_imag'])
            self.assertLess(abs(characteristic_determinant(size, particles, right, left, lam)), 1e-10)
            self.assertGreater(abs(characteristic_determinant(size, particles, right, left, lam+.01)), .01)

    def test_particle_hole_reflection_and_rate_scaling(self):
        for base, other, factor in [('tasep-quarter-l16', 'tasep-holes', 1),
                                    ('tasep-quarter-l16', 'tasep-scaled', 2),
                                    ('asep-bias0.5', 'asep-holes', 1),
                                    ('asep-bias0.5', 'asep-reflected', 1),
                                    ('asep-bias0.5', 'asep-scaled', 2),
                                    ('tasep-quarter-l16', 'asep-bias1', 1)]:
            for key in ('lambda_real', 'lambda_imag', 'gap', 'frequency'):
                self.assertAlmostEqual(factor*self.cases[base][1][key], self.cases[other][1][key], places=12)
        self.assertEqual(len(self.cases['tasep-holes'][1]['roots']), 4)  # Not 12 physical particles.

    def test_analytic_endpoints_and_absent_modes(self):
        for size in SIZES:
            row = self.cases[f'tasep-single-l{size}'][1]
            self.assertAlmostEqual(row['gap'], 2*math.sin(math.pi/size)**2, places=14)
            self.assertAlmostEqual(row['frequency'], math.sin(2*math.pi/size), places=14)
        for name in ('asep-bias0', 'asep-symmetric-half'):
            row = self.cases[name][1]
            self.assertAlmostEqual(row['gap'], 2*math.sin(math.pi/16)**2, places=14)
            self.assertEqual(row['frequency'], 0)
            self.assertEqual(row['roots'], [])
        for name in ('tasep-empty', 'tasep-full', 'asep-empty', 'asep-full'):
            row = self.cases[name][1]
            self.assertFalse(row['has_mode'])
            self.assertNotIn('gap', row)
            self.assertEqual(row['roots'], [])

    def test_scaling_and_drift_claims(self):
        for family in ('half', 'quarter', 'single'):
            gaps = [self.cases[f'tasep-{family}-l{l}'][1]['gap'] for l in SIZES]
            self.assertTrue(all(a > b for a, b in zip(gaps, gaps[1:])))
            z = math.log(gaps[-2]/gaps[-1], 2)
            self.assertLess(abs(z-(2 if family == 'single' else 1.5)), .03)
        for size in SIZES:
            self.assertLess(self.cases[f'tasep-half-l{size}'][1]['frequency'], 1e-12)
        self.assertLess(abs(128*self.cases['tasep-quarter-l128'][1]['frequency']/(2*math.pi)-.5), .01)
        rows = [self.cases[f'asep-bias{b}'][1] for b in BIASES]
        self.assertTrue(all(a['frequency'] > b['frequency'] for a, b in zip(rows, rows[1:])))
        self.assertTrue(all(a['gap'] > b['gap'] > 0 for a, b in zip(rows, rows[1:])))

    def test_bad_modes_and_roots_rejected(self):
        for name in ('tasep-quarter-l16', 'asep-bias0.5'):
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in ('relaxation', 'roots')}
            changes = [('relaxation', '# Outcome: success', '# Outcome: partial'),
                       ('relaxation', 'true,converged', 'false,converged'),
                       ('roots', '# Sites: 16', '# Sites: 15'),
                       ('roots', '\n0,', '\n9,')]
            for table, old, new in changes:
                with self.subTest(name=name, old=old):
                    self.assertIn(old, texts[table])
                    corrupt = dict(texts); corrupt[table] = corrupt[table].replace(old, new, 1)
                    with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            # Change individual cells without depending on platform-specific printed digits.
            for table, column, value in [('relaxation', 'gap', '0'), ('relaxation', 'frequency', '-1'),
                                         ('relaxation', 'residual', '1'), ('roots', 'index', 'nan'),
                                         ('roots', 'Z_real' if name.startswith('tasep') else 'v_real', '42')]:
                lines = texts[table].splitlines()
                header = next(i for i, line in enumerate(lines) if not line.startswith('#'))
                cells = lines[header+1].split(',')
                cells[lines[header].split(',').index(column)] = value
                lines[header+1] = ','.join(cells)
                corrupt = dict(texts); corrupt[table] = '\n'.join(lines)+'\n'
                with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])

    def test_stationary_only_is_not_zero_gap(self):
        for name in ('tasep-empty', 'asep-full'):
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in ('relaxation', 'roots')}
            self.assertIn('false,,,,,', texts['relaxation'])
            texts['relaxation'] = texts['relaxation'].replace('false,,,,,', 'false,0,0,0,0,')
            with self.assertRaises(ValueError): read_exports(texts, CASES[name])


if __name__ == '__main__':
    unittest.main()
