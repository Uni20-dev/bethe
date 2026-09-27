import cmath
import csv
import io
import math
import unittest

from plot_gaudin_tutorial import CASES, DATA, FIELDS, load_cases, read_exports, schema


class GaudinTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.cases = load_cases()

    def energy(self, name): return self.cases[name][1]['energy']
    def values(self, name): return self.cases[name][1]['values']

    def test_pairing_exact_cases_and_diagonal_interaction(self):
        self.assertAlmostEqual(self.energy('rich-two'), -math.sqrt(2), places=12)
        self.assertEqual(self.energy('rich-full'), 8)  # 2 sum epsilon - g M, not 12.
        self.assertEqual(self.energy('rich-empty'), 1)  # A blocked fermion is still present.
        self.assertEqual(self.energy('rich-free'), 1)
        self.assertAlmostEqual((self.energy('rich-weak')-2)/1e-6, -2, delta=1e-5)
        for name, (mode, levels, g, m, blocked) in CASES.items():
            if name.startswith('rich-blocked-g'):
                self.assertAlmostEqual(self.energy(name), 6-g-math.sqrt(9+g*g), places=12)
                self.assertGreater(self.energy(name)-self.energy(name.replace('blocked', 'paired')), 0)

    def test_regular_variables_and_reconstructed_pair_roots(self):
        self.assertAlmostEqual(self.energy('rich-paired-gcollision'), 0, places=12)
        for actual, exact in zip(self.values('rich-paired-gcollision'), (7/9, 2/3, 1/3, 2/9)):
            self.assertAlmostEqual(actual, exact, places=12)
        # Test only away from the collision: this is not public rapidity output.
        for tag, g in [('0.1', .1), ('0.5', .5), ('1', 1), ('2', 2)]:
            name = f'rich-paired-g{tag}'
            total, y = self.energy(name), self.values(name)
            product = -g*total/y[0]
            d = cmath.sqrt(total*total-4*product)
            roots = [(total+d)/2, (total-d)/2]
            for i, root in enumerate(roots):
                residual = 1/g+sum(1/(root-2*j) for j in range(4))-2/(root-roots[1-i])
                self.assertLess(abs(residual), 1e-9)
            for i in range(4):
                self.assertLess(abs(g*sum(1/(2*i-r) for r in roots)-y[i]), 1e-10)
            self.assertEqual(abs(d.imag) > 0, g > 2/3)

    def test_pairing_scale_shift_and_block_order(self):
        base = self.energy('rich-paired-g1')
        self.assertAlmostEqual(self.energy('rich-shifted'), base+20, places=12)
        self.assertAlmostEqual(self.energy('rich-scaled'), 2*base, places=12)
        for name in ('rich-shifted', 'rich-scaled'):
            for a, b in zip(self.values(name), self.values('rich-paired-g1')):
                self.assertAlmostEqual(a, b, places=12)
        self.assertEqual(self.energy('rich-block-order'), self.energy('rich-blocked-g1'))

    def test_central_analytic_spin_reversal_and_permutation(self):
        self.assertAlmostEqual(self.energy('central-two'), -.25-math.sqrt(2)/2, places=12)
        self.assertAlmostEqual(self.energy('central-two-zero'), -.75, places=12)
        for name, exact in [('central-up', 1), ('central-down', 0), ('central-isolated', -.5)]:
            self.assertEqual(self.energy(name), exact)
            self.assertEqual(self.values(name), [])
        self.assertEqual(self.energy('central-reversed'), self.energy('central-s1-b1'))
        self.assertEqual(self.values('central-reversed'), self.values('central-s1-b1'))
        self.assertAlmostEqual(self.energy('central-scaled'), 2*self.energy('central-s0-b1'), places=12)
        self.assertAlmostEqual(self.energy('central-permuted'), self.energy('central-s0-b1'), places=12)
        for a, i in zip(self.values('central-permuted'), (0, 3, 1, 2)):
            self.assertAlmostEqual(a, self.values('central-s0-b1')[i], places=12)
        self.assertLess(self.values('central-s0-b1')[0], 0)  # Not an occupation.

    def test_central_zero_field_multiplet_concavity_and_slopes(self):
        for sz in (-1, 1):
            self.assertAlmostEqual(self.energy(f'central-s{sz}-b0'), self.energy('central-s0-b0'), places=12)
        for sz in (-1, 0, 1):
            slopes = [(self.energy(f'central-s{sz}-b{b}')-self.energy(f'central-s{sz}-b{a}'))/(b-a)
                      for a, b in zip(FIELDS, FIELDS[1:])]
            self.assertTrue(all(-.5 <= s <= .5 for s in slopes))
            self.assertTrue(all(a >= b for a, b in zip(slopes, slopes[1:])))
        coarse = (self.energy('central-deriv-plus')-self.energy('central-deriv-minus'))/.02
        fine = (self.energy('central-deriv-half-plus')-self.energy('central-deriv-half-minus'))/.01
        self.assertLess(abs(coarse-fine), 5e-6)
        self.assertTrue(-.5 < fine < 0)
        # High-field central-down product seed, M=2 largest bath couplings up.
        self.assertLess(abs(self.energy('central-high')-(-50+.5-(1+.7)/2)), .01)

    def test_corrupt_variables_and_unfinished_continuation_rejected(self):
        edits = {
            'rich-blocked-g1': [('states', 'reached_g', '.5'), ('states', 'target_residual', '1'),
                                ('variables', 'blocked', 'true'), ('variables', 'y', 'nan')],
            'central-s0-b1': [('states', 'reached_field', ''), ('states', 'sz', '1'),
                              ('variables', 'v', '0'), ('variables', 'spin', '1')],
            'central-up': [('variables', 'v', '0')],
        }
        for name, mutations in edits.items():
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schema(CASES[name])}
            for table, key, value in mutations:
                lines = texts[table].splitlines()
                reader = csv.DictReader(line for line in lines if not line.startswith('#'))
                rows = list(reader); rows[0][key] = value
                out = io.StringIO(); writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                writer.writeheader(); writer.writerows(rows)
                bad = dict(texts)
                bad[table] = '\n'.join(line for line in lines if line.startswith('#'))+'\n'+out.getvalue()
                with self.subTest(name=name, column=key):
                    with self.assertRaises(ValueError): read_exports(bad, CASES[name])
            for old, new in [('# Outcome: success', '# Outcome: partial'), ('true,converged', 'false,converged')]:
                bad = dict(texts); bad['states'] = texts['states'].replace(old, new, 1)
                with self.assertRaises(ValueError): read_exports(bad, CASES[name])


if __name__ == '__main__': unittest.main()
