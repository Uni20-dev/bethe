import math
import unittest

from plot_multicomponent_tutorial import CASES, DATA, COUPLINGS, load_cases, read_exports, schemas


class MulticomponentTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.cases = load_cases()

    def state(self, name): return self.cases[name][1]['states'][0]
    def energy(self, name): return self.state(name)['energy']

    def test_independent_model_reductions(self):
        for a, b in [('sun-two', 'ref-gy6'), ('sun-singlet', 'ref-ll3'),
                     ('bf-one-fermion', 'ref-ll5'), ('bf-bosons', 'ref-ll5')]:
            self.assertAlmostEqual(self.energy(a), self.energy(b), places=12)
        self.assertAlmostEqual(self.energy('bf-contact'), math.pi**2/2, places=12)
        # Check the spatial roots too, not merely one matching total energy.
        for a, b in zip(self.cases['sun-singlet'][1]['levels'][0], self.cases['ref-ll3'][1]['levels'][0]):
            self.assertAlmostEqual(a, b, places=12)
        for a, b in zip(self.cases['bf-one-fermion'][1]['levels'][0], self.cases['ref-ll5'][1]['levels'][0]):
            self.assertAlmostEqual(a, b, places=12)

    def test_weak_and_strong_limits(self):
        for prefix, slope, limit in [('sun', 6, 80*math.pi**2/27), ('bf', 2.8, 8*math.pi**2/5)]:
            energies = [self.energy(f'{prefix}-c{c}') for c in COUPLINGS]
            self.assertTrue(all(a < b for a, b in zip(energies, energies[1:])))
            self.assertLess(abs((energies[0]-self.energy(prefix+'-free'))/.001-slope), .002)
            self.assertLess(energies[-1], limit)
            self.assertLess(1-energies[-1]/limit, .004)
        self.assertAlmostEqual(self.energy('sun-polarized'), 80*math.pi**2/27, places=12)
        self.assertAlmostEqual(self.energy('bf-fermions'), 8*math.pi**2/5, places=12)
        self.assertGreater(self.energy('bf-c1'), self.energy('bf-one-fermion'))

    def test_component_mapping_nesting_and_permutation(self):
        self.assertEqual([r['nesting_rank'] for r in self.cases['sun-permuted'][1]['components']], [2, None, 0, 1])
        self.assertAlmostEqual(self.energy('sun-imbalanced'), self.energy('sun-permuted'), places=13)
        self.assertEqual(self.cases['sun-imbalanced'][1]['levels'], self.cases['sun-permuted'][1]['levels'])
        self.assertEqual(list(map(len, self.cases['sun-c1'][1]['levels'])), [9, 6, 3])
        self.assertEqual(list(map(len, self.cases['sun-four'][1]['levels'])), [6, 3, 2, 1])
        self.assertEqual(list(map(len, self.cases['bf-c1'][1]['levels'])), [5, 2])
        self.assertEqual(self.cases['bf-bosons'][1]['auxiliary_roots'], [])

    def test_length_scaling_and_free_shells(self):
        for prefix in ('sun', 'bf'):
            self.assertAlmostEqual(self.energy(prefix+'-scaled'), self.energy(prefix+'-c1')/4, places=12)
            for a, b in zip(self.cases[prefix+'-c1'][1]['levels'], self.cases[prefix+'-scaled'][1]['levels']):
                for x, y in zip(a, b): self.assertAlmostEqual(x/2, y, places=12)
            self.assertEqual(self.energy(prefix+'-vacuum'), 0)
            self.assertEqual(self.cases[prefix+'-vacuum'][1]['free_modes'], [])
        self.assertAlmostEqual(self.state('sun-free-shell')['p'], 4*math.pi, places=13)
        self.assertAlmostEqual(self.state('bf-free-shell')['momentum'], 2*math.pi, places=13)
        self.assertAlmostEqual(self.energy('bf-free-shell'), 4*math.pi**2, places=12)
        self.assertEqual([r['mode'] for r in self.cases['bf-free-shell'][1]['free_modes']], [0, 0, 0, 1])

    def test_corrupted_tables_and_incomplete_continuation_rejected(self):
        for name, mutations in [
                ('sun-c1', [('states', 'reached_c', '0.5'), ('states', 'target_residual', '1'),
                            ('roots', 'level', '1'), ('roots', 'rapidity', '42'),
                            ('components', 'nesting_rank', '1'), ('components', 'particles', '4')]),
                ('bf-c1', [('states', 'energy', 'nan'), ('states', 'residual', '1'),
                           ('auxiliary_roots', 'lambda', '42'), ('charge_roots', 'I', '42')])]:
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schemas(CASES[name])[1]}
            for table, column, value in mutations:
                lines = texts[table].splitlines()
                header = next(i for i, line in enumerate(lines) if not line.startswith('#'))
                cells = lines[header+1].split(','); cells[lines[header].split(',').index(column)] = value
                lines[header+1] = ','.join(cells)
                corrupt = dict(texts); corrupt[table] = '\n'.join(lines)+'\n'
                with self.subTest(name=name, column=column):
                    with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            for old, new in [('# Outcome: success', '# Outcome: partial'), ('true,converged', 'false,converged')]:
                corrupt = dict(texts); corrupt['states'] = texts['states'].replace(old, new, 1)
                with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])


if __name__ == '__main__': unittest.main()
