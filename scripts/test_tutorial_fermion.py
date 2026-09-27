import math
import unittest

from plot_fermion_tutorial import CASES, DATA, COUPLINGS, FILLINGS, load_cases, read_exports, schemas


class FermionTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def state(self, name):
        return self.cases[name][1]['states'][0]

    def test_gaudin_yang_contact_and_coupling_limits(self):
        self.assertAlmostEqual(self.state('gy-two-body')['energy'], math.pi**2/2, places=12)
        root = self.cases['gy-two-body'][1]['charge_roots'][-1]['k']
        self.assertAlmostEqual(root*math.tan(root/2), math.pi/2, places=12)
        energies = [self.state(f'gy-c{c}')['energy'] for c in COUPLINGS]
        self.assertTrue(all(a < b for a, b in zip(energies, energies[1:])))
        self.assertAlmostEqual(energies[0], 4*math.pi**2/9, places=12)
        self.assertLess(abs((energies[1]-energies[0])/.001-3), .001)
        limit = 35*math.pi**2/18
        self.assertLess(energies[-1], limit)
        self.assertLess(1-energies[-1]/limit, .003)
        # Different even-particle periodic shell: not the same finite-ring limit.
        self.assertGreater(self.state('gy-polarized')['energy'], limit)

    def test_rescaling_and_spin_reversal(self):
        self.assertAlmostEqual(self.state('gy-scaled')['energy'], self.state('gy-c1')['energy']/4, places=12)
        for table, key in [('charge_roots', 'k'), ('spin_roots', 'rapidity')]:
            for a, b in zip(self.cases['gy-c1'][1][table], self.cases['gy-scaled'][1][table]):
                self.assertAlmostEqual(a[key]/2, b[key], places=12)
        for prefix in ('gy', 'tj'):
            a, b = self.state(prefix+'-imbalanced'), self.state(prefix+'-reversed')
            self.assertAlmostEqual(a['energy'], b['energy'], places=12)
            self.assertEqual(a['sz'], -b['sz'])
            self.assertEqual(a['p'], b['p'])

    def test_free_shell_current_is_not_reduced_modulo_two_pi(self):
        row = self.state('gy-free-shell')
        self.assertEqual(row['momentum_index'], 2)
        self.assertAlmostEqual(row['p'], 4*math.pi, places=14)
        for table in ('free_up', 'free_down'):
            self.assertEqual([r['mode'] for r in self.cases['gy-free-shell'][1][table]], [0, 1])
        self.assertEqual(self.state('gy-polarized')['momentum_index'], 3)

    def test_tj_analytic_states_and_fermionic_translation(self):
        self.assertAlmostEqual(self.state('tj-imbalanced')['energy'], -3-math.sqrt(5), places=12)
        self.assertAlmostEqual(self.state('tj-balanced-n2')['energy'], -4, places=12)
        self.assertAlmostEqual(self.state('tj-four-noholes')['energy'], -6, places=12)
        self.assertAlmostEqual(self.state('tj-six-noholes')['energy'], -5-math.sqrt(13), places=12)
        self.assertAlmostEqual(self.state('tj-four-noholes')['p'], math.pi, places=14)
        self.assertEqual(self.state('tj-six-noholes')['p'], 0)
        self.assertEqual(self.state('tj-balanced-n0')['energy'], 0)
        self.assertAlmostEqual(self.state('tj-polarized-n16')['energy'], 0, places=12)
        for n in FILLINGS[1:]:
            self.assertLess(self.state(f'tj-balanced-n{n}')['energy'], self.state(f'tj-polarized-n{n}')['energy'])
        tables = self.cases['tj-balanced-n14'][1]
        self.assertEqual(len(tables['first_roots']), 9)
        self.assertEqual(len(tables['second_roots']), 2)
        self.assertEqual(self.cases['tj-four-noholes'][1]['second_roots'], [])

    def test_incomplete_or_corrupt_nested_tables_rejected(self):
        for name, root_table, key in [('gy-c1', 'charge_roots', 'k'), ('tj-balanced-n14', 'first_roots', 'rapidity')]:
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schemas(CASES[name])[1]}
            for table, old, new in [('states', '# Outcome: success', '# Outcome: partial'),
                                    ('states', 'true,converged', 'false,converged'),
                                    (root_table, '# N_up:', '# Wrong population:')]:
                corrupt = dict(texts); corrupt[table] = corrupt[table].replace(old, new, 1)
                with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            for table, column, value in [('states', 'energy', 'nan'), ('states', 'residual', '1'),
                                         (root_table, key, '42'), (root_table, 'quantum_number', '42')]:
                lines = texts[table].splitlines()
                header = next(i for i, line in enumerate(lines) if not line.startswith('#'))
                cells = lines[header+1].split(',')
                cells[lines[header].split(',').index(column)] = value
                lines[header+1] = ','.join(cells)
                corrupt = dict(texts); corrupt[table] = '\n'.join(lines)+'\n'
                with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            lines = texts[root_table].splitlines()
            header = next(i for i, line in enumerate(lines) if not line.startswith('#'))
            del lines[header+1]
            corrupt = dict(texts); corrupt[root_table] = '\n'.join(lines)+'\n'
            with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])


if __name__ == '__main__': unittest.main()
