import math
import unittest

from plot_lieb_liniger_tutorial import DATA, load_cases, read_export


class LiebLinigerTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_density_scaling(self):
        a, branches = self.cases['ll-c1-n1']
        b, scaled = self.cases['ll-c2-n2']
        for key, factor in [('Ground energy/length', 8), ('Chemical potential', 4),
                            ('Fermi rapidity Q', 2)]:
            self.assertAlmostEqual(float(b[key]), factor*float(a[key]), places=11)
        for branch in branches:
            for x, y in zip(branches[branch], scaled[branch]):
                self.assertAlmostEqual(y['p'], 2*x['p'], places=12)
                self.assertAlmostEqual(y['energy'], 4*x['energy'], places=10)

    def test_endpoints_symmetry_and_tonks_approach(self):
        for _, branches in self.cases.values():
            self.assertEqual(branches['type-i'][0]['energy'], 0)
            holes = branches['type-ii']
            self.assertEqual(holes[0]['energy'], 0)
            self.assertEqual(holes[-1]['energy'], 0)
            for a, b in zip(holes, reversed(holes)):
                self.assertAlmostEqual(a['energy'], b['energy'], places=10)
        for branch in ('type-i', 'type-ii'):
            curves = [self.cases[f'll-c{c}-n1'][1][branch] for c in (1, 10, 100)]
            for a, b, c in zip(*(rows[1:-1] for rows in curves)):
                p = a['p']
                limit = p*(2*math.pi + (p if branch == 'type-i' else -p))
                self.assertLess(a['energy'], b['energy'])
                self.assertLess(b['energy'], c['energy'])
                self.assertLess(c['energy'], limit)

    def test_invalid_data_rejected(self):
        text = (DATA / 'll-c1-n1.csv').read_text()
        lines = text.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('type-i,0.'))
        bad = [text.replace('# Outcome: success', '# Outcome: failure'),
               text.replace(',converged', ',precision_limit', 1),
               text.replace('# Background status: converged', '# Background status: failed'),
               text.replace('physical inverse length; not reduced modulo 2*pi', 'folded'),
               '\n'.join(lines[:index] + lines[index+1:])]
        for column, value in [(0, 'unknown'), (1, '1'), (3, 'nan'), (6, '-1')]:
            copy = lines.copy()
            row = copy[index].split(','); row[column] = value
            copy[index] = ','.join(row)
            bad.append('\n'.join(copy))
        for candidate in bad:
            with self.assertRaises(ValueError):
                read_export(candidate, 1, 1)


if __name__ == '__main__':
    unittest.main()
