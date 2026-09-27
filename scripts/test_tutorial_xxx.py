import math
import unittest

from plot_xxx_tutorial import DATA, SIZES, load_cases, read_scan, read_spinons


class XXXTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_spinon_coordinates_and_finite_size_trend(self):
        corrections = []
        for n in SIZES:
            rows = self.cases[f'spinons-{n}']
            self.assertAlmostEqual(rows[0]['k'], math.pi/(2*n), places=14)
            self.assertAlmostEqual(rows[-1]['k'], math.pi-math.pi/(2*n), places=14)
            for a, b in zip(rows, reversed(rows)):
                self.assertAlmostEqual(a['energy'], b['energy'], places=11)
                self.assertAlmostEqual(a['epsilon_inf'], math.pi/2*math.sin(a['k']), places=13)
                expected_p = (math.pi*((n-1)//2)+math.pi/2-a['k']) % (2*math.pi)
                self.assertAlmostEqual(math.cos(a['p']), math.cos(expected_p), places=13)
            corrections.append(max(abs(r['bulk_subtracted_energy']-r['epsilon_inf']) for r in rows))
        self.assertGreater(corrections[0], corrections[1])
        self.assertGreater(corrections[1], corrections[2])

    def test_four_site_reference_and_magnons(self):
        for boundary in ('pbc', 'obc'):
            ground, rows = self.cases[f'{boundary}-4']
            expected = (-2 if boundary == 'pbc' else -.75-math.sqrt(3)/2)
            self.assertAlmostEqual(ground, expected, places=12)
            energies = [-1, 0, 0] if boundary == 'pbc' else sorted(-.25+math.cos(math.pi*m/4) for m in (1, 2, 3))
            for row, e in zip(rows, energies):
                self.assertAlmostEqual(row['energy'], e, places=12)
        self.assertEqual(len(self.cases['pbc-12'][1]), 21)
        self.assertEqual(len(self.cases['obc-12'][1]), 21)
        self.assertEqual(math.comb(12, 5)-math.comb(12, 4), 297)

    def test_failed_or_misreferenced_scan_rejected(self):
        text = (DATA / 'xxx-n12-pbc.csv').read_text()
        ground_line = next(line for line in text.splitlines() if line.startswith('# Ground energy: '))
        for old, new in [('# Outcome: success', '# Outcome: partial'),
                         ('# Ground status: converged', '# Ground status: failed'),
                         ('# Converged candidates: 21', '# Converged candidates: 20'),
                         (ground_line, '# Ground energy: 0'),
                         (',true,converged', ',false,failed')]:
            with self.assertRaises(ValueError):
                read_scan(text.replace(old, new, 1), 12, 'pbc')

    def test_spinon_table_join_and_reference(self):
        texts = {name: (DATA / f'xxx-n15-{name}.csv').read_text() for name in ('states', 'spinons')}
        for column, value in [(0, '999'), (1, '0'), (2, 'nan'), (3, '0')]:
            changed = texts.copy(); lines = changed['spinons'].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            row = lines[i].split(','); row[column] = value; lines[i] = ','.join(row)
            changed['spinons'] = '\n'.join(lines)
            with self.assertRaises(ValueError):
                read_spinons(changed, 15)


if __name__ == '__main__':
    unittest.main()
