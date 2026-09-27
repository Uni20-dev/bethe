import math
import unittest

from plot_xyz_tutorial import DATA, POINTS, load_cases, read_export


class XYZTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_xy_limit_and_spinon_band(self):
        meta, _ = self.cases['xyz-xy']
        self.assertAlmostEqual(float(meta['Jz']), 0, places=14)
        self.assertAlmostEqual(float(meta['Single-spinon gap']), (float(meta['Jx'])-float(meta['Jy']))/2, places=13)
        self.assertAlmostEqual(float(meta['Spinon maximum energy']), (float(meta['Jx'])+float(meta['Jy']))/2, places=13)
        for meta, branches in self.cases.values():
            m, maximum = (float(meta[k]) for k in ('Single-spinon gap', 'Spinon maximum energy'))
            self.assertGreater(maximum, m)
            self.assertGreater(m, 0)
            for r in branches['spinon', 0, 0]:
                expected = math.hypot(m*math.cos(r['p']), maximum*math.sin(r['p']))
                self.assertAlmostEqual(r['energy'], expected, places=13)
            envelope = branches['two-spinon', 0, 0]
            self.assertAlmostEqual(envelope[0]['lower'], 2*m, places=13)
            self.assertAlmostEqual(envelope[0]['upper'], 2*maximum, places=13)
            for i, row in enumerate(envelope):
                partner = envelope[(i+(POINTS-1)//2) % (POINTS-1)]
                for key in ('lower', 'upper'):
                    self.assertAlmostEqual(row[key], partner[key], places=12)

    def test_bound_copies_and_strict_merger(self):
        meta, branches = self.cases['xyz-eta075']
        self.assertEqual({s for _, s, _ in branches if s}, {1, 2})
        self.assertEqual(len(self.cases['xyz-xy'][1]), 2)  # s=1 exactly at a merger.
        for s in (1, 2):
            plus, minus = branches['bound', s, 1], branches['bound', s, -1]
            for i, row in enumerate(minus):
                self.assertAlmostEqual(row['energy'], plus[(i+(POINTS-1)//2) % (POINTS-1)]['energy'], places=12)
            self.assertLess(plus[0]['energy'], 2*float(meta['Single-spinon gap']))
        self.assertAlmostEqual(branches['bound', 1, 1][0]['energy'], .2437867651629989, places=12)
        self.assertGreater(branches['bound', 1, 1][(POINTS-1)//2]['energy'],
                           branches['two-spinon', 0, 0][(POINTS-1)//2]['upper'])

    def test_exchange_scaling(self):
        original, scaled = (self.cases[name] for name in ('xyz-eta075', 'xyz-exchange2'))
        for key in ('Jx', 'Jy', 'Jz'):
            self.assertEqual(original[0][key], scaled[0][key])
        for branch, rows in original[1].items():
            for a, b in zip(rows, scaled[1][branch]):
                for key in a:
                    factor = 2 if key in ('energy', 'lower', 'upper') else 1
                    self.assertAlmostEqual(b[key], factor*a[key], places=12)

    def test_invalid_rows_rejected(self):
        text = (DATA / 'xyz-eta075.csv').read_text()
        bad = [text.replace('# Outcome: success', '# Outcome: partial'),
               text.replace(',converged', ',precision_limit', 1),
               text.replace('# Selected bound branches: 2', '# Selected bound branches: 3')]
        lines = text.splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith('bound,1,-1,1,'))
        for column, value in [(1, '3'), (2, '1'), (6, '1'), (7, '1'), (8, 'nan'), (9, '0')]:
            changed = lines.copy(); row = changed[i].split(','); row[column] = value
            changed[i] = ','.join(row); bad.append('\n'.join(changed))
        bad.append('\n'.join(lines[:i]+lines[i+1:]))
        for candidate in bad:
            with self.assertRaises(ValueError):
                read_export(candidate, .75, 1)


if __name__ == '__main__':
    unittest.main()
