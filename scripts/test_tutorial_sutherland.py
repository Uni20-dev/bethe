import math
import unittest

from plot_sutherland_tutorial import DATA, load_cases, read_export


class SutherlandTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_collision_branch_and_boost(self):
        q2 = (math.pi/2)**2
        for exponent in (0, 1, 2):
            metadata, rows = self.cases[f'sutherland-lambda{exponent}']
            states = {r['labels']: r for r in rows}
            self.assertAlmostEqual(float(metadata['Ground energy']), 5*exponent**2*q2, places=12)
            self.assertAlmostEqual(states[(0, 0, 0, 1)]['gap']/q2, 1+3*exponent, places=12)
            self.assertAlmostEqual(states[(1, 1, 1, 1)]['gap']/q2, 4, places=12)
            for labels, row in states.items():
                reflected = states[tuple(-n for n in reversed(labels))]
                self.assertAlmostEqual(row['gap'], reflected['gap'], places=12)
                self.assertEqual(row['momentum_index'], -reflected['momentum_index'])

    def test_window_nesting_and_length_scaling(self):
        small = {r['labels']: r for r in self.cases['sutherland-lambda1'][1]}
        large = {r['labels']: r for r in self.cases['sutherland-window2'][1]}
        long_ring = {r['labels']: r for r in self.cases['sutherland-length8'][1]}
        self.assertEqual(len(small), 15)
        self.assertEqual(len(large), 70)
        self.assertLess(small.keys(), large.keys())
        for labels, row in small.items():
            for key in ('gap', 'energy'):
                self.assertEqual(row[key], large[labels][key])
                self.assertAlmostEqual(row[key], 4*long_ring[labels][key], places=12)
            self.assertAlmostEqual(row['p'], 2*long_ring[labels]['p'], places=12)

    def test_invalid_exact_exports_rejected(self):
        text = (DATA / 'sutherland-lambda1.csv').read_text()
        lines = text.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('0,0 0 0 0,'))
        bad = [text.replace('# Outcome: success', '# Outcome: failure'),
               text.replace('# Status: exact spectral rules', '# Status: failed'),
               text.replace('# States enumerated: 15', '# States enumerated: 14'),
               '\n'.join(lines[:index] + lines[index+1:])]
        for column, value in [(1, '1 0 0 0'), (1, '0 0 0'), (2, 'nan'), (3, '1'), (5, '1')]:
            copy = lines.copy()
            row = copy[index].split(','); row[column] = value
            copy[index] = ','.join(row)
            bad.append('\n'.join(copy))
        # Duplicate a valid state while retaining sequential IDs and row count.
        copy = lines.copy()
        copy[index+1] = '1,' + copy[index].split(',', 1)[1]
        bad.append('\n'.join(copy))
        # A truncated exact-rule row must fail even without a row-status column.
        copy = lines.copy(); copy[index] = copy[index].rsplit(',', 1)[0]
        bad.append('\n'.join(copy))
        for candidate in bad:
            with self.assertRaises(ValueError):
                read_export(candidate, 1, 4, 1)


if __name__ == '__main__':
    unittest.main()
